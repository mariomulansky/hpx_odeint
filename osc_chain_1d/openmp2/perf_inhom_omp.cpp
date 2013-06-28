// Copyright 2013 Mario Mulansky

#include <iostream>
#include <vector>
#include <random>

#include <omp.h>

#include <boost/numeric/odeint.hpp>
#include <boost/timer/timer.hpp>
#include <boost/foreach.hpp>

#include "system.hpp"
#include "nested_omp_algebra.hpp"
#include "resize.hpp"

using boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan;
using boost::numeric::odeint::range_algebra;

using boost::timer::cpu_timer;
using boost::timer::cpu_times;

typedef std::vector< double > dvec;
typedef std::vector< dvec > state_type;

typedef symplectic_rkn_sb3a_mclachlan< state_type ,
                                       state_type ,
                                       double ,
                                       state_type ,
                                       state_type , 
                                       double ,
                                       nested_omp_algebra< range_algebra > > stepper_type;

const double KAPPA = 3.3;
const double LAMBDA = 4.7;
const double beta = 1.0;

int main( int argc , char* argv[] )
{
    int N = 1024;
    int steps = 100;
    double dt = 0.01;
    if( argc > 1 )
        N = atoi( argv[1] );
    int block_size = N/4;
    if( argc > 2 )
        block_size = atoi( argv[2] );
    if( argc > 3 )
        steps = atoi( argv[3] );
    if( argc > 4 )
        dt = atof( argv[4] );
    int N_init = N;
    if( argc > 5 )
        N_init = atoi( argv[5] );
    

    int M = N/block_size;
    int G = block_size;

    std::cout << "Size: " << N << " with " << block_size << " elements per thread" << " and " << steps << " steps." << std::endl;

    //omp_set_schedule( omp_sched_dynamic , block_size );
    omp_set_schedule( omp_sched_static , 1 );

    double avrg_time = 0.0;
    double min_time = 1000000.0;
        
    for( size_t n=0 ; n<12 ; ++n )
    {

        osc_chain system( KAPPA , LAMBDA , beta );

        // initialize
        state_type p_init( M );

        // fully random
        for( size_t i=0 ; i<M ; i++ )
        {
            int G_ = (i%2 == 0) ? G/2 : 3*G/2;
            p_init[i].resize( G_ );
            std::uniform_real_distribution<double> distribution( 0.0 );
            std::mt19937 engine( i ); // Mersenne twister MT19937
            auto generator = std::bind( distribution , engine );
            std::generate( p_init[i].begin() , p_init[i].end() , generator );
        }
    
        state_type q( M );
        state_type p( M );
    
#pragma omp parallel for schedule( runtime )
        for( size_t i=0 ; i<M ; i++ )
        {
            q[i] = dvec( p_init[i].size() , 0.0 );
            p[i] = p_init[i];
        }

        //std::clog << "# Initial energy: " << system.energy( q , p ) << std::endl;
    
        cpu_timer timer;

        integrate_n_steps( stepper_type() , 
                           system , 
                           std::make_pair( boost::ref(q) , boost::ref(p) ) , 
                           0.0 , dt , steps );

        double run_time = static_cast<double>(timer.elapsed().wall)/(1000*1000*1000);

        //std::clog << "# Final energy: " << system.energy( q , p ) << std::endl;

        if( n > 1 )
        {
            min_time = std::min( min_time , run_time );
            avrg_time += run_time;
        }

        std::clog << "G: " << block_size << ", run " << n << ": " << run_time << std::endl;

    }

    std::cout << block_size << '\t' << min_time << '\t' << avrg_time/(10) << std::endl;

    return 0;
}
