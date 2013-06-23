// Copyright 2013 Mario Mulansky

#include <iostream>
#include <vector>
#include <random>

#include <boost/numeric/odeint.hpp>
#include <boost/timer/timer.hpp>
#include <boost/foreach.hpp>

#include "system.hpp"
#include "nested_omp_algebra.hpp"
#include "resize.hpp"

using boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan;
using boost::numeric::odeint::range_algebra;

using boost::timer::auto_cpu_timer;
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

const double KAPPA = 3.5;
const double LAMBDA = 4.5;
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

    int M = N/block_size;
    int G = block_size;

    std::cout << "Size: " << N << " with " << block_size << " elements per thread" << " and " << steps << " steps." << std::endl;

    //omp_set_schedule( omp_sched_dynamic , block_size );
    omp_set_schedule( omp_sched_static , block_size );

    // initialize
    state_type p_init( M , dvec( G , 0.0 ) );

    // fully random
    for( size_t i=0 ; i<M ; i++ )
    {
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
        q[i] = dvec( G , 0.0 );
        p[i] = p_init[i];
    }

    osc_chain sys( KAPPA , LAMBDA , beta );

    std::cout << "Initial energy: " << sys.energy( q , p ) << std::endl;

    {
    auto_cpu_timer timer( 3 , "%w sec\n");

    integrate_n_steps( stepper_type() , 
                       sys , 
                       std::make_pair( std::ref(q) , std::ref(p) ) , 
                       0.0 , dt , steps );
                       //std::ref(obs) );
    }

    // std::cout << "Time: " << elapsed << std::endl;
    std::cout << "Final energy: " << sys.energy( q , p ) << std::endl;

    return 0;
}
