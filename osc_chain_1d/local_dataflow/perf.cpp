// Copyright 2013 Mario Mulansky
//

#include <iostream>
#include <vector>
#include <memory>

#define HPX_LIMIT 6

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/lcos/async.hpp>
#include <hpx/util/unwrapped.hpp>
#include <hpx/include/iostreams.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/foreach.hpp>

#include "local_dataflow_shared_resize.hpp"
#include "local_dataflow_algebra.hpp"
#include "local_dataflow_shared_operations.hpp"
#include "initialize.hpp"
#include "system.hpp"

using hpx::lcos::shared_future;
using hpx::find_here;
using hpx::lcos::wait_all;
using hpx::make_ready_future;
using hpx::lcos::local::dataflow;
using hpx::async;
using hpx::util::unwrapped;

using boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan;
using boost::numeric::odeint::integrate_n_steps;

typedef std::vector< double > dvec;
typedef std::shared_ptr< dvec > shared_vec;
typedef std::vector< shared_future< shared_vec > > state_type;

typedef symplectic_rkn_sb3a_mclachlan< state_type ,
                                       state_type ,
                                       double ,
                                       state_type ,
                                       state_type , 
                                       double ,
                                       local_dataflow_algebra ,
                                       local_dataflow_shared_operations > stepper_type;

int hpx_main(boost::program_options::variables_map& vm)
{


    const std::size_t N = vm["N"].as<std::size_t>();
    const std::size_t G = vm["G"].as<std::size_t>();
    const std::size_t steps = vm["steps"].as<std::size_t>();
    const double dt = vm["dt"].as<double>();
    const std::size_t M = N/G;

    double avrg_time = 0.0;
    double min_time = 1000000.0;

    for( size_t n=0 ; n<12 ; ++n )
    {

        dvec p_init( N );

        std::uniform_real_distribution<double> distribution( -1.0 , 1.0 );
        std::mt19937 engine( 0 ); // Mersenne twister MT19937
        auto generator = std::bind(distribution, engine);

        std::generate( p_init.begin() , 
                       p_init.end() , 
                       std::ref(generator) );

        state_type q( M );
        state_type p( M );

        for( size_t i=0 ; i<M ; ++i )
        {
            q[i] = make_ready_future( std::make_shared<dvec>( ) );
            q[i] = dataflow( unwrapped(initialize_zero( G )) , q[i] );
            p[i] = make_ready_future( std::make_shared<dvec>( ) );
            p[i] = dataflow( unwrapped(initialize_copy( p_init , i*G , G )) , p[i] );
        }

        wait_all( q );
        wait_all( p );

        hpx::util::high_resolution_timer timer;

        integrate_n_steps( stepper_type() , osc_chain , 
                           std::make_pair( boost::ref(q) , boost::ref(p) ) ,
                           0.0 , dt , steps );

        //hpx::cout << "dataflow generation ready\n" << hpx::flush;

        wait_all( q );
        wait_all( p );

        double run_time = timer.elapsed();

        if( n > 1 )
        {
            avrg_time += run_time;
            min_time = std::min( run_time , min_time );
        }

        std::clog << G << ", run: " << n << " run time: " << run_time << std::endl;

    }

    hpx::cout << (boost::format("%d\t%f\t%f\n") % G % min_time % (avrg_time/10)) << hpx::flush;

    return hpx::finalize();
}


int main( int argc , char* argv[] )
{
    boost::program_options::options_description
       desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

    desc_commandline.add_options()
        ( "N",
          boost::program_options::value<std::size_t>()->default_value(1024),
          "Dimension (1024)")
        ;
    desc_commandline.add_options()
        ( "G",
          boost::program_options::value<std::size_t>()->default_value(128),
          "Block size (128)")
        ;
    desc_commandline.add_options()
        ( "steps",
          boost::program_options::value<std::size_t>()->default_value(100),
          "time steps (100)")
        ;
    desc_commandline.add_options()
        ( "dt",
          boost::program_options::value<double>()->default_value(0.01),
          "step size (0.01)")
        ;

    // Initialize and run HPX
    return hpx::init(desc_commandline, argc, argv);
}
