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

using hpx::lcos::future;
using hpx::find_here;
using hpx::lcos::wait;
using hpx::make_ready_future;
using hpx::lcos::local::dataflow;
using hpx::async;
using hpx::util::unwrapped;

using boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan;

typedef std::vector< double > dvec;
typedef std::shared_ptr< dvec > shared_vec;
typedef std::vector< future< shared_vec > > state_type;

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


    std::clog << "Dimension: " << N << ", number of elements per dataflow: " << G;
    std::clog << ", number of dataflow: " << M << ", steps: " << steps << ", dt: " << dt << std::endl;

    dvec p_init( N );

    std::uniform_real_distribution<double> distribution( -1.0 , 1.0 );
    std::mt19937 engine( 0 ); // Mersenne twister MT19937
    auto generator = std::bind(distribution, engine);

    std::generate( p_init.begin() , 
                   p_init.end() , 
                   std::ref(generator) );

    state_type q_in( M );
    state_type p_in( M );

    for( size_t i=0 ; i<M ; ++i )
    {
        q_in[i] = make_ready_future( std::make_shared<dvec>( ) );
        q_in[i] = dataflow( unwrapped(initialize_zero( G )) , q_in[i] );
        p_in[i] = make_ready_future( std::make_shared<dvec>( ) );
        p_in[i] = dataflow( unwrapped(initialize_copy( p_init , i*G , G )) , p_in[i] );
    }

    std::clog << "init dataflows ready" << std::endl;

    wait( q_in );
    wait( p_in );
    std::clog.precision(10);
    std::clog << "Initialization complete, energy: " << energy( q_in , p_in ) << std::endl;

    hpx::util::high_resolution_timer timer;

    integrate_n_steps( stepper_type() , osc_chain , 
                       std::make_pair( boost::ref(q_in) , boost::ref(p_in) ) ,
                       0.0 , dt , steps );

    hpx::cout << "dataflow generation ready\n" << hpx::flush;

    wait( q_in );
    wait( p_in );

    hpx::cout << (boost::format("runtime: %fs\n") %timer.elapsed()) << hpx::flush;

    std::clog << "Integration complete, energy: " << energy( q_in , p_in ) << std::endl;

    std::cout.precision(10);

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
