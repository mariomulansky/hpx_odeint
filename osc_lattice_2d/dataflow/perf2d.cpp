// Copyright 2013 Mario Mulansky
//

#include <iostream>
#include <vector>
#include <memory>

#define HPX_LIMIT 6

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/components/dataflow/dataflow.hpp>
#include <hpx/lcos/async.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/foreach.hpp>

#include "dataflow_shared_resize.hpp"
#include "dataflow_shared_algebra.hpp"
#include "dataflow_shared_operations.hpp"
#include "serialization.hpp"
#include "2d_system.hpp"
#include "hpx_odeint_actions.hpp"
#include "spreading_observer.hpp"

using hpx::async;
using hpx::lcos::future;
using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;
using hpx::lcos::wait;

using boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan;

typedef dataflow_base< double > df_base;
typedef std::vector< double > dvec;
typedef std::vector< dvec > dvecvec;
typedef std::shared_ptr< dvecvec > shared_vec;
typedef std::vector< dataflow_base< shared_vec > > state_type;

typedef symplectic_rkn_sb3a_mclachlan< state_type ,
                                       state_type ,
                                       double ,
                                       state_type ,
                                       state_type , 
                                       double ,
                                       dataflow_shared_algebra_2d ,
                                       dataflow_operations > stepper_type;

int hpx_main(boost::program_options::variables_map& vm)
{

    const std::size_t N1 = vm["N1"].as<std::size_t>();
    const std::size_t N2 = vm["N2"].as<std::size_t>();
    const std::size_t G =  vm["G"].as<std::size_t>();
    const bool fully_random = vm["fully_random"].as<bool>();
    const std::size_t init_length = vm["init_length"].as<std::size_t>();
    const std::size_t steps = vm["steps"].as<std::size_t>();
    const double dt = vm["dt"].as<double>();

    const std::size_t M = N1/G;

    double min_time = 1000000.0;
    double mean_time = 0.0;
    
    for( size_t n=0 ; n<10 ; ++n )
    {

        dvecvec p_init( N1 , dvec( N2 , 0.0 ) );

        std::uniform_real_distribution<double> distribution( -1.0 , 1.0 );
        std::mt19937 engine( 0 ); // Mersenne twister MT19937
        auto generator = std::bind(distribution, engine);

        if( fully_random )
        {
            for( size_t j=0 ; j<N1 ; j++ )
                std::generate( p_init[j].begin() , 
                               p_init[j].end() , 
                               std::ref(generator) );
        } else
        {
            for( size_t j=N1/2-init_length/2 ; j<N1/2+init_length/2 ; j++ )
                std::generate( p_init[j].begin()+N2/2-init_length/2 , 
                               p_init[j].begin()+N2/2+init_length/2 , 
                               std::ref(generator) );
        }

        state_type q( M );
        state_type p( M );

        for( size_t i=0 ; i<M ; ++i )
        {
            q[i] = dataflow< initialize_2d_action >( find_here() , 
                                                     std::allocate_shared< dvecvec >( std::allocator<dvecvec>() ) ,
                                                     G ,
                                                     N2 ,
                                                     0.0 );
            p[i] = dataflow< initialize_2d_from_data_action >( find_here() , 
                                                               std::allocate_shared< dvecvec >( std::allocator<dvecvec>() ) ,
                                                               p_init ,
                                                               i*G ,
                                                               G
                                                               );
        }

        std::vector< future<shared_vec> > futures_q( M );
        std::vector< future<shared_vec> > futures_p( M );
        for( size_t i=0 ; i<M ; ++i )
        {
            futures_q[i] = q[i].get_future();
            futures_p[i] = p[i].get_future();
        }

        wait( futures_q );
        wait( futures_p );

        hpx::util::high_resolution_timer timer;

        stepper_type stepper;
        spreading_observer obs;

        for( size_t t=0 ; t<steps ; ++t )
        {
            auto x = std::make_pair( boost::ref(q) , boost::ref(p) );
            stepper.do_step( system_2d , 
                             x ,
                             t*dt , 
                             dt );

        }

        for( size_t i=0 ; i<M ; ++i )
        {
            futures_q[i] = q[i].get_future();
            futures_p[i] = p[i].get_future();
        }
        wait( futures_q );
        wait( futures_p );

        double run_time = timer.elapsed();

        std::clog << "run " << n << ": " << run_time << std::endl;

        min_time = std::min( min_time , run_time );
        mean_time += run_time;
    }
    
    hpx::cout << (boost::format("%d\t%f\t%f\n") % G % (min_time) % (mean_time/10)) << hpx::flush;

    return hpx::finalize();
}


int main( int argc , char* argv[] )
{
    boost::program_options::options_description
       desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

    desc_commandline.add_options()
        ( "N1",
          boost::program_options::value<std::size_t>()->default_value(1024),
          "Dimension 1 (1024)")
        ;
    desc_commandline.add_options()
        ( "N2",
          boost::program_options::value<std::size_t>()->default_value(1024),
          "Dimension 2 (1024)")
        ;
    desc_commandline.add_options()
        ( "G",
          boost::program_options::value<std::size_t>()->default_value(64),
          "granulartity (64)")
        ;
    desc_commandline.add_options()
        ( "fully_random",
          boost::program_options::value<bool>()->default_value(true),
          "Fully random initial condition (true)")
        ;
    desc_commandline.add_options()
        ( "init_length",
          boost::program_options::value<std::size_t>()->default_value(32),
          "Initial excitation length (32)")
        ;
    desc_commandline.add_options()
        ( "steps",
          boost::program_options::value<std::size_t>()->default_value(100),
          "time steps (100)")
        ;
    desc_commandline.add_options()
        ( "dt",
          boost::program_options::value<double>()->default_value(0.1),
          "step size (0.1)")
        ;

    // Initialize and run HPX
    return hpx::init(desc_commandline, argc, argv);
}
