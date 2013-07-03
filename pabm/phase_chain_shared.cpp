// Copyright Mario Mulansky 2013

#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>

#define HPX_LIMIT 10

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/lcos/wait_all.hpp>
#include <hpx/lcos/async.hpp>
#include <hpx/util/unwrapped.hpp>
#include <hpx/include/iostreams.hpp>

#include <boost/foreach.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread/thread.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <boost/numeric/odeint.hpp>

#include "parallel_adams_bashforth_stepper.hpp"
#include "local_dataflow_algebra_shared.hpp"
#include "future_resize_shared.hpp"

using hpx::lcos::future;
using hpx::find_here;
using hpx::lcos::wait;
using hpx::make_ready_future;
using hpx::lcos::local::dataflow;
using hpx::when_all;
using hpx::async;
using hpx::util::unwrapped;

typedef std::vector<double> dvec;
typedef std::shared_ptr< dvec > shared_vec;
typedef future< shared_vec> state_type;

using boost::numeric::odeint::parallel_adams_bashforth_stepper;
using boost::numeric::odeint::runge_kutta4;
using boost::numeric::odeint::range_algebra;
using boost::numeric::odeint::integrate_n_steps;

const double GAMMA = 1.2;

typedef parallel_adams_bashforth_stepper< 4 , 
                                          state_type , double , state_type , double , 
                                          local_dataflow_algebra< range_algebra >
                                          > pab_stepper_type;

typedef runge_kutta4< state_type , double , state_type , double , 
                      local_dataflow_algebra< range_algebra > > rk_stepper_type;

inline double coupling( const double x )
{
    return sin( x ) - GAMMA * ( 1.0 - cos( x ) );
}


struct rhs_func 
{
    shared_vec operator()( const shared_vec x_ , shared_vec dxdt_ )
    {
        dvec &x = *x_;
        dvec &dxdt = *dxdt_;
        const size_t N = x.size();
        //hpx::cout << boost::format("rhs start %d , %d \n") % x.size() % dxdt.size() << hpx::flush;
        dxdt[0] = coupling( x[1]-x[0] );
        for( size_t i=1 ; i<N-1 ; ++i )
        {
            dxdt[i] = coupling( x[i+1]-x[i] ) + coupling( x[i-1]-x[i] );
        }
        dxdt[N-1] = coupling( x[N-2] - x[N-1] );
        //hpx::cout << "rhs end\n" << hpx::flush;
        return dxdt_;
    }
};

void rhs( const state_type &x , state_type &dxdt , double t )
{
    dxdt = dataflow( hpx::launch::async , unwrapped( rhs_func() ) ,  x , dxdt );
}

int hpx_main(boost::program_options::variables_map& vm)
{
    const size_t N = vm["N"].as<size_t>();
    const std::size_t steps = vm["steps"].as<std::size_t>();
    const double dt = 0.1;

    shared_vec x_init = std::make_shared<dvec>( N );
    std::uniform_real_distribution<double> distribution( 0.0 , 2*3.14159 );
    std::mt19937 engine( 0 ); // Mersenne twister MT19937
    auto generator = std::bind(distribution, engine);
    std::generate( x_init->begin() , x_init->end() , std::ref(generator) );

    state_type x = make_ready_future( x_init );
    state_type x_out = make_ready_future( std::make_shared<dvec>( N ) );

    pab_stepper_type stepper;
    //rk_stepper_type stepper;
    // for some reason this is necessary
    wait( x );
    wait( x_out );

    hpx::cout << (boost::format("%f\n") % ((*(x.get()))[0])) << hpx::flush;

    stepper.do_step( rhs , x , 0.0 , dt );
    wait( x );

    hpx::util::high_resolution_timer timer;
    
    for( size_t t=0 ; t<steps ; ++t )
    {
        stepper.do_step( rhs , x , 0.0 , x_out , dt );
        wait( x );
        wait( x_out );
        for( size_t n=0 ; n<stepper.m_states.size() ; ++n )
            wait( stepper.m_states[n].m_v );
        std::swap( x , x_out );
        if( t%10 == 9 )
            hpx::cout << boost::format( "step %d done\n" ) % (t+1) << hpx::flush;
    }

    wait( x );

    hpx::cout << (boost::format("runtime: %fs\n") %timer.elapsed()) << hpx::flush;

    hpx::cout << (boost::format("%f\n") % ((*(x.get()))[0]) ) << hpx::flush;
    return hpx::finalize();
}


int main( int argc , char* argv[] )
{
    boost::program_options::options_description
       desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

    desc_commandline.add_options()
        ( "N",
          boost::program_options::value<std::size_t>()->default_value(1024),
          "N (1024)")
        ;

    desc_commandline.add_options()
        ( "steps",
          boost::program_options::value<std::size_t>()->default_value(100),
          "Steps (100)")
        ;

    // Initialize and run HPX
    return hpx::init(desc_commandline, argc, argv);
}
