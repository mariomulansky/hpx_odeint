// Copyright 2013 Mario Mulansky

#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <vector>
#include <memory>
#include <cmath>

#include <boost/math/special_functions/pow.hpp>
#include <boost/thread/thread.hpp>

#include <hpx/lcos/future.hpp>
#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/util/unwrapped.hpp>

using hpx::lcos::local::dataflow;
using hpx::lcos::shared_future;
using hpx::lcos::wait_all;
using hpx::util::unwrapped;

using boost::math::pow;

const int KAPPA = 4;
const int LAMBDA = 6;

typedef std::vector< double > dvec;
typedef std::shared_ptr< dvec > shared_vec;
typedef std::vector< shared_future< shared_vec > > state_type;

struct system_first_block
{
    shared_vec operator()( shared_vec q , const double q_r , shared_vec dpdt ) const
    {
        //hpx::cout << (boost::format("first block\n") ) << hpx::flush;

        const size_t N = q->size();

        double coupling_lr = -pow<LAMBDA-1>( (*q)[0] );
        for( size_t i=0 ; i<N-1 ; ++i )
        {
            (*dpdt)[i] = -pow<KAPA-1>( (*q)[i] ) + coupling_lr;
            coupling_lr = pow<LAMBDA-1>( (*q)[i]-(*q)[i+1] );
            (*dpdt)[i] -= coupling_lr;
        }
        (*dpdt)[N-1] = -pow<KAPPA-1>( (*q)[N-1] ) 
            + coupling_lr - pow<LAMBDA-1>( (*q)[N-1] - q_r );
    
        //hpx::cout << (boost::format("first block done\n") ) << hpx::flush;
        return dpdt;
    }
};


struct system_center_block
{

    shared_vec operator() ( shared_vec q , const double q_l , 
                               const double q_r , shared_vec dpdt ) const
    {
        //hpx::cout << (boost::format("center block\n") ) << hpx::flush;
 
        const size_t N = q->size();
        double coupling_lr = -pow<LAMBDA-1>( (*q)[0] - q_l );
        for( size_t i=0 ; i<N-1 ; ++i )
        {
            (*dpdt)[i] = -pow<KAPPA-1>( (*q)[i] ) + coupling_lr;
            coupling_lr = pow<LAMBDA-1>( (*q)[i]-(*q)[i+1] );
            (*dpdt)[i] -= coupling_lr;
        }
        (*dpdt)[N-1] = -pow<KAPPA-1>( (*q)[N-1] ) 
            + coupling_lr - pow<LAMBDA-1>( (*q)[N-1] - q_r );

        return dpdt;
    }
};


struct system_last_block
{

    shared_vec operator()( shared_vec q , const double q_l , shared_vec dpdt ) const
    {
        //hpx::cout << (boost::format("last block\n") ) << hpx::flush;

        const size_t N = q->size();

        double coupling_lr = -pow<LAMBDA-1>( (*q)[0] - q_l );
        for( size_t i=0 ; i<N-1 ; ++i )
        {
            (*dpdt)[i] = -pow<KAPPA-1>( (*q)[i] ) + coupling_lr;
            coupling_lr = pow<LAMBDA-1>( (*q)[i]-(*q)[i+1] );
            (*dpdt)[i] -= coupling_lr;
        }
        (*dpdt)[N-1] = -pow<KAPPA-1>( (*q)[N-1] ) 
            + coupling_lr - pow<LAMBDA-1>( (*q)[N-1] );

        return dpdt;
    }
};

void osc_chain( state_type &q , state_type &dpdt )
{
    // works on shared data, but coupling data is provided as copy
    const size_t N = q.size();

    //hpx::cout << boost::format("system call size: %d , %d ...\n") % (q.size()) % (dpdt.size()) << hpx::flush;

    //state_type dpdt_(N);
    // first row
    dpdt[0] = dataflow( hpx::launch::async , unwrapped(system_first_block()) , q[0] , 
                        dataflow( hpx::launch::sync , unwrapped([](shared_vec v) 
        { return (*v)[0]; }) , q[1] ) , 
                        dpdt[0] );
    // middle rows
    for( size_t i=1 ; i<N-1 ; i++ )
    {
        dpdt[i] = dataflow( hpx::launch::async , unwrapped(system_center_block()) , q[i] , 
                            dataflow( hpx::launch::sync , unwrapped([](shared_vec v) 
            { return (*v)[v->size()-1]; }) ,
                                      q[i-1] ) , 
                            dataflow( hpx::launch::sync , unwrapped([](shared_vec v) 
            { return (*v)[0]; }) , q[i+1] ) ,
                            dpdt[i] );
    }
    dpdt[N-1] = dataflow( hpx::launch::async , unwrapped(system_last_block()) , q[N-1] , 
                          dataflow( hpx::launch::sync , unwrapped([](shared_vec v) 
        { return (*v)[v->size()-1]; }),
                                    q[N-2] ) ,
                          dpdt[N-1] );

    /*
    // synchronization to make sure q doesn't get changed while dpdt isn't ready
    for( size_t i=1 ; i<N-1 ; i++ )
    {
        q[i] = dataflow( hpx::launch::async , 
                         unwrap([]( shared_vec x , shared_vec sync){ return x; } ) ,
                         q[i] , 
                         dpdt[i] );
    }
    */
}

void osc_chain_gb( state_type &q , state_type &dpdt )
{
    // works on shared data, but coupling data is provided as copy
    const size_t N = q.size();

    // first row
    dpdt[0] = dataflow( hpx::launch::async , unwrapped(system_first_block()) , q[0] , 
                        dataflow( hpx::launch::sync , unwrapped([](shared_vec v) 
        { return (*v)[0]; }) , q[1] ) , 
                        dpdt[0] );
    // middle rows
    for( size_t i=1 ; i<N-1 ; i++ )
    {
        dpdt[i] = dataflow( hpx::launch::async , unwrapped(system_center_block()) , q[i] , 
                            dataflow( hpx::launch::sync , unwrapped([](shared_vec v) 
            { return (*v)[v->size()-1]; }) ,
                                      q[i-1] ) , 
                            dataflow( hpx::launch::sync , unwrapped([](shared_vec v) 
            { return (*v)[0]; }) , q[i+1] ) ,
                            dpdt[i] );
    }
    dpdt[N-1] = dataflow( hpx::launch::async , unwrapped(system_last_block()) , q[N-1] , 
                          dataflow( hpx::launch::sync , unwrapped([](shared_vec v) 
        { return (*v)[v->size()-1]; }), 
                                    q[N-2] ) , 
                          dpdt[N-1] );
    // global barrier
    wait_all( dpdt );

}


double energy( const dvec &q , const dvec &p )
{
    const size_t N = q.size();
    double energy = 0.5*pow<LAMBDA>( q[0] ) / LAMBDA;
    for( size_t i=0 ; i<N-1 ; ++i )
    {
        energy += 0.5*p[i]*p[i] + pow<KAPPA>( q[i] ) / KAPPA
            + pow<LAMBDA>( q[i]-q[i+1] ) / LAMBDA;
    }
    energy += 0.5*p[N-1]*p[N-1] + pow<KAPPA>( q[N-1] ) / KAPPA
        + 0.5*pow<LAMBDA>( q[N-1] ) / LAMBDA;
    return energy;
}

template< typename S >
double energy( const S &q_fut , const S &p_fut )
{
    dvec q,p;
    for( size_t i=0 ; i<q_fut.size() ; ++i )
    {
        q.insert( q.end() , q_fut[i].get()->begin() , q_fut[i].get()->end() );
        p.insert( p.end() , p_fut[i].get()->begin() , p_fut[i].get()->end() );
    }
    return energy( q , p );
}

#endif
