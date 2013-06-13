// Copyright 2013 Mario Mulansky

#ifndef SYSTEM_2D_HPP
#define SYSTEM_2D_HPP

#include <vector>
#include <memory>
#include <cmath>

#include <boost/math/special_functions/pow.hpp>
#include <boost/thread/thread.hpp>

#include <hpx/lcos/future.hpp>
#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/util/unwrap.hpp>

using hpx::lcos::local::dataflow;
using hpx::lcos::future;
using hpx::lcos::wait;
using hpx::util::unwrap;

using boost::math::pow;

typedef std::vector< double > dvec;
typedef std::vector< dvec > dvecvec;
typedef std::shared_ptr< dvecvec > shared_vecvec;
typedef std::vector< future< shared_vec > > state_type;

template< int Kappa , int Lambda >
struct system_first_block
{
    shared_vecvec operator()( shared_vecvec q , const dvec q_d , shared_vecvec dpdt ) const
    {
        //hpx::cout << (boost::format("first block\n") ) << hpx::flush;

        const size_t N = q->size();

        double coupling_lr = 0.0;
        const size_t M = (*q)[0].size();
        dvec coupling_ud( M , 0.0 );
        for( size_t i=0 ; i<N-1 ; ++i )
        {
            for( size_t j=0 ; j<M-1 ; ++j )
            {
                (*dpdt)[i][j] = -pow<Kappa-1>( (*q)[i][j] ) 
                    + coupling_lr + coupling_ud[j];
                coupling_lr = pow<Lambda-1>( (*q)[i][j]-(*q)[i][j+1] );
                coupling_ud[j] = pow<Lambda-1>( (*q)[i][j]-(*q)[i+1][j] );
                (*dpdt)[i][j] -= coupling_lr + coupling_ud[j];
            }
            (*dpdt)[i][M-1] = -pow<Kappa-1>( (*q)[i][M-1] ) 
                + coupling_lr + coupling_ud[M-1];
            coupling_ud[M-1] = pow<Lambda-1>( (*q)[i][M-1]-(*q)[i+1][M-1] );
            coupling_lr = 0.0;
            (*dpdt)[i][M-1] -= coupling_ud[M-1];
        }
        for( size_t j=0 ; j<M-1 ; ++j )
        {
            (*dpdt)[N-1][j] = -pow<Kappa-1>( (*q)[N-1][j] ) 
                + coupling_lr + coupling_ud[j];
            coupling_lr = pow<Lambda-1>( (*q)[N-1][j]-(*q)[N-1][j+1] );
            (*dpdt)[N-1][j] -= coupling_lr + pow<Lambda-1>( (*q)[N-1][j] - q_d[j] );
        }
        (*dpdt)[N-1][M-1] = -pow<Kappa-1>( (*q)[N-1][M-1] ) 
            + coupling_lr + coupling_ud[M-1]
            - pow<Lambda-1>( (*q)[N-1][M-1] - q_d[M-1] );
    
        //hpx::cout << (boost::format("first block done\n") ) << hpx::flush;
        return dpdt;
    }
};

template< int Kappa , int Lambda >
struct system_center_block
{

    shared_vecvec operator() ( shared_vecvec q , const dvec q_u , 
                               const dvec q_d , shared_vecvec dpdt ) const
    {
        //hpx::cout << (boost::format("center block\n") ) << hpx::flush;
 
        const size_t N = q->size();
        const size_t M = (*q)[0].size();

        double coupling_lr = 0.0;
        dvec coupling_ud( M );

        for( size_t j=0 ; j<M-1 ; ++j )
        {
            (*dpdt)[0][j] = -pow<Kappa-1>( (*q)[0][j] ) + coupling_lr
                - pow<Lambda-1>( (*q)[0][j] - q_u[j] );
            coupling_lr = pow<Lambda-1>( (*q)[0][j]-(*q)[0][j+1] );
            coupling_ud[j] = pow<Lambda-1>( (*q)[0][j]-(*q)[1][j] );
            (*dpdt)[0][j] -= coupling_lr + coupling_ud[j];
        }
        (*dpdt)[0][M-1] = -pow<Kappa-1>( (*q)[0][M-1] ) + coupling_lr
            - pow<Lambda-1>( (*q)[0][M-1]-q_u[M-1] );
        coupling_ud[M-1] = pow<Lambda-1>( (*q)[0][M-1]-(*q)[1][M-1] );
        coupling_lr = 0.0;
        (*dpdt)[0][M-1] -= coupling_ud[M-1];


        for( size_t i=1 ; i<N-1 ; ++i )
        {
            for( size_t j=0 ; j<M-1 ; ++j )
            {
                (*dpdt)[i][j] = -pow<Kappa-1>( (*q)[i][j] ) 
                    + coupling_lr + coupling_ud[j];
                coupling_lr = pow<Lambda-1>( (*q)[i][j]-(*q)[i][j+1] );
                coupling_ud[j] = pow<Lambda-1>( (*q)[i][j]-(*q)[i+1][j] );
                (*dpdt)[i][j] -= coupling_lr + coupling_ud[j];
            }
            (*dpdt)[i][M-1] = -pow<Kappa-1>( (*q)[i][M-1] ) 
                + coupling_lr + coupling_ud[M-1];
            coupling_ud[M-1] = pow<Lambda-1>( (*q)[i][M-1]-(*q)[i+1][M-1] );
            coupling_lr = 0.0;
            (*dpdt)[i][M-1] -= coupling_ud[M-1];
        }

        for( size_t j=0 ; j<M-1 ; ++j )
        {
            (*dpdt)[N-1][j] = -pow<Kappa-1>( (*q)[N-1][j] ) 
                + coupling_lr + coupling_ud[j];
            coupling_lr = pow<Lambda-1>( (*q)[N-1][j]-(*q)[N-1][j+1] );
            (*dpdt)[N-1][j] -= coupling_lr + pow<Lambda-1>( (*q)[N-1][j] - q_d[j] );
        }
        (*dpdt)[N-1][M-1] = -pow<Kappa-1>( (*q)[N-1][M-1] ) 
            + coupling_lr + coupling_ud[M-1]
            - pow<Lambda-1>( (*q)[N-1][M-1] - q_d[M-1] );
    
        //hpx::cout << (boost::format("center block done\n") ) << hpx::flush;

        return dpdt;
    }
};

template< int Kappa , int Lambda >
struct system_last_block
{

    typedef shared_vec result_type;

    shared_vecvec operator()( shared_vecvec q , const dvec q_u , shared_vecvec dpdt ) const
    {
        //hpx::cout << (boost::format("last block\n") ) << hpx::flush;

        const size_t N = q->size();
        const size_t M = (*q)[0].size();

        double coupling_lr = 0.0;
        dvec coupling_ud( M );

        //hpx::cout << (boost::format("last block iterating...\n") ) << hpx::flush;

        for( size_t j=0 ; j<M-1 ; ++j )
        {
            (*dpdt)[0][j] = -pow<Kappa-1>( (*q)[0][j] ) + coupling_lr
                - pow<Lambda-1>( (*q)[0][j] - q_u[j] );
            coupling_lr = pow<Lambda-1>( (*q)[0][j]-(*q)[0][j+1] );
            coupling_ud[j] = pow<Lambda-1>( (*q)[0][j]-(*q)[1][j] );
            (*dpdt)[0][j] -= coupling_lr + coupling_ud[j];
        }
        (*dpdt)[0][M-1] = -pow<Kappa-1>( (*q)[0][M-1] ) + coupling_lr
            - pow<Lambda-1>( (*q)[0][M-1]-q_u[M-1] );
        coupling_ud[M-1] = pow<Lambda-1>( (*q)[0][M-1]-(*q)[1][M-1] );
        coupling_lr = 0.0;
        (*dpdt)[0][M-1] -= coupling_ud[M-1];

        for( size_t i=1 ; i<N-1 ; ++i )
        {
            //hpx::cout << (boost::format("row %d\n") % i ) << hpx::flush;
            for( size_t j=0 ; j<M-1 ; ++j )
            {
                //hpx::cout << (boost::format("row %d , col %d \n") % i %j ) << hpx::flush;
                (*dpdt)[i][j] = -pow<Kappa-1>( (*q)[i][j] ) 
                    + coupling_lr + coupling_ud[j];
                coupling_lr = pow<Lambda-1>( (*q)[i][j]-(*q)[i][j+1] );
                coupling_ud[j] = pow<Lambda-1>( (*q)[i][j]-(*q)[i+1][j] );
                (*dpdt)[i][j] -= coupling_lr + coupling_ud[j];
            }
            //hpx::cout << (boost::format("row %d , col %d \n") % i % (M-1) ) << hpx::flush;
            (*dpdt)[i][M-1] = -pow<Kappa-1>( (*q)[i][M-1] ) 
                + coupling_lr + coupling_ud[M-1];
            coupling_ud[M-1] = pow<Lambda-1>( (*q)[i][M-1]-(*q)[i+1][M-1] );
            coupling_lr = 0.0;
            (*dpdt)[i][M-1] -= coupling_ud[M-1];
        }

        for( size_t j=0 ; j<M-1 ; ++j )
        {
            (*dpdt)[N-1][j] = -pow<Kappa-1>( (*q)[N-1][j] ) 
                + coupling_lr + coupling_ud[j];
            coupling_lr = pow<Lambda-1>( (*q)[N-1][j]-(*q)[N-1][j+1] );
            (*dpdt)[N-1][j] -= coupling_lr;
        }
        (*dpdt)[N-1][M-1] = -pow<Kappa-1>( (*q)[N-1][M-1] ) 
            + coupling_lr + coupling_ud[M-1];
    
        //hpx::cout << (boost::format("last block done\n") ) << hpx::flush;

        return dpdt;
    }
};

template< int Kappa , int Lambda >
struct system_2d
{

    void operator()( state_type &q , state_type &dpdt )
    {
        // works on shared data, but coupling data is provided as copy
        const size_t N = q.size();

        //hpx::cout << boost::format("system call size: %d , %d ...\n") % (q.size()) % (dpdt.size()) << hpx::flush;

        //state_type dpdt_(N);
        // first row
        dpdt[0] = dataflow( hpx::launch::async , 
                            unwrap(system_first_block<Kappa,Lambda>()) , 
                            q[0] , 
                            dataflow( hpx::launch::sync , unwrap([](shared_vecvec v) 
            { return (*v)[0]; }) , q[1] ) , 
                            dpdt[0] );
        // middle rows
        for( size_t i=1 ; i<N-1 ; i++ )
            {
                dpdt[i] = dataflow( hpx::launch::async , 
                                    unwrap(system_center_block<Kappa,Lambda>()) , 
                                    q[i] , 
                                    dataflow( hpx::launch::sync , unwrap([](shared_vecvec v) 
                    { return (*v)[v->size()-1]; }) ,
                                              q[i-1] ) , 
                                    dataflow( hpx::launch::sync , unwrap([](shared_vecvec v) 
                    { return (*v)[0]; }) , q[i+1] ) ,
                                    dpdt[i] );
            }
        dpdt[N-1] = dataflow( hpx::launch::async , 
                              unwrap(system_last_block<Kappa,Lambda>()) , 
                              q[N-1] , 
                              dataflow( hpx::launch::sync , unwrap([](shared_vecvec v) 
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


    double energy( const dvecvec &q , const dvecvec &p )
    {
        using std::abs;
        const size_t N = q.size();
        double energy = 0.0;
        for( size_t i=0 ; i<N-1 ; ++i )
        {
            const size_t M = q[i].size();
            for( size_t j=0 ; j<M-1 ; ++j )
            {
                energy += 0.5*p[i][j]*p[i][j] + pow<Kappa>( q[i][j] ) / Kappa
                    + pow<Lambda>( abs(q[i][j]-q[i][j+1]) ) / Lambda
                    + pow<Lambda>( abs(q[i][j]-q[i+1][j]) ) / Lambda;
            }
            energy += 0.5*p[i][M-1]*p[i][M-1] + pow<Kappa>( q[i][M-1] ) / Kappa
                + pow<Lambda>( abs(q[i][M-1]-q[i+1][M-1]) ) / Lambda;
        }
        const size_t M = q[N-1].size();
        for( size_t j=0 ; j<M-1 ; ++j )
        {
            energy += 0.5*p[N-1][j]*p[N-1][j] + pow<Kappa>( q[N-1][j] ) / Kappa
                + pow<Lambda>( abs(q[N-1][j]-q[N-1][j+1]) ) / Lambda;
        }
        energy += 0.5*p[N-1][M-1]*p[N-1][M-1] + pow<Kappa>( q[N-1][M-1] ) / Kappa;
        return energy;
    }


    template< typename S >
    double energy( const S &q_fut , const S &p_fut )
    {
        dvecvec q,p;
        for( size_t i=0 ; i<q_fut.size() ; ++i )
        {
            for( size_t j=0 ; j<(*(q_fut[i].get())).size() ; ++j )
            {
                q.push_back( (*(q_fut[i].get()))[j] );
                p.push_back( (*(p_fut[i].get()))[j] );
            }
        }
        return energy( q , p );
    }
};

#endif
