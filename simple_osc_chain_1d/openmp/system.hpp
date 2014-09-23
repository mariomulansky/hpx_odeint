/* stronlgy nonlinear hamiltonian chain */

#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <vector>
#include <cmath>
#include <iostream>

#include <omp.h>

#include <boost/math/special_functions/pow.hpp>

typedef std::vector< double > dvec;

using boost::math::pow;

template<int KAPPA, int LAMBDA>
struct rhs_func {
    
    void operator()( dvec &dpdt , const dvec &q , double q_l , double q_r )
    {
        const size_t N = q.size();
        double coupling_lr = pow<LAMBDA-1>( q_l - q[0] );
        for( size_t i=0 ; i<N-1 ; ++i )
        {
            dpdt[i] = dpdt[i] = -pow<KAPPA-1>( q[i] )
                + coupling_lr;
            coupling_lr = pow<LAMBDA-1>( q[i] - q[i+1] );
            dpdt[i] -= coupling_lr;
        }
        dpdt[N-1] = dpdt[N-1] = -pow<KAPPA-1>( q[N-1] )
	    + coupling_lr - pow<LAMBDA-1>( q[N-1] - q_r );
    }
};

template< int KAPPA, int LAMBDA >
struct osc_chain {

    const double m_beta;
    int m_threads;

    osc_chain( const double beta )
        : m_beta( beta ) ,  m_threads(0)
    { }

    template< class StateIn , class StateOut >
    void operator()( const StateIn &q , StateOut &dpdt )
    {
        // std::cout << "system" << std::endl;
        // q and dpdt are 2d
        const int N = q.size();

#ifndef NO_OMP
#pragma omp parallel for schedule(runtime)
#endif	
        for( int i=0 ; i<N ; ++i )
        {
            rhs_func<KAPPA, LAMBDA> f;
            if( i==0 )
                f( dpdt[i] , q[i] , 0.0 , q[i+1][0] );
            else if ( i<N-1 )
                f( dpdt[i] , q[i] , q[i-1][q[i-1].size()-1] , q[i+1][0] );
            else
                f( dpdt[i] , q[i] , q[i-1][q[i-1].size()-1] , 0.0 );
        }
    }

    template< class StateIn >
    double energy( const StateIn &q , const StateIn &p )
    {
        // q and dpdt are 2d
        const size_t N = q.size();
        double energy = 0.5*pow<LAMBDA>( q[0][0] ) / LAMBDA;
#ifndef NO_OMP
#pragma omp parallel
        {
#pragma omp master
            {
                if( m_threads == 0 )
                    m_threads = omp_get_num_threads();
            }

#pragma omp for reduction(+:energy) schedule(runtime)
#endif //NO_OMP
            for( size_t i=0 ; i<N ; ++i )
            {
                const size_t M=q[i].size();
                for( size_t j=0 ; j<M-1 ; ++j )
                {
                    energy += p[i][j]*p[i][j] / 2.0
                        + pow<KAPPA>( q[i][j] ) / KAPPA
                        + pow<LAMBDA>( q[i][j]-q[i][j+1] ) / LAMBDA;
                }
                energy += p[i][M-1]*p[i][M-1] / 2.0
                    + pow<KAPPA>( q[i][M-1] ) / KAPPA;
                if( i<N-1 )
                    energy += pow<LAMBDA>( q[i][M-1]-q[i+1][0] ) / LAMBDA;
                else
                    energy += 0.5*pow<LAMBDA>( q[i][M-1] ) / LAMBDA;
            }
        }
        return energy;
    }
};

#endif
