/* stronlgy nonlinear hamiltonian lattice in 2d */

#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <vector>
#include <cmath>
#include <iostream>

#include <omp.h>

#include <boost/math/special_functions/sign.hpp>

typedef std::vector< double > dvec;

namespace checked_math {
    inline double pow( double x , double y )
    {
        if( x==0.0 )
            // 0**y = 0, don't care for y = 0 or NaN
            return 0.0;
        using std::pow;
        using std::abs;
        return pow( abs(x) , y );
    }
}

double signed_pow( double x , double k )
{
    using boost::math::sign;
    return checked_math::pow( x , k ) * sign(x);
}

struct rhs_func {
    const double m_kap;
    const double m_lam;
    
    rhs_func( const double kap , const double lam )
        : m_kap( kap ) , m_lam( lam ) 
    { }

    void operator()( dvec &dpdt , const dvec &q , double q_l , double q_r )
    {
        const size_t N = q.size();
        double coupling_lr = signed_pow( q_l - q[0] , m_lam-1 );
        for( size_t i=0 ; i<N-1 ; ++i )
        {
            dpdt[i] = dpdt[i] = -signed_pow( q[i] , m_kap-1 )
                + coupling_lr;
            coupling_lr = signed_pow( q[i] - q[i+1] , m_lam-1 );
            dpdt[i] -= coupling_lr;
        }
        dpdt[N-1] = dpdt[N-1] = -signed_pow( q[N-1] , m_kap-1 )
                + coupling_lr - signed_pow( q[N-1] - q_r , m_lam-1 );
    }
};

struct osc_chain {

    const double m_beta;
    const double m_kap;
    const double m_lam;
    int m_threads;

    osc_chain( const double kap , const double lam , 
            const double beta )
        : m_kap( kap ) , m_lam( lam ) , m_beta( beta ) , 
          m_threads(0)
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
            rhs_func f( m_kap , m_lam );
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
        using checked_math::pow;
        // q and dpdt are 2d
        const size_t N = q.size();
        double energy = 0.5*pow( q[0][0] , m_lam ) / m_lam;;
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
                        + pow( q[i][j] , m_kap ) / m_kap
                        + pow( q[i][j]-q[i][j+1] , m_lam ) / m_lam;
                }
                energy += p[i][M-1]*p[i][M-1] / 2.0
                    + pow( q[i][M-1] , m_kap ) / m_kap;
                if( i<N-1 )
                    energy += pow( q[i][M-1]-q[i+1][0] , m_lam ) / m_lam;
                else
                    energy += 0.5*pow( q[i][M-1] , m_lam ) / m_lam;
            }
        }
        return energy;
    }
};

#endif
