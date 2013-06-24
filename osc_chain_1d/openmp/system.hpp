/* stronlgy nonlinear hamiltonian lattice in 1d */

#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <vector>
#include <cmath>
#include <iostream>

#include <omp.h>

#include <boost/math/special_functions/sign.hpp>


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

        double coupling_lr( 0.0 );
        int last_i = -1;

#ifndef NO_OMP
#pragma omp parallel for firstprivate( coupling_lr , last_i ) schedule(runtime)
#endif	
        for( int i=0 ; i<N-1 ; ++i )
        {
            // non-continuous execution
            if( i != (last_i+1) )
            {
                if( i>0 )
                    coupling_lr = signed_pow( q[i-1]-q[i] , m_lam-1 );
                else
                    coupling_lr = 0.0;
            }
            dpdt[i] = -signed_pow( q[i] , m_kap-1 )
                + coupling_lr;
            coupling_lr = signed_pow( q[i]-q[i+1] , m_lam-1 );
            dpdt[i] -= coupling_lr;
            last_i = i;
        }
        dpdt[N-1] = -signed_pow( q[N-1] , m_kap-1 )
            + coupling_lr;
    }

    template< class StateIn >
    double energy( const StateIn &q , const StateIn &p )
    {
        using checked_math::pow;
        // q and dpdt are 2d
        const size_t N = q.size();
        double energy = 0.0;
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
            for( size_t i=0 ; i<N-1 ; ++i )
            {
                energy += p[i]*p[i] / 2.0
                    + pow( q[i] , m_kap ) / m_kap
                    + pow( q[i]-q[i+1] , m_lam ) / m_lam;
            }
        }
        return energy;
    }
};

#endif
