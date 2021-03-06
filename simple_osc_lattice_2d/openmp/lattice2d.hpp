/* stronlgy nonlinear hamiltonian lattice in 2d */

#ifndef LATTICE2D_HPP
#define LATTICE2D_HPP

#include <vector>
#include <cmath>
#include <iostream>

#include <omp.h>
#include <boost/math/special_functions/sign.hpp>
#include <boost/math/special_functions/pow.hpp>

using boost::math::pow;

template< size_t Kappa , size_t Lambda >
struct lattice2d {

    const double m_beta;
    int m_threads;

    lattice2d( const double beta = 1.0 )
        : m_beta( beta ) , m_threads(0)
    { }

    template< class StateIn , class StateOut >
    void operator()( const StateIn &q , StateOut &dpdt )
    {
        // std::cout << "system" << std::endl;
        // q and dpdt are 2d
        const int N = q.size();
        const int M = q[0].size();

        double coupling_lr( 0.0 );
        std::vector<double> coupling_ud( M , 0.0 );
        int last_i = -1;

#ifndef NO_OMP
#pragma omp parallel for firstprivate( coupling_lr , coupling_ud , last_i ) schedule(runtime)
#endif	
        for( int i=0 ; i<N ; ++i )
        {
            // non-continuous execution
            if( i != (last_i+1) )
            {
                // initialize temporaries in each thread
                for( size_t j=0 ; j<M ; ++j )
                {
                    if( i > 0 )
                        coupling_ud[j] = pow<Lambda-1>( q[i-1][j]-q[i][j] );
                    else
                        coupling_ud[j] = 0.0;
                    //std::cout << coupling_ud[j] << std::endl;
                }
                coupling_lr = 0.0;
            }

            // actual work
            for( size_t j=0 ; j<M-1 ; ++j )
            {
                dpdt[i][j] = -pow<Kappa-1>( q[i][j] )
                    + coupling_lr + coupling_ud[j];
                coupling_lr = pow<Lambda-1>( q[i][j]-q[i][j+1] );
                if( i<N-1 )
                    coupling_ud[j] = pow<Lambda-1>( q[i][j]-q[i+1][j] );
                else
                    coupling_ud[j] = 0.0;
                dpdt[i][j] -= coupling_lr + coupling_ud[j];
                //std::cout << dpdt[i][j] << ": " << q[i][j] << std::endl;
            }
            dpdt[i][M-1] = -pow<Lambda-1>( q[i][M-1] )
                + coupling_lr + coupling_ud[M-1];
            coupling_lr = 0.0;
            if( i<N-1 )
                coupling_ud[M-1] = pow<Lambda-1>( q[i][M-1]-q[i+1][M-1] );
            else
                coupling_ud[M-1] = 0.0;
            dpdt[i][M-1] -= coupling_ud[M-1];
            last_i = i;
        }

        // for( int i=0 ; i<N ; ++i )
        // {
        //     for( int j=0 ; j<M ; ++j )
        //     {
        //         std::cout << dpdt[i][j] << '\t';
        //     }
        //     std::cout << std::endl;
        // }

    }

    template< class StateIn >
    double energy( const StateIn &q , const StateIn &p )
    {
        // q and dpdt are 2d
        const size_t N = q.size();
        const size_t M = q[0].size();
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
                for( size_t j=0 ; j<M-1 ; ++j )
                {
                    energy += p[i][j]*p[i][j] / 2.0
                        + pow<Kappa>( q[i][j] ) / Kappa
                        + pow<Lambda>( q[i][j]-q[i][j+1] ) / Lambda
                        + pow<Lambda>( q[i][j]-q[i+1][j] ) / Lambda;
                }
                energy += p[i][M-1]*p[i][M-1] / 2.0
                    + pow<Kappa>( q[i][M-1] ) / Kappa
                    + pow<Lambda>( q[i][M-1]-q[i+1][M-1] ) / Lambda;
            }
#ifndef NO_OMP
        }
#endif
        for( size_t j=0 ; j<M-1 ; ++j )
        {
            energy += p[N-1][j]*p[N-1][j] / 2.0
                + pow<Kappa>( q[N-1][j] ) / Kappa
                + pow<Lambda>( q[N-1][j]-q[N-1][j+1] ) / Lambda;
        }
        energy += p[N-1][M-1]*p[N-1][M-1] / 2.0
            + pow<Kappa>( q[N-1][M-1] ) / Kappa;
        return energy;
    }


//     template< class StateIn , class StateOut >
//     double local_energy( const StateIn &q , const StateIn &p , StateOut &energy )
//     {
//         // q and dpdt are 2d
// 	const int N = q.size();
//         double e = 0.0;
//         int i;
// #ifndef NO_OMP
// #pragma omp parallel for private(i) reduction( + : e )
// #endif
// 	for( i = 0 ; i < N ; ++i )
// 	{
// 	    const int i_l = (i-1+N) % N;
// 	    const int i_r = (i+1) % N;
// 	    for( size_t j = 0 ; j < N ; ++j )
// 	    {
// 		const int j_l = (j-1+N) % N;
// 		const int j_r = (j+1) % N;
// 		energy[i][j] = p[i][j]*p[i][j] / 2.0
//                     + m_omega[i+m_start][j+m_start] * pow<Kappa>( q[i][j] ) / Kappa
// 		    + m_beta * pow<Lambda>( q[i][j] - q[i][j_l] ) / Lambda / 2
// 		    + m_beta * pow<Lambda>( q[i][j] - q[i][j_r] ) / Lambda / 2
// 		    + m_beta * pow<Lambda>( q[i][j] - q[i_l][j] ) / Lambda / 2
// 		    + m_beta * pow<Lambda>( q[i][j] - q[i_r][j] ) / Lambda / 2;
//                 e += energy[i][j];
//             }
//         }
//         //rescale
//         e = 1.0/e;
// #ifndef NO_OMP
// #pragma omp parallel for private(i)
// #endif
// 	for( i = 0 ; i < N ; ++i )
//             for( size_t j = 0 ; j < N ; ++j )
//                 energy[i][j] *= e;
//         return 1.0/e;
//     }


};

#endif
