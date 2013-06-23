// Copyright Mario Mulansky 2013

#include <vector>
#include <iostream>

#include <boost/numeric/odeint/util/resize.hpp>

namespace boost { namespace numeric { namespace odeint {

typedef std::vector< double > state_type;

template<>
struct resize_impl< state_type , state_type >
{
    static void resize( state_type &out , const state_type &in )
    {
        size_t N = boost::size( in );
        out.resize( N );
#pragma omp parallel for schedule( runtime )
        for( size_t n=0 ; n<N ; ++n )
        {
            out[n] = 0.0;
        }
    }
};

} } } 
