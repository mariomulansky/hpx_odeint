// Copyright 2013 Mario Mulansky
// resizing functionality for odeint
#ifndef FUTURE_RESIZE_HPP
#define FUTURE_RESIZE_HPP

#include <iostream>
#include <vector>

#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resize.hpp>
#include <boost/numeric/odeint/util/same_size.hpp>

#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/util/unwrapped.hpp>

using hpx::lcos::local::dataflow;
using hpx::lcos::future;
using hpx::util::unwrapped;
using hpx::make_ready_future;

namespace boost { 
namespace numeric { 
namespace odeint {

template< typename T >
struct state_wrapper< future<T> >
{
    typedef state_wrapper< future<T> > state_wrapper_type;

    state_wrapper()
    {
        m_v = make_ready_future( T() );
    };

    future<T> m_v;

};

template< typename T , typename A >
struct is_resizeable< future< std::vector< T , A > > >
{
    typedef boost::true_type type;
    const static bool value = type::value;
};

template< typename T , typename A >
struct same_size_impl< future< std::vector< T , A > > , future< std::vector< T , A > > >
{
    typedef future< std::vector< T , A > > future_type;
    static bool same_size( const future_type &x1 ,
                           const future_type &x2 )
    {
        //hpx::cout << "same size ?" << hpx::endl << hpx::flush;
        return ( ( x1.get().size() == x2.get().size() ) );
    }
};

template< typename T , typename A >
struct resize_impl< future< std::vector< T , A > > , future< std::vector< T , A > > >
{
    typedef future< std::vector< T , A > > future_type;
    typedef std::vector< T , A > state_type;
    static void resize( future_type &x1 ,
                        const future_type &x2 )
    {
        //hpx::cout << "resizing..." << hpx::endl;
        // allocate required memory
        x1 = dataflow( hpx::launch::async , 
                       unwrapped( []( state_type v1 , const state_type &v2 ) -> state_type
                                  {
                                      v1.resize( v2.size() );
                                      return v1;
                                  } ) ,
                       x1 , x2 );
        //hpx::cout << (boost::format("resized %d") % (x1.get().size())) << hpx::endl << hpx::flush;
    }
};

} } }

#endif
