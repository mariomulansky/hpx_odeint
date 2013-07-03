// Copyright 2013 Mario Mulansky
#ifndef DATAFLOW_SHARED_ALGEBRA_HPP
#define DATAFLOW_SHARED_ALGEBRA_HPP

#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/util/unwrapped.hpp>

using hpx::lcos::local::dataflow;
using hpx::lcos::future;
using hpx::lcos::future_traits;
using hpx::util::unwrapped;

template< typename Algebra >
struct local_dataflow_algebra
{
    Algebra m_algebra;

    local_dataflow_algebra( Algebra a = Algebra() )
        : m_algebra( a )
    { }

    //template< class S1 , class S2 , class S3 , class Op >
    // for now just a single state  type
    template< typename S , typename Op >
    void for_each3( S &s1 , const S &s2 , const S &s3 , Op op )
    {
        typedef typename future_traits<S>::value_type state_type;
        s1 = dataflow( hpx::launch::sync , 
                       unwrapped( [this,op]( state_type x1 , state_type x2 , state_type x3 ) -> state_type  
                                  {
                                      //hpx::cout << "for_each3\n" << hpx::flush;
                                      this->m_algebra.for_each3( x1 , x2 , x3 , op );
                                      return x1;
                                  } ) ,
                       s1 , s2 , s3 );
    }

    template< typename S , typename Op >
    void for_each4( S &s1 , const S &s2 , const S &s3 , const S &s4 , Op op )
    {
        typedef typename future_traits<S>::value_type state_type;
        s1 = dataflow( hpx::launch::sync , 
                       unwrapped( [this,op]( state_type x1 , state_type x2 , state_type x3 , state_type x4 ) -> state_type  
                                  {
                                      //hpx::cout << "for_each4\n" << hpx::flush;
                                      this->m_algebra.for_each4( x1 , x2 , x3 , x4 , op );
                                      return x1;
                                  } ) ,
                       s1 , s2 , s3 , s4 );
    }

    template< typename S , typename Op >
    void for_each5( S &s1 , const S &s2 , const S &s3 , const S &s4 , const S &s5 , Op op )
    {
        typedef typename future_traits<S>::value_type state_type;
        s1 = dataflow( hpx::launch::sync , 
                       unwrapped( [this,op]( state_type x1 , state_type x2 , 
                                             state_type x3 , state_type x4 ,
                                             state_type x5 ) -> state_type  
                                  {
                                      //hpx::cout << "for_each5\n" << hpx::flush;
                                      this->m_algebra.for_each5( x1 , x2 , x3 , x4 , x5 , op );
                                      return x1;
                                  } ) ,
                       s1 , s2 , s3 , s4 , s5 );
    }

    template< typename S , typename Op >
    void for_each6( S &s1 , const S &s2 , const S &s3 , const S &s4 , const S &s5 , const S &s6 , Op op )
    {
        typedef typename future_traits<S>::value_type state_type;
        s1 = dataflow( hpx::launch::sync , 
                       unwrapped( [this,op]( state_type x1 , state_type x2 , 
                                             state_type x3 , state_type x4 ,
                                             state_type x5 , state_type x6 ) -> state_type  
                                  {
                                      //hpx::cout << "for_each6\n" << hpx::flush;
                                      this->m_algebra.for_each6( x1 , x2 , x3 , x4 , x5 , x6 , op );
                                      return x1;
                                  } ) ,
                       s1 , s2 , s3 , s4 , s5 , s6 );
    }
};

#endif
