//  Copyright (c) 2013 Mario Mulansky
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// This test case demonstrates the issue described in #771: Setting HPX_LIMIT
// gives huge boost MPL errors.

// anything > 10 triggers #771
#define HPX_LIMIT 15

#include <hpx/hpx.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/components/dataflow/dataflow.hpp>
#include <hpx/util/lightweight_test.hpp>
#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/util/unwrapped.hpp>


using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;
//using hpx::lcos::local::dataflow;
using hpx::lcos::future;
using hpx::util::unwrapped;


// define large action
double func(double x1 , double x2 , double x3 , double x4 , double x5 , double x6 , double x7)
{
    return x1;
}

HPX_PLAIN_ACTION(func , large_action);

int main()
{
    // action too big to compile...
    dataflow_base<double> df = dataflow<large_action>(
                                                      find_here() , 1 , 1 , 1 , 1 , 1 , 1 , 1);
    HPX_TEST_EQ(df.get_future().get(), 1);


    future< double > f = make_ready_future( 1.0 );
    f = hpx::lcos::local::dataflow( hpx::launch::sync , 
                                    unwrapped( [](double x1 , double x2 , double x3 , double x4 ,
                                                  double x5 , double x6 , double x7 ) -> double 
                                               { return x1; } ) , f , f , f , f , f , f , f );
    return hpx::util::report_errors();
}
