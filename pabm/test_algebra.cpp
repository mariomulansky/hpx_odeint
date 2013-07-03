// Copyright 2013 Mario Mulansky
//

#include <iostream>
#include <vector>
#include <memory>

#define HPX_LIMIT 6

#include <hpx/hpx.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/lcos/async.hpp>
#include <hpx/util/unwrapped.hpp>
#include <hpx/include/iostreams.hpp>

#include <boost/numeric/odeint.hpp>

#include "local_dataflow_algebra.hpp"

using hpx::lcos::future;
using hpx::find_here;
using hpx::lcos::wait;
using hpx::make_ready_future;
using hpx::lcos::local::dataflow;
using hpx::async;
using hpx::util::unwrapped;

typedef future< double > state_type;

using boost::numeric::odeint::vector_space_algebra;
using boost::numeric::odeint::default_operations;

int main()
{

    state_type x = make_ready_future( 1.0 );
    state_type y = make_ready_future( 2.0 );
    state_type z = make_ready_future( 1.0 );

    local_dataflow_algebra< vector_space_algebra > algebra;
    
    x.wait();
    y.wait();
    z.wait();
    
    algebra.for_each3( z , x , y, default_operations::scale_sum2<double>( 1.0 , 1.0 ) );

    std::cout << z.get() << std::endl;

    return 0;
}
