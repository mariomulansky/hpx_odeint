// Copyright 2013 Mario Mulansky

#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <cmath>

typedef std::vector< double > dvec;

const double KAPPA = 2.5;
const double LAMBDA = 3.5;

int main( int argc , char **argv ) 
{
    int N = 1024*1024;
    if( argc > 1 )
        N = atoi( argv[1] );

    dvec x( N );
    dvec y( N );
    std::uniform_real_distribution<double> distribution( 0.0 );
    std::mt19937 engine( 0 ); // Mersenne twister MT19937
    auto generator = std::bind( distribution , engine );
    std::generate( x.begin() , x.end() , generator );

    for( int n=0 ; n<10 ; n++ )
    {
        std::transform( x.begin() , x.end() , y.begin() ,
                        []( const double x ) -> double 
                        { 
                            const double y = std::pow( x , KAPPA );
                            return std::pow( y , LAMBDA );
                        } );
    }
}
