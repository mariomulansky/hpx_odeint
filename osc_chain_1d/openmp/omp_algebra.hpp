#ifndef OMP_ALGEBRA_HPP
#define OMP_ALGEBRA_HPP

struct omp_algebra
{

    template< class S1 , class S2 , class S3 , class Op >
    void for_each3( S1 &s1 , S2 &s2 , S3 &s3 , Op op )
    {
#pragma omp parallel for schedule(runtime)
        for( size_t i=0 ; i<boost::size(s1) ; ++i )
            op( s1[i] , s2[i] , s3[i] );
    }

};

#endif
