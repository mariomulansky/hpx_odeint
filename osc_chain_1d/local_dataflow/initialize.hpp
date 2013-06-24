#ifndef INITIALIZE_HPP
#define INITIALIZE_HPP

#include <memory>

typedef std::vector< double > dvec;
typedef std::shared_ptr< dvec > shared_vec;

struct initialize_zero
{
    const size_t m_N;

    initialize_zero( const size_t N )
        : m_N( N )
    { }

    shared_vec operator()( shared_vec v ) const
    {
        //hpx::cout << "initializing vector with zero ...\n" << hpx::flush;
        v->resize( m_N );
        std::fill( v->begin() , v->end() , 0.0 );
        return v;
    }
};


struct initialize_copy
{
    const dvec &m_data; // why no reference here?
    const size_t m_index;
    const size_t m_len;

    initialize_copy( const dvec &data , const size_t index , const size_t len )
        : m_data( data ) , m_index( index ) , m_len( len )
    { }

    shared_vec operator()( shared_vec v ) const
    {
        //hpx::cout << boost::format("initializing vector from data at index %d ...\n") % m_index << hpx::flush;
        v->resize( m_len );
        std::copy( &(m_data[m_index]) , &(m_data[m_index+m_len]) , v->begin() );
        return v;
    }
};

#endif
