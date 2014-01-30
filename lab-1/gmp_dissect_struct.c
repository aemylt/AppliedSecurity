#include "gmp_dissect_struct.h"

int main( int argc, char* argv[] ) {
  mpz_t x;

  mpz_init( x );

  gmp_scanf( "%Zd", x );

  size_t n = abs( x->_mp_size );

  mp_limb_t* t = x->_mp_d;



  for( int i = 0; i < n; i++ ) {
    if( i != 0 ) {
      printf( "+" );
    }

    gmp_printf( "%llu*(2^(64))^(%d)", t[ i ], i );
  }

  printf( "\n" );

  mpz_clear( x );

  return 0;
}

