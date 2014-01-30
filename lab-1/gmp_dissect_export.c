#include "gmp_dissect_export.h"

int main( int argc, char* argv[] ) {
  mpz_t x;

  mpz_init( x );

  gmp_scanf( "%Zd", x );

  size_t n = mpz_size( x );

  mp_limb_t t[ n ];

  mpz_export( t, NULL, -1, sizeof( mp_limb_t ), -1, 0, x );

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
