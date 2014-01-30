#include "gmp_helloworld.h"

int main( int argc, char* argv[] ) {
  mpz_t r, x, y;

  mpz_init( r );
  mpz_init( x );
  mpz_init( y );

  gmp_scanf( "%Zd",  x );
  gmp_scanf( "%Zd",  y );

  mp_size_t size_x = x->_mp_size;
  mp_size_t size_y = y->_mp_size;
  mp_size_t size_r = (size_x > size_y ? size_x : size_y) + 1;

  r->_mp_d[size_r - 1] = mpn_add( r->_mp_d, x->_mp_d, x->_mp_size, y->_mp_d, y->_mp_size );
  if (r->_mp_d[size_r - 1] == 0) size_r--;
  r->_mp_size = size_r;

  gmp_printf( "%Zd\n", r );

  mpz_clear( r );
  mpz_clear( x );
  mpz_clear( y );

  return 0;
}
