#include "gmp_sum.h"

int main( int argc, char* argv[] ) {
  mpz_t x, y;
  mp_size_t size_x, size_y;
  int i;

  mpz_init( x );
  mpz_init( y );

  FILE *f = fopen(argv[1], "r");

  size_x = x->_mp_size;

	gmp_fscanf( f, "%Zd",  x );

  for (i = 0; i < 99; i++) {
		gmp_fscanf( f, "%Zd",  y );
		size_y = y->_mp_size;
		size_x = (size_x > size_y ? size_x : size_y) + 1;

		x->_mp_d[size_x - 1] = mpn_add( x->_mp_d, x->_mp_d, x->_mp_size, y->_mp_d, y->_mp_size );
		if (x->_mp_d[size_x - 1] == 0) size_x--;
		x->_mp_size = size_x;
  }

  gmp_printf( "%Zd\n", x );

  fclose(f);

  mpz_clear( x );
  mpz_clear( y );

  return 0;
}
