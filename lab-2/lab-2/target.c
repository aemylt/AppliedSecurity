#include "target.h"

int t;

bool match( char* x, char* y ) {
  int l_x = strlen( x );
  int l_y = strlen( y );

  t = 0;

  if( l_x != l_y ) {
    return false;
  }

  t = 1;

  for( int i = 0; i < l_x; i++, t++ ) {
    if( x[ i ] != y[ i ] ) {
      return false;
    }
  }

  return true;
}

int main( int argc, char* argv[] ) {
  char G[ 80 ], P[] = "password";

  while( true ) {
    fscanf( stdin, "%80s", G );

    if( feof( stdin ) ) {
      break;
    }

    int r = match( G, P );

    fprintf( stdout, "%d\n", t );
    fprintf( stdout, "%d\n", r );

    fflush( stdout );
    fflush( stderr );
  }

  return 0;
}
