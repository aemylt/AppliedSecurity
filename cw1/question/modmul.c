#include "modmul.h"

void exp_mod(mpz_t r, mpz_t x, mpz_t y, mpz_t N) {
  mpz_t tmp;
  mpz_init(tmp);
  mpz_set_ui(tmp, 1);
  int i;

  for (i = 1023; i >= 0; i--) {
    mpz_mul(tmp, tmp, tmp);
    mpz_mod(tmp, tmp, N);
    if (mpz_tstbit(y, i)) {
      mpz_mul(tmp, tmp, x);
      mpz_mod(tmp, tmp, N);
    }
  }
  mpz_set(r, tmp);
}

/*
Perform stage 1:

- read each 3-tuple of N, e and m from stdin,
- compute the RSA encryption c,
- then write the ciphertext c to stdout.
*/

void stage1() {

  mpz_t N, e, m, c;

  // Initialise integers
  mpz_init(N);
  mpz_init(e);
  mpz_init(m);
  mpz_init(c);

  // Repeat until we reach end of stream
  while(gmp_scanf("%ZX", N) > 0) {

    gmp_scanf("%ZX", e);
    gmp_scanf("%ZX", m);

    // c = m^e (mod N)
    exp_mod(c, m, e, N);

    // Print to stdout
    gmp_printf("%ZX\n", c);

  }

}

/*
Perform stage 2:

- read each 9-tuple of N, d, p, q, d_p, d_q, i_p, i_q and c from stdin,
- compute the RSA decryption m,
- then write the plaintext m to stdout.
*/

void stage2() {

  mpz_t N, d, p, q, d_p, d_q, i_p, i_q, c, m_p, m_q, m;

  // Initialise integers
  mpz_init(N);
  mpz_init(d);
  mpz_init(p);
  mpz_init(q);
  mpz_init(d_p);
  mpz_init(d_q);
  mpz_init(i_p);
  mpz_init(i_q);
  mpz_init(c);
  mpz_init(m_p);
  mpz_init(m_q);
  mpz_init(m);

  // Repeat until we reach end of stream
  while (gmp_scanf("%ZX", N) > 0) {

    gmp_scanf("%ZX", d);
    gmp_scanf("%ZX", p);
    gmp_scanf("%ZX", q);
    gmp_scanf("%ZX", d_p);
    gmp_scanf("%ZX", d_q);
    gmp_scanf("%ZX", i_p);
    gmp_scanf("%ZX", i_q);
    gmp_scanf("%ZX", c);

    // m_p = c^d (mod p)
    exp_mod(m_p, c, d_p, p);
    // m_q = c^d (mod q)
    exp_mod(m_q, c, d_q, q);

    // Compute chinese remainder theorem:
    // m = (m_p * q * q^-1 (mod p)) + (m_q * p * p^-1 (mod q)) (mod N)
    mpz_mul(m_p, m_p, q);
    mpz_mul(m_p, m_p, i_q);
    mpz_mul(m_q, m_q, p);
    mpz_mul(m_q, m_q, i_p);
    mpz_add(m, m_p, m_q);
    mpz_mod(m, m, N);

    // Print to stdout
    gmp_printf("%ZX\n", m);

  }

}

/*
Perform stage 3:

- read each 5-tuple of p, q, g, h and m from stdin,
- compute the ElGamal encryption c = (c_1,c_2),
- then write the ciphertext c to stdout.
*/

void stage3() {

  mpz_t p, q, g, h, m, c_1, c_2, w;
#ifndef DEBUG
  gmp_randstate_t state;
#endif

  // Initialise integers
  mpz_init(p);
  mpz_init(q);
  mpz_init(g);
  mpz_init(h);
  mpz_init(m);
  mpz_init(c_1);
  mpz_init(c_2);
  mpz_init(w);

#ifndef DEBUG
  // Initialise random state with a Mersenne Twister algorithm
  gmp_randinit_mt(state);
#else
  // Set random key w to 1
  mpz_set_ui(w, 1);
#endif

  // Repeat until we reach end of stream
  while (gmp_scanf("%ZX", p) > 0) {

    gmp_scanf("%ZX", q);
    gmp_scanf("%ZX", g);
    gmp_scanf("%ZX", h);
    gmp_scanf("%ZX", m);

#ifndef DEBUG
    // Get random value for w
    mpz_urandomm(w, state, q);
#endif

    // c_1 = g^w (mod p)
    // c_2 = m * h^w (mod p)
    exp_mod(c_1, g, w, p);
    exp_mod(c_2, h, w, p);
    mpz_mul(c_2, m, c_2);
    mpz_mod(c_2, c_2, p);

    // Print to stdout
    gmp_printf("%ZX\n", c_1);
    gmp_printf("%ZX\n", c_2);

  }

}

/*
Perform stage 4:

- read each 5-tuple of p, q, g, x and c = (c_1,c_2) from stdin,
- compute the ElGamal decryption m,
- then write the plaintext m to stdout.
*/

void stage4() {

  mpz_t p, q, g, x, c_1, c_2, m, tmp;

  // Initialise integers
  mpz_init(p);
  mpz_init(q);
  mpz_init(g);
  mpz_init(x);
  mpz_init(c_1);
  mpz_init(c_2);
  mpz_init(m);
  mpz_init(tmp);

  // Repeat until we reach end of stream
  while (gmp_scanf("%ZX", p) > 0) {

    gmp_scanf("%ZX", q);
    gmp_scanf("%ZX", g);
    gmp_scanf("%ZX", x);
    gmp_scanf("%ZX", c_1);
    gmp_scanf("%ZX", c_2);

    // m = c_1^(q-x) * c_2 (mod p)
    mpz_sub(tmp, q, x);
    exp_mod(c_1, c_1, tmp, p);
    mpz_mul(m, c_1, c_2);
    mpz_mod(m, m, p);

    // Print to stdout
    gmp_printf("%ZX\n", m);

  }

}

/*
The main function acts as a driver for the assignment by simply invoking
the correct function for the requested stage.
*/

int main( int argc, char* argv[] ) {
  if     ( !strcmp( argv[ 1 ], "stage1" ) ) {
    stage1();
  }
  else if( !strcmp( argv[ 1 ], "stage2" ) ) {
    stage2();
  }
  else if( !strcmp( argv[ 1 ], "stage3" ) ) {
    stage3();
  }
  else if( !strcmp( argv[ 1 ], "stage4" ) ) {
    stage4();
  }

  return 0;
}
