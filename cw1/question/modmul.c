#include "modmul.h"

// Inline function to test if the i-th bit of y is 1 or 0.
// Does this by finding the (i / 64)-th limb of y and selecting the (i mod 64)-th bit from there.
inline int test_bit(mpz_t y, int64_t i) {
  return (y->_mp_d[i >> 6] >> (i & 63)) & 1;
}

// Create variables rho and omega for Montgomery methods
void mont_init_cube(mpz_t rho2, mpz_t rho3, mpz_t omega, mpz_t N) {
  uint64_t i;
  mpz_t sub_op;
  mpz_init_set_ui(sub_op, 1);
  mpz_mul_2exp(sub_op, sub_op, 64);
  mpz_set_ui(omega, 1);
  for (i = 1; i < 64; i++) {
    mpz_mul(omega, omega, omega);
    mpz_mul(omega, omega, N);
    mpz_set_ui(omega, omega->_mp_d[0]);
  }
  mpz_sub(omega, sub_op, omega);
  mpz_set_ui(rho2, 1);
  for (i = 1; i <= N->_mp_size << 7; i++) {
    mpz_add(rho2, rho2, rho2);
    if (mpz_cmp(rho2, N) > -1) mpz_sub(rho2, rho2, N);
  }
  mpz_set(rho3, rho2);
  for (i = 1; i <= N->_mp_size << 6; i++) {
    mpz_add(rho3, rho3, rho3);
    if (mpz_cmp(rho3, N) > -1) mpz_sub(rho3, rho3, N);
  }
}

// Create variables rho and omega for Montgomery methods
void mont_init(mpz_t rho2, mpz_t omega, mpz_t N) {
  uint64_t i;
  mpz_t sub_op;
  mpz_init_set_ui(sub_op, 1);
  mpz_mul_2exp(sub_op, sub_op, 64);
  mpz_set_ui(omega, 1);
  for (i = 1; i < 64; i++) {
    mpz_mul(omega, omega, omega);
    mpz_mul(omega, omega, N);
    mpz_set_ui(omega, omega->_mp_d[0]);
  }
  mpz_sub(omega, sub_op, omega);
  mpz_set_ui(rho2, 1);
  for (i = 1; i <= N->_mp_size << 7; i++) {
    mpz_add(rho2, rho2, rho2);
    if (mpz_cmp(rho2, N) > -1) mpz_sub(rho2, rho2, N);
  }
}

// Perform modular multiplication via Montgomery reduction
void mont_mul(mpz_t r, mpz_t x, mpz_t y, mpz_t N, mpz_t omega) {
  mpz_t tmp, yi_x, u;
  uint64_t i, j;
  mpz_init_set_ui(tmp, 0);
  mpz_init(u);
  mpz_init(yi_x);
  for (i = 0; i < N->_mp_size; i++) {
    mpz_set_ui(u, x->_mp_d[0]);
    mpz_mul_ui(u, u, (y->_mp_size > i) ? y->_mp_d[i] : 0);
    mpz_add_ui(u, u, tmp->_mp_d[0]);
    mpz_mul(u, u, omega);
    mpz_set_ui(u, u->_mp_d[0]);
    mpz_mul(u, u, N);
    mpz_mul_ui(yi_x, x, (y->_mp_size > i) ? y->_mp_d[i] : 0);
    mpz_add(tmp, tmp, yi_x);
    mpz_add(tmp, tmp, u);
    for (j = 0; j < tmp->_mp_size - 1; j++) {
      tmp->_mp_d[j] = tmp->_mp_d[j + 1];
    }
    tmp->_mp_size -= 1;
  }
  if (mpz_cmp(tmp, N) > -1) mpz_sub(tmp, tmp, N);
  mpz_set(r, tmp);
}

// Perform a modular operation via Montgomery Reduction
void mont_red(mpz_t x, mpz_t N, mpz_t omega) {
  mpz_t u;
  uint64_t i, j;
  mpz_init(u);
  for (i = 0; i < N->_mp_size; i++) {
    mpz_mul_ui(u, omega, x->_mp_d[0]);
    mpz_set_ui(u, u->_mp_d[0]);
    mpz_mul(u, u, N);
    mpz_add(x, x, u);
    for (j = 0; j < x->_mp_size - 1; j++) {
      x->_mp_d[j] = x->_mp_d[j + 1];
    }
    x->_mp_size -= 1;
  }
  if (mpz_cmp(x, N) > -1) mpz_sub(x, x, N);
}

void exp_mod_mont(mpz_t r, mpz_t x, mpz_t y, mpz_t N, mpz_t rho2, mpz_t omega) {
  mpz_t tmp, mp_u, and_op;
  int64_t i, j, l;
  uint64_t u;
  mpz_t *T;

  // Preprocess results for y = 1,3,5..2^k - 1
  T = malloc(sizeof(mpz_t) << (K_BITS - 1));
  mpz_init_set(T[0], x);
  mpz_init(tmp);
  mont_mul(tmp, x, x, N, omega);
  for (i = 1; i < 1 << (K_BITS - 1); i++) {
    mpz_init(T[i]);
    mont_mul(T[i], T[i - 1], tmp, N, omega);
  }

  mpz_init(mp_u);
  mpz_init(and_op);
  mpz_set_ui(tmp, 1);
  mont_mul(tmp, tmp, rho2, N, omega);
  // Set i to the size of y for 64-bit processors.
  i = y->_mp_size << 6;

  while (i >= 0) {
    // If y_i = 1 then find l and u, otherwise set l = i and u = 0.
    if (test_bit(y, i)) {
      l = i - K_BITS + 1;
      if (l < 0) l = 0;
      while (!test_bit(y, l)) l++;
      // Calculate u by shifting right l bits and performing a bitwise AND with 2^(i - l + 1) - 1.
      mpz_fdiv_q_2exp(mp_u, y, l);
      mpz_set_ui(and_op, (1 << (i - l + 1)) - 1);
      mpz_and(mp_u, mp_u, and_op);
      u = mpz_get_ui(mp_u);
    } else {
      l = i;
      u = 0;
    }
    // t = t^(2^(i - l + 1))
    for (j = 0; j < i - l + 1; j++) {
      mont_mul(tmp, tmp, tmp, N, omega);
    }
    if (u != 0) {
      // Multiply by x^((u - 1)/2) (mod N)
      mont_mul(tmp, tmp, T[(u - 1) >> 1], N, omega);
    }
    i = l - 1;
  }
  mpz_set(r, tmp);
}

// Compute r = x^y (mod N) via sliding window.
void exp_mod(mpz_t r, mpz_t x, mpz_t y, mpz_t N) {
  mpz_t rho2, rho3, omega, x_tmp;

  mpz_init(x_tmp);
  mpz_set(x_tmp, x);

  mpz_init(rho2);
  mpz_init(rho3);
  mpz_init(omega);

  mont_init_cube(rho2, rho3, omega, N);
  mont_mul(x_tmp, x_tmp, rho3, N, omega);
  mont_red(x_tmp, N, omega);

  exp_mod_mont(r, x_tmp, y, N, rho2, omega);
  mont_red(r, N, omega);
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

//    exp_mod(m, c_p, d, N);
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

  mpz_t p, q, g, h, m, c_1, c_2, w, rho2, omega;
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
  mpz_init(rho2);
  mpz_init(omega);

#ifndef DEBUG
  mpz_init(w);
  // Initialise random state with a Mersenne Twister algorithm
  gmp_randinit_mt(state);
#else
  // Set random key w to 1
  mpz_init_set_ui(w, 1);
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
    mont_init(rho2, omega, p);
    mont_mul(g, g, rho2, p, omega);
    exp_mod_mont(c_1, g, w, p, rho2, omega);
    mont_red(c_1, p, omega);

    mont_mul(h, h, rho2, p, omega);
    exp_mod_mont(c_2, h, w, p, rho2, omega);
    mont_mul(m, m, rho2, p, omega);
    mont_mul(c_2, c_2, m, p, omega);
    mont_red(c_2, p, omega);

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

  mpz_t p, q, g, x, c_1, c_2, m, tmp, rho2, omega;

  // Initialise integers
  mpz_init(p);
  mpz_init(q);
  mpz_init(g);
  mpz_init(x);
  mpz_init(c_1);
  mpz_init(c_2);
  mpz_init(m);
  mpz_init(tmp);
  mpz_init(rho2);
  mpz_init(omega);

  // Repeat until we reach end of stream
  while (gmp_scanf("%ZX", p) > 0) {

    gmp_scanf("%ZX", q);
    gmp_scanf("%ZX", g);
    gmp_scanf("%ZX", x);
    gmp_scanf("%ZX", c_1);
    gmp_scanf("%ZX", c_2);

    // m = c_1^(q-x) * c_2 (mod p)
    mont_init(rho2, omega, p);
    mpz_sub(tmp, q, x);

    mont_mul(c_1, c_1, rho2, p, omega);
    exp_mod_mont(m, c_1, tmp, p, rho2, omega);
    mont_mul(c_2, c_2, rho2, p, omega);
    mont_mul(m, m, c_2, p, omega);
    mont_red(m, p, omega);

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
