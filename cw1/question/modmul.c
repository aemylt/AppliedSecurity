#include "modmul.h"

// Inline function to test if the i-th bit of y is 1 or 0.
// Does this by finding the (i / 64)-th limb of y and selecting the (i mod 64)-th bit from there.
inline int test_bit(mpz_t y, int64_t i) {
  return (y->_mp_size > i >> 6) ? (y->_mp_d[i >> 6] >> (i & 63)) & 1 : 0;
}

inline int compare(mpz_t x, mpz_t y) {
  if (x->_mp_size > y->_mp_size) return 1;
  else if (y->_mp_size > x->_mp_size) return -1;
  else return mpn_cmp(x->_mp_d, y->_mp_d, x->_mp_size);
}

// Create variables rho and omega for Montgomery methods
void mont_init_cube(mpz_t rho2, mpz_t rho3, uint64_t *omega, mpz_t N) {
  uint64_t i, tmp, result;
  tmp = 1;
  for (i = 1; i < 64; i++) {
    tmp *= tmp;
    tmp *= N->_mp_d[0];
  }
  *omega = 0 - tmp;
  mpz_set_ui(rho2, 1);
  _mpz_realloc(rho2, (N->_mp_size << 1) + 1);
  for (i = 1; i <= N->_mp_size << 7; i++) {
    result = mpn_lshift(rho2->_mp_d, rho2->_mp_d, rho2->_mp_size, 1);
    if (result) {
      rho2->_mp_d[rho2->_mp_size] = result;
      rho2->_mp_size++;
    }
    if (compare(rho2, N) > -1) mpz_sub(rho2, rho2, N);
  }
  mpz_set(rho3, rho2);
  for (i = 1; i <= N->_mp_size << 6; i++) {
    result = mpn_lshift(rho3->_mp_d, rho3->_mp_d, rho3->_mp_size, 1);
    if (result) {
      rho3->_mp_d[rho3->_mp_size] = result;
      rho3->_mp_size++;
    }
    if (compare(rho3, N) > -1) mpz_sub(rho3, rho3, N);
  }
}

// Create variables rho and omega for Montgomery methods
void mont_init(mpz_t rho2, uint64_t *omega, mpz_t N) {
  uint64_t i, tmp, result;
  tmp = 1;
  for (i = 1; i < 64; i++) {
    tmp *= tmp;
    tmp *= N->_mp_d[0];
  }
  *omega = 0 - tmp;
  mpz_set_ui(rho2, 1);
  _mpz_realloc(rho2, (N->_mp_size << 1) + 1);
  for (i = 1; i <= N->_mp_size << 7; i++) {
    result = mpn_lshift(rho2->_mp_d, rho2->_mp_d, rho2->_mp_size, 1);
    if (result) {
      rho2->_mp_d[rho2->_mp_size] = result;
      rho2->_mp_size++;
    }
    if (compare(rho2, N) > -1) mpz_sub(rho2, rho2, N);
  }
}

// Perform modular multiplication via Montgomery reduction
void mont_mul(mpz_t r, mpz_t x, mpz_t y, mpz_t N, uint64_t omega) {
  mpz_t tmp, t, yi_x;
  uint64_t i, u;
  mpz_init_set_ui(tmp, 0);
  mpz_init(yi_x);
  mpz_init(t);
  for (i = 0; i < N->_mp_size; i++) {
    u = (x->_mp_d[0] * ((y->_mp_size > i) ? y->_mp_d[i] : 0) + tmp->_mp_d[0]) * omega;
    mpz_mul_ui(t, N, u);
    mpz_mul_ui(yi_x, x, (y->_mp_size > i) ? y->_mp_d[i] : 0);
    mpz_add(tmp, tmp, yi_x);
    mpz_add(tmp, tmp, t);
    mpn_rshift(tmp->_mp_d, tmp->_mp_d, tmp->_mp_size, 32);
    mpn_rshift(tmp->_mp_d, tmp->_mp_d, tmp->_mp_size, 32);
    tmp->_mp_size -= 1;
  }
  if (compare(tmp, N) > -1) mpz_sub(tmp, tmp, N);
  mpz_set(r, tmp);
}

// Perform a modular operation via Montgomery Reduction
void mont_red(mpz_t x, mpz_t N, uint64_t omega) {
  uint64_t i, u;
  mpz_t tmp;
  mpz_init(tmp);
  for (i = 0; i < N->_mp_size; i++) {
    u = x->_mp_d[0] * omega;
    mpz_mul_ui(tmp, N, u);
    mpz_add(x, x, tmp);
    mpn_rshift(x->_mp_d, x->_mp_d, x->_mp_size, 32);
    mpn_rshift(x->_mp_d, x->_mp_d, x->_mp_size, 32);
    x->_mp_size -= 1;
  }
  if (compare(x, N) > -1) mpz_sub(x, x, N);
}

void exp_mod_mont(mpz_t r, mpz_t x, mpz_t y, mpz_t N, mpz_t rho2, uint64_t omega) {
  mpz_t tmp;
  int64_t i, j, l, i_digit, i_bit, l_digit, l_bit;
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
      i_digit = i >> 6;
      i_bit = i & 63;
      l_digit = l >> 6;
      l_bit = l & 63;
      if (i_digit == l_digit) u = (y->_mp_d[i_digit] << (63 - i_bit)) >> (63 - i_bit + l_bit);
      else u = ( y->_mp_d[l_digit] >> l_bit) | (((y->_mp_d[i_digit] << (63 - i_bit)) >> (63 - i_bit)) << (64 - l_bit));
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
void exp_mod_cube(mpz_t r, mpz_t x, mpz_t y, mpz_t N) {
  mpz_t rho2, rho3, x_tmp;
  uint64_t omega;

  mpz_init_set(x_tmp, x);

  mpz_init(rho2);
  mpz_init(rho3);

  mont_init_cube(rho2, rho3, &omega, N);
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

  mpz_t N, e, m, c, rho2;
  uint64_t omega;

  // Initialise integers
  mpz_init(N);
  mpz_init(e);
  mpz_init(m);
  mpz_init(c);
  mpz_init(rho2);

  // Repeat until we reach end of stream
  while(gmp_scanf("%ZX", N) > 0) {

    gmp_scanf("%ZX", e);
    gmp_scanf("%ZX", m);

    // c = m^e (mod N)
    mont_init(rho2, &omega, N);
    mont_mul(m, m, rho2, N, omega);
    exp_mod_mont(c, m, e, N, rho2, omega);
    mont_red(c, N, omega);

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
    exp_mod_cube(m_p, c, d_p, p);
    // m_q = c^d (mod q)
    exp_mod_cube(m_q, c, d_q, q);

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

  mpz_t p, q, g, h, m, c_1, c_2, w, rho2;
  uint64_t omega;
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
    mont_init(rho2, &omega, p);
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

  mpz_t p, q, g, x, c_1, c_2, m, tmp, rho2;
  uint64_t omega;

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

  // Repeat until we reach end of stream
  while (gmp_scanf("%ZX", p) > 0) {

    gmp_scanf("%ZX", q);
    gmp_scanf("%ZX", g);
    gmp_scanf("%ZX", x);
    gmp_scanf("%ZX", c_1);
    gmp_scanf("%ZX", c_2);

    // m = c_1^(q-x) * c_2 (mod p)
    mont_init(rho2, &omega, p);
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
