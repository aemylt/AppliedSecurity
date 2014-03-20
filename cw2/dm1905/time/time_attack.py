import sys
import subprocess
import random

# Defining the word length
w = 64
b = 1 << w

# Helper function for ceiling division
# Works by upside-down floor division
def ceildiv(x, y):
    return -(-x//y)

# Helper function for interacting with the system
# Converts an integer into a hexadecimal string, does any necessary padding
# and returns the result as a string.
def interact(c):
    target_in.write("%X\n" % c)
    target_in.flush()

    l = int(target_out.readline().strip())
    m = int(target_out.readline().strip(), 16)
    return l, m

# Helper function to initialise Montgomery variables
def mont_init(N):
    l_N = ceildiv(len("%X" % N), 16)
    rho_2 = 1
    for i in range(2 * l_N * w):
        rho_2 += rho_2
        if rho_2 >= N:
            rho_2 -= N
    omega = 1
    for i in range(1, w):
        omega = (omega * omega * N) & (b - 1)
    omega = b - omega
    return l_N, rho_2, omega

# Perform Montgomery Multiplication.
# Return the computed result and whether a reduction was required or not
def mont_mul(x, y, N, l_N, omega):
    r = 0
    for i in range(l_N):
        u = (((r & (b - 1)) + ((y >> (i * w)) & (b - 1)) * (x & (b - 1))) * omega) & (b - 1)
        r = (r + ((y >> (i * w)) & (b - 1)) * x + u * N) >> w
    if r > N:
        return r - N, True
    else:
        return r, False

# Perform Montgomery Reduction
def mont_red(t, N, l_N, omega):
    for i in range(l_N):
        u = (((t >> (i * w)) & (b - 1)) * omega) & (b - 1)
        t += (u * N) << (i * w)
    t >>= (w * l_N)
    if t > N:
        return t - N
    else:
        return t

# Montgomery exponentiation within a L2R binary loop
def mont_exp_loop(x, y, N, rho_2, l_N, omega):
    t = mont_red(rho_2, N, l_N, omega)

    for i in range(len("%X" % y) * 4 - 1, -1, -1):
        t, _ = mont_mul(t, t, N, l_N, omega)
        if (y >> i) & 1:
            t, _ = mont_mul(t, x, N, l_N, omega)
    return t

# Generate test ciphertexts
def generate_cs(N, rho_2, l_N, omega, num_interactions):
    c_time = []
    c_message = []
    c_p = []
    c_cur = []

    print "Generating messages"
    for i in range(13000):
        c = random.randint(2, N)
        c_mont, _ = mont_mul(c, rho_2, N, l_N, omega)
        c_p.append(c_mont)
        c_cur.append(c_p[-1])
        time, test_message = interact(c)
        num_interactions += 1
        c_time.append(time)
    return c_time, c_p, c_cur, num_interactions, c, test_message

# Get N and e
public = open(sys.argv[2], 'r')
N_hex = public.readline()
e_hex = public.readline()
public.close()

N = int(N_hex, 16)
e = int(e_hex, 16)

# Set up attack target
target = subprocess.Popen(args   = sys.argv[ 1 ],
                          stdout = subprocess.PIPE, 
                          stdin  = subprocess.PIPE)

target_out = target.stdout
target_in  = target.stdin

l_N, rho_2, omega = mont_init(N)

test_c = 1
d = 1
test_message = 0
num_interactions = 0

# Repeat until we manage to successfully decrypt a message
while(pow(test_c, d, N) != test_message):

    c_time, c_p, c_cur, num_interactions, test_c, test_message = generate_cs(N, rho_2, l_N, omega, num_interactions)

    d = 1
    n = 1
    print "Testing key bits"
    # Test 64 bits
    while n < 64:
        print "Guessing bit %d" % n
        c_0 = []
        c_1 = []
        # Fill in the reduction sets
        bit0_red = []
        bit0_nored = []
        bit1_red = []
        bit1_nored = []
        for i, c in enumerate(c_p):
            ci_0 = c_cur[i]
            ci_0, _ = mont_mul(ci_0, ci_0, N, l_N, omega)
            c_0.append(ci_0)
            ci_0, red0 = mont_mul(ci_0, ci_0, N, l_N, omega)
            if red0:
                bit0_red.append(c_time[i])
            else:
                bit0_nored.append(c_time[i])
            ci_1 = c_cur[i]
            ci_1, _ = mont_mul(ci_1, ci_1, N, l_N, omega)
            ci_1, _ = mont_mul(ci_1, c, N, l_N, omega)
            c_1.append(ci_1)
            ci_1, red1 = mont_mul(ci_1, ci_1, N, l_N, omega)
            if red1:
                bit1_red.append(c_time[i])
            else:
                bit1_nored.append(c_time[i])

        # Calculate the distinguisher
        mean_bit0_red = sum(bit0_red) // len(bit0_red)
        mean_bit0_nored = sum(bit0_nored) // len(bit0_nored)
        mean_bit1_red = sum(bit1_red) // len(bit1_red)
        mean_bit1_nored = sum(bit1_nored) // len(bit1_nored)
        diff_0 = mean_bit0_red - mean_bit0_nored
        diff_1 = mean_bit1_red - mean_bit1_nored

        # Statistically guess a bit
        # If there is no difference, the default currently is 1
        if diff_0 > diff_1:
            d <<= 1
            c_cur = c_0
        else:
            d = (d << 1) + 1
            c_cur = c_1
        n += 1
        # Test if we have recovered the key by decrypting a ciphertext
        # This is done by bruteforcing the final bit, as we can't exploit
        # the square in the next round due to there not being another round
        if pow(test_c, d << 1, N) == test_message:
            d <<= 1
            break
        elif pow(test_c, (d << 1) + 1, N) == test_message:
            d = (d << 1) + 1
            break

# Print results
print "Key Recovered"
print "Key: %X" % d
print "Message decrypted by key: %X" % pow(test_c, d, N)
print "Message decrypted by system: %X" % test_message
print "Number of interactions with system: %d" % num_interactions
