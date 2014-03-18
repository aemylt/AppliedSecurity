import Crypto.Util.number as number
import sys
import subprocess
import random

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

def mont_mul(x, y, N, l_N, omega):
    r = 0
    for i in range(l_N):
        u = (((r & (b - 1)) + ((y >> (i * w)) & (b - 1)) * (x & (b - 1))) * omega) & (b - 1)
        r = (r + ((y >> (i * w)) & (b - 1)) * x + u * N) >> w
    if r > N:
        return r - N, True
    else:
        return r, False

def mont_red(t, N, l_N, omega):
    for i in range(l_N):
        u = (((t >> (i * w)) & (b - 1)) * omega) & (b - 1)
        t += (u * N) << (i * w)
    t >>= (w * l_N)
    if t > N:
        return t - N
    else:
        return t

def mont_exp_loop(x, y, N, rho_2, l_N, omega):
    t = mont_red(rho_2, N, l_N, omega)

    for i in range(len("%X" % y) * 4 - 1, -1, -1):
        t, _ = mont_mul(t, t, N, l_N, omega)
        if (y >> i) & 1:
            t, _ = mont_mul(t, x, N, l_N, omega)
    return t

def mont_exp(x, y, N, rho_2, l_N, omega):
    x_p, _ = mont_mul(x, rho_2, N, l_N, omega)
    return mont_red(mont_exp_loop(x_p, y, N, rho_2, l_N, omega), N, l_N, omega)

def test_message(c, d, N, rho_2, l_N, omega):
    bit0_red = []
    bit0_nored = []
    bit1_red = []
    bit1_nored = []
    c_p, _ = mont_mul(c, rho_2, N, l_N, omega)
    t = mont_exp_loop(c_p, d, N, rho_2, l_N, omega)
    t, _ = mont_mul(t, t, N, l_N, omega)
    t_0 = t
    t_1, _ = mont_mul(t, c_p, N, l_N, omega)
    _, red_0 = mont_mul(t_0, t_0, N, l_N, omega)
    _, red_1 = mont_mul(t_1, t_1, N, l_N, omega)
    if red_0:
        time, _ = interact(c)
        bit0_red.append(time)
    else:
        time, _ = interact(c)
        bit0_nored.append(time)

    if red_1:
        time, _ = interact(c)
        bit1_red.append(time)
    else:
        time, _ = interact(c)
        bit1_nored.append(time)
    return bit0_red, bit0_nored, bit1_red, bit1_nored

# Get N and e
public = open(sys.argv[2], 'r')
N_hex = public.readline()
e_hex = public.readline()
public.close()

N = int(N_hex, 16)
e = int(e_hex, 16)

target = subprocess.Popen(args   = sys.argv[ 1 ],
                          stdout = subprocess.PIPE, 
                          stdin  = subprocess.PIPE)

target_out = target.stdout
target_in  = target.stdin

l_N, rho_2, omega = mont_init(N)

c = random.randint(0, N)

print test_message(c, 1, N, rho_2, l_N, omega)
