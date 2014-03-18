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
        t += u * N * (1 << (i * w))
    t >>= (w * l_N)
    if t > N:
        return t - N
    else:
        return t

def mont_exp(x, y, N, rho_2, l_N, omega):
    x_p, _ = mont_mul(x, rho_2, N, l_N, omega)
    t = mont_red(rho_2, N, l_N, omega)

    for i in range(len("%X" % y) * 4 - 1, -1, -1):
        t, _ = mont_mul(t, t, N, l_N, omega)
        if (y >> i) & 1:
            t, _ = mont_mul(t, x_p, N, l_N, omega)
    r = mont_red(t, N, l_N, omega)
    return r

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
test = mont_exp(c, e, N, rho_2, l_N, omega)

print test
print pow(c, e, N)

c = random.randint(0, N)
l, m = interact(c)
