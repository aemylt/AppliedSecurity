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

cs = []
c_time = []
c_message = []
c_p = []
c_cur = []

for i in range(10000):
    cs.append(random.randint(2, N))
    c_mont, _ = mont_mul(cs[-1], rho_2, N, l_N, omega)
    c_p.append(c_mont)
    c_cur.append(c_p[-1])
    time, message = interact(cs[-1])
    c_time.append(time)
    c_message.append(message)

d = 1
n = 1
while n < 64:
    c_0 = []
    c_1 = []
    bit0_red = []
    bit0_nored = []
    bit1_red = []
    bit1_nored = []
    for i, c in enumerate(cs):
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
        ci_1, _ = mont_mul(ci_1, c_p[i], N, l_N, omega)
        c_1.append(ci_1)
        ci_1, red1 = mont_mul(ci_1, ci_1, N, l_N, omega)
        if red1:
            bit1_red.append(c_time[i])
        else:
            bit1_nored.append(c_time[i])

    mean_bit0_red = sum(bit0_red) // len(bit0_red)
    mean_bit0_nored = sum(bit0_nored) // len(bit0_nored)
    mean_bit1_red = sum(bit1_red) // len(bit1_red)
    mean_bit1_nored = sum(bit1_nored) // len(bit1_nored)

    if mean_bit0_red - mean_bit0_nored > mean_bit1_red - mean_bit1_nored:
        d <<= 1
        c_cur = c_0
    else:
        d = (d << 1) + 1
        c_cur = c_1
    n += 1
    print d
    test_c = random.randint(2, N)
    _, interaction = interact(test_c)
    if pow(test_c, d << 1, N) == interaction:
        d <<= 1
        break
    elif pow(test_c, (d << 1) + 1, N) == interaction:
        d = (d << 1) + 1
        break

print "%X" % d
print pow(test_c, d, N)
print interaction
