import Crypto.Util.number as number
import sys
import subprocess
import math

def ceildiv(x, y):
    q, r = divmod(x, y)
    if r != 0:
        q += 1
    return q

def interact(f):
    c_prime = (c * pow(f, e, N)) % N
    c_primeh = "%X" % c_prime
    while len(c_primeh) < k * 2:
        c_primeh = "0" + c_primeh

    target_in.write("%s\n" % c_primeh)
    target_in.flush()

    l = int(target_out.readline().strip(), 16)
    if (l > 2):
        print l
        sys.exit(0)
    return l

public = open(sys.argv[2], 'r')
N_hex = public.readline()
e_hex = public.readline()
c_hex = public.readline()

N = int(N_hex, 16)
e = int(e_hex, 16)
c = int(c_hex, 16)

k = len(N_hex) / 2
B = pow(2, 8 * (k - 1))

target = subprocess.Popen(args   = sys.argv[ 1 ],
                          stdout = subprocess.PIPE, 
                          stdin  = subprocess.PIPE)

target_out = target.stdout
target_in  = target.stdin

f_1 = 1

l = 0
counter = 0

while l != 1:
    f_1 *= 2
    l = interact(f_1)
    counter += 1

f_2 = int((N + B) / B) * (f_1 / 2)

while l == 1:
    l = interact(f_2)
    if l == 1:
        f_2 = f_2 + f_1 / 2
    counter += 1

m_min = ceildiv(N, f_2)
m_max = (N + B) / f_2

while m_min != m_max:
    f_tmp = 2 * B / (m_max - m_min)
    i = (f_tmp * m_min) / N
    f_3 = ceildiv(i * N, m_min)
    l = interact(f_3)
    if l == 1:
        m_min = ceildiv(i * N + B, f_3)
    else:
        m_max = (i * N + B) / f_3
    counter += 1

m = m_min
print m
print "%X" % m
print c
print pow(m, e, N)
print counter
