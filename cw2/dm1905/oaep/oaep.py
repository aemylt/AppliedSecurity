import Crypto.Util.number as number
import sys
import subprocess
import math

def ceildiv(x, y):
    q, r = divmod(x, y)
    if r != 0:
        q += 1
    return q

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

while l != 1:
    f_1 *= 2
    c_prime = (c * pow(f_1, e, N)) % N
    c_primeh = "%X" % c_prime
    if len(c_primeh) % 2:
        c_primeh = "0" + c_primeh

    target_in.write("%s\n" % c_primeh)
    target_in.flush()

    l = int(target_out.readline().strip(), 16)
    if (l > 2):
        print l
        sys.exit(0)

f_2 = int((N + B) / B) * (f_1 / 2)

while l == 1:
    c_prime = (c * pow(f_2, e, N)) % N
    c_primeh = "%X" % c_prime
    if len(c_primeh) % 2:
        c_primeh = "0" + c_primeh

    target_in.write("%s\n" % c_primeh)
    target_in.flush()

    l = int(target_out.readline().strip(), 16)
    if l == 1:
        f_2 = f_2 + f_1 / 2
    if l > 2:
        print l
        sys.exit(0)

m_min = int(ceildiv(N, f_2))
m_max = int((N + B) / f_2)

while m_min != m_max:
    f_tmp = int((2 * B) / (m_max - m_min))
    i = int((f_tmp * m_min) / N) 
    f_3 = int(ceildiv(i * N, m_min))
    c_prime = (c * pow(f_3, e, N)) % N
    c_primeh = "%X" % c_prime
    while len(c_primeh) < k * 2:
        c_primeh = "0" + c_primeh

    target_in.write("%s\n" % c_primeh)
    target_in.flush()

    l = int(target_out.readline().strip(), 16)
    if l == 1:
        m_min = int(ceildiv((i * N + B), f_3))
    else:
        m_max = int((i * N + B) / f_3)
    if l > 2:
        print l
        sys.exit(0)

m = m_min
print m
print "%X" % m
