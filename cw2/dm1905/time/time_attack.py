import Crypto.Util.number as number
import sys
import subprocess
import random

# Helper function for interacting with the system
# Converts an integer into a hexadecimal string, does any necessary padding
# and returns the result as a string.
def interact(c):
    target_in.write("%X\n" % c)
    target_in.flush()

    l = int(target_out.readline().strip())
    m = int(target_out.readline().strip(), 16)
    return l, m

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

c = random.randint(0, N)
l, m = interact(c)

print c

print l
print m
