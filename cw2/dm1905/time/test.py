import sys
import subprocess
import random
import math
import os

def interact(N, c):
    target_in.write("%X\n" % c)
    target_in.write("%X\n" % N)
    target_in.write("2\n")
    target_in.flush()

    l = int(target_out.readline().strip())
    m = int(target_out.readline().strip(), 16)
    return l

target = subprocess.Popen(args   = os.path.realpath(sys.argv[ 1 ]),
                          stdout = subprocess.PIPE, 
                          stdin  = subprocess.PIPE)

target_out = target.stdout
target_in  = target.stdin

timings = {}

for _ in range(100000):
    N = random.getrandbits(128)
    c = random.randint(0, N)
    time = interact(N, c)
    if time not in timings:
        timings[time] = True

print timings
