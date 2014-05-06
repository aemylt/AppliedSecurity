import scipy.io as sio
import numpy

data = sio.loadmat('WS1.mat')

traces = data['traces']
byte_hamming_weight = data['byte_Hamming_weight'][0]
sub_bytes = data['SubBytes'][0]
inputs = [input_data[0] for input_data in data['inputs']]

n = len(traces[0])

keys = [i for i in range(256)]

after_sbox = [[sub_bytes[input_data ^ key] for key in keys] for input_data in inputs]

power_consumption = numpy.array([byte_hamming_weight[after] for after in after_sbox])

key_trace = [[0 for i in range(n)] for j in range(256)]

chunksize = 50
chunks = n / 50

for i in range(256):
    for j in range(1, chunks + 1):
        cmatrix = numpy.corrcoef(traces[:, (j - 1) * chunksize:j * chunksize].T, power_consumption[:, i])
        key_trace[i][(j - 1) * chunksize:j * chunksize] = cmatrix[chunksize][:chunksize]

f = file('test.txt', 'w')

f.write(str([max(byte) for byte in key_trace]))

f.close()
