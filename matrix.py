import numpy as np
import os

def matrix(uniq, sion, dion, sig, b):
    # ;CONSTRUCTS COEFFICIENT MATRIX
    # ;ASSUMES THE USE OF COLUMN VECTORS FOR dN/dt = A*N
    # ;IF DESIRE ROW VECTORS (ARRAYS) TAKE TRANSPOSE

    a = np.zeros((uniq.size, uniq.size))

    for i in range(0, sion.size):
        m = np.where(uniq == sion[i])[0][0]
        # print(uniq[m])
        a[m, m] = -1*sig[i]*b[i] + (-1*sig[i]*(1-b[i]))
        # print(m)
        # print(i)
    for i in range(0, dion.size):
        c = np.where(uniq == dion[i])[0][0]
        # print(c)
        # print(uniq[c])
        r = np.where(uniq == sion[i])[0][0]
        # print(uniq[r])
        # print("++")
        # if c != r:
        #     if i != 0:
        a[r,c] = sig[i] * b[i]
        # a[c,r] = sig[i] * b[i]
        # print(a[c, r])
    dirname = os.path.dirname(__file__)
    matrixFile = os.path.join(dirname, "matrix.csv")
    np.savetxt(matrixFile, a, delimiter = ',' )

    return a
