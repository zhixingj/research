import numpy as np
import os

def sn1a(one_a, z, z_max, b, numiso, a):
    iso1a = np.zeros((z_max, numiso))
    f = np.tanh(b)
    g = np.tanh(a-b)
    for i in range(0, z_max):
        for j in range(0, numiso):
            iso1a[i,j] = (one_a[j]) * z[i] * ((np.tanh(a*z[i]-b))+f) / (g+f)
    return iso1a