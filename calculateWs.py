import numpy as np
import sympy as sp
def calculateYs(minF, minK, final, solarM, wsIndex, N):
    Y = []
    Num = N+1
    fraction = sp.symbols('f1:%d'%Num)
    tauValue = sp.symbols('k1:%d'%Num)
    minFAndK = minF + minK
    for iso in wsIndex:
        finalValue = sp.lambdify(list(fraction)+list(tauValue), final[iso])
        finalValue = finalValue(*minFAndK)
        ratio = finalValue/solarM[iso]
        Y.append(np.log10(float(ratio)))
    return Y