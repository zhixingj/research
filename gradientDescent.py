import tensorflow as tf
import numpy as np
import os
import cProfile
from sympy.utilities.lambdify import implemented_function
from sympy import *
from matplotlib import pyplot as plt
import random
from scipy.stats import chisquare

# from tensorflow.python.framework.ops import disable_eager_execution

def gradDescent(exp, solarM, index, N):
    Num = N+1
    fraction = symbols('f1:%d'%Num)
    tauValue = symbols('k1:%d'%Num)
    fs = []
    ks = []
    for i in range(N):
        fs.append(tf.Variable(random.random()*0.05))
        ks.append(tf.Variable(random.random()*0.05))

    # f11 = tf.Variable(0.05)
    # f22 = tf.Variable(0.3)
    # k11 = tf.Variable(0.6)
    # k22 = tf.Variable(0.1)

    # f1,k1,f2,k2 = symbols('f1 k1 f2 k2')
    wsFinals = []
    wsSolarM = []
    for iso in index:
        func = exp[iso]
        func = lambdify(list(fraction) + list(tauValue), func)
        wsFinals.append(func)
        wsSolarM.append(solarM[iso])
    optimizer = tf.compat.v1.train.GradientDescentOptimizer(10)
    lossAll = []
    loss = 100
    step = 0
    minLoss = 5e-8
    minF = []
    minK = []
    minStep = 1
    for i in range(2000):
        fOutBound = [x for x in range(len(fs)) if fs[x] <= 0 or fs[x] >= 1]
        kOutBound = [y for y in range(len(ks)) if ks[y] <= 0 or ks[y] >= 1]
        for i in fOutBound:
            fs[i].assign(random.random()*0.1)
        for i in kOutBound:
            ks[i].assign(random.random()*0.1)

        with tf.GradientTape() as gt:
            step +=1
            gt.watch(fs)
            gt.watch(ks)
            loss = 0
            fAndK =fs+ks
            abu = []
            for iso in range(len(index)):
                wsValue = wsFinals[iso](*fAndK)
                abu.append(wsValue)
                loss += abs(wsValue-wsSolarM[iso])
            # print('Loss:', loss)
            if loss < minLoss:
                minAbu = abu
                minLoss = loss
                minF, minK = [], []
                for f in fs:
                    minF.append(f.read_value())
                for k in ks:
                    minK.append(k.read_value())
                minStep = step
            lossAll.append(loss)
            gradients = gt.gradient(loss, fs+ks)
            optimizer.apply_gradients(zip(gradients, fs+ks))
         
    print('--Final MinVariables: \n', 'F: ',minF,'\nK: ', minK)
    steps = np.arange(step)
    print('-Min Loss:', minLoss)
    plt.scatter(steps,np.array(lossAll), marker = '.', s = 1)
    plt.scatter(minStep, minLoss, marker = '*', s = 10, c='r')
    plt.axhline(y=0, color='b', linestyle='--')
    plt.ylim(bottom = 0)
    plt.show(block=True)
    if minF == [] or minK ==[]:
        print('No ideal min loss')
        exit(0)
    else:
        chiSq = chisquare(abu, wsSolarM)
        print('-------Chi Square:', chiSq)
        return minF, minK
