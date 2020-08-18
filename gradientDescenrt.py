import tensorflow as tf
import numpy as np
import os
import cProfile
from sympy.utilities.lambdify import implemented_function
from sympy import *
from matplotlib import pyplot as plt

# from tensorflow.python.framework.ops import disable_eager_execution

def gradDescent(exp, solarM, se76,  kr82, sr86):

    f11 = tf.Variable(0.05)
    f22 = tf.Variable(0.3)
    k11 = tf.Variable(0.6)
    k22 = tf.Variable(0.1)

    f1,k1,f2,k2 = symbols('f1 k1 f2 k2')

    trainables = [f11, f22, k11, k22]

    # ge70Exp = exp[ge70]
    se76Exp = exp[se76]
    # kr80Exp = exp[kr80]
    kr82Exp = exp[kr82]
    sr86Exp = exp[sr86]
    # sr87Exp = exp[sr87]

    # ge70Ws = solarM[ge70]
    se76Ws = solarM[se76]
    # kr80Ws = solarM[kr80]
    kr82Ws = solarM[kr82]
    sr86Ws = solarM[sr86]
    # sr87Ws = solarM[sr87]

    optimizer = tf.compat.v1.train.GradientDescentOptimizer(100)
    lossAll = []
    loss = 100
    step = 0
    minLoss = 1e-7
    minVar = []
    minStep = 1
    for i in range(1000):
        if f11 <= 0 or f11>=1:
            print(f11.value(), 'reassigned!!!!!!! f11')
            f11.assign(0.4)
        if f22 <= 0 or f22>=1:
            print(f22.value(), 'reassigned!!!!!!! f22')
            f22.assign(0.5)

        if k11 <= 0 or k11>=1:
            print(k11.value(), 'reassigned!!!!!!! k11')
            k11.assign(0.2)

        if k22 <= 0 or k22>=1:
            print(k22.value(), 'reassigned!!!!!!! k22')
            k22.assign(0.3)

        if f11+f22>=1:
            print(f11.value(), f22.value(), 'reassigned!!!!!!! f11 & f22')
            f11.assign(0.5)
            f22.assign(0.3)
        with tf.GradientTape() as gt:
            step +=1
            gt.watch(trainables)
            # lossGe70 = abs(ge70Exp - ge70Ws)
            lossSe76 = abs(se76Exp - se76Ws)
            # lossKr80 = abs(kr80Exp - kr80Ws)
            lossKr82 = abs(kr82Exp - kr82Ws)
            lossSr86 = abs(sr86Exp - sr86Ws)
            # lossSr87 = abs(sr87Exp - sr87Ws)
            # loss = lossGe70 + lossSe76 + lossKr80 + lossKr82 + lossSr86 + lossSr87
            loss = lossSe76 + lossKr82 + lossSr86
            loss = lambdify([f1,f2,k1,k2], loss)(f11,f22,k11,k22)
            if loss < minLoss:
                minLoss = loss
                minVar = []
                for v in gt.watched_variables():
                    minVar.append(v.read_value())
                minStep = step

            lossAll.append(loss)
            gradients = gt.gradient(loss, trainables)
            # print('---G:', gradients)
            optimizer.apply_gradients(zip(gradients, trainables))
            # print('loss:', loss)
    print('--Last Var:', trainables)
    print('--Final MinVariables: ', minVar)
    steps = np.arange(step)
    print('-Min Loss:', minLoss)
    plt.scatter(steps,np.array(lossAll), marker = '.', s=0.1)
    plt.scatter(minStep, minLoss, marker = '*', s = 3, c='r')
    plt.axhline(y=0, color='b', linestyle='--')
    plt.show(block=True)
    if minVar == []:
        print('No ideal min loss')
        exit(0)
    else:
        return minVar
        # optimizer.apply_gradients(zip(gradient, trainables))
