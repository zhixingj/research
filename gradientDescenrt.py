import tensorflow as tf
import numpy as np
import os
import cProfile
from sympy.utilities.lambdify import implemented_function
from sympy import *
from matplotlib import pyplot as plt

# from tensorflow.python.framework.ops import disable_eager_execution

def gradDescent(exp, solarM, kr82, sr86, sr87):

    f11 = tf.Variable(0.13866)
    f22 = tf.Variable(0.2654)
    k11 = tf.Variable(0.03937)
    k22 = tf.Variable(0.03665)

    f1,k1,f2,k2 = symbols('f1 k1 f2 k2')

    trainables = [f11, f22, k11, k22]

    # ge70Exp = exp[ge70]
    # se76Exp = exp[se76]
    # kr80Exp = exp[kr80]
    kr82Exp = exp[kr82]
    sr86Exp = exp[sr86]
    sr87Exp = exp[sr87]

    # ge70Ws = solarM[ge70]
    # se76Ws = solarM[se76]
    # kr80Ws = solarM[kr80]
    kr82Ws = solarM[kr82]
    sr86Ws = solarM[sr86]
    sr87Ws = solarM[sr87]

    optimizer = tf.compat.v1.train.GradientDescentOptimizer(40)
    lossAll = []
    loss = 100
    step = 0
    minLoss = 100
    minVar = []
    minStep = 1
    for i in range(700):
        if f11 <= 0 or f11>=1:
            f11.assign(0.15)
            print(f11, 'reassigned!!!!!!! f11')
        if f22 <= 0 or f22>=1:
            f22.assign(0.25)
            print(f22, 'reassigned!!!!!!! f22')

        if k11 <= 0 or k11>=1:
            k11.assign(0.04)
            print(k11, 'reassigned!!!!!!! k11')

        if k22 <= 0 or k22>=1:
            k22.assign(0.04)
            print(k22, 'reassigned!!!!!!! k22')

        if f11+f22>=1:
            f11.assign(0.3)
            f22.assign(0.4)
            print(f11, f22, 'reassigned!!!!!!! f11 & f22')

        with tf.GradientTape() as gt:
            step +=1
            gt.watch(trainables)
            # predictGe70 = ge70Exp.subs([(f1, f11), (f2, f22), (k1,k11), (k2,k22)])
            # predictSe76 = se76Exp.subs([(f1, f11), (f2, f22), (k1,k11), (k2,k22)])
            lossKr82 = abs(kr82Exp - kr82Ws)
            lossSr86 = abs(sr86Exp - sr86Ws)
            lossSr87 = abs(sr87Exp - sr87Ws)
            # loss = lossGe70 + lossSe76 + lossKr80
            loss = lossKr82 + lossSr86 + lossSr87
            loss = lambdify([f1,f2,k1,k2], loss)(f11,f22,k11,k22)
            if loss<minLoss:
                minLoss = loss
                minVar = gt.watched_variables()
                minStep = step
            lossAll.append(loss)
            gradients = gt.gradient(loss, trainables)
            # print('---G:', gradients)
            optimizer.apply_gradients(zip(gradients, trainables))
            # print('loss:', loss)
    print('MMMinVariables: ', minVar)
    steps = np.arange(step)
    print('Loss:', loss)
    # print('SSSteps: ', steps)
    plt.scatter(steps,np.array(lossAll), marker = '.', s=0.1)
    plt.scatter(minStep, minLoss, marker = '*', s = 3, c='g')
    plt.axhline(y=0, color='r', linestyle='--')
    plt.show()
    return minVar
        # optimizer.apply_gradients(zip(gradient, trainables))
