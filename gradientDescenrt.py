import tensorflow as tf
import numpy as np
import os
import cProfile
from sympy.utilities.lambdify import implemented_function
from sympy import *
from matplotlib import pyplot as plt

# from tensorflow.python.framework.ops import disable_eager_execution

def gradDescent(exp, solarM, ge70, se76, kr80):

    f11 = tf.Variable(0.5)
    f22 = tf.Variable(0.4)
    k11 = tf.Variable(0.6)
    k22 = tf.Variable(0.5)

    f1,k1,f2,k2 = symbols('f1 k1 f2 k2')

    trainables = [f11, f22, k11, k22]

    ge70Exp = exp[ge70]
    se76Exp = exp[se76]
    kr80Exp = exp[kr80]

    ge70Ws = solarM[ge70]
    se76Ws = solarM[se76]
    kr80Ws = solarM[kr80]

    optimizer = tf.compat.v1.train.GradientDescentOptimizer(3)
    lossAll = []
    loss = 100
    step = 0
    while loss>0.000001:
        if f11 <= 0 or f11>=1:
            f11.assign(0.5)
            print('reassigned!!!!!!!')
        if f22 <= 0 or f22>=1:
            f22.assign(0.5)
            print('reassigned!!!!!!!')

        if k11 <= 0 or k11>=1:
            k11.assign(0.5)
            print('reassigned!!!!!!!')

        if k22 <= 0 or k22>=1:
            k22.assign(0.5)
            print('reassigned!!!!!!!')

        if f11+f22>=1:
            f11.assign(0.4)
            f22.assign(0.5)
            print('reassigned!!!!!!!')

        with tf.GradientTape() as gt:
            step +=1
            gt.watch(trainables)
            # predictGe70 = ge70Exp.subs([(f1, f11), (f2, f22), (k1,k11), (k2,k22)])
            # predictSe76 = se76Exp.subs([(f1, f11), (f2, f22), (k1,k11), (k2,k22)])
            lossGe70 = abs(ge70Exp - ge70Ws)
            lossSe76 = abs(se76Exp - se76Ws)
            lossKr80 = abs(kr80Exp - kr80Ws)
            loss = lossGe70 + lossSe76 + lossKr80
            loss = lambdify([f1,f2,k1,k2], loss)(f11,f22,k11,k22)
            lossAll.append(loss)
            gradients = gt.gradient(loss, trainables)
            # print('---G:', gradients)
            optimizer.apply_gradients(zip(gradients, trainables))
            # print('loss:', loss)
    print('Variables: ',gt.watched_variables())
    steps = np.arange(step)
    print('Loss:', loss)
    # print('SSSteps: ', steps)
    plt.scatter(steps,lossAll, marker = '.', s=0.1)
    plt.show()
    return gt.watched_variables()
        # optimizer.apply_gradients(zip(gradient, trainables))
