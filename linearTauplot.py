#f 0.4 0.49 t 0.02 0.034

import numpy as np
import os
from numpy import linalg as LA
from scipy.linalg import expm, sinm, cosm
from matplotlib import pyplot as plt
from weak_s_target import weak_s_target
from weak_s_daughter import weak_s_daughter
from weak_s_index import weak_s_index
from w_cs import w_cs
from weak_source import weak_source
from matrix import matrix
from branch import  branch
from sympy import *
from gradientDescent import gradDescent
from multipletau import multipleTau
from calculateWs import calculateYs

#GET SOLAR ABUN
# def linearTau(a):

numiso = 287

dirname = os.path.dirname(__file__)
sollo20File = os.path.join(dirname, 'sollo20.dat')
abu_array = np.loadtxt(sollo20File, usecols = 1, skiprows = 2, dtype = 'float')

abu_array2 = np.empty(287*2, dtype = 'object')

# Make Ion Array
isotopeAsFile = os.path.join(dirname, 'isotope_As.dat')
isomass = np.loadtxt(isotopeAsFile, dtype = 'str')

isoNameFile = os.path.join(dirname, 'isotope_names.dat')
isoname = np.loadtxt(isoNameFile, dtype = 'str')
isoname = np.char.lower(isoname)

ion = isomass.astype(str)
ion = np.char.add(isoname, ion)


sr86 = np.where(ion =='sr86')
sr87 = np.where(ion == 'sr87')

#Get main s process abundances
sfracFile = os.path.join(dirname, 'sfrac_sollo09_gallino.dat')
sfrac = np.loadtxt(sfracFile, skiprows = 1, dtype = 'float')
mains = sfrac*abu_array

# ;get weak s ion targets in an array
sion = weak_s_target(ion)
sion = np.char.lower(sion)

# ;get weak s ion daughters in an array
dion = weak_s_daughter(ion)
dion = np.char.lower(dion)

nion = weak_source(ion)
nion = np.char.lower(nion)

# ;get unique iso names for defining matrix rows and columns
u,indices = np.unique(sion, return_index = True)
array = sion[np.sort(indices)]
uniq = np.append(array,dion[dion.size-1])

# ;get indicies for weak s path relative to full ion array
# ;these are used for:
# ;reconstituting final weak s into ion format
# ;parsing solar abun and isomass to compare weak s with them
# ;note s_index will NOT contain unstable isos
s_index = weak_s_index(uniq, ion)
# print(s_index)


s_index_new = weak_s_index(nion, ion)



# ;TAKE SOLAR TO BE SOLAR - MAIN S-PROCESS, FIT THE WEAK S-PROCESS TO THIS!
# ;TO UNDO, JUST REMOVE -MAINS(XXX) TERMS FROM BELOW

Ge70=np.where(ion == 'ge70')
Se76=np.where(ion == 'se76')
Kr80=np.where(ion == 'kr80')
Kr82=np.where(ion == 'kr82')
Sr86=np.where(ion == 'sr86')
Sr87=np.where(ion == 'sr87')

Ge70=np.where(uniq == 'ge70')
Se76=np.where(uniq == 'se76')
Kr80=np.where(uniq == 'kr80')
Kr82=np.where(uniq == 'kr82')
Sr86=np.where(uniq == 'sr86')
Sr87=np.where(uniq == 'sr87')

# print(ion)
# print(uniq)

# ;get solar/mass arrays (careful: must use abundances after decays)
# ; Se79 -> Br89
# ; Zr93 -> Nb93
# ; Tc99 -> Ru99
solar=abu_array[s_index_new] #-mains(s_index)
mass=isomass[s_index_new]
# name = ion[s_index]
# ;Se79=where(uniq eq 'Se79')
Zr93=np.where(nion == 'zr93')
Tc99=np.where(nion == 'tc99')
# ;Br79=where(ion eq 'Br79')
Nb93=np.where(ion == 'nb93')
Ru99=np.where(ion == 'ru99')
# ;solar[Se79]=abu_array[Br79];-mains(Br79)
solar[Zr93]=abu_array[Nb93] #-mains(Nb93)
solar[Tc99]=abu_array[Ru99] #-mains(Ru99)
# ;mass[Se79]=79d0
mass[Zr93]=93
mass[Tc99]=99

solarM=abu_array[s_index_new]-mains[s_index_new]
mass=isomass[s_index_new]
# ;Se79=where(uniq eq 'Se79')
Zr93=np.where(nion == 'zr93')
Tc99=np.where(nion == 'tc99')
# ;Br79=where(ion eq 'br79')
Nb93=np.where(ion == 'nb93')
Ru99=np.where(ion == 'ru99')
# ;solarM[Se79]=abu_array[Br79]-mains(Br79)
solarM[Zr93]=abu_array[Nb93]-mains[Nb93]
solarM[Tc99]=abu_array[Ru99]-mains[Ru99]
# ;mass[Se79]=79
mass[Zr93]=93
mass[Tc99]=99

# ;15% Kr80 and 3% Kr82 made by p-process
solarM[Kr80] = solarM[Kr80]-0.15 * solar[Kr80]
solarM[Kr82] = solarM[Kr82]-0.03 * solar[Kr82]



submain=1
zero = np.where(solar == 0)

# ;get weak s cross sections: they match sion array

mains_process = mains[s_index_new]

sig=w_cs(ion,nion)
# b=branch(sion)

# f1,k1,f2,k2 = symbols('f1 k1 f2 k2')
#
f1 = 0.37
k1 = 0.03
f2 = 0.29
k2 = 0.03

# min = 1000
# minF = 1000
# minK = 1000
# f: 0.011-0.111
# k: 0.05-0.15

# for i in range(1,500):
#     f = 0.011+ i * 1/5000
#     for q in range(1,500):
#         k = 0.05 + q * 1/10000

Fe56=np.where(nion == 'fe56')
factorOneA = (f1 * solar[Fe56][0])/(k1)
factorOneB = (f2 * solar[Fe56][0])/(k2)

firstBranch = np.where(nion == "cu64")
firstPosition = firstBranch[0][0]

init_A = 1
init_B = 1
final_A = zeros(nion.size, 1)
final_B = zeros(nion.size, 1)


for i in range(0,firstPosition):
    init_A = init_A * (1/(1+(1/(k1*sig[i]))))
    final_A[i] = init_A
for i in range(0,firstPosition):
    init_B = init_B * (1/(1+(1/(k2*sig[i]))))
    final_B[i] = init_B


cu641_A = final_A[firstPosition-1] * (1/(1-0.3311)+(1/(k1*sig[firstPosition])))**(-1)
cu642_A = final_A[firstPosition-1] * (1/(0.3311)+(1/(k1*sig[firstPosition])))**(-1)
cu64_A = cu641_A + cu642_A
final_A[firstPosition] = cu64_A

cu641_B = final_B[firstPosition-1] * (1/(1-0.3311)+(1/(k2*sig[firstPosition])))**(-1)
cu642_B = final_B[firstPosition-1] * (1/(0.3311)+(1/(k2*sig[firstPosition])))**(-1)
cu64_B = cu641_B + cu642_B
final_B[firstPosition] = cu64_B



niPo = np.where(nion == "ni64")
niPosition = niPo[0][0]
ni64_A = cu641_A * 0.4949912 * (1+(1/(k1*sig[niPosition])))**(-1)
ni64_B = cu641_B * 0.4949912 * (1+(1/(k2*sig[niPosition])))**(-1)

znPo = np.where(nion == "zn64")
znPosition = znPo[0][0]
zn64_A = cu642_A * 2.023255 * (1+(1/(k1*sig[znPosition])))**(-1)
zn64_B = cu642_B * 2.023255 * (1+(1/(k2*sig[znPosition])))**(-1)

final_A[niPosition] = ni64_A
final_A[znPosition] = zn64_A
final_B[niPosition] = ni64_B
final_B[znPosition] = zn64_B

# print(firstPosition)
# print(znPosition)

# ni65 = ni64 * (1+(1/(k*sig[znPosition+1])))**(-1)
# zn65 = zn64 * (1+(1/(k*sig[znPosition+2])))**(-1)
#
# final[znPosition+1] = ni65
# final[znPosition+2] = zn65

cu65Po = np.where(nion == "cu65")
cu65Position = cu65Po[0][0]
# cu65 = ni64 * (1+(1/(k*sig[cu65Position])))**(-1) + zn64 * (1+(1/(k*sig[cu65Position])))**(-1)
cu65_A = zn64_A * (1+(1/(k1*sig[cu65Position])))**(-1)
cu65_B = zn64_B * (1+(1/(k2*sig[cu65Position])))**(-1)

final_A[cu65Position] = cu65_A
final_B[cu65Position] = cu65_B



secondBranch  = np.where(nion == "br80")
secondPosition = secondBranch[0][0]
init2_A = cu65_A
init2_B = cu65_B

for i in range(cu65Position+1,secondPosition):
    init2_A = init2_A * (1/(1+(1/(k1*sig[i]))))
    final_A[i] = init2_A
for i in range(cu65Position+1,secondPosition):
    init2_B = init2_B * (1/(1+(1/(k2*sig[i]))))
    final_B[i] = init2_B

br801_A = final_A[secondPosition-1] * (1/(1-0.0325)+(1/(k1*sig[secondPosition])))**(-1)
br802_A = final_A[secondPosition-1] * (1/(0.0325)+(1/(k1*sig[secondPosition])))**(-1)
br80_A = br801_A+br802_A
final_A[secondPosition] = br80_A

br801_B = final_B[secondPosition-1] * (1/(1-0.0325)+(1/(k2*sig[secondPosition])))**(-1)
br802_B = final_B[secondPosition-1] * (1/(0.0325)+(1/(k2*sig[secondPosition])))**(-1)
br80_B = br801_B + br802_B
final_B[secondPosition] = br80_B


sePo = np.where(nion == "se80")
sePosition = sePo[0][0]
se80_A = br801_A * 0.03359173 * (1+(1/(k1*sig[sePosition])))**(-1)
se80_B = br801_B * 0.03359173 * (1+(1/(k2*sig[sePosition])))**(-1)

krPo = np.where(nion == "kr80")
krPosition = krPo[0][0]
kr80_A = br802_A * 29.76923 * (1+(1/(k1*sig[krPosition])))**(-1)
kr80_B = br802_B * 29.76923 * (1+(1/(k2*sig[krPosition])))**(-1)

final_A[sePosition] = se80_A
final_A[krPosition] = kr80_A

final_B[sePosition] = se80_B
final_B[krPosition] = kr80_B

# se81 = se80 * (1+(1/(k*sig[krPosition+1])))**(-1)
# kr81 = kr80 * (1+(1/(k*sig[krPosition+2])))**(-1)
#
# final[krPosition+1] = se81
# final[krPosition+2] = kr81


br81Po = np.where(nion == "br81")
br81Position = br81Po[0][0]
# br81 = se80 * (1+(1/(k*sig[br81Position])))**(-1) + kr80 * (1+(1/(k*sig[br81Position])))**(-1)
br81_A = kr80_A * (1+(1/(k1*sig[br81Position])))**(-1)
br81_B = kr80_B * (1+(1/(k2*sig[br81Position])))**(-1)

final_A[br81Position] = br81_A
final_B[br81Position] = br81_B

ruPo = np.where(nion == "ru102")
ruPosition = ruPo[0][0]


init3_A = br81_A
for i in range(br81Position+1, ruPosition+1):
    init3_A = init3_A * (1/(1+(1/(k1*sig[i]))))
    final_A[i] = init3_A

init3_B = br81_B
for i in range(br81Position+1, ruPosition+1):
    init3_B = init3_B * (1/(1+(1/(k2*sig[i]))))
    final_B[i] = init3_B


final_A = final_A * factorOneA
final_B = final_B * factorOneB
for i in range(sig.shape[0]):
    final_A[i] = final_A[i]/sig[i]
    final_B[i] = final_B[i]/sig[i]


final = final_A + final_B
final_type = np.asarray(final).astype(np.float64)
finalRatio = []
for i in range(solar.shape[0]):
    finalRatio.append(final_type[i]/solar[i])



Ge70=np.where(nion == 'ge70')
Se76=np.where(nion == 'se76')
Kr80=np.where(nion == 'kr80')
Kr82=np.where(nion == 'kr82')
Sr86=np.where(nion == 'sr86')
Sr87=np.where(nion == 'sr87')

N = 10
final = multipleTau(N)

allWsIso = [Ge70[0][0], Se76[0][0], Kr80[0][0], Kr82[0][0], Sr86[0][0], Sr87[0][0]]
minF, minK = gradDescent(final, solarM, allWsIso, N)

Ys = calculateYs(minF, minK, final, solarM, allWsIso, N)

ge70Value = final[Ge70[0][0]]
se76Value = final[Se76[0][0]]
kr80Value = final[Kr80[0][0]]
kr82Value = final[Kr82[0][0]]
sr86Value = final[Sr86[0][0]]
sr87Value = final[Sr87[0][0]]


sOnlyM = np.zeros(6)
sOnlyM[0] = mass[Ge70]
sOnlyM[1] = mass[Se76]
sOnlyM[2] = mass[Kr80]
sOnlyM[3] = mass[Kr82]
sOnlyM[4] = mass[Sr86]
sOnlyM[5] = mass[Sr87]

# a3 = np.log10(float(ge70Value / solar[Ge70[0]]))
# b3 = np.log10(float(se76Value / solar[Se76[0]]))
# c3 = np.log10(float(kr80Value / solar[Kr80[0]]))
# d3 = np.log10(float(kr82Value / solar[Kr82[0]]))
# w3 = np.log10(float(sr86Value / solar[Sr86[0]]))
# q3 = np.log10(float(sr87Value / solar[Sr87[0]]))

aa3 = np.log10((solar[Ge70][0]-mains_process[Ge70][0])/solar[Ge70][0])
bb3 = np.log10((solar[Se76][0]-mains_process[Se76][0])/solar[Se76][0])
cc3 = np.log10((solar[Kr80][0]-mains_process[Kr80][0])/solar[Kr80][0])
dd3 = np.log10((solar[Kr82][0]-mains_process[Kr82][0])/solar[Kr82][0])
ww3 = np.log10((solar[Sr86][0]-mains_process[Sr86][0])/solar[Sr86][0])
qq3 = np.log10((solar[Sr87][0]-mains_process[Sr87][0])/solar[Sr87][0])

aa3 = np.log10((solarM[Ge70][0])/solar[Ge70][0])
bb3 = np.log10((solarM[Se76][0])/solar[Se76][0])
cc3 = np.log10((solarM[Kr80][0])/solar[Kr80][0])
dd3 = np.log10((solarM[Kr82][0])/solar[Kr82][0])
ww3 = np.log10((solarM[Sr86][0])/solar[Sr86][0])
qq3 = np.log10((solarM[Sr87][0])/solar[Sr87][0])



# array = np.array([a3,b3,c3,d3,w3,q3])
array2 = np.array([aa3,bb3,cc3,dd3,ww3,qq3])


mass = mass.astype(int)
ppplot= plt.figure()
axes= ppplot.add_axes([0.1,0.1,0.8,0.8])
plt.xlabel('mass')
plt.ylabel("log10(weak/solar)")
plt.title("weak tau,f")
# adding axes
axes.scatter(mass,np.log10(finalRatio), marker='.')
# axes.scatter(sOnlyM, array, marker = 's')
axes.scatter(sOnlyM, Ys, marker = 's')
axes.scatter(sOnlyM, array2, marker = '*')
axes.set_xlim([56,105])
axes.set_ylim([-5,1])
plt.axhline(y=0, color='r', linestyle='--')
plt.show()
