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
from sympy import symbols, Eq, solve, nonlinsolve, zeros
from linearTau import getY, getIndex, gradDescent

numiso = 287

#GET SOLAR ABUN

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

Ge70=np.where(nion == 'ge70')
Se76=np.where(nion == 'se76')
Kr80=np.where(nion == 'kr80')
Kr82=np.where(nion == 'kr82')
Sr86=np.where(nion == 'sr86')
Sr87=np.where(nion == 'sr87')

# Ge70=np.where(uniq == 'ge70')
# Se76=np.where(uniq == 'se76')
# Kr80=np.where(uniq == 'kr80')
# Kr82=np.where(uniq == 'kr82')
# Sr86=np.where(uniq == 'sr86')
# Sr87=np.where(uniq == 'sr87')

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
# print(solarM)
# print(solar)
# print(mains[s_index_new])
mains_process = mains[s_index_new]

sig=w_cs(ion,nion)
# b=branch(sion)

# f,k = symbols('f k')
#
# f = 0.1108
# k = 0.0482

# f = 0.016
# k = 0.07
# f = 0.09519999999999999
# k = 0.0502
# min = 1000
# minF = 1000
# minK = 1000

# f = 0.0878
# k = 0.0428

f = 0.096
k = 0.050100000000000006

# f = 0.089
# k = 0.0465


Fe56=np.where(nion == 'fe56')
factorOne = (f * solar[Fe56][0])/(k)
firstBranch = np.where(nion == "cu64")
firstPosition = firstBranch[0][0]

init = 1
final = zeros(nion.size,1)


for i in range(0,firstPosition):
    init = init * 1/(1+(1/(k*sig[i])))
    final[i] = init

cu641 = final[firstPosition-1] * (1/(1-0.3311)+(1/(k*sig[firstPosition])))**(-1)
cu642 = final[firstPosition-1] * (1/(0.3311)+(1/(k*sig[firstPosition])))**(-1)
cu64 = cu641 + cu642
final[firstPosition] = cu64



niPo = np.where(nion == "ni64")
niPosition = niPo[0][0]
ni64 = cu641 * 0.4949912 * (1+(1/(k*sig[niPosition])))**(-1)

znPo = np.where(nion == "zn64")
znPosition = znPo[0][0]
zn64 = cu642 * 2.023255 * (1+(1/(k*sig[znPosition])))**(-1)

final[niPosition] = ni64
final[znPosition] = zn64

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
cu65 = zn64 * (1+(1/(k*sig[cu65Position])))**(-1)

final[cu65Position] = cu65



secondBranch  = np.where(nion == "br80")
secondPosition = secondBranch[0][0]
init2 = cu65

for i in range(cu65Position+1,secondPosition):
    init2 = init2 * 1/(1+(1/(k*sig[i])))
    final[i] = init2

br801 = final[secondPosition-1] * (1/(1-0.0325)+(1/(k*sig[secondPosition])))**(-1)
br802 = final[secondPosition-1] * (1/(0.0325)+(1/(k*sig[secondPosition])))**(-1)
br80 = br801+br802
final[secondPosition] = br80

sePo = np.where(nion == "se80")
sePosition = sePo[0][0]
se80 = br801 * 0.03359173 * (1+(1/(k*sig[sePosition])))**(-1)

krPo = np.where(nion == "kr80")
krPosition = krPo[0][0]
kr80 = br802 * 29.76923 * (1+(1/(k*sig[krPosition])))**(-1)

final[sePosition] = se80
final[krPosition] = kr80

# se81 = se80 * (1+(1/(k*sig[krPosition+1])))**(-1)
# kr81 = kr80 * (1+(1/(k*sig[krPosition+2])))**(-1)
#
# final[krPosition+1] = se81
# final[krPosition+2] = kr81


br81Po = np.where(nion == "br81")
br81Position = br81Po[0][0]
# br81 = se80 * (1+(1/(k*sig[br81Position])))**(-1) + kr80 * (1+(1/(k*sig[br81Position])))**(-1)
br81 = kr80 * (1+(1/(k*sig[br81Position])))**(-1)

final[br81Position] = br81

ruPo = np.where(nion == "ru102")
ruPosition = ruPo[0][0]


init3 = br81
for i in range(br81Position+1, ruPosition+1):
    init3 = init3 * 1/(1+(1/(k*sig[i])))
    final[i] = init3

# abudance * sig (weak-s)
final = final * factorOne
# abudance * sig
final2 = final.copy()
final2_type = np.asarray(final2).astype(np.float64)

# abudance (weak-s)
for i in range(sig.shape[0]):
    final[i] = final[i]/sig[i]

final_type = np.asarray(final).astype(np.float64)

# abundance (s)
s_process = final_type + mains[s_index_new]
finalRatio = final_type.copy()
finalRatio_Stotal = s_process.copy()
for i in range(final_type.shape[0]):

    finalRatio[i] = final_type[i]/solar[i]
    finalRatio_Stotal[i] = s_process[i]/solar[i]
print(finalRatio.size)

# calculate nC
arrayIndex = mass.copy()
arrayIndex = arrayIndex.astype(int)
unstablePosi = np.where(arrayIndex == 1)
arrayIndex[unstablePosi[0][0]] = 64
arrayIndex[unstablePosi[0][1]] = 80
arrayIndex = arrayIndex-56
sumA = arrayIndex * final_type
sumF = np.sum(sumA)
nC = sumF / (f*solar[Fe56][0])



Ge70=np.where(nion == 'ge70')
Se76=np.where(nion == 'se76')
Kr80=np.where(nion == 'kr80')
Kr82=np.where(nion == 'kr82')
Sr86=np.where(nion == 'sr86')
Sr87=np.where(nion == 'sr87')

ge70Value = final[Ge70[0][0]]
se76Value = final[Se76[0][0]]
kr80Value = final[Kr80[0][0]]
kr82Value = final[Kr82[0][0]]
sr86Value = final[Sr86[0][0]]
sr87Value = final[Sr87[0][0]]

# d1 = solarM[Ge70][0]-ge70Value[0]
# d2 = se76Value[0]-solarM[Se76][0]
# d3 = kr80Value[0]-solarM[Kr80][0]
# d4 = kr82Value[0]-solarM[Kr82][0]
# d5 = sr86Value[0]-solarM[Sr86][0]
# d6 = sr87Value[0]-solarM[Sr87][0]
# sum = d1**2 + d2**2 + d3**2 + d4**2 + d5**2 + d6**2



sOnlyM = np.zeros(6)
sOnlyM[0] = mass[Ge70]
sOnlyM[1] = mass[Se76]
sOnlyM[2] = mass[Kr80]
sOnlyM[3] = mass[Kr82]
sOnlyM[4] = mass[Sr86]
sOnlyM[5] = mass[Sr87]


        # print(ge70Value)
        # print(se76Value)
# print(solar[Ge70])
# print(solar[Se76])
# print(solar[Kr80])
# print(solar[Kr82])
# print(solar[Sr86])
# print(solar[Sr87])

#
# print(solarM[Ge70])
# print(solarM[Se76])
# print(solarM[Kr80])
# print(solarM[Kr82])
# print(solarM[Sr86])
# print(solarM[Sr87])
#
# print(ge70Value)
# print(se76Value)
# print(kr80Value)
# print(kr82Value)
# print(sr86Value)
# print(sr87Value)

# print(np.log10(ge70Value[0]/solarM[Ge70][0]))
# print(np.log10(se76Value[0]/solarM[Se76][0]))
# print(np.log10(kr80Value[0]/solarM[Kr80][0]))
# print(np.log10(kr82Value[0]/solarM[Kr82][0]))
# print(np.log10(sr86Value[0]/solarM[Sr86][0]))
# print(np.log10(sr87Value[0]/solarM[Sr87][0]))

a3 = np.log10(float(ge70Value/solar[Ge70[0][0]]))
b3 = np.log10(float(se76Value/solar[Se76[0][0]]))
c3 = np.log10(float(kr80Value/solar[Kr80[0][0]]))
d3 = np.log10(float(kr82Value/solar[Kr82[0][0]]))
w3 = np.log10(float(sr86Value/solar[Sr86[0][0]]))
q3 = np.log10(float(sr87Value/solar[Sr87[0][0]]))

aa3 = np.log10((solar[Ge70][0]-mains_process[Ge70][0])/solar[Ge70][0])
bb3 = np.log10((solar[Se76][0]-mains_process[Se76][0])/solar[Se76][0])
cc3 = np.log10((solar[Kr80][0]-mains_process[Kr80][0])/solar[Kr80][0])
dd3 = np.log10((solar[Kr82][0]-mains_process[Kr82][0])/solar[Kr82][0])
ww3 = np.log10((solar[Sr86][0]-mains_process[Sr86][0])/solar[Sr86][0])
qq3 = np.log10((solar[Sr87][0]-mains_process[Sr87][0])/solar[Sr87][0])



# a3 = np.log10(finalRatio_Stotal[Ge70][0])
# b3 = np.log10(finalRatio_Stotal[Se76][0])
# c3 = np.log10(finalRatio_Stotal[Kr80][0])
# d3 = np.log10(finalRatio_Stotal[Kr82][0])
# w3 = np.log10(finalRatio_Stotal[Sr86][0])
# q3 = np.log10(finalRatio_Stotal[Sr87][0])

array = np.array([a3,b3,c3,d3,w3,q3])
array2 = np.array([aa3,bb3,cc3,dd3,ww3,qq3])


mass = mass.astype(int)
ppplot= plt.figure()
axes= ppplot.add_axes([0.1,0.1,0.8,0.8])
plt.xlabel('mass')
plt.ylabel("log10(weak/solar)")
plt.title('weak tau = 0.0465, f = 0.089, nC = ' + str(nC))
# adding axes
Ys = getY()
print('Got Ys:', Ys)
# print(mass.size, np.log10(finalRatio).size)
axes.scatter(mass,np.log10(finalRatio), marker='.')
axes.scatter(sOnlyM, Ys, marker = 's')
axes.scatter(sOnlyM, array2, marker = '*')
axes.set_xlim([56,105])
axes.set_ylim([-5,1])
plt.axhline(y=0, color='r', linestyle='--')
plt.show(block = True)


# ppplot= plt.figure()
# axes= ppplot.add_axes([0.1,0.1,0.8,0.8])
# plt.xlabel('mass')
# plt.ylabel("weak")
# plt.title('weak')
# # adding axes
# axes.scatter(mass,np.log10(final2_type), marker='.')
# # axes.scatter(sOnlyM, array, marker = 's')
# axes.set_xlim([56,152])
# axes.set_ylim([-10,5])
# plt.axhline(y=0, color='r', linestyle='--')
# plt.show()


# realequ1 = solarM[Ge70]/solarM[Se76]
# real1 = realequ1[0]
# eq1 = Eq(ge70Value[0]/se76Value[0]-real1,0)
#
# realequ2 = solarM[Kr82]/solarM[Sr86]
# real2 = realequ2[0]
# eq2 = Eq(kr82Value[0]/sr86Value[0]-real2,0)
#
# realequ3 = solarM[Sr86]/solarM[Sr87]
# real3 = realequ3[0]
# eq3 = Eq(sr86Value[0]/sr87Value[0]-real3,0)
# # print(final[Sr86])
# # sol = solve((eq1,eq2), (f,k))
# sol2 = nonlinsolve([eq1],[k])
# sol3 = nonlinsolve([eq2],[k])
# print(eq1)
# print(eq2)
#
# sol4 = solve(eq3)
# # print(sol3)


# eq1 = Eq(ge70Value[0],0)
# eq2 = Eq(se76Value[0],0)
# eq3 = kr80Value = final[Kr80]
