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



sig=w_cs(ion,nion)
# b=branch(sion)

# f1,k1,f2,k2 = symbols('f1 k1 f2 k2')
#
# f1 = 1/100
# k1 = 0.03
# f2 = 2/100
# k2 = 0.01

# min = 1000
# minF = 1000
# minK = 1000
# f: 0.011-0.111
# k: 0.05-0.15

# for i in range(1,500):
#     f = 0.011+ i * 1/5000
#     for q in range(1,500):
#         k = 0.05 + q * 1/10000



def multipleTau(numberE):

    numEquations = numberE+1
    fraction = symbols('f1:%d'%numEquations)
    tauValue = symbols('k1:%d'%numEquations)
    finalSum = np.empty(numberE, dtype = "object")

    for j in range(numberE):
        Fe56=np.where(nion == 'fe56')
        factorOne = (fraction[j] * solar[Fe56][0])/(tauValue[j])

        firstBranch = np.where(nion == "cu64")
        firstPosition = firstBranch[0][0]

        init = 1
        final = zeros(nion.size, 1)

        for i in range(0,firstPosition):
            init = init * 1/(1+(1/(tauValue[j]*sig[i])))
            final[i] = init

        cu641 = final[firstPosition-1] * (1/(1-0.3311)+(1/(tauValue[j]*sig[firstPosition])))**(-1)
        cu642 = final[firstPosition-1] * (1/(0.3311)+(1/(tauValue[j]*sig[firstPosition])))**(-1)
        cu64 = cu641 + cu642
        final[firstPosition] = cu64

        niPo = np.where(nion == "ni64")
        niPosition = niPo[0][0]
        ni64 = cu641 * 0.4949912 * (1+(1/(tauValue[j]*sig[niPosition])))**(-1)

        znPo = np.where(nion == "zn64")
        znPosition = znPo[0][0]
        zn64 = cu642 * 2.023255 * (1+(1/(tauValue[j]*sig[znPosition])))**(-1)

        final[niPosition] = ni64
        final[znPosition] = zn64


        cu65Po = np.where(nion == "cu65")
        cu65Position = cu65Po[0][0]
        # cu65 = ni64 * (1+(1/(k*sig[cu65Position])))**(-1) + zn64 * (1+(1/(k*sig[cu65Position])))**(-1)
        cu65 = zn64 * (1+(1/(tauValue[j]*sig[cu65Position])))**(-1)

        final[cu65Position] = cu65

        secondBranch  = np.where(nion == "br80")
        secondPosition = secondBranch[0][0]
        init2 = cu65

        for i in range(cu65Position+1,secondPosition):
            init2 = init2 * 1/(1+(1/(tauValue[j]*sig[i])))
            final[i] = init2

        br801 = final[secondPosition-1] * (1/(1-0.0325)+(1/(tauValue[j]*sig[secondPosition])))**(-1)
        br802 = final[secondPosition-1] * (1/(0.0325)+(1/(tauValue[j]*sig[secondPosition])))**(-1)
        br80 = br801+br802
        final[secondPosition] = br80


        sePo = np.where(nion == "se80")
        sePosition = sePo[0][0]
        se80 = br801 * 0.03359173 * (1+(1/(tauValue[j]*sig[sePosition])))**(-1)

        krPo = np.where(nion == "kr80")
        krPosition = krPo[0][0]
        kr80 = br802 * 29.76923 * (1+(1/(tauValue[j]*sig[krPosition])))**(-1)

        final[sePosition] = se80
        final[krPosition] = kr80

        br81Po = np.where(nion == "br81")
        br81Position = br81Po[0][0]
        br81 = kr80 * (1+(1/(tauValue[j]*sig[br81Position])))**(-1)

        final[br81Position] = br81

        ruPo = np.where(nion == "ru102")
        ruPosition = ruPo[0][0]

        init3 = br81
        for i in range(br81Position+1, ruPosition+1):
            init3 = init3 * 1/(1+(1/(tauValue[j]*sig[i])))
            final[i] = init3

        final = final * factorOne
        for i in range(sig.shape[0]):
            final[9] = final[i]/sig[i]

        finalSum[j] = final

    finalTotal = np.sum(finalSum)
    print('Got final total!')
    return finalTotal


Ge70=np.where(nion == 'ge70')
Se76=np.where(nion == 'se76')
Kr80=np.where(nion == 'kr80')
Kr82=np.where(nion == 'kr82')
Sr86=np.where(nion == 'sr86')
Sr87=np.where(nion == 'sr87')
print('solarM',solarM[Ge70[0][0]])
#
# ge70Value = final[Ge70][0]
# se76Value = final[Se76][0]
# kr80Value = final[Kr80][0]
# kr82Value = final[Kr82][0]
# sr86Value = final[Sr86][0]
# sr87Value = final[Sr87][0]
# multipleFinal= multipleTau(3)
# print(multipleFinal[Ge70])
