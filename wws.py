import numpy as np
import os
from numpy import linalg as LA
from scipy.linalg import expm, sinm, cosm
from matplotlib import pyplot as plt
from weak_s_target import weak_s_target
from weak_s_daughter import weak_s_daughter
from weak_s_index import weak_s_index
from w_cs import w_cs
from matrix import matrix
from branch import  branch
# ;CALCULATES FULL S PROCESS NUCLEOSYNTHESIS VIA MATRIX METHOD
# ;USES b+/b- RATIA FROM KEPLER DECAY.DAT
# ;USES N CAPTURE CROSS SECTIONS FROM KEPLER
# ;USES n/b- RATIA CALCULATED FROM HALF LIVES

# erase
# !P.BACKGROUND = 'FFFFFF'XL   ;Set background color to white
# !P.COLOR = 'OOOOOO'XL
# erase

# ;DO WEAK S PATH FIRST

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

# ;get unique iso names for defining matrix rows and columns
u,indices = np.unique(sion, return_index = True)
array = sion[np.sort(indices)]
uniq = np.append(array,dion[dion.size-1])
# uniq = [array, dion[dion.size-1]]

# ;get indicies for weak s path relative to full ion array
# ;these are used for:
# ;reconstituting final weak s into ion format
# ;parsing solar abun and isomass to compare weak s with them
# ;note s_index will NOT contain unstable isos
s_index = weak_s_index(uniq, ion)
# print(s_index)

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
solar=abu_array[s_index] #-mains(s_index)
mass=isomass[s_index]
# name = ion[s_index]
# ;Se79=where(uniq eq 'Se79')
Zr93=np.where(uniq == 'zr93')
Tc99=np.where(uniq == 'tc99')
# ;Br79=where(ion eq 'Br79')
Nb93=np.where(ion == 'nb93')
Ru99=np.where(ion == 'ru99')
# ;solar[Se79]=abu_array[Br79];-mains(Br79)
solar[Zr93]=abu_array[Nb93] #-mains(Nb93)
solar[Tc99]=abu_array[Ru99] #-mains(Ru99)
# ;mass[Se79]=79d0
mass[Zr93]=93
mass[Tc99]=99

solarM=abu_array[s_index]-mains[s_index]
mass=isomass[s_index]
# ;Se79=where(uniq eq 'Se79')
Zr93=np.where(uniq == 'zr93')
Tc99=np.where(uniq == 'tc99')
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
sig=w_cs(ion,sion)


# ;get branching ratios
b=branch(sion)

# print(uniq.size) 49
# print(sig.size) 50
# print(b.size) 50
# print(dion.size)50
# print(sion.size)50
# ;get coefficient matrix
a=matrix(uniq,sion,dion,sig,b)

# ;q=strarr(1,n_elements(sion))
# ;q[0,*]=sion[*]
# ;print,q

# ;USING dN/dt = A*N, so N is column matrix

# ;0000000000000000000000000000
# ;BEGIN LOOP FOR FINDING TAU
# ;0000000000000000000000000000

# ;tau=dindgen(1000)/1000+0.01
# ;tau[0]=tau[1]/2d0
# ;i=0

tt=np.arange(100)/100+0.01
ttgood=np.zeros(tt.size)
# ;FOR x=0,n_elements(tt)-1 do begin
# ;tau=[.4,.3,tt[x]]
# ;tau=[.1,.2,.25,.3,.35,.4]
# ;tau=[.1,.15,.18,.22,.25,.3,.33,.35,.37,.4]
tau=np.array([.2,.25,.3,.35,.4])
# taunew = np.array([0.29 0.27 0.32 0.38 0.45])
# ;tau=[.1,.2,.3,.4]
# ;tau=[.2,.3]
# for i1 in range(0,10):
#     tau[0] = tau[0] + 0.01
#     for j1 in range(0,10):
#         tau[1] = tau[1] + 0.01
#         for q1 in range(0,10):
#             tau[2] = tau[2] + 0.01
#             for t1 in range(0,10):
        #         for r1 in range(0,10) :
        #
errormin = 10000
taumin=0
minW = 0
count = 0

Ge70=np.where(uniq == 'ge70')
Se76=np.where(uniq == 'se76')
Kr80=np.where(uniq == 'kr80')
Kr82=np.where(uniq == 'kr82')
Sr86=np.where(uniq == 'sr86')
Sr87=np.where(uniq == 'sr87')

for i1 in np.arange(0.001,2,0.02):
    tau[4] = 0.32*i1
    tau[3] = 0.28*i1
    tau[2] = 0.26*i1
    tau[1] = 0.2*i1
    tau[0] = 0.16*i1

    sr=np.empty((tau.size, 2), dtype = 'object')
    # ;a=transpose(a)
    piece=0

    if piece == 1:
        upp = 4
        a = a[0:upp,0:upp]
        uniq=uniq[0:upp]
        solar=solar[0:upp]
        abu_array=abu_array[0:upp]
        sion = sion[0:upp]
        dion = dion[0:upp]

    aa = np.zeros((6, tau.size))
    bb=np.array([[1],[1],[1],[1],[1],[1]])
    sfits=np.zeros((solar.size, tau.size))

    for i in range(0, tau.size):
    # for i in range(0, 1):
        m=a*tau[i]
        eval, evec = LA.eig(m)
        #; Check the results using the eigenvalue equation:
        # maxErr = 0
        # for j in range(0, uniq.size):
        #     # ; A*z = lambda*z
        #     alhs = m.dot(evec[j,:])
        #     arhs = evaluate[j] * evec[j,:]
        #     maxErr = maxErr > MAX(ABS(alhs - arhs))
        # expA1 = np.dot(evec.T, (np.diag(np.exp(eval))).T).T
        expA = expm(m)
        # dirname = os.path.dirname(__file__)
        matrix_expmFile = os.path.join(dirname, "expm_matrix.csv")
        np.savetxt(matrix_expmFile, expA, delimiter = ',' )
        # print(expA[5,0])

        # ;n=expA[0,*]
        n=(expA[0,:])
        # n=(expA[:,0])

        if submain == 1:
            if zero != -1:
                n[zero] = 0

        if piece == 2:
            # ;COMPARE WITH SOLAR
            # ;PICK OUT S-ONLY ISOTOPES TO MATCH
            if max(n/solar) == n[Ge70]/solar[Ge70]:
                print(tau[i],'Ge70')
                sr[i,0]='ge70'
                sr[i,1]=str(tau[i])
            if max(n/solar) == n[Se76]/solar[Se76]:
                print(tau[i],'Se76')
                sr[i,0]='se76'
                sr[i, 1]=str(tau[i])

            if max(n/solar) == n[Kr80]/solar[Kr80]:
                # ccc = max(n/solar)
                # maaa=np.where(n/solar == ccc)
                # print(uniq[maaa])
                print(tau[i],'Kr80')
                sr[i,0]='kr80'
                sr[i,1]=str(tau[i])

            if max(n/solar) == n[Kr82]/solar[Kr82]:
                print(tau[i],'Kr82')
                sr[i, 0]='kr82'
                sr[i, 1]=str(tau[i])

            if max(n/solar) == n[Sr86]/solar[Sr86]:
                print,tau[i],'Sr86'
                sr[i, 0]='sr86'
                sr[i, 1]=str(tau[i])

            if max(n/solar) == n[Sr87]/solar[Sr87]:
                print,tau[i],'Sr87'
                sr[i, 0]='sr87'
                sr[i, 1]=str(tau[i])
            #line 257

        # ;fac=n[Sr86]/solar[Sr86]

        # ;convert n to mass fractions
        # ;n=atoms2massfrac(n,solar,abu_array)


        if submain == 0:
            fac=max(n/solar)
            maxx=np.where(n/solar == fac)
            for k in range(0,uniq.size):
                n[k]=n[k]/fac
            ratio=n/solar
            # plot,mass,alog10(n/solar),psym=2,xrange=[mass[0],mass[n_elements(sion)-1]]

        if submain == 1:
            nonzero= np.where((n != 0) & (solarM != 0))
            fac=max(n[nonzero]/solarM[nonzero])
            # print(fac)
            maxx = np.where(n/solarM == fac)
            for k in range(0, uniq.size):
                n[k]=n[k]/fac
            ratio=n/solarM #something is zero
            ifn = np.where(ratio == float('inf'))
            ratio[ifn] = 100
            # print(ratio, 'ratio')
            # print(ratio)
            # print(n)
            mass = mass.astype(int)

           
        aa[:,i]=[ratio[Ge70],ratio[Se76],ratio[Kr80],ratio[Kr82],ratio[Sr86],ratio[Sr87]]
        sfits[:, i]=ratio
    #
    # print(aa)
    # print(sfits)

    # ;a2=0.96959363d0 ;.2
    # ;a3=0.0d0 ;.25
    # ;a4=1.3878442d0 ;.3

    ttee=1
    if ttee == 1:
        # ;a1=0.0 ;.1
        # ;a2=1.5 ;.2
        # ;a3=0.2 ;.25
        # ;a4=0.4 ;.3
        # ;a5=0.6 ;.35
        # ;a6=0.7 ;.4
            
        a1=1.5 #.2
        a2=0.2 #.25
        a3=0.4 #.3
        a4=0.6 #.35
        a5=0.7 #.4
        for w1 in np.arange(0.001,2,0.02):
            a1=(w1*1.5)
            a2=(w1*0.2)
            a3=(w1*0.4)
            a4=(w1*0.6)
            a5=(w1*0.7)
            t1=sfits[:,0]*(a1)
            t2=sfits[:,1]*(a2)
            t3=sfits[:,2]*(a3)
            t4=sfits[:,3]*(a4)
            t5=sfits[:,4]*(a5)

            array=t1+t2+t3+t4+t5

            sOnly = np.zeros(6)
            sOnly[0] = array[Ge70]
            sOnly[1] = array[Se76]
            sOnly[2] = array[Kr80]
            sOnly[3] = array[Kr82]
            sOnly[4] = array[Sr86]
            sOnly[5] = array[Sr87]

            overP = np.where(array > 1)
            # print('overproductions:', uniq[overP])

            sum = 0

            for i in range(0,6):
                sum = sum + abs(sOnly[i]-1)
            sOnlyE = sum/6
            error = np.sum(array[overP]-1)
            
            # print('overp', array[overP])
            # print(error, 'error')
            if sOnlyE < 0.25 and errormin > error:
                errormin = error
                taumin = tau.copy()
                minW = [a1,a2,a3,a4,a5]
                # print(errormin)
                # print(taumin)
                sOnlyMin = sOnly.copy()
                sOnlyM = np.zeros(6)
                sOnlyM[0] = mass[Ge70]
                sOnlyM[1] = mass[Se76]
                sOnlyM[2] = mass[Kr80]
                sOnlyM[3] = mass[Kr82]
                sOnlyM[4] = mass[Sr86]
                sOnlyM[5] = mass[Sr87]
                print('overproductions:', uniq[overP])
                print('-------error',error)

            # print(count)

ppplot= plt.figure()
axes= ppplot.add_axes([0.1,0.1,0.8,0.8])
plt.xlabel('mass')
plt.ylabel("log10(weighted)")
plt.title('New Weighted')
# adding axes
axes.scatter(mass,np.log10(array), marker='.')
axes.scatter(sOnlyM, np.log10(sOnlyMin), marker = 's')
axes.set_xlim([56,102])
axes.set_ylim([-5,1])
plt.axhline(y=0, color='r', linestyle='--')
plt.show()


# ppplot2= plt.figure()
# axes= ppplot2.add_axes([0.1,0.1,0.8,0.8])
# plt.title("weighted from paper")
# plt.xlabel("mass")
# plt.ylabel("alog10(t1+t2+t3+t4+t5)")
# plt.axhline(y=0, color='r', linestyle='-')

# # adding axes
# axes.plot(mass,np.log10(array), marker='.')
# # aa[:,i]=[ratio[Ge70],ratio[Se76],ratio[Kr80],ratio[Kr82],ratio[Sr86],ratio[Sr87]]
# # sfits[:, i]=ratio
# specific = np.array([70,76,80,82,86,87])

# axes.plot(specific,np.log10(aa), marker='.')
# axes.set_xlim([56,105])
# axes.set_ylim([-5,0.5])
# plt.show()



'''
togg=0
for p in range(0, cc):
  if cc[p] > 0.0:
    togg+=1
if togg == 3:
    for x in range(0,tt.size):
        ttgood[x]=tt[x]

# ;ENDFOR ;TT LOOP

toobig=np.where(array > 1)
array[toobig]=1

# ;MUST RENAME SION SO THEY MATCH THEIR DECAY PRODUCTS
# ; Se79 -> Br89
# ; Zr93 -> Nb93
# ; Tc99 -> Ru99
# ;uniq[Se79]='Br79'
uniq[Zr93]='nb93'
uniq[Tc99]='ru99'

array=array*solarM/solar

s_index=weak_s_index(uniq,ion)
wsfrac=np.zeros(numiso)
for i in range(0, s_index.size):
  wsfrac[s_index[i]]=array[i]
# ;  wsfrac[s_index[i]]=wsabun[i]

ga69=np.where(ion == 'ga69')
wsabun=wsfrac
p_in=get_p(ion)
rfrac=np.zeros(numiso)
for i in range(ga69[0],numiso):
  rfrac[i]=1-(sfrac[i]+wsfrac[i])
rfrac[p_in]=0

# ;openpsfl,'~/Desktop/gce/sprocess/MAIN_WEAK_COLOR.ps'
# xr=[0,240]
# yr=[-5,1]
# plot,xr,yr,/nodata,TITLE='Red: Weak   Blue: Main',XTITLE='Mass Number',YTITLE='Log(abun/solar)',/YS,/XS
# hline,0
# oplot,mass,alog10(array),color=rgb(1d0,0d0,0d0)
# nozero=where(sfrac gt 0.0)
# oplot,isomass[nozero],alog10(sfrac[nozero]),color=rgb(0,0,1d0)
# nozero=where(rfrac gt 0.0)
# oplot,isomass[nozero],alog10(rfrac[nozero]),color=rgb(0,1d0,0)

# ;closeps
# stop

# plot,isomass,alog10(wsfrac),psym=2,yrange=[-8,1]
# oplot,isomass,alog10(sfrac),psym=6

stot=wsfrac+sfrac
big=np.where(stot > 1.0)
if big != -1:
  print('[WS.PRO] overproductions!')
  print(ion[big])
  overprod=stot[big] - 1
  wsfrac[big]-=overprod


# plot,isomass,alog10(wsfrac+sfrac),psym=2,yrange=[-8,1]


wsfrac=str(wsfrac)
# for i in range(0,numiso):
#   wsfrac[i]=strsplit(wsfrac[i],/extract)

wsFile = os.path.join(dirname, 'wsfrac_gallino.dat')
header = ';WS abundance fractions relative to solar, fit with tau=[.2,.25,.3,.35,.4] & cc=[1.5,0.2,0.4,0.6,0.7], exp(matrix) method employed. 7/5/2011'
np.savetxt(wsFile, wsfrac, header = header, delimiter = ',', fmt='%s')

good=np.where(ttgood != 0)
print('--------------')
print(ttgood[good])
print('--------------')

kk=sr[:,1])
pgood=np.where(kk != 0)
srgood=np.empty(pgood.size,2)
srgood[:,0]=sr[pgood, 0]
srgood[:, 1]=sr[pgood, 1]
# ;print,'-----------------'
# ;print,srgood
# ;print,'-----------------'

# ;END LOOP FOR FINDING TAU

# ;EXPORT TO DATA FILE
# ;MUST RENAME SION SO THEY MATCH THEIR DECAY PRODUCTS
# ; Se79 -> Br89
# ; Zr93 -> Nb93
# ; Tc99 -> Ru99
uniq[Se79]='br79'
uniq[Zr93]='nb93'
uniq[Tc99]='ru99'

# ;FIX s_index TO REFLECT UNSTABLE -> STABLE

# ;----finish export here
resultArray = np.stack(uniq, n, axis = 1)
wsabunFile = os.path.join(dirname, 'weaksabun.dat')
np.savetxt(wsabunFile, resultArray, delimiter=',')'''
