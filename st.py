import numpy as np
import datetime
from datetime import date
import os
from Solar_decompJuly1 import solar_decomp
'''!PATH = '~/dhome/Desktop/Isotopes/GCE'
#DEFINE VARIABLES AND ARRAYS#'''

numiso = 287
num_el = 83
z_max = 300
logz = np.arange(z_max, dtype = 'float') / ((z_max - 1) / 5.5) - 5
#logz = dindgen(z_max)/((z_max-1)/6.9) - 6.9
#logz=x
z = 10**logz
zero = np.where(abs(z - 1) == min(abs(z - 1)))
z[zero] = 1.0

iso_ws = np.zeros((z_max, numiso))
iso_s = np.zeros((z_max, numiso))
iso_ss = np.zeros((z_max, numiso))
iso_ls = np.zeros((z_max, numiso))
iso_hs = np.zeros((z_max, numiso))
iso_r = np.zeros((z_max, numiso))
iso_nup = np.zeros((z_max, numiso))
iso_gam = np.zeros((z_max, numiso))

sollo09 = np.zeros(numiso)
data = np.zeros((numiso, 13))
iso = np.zeros((z_max, numiso))

element_model = np.zeros((z_max, num_el))
element_model_snii = np.zeros((z_max, num_el))
element_model_sn1a = np.zeros((z_max, num_el))
element_lodders = np.zeros((num_el))
hsproc = np.zeros((z_max, num_el))
lsproc = np.zeros((z_max, num_el))
ssproc = np.zeros((z_max, num_el))
sproc = np.zeros((z_max, num_el))
wsproc = np.zeros((z_max, num_el))
rproc = np.zeros((z_max, num_el))
nproc = np.zeros((z_max, num_el))
gproc = np.zeros((z_max, num_el))
IIproc = np.zeros((z_max, num_el))
Iproc = np.zeros((z_max, num_el))


#DEFINE VARIABLES AND ARRAYS#END

#MAKE ARRAYS#
dirname = os.path.dirname(__file__)
isoNameFile = os.path.join(dirname, 'isotope_names.dat')
iso_names = np.loadtxt(isoNameFile, dtype = 'str')
iso_names = np.char.lower(iso_names)

dirname = os.path.dirname(__file__)
isoAsFile = os.path.join(dirname, 'isotope_As.dat')
iso_as = np.loadtxt(isoAsFile, dtype = 'str')

isoListIni = np.char.add(iso_as, "_")
isoList = np.char.add(isoListIni, iso_names)

dirname = os.path.dirname(__file__)
isoZsFile = os.path.join(dirname, 'isotope_Zs.dat')
iso_zs = np.loadtxt(isoZsFile, dtype = 'str')

dirname = os.path.dirname(__file__)
isoNameFile = os.path.join(dirname, 'isotope_names.dat')
isoname = np.loadtxt(isoNameFile, dtype = 'str')
isoname = np.char.lower(isoname)

ion = iso_as.astype(str)
ion = np.char.add(isoname, ion)

# for i in range(0, numiso):
  
#     ion[i] = np.char.add(isoname[i].lower(), ion[i].split()[0])
#     print(ion[i])

dirname = os.path.dirname(__file__)
elementNameFile = os.path.join(dirname, 'element_names.dat')
element_names = np.loadtxt(elementNameFile, dtype = 'str', skiprows = 1)
element_names = np.char.lower(element_names)

dirname = os.path.dirname(__file__)
sollo20File = os.path.join(dirname, 'sollo20.dat')
data4 = np.loadtxt(sollo20File, dtype = 'float', usecols = 1, skiprows = 2)
sollo20 = data4.copy()

dirname = os.path.dirname(__file__)
amuFile = os.path.join(dirname, 'amu.dat')
amu = np.loadtxt(amuFile, dtype = 'float', skiprows = 5)

sollo20 =sollo20 / amu

fe56 = np.where(ion == 'fe56')
mg24 = np.where(ion == 'mg24')
o16 = np.where(ion == 'o16')

dirname = os.path.dirname(__file__)
starfit_613File = os.path.join(dirname, 'starfit_274823_6_13_2013.dat')
data3 = np.loadtxt(starfit_613File)
massive = np.zeros(numiso)
massive[0:286] = data3[0:286]

massive[0:9] = 0.0
massive[76:285] = 0.0


h = np.where(element_names == 'h')
c = np.where(element_names == 'c')
nit = np.where(element_names == 'n')
o = np.where(element_names == 'o')
fl = np.where(element_names == 'f')
neon = np.where(element_names == 'ne')
na = np.where(element_names == 'na')
mg = np.where(element_names == 'mg')
al = np.where(element_names == 'al')
si = np.where(element_names == 'si')
ca = np.where(element_names == 'ca')
pp = np.where(element_names == 'p')
ss = np.where(element_names == 's')
cl = np.where(element_names == 'cl')
ar = np.where(element_names == 'ar')
pot = np.where(element_names == 'k')
sc = np.where(element_names == 'sc')
ti = np.where(element_names == 'ti')
va = np.where(element_names == 'v')
cr = np.where(element_names == 'cr')
mn = np.where(element_names == 'mn')
fe = np.where(element_names == 'fe')
co = np.where(element_names == 'co')
ni = np.where(element_names == 'ni')
cu = np.where(element_names == 'cu')
zn = np.where(element_names == 'zn')
ga = np.where(element_names == 'ga')
ger = np.where(element_names == 'ge')
arsenic = np.where(element_names == 'as')
se = np.where(element_names == 'se')
br = np.where(element_names == 'br')
kr = np.where(element_names == 'kr')
rb = np.where(element_names == 'rb')
sr = np.where(element_names == 'sr')
yt = np.where(element_names == 'y')
zr = np.where(element_names == 'zr')
nb = np.where(element_names == 'nb')
pb = np.where(element_names == 'pb')
mo = np.where(element_names == 'mo')
pd = np.where(element_names == 'pd')
ag = np.where(element_names == 'ag')
uran = np.where(element_names == 'u')
th = np.where(element_names == 'th')
au = np.where(element_names == 'au')
pt = np.where(element_names == 'pt')
gd = np.where(element_names == 'gd')
hf = np.where(element_names == 'hf')
ho = np.where(element_names == 'ho')
ir = np.where(element_names == 'ir')
lu = np.where(element_names == 'lu')
dy = np.where(element_names == 'dy')
os = np.where(element_names == 'os')

ba = np.where(element_names == 'ba')
eu = np.where(element_names == 'eu')
er = np.where(element_names == 'er')



space = 25
a = np.zeros(space)
a = np.arange(space + 1, dtype = 'float') / ((space+1) / 0.04)+5.0


space1 = 25
b = np.zeros(space1)
b = np.arange(space1 + 1, dtype = 'float') / ((space1 + 1) / 0.04) + 2.7


f_space = 10
fp = np.zeros(f_space)
fp = np.arange(f_space + 1, dtype = 'float') / ((f_space + 1) / 0.002)+0.692 


space2 = 1000
hs = np.zeros(space2)
hs = np.arange(space2+1, dtype = 'float')/((space2+1)/0.07)+1.47

space3 = 1000
rp = np.zeros(space3)
rp = np.arange(space3+1, dtype = 'float')/((space3+1)/0.004)+0.935

space4 = 30
wp = np.zeros(space4)
wp = np.arange(space4+1, dtype = 'float')/((space4+1)/0.02)+2.87

space6 = 30
ls = np.zeros(space6)
ls = np.arange(space6+1, dtype = 'float')/((space6+1)/0.01)+1.20

space7 = 26
agba = np.zeros(space7)
agba = np.arange(space7+1, dtype = 'float')/((space7+1)/100.0)+1.0

space8 = 10
agbb = np.zeros(space8)
agbb=[-.00001,-.00005,-.0001,-.0005,-.001,-.005,-.01,-.05,-.1,-.5]


space9 = 10
cof = np.zeros(space9)
cof = np.zeros(space9+1, dtype = 'float')/((space9+1)/100.0)+1.0

#red_chi = DBLARR(space*space1)
N = 0
a_array = np.zeros(space)
b_array = np.zeros(space1)
hs_array = np.zeros(space2)
rp_array = np.zeros(space3)
#param = np.zeros(3,space*space1)
#if ggg == 1 then param = np.zeros(3,space2*space3)
'''if ggg == 'ba':
    param = np.zeros((space2, 2))
if ggg == 'eu':
    param = np.zeros((space3, 2))
#param = np.zeros(f_space, 2)
if ggg == 'mg' or ggg == 'o':
    param = np.zeros((space*space1*f_space, 4))
if ggg == 'sr':
    param = np.zeros((space4*space6, 3))
if ggg == 'zr':
    param = np.zeros((space4, 2))

#if ggg == 'pb' then param = np.zeros(2,space5)
if ggg == 'pb':
    param = np.zeros((space7*space8*space9, 4))
'''
chimin = np.zeros((space, 3))
kk=0
'''MASTER LOOP FOR FINDING PARAMETER SPACES
#FOR D = 0,space-1 DO BEGIN
#FOR U = 0,space1-1 DO BEGIN
#FOR F = 0,f_space -1 DO BEGIN
#FOR V = 0,space2-1 DO BEGIN
#FOR M = 0,space3-1 DO BEGIN
#FOR G = 0,space4-1 DO BEGIN
#FOR ST = 0,space5-1 DO BEGIN (obsolete)
#FOR SL = 0,space6-1 DO BEGIN
#FOR AGA = 0,space7-1 DO BEGIN
#FOR AGB = 0,space8-1 DO BEGIN
#FOR coo = 0,space9-1 DO BEGIN'''

#BEST FIT VALUES FOUND BY SCANNING PARAMETER SPACE (AFTER 6/11/2013)

SL=0
ls[SL] = 1.227
G = 0
wp[G] = 1.230

#3 PARAM
AGA=0
AGB=0
coo=0
agba[AGA]=200.0
agbb[AGB]=-.23
cof[coo]=-2.e-11

F=0
fp[F]=0.6935

D=0
U=0
a[D]=5.024
b[U]=2.722

V=0
M=0
hs[V]=1.509
#hs[V]=2.0
rp[M]=0.938



chimin = np.zeros((space, 3))
kk=0

SL=0
ls[SL] = 1.227
G = 0
wp[G] = 1.230

#3 PARAM
AGA=0
AGB=0
coo=0
agba[AGA]=200.0
agbb[AGB]=-.23
cof[coo]=-2.e-11

F=0
fp[F]=0.6935

D=0
U=0
a[D]=5.024
b[U]=2.722

V=0
M=0
hs[V]=1.509
#hs[V]=2.0
rp[M]=0.938


data = solar_decomp(sollo20, amu, iso_zs, fp[F], ion, massive)
isotopes = data
isotopes = isotopes.astype(float)



# s = data[:, 2]
# ws = data[:, 3]
# r = data[:, 4]
# nu_p = data[:, 5]
# gamma = data[:, 6]
# one_a = data[:, 7]
# ii = data[:, 8]
# #lightest = data[:, 9]
# bbn=data[:, 9]
# gcr=data[:, 10]
# novae_priGCR_nu=data[:,11]
# h_burn=data[:, 12]

# data[:, 1]=iso_as

# snii_z = data[:, 0]
# snii_a = data[:, 1]

# isotopes = isotopes.astype(float)

# fraction = np.zeros((287, 11))
# for i in range(0,numiso):
#   ticks=0
#   for j in range(0, 11):
#     fraction[i, j] = isotopes[i, 2 + j] / sollo20[i]
#     if fraction[i, j] < 0.0:
#         fraction[i, j] = 0.0
#     if fraction[i, j] == 1.0:
#       ticks=1
#       col=j
#   if ticks == 1:
#     fraction[i,:]=0.0
#     fraction[i, col]=1.0

# #zero out tiny negligible contributions
# for i in range(0, numiso):
#   for j in range(0, 11):
#     if fraction[i, j] < 1.e-6:
#         fraction[i, j]=0.0

# #RE-MAKE FRACTION INTO STRING ARRAYS

# fraction2 = fraction
# fraction = np.empty((287, 11), dtype = 'object')

# # mmm = str(fraction2[282,0]).split()
# # print(mmm)
# for i in range(0, numiso):
#     for j in range(0,11):
#         fraction[i, j] = str(fraction2[i, j]).split()[0]
#         if float(fraction[i, j]) == 0.0000000:
#             fraction[i, j]= '0.0'
# # print(fraction)
# fraction = fraction.astype(float)
# print(fraction)

# dirname = os.path.dirname(__file__)
# resultFile = os.path.join(dirname, "df.csv")
# header = "fff"
# faa = np.array([ 1.2, 2.3, 3.2, 4.4, 6.6, 2.3]
# np.save
# np.savetxt(resultFile, fa, header=header, delimiter=",")

# # dirname = os.path.dirname(__file__)
# # resultFile = os.path.join(dirname, 'result.csv')
# #header = 'iso_Zs, amu, abu_array, sproc, wsproc, rproc, nproc, gproc, w7, massive2, bbn, gcr, novae_priGCR_nu, h_burn'




# # f = open("myaaa.txt", "w")
# # f.write("fff")
# # # for i in range(0,numiso):
# # #     for j in range(0,11):
# # #         f.write(fraction[i,j])
# # #     f.write("\n")

# # f.close()
# # print(f)



# # b = np.array(["a", "b"])
# # dirname = os.path.dirname(__file__)
# # name = 'fraction_' + str(date.today()) + '.dat'
# # fracFile = os.path.join(dirname, "a.txt")
# # # fracContent = ''
# # # for i in range(0, numiso):
# # #     fracContent =  fracContent + '\n' + ion[i] +'\s' + sollo20[i] + '\s',fraction[i, 0] +'\s' + fraction[i, 1] + '\s' + fraction[i, 2] + '\s' + fraction[i, 3] + '\s' + fraction[i, 4] + '\s' + fraction[i，5] + '\s' + fraction[i, 6] + '\s' + fraction[i, 7] + '\s' +fraction[i, 8] + '\s' + fraction[i, 9] + '\s' + fraction[i, 10]
# # # header = 'isotope, &, Solar Abun, &, Main S, &, Weak S, &, R, &, Nu-P, &, Gamma, &, SNIa, &, Massive, &, BBN, &, GCR Spallation, &, Novae/Pri GCR/Nu, &, H-Burn'
# # np.savetxt(fracFile, b)







