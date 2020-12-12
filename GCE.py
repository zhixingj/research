import numpy as np
from datetime import date
import os
from Solar_decompJuly1 import solar_decomp
from snii import snii
from sn1a import sn1a
from light import light
from gce_data import gce_data
import math
# np.set_printoptions(threshold=np.inf)

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

data = np.zeros((numiso, 13))
iso = np.zeros((z_max, numiso))

element_model = np.zeros((z_max, num_el))
element_model_snii = np.zeros((z_max, num_el))
element_model_sn1a = np.zeros((z_max, num_el))
element_lodders = np.zeros((num_el))
hsproc = np.zeros((z_max, num_el))
ssproc = np.zeros((z_max, num_el))
lsproc = np.zeros((z_max, num_el))
wsproc = np.zeros((z_max, num_el))
rproc = np.zeros((z_max, num_el))
nproc = np.zeros((z_max, num_el))
sproc = np.zeros((z_max, num_el))
gproc = np.zeros((z_max, num_el))
IIproc = np.zeros((z_max, num_el))
Iproc = np.zeros((z_max, num_el))


#DEFINE VARIABLES AND ARRAYS#END

#MAKE ARRAYS#
dirname = os.path.dirname(__file__)
isoNameFile = os.path.join(dirname, 'isotope_names.dat')
iso_names = np.loadtxt(isoNameFile, dtype = 'str')
np.char.lower(iso_names)

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
np.char.lower(isoname)

ion = iso_as.astype(str)
ion = np.char.add(isoname, ion)
ion = np.asarray(ion)

'''OPENW,lun,'~/Desktop/gce/ions.dat',/GET_LUN
#for i=0,n_elements(ion)-1 do begin 
#PRINTF,lun,ion[i]
#endfor
#close,lun
#free_lun,lun
#stop'''

dirname = os.path.dirname(__file__)
elementNameFile = os.path.join(dirname, 'element_names.dat')
element_names = np.loadtxt(elementNameFile, dtype = 'str', skiprows = 1)
np.char.lower(element_names)
'''OPENR,lun, path + 'solar_abun.dat', /GET_LUN
#header2 = STRARR(1)
#READF,lun, header2
#data4 = DBLARR(3,numiso)
#READF,lun,data4

#sollo09 = data4(2,*)'''

dirname = os.path.dirname(__file__)
sollo20File = os.path.join(dirname, 'sollo20.dat')
data4 = np.loadtxt(sollo20File, usecols = 1, skiprows = 2)

#sollo20 is in mass fractions
sollo20 = data4.copy()

#convert to mol frac
dirname = os.path.dirname(__file__)
amuFile = os.path.join(dirname, 'isomasses_numbers.dat')
amu = np.loadtxt(amuFile, skiprows = 5)

sollo20 =sollo20 / amu

fe56 = np.where(ion == 'fe56')
mg24 = np.where(ion == 'mg24')
o16 = np.where(ion == 'o16')

#GET HEGER DATA#

#NEW starfit DATA 6/13/2012 (IS IN MOL FRACTIONS!!!) 
dirname = os.path.dirname(__file__)
starfit_613File = os.path.join(dirname, 'starfit_274823_6_13_2013.dat')
data3 = np.loadtxt(starfit_613File)
massive = np.zeros(numiso)
massive[0:286] = data3[0:286]

#ONLY USE C12 TO GA69 HEGERDATA
massive[0:9] = 0.0
massive[76:285] = 0.0

#END HEGER DATA#   To line296
element_names = [x.lower() for x in element_names]
element_names = np.asarray(element_names)
h = np.where(element_names == 'h')[0][0]
c = np.where(element_names == 'c')[0][0]
nit = np.where(element_names == 'n')[0][0]
o = np.where(element_names == 'o')[0][0]
fl = np.where(element_names == 'f')[0][0]
neon = np.where(element_names == 'ne')[0][0]
na = np.where(element_names == 'na')[0][0]
mg = np.where(element_names == 'mg')[0][0]
al = np.where(element_names == 'al')[0][0]
si = np.where(element_names == 'si')[0][0]
ca = np.where(element_names == 'ca')[0][0]
pp = np.where(element_names == 'p')[0][0]
ss = np.where(element_names == 's')[0][0]
cl = np.where(element_names == 'cl')[0][0]
ar = np.where(element_names == 'ar')[0][0]
pot = np.where(element_names == 'k')[0][0]
sc = np.where(element_names == 'sc')[0][0]
ti = np.where(element_names == 'ti')[0][0]
va = np.where(element_names == 'v')[0][0]
cr = np.where(element_names == 'cr')[0][0]
mn = np.where(element_names == 'mn')[0][0]
fe = np.where(element_names == 'fe')[0][0]
co = np.where(element_names == 'co')[0][0]
ni = np.where(element_names == 'ni')[0][0]
cu = np.where(element_names == 'cu')[0][0]
zn = np.where(element_names == 'zn')[0][0]
ga = np.where(element_names == 'ga')[0][0]
ger = np.where(element_names == 'ge')[0][0]
arsenic = np.where(element_names == 'as')[0][0]
se = np.where(element_names == 'se')[0][0]
br = np.where(element_names == 'br')[0][0]
kr = np.where(element_names == 'kr')[0][0]
rb = np.where(element_names == 'rb')[0][0]
sr = np.where(element_names == 'sr')[0][0]
yt = np.where(element_names == 'y')[0][0]
zr = np.where(element_names == 'zr')[0][0]
nb = np.where(element_names == 'nb')[0][0]
pb = np.where(element_names == 'pb')[0][0]
mo = np.where(element_names == 'mo')[0][0]
pd = np.where(element_names == 'pd')[0][0]
ag = np.where(element_names == 'ag')[0][0]
uran = np.where(element_names == 'u')[0][0]
th = np.where(element_names == 'th')[0][0]
au = np.where(element_names == 'au')[0][0]
pt = np.where(element_names == 'pt')[0][0]
gd = np.where(element_names == 'gd')[0][0]
hf = np.where(element_names == 'hf')[0][0]
ho = np.where(element_names == 'ho')[0][0]
ir = np.where(element_names == 'ir')[0][0]
lu = np.where(element_names == 'lu')[0][0]
dy = np.where(element_names == 'dy')[0][0]
osmium = np.where(element_names == 'os')[0][0]

ba = np.where(element_names == 'ba')[0][0]
eu = np.where(element_names == 'eu')[0][0]
er = np.where(element_names == 'er')[0][0]

#DEFINE PARAMETERS#

space = 25
a = np.zeros(space)
#a(0) = 4.245D0
#a(0)=24.D0
a = np.arange(space + 1, dtype = 'float') / ((space+1) / 0.04)+5.0
#a = dindgen(space+1)/((space+1)/0.002D0)+4.639D0

space1 = 25
b = np.zeros(space1)
#b(0) = 4.814
#b(0)=16.D0
b = np.arange(space1 + 1, dtype = 'float') / ((space1 + 1) / 0.04) + 2.7
#b = dindgen(space1+1)/((space1+1)/0.001D0)+2.399D0

f_space = 10
fp = np.zeros(f_space)
fp = np.arange(f_space + 1, dtype = 'float') / ((f_space + 1) / 0.002)+0.692 
#fp = DINDGEN(f_space+1)/((f_space+1)/0.0002D0)+0.7082D0 
#fp = 0.7083D0

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
#agbb = np.arange(space8+1, dtype = 'float')/((space8+1)/-0.1)-0.001
agbb=[-.00001,-.00005,-.0001,-.0005,-.001,-.005,-.01,-.05,-.1,-.5]
#agbb*=100.0

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

#END BEST FIT VALUES (AFTER 3/8/2012)

'''SL=0
#ls(SL)=1.2276d0
#G=0
#wp(G) = 1.2287D0


#ST=0
#ss(ST) = 0.325
#AGA=0
#AGB=0
#BEST FIT
#agba(AGA)=307.7D0
#agbb(AGB)=1.623D0
#BEST BEHAVIOUR
#agba(AGA)=500.D0
#agbb(AGB)=3.397D0

#3 PARAM
#AGA=0
#AGB=0
#coo=0
#agba(AGA)=200.D0
#agbb(AGB)=-.29d0
#cof(coo)=-2.e-11

#F=0
#fp(F) = 0.709D0

#D=0
#U=0
#a(D) = 4.635D0
#b(U) = 2.397D0

#V=0
#M=0
#hs(V) = 1.487D0
#rp(M) = 0.968D0

#G=0
#wp(G) = 0.8175D0

#sp(V) = 1.52D0
#rp(M) = 1.061D0

#END BEST FIT VALUES (PRIOR TO 3/8/2012)'''

#CALL ALL SUBROUTINES#

#DO SOLAR DECOMPOSITION
ion = np.asarray([x.lower() for x in ion])
#print('sollo 20', sollo20, '\namu', amu,'\niso_zs', iso_zs, '\nfp[F]', fp[F], '\nion', ion, '\nmassive', massive)
data = solar_decomp(sollo20, amu, iso_zs, fp[F], ion, massive)
isotopes = data.astype(float)
# print(isotopes)
#aqw=dblarr(numiso)
#aqw[*]=isotopes[7,*]*amu[*]
#binwrite,array=aqw[*],filename='data/w_hybridsolarfit1'
#stop

#these arrays are for isobox plots
s = data[:, 2]
ws = data[:, 3]
r = data[:, 4]
nu_p = data[:, 5]
gamma = data[:, 6]
one_a = data[:, 7]
ii = data[:, 8]
#lightest = data[:, 9]
bbn=data[:, 9]
gcr=data[:, 10]
novae_priGCR_nu=data[:,11]
h_burn=data[:, 12]

data[:, 1]=iso_as

#FINISH HEGER DATA
snii_z = data[:, 0]
snii_a = data[:, 1]

fraction = np.zeros((287, 11))
for i in range(0,numiso):
  ticks = 0
  for j in range(0, 11):
    fraction[i, j] = isotopes[i, 2 + j] / sollo20[i]
    if fraction[i, j] < 0.0:
        fraction[i, j] = 0.0
    if fraction[i, j] == 1.0:
      ticks=1
      col=j
  if ticks == 1:
    fraction[i,:]=0.0
    fraction[i, col]=1.0
#zero out tiny negligible contributions
for i in range(0, numiso):
  for j in range(0, 11):
    if fraction[i, j] < 1.e-6:
        fraction[i, j]=0.0

#RE-MAKE FRACTION INTO STRING ARRAYS
fraction2 = fraction.copy()
fraction = np.empty((287, 11), dtype = 'object')
for i in range(0, numiso):
    for j in range(0,11):
        fraction[i, j] = str(fraction2[i, j]).split()[0]
        if float(fraction[i, j]) == 0:
            fraction[i, j] = '\\nodata'

# units are mol frac
dirname = os.path.dirname(__file__)
fracFile = os.path.join(dirname, 'fraction_' + str(date.today()) + '.csv')
content = np.empty((numiso, 13), dtype = 'object')
for i in range(0, numiso):
    content[i, 0] = ion[i]
    content[i, 1] = sollo20[i]
    content[i, 2] = fraction[i, 0]
    content[i, 3] = fraction[i, 1]
    content[i, 4] = fraction[i, 1]
    content[i, 5] = fraction[i, 3]
    content[i, 6] = fraction[i, 4]
    content[i, 7] = fraction[i, 5]
    content[i, 8] = fraction[i, 6]
    content[i, 9] = fraction[i, 7]
    content[i, 10] = fraction[i, 8]
    content[i, 11] = fraction[i, 9]
    content[i, 12] = fraction[i, 10]


# content = np.array([ion[i], sollo20[i], fraction[i, 0], fraction[i, 1], fraction[i, 2], fraction[i, 3], fraction[i, 4], fraction[i, 5], fraction[i, 6], fraction[i, 7], fraction[i, 8],fraction[i, 9],fraction[i, 10]])
header = 'isotope, Solar Abun, Main S, Weak S, R, Nu-P, Gamma, SNIa, Massive, BBN, GCR Spallation, Novae/Pri GCR/Nu, H-Burn'
np.savetxt(fracFile, content, delimiter = ',', header = header, fmt='%s')


#line 606

exportmodel = 0

if exportmodel == 1:
    #MATRICIES FOR RECASTING MODEL AS X*F
    massive_scaled = np.zeros(numiso)
    dirname = os.path.dirname(__file__)
    massive3File = os.path.join(dirname, 'massivedatascaledforminus3.dat')
    massive_scaled = np.loadtxt(massive3File)

    slopes = np.zeros(numiso)
    dirname = os.path.dirname(__file__)
    slopeslogFile = os.path.join(dirname, 'slopeslog.dat')
    slopes = np.loadtxt(slopeslogFile)

    ba130 = np.where(ion == 'ba130')
    pb204 = np.where(ion == 'pb204')
    bi209 = np.where(ion == 'bi209')

    light_s = np.zeros(numiso)
    light_s[0 : ba130[0][0]] = data[0 : ba130[0][0], 2]

    heavy_s = np.zeros(numiso)
    heavy_s[ba130[0][0] : pb204] = data[ba130[0][0] : pb204[0][0], 2]

    strong_s = np.zeros(numiso)
    strong_s[pb204[0][0] : 287] = data[pb204[0][0] : 287, 2]

    strong_s_func = np.zeros(numiso)
    strong_s_func[pb204[0][0] : bi209[0][0]] = 1.0

    #ADD D TO H_BURN ARRAY, IT EVOLVES NOW THE SAME (~Z)
    #MUST ALSO ADD TERMS FOR H1 SCALING: (D + He3 + He4) only the z-dependence part
    #RENAME ARRAY, OTHERWISE SOLAR DECOMP ISOBOX WILL BE WRONG
    h_burnx = h_burn
    h_burnx[1] = sollo20[1]-bbn[1]
    h_burnx[0] = -1.0 * ((sollo20[1]-bbn[1])*amu[1] + (sollo20[2]-bbn[2])*amu[2] + (sollo20[3]-bbn[3])*amu[3])/amu[0]

    #hyd gets * 1.d0*z (total metallicity for sum)
    hyd = np.zeros(numiso)
    hyd[0] = -1.0 * 0.0153 / amu[0]

    #hterms gets * 1.d0
    #put in -1.d0 (sum of X+Y+Z)
    #put in Y and D the bbn part (y-intercept)
    #subtract BBN for H1 (which in mol frac) -- it comes from x_matrix[14] #& add back BBN for H1 in mass fractions
    hterms = np.zeros(numiso)
    hterms[0] = (1.0 - (bbn[1] * amu[1] + bbn[2] * amu[2] + bbn[3]*amu[3])) / amu[0] - bbn[0]    #+ bbn[0]*amu[0])

    #hyd=dblarr(numiso)
    #hmult=(0.0153d0+h_burn(2)*3.016029d0+h_burn(3)*4.002603d0)
    #hyd(0)=0.99998060d0*hmult/1.007825D0
    #hyd(1)=(1.d0-0.99998060d0)*hmult/2.014102D0*2.d0

    #massfrac=dblarr(numiso)
    #massfrac(0)=(1.d0-3.016029d0*bbn(2)-4.002603d0*bbn(3))*0.99998060d0/1.007825D0-bbn[0]
    #massfrac(1)=(1.d0-3.016029d0*bbn(2)-4.002603d0*bbn(3))*(1.d0-0.99998060d0)/2.014102D0*2.d0-bbn[1]


    x_matrix=np.zeros((numiso, 17))
    #x_matrix[0,*]=ls
    #x_matrix[1,*]=hs
    #x_matrix[2,*]=ss
    #x_matrix[3,*]=ws
    #x_matrix[4,*]=r
    #x_matrix[5,*]=nu_p
    #x_matrix[6,*]=gamma
    #x_matrix[7,*]=one_a
    #x_matrix[8,*]=slopes
    #x_matrix[9,*]=massive_scaled
    #x_matrix[10,*]=gcr
    #x_matrix[11,*]=novae_priGCR_nu
    #x_matrix[12,*]=bbn
    #x_matrix[13,*]=h_burn
    x_matrix[:, 0]=light_s
    x_matrix[:, 1]=heavy_s
    x_matrix[:, 2]=strong_s
    x_matrix[:, 3]=strong_s_func
    x_matrix[:, 4]=ws
    x_matrix[:, 5]=r
    x_matrix[:, 6]=nu_p
    x_matrix[:, 7]=gamma
    x_matrix[:, 8]=one_a
    x_matrix[:, 9]=slopes
    x_matrix[:, 10]=massive_scaled
    x_matrix[:, 11]=gcr
    x_matrix[:, 12]=novae_priGCR_nu
    x_matrix[:, 13]=bbn
    x_matrix[:, 14]=h_burnx
    x_matrix[:, 15]=hterms
    x_matrix[:, 16]=hyd


    dirname6 = os.path.dirname(__file__)
    fileName = os.path.join(dirname6, 'x_matrix_'+ str(date.today()) + '.dat')
    np.savetxt(fileName, x_matrix)
    #format='(17(E17.10))' 

    #if z >= -2.5 use this
    #f_matrix=DBLARR(17,1)
    #f_matrix[0]=z**1.227d0
    #f_matrix[1]=z**1.509d0
    #f_matrix[2]=1.d0
    #f_matrix[3]=-2.e-11*(1-(tanh(200d0*z-0.29)/tanh(200d0*1.d0-0.29d0)))
    #f_matrix[4]=z**1.230d0
    #f_matrix[5]=z**.938d0
    #f_matrix[6]=z**.938d0
    #f_matrix[7]=z**((1.509d0+.938d0)/2)
    #f_matrix[8]=z*((tanh(5.024d0*z-2.722d0))+tanh(2.722d0))/(tanh(5.024d0-2.722d0)+tanh(2.722d0))
    #f_matrix[9]=(z-0.00317)
    #f_matrix[10]=1.d0
    #f_matrix[11]=z**1.509d0
    #f_matrix[12]=z**.938d0
    #f_matrix[13]=1.d0
    #f_matrix[14]=z**.938d0
    #f_matrix[15]=1.d0
    #f_matrix[16]=1.d0*z

    #if z < -2.5 use this
    #f_matrix=DBLARR(17,1)
    #f_matrix[0]=z**1.227d0
    #f_matrix[1]=z**1.509d0
    #f_matrix[2]=1.d0
    #f_matrix[3]=-2.e-11*(1-(tanh(200d0*z-0.29)/tanh(200d0*1.d0-0.29d0)))
    #f_matrix[4]=z**1.230d0
    #f_matrix[5]=z**.939d0
    #f_matrix[6]=z**.939d0
    #f_matrix[7]=z**((1.509d0+.938d0)/2)
    #f_matrix[8]=z*((tanh(5.024d0*z-2.722d0))+tanh(2.722d0))/(tanh(5.024d0-2.722d0)+tanh(2.722d0))
    #f_matrix[9]=0.d0
    #f_matrix[10]=z/0.00317
    #f_matrix[11]=z**1.509d0
    #f_matrix[12]=z**.938d0
    #f_matrix[13]=1.d0
    #f_matrix[14]=z**.938d0
    #f_matrix[15]=1.d0
    #f_matrix[16]=1.d0*z


#forward_function sn1a #?
one_a = one_a.astype(np.float)
iso1A = sn1a(one_a, z, z_max, b[U], numiso, a[D])

#forward_function snii
ii = ii.astype(np.float)
snii_z = snii_z.astype(np.float)
snii_a = snii_a.astype(np.float)
isoii = snii(ii, z, z_max, snii_z, snii_a, numiso, massive, ion, fp[F])

#m=dblarr(numiso)
#m[*]=massive_scaled[*]*amu[*]
#sollo09[*]*=amu[*]
#print(np.log10(total(m[60:63])/total(sollo09[60:63]))-np.log10(total(m[13:15])/total(sollo09[13:15]))

zero = np.where(abs(z - 1.0) == min(abs(z - 1.0)))
zm3 = np.where(abs(np.log10(z)+2.5) == min(abs(np.log10(z)+2.5)))
zm1 = np.where(abs(np.log10(z)+0.5) == min(abs(np.log10(z)+0.5)))
zm2 = np.where(abs(np.log10(z)+1.5) == min(abs(np.log10(z)+1.5)))
zm4 = np.where(abs(np.log10(z)+3.5) == min(abs(np.log10(z)+3.5)))

#forward_function light #?
#make lightest array for evolving 1<=A<=11
lightest = (bbn.astype(np.float) + gcr.astype(np.float) + novae_priGCR_nu.astype(np.float) + h_burn.astype(np.float))
lightest[0:2]=sollo20[0:2]
iso_lt = light(lightest, z, z_max, numiso, hs[V], rp[M])

ba130 = np.where(ion == 'ba130')
pb204 = np.where(ion == 'pb204')

#data_str=DBLARR(numiso)
#data_str(pb204:286)=data(2,pb204:286)
#data_ss=agbstr(data_str, z, z_max, agbb(AGB),numiso,agba(AGA))


#MAKE S R NU_P GAMMA ARRAYS
for i in range(0, z_max):
    data = data.astype(np.float)
    z = z.astype(np.float)
    ls = ls.astype(np.float)
    #Pb204-208 are i=278-281
    iso_ls[i, 0:ba130[0][0]] = data[0:ba130[0][0], 2] * (z[i])**(ls[SL])
    iso_hs[i, ba130[0][0]:pb204[0][0]] = data[ba130[0][0] : pb204[0][0], 2]*(z[i])**(hs[V])
    #iso_ss(pb204:286,i) = data(2,pb204:286)*(z(i))**(ss(ST)) 
    #iso_s(0:277,i) = data(2,0:277)*(z(i))**(sp(V))
    #iso_s(282:286,i) = data(2,282:286)*(z(i))**(sp(V))
    iso_ws[i, :] = data[:, 3] * (z[i])**(wp[G])
    iso_r[i, :] = data[:, 4]*(z[i])**(rp[M])
    iso_nup[i, :] = data[:, 5]*(z[i])**(rp[M])
    iso_gam[i, :] = data[:, 6]*(z[i])**((hs[V]+rp[M])/2.0+1.0)

#z_str_cutoff = closest(logz, -2.255)

#agba=10.
#agbb=2.
#f=tanh(agbb)
#g=tanh(agba-agbb)
#plot,logz,np.log10((((tanh(agba*z-agbb))+f)/(g+f)))

#AGA=0
#AGB=0
#agba(AGA)=10.D0
#agbb(AGB)=2.397D0

testing = 0
if testing == 1:
    AGA=0
    AGB=0
    coo=0
    agba[AGA]=200.0
    agbb[AGB]=-.23
    cof[coo]=-2.e-11

z_solar = min(abs(logz))
#f=tanh(agbb(AGB))
#g=tanh(agba(AGA)*z(z_solar)+agbb(AGB))
#plot,x_model,np.log10(iso_ss(pb204,*)/element_model(fe,*))

#zl=closest(logz,-4.5d0)
#z_solar=closest(logz,0.d0)
##yl=data(2,pb204)/10000.d0
#f=tanh(agba(AGA)*z(z_solar)+agbb(AGB))
#g=tanh(agba(AGA)*z(zl)+agbb(AGB))
#mult=1.e2

#redo strong, cutoff at logz=-2.5
for i in range(0, z_max):
    #iso_ss(pb204:286,i) = (data(2,pb204:286)/mult-data(2,pb204:286))*((tanh(agba(AGA)*z(i)-agbb(AGB)))-f)/(g-f)+data(2,pb204:286)
    #iso_ss(pb204:286,i) = data(2,pb204:286)*(tanh(agba(AGA)*z(i)+agbb(AGB))-f)/(g-f)
    #iso_ss(pb204:286,i) = cof(coo)*(tanh(agba(AGA)*z(i)+agbb(AGB))-tanh(agba(AGA)*z(z_solar)+agbb(AGB)))+data(2,pb204:286)

    #THESIS UP TO DEFENSE VALUES
    #iso_ss(pb204:286,i) = cof(coo)*(1-(tanh(agba(AGA)*z(i)+agbb(AGB))/tanh(agba(AGA)*z(z_solar)+agbb(AGB))))+data(2,pb204:286)
    agba[AGA]=140
    agbb[AGB]=-0.05
    #200,-0.8,-0.23

    #BETTER FIT VALUES
    #190,-0.1,-0.05

    #BEST FIT 2 PARAMETERS 7 TYPE IA BASED FUNCTION
    agba[AGA]=140.
    agbb[AGB]=-0.05
    iso_ss[i, pb204[0][0]:287] = data[pb204[0][0]:287, 2]*(np.tanh(agba[AGA]*z[i]+agbb[AGB])+np.tanh(agbb[AGB]))/(np.tanh(agba[AGA]+agbb[AGB])+np.tanh(agbb[AGB]))

    #iso_ss(pb204:286,i) = data(2,pb204:286)*(tanh(agba(AGA)*z(i)-agbb(AGB))+tanh(agbb(AGB)))/(tanh(agba(AGA)+agbb(AGB)))
    #iso_ss(pb204:286,i) = data(2,pb204:286)*(tanh(agba(AGA)*z(i)-agbb(AGB))+tanh(-.05))/(tanh(agba(AGA))-tanh(agbb(AGB)))

#using zl and yl (AND 3 PARAM) technique, we now do not have BBN boundary conditions, so reassign negative values
for i in range(0, z_max):
  for j in range(pb204[0][0], 287):
    if iso_ss[i,j] < 0.0:
        iso_ss[i,j] = 0.0

#END SUBROUTINES#

#MAKE ISO/ELEMENT MATRICIES#


iso = np.zeros((z_max, numiso))
iso = isoii + iso1A + iso_ws + iso_r + iso_nup + iso_gam + iso_lt + iso_ss + iso_ls + iso_hs
model_massfrac = np.zeros((z_max, numiso))
lodders_massfrac = np.zeros(numiso)

for i in range (0,numiso):
    #model_massfrac(i,*) = iso(i,*)*DOUBLE(iso_As(i))
    #lodders_massfrac(i) = sollo20(i)*DOUBLE(iso_As(i)) 
    model_massfrac[:, i] = iso[:, i]*amu[i]
    lodders_massfrac[i] = sollo20[i]*amu[i]
# print(min(z+1), min(z+3))
# print(logz)
z_neg1_index = np.where(abs(logz+1) == min(abs(logz+1)))
z_neg3_index = np.where(abs(logz+3)  == min(abs(logz+3)))
# print(z)
neg1_iso = iso[z_neg1_index,:]
neg3_iso = iso[z_neg3_index, :]
#scale linearly between bbn and solar
#s=dblarr(numiso)
#scaled=dblarr(300,numiso)
#sc_em=dblarr(300,83)
#for i=0,numiso-1 do begin
#   s[i]=(sollo09[i]-bbn[i])/(1.d0-0.d0)
#endfor 

'''#for i=0,n_elements(logz)-1 do begin
#for j=0,numiso-1 do begin
#   scaled[i,j]=s[j]*(z[i]-1.d0)+sollo09[j]
#endfor
#endfor'''

i = 0
j = 0
while i < numiso:
    index = np.where(iso_names[i] == iso_names)
    element_model[:,j] = np.sum(iso[:, index[0]], axis = 1)
    model_massfrac[:, j] = np.sum(model_massfrac[:, index[0]], axis = 1)
    lodders_massfrac[j] = np.sum(lodders_massfrac[index[0]])
    lsproc[:, j] = np.sum(iso_ls[:, index[0]], axis = 1)
    hsproc[:, j] = np.sum(iso_hs[:, index[0]], axis = 1)
    ssproc[:, j] = np.sum(iso_ss[:, index[0]], axis = 1)
    sproc[:, j] = np.sum(iso_s[:, index[0]], axis = 1)
    rproc[:, j] = np.sum(iso_r[:, index[0]], axis = 1)
    nproc[:, j] = np.sum(iso_nup[:, index[0]], axis = 1)
    gproc[:, j] = np.sum(iso_gam[:, index[0]], axis = 1)
    wsproc[:, j] = np.sum(iso_ws[:, index[0]], axis = 1)
    IIproc[:, j] = np.sum(isoii[:, index[0]], axis = 1)
    Iproc[:, j] = np.sum(iso1A[:, index[0]], axis = 1)
    element_lodders[j] = np.sum(sollo20[index[0]])
    element_model_snii[:, j] = np.sum(isoii[:, index[0]], axis = 1)
    element_model_sn1a[:, j] = np.sum(iso1A[:, index[0]], axis = 1)

    #sc_em(*,j) = total(scaled(*,index))

    j = j + 1
    #i = i + n_elements(index)
    i += index[0].size
# line 927

#MAKE ELEMENT MODEL FOR SCALED SOLAR
element_lodders = element_lodders.astype(np.float)
bbn = bbn.astype(np.float)
scaled_solar = np.zeros((83, 300))
scaled_feh = np.zeros(300)
for j in range(0,element_model[0, :].size): #elements
    for i in range(0, element_model[:, 0].size): #z
        #if j eq 0 then scaled_solar[i,j]=-element_lodders[j]*z[i]+bbn[0]+bbn[1]
        if j == 0:
            scaled_solar[j,i]=(element_lodders[j]-bbn[0]-bbn[1])*z[i] #+bbn[0]+bbn[1]
        if j == 1:
            scaled_solar[j,i] = (element_lodders[j]-bbn[2]-bbn[3])*z[i] #+bbn[2]+bbn[3]
        if j == 2:
            scaled_solar[j,i]=(element_lodders[j]-bbn[4]-bbn[5])*z[i] #+bbn[4]+bbn[5]
        if j > 2:
            scaled_solar[j,i]=element_lodders[j]*z[i]

for i in range(0, element_model[:, 0].size): #z
    scaled_feh[i]=np.log10(scaled_solar[25, i]/element_lodders[25]) - np.log10(scaled_solar[0, i]/element_lodders[0])

#END ISO/ELEMENT MATRICIES#

#MAKE ARRAYS AND CALL CHI COMPUTE

fe_h_model = np.zeros(z_max)
mg_fe_model = np.zeros(z_max)
mn_fe_model = np.zeros(z_max)
na_fe_model = np.zeros(z_max)
ba_fe_model = np.zeros(z_max)
ni_fe_model = np.zeros(z_max)
eu_fe_model = np.zeros(z_max)

sr_fe_model = np.zeros(z_max)
zr_fe_model = np.zeros(z_max)
metallicity = np.zeros(z_max)
mg_h_model = np.zeros(z_max)
o_h_model = np.zeros(z_max)
ba_h_model = np.zeros(z_max)
sba_h_model = np.zeros(z_max)
rba_h_model = np.zeros(z_max)
gba_h_model = np.zeros(z_max)
eu_h_model = np.zeros(z_max)
mn_h_model = np.zeros(z_max)
ni_h_model = np.zeros(z_max)
solarz = np.sum(lodders_massfrac[2:83])
ba_fe_model_hs = np.zeros(z_max)
ba_fe_model_r = np.zeros(z_max)
ba_fe_model_g = np.zeros(z_max)
sr_fe_model_ls = np.zeros(z_max)
sr_fe_model_r = np.zeros(z_max)
sr_fe_model_ws = np.zeros(z_max)
sr_ba_model = np.zeros(z_max)
eu_fe_model_s = np.zeros(z_max)
eu_fe_model_r = np.zeros(z_max)
eu_fe_model_g = np.zeros(z_max)
mg_fe_model_ii = np.zeros(z_max)
mg_fe_model_1a = np.zeros(z_max)
o_fe_model_ii = np.zeros(z_max)
o_fe_model_1a = np.zeros(z_max)
mn_fe_model_ii = np.zeros(z_max)
mn_fe_model_1a = np.zeros(z_max)

c_fe_model = np.zeros(z_max)
n_fe_model = np.zeros(z_max)
o_fe_model = np.zeros(z_max)
f_fe_model = np.zeros(z_max)
ne_fe_model = np.zeros(z_max)
na_fe_model = np.zeros(z_max)
mg_fe_model = np.zeros(z_max)
al_fe_model = np.zeros(z_max)
si_fe_model = np.zeros(z_max)
p_fe_model = np.zeros(z_max)
s_fe_model = np.zeros(z_max)
cl_fe_model = np.zeros(z_max)
ar_fe_model = np.zeros(z_max)
k_fe_model = np.zeros(z_max)
ca_fe_model = np.zeros(z_max)
ir_fe_model = np.zeros(z_max)
er_fe_model = np.zeros(z_max)
sc_fe_model = np.zeros(z_max)
ti_fe_model = np.zeros(z_max)
v_fe_model = np.zeros(z_max)
cr_fe_model = np.zeros(z_max)
mn_fe_model = np.zeros(z_max)
co_fe_model = np.zeros(z_max)
os_fe_model = np.zeros(z_max)
ni_fe_model = np.zeros(z_max)
cu_fe_model = np.zeros(z_max)
zn_fe_model = np.zeros(z_max)
ga_fe_model = np.zeros(z_max)
ge_fe_model = np.zeros(z_max)
as_fe_model = np.zeros(z_max)
se_fe_model = np.zeros(z_max)
br_fe_model = np.zeros(z_max)
kr_fe_model = np.zeros(z_max)
rb_fe_model = np.zeros(z_max)
sr_fe_model = np.zeros(z_max)
y_fe_model = np.zeros(z_max)
zr_fe_model = np.zeros(z_max)
pb_fe_model = np.zeros(z_max)
pb_fe_model_ss = np.zeros(z_max)
pb_fe_model_r = np.zeros(z_max)
ag_fe_model = np.zeros(z_max)
al_fe_model = np.zeros(z_max)
au_fe_model = np.zeros(z_max)
cr_fe_model = np.zeros(z_max)
dy_fe_model = np.zeros(z_max)
gd_fe_model = np.zeros(z_max)
hf_fe_model = np.zeros(z_max)
ho_fe_model = np.zeros(z_max)
lu_fe_model = np.zeros(z_max)
mo_fe_model = np.zeros(z_max)
pd_fe_model = np.zeros(z_max)
pt_fe_model = np.zeros(z_max)
u_fe_model = np.zeros(z_max)
cr_fe_model_ii = np.zeros(z_max)
cr_fe_model_ia = np.zeros(z_max)
mn_fe_model_ii = np.zeros(z_max)
mn_fe_model_ia = np.zeros(z_max)

fe_model=np.zeros(z_max)
sr_model=np.zeros(z_max)
pb_model=np.zeros(z_max)
mg_model=np.zeros(z_max)
o_model=np.zeros(z_max)
eu_model=np.zeros(z_max)

fe_modelII=np.zeros(z_max)
sr_modelII=np.zeros(z_max)
pb_modelII=np.zeros(z_max)
mg_modelII=np.zeros(z_max)
o_modelII=np.zeros(z_max)
eu_modelII=np.zeros(z_max)
eu_modelR=np.zeros(z_max)

sr_modelW=np.zeros(z_max)

ba_mg_model=np.zeros(z_max)
eu_mg_model=np.zeros(z_max)

feP=np.zeros(z_max)

for i in range(0, z_max):
    c_fe_model[i] = np.log10((element_model[i, c] / element_lodders[c])/(element_model[i, fe]/element_lodders[fe]))
    n_fe_model[i] = np.log10((element_model[i, nit] / element_lodders[nit])/(element_model[i, fe]/element_lodders[fe]))
    o_fe_model[i] = np.log10((element_model[i, o]/element_lodders[o])/(element_model[i, fe]/element_lodders[fe]))
    f_fe_model[i] = np.log10((element_model[i, fl]/element_lodders[fl])/(element_model[i, fe]/element_lodders[fe]))
    ne_fe_model[i] = np.log10((element_model[i, neon]/element_lodders[neon])/(element_model[i, fe]/element_lodders[fe]))
    na_fe_model[i] = np.log10((element_model[i, na]/element_lodders[na])/(element_model[i, fe]/element_lodders[fe]))
    mg_fe_model[i] = np.log10((element_model[i, mg]/element_lodders[mg])/(element_model[i, fe]/element_lodders[fe]))
    al_fe_model[i] = np.log10((element_model[i, al]/element_lodders[al])/(element_model[i, fe]/element_lodders[fe]))
    si_fe_model[i] = np.log10((element_model[i, si]/element_lodders[si])/(element_model[i, fe]/element_lodders[fe]))
    p_fe_model[i] = np.log10((element_model[i, pp]/element_lodders[pp])/(element_model[i, fe]/element_lodders[fe]))
    s_fe_model[i] = np.log10((element_model[i, ss]/element_lodders[ss])/(element_model[i, fe]/element_lodders[fe]))
    cl_fe_model[i] = np.log10((element_model[i, cl]/element_lodders[cl])/(element_model[i, fe]/element_lodders[fe]))
    ar_fe_model[i] = np.log10((element_model[i, ar]/element_lodders[ar])/(element_model[i, fe]/element_lodders[fe]))
    k_fe_model[i] = np.log10((element_model[i, pot]/element_lodders[pot])/(element_model[i, fe]/element_lodders[fe]))
    ca_fe_model[i] = np.log10((element_model[i, ca]/element_lodders[ca])/(element_model[i, fe]/element_lodders[fe]))
    ir_fe_model[i] = np.log10((element_model[i, ir]/element_lodders[ir])/(element_model[i, fe]/element_lodders[fe]))

    sc_fe_model[i] = np.log10((element_model[i, sc]/element_lodders[sc])/(element_model[i, fe]/element_lodders[fe]))
    ti_fe_model[i] = np.log10((element_model[i, ti]/element_lodders[ti])/(element_model[i, fe]/element_lodders[fe]))
    v_fe_model[i] = np.log10((element_model[i, va]/element_lodders[va])/(element_model[i, fe]/element_lodders[fe]))
    cr_fe_model[i] = np.log10((element_model[i, cr]/element_lodders[cr])/(element_model[i, fe]/element_lodders[fe]))
    mn_fe_model[i] = np.log10((element_model[i, mn]/element_lodders[mn])/(element_model[i, fe]/element_lodders[fe]))
    co_fe_model[i] = np.log10((element_model[i, co]/element_lodders[co])/(element_model[i, fe]/element_lodders[fe]))
    ni_fe_model[i] = np.log10((element_model[i, ni]/element_lodders[ni])/(element_model[i, fe]/element_lodders[fe]))
    cu_fe_model[i] = np.log10((element_model[i, cu]/element_lodders[cu])/(element_model[i, fe]/element_lodders[fe]))
    zn_fe_model[i] = np.log10((element_model[i, zn]/element_lodders[zn])/(element_model[i, fe]/element_lodders[fe]))
    ga_fe_model[i] = np.log10((element_model[i, ga]/element_lodders[ga])/(element_model[i, fe]/element_lodders[fe]))
    ge_fe_model[i] = np.log10((element_model[i, ger]/element_lodders[ger])/(element_model[i, fe]/element_lodders[fe]))
    as_fe_model[i] = np.log10((element_model[i, arsenic]/element_lodders[arsenic])/(element_model[i, fe]/element_lodders[fe]))
    se_fe_model[i] = np.log10((element_model[i, se]/element_lodders[se])/(element_model[i, fe]/element_lodders[fe]))
    br_fe_model[i] = np.log10((element_model[i, br]/element_lodders[br])/(element_model[i, fe]/element_lodders[fe]))
    kr_fe_model[i] = np.log10((element_model[i, kr]/element_lodders[kr])/(element_model[i, fe]/element_lodders[fe]))
    rb_fe_model[i] = np.log10((element_model[i, rb]/element_lodders[rb])/(element_model[i, fe]/element_lodders[fe]))
    sr_fe_model[i] = np.log10((element_model[i, sr]/element_lodders[sr])/(element_model[i, fe]/element_lodders[fe]))
    y_fe_model[i] = np.log10((element_model[i, yt]/element_lodders[yt])/(element_model[i, fe]/element_lodders[fe]))
    zr_fe_model[i] = np.log10((element_model[i, zr]/element_lodders[zr])/(element_model[i, fe]/element_lodders[fe]))
    pb_fe_model[i] = np.log10((element_model[i, pb]/element_lodders[pb])/(element_model[i, fe]/element_lodders[fe]))
    os_fe_model[i] = np.log10((element_model[i, osmium]/element_lodders[osmium])/(element_model[i, fe]/element_lodders[fe]))

    ag_fe_model[i] = np.log10((element_model[i, ag]/element_lodders[ag])/(element_model[i, fe]/element_lodders[fe]))
    al_fe_model[i] = np.log10((element_model[i, al]/element_lodders[al])/(element_model[i, fe]/element_lodders[fe]))
    au_fe_model[i] = np.log10((element_model[i, au]/element_lodders[au])/(element_model[i, fe]/element_lodders[fe]))
    cr_fe_model[i] = np.log10((element_model[i, cr]/element_lodders[cr])/(element_model[i, fe]/element_lodders[fe]))
    dy_fe_model[i] = np.log10((element_model[i, dy]/element_lodders[dy])/(element_model[i, fe]/element_lodders[fe]))
    gd_fe_model[i] = np.log10((element_model[i, gd]/element_lodders[gd])/(element_model[i, fe]/element_lodders[fe]))
    hf_fe_model[i] = np.log10((element_model[i, hf]/element_lodders[hf])/(element_model[i, fe]/element_lodders[fe]))
    ho_fe_model[i] = np.log10((element_model[i, ho]/element_lodders[ho])/(element_model[i, fe]/element_lodders[fe]))
    lu_fe_model[i] = np.log10((element_model[i, lu]/element_lodders[lu])/(element_model[i, fe]/element_lodders[fe]))
    mo_fe_model[i] = np.log10((element_model[i, mo]/element_lodders[mo])/(element_model[i, fe]/element_lodders[fe]))
    pd_fe_model[i] = np.log10((element_model[i, pd]/element_lodders[pd])/(element_model[i, fe]/element_lodders[fe]))
    pt_fe_model[i] = np.log10((element_model[i, pt]/element_lodders[pt])/(element_model[i, fe]/element_lodders[fe]))
    u_fe_model[i] = np.log10((element_model[i, uran]/element_lodders[uran])/(element_model[i, fe]/element_lodders[fe]))

    fe_h_model[i] = np.log10((element_model[i, fe]/element_lodders[fe])/(element_model[i, h]/element_lodders[h]))
    mn_fe_model[i] = np.log10((element_model[i, mn]/element_lodders[mn])/(element_model[i, fe]/element_lodders[fe]))
    na_fe_model[i] = np.log10((element_model[i, na]/element_lodders[na])/(element_model[i, fe]/element_lodders[fe]))
    eu_fe_model[i] = np.log10((element_model[i, eu]/element_lodders[eu])/(element_model[i, fe]/element_lodders[fe]))
    er_fe_model[i] = np.log10((element_model[i, er]/element_lodders[er])/(element_model[i, fe]/element_lodders[fe]))
    ba_fe_model[i] = np.log10((element_model[i, ba]/element_lodders[ba])/(element_model[i, fe]/element_lodders[fe]))
    sr_fe_model[i] = np.log10((element_model[i, sr]/element_lodders[sr])/(element_model[i, fe]/element_lodders[fe]))
    zr_fe_model[i] = np.log10((element_model[i, zr]/element_lodders[zr])/(element_model[i, fe]/element_lodders[fe]))
    metallicity[i] = np.log10(np.sum(model_massfrac[i, 2:83] / solarz))

    ni_fe_model[i] = np.log10((element_model[i, ni] / element_lodders[ni])/(element_model[i, fe]/element_lodders[fe]))

    pb_fe_model_ss[i] = np.log10((ssproc[i, pb] / element_lodders[pb])/(element_model[i, fe]/element_lodders[fe]))
    pb_fe_model_r[i] = np.log10((rproc[i, pb] / element_lodders[pb])/(element_model[i, fe]/element_lodders[fe]))

    ba_fe_model_hs[i] = np.log10((hsproc[i, ba]/element_lodders[ba])/(element_model[i, fe]/element_lodders[fe]))
    ba_fe_model_r[i] = np.log10((rproc[i, ba]/element_lodders[ba])/(element_model[i, fe]/element_lodders[fe]))
    ba_fe_model_g[i] = np.log10((gproc[i, ba]/element_lodders[ba])/(element_model[i, fe]/element_lodders[fe]))

    sr_fe_model_ls[i] = np.log10((lsproc[i, sr]/element_lodders[sr])/(element_model[i, fe]/element_lodders[fe]))
    sr_fe_model_r[i] = np.log10((rproc[i, sr]/element_lodders[sr])/(element_model[i, fe]/element_lodders[fe]))
    sr_fe_model_ws[i] = np.log10((wsproc[i, sr]/element_lodders[sr])/(element_model[i, fe]/element_lodders[fe]))

    sr_ba_model[i] = np.log10((element_model[i, sr]/element_lodders[sr])/(element_model[i, ba]/element_lodders[ba]))

    eu_fe_model_s[i] = np.log10((sproc[i,eu]/element_lodders[eu])/(element_model[i, fe]/element_lodders[fe]))
    eu_fe_model_r[i] = np.log10((rproc[i,eu]/element_lodders[eu])/(element_model[i, fe]/element_lodders[fe]))
    eu_fe_model_g[i] = np.log10((gproc[i,eu]/element_lodders[eu])/(element_model[i, fe]/element_lodders[fe]))

    sba_h_model[i] = np.log10((sproc[i, ba]/element_lodders[ba])/(element_model[i, h]/element_lodders[h]))
    rba_h_model[i] = np.log10((rproc[i, ba]/element_lodders[ba])/(element_model[i, h]/element_lodders[h]))
    gba_h_model[i] = np.log10((gproc[i, ba]/element_lodders[ba])/(element_model[i, h]/element_lodders[h]))

    eu_fe_model_s[i] = np.log10((sproc[i,eu]/element_lodders[eu])/(element_model[i, fe]/element_lodders[fe]))
    eu_fe_model_r[i] = np.log10((rproc[i,eu]/element_lodders[eu])/(element_model[i, fe]/element_lodders[fe]))
    eu_fe_model_g[i] = np.log10((gproc[i,eu]/element_lodders[eu])/(element_model[i, fe]/element_lodders[fe]))

    mg_fe_model_ii[i] = np.log10((IIproc[i,mg]/element_lodders[mg])/(element_model[i, fe]/element_lodders[fe]))
    mg_fe_model_1a[i] = np.log10((Iproc[i,mg]/element_lodders[mg])/(element_model[i, fe]/element_lodders[fe]))
    mn_fe_model_ii[i] = np.log10((IIproc[i,mn]/element_lodders[mn])/(element_model[i, fe]/element_lodders[fe]))
    mn_fe_model_1a[i] = np.log10((Iproc[i,mn]/element_lodders[mn])/(element_model[i, fe]/element_lodders[fe]))
    #mg_fe_model_ii[i] = np.log10((IIproc[i,mg]/element_lodders[mg])/(IIproc[i, fe]/element_lodders[fe]))
    #mg_fe_model_1a[i] = np.log10((Iproc[i,mg]/element_lodders[mg])/(Iproc[i, fe]/element_lodders[fe]))

    cr_fe_model_ii[i] = np.log10((IIproc[i,cr]/element_lodders[cr])/(element_model[i, fe]/element_lodders[fe]))
    cr_fe_model_ia[i] = np.log10((Iproc[i,cr]/element_lodders[cr])/(element_model[i, fe]/element_lodders[fe]))

    o_fe_model_ii[i] = np.log10((IIproc[i,o]/element_lodders[o])/(element_model[i, fe]/element_lodders[fe]))
    o_fe_model_1a[i] = np.log10((Iproc[i,o]/element_lodders[o])/(element_model[i, fe]/element_lodders[fe]))

    mn_fe_model_ii[i] = np.log10((IIproc[i,mn]/element_lodders[mn])/(element_model[i, fe]/element_lodders[fe]))
    mn_fe_model_1a[i] = np.log10((Iproc[i,mn]/element_lodders[mn])/(element_model[i, fe]/element_lodders[fe]))

    mg_h_model[i] = np.log10((element_model[i,mg]/element_lodders[mg])/(element_model[i, h]/element_lodders[h]))
    mn_h_model[i] = np.log10((element_model[i,mn]/element_lodders[mn])/(element_model[i, h]/element_lodders[h]))
    o_h_model[i] = np.log10((element_model[i,o]/element_lodders[o])/(element_model[i, h]/element_lodders[h]))
    ni_h_model[i] = np.log10((element_model[i, ni]/element_lodders[ni])/(element_model[i, h]/element_lodders[h]))
    ba_h_model[i] = np.log10((element_model[i, ba]/element_lodders[ba])/(element_model[i, h]/element_lodders[h]))
    eu_h_model[i] = np.log10((element_model[i,eu]/element_lodders[eu])/(element_model[i, h]/element_lodders[h]))

    fe_model[i]=np.log10(element_model[i, fe]/element_lodders[fe])
    sr_model[i]=np.log10(element_model[i, sr]/element_lodders[sr])
    pb_model[i]=np.log10(element_model[i,pb]/element_lodders[pb])
    mg_model[i]=np.log10(element_model[i,mg]/element_lodders[mg])
    o_model[i]=np.log10(element_model[i,o]/element_lodders[o])
    eu_model[i]=np.log10(element_model[i,eu]/element_lodders[eu])

    fe_modelII[i]=np.log10(IIproc[i, fe]/element_lodders[fe])
    sr_modelII[i]=np.log10(IIproc[i, sr]/element_lodders[sr])
    pb_modelII[i]=np.log10(IIproc[i, pb]/element_lodders[pb])
    mg_modelII[i]=np.log10(IIproc[i,mg]/element_lodders[mg])
    o_modelII[i]=np.log10(IIproc[i,o]/element_lodders[o])
    eu_modelR[i]=np.log10(rproc[i,eu]/element_lodders[eu])

    sr_modelW[i]=np.log10(wsproc[i, sr]/element_lodders[sr])

    ba_mg_model[i]=np.log10((element_model[i, ba]/element_lodders[ba])/(element_model[i,mg]/element_lodders[mg]))
    eu_mg_model[i]=np.log10((element_model[i,eu]/element_lodders[eu])/(element_model[i,mg]/element_lodders[mg]))
    
    feP[i]= (element_model[i, fe] / element_lodders[fe])


'''#COMPUTE FIT PARAMETERS IN TERMS OF Z
varswitch = 0
if varswitch == 1:
    #z := xi , z(zero)=1 - solar
    xi=z
    z=10**(metallicity)

    near = min(abs(metallicity - 0.5), p05)  #p05 is the index where feh[*] is closest to +0.5
    near = min(abs(metallicity + 0.0), z0)
    near = min(abs(metallicity + 0.5), m05)
    near = min(abs(metallicity + 1.0), m1)
    near = min(abs(metallicity + 2.0), m2)
    near = min(abs(metallicity + 3.0), m3)
    near = min(abs(metallicity + 4.0), m4)

    near = min(abs(fe_model - 0.5), fp05)  #p05 is the index where feh[*] is closest to +0.5
    near = min(abs(fe_model + 0.0), fz0)
    near = min(abs(fe_model + 0.5), fm05)
    near = min(abs(fe_model + 1.0), fm1)
    near = min(abs(fe_model + 2.0), fm2)
    near = min(abs(fe_model + 3.0), fm3)
    near = min(abs(fe_model + 4.0), fm4)

    #x=dln(s)/dln(z)=(z/s)*(ds/dz)

    #l=(z)**ls(SL)
    #s=(z)**hs(V)
    #r=(z)**rp(M)
    #w=(z)**wp(G)

    #dl=np.gradient(np.log(l))/np.gradient(np.log(xi))+np.gradient(np.log(xi))/np.gradient(np.log(z))
    #ds=np.gradient(np.log(s))/np.gradient(np.log(xi))+np.gradient(np.log(xi))/np.gradient(np.log(z))
    #dr=np.gradient(np.log(r))/np.gradient(np.log(xi))+np.gradient(np.log(xi))/np.gradient(np.log(z))
    #dw=np.gradient(np.log(w))/np.gradient(np.log(xi))+np.gradient(np.log(xi))/np.gradient(np.log(z))

    l=(xi)**ls[SL]
    s=(xi)**hs[V]
    r=(xi)**rp[M]
    w=(xi)**wp[G]
    II=isoii[:, 61]
    Ia=iso1A[:, 61]
    g=iso_gam[:, 176]

    dl=np.gradient(np.log(l))/np.gradient(np.log(z))
    ds=np.gradient(np.log(s))/np.gradient(np.log(z))
    dr=np.gradient(np.log(r))/np.gradient(np.log(z))
    dw=np.gradient(np.log(w))/np.gradient(np.log(z))
    dII=np.gradient(np.log(II))/np.gradient(np.log(z))
    dIa=np.gradient(np.log(Ia))/np.gradient(np.log(z))
    dg=np.gradient(np.log(g))/np.gradient(np.log(z))

    ind = z0
    find = fz0

    print('ls=',dl[ind])
    print('hs=',ds[ind])
    print('r=',dr[ind])
    print('w=',dw[ind])
    print('II=',dII[ind])
    print('Ia=',dIa[ind])
    print('g=',dg[ind])

    dl_fe=np.gradient(np.log(l))/np.gradient(np.log(feP))
    ds_fe=np.gradient(np.log(s))/np.gradient(np.log(feP))
    dr_fe=np.gradient(np.log(r))/np.gradient(np.log(feP))
    dw_fe=np.gradient(np.log(w))/np.gradient(np.log(feP))
    dII_fe=np.gradient(np.log(II))/np.gradient(np.log(feP))
    dIa_fe=np.gradient(np.log(Ia))/np.gradient(np.log(feP))
    dg_fe=np.gradient(np.log(g))/np.gradient(np.log(feP))

    print('ls_fe=',dl_fe[find])
    print('hs_fe=',ds_fe[find])
    print('r_fe=',dr_fe[find])
    print('w_fe=',dw_fe[find])
    print('II_fe=',dII_fe[find])
    print('Ia_fe=',dIa_fe[find])
    print('g_fe=',dg_fe[find])

    if psmake == 1:
        closeps

    print('[GCE] line 1159: Do not continue, yr is now wrong') 
#re-call gce_yr function, since yrange was changed for process plot
#yr=gce_yr(ggg)
# END PROCESS PLOT

#CHOOSE VARIABLE FOR PLOTS

IF ggg EQ 'pb' THEN BEGIN
y_model = pb_fe_model
y_title = '[Pb/Fe]'
ENDIF

IF ggg EQ 'mg' THEN BEGIN
y_model = mg_fe_model
y_title = '[Mg/Fe]'
ENDIF

IF ggg EQ 'ba' THEN BEGIN
y_model = ba_fe_model
y_title = '[Ba/Fe]'
ENDIF

IF ggg EQ 'o' THEN BEGIN
y_model = o_fe_model
y_title = '[O/Fe]'
ENDIF

IF ggg EQ 'eu' THEN BEGIN
y_model = eu_fe_model
y_title = '[Eu/Fe]'
ENDIF

IF ggg EQ 'mn' THEN BEGIN
y_model = mn_fe_model
y_title = '[Mn/Fe]'
ENDIF

IF ggg EQ 'sr' THEN BEGIN
y_model = sr_fe_model
y_title = '[Sr/Fe]'
ENDIF

IF ggg EQ 'zr' THEN BEGIN
y_model = zr_fe_model
y_title = '[Zr/Fe]'
ENDIF'''


'''chi_squared = dblarr(nx)
;for i =0,nx-1 do begin
chi_squared(*) = 0.D0
;endfor

y_model = y_model
x_model = fe_h_model
x_data = x
y_data = av
y_err_data = sig

;w1=tanh(w)
w1=w

forward_function chi_compute
chi_squared = chi_compute(x_data, y_data, y_err_data, x_model, y_model, nx, w)
red_chi = total(chi_squared)/(nx - 2 - 1 )

;if ggg eq 0 then param(*,N) = [[a(D)], [b(U)], [red_chi]]
;if ggg eq 1 then param(*,N) = [[sp(V)], [rp(M)], [red_chi]]
if ggg eq 'ba' then param(*,N) = [[hs(V)], [red_chi]]

if ggg eq 'eu' then param(*,N) = [[rp(M)], [red_chi]]

;f_space param
;param(*,N) = [[fp(F)],[red_chi]]
if ggg eq 'mg' or ggg eq 'o' then param(*,N) = [[a(D)], [b(U)],[fp(F)],[red_chi]]

if ggg eq 'sr' then param(*,N) = [[wp(G)],[ls(SL)],[red_chi]]
if ggg eq 'zr' then param(*,N) = [[wp(G)],[red_chi]]

;if ggg eq 'pb' then param(*,N) = [[ss(ST)],[red_chi]]
if ggg eq 'pb' then param(*,N) = [[agba(AGA)],[agbb(AGB)],[cof(coo)],[red_chi]]


astr = string(a(D))
bstr = string(b(U))
cstr = string(red_chi)
sstr = string(hs(V))
rstr = string(rp(M))
name = sstr+''+rstr+''+cstr

;PRINT,'-------------------'
;if N MOD 50 eq 0 then PRINT,N, param[*,N] 
;PRINT,param(*,N)
;print,red_chi
;PRINT,'-------------------'

N=N+1


if N eq 1 then begin
erase
;plot,x_model,alog10(element_model(pb,*))
;plot,x_model,alog10(iso_ss(pb204,*)/1.8872500e-12)-x_model
plot,x_model,y_model,yrange=yr,xrange=[-5.0,0.5]
oplot,feh,mgfe,psym=2
;plot,x_model,alog10(element_model(pb,*))-x_model,yrange=yr,xrange=[-5.0,0.5]
;oplot,feh,mgfe,psym=2
endif 

if N MOD 50 eq 0 then begin
;oplot,x_model,y_model,line=2
endif

if N gt 1 then begin
;plot,x_model,alog10(iso_ss(pb204,*)/1.8872500e-12)-x_model
;plot,x_model,alog10(iso_ss(pb204,*)/1.8872500e-12)-x_model
;plot,x_model,alog10(element_model(pb,*)),line=2 
oplot,x_model,y_model,line=2
endif

;if testing eq 1 then stop

;ENDFOR
;ENDFOR
;ENDFOR

;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ftmax=max(ftot)

yy=REPLICATE(1,N_ELEMENTS(x))#y
ya=av#REPLICATE(1,N_ELEMENTS(y))
yv=var#REPLICATE(1,N_ELEMENTS(y))
ww=w#REPLICATE(1,N_ELEMENTS(y))

;ff=EXP(-(yy-ya)^2/yv)*ww
ff=EXP(-0.5*(yy-ya)^2/yv)*ww
fmax=max(ff)
ff/=fmax
; ff/=ftmax

ff=ALOG10(ff>10.D0^mrange[1])/(-mrange[1])+1
colorfunc,colorscheme,ff,c1

xbin=[x[0]-xres,x]+0.5*xres
ybin=[y[0]-yres,y]+0.5*yres

xoverlap=0.05
yoverlap=0.05

ftot/=ftmax
;ftot/=fmax
ftot=ALOG10(ftot>10.D0^mrange[1])/(-mrange[1])+1
colorfunc,colorscheme,ftot,c


x=xr[0]+(dindgen((xr[1]-xr[0])/xres)+0.5)*xres
y=yr[0]+(dindgen((yr[1]-yr[0])/yres)+0.5)*yres

!P.BACKGROUND = 'FFFFFF'XL   ;Set background color to white
!P.COLOR = 'OOOOOO'XL  

;!P.BACKGROUND = rgb(1.)     ;Set background color to white
;!P.COLOR = rgb(0)           ;Set plotline to black
; erase'''
gce_data(1)