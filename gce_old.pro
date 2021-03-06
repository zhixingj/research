PRO gce

;erase
;!P.BACKGROUND = 'FFFFFF'XL   ;Set background color to white
;!P.COLOR = 'OOOOOO'XL  
;erase

PATH='~/Desktop/gce/'
ggg=''
;read,ggg,prompt='Enter: 0-Mg, 1-Ba, 2-O, 3-Eu, 4-Mn, 5-Sr, 6-Zr, 7-process plot, 8-colorfolder (for pplot), 9-Pb: '
read,ggg,prompt='Enter element (or "color", "process", or "multi"): '
read,psmake,prompt='Enter: 0-no ps file, 1-generate ps file: '
;read,pcompare,prompt='Enter: 0-no prantzos comparison, 1-compare with prantzos model: '
;psmake=0
pcompare=0

caldat,systime(/julian),m,d,y
m=strsplit(string(m),/extract)
d=strsplit(string(d),/extract)
y=strsplit(string(y),/extract)
date=m+'_'+d+'_'+y

if psmake eq 1 then begin
  if ggg eq 'mg' then tit='Mg'
  if ggg eq 'ba' then tit='Ba'
  if ggg eq 'o' then tit='O'
  if ggg eq 'eu' then tit='Eu'
  if ggg eq 'mn' then tit='Mn'
  if ggg eq 'sr' then tit='Sr'
  if ggg eq 'zr' then tit='Zr'
  if ggg eq 'pb' then tit='Pb'
endif
  
;ggg=1
freb=1
;psmake=0

if pcompare eq 1 then begin
if ggg ne 1 and ggg ne 3 then pran=prantzos(ggg)
endif

;if ggg lt 7.0 then data_el=gce_data(ggg,freb)
;if ggg ge 7.00000 then data_el=gce_data(0,freb)
;if ggg eq 9 then data_el=gce_data(ggg,freb)

if freb eq 1 then begin
if ggg eq 'mg' then data_el=gce_data(ggg,freb)
if ggg ne 'mg' then data_el=gce_data_alt(ggg,freb)
endif

if freb eq 0 then data_el=gce_data_old(ggg,freb)

mgfe=data_el[0,*]
feh=data_el[1,*]
mgfe_sig=data_el[2,*]
feh_sig=data_el[3,*]

;zone=where(feh[*] ge -3.1 and feh[*] le -2.9 and mgfe[*] lt 1.0)
;print,'mean: ',mean(mgfe[zone])
;plot,feh,mgfe,psym=2
;stop

yr=gce_yr(ggg)
;yr=[-3.0,3.0]


x=where(mgfe_sig EQ 0)
IF x[0] NE -1 THEN BEGIN
    mgfe_sig[x]=feh_sig[x]
ENDIF

x=WHERE((mgfe_sig NE 0) OR (feh_sig NE 0))
feh     =feh     [x]
mgfe    =mgfe    [x]
feh_sig =feh_sig [x]
mgfe_sig=mgfe_sig[x]

nstar=N_ELEMENTS(feh)

xres=0.1D0/5
yres=0.02D0/5

xr=[-5,1]

x=xr[0]+(dindgen((xr[1]-xr[0])/xres)+0.5D0)*xres
y=yr[0]+(dindgen((yr[1]-yr[0])/yres)+0.5D0)*yres

nx = N_ELEMENTS(x)

ftot=x#y
ftot[*]=0.D0

FOR i=0,nstar-1 DO BEGIN

    x0=feh[i]
    y0=mgfe[i]

    sigx=feh_sig[i]
    sigy=mgfe_sig[i]

     fx=exp(-0.5D*((x-x0)/sigx)^2)
     fy=exp(-0.5D*((y-y0)/sigy)^2)

    f=fx#fy
    f/=TOTAL(f)

    ftot+=f
ENDFOR

colbarwidth=6
colbarmargin=7
colbardist=3

;backup
;xmargin=[6,.2+colbarwidth+colbarmargin+colbardist]
;ymargin=[4.2,0.2]
;xmargin=[8,.2+colbarwidth+colbarmargin+colbardist]
;ymargin=[4.2,0.2]
;xmargin=[0,.2+colbarwidth+colbarmargin+colbardist]
xmargin=[2.5,0]
ymargin=[2.5,0]
colorscheme='test26'
mrange=[0,-3]

w=total(ftot,2)
s1=y##ftot
s2=y^2##ftot
av=s1/w
var=s2/w-av^2
sig=SQRT(var)
w=tanh(w)

;!PATH = '~/dhome/Desktop/Isotopes/GCE'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE VARIABLES AND ARRAYS;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

numiso = 287
num_el = 83
z_max = 300
logz = dindgen(z_max)/((z_max-1)/5.5) - 5
;logz = dindgen(z_max)/((z_max-1)/6.9) - 6.9
;logz=x
z = 10^logz
zero = where(abs(z - 1) eq min(abs(z - 1)))
z(zero)=1.D0

iso_ws = DBLARR(numiso,z_max)
iso_s = DBLARR(numiso,z_max)
iso_ss = DBLARR(numiso,z_max)
iso_ls = DBLARR(numiso,z_max)
iso_hs = DBLARR(numiso,z_max)
iso_r = DBLARR(numiso,z_max)
iso_nup = DBLARR(numiso,z_max)
iso_gam = DBLARR(numiso,z_max)

iso_names = STRARR(1,numiso)
iso_as = STRARR(1,numiso)
iso_zs = STRARR(1,numiso)
sollo09 = DBLARR(numiso)
element_names = STRARR(1,num_el)
data = DBLARR(13,numiso)
iso = DBLARR(numiso,z_max)
element_model = DBLARR(num_el, z_max)
element_model_snii = DBLARR(num_el, z_max)
element_model_sn1a = DBLARR(num_el, z_max)
element_lodders = DBLARR(num_el)
hsproc = DBLARR(num_el, z_max)
lsproc = DBLARR(num_el, z_max)
ssproc = DBLARR(num_el, z_max)
sproc = DBLARR(num_el, z_max)
wsproc = DBLARR(num_el, z_max)
rproc = DBLARR(num_el, z_max)
nproc = DBLARR(num_el, z_max)
gproc = DBLARR(num_el, z_max)
IIproc = DBLARR(num_el, z_max)
Iproc = DBLARR(num_el, z_max)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFINE VARIABLES AND ARRAYS;END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;
;MAKE ARRAYS;
;;;;;;;;;;;;;

path = '~/Desktop/gce/'

OPENR, lun, path+'isotope_names.dat', /GET_LUN
READF, lun, iso_names

close,lun
free_lun,lun

OPENR, lun, path + 'isotope_As.dat', /GET_LUN
READF, lun, iso_as
close,lun
free_lun,lun

iso_list = iso_As + '_' + iso_names

OPENR, lun, path + 'isotope_Zs.dat', /GET_LUN
READF, lun, iso_zs
close,lun
free_lun,lun

OPENR,lun,'~/Desktop/gce/isotope_names.dat',/GET_LUN
isoname=STRARR(287)
READF,lun,isoname
close,lun
free_lun,lun

ion=string(uint(iso_as))
for i=0,numiso-1 do begin
ion[i]=isoname[i]+strsplit(ion[i],/extract)
endfor

close,lun
free_lun,lun

;OPENW,lun,'~/Desktop/gce/ions.dat',/GET_LUN
;for i=0,n_elements(ion)-1 do begin 
;PRINTF,lun,ion[i]
;endfor
;close,lun
;free_lun,lun
;stop

OPENR, lun, path + 'element_names.dat', /GET_LUN
header3 = STRARR(1)
READF,lun, header3
READF, lun, element_names
close,lun
free_lun,lun

;OPENR,lun, path + 'solar_abun.dat', /GET_LUN
;header2 = STRARR(1)
;READF,lun, header2
;data4 = DBLARR(3,numiso)
;READF,lun,data4

;sollo09 = data4(2,*)

OPENR,lun, path + 'sollo09.dat', /GET_LUN
header2 = STRARR(28)
READF,lun, header2
data4 = STRARR(numiso)
READF,lun,data4

;sollo09 is in mass fractions
sollo09 = strmid(data4,6,12)
sollo09=double(sollo09)

close,lun
free_lun,lun

;convert to mol frac
OPENR,lun, path + 'isomasses_numbers.dat', /GET_LUN
header = STRARR(5)
READF,lun, header
amu = dblarr(numiso)
READF,lun,amu
close,lun
free_lun,lun

sollo09[*]=sollo09[*]/amu[*]

fe56=where(ion eq 'Fe56')
mg24=where(ion eq 'Mg24')
o16=where(ion eq 'O16')

;;;;;;;;;;;;;;;;
;GET HEGER DATA;
;;;;;;;;;;;;;;;;

;NEW starfit DATA 6/13/2012 (IS IN MOL FRACTIONS!!!) 
OPENR,lun,'~/Desktop/gce/starfit_274823_6_13_2013.dat',/GET_LUN
;OPENR,lun,'~/Desktop/gce/cayrelfit_3_8_2012.dat',/get_lun
;heads=strarr(4)
;readf,lun,heads
data3 = DBLARR(286)
READF,lun,data3
massive = DBLARR(numiso)
massive(0:285) = data3(0:285)
close,lun
free_lun,lun

;ONLY USE C12 TO GA69 HEGERDATA
massive(0:8) = 0.0
massive(76:284) = 0.0

;plot,iso_as[*],alog10(massive[*]/sollo09[*]),psym=2,xrange=[8,80],xs=1
;hline,alog10(massive[fe56]/sollo09[fe56])
;hline,alog10(massive[o16]/sollo09[o16])
;vline,56
;vline,16
;stop
;binwrite,array=massive,filename='massivePRE'
;stop

;;;;;;;;;;;;;;;;
;END HEGER DATA;
;;;;;;;;;;;;;;;;

h = where(element_names eq 'H')
c = where(element_names eq 'C')
nit = where(element_names eq 'N')
o = where(element_names eq 'O')
fl = where(element_names eq 'F')
neon = where(element_names eq 'Ne')
na = where(element_names eq 'Na')
mg = where(element_names eq 'Mg')
al = where(element_names eq 'Al')
si = where(element_names eq 'Si')
ca = where(element_names eq 'Ca')
pp = where(element_names eq 'P')
ss = where(element_names eq 'S')
cl = where(element_names eq 'Cl')
ar = where(element_names eq 'Ar')
pot = where(element_names eq 'K')
sc = where(element_names eq 'Sc')
ti = where(element_names eq 'Ti')
va = where(element_names eq 'V')
cr = where(element_names eq 'Cr')
mn = where(element_names eq 'Mn')
fe = where(element_names eq 'Fe')
co = where(element_names eq 'Co')
ni = where(element_names eq 'Ni')
cu = where(element_names eq 'Cu')
zn = where(element_names eq 'Zn')
ga = where(element_names eq 'Ga')
ger = where(element_names eq 'Ge')
as = where(element_names eq 'As')
se = where(element_names eq 'Se')
Br = where(element_names eq 'Br')
Kr = where(element_names eq 'Kr')
rb = where(element_names eq 'Rb')
sr = where(element_names eq 'Sr')
yt = where(element_names eq 'Y')
zr = where(element_names eq 'Zr')
nb = where(element_names eq 'Nb')
pb = where(element_names eq 'Pb')
mo = where(element_names eq 'Mo')
pd = where(element_names eq 'Pd')
ag = where(element_names eq 'Ag')
uran = where(element_names eq 'U')
th = where(element_names eq 'Th')
au = where(element_names eq 'Au')
pt = where(element_names eq 'Pt')
gd = where(element_names eq 'Gd')
hf = where(element_names eq 'Hf')
ho = where(element_names eq 'Ho')
ir = where(element_names eq 'Ir')
lu = where(element_names eq 'Lu')
dy = where(element_names eq 'Dy')
os = where(element_names eq 'Os')

ba = where(element_names eq 'Ba')
eu = where(element_names eq 'Eu')
er = where(element_names eq 'Er')

;;;;;;;;;;;;;;;;;;;
;DEFINE PARAMETERS;
;;;;;;;;;;;;;;;;;;;

space = 25
a = DBLARR(space)
;a(0) = 4.245D0
;a(0)=24.D0
a = dindgen(space+1)/((space+1)/0.04D0)+5.0D0
;a = dindgen(space+1)/((space+1)/0.002D0)+4.639D0

space1 = 25
b = DBLARR(space1)
;b(0) = 4.814D0
;b(0)=16.D0
b = dindgen(space1+1)/((space1+1)/0.04D0)+2.7D0
;b = dindgen(space1+1)/((space1+1)/0.001D0)+2.399D0

f_space = 10
fp = DBLARR(f_space)
fp = DINDGEN(f_space+1)/((f_space+1)/.002D0)+0.692D0 
;fp = DINDGEN(f_space+1)/((f_space+1)/0.0002D0)+0.7082D0 
;fp = 0.7083D0

space2 = 1000
hs = DBLARR(space2)
hs = dindgen(space2+1)/((space2+1)/0.07D0)+1.47D0

space3 = 1000
rp = DBLARR(space3)
rp = dindgen(space3+1)/((space3+1)/0.004D0)+0.935D0

space4 = 30
wp = DBLARR(space4)
wp = dindgen(space4+1)/((space4+1)/.02D0)+2.87D0

space6 = 30
ls = DBLARR(space6)
ls = dindgen(space6+1)/((space6+1)/0.01D0)+1.20D0

space7 = 26
agba = DBLARR(space7)
agba = dindgen(space7+1)/((space7+1)/100.D0)+1.0D0

space8 = 10
agbb = DBLARR(space8)
;agbb = dindgen(space8+1)/((space8+1)/-0.1D0)-0.001D0
agbb=[-.00001,-.00005,-.0001,-.0005,-.001,-.005,-.01,-.05,-.1,-.5]
;agbb*=100.d0

space9 = 10
cof = DBLARR(space9)
cof = dindgen(space9+1)/((space9+1)/100.D0)+1.D0

;red_chi = DBLARR(space*space1)
N = 0
a_array = DBLARR(space)
b_array = DBLARR(space1)
hs_array = DBLARR(space2)
rp_array = DBLARR(space3)
;param = DBLARR(3,space*space1)
;if ggg eq 1 then param = DBLARR(3,space2*space3)
if ggg eq 'ba' then param = DBLARR(2,space2)
if ggg eq 'eu' then param = DBLARR(2,space3)
;param = DBLARR(2,f_space)
if ggg eq 'mg' or ggg eq 'o' then param = DBLARR(4,space*space1*f_space)
if ggg eq 'sr' then param = DBLARR(3,space4*space6)
if ggg eq 'zr' then param = DBLARR(2,space4)

;if ggg eq 'pb' then param = DBLARR(2,space5)
if ggg eq 'pb' then param = DBLARR(4,space7*space8*space9)



close,lun
free_lun,lun

chimin = dblarr(3,space)
kk=0
;MASTER LOOP FOR FINDING PARAMETER SPACES
;FOR D = 0,space-1 DO BEGIN
;FOR U = 0,space1-1 DO BEGIN
;FOR F = 0,f_space -1 DO BEGIN
;FOR V = 0,space2-1 DO BEGIN
;FOR M = 0,space3-1 DO BEGIN
;FOR G = 0,space4-1 DO BEGIN
;FOR ST = 0,space5-1 DO BEGIN (obsolete)
;FOR SL = 0,space6-1 DO BEGIN
;FOR AGA = 0,space7-1 DO BEGIN
;FOR AGB = 0,space8-1 DO BEGIN
;FOR coo = 0,space9-1 DO BEGIN

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;BEST FIT VALUES FOUND BY SCANNING PARAMETER SPACE (AFTER 6/11/2013)

SL=0
ls(SL)=1.227d0
G=0
wp(G)=1.230D0

;3 PARAM
AGA=0
AGB=0
coo=0
agba(AGA)=200.D0
agbb(AGB)=-.23d0
cof(coo)=-2.e-11

F=0
fp(F)=0.6935D0

D=0
U=0
a(D)=5.024D0
b(U)=2.722D0

V=0
M=0
hs(V)=1.509D0
;hs(V)=2.0D0
rp(M)=0.938D0

;END BEST FIT VALUES (AFTER 3/8/2012)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;BEST FIT VALUES FOUND BY SCANNING PARAMETER SPACE (PRIOR TO 3/8/2012)

;SL=0
;ls(SL)=1.2276d0
;G=0
;wp(G) = 1.2287D0


;ST=0
;ss(ST) = 0.325
;AGA=0
;AGB=0
;BEST FIT
;agba(AGA)=307.7D0
;agbb(AGB)=1.623D0
;BEST BEHAVIOUR
;agba(AGA)=500.D0
;agbb(AGB)=3.397D0

;3 PARAM
;AGA=0
;AGB=0
;coo=0
;agba(AGA)=200.D0
;agbb(AGB)=-.29d0
;cof(coo)=-2.e-11

;F=0
;fp(F) = 0.709D0

;D=0
;U=0
;a(D) = 4.635D0
;b(U) = 2.397D0

;V=0
;M=0
;hs(V) = 1.487D0
;rp(M) = 0.968D0

;G=0
;wp(G) = 0.8175D0

;sp(V) = 1.52D0
;rp(M) = 1.061D0

;END BEST FIT VALUES (PRIOR TO 3/8/2012)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;
;CALL ALL SUBROUTINES;
;;;;;;;;;;;;;;;;;;;;;;

;DO SOLAR DECOMPOSITION
data=solar_decomp(sollo09,amu,iso_zs,fp(F),ion,massive)
isotopes=data

;aqw=dblarr(numiso)
;aqw[*]=isotopes[7,*]*amu[*]
;binwrite,array=aqw[*],filename='data/w_hybridsolarfit1'
;stop

;these arrays are for isobox plots
s = data(2,*)
ws = data(3,*)
r = data(4,*)
nu_p = data(5,*)
gamma = data(6,*)
one_a = data(7,*)
ii = data(8,*)
;lightest = data(9,*)
bbn=data(9,*)
gcr=data(10,*)
novae_priGCR_nu=data(11,*)
h_burn=data(12,*)

data(1,*)=iso_as(*)

;FINISH HEGER DATA
snii_z = data(0,*)
snii_a = data(1,*)

fraction=dblarr(11,287)
for i=0,numiso-1 do begin
  ticks=0
  for j=0,10 do begin
    fraction[j,i]=isotopes[2+j,i]/sollo09[i]
    if fraction[j,i] lt 0.d0 then fraction[j,i]=0.d0
    if fraction[j,i] eq 1.d0 then begin
      ticks=1
      col=j
    endif
  endfor
  if ticks eq 1 then begin
    fraction[*,i]=0.d0
    fraction[col,i]=1.d0
  endif
endfor

;zero out tiny negligible contributions
for i=0,numiso-1 do begin
  for j=0,10 do begin
    if fraction[j,i] lt 1.e-6 then fraction[j,i]=0.d0
  endfor
endfor

;RE-MAKE FRACTION INTO STRING ARRAYS
fraction2=fraction
fraction=strarr(11,287)
for i=0,numiso-1 do begin
  for j=0,10 do begin
    fraction[j,i]=strsplit(string(fraction2[j,i]),/extract)
    if fraction[j,i] eq '0.0000000' then fraction[j,i]='\nodata'
  endfor
endfor

;units are mol frac
openw,lun,'fraction_'+date+'.dat',/get_lun
printf,lun,'isotope, &, Solar Abun, &, Main S, &, Weak S, &, R, &, Nu-P, &, Gamma, &, SNIa, &, Massive, &, BBN, &, GCR Spallation, &, Novae/Pri GCR/Nu, &, H-Burn'
for i=0,numiso-1 do begin
printf,lun,ion[i],'&',sollo09[i],'&',fraction[0,i],'&',fraction[1,i],'&',fraction[2,i],'&',fraction[3,i],'&',fraction[4,i],'&',fraction[5,i],'&',fraction[6,i],'&',fraction[7,i],'&',fraction[8,i],'&',fraction[9,i],'&',fraction[10,i],format='(25(A17))'
endfor
close,lun
free_lun,lun

exportmodel=0

IF exportmodel EQ 1 THEN BEGIN

;MATRICIES FOR RECASTING MODEL AS X*F
massive_scaled=DBLARR(numiso)
openr,lun,'~/gce/massivedatascaledforminus3.dat',/get_lun
readf,lun,massive_scaled
close,lun
free_lun,lun

slopes=DBLARR(numiso)
openr,lun,'~/gce/slopeslog.dat',/get_lun
readf,lun,slopes
close,lun
free_lun,lun

ba130=where(ion eq 'Ba130')
pb204=where(ion eq 'Pb204')
bi209=where(ion eq 'Bi209')

light_s=dblarr(numiso)
light_s(0:ba130-1)=data(2,0:ba130-1)

heavy_s=dblarr(numiso)
heavy_s(ba130:pb204-1)=data(2,ba130:pb204-1)

strong_s=dblarr(numiso)
strong_s(pb204:286)=data(2,pb204:286)

strong_s_func=dblarr(numiso)
strong_s_func(pb204:bi209)=1.d0

;ADD D TO H_BURN ARRAY, IT EVOLVES NOW THE SAME (~Z)
;MUST ALSO ADD TERMS FOR H1 SCALING: (D + He3 + He4) only the z-dependence part
;RENAME ARRAY, OTHERWISE SOLAR DECOMP ISOBOX WILL BE WRONG
h_burnx=h_burn
h_burnx[1]=sollo09[1]-bbn[1]
h_burnx[0]=-1.d0*((sollo09[1]-bbn[1])*amu[1] + (sollo09[2]-bbn[2])*amu[2] + (sollo09[3]-bbn[3])*amu[3])/amu[0]

;hyd gets * 1.d0*z (total metallicity for sum)
hyd=dblarr(numiso)
hyd[0]= -1.d0*0.0153d0/amu[0]

;hterms gets * 1.d0
;put in -1.d0 (sum of X+Y+Z)
;put in Y and D the bbn part (y-intercept)
;subtract BBN for H1 (which in mol frac) -- it comes from x_matrix[14] ;& add back BBN for H1 in mass fractions
hterms=dblarr(numiso)
hterms[0]=(1.d0 - (bbn[1]*amu[1] + bbn[2]*amu[2] + bbn[3]*amu[3]))/amu[0] - bbn[0]    ;+ bbn[0]*amu[0])

;hyd=dblarr(numiso)
;hmult=(0.0153d0+h_burn(2)*3.016029d0+h_burn(3)*4.002603d0)
;hyd(0)=0.99998060d0*hmult/1.007825D0
;hyd(1)=(1.d0-0.99998060d0)*hmult/2.014102D0*2.d0

;massfrac=dblarr(numiso)
;massfrac(0)=(1.d0-3.016029d0*bbn(2)-4.002603d0*bbn(3))*0.99998060d0/1.007825D0-bbn[0]
;massfrac(1)=(1.d0-3.016029d0*bbn(2)-4.002603d0*bbn(3))*(1.d0-0.99998060d0)/2.014102D0*2.d0-bbn[1]


x_matrix=DBLARR(17,numiso)
;x_matrix[0,*]=ls
;x_matrix[1,*]=hs
;x_matrix[2,*]=ss
;x_matrix[3,*]=ws
;x_matrix[4,*]=r
;x_matrix[5,*]=nu_p
;x_matrix[6,*]=gamma
;x_matrix[7,*]=one_a
;x_matrix[8,*]=slopes
;x_matrix[9,*]=massive_scaled
;x_matrix[10,*]=gcr
;x_matrix[11,*]=novae_priGCR_nu
;x_matrix[12,*]=bbn
;x_matrix[13,*]=h_burn
x_matrix[0,*]=light_s
x_matrix[1,*]=heavy_s
x_matrix[2,*]=strong_s
x_matrix[3,*]=strong_s_func
x_matrix[4,*]=ws
x_matrix[5,*]=r
x_matrix[6,*]=nu_p
x_matrix[7,*]=gamma
x_matrix[8,*]=one_a
x_matrix[9,*]=slopes
x_matrix[10,*]=massive_scaled
x_matrix[11,*]=gcr
x_matrix[12,*]=novae_priGCR_nu
x_matrix[13,*]=bbn
x_matrix[14,*]=h_burnx
x_matrix[15,*]=hterms
x_matrix[16,*]=hyd

openw,lun,'x_matrix_'+date+'.dat',/get_lun
for i=0,0 do begin
printf,lun,x_matrix[*,*],format='(17(A17))' 
;format='(17(E17.10))' 
endfor
close,lun
free_lun,lun
stop
;if z >= -2.5 use this
;f_matrix=DBLARR(17,1)
;f_matrix[0]=z^1.227d0
;f_matrix[1]=z^1.509d0
;f_matrix[2]=1.d0
;f_matrix[3]=-2.e-11*(1-(tanh(200d0*z-0.29)/tanh(200d0*1.d0-0.29d0)))
;f_matrix[4]=z^1.230d0
;f_matrix[5]=z^.938d0
;f_matrix[6]=z^.938d0
;f_matrix[7]=z^((1.509d0+.938d0)/2)
;f_matrix[8]=z*((tanh(5.024d0*z-2.722d0))+tanh(2.722d0))/(tanh(5.024d0-2.722d0)+tanh(2.722d0))
;f_matrix[9]=(z-0.00317)
;f_matrix[10]=1.d0
;f_matrix[11]=z^1.509d0
;f_matrix[12]=z^.938d0
;f_matrix[13]=1.d0
;f_matrix[14]=z^.938d0
;f_matrix[15]=1.d0
;f_matrix[16]=1.d0*z

;if z < -2.5 use this
;f_matrix=DBLARR(17,1)
;f_matrix[0]=z^1.227d0
;f_matrix[1]=z^1.509d0
;f_matrix[2]=1.d0
;f_matrix[3]=-2.e-11*(1-(tanh(200d0*z-0.29)/tanh(200d0*1.d0-0.29d0)))
;f_matrix[4]=z^1.230d0
;f_matrix[5]=z^.939d0
;f_matrix[6]=z^.939d0
;f_matrix[7]=z^((1.509d0+.938d0)/2)
;f_matrix[8]=z*((tanh(5.024d0*z-2.722d0))+tanh(2.722d0))/(tanh(5.024d0-2.722d0)+tanh(2.722d0))
;f_matrix[9]=0.d0
;f_matrix[10]=z/0.00317
;f_matrix[11]=z^1.509d0
;f_matrix[12]=z^.938d0
;f_matrix[13]=1.d0
;f_matrix[14]=z^.938d0
;f_matrix[15]=1.d0
;f_matrix[16]=1.d0*z

ENDIF

forward_function sn1a
iso1A = sn1a(one_a, z, z_max, b(U), numiso, a(D))

;forward_function snii
isoii = snii(ii, z, z_max, snii_z, snii_a, numiso, massive, ion, fp(F))

;m=dblarr(numiso)
;m[*]=massive_scaled[*]*amu[*]
;sollo09[*]*=amu[*]
;print,alog10(total(m[60:63])/total(sollo09[60:63]))-alog10(total(m[13:15])/total(sollo09[13:15]))

;stop

zero = where(abs(z - 1.0) eq min(abs(z - 1.0)))
zm3 = where(abs(alog10(z)+2.5) eq min(abs(alog10(z)+2.5)))
zm1 = where(abs(alog10(z)+0.5) eq min(abs(alog10(z)+0.5)))
zm2 = where(abs(alog10(z)+1.5) eq min(abs(alog10(z)+1.5)))
zm4 = where(abs(alog10(z)+3.5) eq min(abs(alog10(z)+3.5)))

forward_function light
;make lightest array for evolving 1<=A<=11
lightest=(bbn+gcr+novae_priGCR_nu+h_burn)
lightest[0:1]=sollo09[0:1]
iso_lt = light(lightest, z, z_max, numiso, hs(V), rp(M))

ba130=where(ion eq 'Ba130')
pb204=where(ion eq 'Pb204')

;data_str=DBLARR(numiso)
;data_str(pb204:286)=data(2,pb204:286)
;data_ss=agbstr(data_str, z, z_max, agbb(AGB),numiso,agba(AGA))

; MAKE S R NU_P GAMMA ARRAYS
for i = 0l, z_max-1 do begin
   ;Pb204-208 are i=278-281
   iso_ls(0:ba130-1,i) = data(2,0:ba130-1)*(z(i))^(ls(SL))
   iso_hs(ba130:pb204-1,i) = data(2,ba130:pb204-1)*(z(i))^(hs(V))
;   iso_ss(pb204:286,i) = data(2,pb204:286)*(z(i))^(ss(ST)) 
   ;iso_s(0:277,i) = data(2,0:277)*(z(i))^(sp(V))
   ;iso_s(282:286,i) = data(2,282:286)*(z(i))^(sp(V))
   iso_ws(*,i) = data(3,*)*(z(i))^(wp(G))
   iso_r(*,i) = data(4,*)*(z(i))^(rp(M))
   iso_nup(*,i) = data(5,*)*(z(i))^(rp(M))
   iso_gam(*,i) = data(6,*)*(z(i))^((hs(V)+rp(M))/2.D0+1.D0)
endfor

z_str_cutoff=closest(logz,-2.255d0)

;agba=10.
;agbb=2.
;f=tanh(agbb)
;g=tanh(agba-agbb)
;plot,logz,alog10((((tanh(agba*z-agbb))+f)/(g+f)))

;AGA=0
;AGB=0
;agba(AGA)=10.D0
;agbb(AGB)=2.397D0

testing = 0
if testing eq 1 then begin
AGA=0
AGB=0
coo=0
agba(AGA)=200.D0
agbb(AGB)=-.23d0
cof(coo)=-2.e-11
endif

z_solar=closest(logz,0.d0)
;f=tanh(agbb(AGB))
;g=tanh(agba(AGA)*z(z_solar)+agbb(AGB))
;plot,x_model,alog10(iso_ss(pb204,*)/element_model(fe,*))

;zl=closest(logz,-4.5d0)
;z_solar=closest(logz,0.d0)
;;yl=data(2,pb204)/10000.d0
;f=tanh(agba(AGA)*z(z_solar)+agbb(AGB))
;g=tanh(agba(AGA)*z(zl)+agbb(AGB))
;mult=1.e2

;redo strong, cutoff at logz=-2.5
for i = 0l, z_max-1 do begin
;iso_ss(pb204:286,i) = (data(2,pb204:286)/mult-data(2,pb204:286))*((tanh(agba(AGA)*z(i)-agbb(AGB)))-f)/(g-f)+data(2,pb204:286)
;iso_ss(pb204:286,i) = data(2,pb204:286)*(tanh(agba(AGA)*z(i)+agbb(AGB))-f)/(g-f)
;iso_ss(pb204:286,i) = cof(coo)*(tanh(agba(AGA)*z(i)+agbb(AGB))-tanh(agba(AGA)*z(z_solar)+agbb(AGB)))+data(2,pb204:286)

;THESIS UP TO DEFENSE VALUES
;iso_ss(pb204:286,i) = cof(coo)*(1-(tanh(agba(AGA)*z(i)+agbb(AGB))/tanh(agba(AGA)*z(z_solar)+agbb(AGB))))+data(2,pb204:286)

agba(AGA)=140.
agbb(AGB)=-0.05
;200,-0.8,-0.23

;BETTER FIT VALUES
;190,-0.1,-0.05

;BEST FIT 2 PARAMETERS 7 TYPE IA BASED FUNCTION
agba(AGA)=140.
agbb(AGB)=-0.05
iso_ss(pb204:286,i) = data(2,pb204:286)*(tanh(agba(AGA)*z(i)+agbb(AGB))+tanh(agbb(AGB)))/(tanh(agba(AGA)+agbb(AGB))+tanh(agbb(AGB)))

;iso_ss(pb204:286,i) = data(2,pb204:286)*(tanh(agba(AGA)*z(i)-agbb(AGB))+tanh(agbb(AGB)))/(tanh(agba(AGA)+agbb(AGB)))
;iso_ss(pb204:286,i) = data(2,pb204:286)*(tanh(agba(AGA)*z(i)-agbb(AGB))+tanh(-.05))/(tanh(agba(AGA))-tanh(agbb(AGB)))

endfor

;using zl and yl (AND 3 PARAM) technique, we now do not have BBN boundary conditions, so reassign negative values
for i = 0l, z_max-1 do begin
  for j=pb204[0], 286 do begin
    if iso_ss(j,i) lt 0.d0 then iso_ss(j,i) = 0.d0
  endfor
endfor

isomake=0
if isomake eq 1 then begin
boxx=isotopes
boxx[1,*]=iso_as[*]
openw,lun,'~/Desktop/gce/isobox/boxmatrix'+date+'.dat',/GET_LUN
for i=0,1 do begin
  printf,lun,boxx[*,*],format='(13(A17))'
;THIS PRINTS TWO COPIES OF ISOTOPES INTO BOXMATRIX.DAT
;SO YOU NEED TO GO INTO FILE AND MANUALLY ERASE THE SECOND!!
endfor
print,'[GCE] EXPORTED TO ISOBOX FOR BOXMATRIX.DAT'
STOP
endif

;;;;;;;;;;;;;;;;;
;END SUBROUTINES;
;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;MAKE ISO/ELEMENT MATRICIES;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

iso = DBLARR(numiso,z_max)
iso = isoii + iso1A + iso_ws + iso_r + iso_nup + iso_gam + iso_lt + iso_ss + iso_ls + iso_hs
model_massfrac = DBLARR(numiso,z_max)
lodders_massfrac = DBLARR(numiso)

for i=0,numiso-1 do begin
;model_massfrac(i,*) = iso(i,*)*DOUBLE(iso_As(i))
;lodders_massfrac(i) = sollo09(i)*DOUBLE(iso_As(i)) 
model_massfrac(i,*) = iso(i,*)*amu[i]
lodders_massfrac(i) = sollo09(i)*amu[i]
endfor

;;scale linearly between bbn and solar
;s=dblarr(numiso)
;scaled=dblarr(300,numiso)
;sc_em=dblarr(300,83)
;for i=0,numiso-1 do begin
;   s[i]=(sollo09[i]-bbn[i])/(1.d0-0.d0)
;endfor 

;for i=0,n_elements(logz)-1 do begin
;for j=0,numiso-1 do begin
;   scaled[i,j]=s[j]*(z[i]-1.d0)+sollo09[j]
;endfor
;endfor

i = 0
j = 0
while i lt numiso do begin
   index = where(iso_names(i) eq iso_names)
   
   element_model(j,*) = total(iso(index,*),1)
   model_massfrac(j,*) = total(model_massfrac(index,*),1)
   lodders_massfrac(j) = total(lodders_massfrac(index),1)
   lsproc(j,*) = total(iso_ls(index,*),1)
   hsproc(j,*) = total(iso_hs(index,*),1)
   ssproc(j,*) = total(iso_ss(index,*),1)
   sproc(j,*) = total(iso_s(index,*),1)
   rproc(j,*) = total(iso_r(index,*),1)
   nproc(j,*) = total(iso_nup(index,*),1)
   gproc(j,*) = total(iso_gam(index,*),1)
   wsproc(j,*) = total(iso_ws(index,*),1)
   IIproc(j,*) = total(isoII(index,*),1)
   Iproc(j,*) = total(iso1A(index,*),1)
   element_lodders(j) = total(sollo09(index),1)
   element_model_snii(j,*) = total(isoII(index,*),1)
   element_model_sn1a(j,*) = total(iso1A(index,*),1)

;   sc_em(*,j) = total(scaled(*,index))

   j = j + 1
   ;i = i + n_elements(index)
   i+=n_elements(index)
endwhile

;stop

;MAKE ELEMENT MODEL FOR SCALED SOLAR
scaled_solar=dblarr(300,83)
scaled_feh=dblarr(300)
for j=0,n_elements(element_model(*,0))-1 do begin ;elements
for i=0,n_elements(element_model(0,*))-1 do begin ;z
  ;if j eq 0 then scaled_solar[i,j]=-element_lodders[j]*z[i]+bbn[0]+bbn[1]
  if j eq 0 then scaled_solar[i,j]=(element_lodders[j]-bbn[0]-bbn[1])*z[i];+bbn[0]+bbn[1]
  if j eq 1 then scaled_solar[i,j]=(element_lodders[j]-bbn[2]-bbn[3])*z[i];+bbn[2]+bbn[3]
  if j eq 2 then scaled_solar[i,j]=(element_lodders[j]-bbn[4]-bbn[5])*z[i];+bbn[4]+bbn[5]
  if j gt 2 then scaled_solar[i,j]=element_lodders[j]*z[i]
endfor
endfor

for i=0,n_elements(element_model(0,*))-1 do begin ;z
  scaled_feh[i]=alog10(scaled_solar[i,25]/element_lodders[25]) - alog10(scaled_solar[i,0]/element_lodders[0])
endfor

;EXPORT TO COLOR FOLDER
scalecolor=0
if scalecolor eq 1 then begin
data_el = scaled_solar
openw,lun,'~/Desktop/gce/pplot/colordat/smatrix'+date+'.dat'
for i=0,82 do begin
   printf,lun,scaled_solar(*,i)
endfor

close,lun
free_lun,lun

openw,lun,'~/Desktop/gce/pplot/colordat/fe_h_s'+date+'.dat'
for i = 0,z_max-1 do printf,lun,scaled_feh(i)

close,lun
free_lun,lun
stop
endif

close,lun
free_lun,lun

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; END ISO/ELEMENT MATRICIES;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; MAKE ARRAYS AND CALL CHI COMPUTE

fe_h_model = DBLarr(z_max)
mg_fe_model = DBLarr(z_max)
mn_fe_model = DBLarr(z_max)
na_fe_model = DBLarr(z_max)
ba_fe_model = DBLarr(z_max)
ni_fe_model = DBLarr(z_max)
eu_fe_model = DBLarr(z_max)

sr_fe_model = DBLarr(z_max)
zr_fe_model = DBLarr(z_max)
metallicity = DBLarr(z_max)
mg_h_model = DBLarr(z_max)
o_h_model = DBLarr(z_max)
ba_h_model = DBLarr(z_max)
sba_h_model = DBLarr(z_max)
rba_h_model = DBLarr(z_max)
gba_h_model = DBLarr(z_max)
eu_h_model = DBLarr(z_max)
mn_h_model = DBLarr(z_max)
ni_h_model = DBLarr(z_max)
solarz = TOTAL(lodders_massfrac(2:82),1)
ba_fe_model_hs = DBLarr(z_max)
ba_fe_model_r = DBLarr(z_max)
ba_fe_model_g = DBLarr(z_max)
sr_fe_model_ls = DBLarr(z_max)
sr_fe_model_r = DBLarr(z_max)
sr_fe_model_ws = DBLarr(z_max)
sr_ba_model = DBLarr(z_max)
eu_fe_model_s = DBLarr(z_max)
eu_fe_model_r = DBLarr(z_max)
eu_fe_model_g = DBLarr(z_max)
mg_fe_model_ii = DBLarr(z_max)
mg_fe_model_1a = DBLarr(z_max)
o_fe_model_ii = DBLarr(z_max)
o_fe_model_1a = DBLarr(z_max)
mn_fe_model_ii = DBLarr(z_max)
mn_fe_model_1a = DBLarr(z_max)

c_fe_model = DBLARR(z_max)
n_fe_model = DBLARR(z_max)
o_fe_model = DBLARR(z_max)
f_fe_model = DBLARR(z_max)
ne_fe_model = DBLARR(z_max)
na_fe_model = DBLARR(z_max)
mg_fe_model = DBLARR(z_max)
al_fe_model = DBLARR(z_max)
si_fe_model = DBLARR(z_max)
p_fe_model = DBLARR(z_max)
s_fe_model = DBLARR(z_max)
cl_fe_model = DBLARR(z_max)
ar_fe_model = DBLARR(z_max)
k_fe_model = DBLARR(z_max)
ca_fe_model = DBLARR(z_max)
ir_fe_model = DBLARR(z_max)
er_fe_model = DBLARR(z_max)
sc_fe_model = DBLARR(z_max)
ti_fe_model = DBLARR(z_max)
v_fe_model = DBLARR(z_max)
cr_fe_model = DBLARR(z_max)
mn_fe_model = DBLARR(z_max)
co_fe_model = DBLARR(z_max)
os_fe_model = DBLARR(z_max)
ni_fe_model = DBLARR(z_max)
cu_fe_model = DBLARR(z_max)
zn_fe_model = DBLARR(z_max)
ga_fe_model = DBLARR(z_max)
ge_fe_model = DBLARR(z_max)
as_fe_model = DBLARR(z_max)
se_fe_model = DBLARR(z_max)
br_fe_model = DBLARR(z_max)
kr_fe_model = DBLARR(z_max)
rb_fe_model = DBLARR(z_max)
sr_fe_model = DBLARR(z_max)
y_fe_model = DBLARR(z_max)
zr_fe_model = DBLARR(z_max)
pb_fe_model = DBLARR(z_max)
pb_fe_model_ss = DBLARR(z_max)
pb_fe_model_r = DBLARR(z_max)
ag_fe_model = DBLARR(z_max)
al_fe_model = DBLARR(z_max)
au_fe_model = DBLARR(z_max)
cr_fe_model = DBLARR(z_max)
dy_fe_model = DBLARR(z_max)
gd_fe_model = DBLARR(z_max)
hf_fe_model = DBLARR(z_max)
ho_fe_model = DBLARR(z_max)
lu_fe_model = DBLARR(z_max)
mo_fe_model = DBLARR(z_max)
pd_fe_model = DBLARR(z_max)
pt_fe_model = DBLARR(z_max)
u_fe_model = DBLARR(z_max)
cr_fe_model_ii = DBLARR(z_max)
cr_fe_model_ia = DBLARR(z_max)
mn_fe_model_ii = DBLARR(z_max)
mn_fe_model_ia = DBLARR(z_max)

fe_model=DBLARR(z_max)
sr_model=DBLARR(z_max)
pb_model=DBLARR(z_max)
mg_model=DBLARR(z_max)
o_model=DBLARR(z_max)
eu_model=DBLARR(z_max)

fe_modelII=DBLARR(z_max)
sr_modelII=DBLARR(z_max)
pb_modelII=DBLARR(z_max)
mg_modelII=DBLARR(z_max)
o_modelII=DBLARR(z_max)
eu_modelII=DBLARR(z_max)
eu_modelR=DBLARR(z_max)

sr_modelW=DBLARR(z_max)

ba_mg_model=DBLARR(z_max)
eu_mg_model=DBLARR(z_max)

feP=DBLARR(z_max)

for i = 0, z_max-1 do begin
c_fe_model(i) = alog10((element_model(c,i)/element_lodders(c))/(element_model(fe,i)/element_lodders(fe)))
n_fe_model(i) = alog10((element_model(nit,i)/element_lodders(nit))/(element_model(fe,i)/element_lodders(fe)))
o_fe_model(i) = alog10((element_model(o,i)/element_lodders(o))/(element_model(fe,i)/element_lodders(fe)))
f_fe_model(i) = alog10((element_model(fl,i)/element_lodders(fl))/(element_model(fe,i)/element_lodders(fe)))
ne_fe_model(i) = alog10((element_model(neon,i)/element_lodders(neon))/(element_model(fe,i)/element_lodders(fe)))
na_fe_model(i) = alog10((element_model(na,i)/element_lodders(na))/(element_model(fe,i)/element_lodders(fe)))
mg_fe_model(i) = alog10((element_model(mg,i)/element_lodders(mg))/(element_model(fe,i)/element_lodders(fe)))
al_fe_model(i) = alog10((element_model(al,i)/element_lodders(al))/(element_model(fe,i)/element_lodders(fe)))
si_fe_model(i) = alog10((element_model(si,i)/element_lodders(si))/(element_model(fe,i)/element_lodders(fe)))
p_fe_model(i) = alog10((element_model(pp,i)/element_lodders(pp))/(element_model(fe,i)/element_lodders(fe)))
s_fe_model(i) = alog10((element_model(ss,i)/element_lodders(ss))/(element_model(fe,i)/element_lodders(fe)))
cl_fe_model(i) = alog10((element_model(cl,i)/element_lodders(cl))/(element_model(fe,i)/element_lodders(fe)))
ar_fe_model(i) = alog10((element_model(ar,i)/element_lodders(ar))/(element_model(fe,i)/element_lodders(fe)))
k_fe_model(i) = alog10((element_model(pot,i)/element_lodders(pot))/(element_model(fe,i)/element_lodders(fe)))
ca_fe_model(i) = alog10((element_model(ca,i)/element_lodders(ca))/(element_model(fe,i)/element_lodders(fe)))
ir_fe_model(i) = alog10((element_model(ir,i)/element_lodders(ir))/(element_model(fe,i)/element_lodders(fe)))

sc_fe_model(i) = alog10((element_model(sc,i)/element_lodders(sc))/(element_model(fe,i)/element_lodders(fe)))
ti_fe_model(i) = alog10((element_model(ti,i)/element_lodders(ti))/(element_model(fe,i)/element_lodders(fe)))
v_fe_model(i) = alog10((element_model(va,i)/element_lodders(va))/(element_model(fe,i)/element_lodders(fe)))
cr_fe_model(i) = alog10((element_model(cr,i)/element_lodders(cr))/(element_model(fe,i)/element_lodders(fe)))
mn_fe_model(i) = alog10((element_model(mn,i)/element_lodders(mn))/(element_model(fe,i)/element_lodders(fe)))
co_fe_model(i) = alog10((element_model(co,i)/element_lodders(co))/(element_model(fe,i)/element_lodders(fe)))
ni_fe_model(i) = alog10((element_model(ni,i)/element_lodders(ni))/(element_model(fe,i)/element_lodders(fe)))
cu_fe_model(i) = alog10((element_model(cu,i)/element_lodders(cu))/(element_model(fe,i)/element_lodders(fe)))
zn_fe_model(i) = alog10((element_model(zn,i)/element_lodders(zn))/(element_model(fe,i)/element_lodders(fe)))
ga_fe_model(i) = alog10((element_model(ga,i)/element_lodders(ga))/(element_model(fe,i)/element_lodders(fe)))
ge_fe_model(i) = alog10((element_model(ger,i)/element_lodders(ger))/(element_model(fe,i)/element_lodders(fe)))
as_fe_model(i) = alog10((element_model(as,i)/element_lodders(as))/(element_model(fe,i)/element_lodders(fe)))
se_fe_model(i) = alog10((element_model(se,i)/element_lodders(se))/(element_model(fe,i)/element_lodders(fe)))
br_fe_model(i) = alog10((element_model(br,i)/element_lodders(br))/(element_model(fe,i)/element_lodders(fe)))
kr_fe_model(i) = alog10((element_model(kr,i)/element_lodders(kr))/(element_model(fe,i)/element_lodders(fe)))
rb_fe_model(i) = alog10((element_model(rb,i)/element_lodders(rb))/(element_model(fe,i)/element_lodders(fe)))
sr_fe_model(i) = alog10((element_model(sr,i)/element_lodders(sr))/(element_model(fe,i)/element_lodders(fe)))
y_fe_model (i) = alog10((element_model(yt,i)/element_lodders(yt))/(element_model(fe,i)/element_lodders(fe)))
zr_fe_model(i) = alog10((element_model(zr,i)/element_lodders(zr))/(element_model(fe,i)/element_lodders(fe)))
pb_fe_model(i) = alog10((element_model(pb,i)/element_lodders(pb))/(element_model(fe,i)/element_lodders(fe)))
os_fe_model(i) = alog10((element_model(os,i)/element_lodders(os))/(element_model(fe,i)/element_lodders(fe)))

ag_fe_model(i) = alog10((element_model(ag,i)/element_lodders(ag))/(element_model(fe,i)/element_lodders(fe)))
al_fe_model(i) = alog10((element_model(al,i)/element_lodders(al))/(element_model(fe,i)/element_lodders(fe)))
au_fe_model(i) = alog10((element_model(au,i)/element_lodders(au))/(element_model(fe,i)/element_lodders(fe)))
cr_fe_model(i) = alog10((element_model(cr,i)/element_lodders(cr))/(element_model(fe,i)/element_lodders(fe)))
dy_fe_model(i) = alog10((element_model(dy,i)/element_lodders(dy))/(element_model(fe,i)/element_lodders(fe)))
gd_fe_model(i) = alog10((element_model(gd,i)/element_lodders(gd))/(element_model(fe,i)/element_lodders(fe)))
hf_fe_model(i) = alog10((element_model(hf,i)/element_lodders(hf))/(element_model(fe,i)/element_lodders(fe)))
ho_fe_model(i) = alog10((element_model(ho,i)/element_lodders(ho))/(element_model(fe,i)/element_lodders(fe)))
lu_fe_model(i) = alog10((element_model(lu,i)/element_lodders(lu))/(element_model(fe,i)/element_lodders(fe)))
mo_fe_model(i) = alog10((element_model(mo,i)/element_lodders(mo))/(element_model(fe,i)/element_lodders(fe)))
pd_fe_model(i) = alog10((element_model(pd,i)/element_lodders(pd))/(element_model(fe,i)/element_lodders(fe)))
pt_fe_model(i) = alog10((element_model(pt,i)/element_lodders(pt))/(element_model(fe,i)/element_lodders(fe)))
u_fe_model(i) = alog10((element_model(uran,i)/element_lodders(uran))/(element_model(fe,i)/element_lodders(fe)))

fe_h_model(i) = alog10((element_model(fe,i)/element_lodders(fe))/(element_model(h,i)/element_lodders(h)))
mn_fe_model(i) = alog10((element_model(mn,i)/element_lodders(mn))/(element_model(fe,i)/element_lodders(fe)))
na_fe_model(i) = alog10((element_model(na,i)/element_lodders(na))/(element_model(fe,i)/element_lodders(fe)))
eu_fe_model(i) = alog10((element_model(eu,i)/element_lodders(eu))/(element_model(fe,i)/element_lodders(fe)))
er_fe_model(i) = alog10((element_model(er,i)/element_lodders(er))/(element_model(fe,i)/element_lodders(fe)))
ba_fe_model(i) = alog10((element_model(ba,i)/element_lodders(ba))/(element_model(fe,i)/element_lodders(fe)))
sr_fe_model(i) = alog10((element_model(sr,i)/element_lodders(sr))/(element_model(fe,i)/element_lodders(fe)))
zr_fe_model(i) = alog10((element_model(zr,i)/element_lodders(zr))/(element_model(fe,i)/element_lodders(fe)))
metallicity(i) = alog10(TOTAL(model_massfrac(2:82,i)/solarz,1))

ni_fe_model(i) = alog10((element_model(ni,i)/element_lodders(ni))/(element_model(fe,i)/element_lodders(fe)))

pb_fe_model_ss(i) = alog10((ssproc(pb,i)/element_lodders(pb))/(element_model(fe,i)/element_lodders(fe)))
pb_fe_model_r(i) = alog10((rproc(pb,i)/element_lodders(pb))/(element_model(fe,i)/element_lodders(fe)))

ba_fe_model_hs(i) = alog10((hsproc(ba,i)/element_lodders(ba))/(element_model(fe,i)/element_lodders(fe)))
ba_fe_model_r(i) = alog10((rproc(ba,i)/element_lodders(ba))/(element_model(fe,i)/element_lodders(fe)))
ba_fe_model_g(i) = alog10((gproc(ba,i)/element_lodders(ba))/(element_model(fe,i)/element_lodders(fe)))

sr_fe_model_ls(i) = alog10((lsproc(sr,i)/element_lodders(sr))/(element_model(fe,i)/element_lodders(fe)))
sr_fe_model_r(i) = alog10((rproc(sr,i)/element_lodders(sr))/(element_model(fe,i)/element_lodders(fe)))
sr_fe_model_ws(i) = alog10((wsproc(sr,i)/element_lodders(sr))/(element_model(fe,i)/element_lodders(fe)))

sr_ba_model(i) = alog10((element_model(sr,i)/element_lodders(sr))/(element_model(ba,i)/element_lodders(ba)))

eu_fe_model_s(i) = alog10((sproc(eu,i)/element_lodders(eu))/(element_model(fe,i)/element_lodders(fe)))
eu_fe_model_r(i) = alog10((rproc(eu,i)/element_lodders(eu))/(element_model(fe,i)/element_lodders(fe)))
eu_fe_model_g(i) = alog10((gproc(eu,i)/element_lodders(eu))/(element_model(fe,i)/element_lodders(fe)))

sba_h_model(i) = alog10((sproc(ba,i)/element_lodders(ba))/(element_model(h,i)/element_lodders(h)))
rba_h_model(i) = alog10((rproc(ba,i)/element_lodders(ba))/(element_model(h,i)/element_lodders(h)))
gba_h_model(i) = alog10((gproc(ba,i)/element_lodders(ba))/(element_model(h,i)/element_lodders(h)))

eu_fe_model_s(i) = alog10((sproc(eu,i)/element_lodders(eu))/(element_model(fe,i)/element_lodders(fe)))
eu_fe_model_r(i) = alog10((rproc(eu,i)/element_lodders(eu))/(element_model(fe,i)/element_lodders(fe)))
eu_fe_model_g(i) = alog10((gproc(eu,i)/element_lodders(eu))/(element_model(fe,i)/element_lodders(fe)))

mg_fe_model_ii(i) = alog10((IIproc(mg,i)/element_lodders(mg))/(element_model(fe,i)/element_lodders(fe)))
mg_fe_model_1a(i) = alog10((Iproc(mg,i)/element_lodders(mg))/(element_model(fe,i)/element_lodders(fe)))
mn_fe_model_ii(i) = alog10((IIproc(mn,i)/element_lodders(mn))/(element_model(fe,i)/element_lodders(fe)))
mn_fe_model_1a(i) = alog10((Iproc(mn,i)/element_lodders(mn))/(element_model(fe,i)/element_lodders(fe)))
;mg_fe_model_ii(i) = alog10((IIproc(mg,i)/element_lodders(mg))/(IIproc(fe,i)/element_lodders(fe)))
;mg_fe_model_1a(i) = alog10((Iproc(mg,i)/element_lodders(mg))/(Iproc(fe,i)/element_lodders(fe)))

cr_fe_model_ii(i) = alog10((IIproc(cr,i)/element_lodders(cr))/(element_model(fe,i)/element_lodders(fe)))
cr_fe_model_ia(i) = alog10((Iproc(cr,i)/element_lodders(cr))/(element_model(fe,i)/element_lodders(fe)))

o_fe_model_ii(i) = alog10((IIproc(o,i)/element_lodders(o))/(element_model(fe,i)/element_lodders(fe)))
o_fe_model_1a(i) = alog10((Iproc(o,i)/element_lodders(o))/(element_model(fe,i)/element_lodders(fe)))

mn_fe_model_ii(i) = alog10((IIproc(mn,i)/element_lodders(mn))/(element_model(fe,i)/element_lodders(fe)))
mn_fe_model_1a(i) = alog10((Iproc(mn,i)/element_lodders(mn))/(element_model(fe,i)/element_lodders(fe)))

mg_h_model(i) = alog10((element_model(mg,i)/element_lodders(mg))/(element_model(h,i)/element_lodders(h)))
mn_h_model(i) = alog10((element_model(mn,i)/element_lodders(mn))/(element_model(h,i)/element_lodders(h)))
o_h_model(i) = alog10((element_model(o,i)/element_lodders(o))/(element_model(h,i)/element_lodders(h)))
ni_h_model(i) = alog10((element_model(ni,i)/element_lodders(ni))/(element_model(h,i)/element_lodders(h)))
ba_h_model(i) = alog10((element_model(ba,i)/element_lodders(ba))/(element_model(h,i)/element_lodders(h)))
eu_h_model(i) = alog10((element_model(eu,i)/element_lodders(eu))/(element_model(h,i)/element_lodders(h)))

fe_model(i)=alog10(element_model(fe,i)/element_lodders(fe))
sr_model(i)=alog10(element_model(sr,i)/element_lodders(sr))
pb_model(i)=alog10(element_model(pb,i)/element_lodders(pb))
mg_model(i)=alog10(element_model(mg,i)/element_lodders(mg))
o_model(i)=alog10(element_model(o,i)/element_lodders(o))
eu_model(i)=alog10(element_model(eu,i)/element_lodders(eu))

fe_modelII(i)=alog10(IIproc(fe,i)/element_lodders(fe))
sr_modelII(i)=alog10(IIproc(sr,i)/element_lodders(sr))
pb_modelII(i)=alog10(IIproc(pb,i)/element_lodders(pb))
mg_modelII(i)=alog10(IIproc(mg,i)/element_lodders(mg))
o_modelII(i)=alog10(IIproc(o,i)/element_lodders(o))
eu_modelR(i)=alog10(rproc(eu,i)/element_lodders(eu))

sr_modelW(i)=alog10(wsproc(sr,i)/element_lodders(sr))

ba_mg_model[i]=alog10((element_model(ba,i)/element_lodders(ba))/(element_model(mg,i)/element_lodders(mg)))
eu_mg_model[i]=alog10((element_model(eu,i)/element_lodders(eu))/(element_model(mg,i)/element_lodders(mg)))

feP(i)=(element_model(fe,i)/element_lodders(fe))

ENDFOR

;stop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;COMPUTE FIT PARAMETERS IN TERMS OF Z
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
varswitch=0
IF varswitch eq 1 then begin
;z := xi , z(zero)=1 - solar
xi=z
z=10^(metallicity[*])

near = Min(Abs(metallicity - 0.5), p05)  ;p05 is the index where feh[*] is closest to +0.5
near = Min(Abs(metallicity + 0.0), z0)
near = Min(Abs(metallicity + 0.5), m05)
near = Min(Abs(metallicity + 1.0), m1)
near = Min(Abs(metallicity + 2.0), m2)
near = Min(Abs(metallicity + 3.0), m3)
near = Min(Abs(metallicity + 4.0), m4)

near = Min(Abs(fe_model - 0.5), fp05)  ;p05 is the index where feh[*] is closest to +0.5
near = Min(Abs(fe_model + 0.0), fz0)
near = Min(Abs(fe_model + 0.5), fm05)
near = Min(Abs(fe_model + 1.0), fm1)
near = Min(Abs(fe_model + 2.0), fm2)
near = Min(Abs(fe_model + 3.0), fm3)
near = Min(Abs(fe_model + 4.0), fm4)


;x=dln(s)/dln(z)=(z/s)*(ds/dz)

;l=(z)^ls(SL)
;s=(z)^hs(V)
;r=(z)^rp(M)
;w=(z)^wp(G)

;dl=deriv(alog(l))/deriv(alog(xi))+deriv(alog(xi))/deriv(alog(z))
;ds=deriv(alog(s))/deriv(alog(xi))+deriv(alog(xi))/deriv(alog(z))
;dr=deriv(alog(r))/deriv(alog(xi))+deriv(alog(xi))/deriv(alog(z))
;dw=deriv(alog(w))/deriv(alog(xi))+deriv(alog(xi))/deriv(alog(z))

l=(xi)^ls(SL)
s=(xi)^hs(V)
r=(xi)^rp(M)
w=(xi)^wp(G)
II=isoii(61,*)
Ia=iso1A(61,*)
g=iso_gam(176,*)

dl=deriv(alog(l))/deriv(alog(z))
ds=deriv(alog(s))/deriv(alog(z))
dr=deriv(alog(r))/deriv(alog(z))
dw=deriv(alog(w))/deriv(alog(z))
dII=deriv(alog(II))/deriv(alog(z))
dIa=deriv(alog(Ia))/deriv(alog(z))
dg=deriv(alog(g))/deriv(alog(z))

ind=z0
find=fz0

print,'ls=',dl[ind]
print,'hs=',ds[ind]
print,'r=',dr[ind]
print,'w=',dw[ind]
print,'II=',dII[ind]
print,'Ia=',dIa[ind]
print,'g=',dg[ind]

dl_fe=deriv(alog(l))/deriv(alog(feP))
ds_fe=deriv(alog(s))/deriv(alog(feP))
dr_fe=deriv(alog(r))/deriv(alog(feP))
dw_fe=deriv(alog(w))/deriv(alog(feP))
dII_fe=deriv(alog(II))/deriv(alog(feP))
dIa_fe=deriv(alog(Ia))/deriv(alog(feP))
dg_fe=deriv(alog(g))/deriv(alog(feP))

print,'ls_fe=',dl_fe[find]
print,'hs_fe=',ds_fe[find]
print,'r_fe=',dr_fe[find]
print,'w_fe=',dw_fe[find]
print,'II_fe=',dII_fe[find]
print,'Ia_fe=',dIa_fe[find]
print,'g_fe=',dg_fe[find]


;stop

openpsfl,'aaaa2.ps'

c = LONARR(7)
c(6)=rgb(150,0,0)   ;ORANGE	
c(5)=rgb(1.0,0,1.0) ;PINK
c(4)=rgb(100,50,50) ;BROWN
c(3)=rgb(0,0,1.0)   ;BLUE 
c(2)=rgb(1.0,0,0)   ;RED 
c(1)=rgb(0,1.0,1.0) ;TEAL
c(0)=rgb(0,1.0,0)   ;GREEN 

th=4

; yrange=[0,2.0] ;for just: hs, ls, w, r
yrange=[0,5.]
xrange=[-4,0.5]
PLOT,[0,1],/NODATA, $
      XRANGE=xrange , $
      YRANGE=yrange, $
      XS=1, $
      YS=1, $
      XMARGIN=[1.0,0], $
      YMARGIN=[2.5,0], $
      ;ytitle='!Z(00F1)', $
      ytitle='exponent', $
      xtitle='[Z] (solid line), [Fe] (dashed line)', $
      charsize=1.5

oplot,metallicity,dl,color=c[0],thick=th
oplot,metallicity,ds,color=c[1],thick=th
oplot,metallicity,dr,color=c[2],thick=th
;oplot,metallicity,dw,color=c[3],thick=th
oplot,metallicity,dII,color=c[3],thick=th
oplot,metallicity,dIa,color=c[4],thick=th
oplot,metallicity,dg,color=c[5],thick=th


oplot,fe_model,dl_fe,color=c[0],thick=th,line=2
oplot,fe_model,ds_fe,color=c[1],thick=th,line=2
oplot,fe_model,dr_fe,color=c[2],thick=th,line=2
;oplot,fe_model,dw_fe,color=c[3],thick=th,line=2
oplot,metallicity,dII_fe,color=c[3],thick=th,line=2
oplot,metallicity,dIa_fe,color=c[4],thick=th,line=2
oplot,metallicity,dg_fe,color=c[5],thick=th,line=2


vline,0,line=1

name=['ls, ws','hs','r','Massive','SNe Ia','!7c!3']

;FOR JUST: ls, hs, w, r
;legende,/bottom,/left,/DRAWFRAME,FRAMETOP=42,FRAMEBOTTOM=7.7,FRAMERIGHT=6,FRAMELEFT=7

;FOR ALL
legende,/top,/left,/DRAWFRAME,FRAMETOP=12,FRAMEBOTTOM=70,FRAMERIGHT=8,FRAMELEFT=7


for i=0,4 do begin
   legende_line,0,name[i],color=c[i],thick=4
endfor

closeps


print,'[VARSWITCH]: DO NOT CONTINUE. SET varswitch EQ 0 @ LINE 1232'
STOP
ENDIF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;END COMPUTE FIT PARAMETERS IN TERMS OF Z
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;binwrite,array=sr_modelW,filename='~/kepler/zdep/GCH_data/srW'

;binwrite,array=fe_model,filename='~/kepler/zdep/GCH_data/fe'
;binwrite,array=sr_model,filename='~/kepler/zdep/GCH_data/sr'
;binwrite,array=pb_model,filename='~/kepler/zdep/GCH_data/pb'
;binwrite,array=mg_model,filename='~/kepler/zdep/GCH_data/mg'
;binwrite,array=o_model,filename='~/kepler/zdep/GCH_data/o'
;binwrite,array=eu_modelR,filename='~/kepler/zdep/GCH_data/euR'

;binwrite,array=fe_modelII,filename='~/kepler/zdep/GCH_data/feII'
;binwrite,array=sr_modelII,filename='~/kepler/zdep/GCH_data/srII'
;binwrite,array=pb_modelII,filename='~/kepler/zdep/GCH_data/pbII'
;binwrite,array=mg_modelII,filename='~/kepler/zdep/GCH_data/mgII'
;binwrite,array=o_modelII,filename='~/kepler/zdep/GCH_data/oII'
;binwrite,array=eu_modelII,filename='~/kepler/zdep/GCH_data/euII'
;stop

;MORE PLOTS
moreplots = 0
if moreplots eq 1 then begin
openpsfl,'~/gce/plots/xi_vs_z_'+date+'.eps'
xr=[-5.0,0.5]
yr=[-5.0,0.5]
!X.MARGIN=[1,0] 
!Y.MARGIN=[2.5,0] 
PLOT,x,y,/noerase,/nodata,xs=1,ys=1,XR=xr,YR=yr,charsize=1.5,$
      XTITLE='log(!7n!3)', $
      YTITLE='[Z]'
oplot,logz,metallicity
unity=dindgen(100)/15-5
oplot,unity,unity,line=2
vline,0,line=2
hline,0,line=2
closeps

stop

openpsfl,'~/gce/plots/frebel_compare/'+ggg+'.ps'
plot,feh,mgfe,ytitle='['+ggg+'/Fe]',psym=7
oploterr,feh,mgfe,mgfe_sig
oplot,fe_h_model,cr_fe_model
oplot,fe_h_model,cr_fe_model_ii,color=rgb(1.,0,0)
oplot,fe_h_model,cr_fe_model_ia,color=rgb(0,1.,0)
closeps
stop
endif


;MULTI PLOTS
if ggg eq 'multi' then begin
!p.multi=[0,4,4]
stop
endif

;END MULTI PLOTS

;PROCESS PLOTS
process=0
if process eq 1 then begin
hsz=iso_hs(180,zero)
rz=iso_r(204,zero)
gz=iso_gam(176,zero)
Iaz=iso1A(61,zero)
wsz=iso_ws(94,zero)
iiz=isoii(61,zero)
lsz=iso_ls(94,zero)
ssz=iso_ss(pb204,zero)

xr=[-6.0,1.0]
yr=[-6.5,2.0]

if psmake eq 1 then begin
  openpsfl,'~/Desktop/gce/plots/processplot_'+date+'.ps'
endif

el_lod_fe=2.9174709e-06
el_lod_o=0.00041348190

!X.MARGIN=[1,0]
!Y.MARGIN=[2,0]

th=4
plot,xr,yr,/nodata,XTITLE='[Fe/H]',YTITLE='[Abundance/O]',/YS,/XS
oplot,fe_h_model,alog10((iso_ls(94,*)/lsz[0])/(element_model(o,*)/el_lod_o)),color=rgb(0,0,1.D0),thick=th,line=1    ;LIGHT BLUE (LS)
oplot,fe_h_model,alog10((iso_ss(pb204,*)/ssz[0])/(element_model(o,*)/el_lod_o)),color=rgb(0,1d0,1.d0),thick=2*th,line=0 ;BLUE (STRONG)
oplot,fe_h_model,alog10((iso_hs(180,*)/hsz[0])/(element_model(o,*)/el_lod_o)),color=rgb(1.D0,0,0),thick=th,line=3   ;RED (HS)
oplot,fe_h_model,alog10((iso_r(204,*)/rz[0])/(element_model(o,*)/el_lod_o)),color=rgb(0,0,0),thick=2*th,line=1  ;BLACK (RPROC)
oplot,fe_h_model,alog10((iso_gam(176,*)/gz[0])/(element_model(o,*)/el_lod_o)),color=rgb(0,1.D0,0),thick=th,line=4   ;GREEN (GAMMA)
oplot,fe_h_model,alog10((iso1A(61,*)/Iaz[0])/(element_model(o,*)/el_lod_o)),color=rgb(1.d0,1.D0,0),thick=th,line=5  ;YELLOW (Ia)
oplot,fe_h_model,alog10((iso_ws(94,*)/wsz[0])/(element_model(o,*)/el_lod_o)),color=rgb(1.D0,0,1.D0),thick=2*th,line=2 ;PINK (WS) 
oplot,fe_h_model,alog10((isoii(61,*)/iiz[0])/(element_model(o,*)/el_lod_o)),color=rgb(100,50,50),thick=th,line=0    ;BROWN (II) 
onset = where(ABS(alog10(iso1A(61,*)/Iaz[0])+2) eq min(ABS(alog10(iso1A(61,*)/Iaz[0])+2)))
print,fe_h_model(onset)
;1.40948e-5
;vline,0,line=2
;hline,0,line=2

c = LONARR(8)
c(7)=rgb(0,0,1d0) ;LIGHT BLUE (LS)
c(6)=rgb(0,1d0,1.D0)    ;BLUE (STRONG)
c(5)=rgb(100,50,50)   ;BROWN (II) 
c(4)=rgb(0,0,0)       ;BLUE (RPROC)
c(3)=rgb(1.D0,0,1.D0) ;PINK (WS)
c(2)=rgb(1.D0,0,0)    ;RED (HS)
c(1)=rgb(0,1.D0,0)    ;GREEN (GAMMA)
c(0)=rgb(1.D0,1.D0,0) ;YELLOW (Ia)

ll = INTARR(8)
ll(7)=1   ;line (LS) 
ll(6)=0   ;line (SS) 
ll(5)=0   ;line (II) 
ll(4)=1   ;dot (RPROC)
ll(3)=2   ;dash (WS)
ll(2)=3   ;dot-dash (HS)
ll(1)=4   ;dot-dot-dot-dash (GAMMA)
ll(0)=5   ;long-dash (Ia)
lab_proc=['Type Ia SNe','Gamma','Heavy S','Weak S','R-Process','Massive','Strong S','Light S']

thi = INTARR(8)
thi(7)=th
thi(6)=2*th
thi(5)=th
thi(4)=2*th
thi(3)=2*th
thi(2)=th
thi(1)=th
thi(0)=th

LEGENDE,/BOTTOM,/RIGHT,CHARSIZE=1.2,/drawframe,FRAMETOP=120,FRAMEBOTTOM=6,FRAMERIGHT=11,FRAMELEFT=7,xpos=0.98
FOR i=0,N_ELEMENTS(lab_proc)-1 DO BEGIN
;   LEGENDE_COL,c[i],lab_proc[i],/noframe
    LEGENDE_LINE,ll[i],lab_proc[i],color=c[i],thick=thi[i]
ENDFOR

if psmake eq 1 then begin
  closeps
endif

print,'[GCE] line 1159: Do not continue, yr is now wrong' 
stop 
endif
;re-call gce_yr function, since yrange was changed for process plot
yr=gce_yr(ggg)
; END PROCESS PLOT


;EXPORT TO COLOR FOLDER
color=0
if color eq 1 then begin
data_el = element_model
openw,lun,'~/Desktop/gce/pplot/colordat/matrix'+date+'.dat'
for i=0,82 do begin
   printf,lun,data_el(i,*)
endfor

close,lun
free_lun,lun

openw,lun,'~/Desktop/gce/pplot/colordat/fe_h_'+date+'.dat'
for i = 0,z_max-1 do printf,lun,fe_h_model(i)

close,lun
free_lun,lun

openw,lun,'~/Desktop/gce/pplot/colordat/solar_abun_'+date+'.dat'
printf,lun,element_lodders

close,lun
free_lun,lun

print,'[GCE] EXPORTED MODEL TO PPLOT FOR ISOBOX PLOTS.'
stop
endif
;END EXPORT

;;;;;;;;;;;;;;;;;;;;;;;;;;
;CHOOSE VARIABLE FOR PLOTS

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
ENDIF

IF ggg EQ 'n' THEN BEGIN
y_model = n_fe_model
y_title = '[N/Fe]'
ENDIF

IF ggg EQ 'n' THEN BEGIN
y_model = n_fe_model
y_title = '[N/Fe]'
ENDIF

IF ggg EQ 'ca' THEN BEGIN
y_model = ca_fe_model
y_title = '[Ca/Fe]'
ENDIF

IF ggg EQ 'ge' THEN BEGIN
y_model = ge_fe_model
y_title = '[Ge/Fe]'
ENDIF

IF ggg EQ 'al' THEN BEGIN
y_model = al_fe_model
y_title = '[Al/Fe]'
ENDIF

IF ggg EQ 'na' THEN BEGIN
y_model = na_fe_model
y_title = '[Na/Fe]'
ENDIF

;
;;;;;;;;;;;;;;;;;;;;;;;;;;

chi_squared = dblarr(nx)
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

;colorind = (ST+1)/n_elements(space5)
;openpsfl,'plots/pbf_bestlooking.ps'
;plot,x_model,y_model,yrange=yr,xrange=[-5.0,0.5]
;oplot,feh,mgfe,psym=2
;closeps
;stop

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


;if ggg eq 'ba' then print,where(param[1,*] eq min(param[1,*])),param[*,where(param[1,*] eq min(param[1,*]))]
;if ggg eq 'eu' then print,where(param[1,*] eq min(param[1,*])),param[*,where(param[1,*] eq min(param[1,*]))]
;if ggg eq 'mg' or ggg eq 'o' then print,where(param[3,*] eq min(param[3,*])),param[*,where(param[3,*] eq min(param[3,*]))]
;;if ggg eq 5 then print,where(param[1,*] eq min(param[1,*])),param[*,where(param[1,*] eq min(param[1,*]))]
;;if ggg eq 6 then print,where(param[1,*] eq min(param[1,*])),param[*,where(param[1,*] eq min(param[1,*]))]
;if ggg eq 'sr' then print,where(param[2,*] eq min(param[2,*])),param[*,where(param[2,*] eq min(param[2,*]))]
;if ggg eq 'pb' then print,where(param[3,*] eq min(param[3,*])),param[*,where(param[3,*] eq min(param[3,*]))]

;stop

if psmake eq 1 then begin
openpsfl,'~/Desktop/gce/plots/'+tit+'fe_'+date+'.ps'
endif

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

xbin=[x[0]-xres,x]+0.5D0*xres
ybin=[y[0]-yres,y]+0.5D0*yres

xoverlap=0.05D0
yoverlap=0.05D0

ftot/=ftmax
;ftot/=fmax
ftot=ALOG10(ftot>10.D0^mrange[1])/(-mrange[1])+1
colorfunc,colorscheme,ftot,c


x=xr[0]+(dindgen((xr[1]-xr[0])/xres)+0.5D0)*xres
y=yr[0]+(dindgen((yr[1]-yr[0])/yres)+0.5D0)*yres

!P.BACKGROUND = 'FFFFFF'XL   ;Set background color to white
!P.COLOR = 'OOOOOO'XL  

;!P.BACKGROUND = rgb(1.)     ;Set background color to white
;!P.COLOR = rgb(0)           ;Set plotline to black
; erase
;;;

PLOT,x,y,/nodata,xs=5,ys=5,XR=xr,YR=yr,$
      XMARGIN=xmargin, $
      YMARGIN=ymargin

skip=0
if ggg eq 'pb' then skip=1

if skip eq 0 then begin
;;;
ny=N_ELEMENTS(y)
nx=N_ELEMENTS(x)

;;;
;plot both colors
;c=[[[c1[*,*,0]]],[[c[*,*,1]]],[[c[*,*,2]]]]

;plot only red
c=[[[c[*,*,0]]],[[c[*,*,1]]],[[c[*,*,2]]]]
dx=!P.clip[2]-!P.clip[0]+1
dy=!P.clip[3]-!P.clip[1]+1

TVLCT,r_sav,g_sav,b_sav,/GET
TVLCT, INDGEN(256), INDGEN(256), INDGEN(256)
IF !D.FLAGS mod 2 EQ 0 THEN BEGIN
    c=congrid(c,dx,dy,3,/center,/cubic)
ENDIF
tv,c,!P.clip[0],!P.clip[1],true=3,xsize=dx,ysize=dy
TVLCT,r_sav,g_sav,b_sav

;;;
endif

;PLOTS,feh,mgfe,color=rgb(1.D0,0,0),psym=1


;;;;
;MG;
;;;;
if ggg eq 'mg' then begin

if freb eq 1 then begin
  OPLOT,feh[0:852],mgfe[0:852],color=rgb(0,0,1.d0),psym=7
;psym=sym(12)
  OPLOT,feh[853:1576],mgfe[853:1576],color=rgb(1d0,0,0),psym=2
;  OPLOT,feh[853:1029],mgfe[853:1029],color=rgb(1d0,1d0,0),psym=sym(11)
LEGENDE,/TOP,/RIGHT,/DRAWFRAME,FRAMETOP=8,FRAMEBOTTOM=18,FRAMERIGHT=9,FRAMELEFT=5,CHARSIZE=1.2
legende_symbol,7,'Frebel',color=rgb(0,0,1.d0)
legende_symbol,2,'Soubiran',color=rgb(1d0,0,0)
ENDIF
;PLOT,feh[0:852],mgfe[0:852],color=rgb(0,0,1.d0),psym=sym(12)
;OPLOT,feh[853:1576],mgfe[853:1576],color=rgb(1d0,1d0,0),psym=sym(11)
;;  OPLOT,feh[853:1029],mgfe[853:1029],color=rgb(1d0,1d0,0),psym=sym(11)
;LEGENDE,/TOP,/RIGHT,/DRAWFRAME,FRAMETOP=8,FRAMEBOTTOM=18,FRAMERIGHT=9,FRAMELEFT=5,CHARSIZE=1.2
;legende_symbol,7,'Frebel',color=rgb(0,0,1.d0)
;legende_symbol,1,'Soubiran',color=rgb(1d0,1d0,0)

if freb eq 0 then begin
;Zhang
OPLOT,feh[954:985],mgfe[954:985],color=rgb(1.D0,0,0),psym=sym(1)
;Mashonkina
OPLOT,feh[939:953],mgfe[939:953],color=rgb(139,69,19),psym=sym(2)
;Cayrel
OPLOT,feh[179:213],mgfe[179:213],color=rgb(0,0,1.D0),psym=sym(8)
;Fulbright
OPLOT,feh[0:178],mgfe[0:178],color=rgb(0,100,0),psym=sym(4)
;Barklem
OPLOT,feh[986:1227],mgfe[986:1227],color=rgb(0,1.D0,1.D0),psym=sym(10)
;Soubiran
OPLOT,feh[214:938],mgfe[214:938],color=rgb(1.D0,1.D0,0),psym=sym(16)
endif

;oplot,fe_h_model,mg_fe_model_1a,color=rgb(0,1d0,1d0)
;oplot,fe_h_model,mg_fe_model_ii,color=rgb(100,50,50)
;c=isocolors(3,saturation=.05)
;c[0]=rgb(0,1d0,0)
;c[1]=rgb(100,50,50)
;c[2]=rgb(0,1d0,1d0)
;lab_proc=['SNIa+Massive','Massive','SNIa']
;;LEGENDE,/LEFT,/TOP
;LEGENDE,/RIGHT,/BOTTOM
;FOR i=0,N_ELEMENTS(lab_proc)-1 DO BEGIN
;   LEGENDE_COL,c[i],lab_proc[i],/noframe
;ENDFOR

if pcompare eq 1 then oplot,pran[0,*],pran[1,*],thick=2
endif

;;;;
;PB;
;;;;
if ggg eq 'pb' then begin

if freb eq 1 then begin
OPLOT,feh,mgfe,color=rgb(1.d0,0,0),psym=sym(12)

oplot,fe_h_model,pb_fe_model_ss,color=rgb(0,0,1d0),THICK=4,line=5
oplot,fe_h_model,pb_fe_model_r,color=rgb(100,50,50),THICK=4,line=3

c=isocolors(3,saturation=.05)
c[0]=rgb(0,1d0,0)
c[1]=rgb(0,0,1d0)
c[2]=rgb(100,50,50)

ll=intarr(3)
ll[0]=0
ll[1]=5
ll[2]=3

lab_proc=['hs+R','hs','R-Process']
LEGENDE,/RIGHT,/TOP,/DRAWFRAME,FRAMETOP=9.5,FRAMEBOTTOM=39,FRAMERIGHT=9,FRAMELEFT=6.2,CHARSIZE=1.2,XPOS=1.0,YPOS=1.0
;LEGENDE,/RIGHT,/BOTTOM
FOR i=0,N_ELEMENTS(lab_proc)-1 DO BEGIN
   ;LEGENDE_COL,c[i],lab_proc[i],/noframe
   LEGENDE_LINE,ll[i],lab_proc[i],color=c[i],thick=4
 ENDFOR
LEGENDE,/LEFT,/TOP,/DRAWFRAME,FRAMETOP=10,FRAMEBOTTOM=8,FRAMERIGHT=6,FRAMELEFT=4.8,CHARSIZE=1.2,XPOS=-0.02,YPOS=1.0
legende_symbol,sym(12),'Frebel',color=rgb(1.0,0,0)

;now do legend for data points

endif

endif

;;;;
;SR;
;;;;
if ggg eq 'sr' then begin

if freb eq 1 then begin
OPLOT,feh[0:465],mgfe[0:465],color=rgb(0,0,1.d0),psym=sym(12)
OPLOT,feh[466:504],mgfe[466:504],color=rgb(1d0,0,0),psym=sym(12)
OPLOT,feh[505:518],mgfe[505:518],psym=2

oplot,fe_h_model,sr_fe_model_ls,color=rgb(0,1.d0,1.d0),THICK=4,line=3
oplot,fe_h_model,sr_fe_model_ws,color=rgb(100,50,50),THICK=4,line=5

c=isocolors(3,saturation=.05)
c[0]=rgb(0,1d0,0)
c[1]=rgb(0,1.d0,1.d0)
c[2]=rgb(100,50,50)

ll=intarr(3)
ll[0]=0
ll[1]=3
ll[2]=5

lab_proc=['ls+Weak','ls','Weak S']
;LEGENDE,/LEFT,/TOP
LEGENDE,/RIGHT,/BOTTOM,/DRAWFRAME,FRAMETOP=40,FRAMEBOTTOM=7,FRAMERIGHT=8,FRAMELEFT=7,XPOS=1.05,YPOS=0.02
FOR i=0,N_ELEMENTS(lab_proc)-1 DO BEGIN
   ;LEGENDE_COL,c[i],lab_proc[i],/noframe
    LEGENDE_LINE,ll[i],lab_proc[i],color=c[i],thick=4
ENDFOR

LEGENDE,/TOP,/RIGHT,/DRAWFRAME,FRAMETOP=8,FRAMEBOTTOM=36,FRAMERIGHT=10,FRAMELEFT=5,CHARSIZE=1.2,XPOS=1.05,YPOS=1.02
legende_symbol,sym(12),'Frebel',color=rgb(0,0,1.d0)
legende_symbol,sym(12),'Mashonkina',color=rgb(1d0,0,0)
legende_symbol,2,'Jehin'

endif

endif

;;;;
;BA;
;;;;
if ggg eq 'ba' then begin

if freb eq 1 then begin
OPLOT,feh,mgfe,color=rgb(0,0,1.0),psym=sym(12)
LEGENDE,/LEFT,/TOP,/DRAWFRAME,FRAMETOP=10,FRAMEBOTTOM=8,FRAMERIGHT=6,FRAMELEFT=4.8,CHARSIZE=1.2,XPOS=-0.1,YPOS=1.03
legende_symbol,sym(12),'Frebel',color=rgb(0,0,1.d0)
endif

oplot,fe_h_model,ba_fe_model_r,color=rgb(1d0,1d0,0),THICK=4,line=5     ;1
oplot,fe_h_model,ba_fe_model_hs,color=rgb(1d0,0,0),THICK=4,line=3      ;2
oplot,fe_h_model,ba_fe_model_g,color=rgb(0,0,0),THICK=4,line=4   ;3

c=isocolors(4,saturation=.05)
c[0]=rgb(0,1d0,0)
c[1]=rgb(1d0,1d0,0)
c[2]=rgb(1d0,0,0)
c[3]=rgb(0,0,0)

ll=intarr(4)
ll[0]=0
ll[1]=5
ll[2]=3
ll[3]=4

thi=intarr(4)
thi[0]=4
thi[1]=4
thi[2]=4
thi[3]=4

lab_proc=['R+hs+G','R-Process','hs','Gamma']
;LEGENDE,/LEFT,/TOP
LEGENDE,/RIGHT,/TOP,/DRAWFRAME,FRAMETOP=9.5,FRAMEBOTTOM=46,FRAMERIGHT=9,FRAMELEFT=6.2,CHARSIZE=1.2,XPOS=1.09,YPOS=1.03
FOR i=0,N_ELEMENTS(lab_proc)-1 DO BEGIN
   ;LEGENDE_COL,c[i],lab_proc[i],line=ll[i],/noframe
    LEGENDE_LINE,ll[i],lab_proc[i],color=c[i],thick=thi[i]
ENDFOR

if freb eq 0 then begin
;Zhang
OPLOT,feh[0:31],mgfe[0:31],color=rgb(1.D0,0,0),psym=sym(1)
;Mashonkina
OPLOT,feh[32:46],mgfe[32:46],color=rgb(139,69,19),psym=sym(2)
;Andrievsky
OPLOT,feh[47:87],mgfe[47:87],color=rgb(0,0,1.D0),psym=sym(2)
;Fulbright
OPLOT,feh[88:266],mgfe[88:266],color=rgb(255,255,0),psym=sym(9)
;Barklem
OPLOT,feh[267:485],mgfe[267:485],color=rgb(0,100,0),psym=sym(10)
endif
endif
;;;

;;;;
;O;
;;;;
if ggg eq 'o' then begin

if freb eq 1 then OPLOT,feh,mgfe,color=rgb(0,0,1d0),psym=sym(12)

if freb eq 0 then begin
;Zhang
OPLOT,feh[437:467],mgfe[437:467],color=rgb(1.D0,0,0),psym=sym(1)
;Soubiran
OPLOT,feh[0:414],mgfe[0:414],color=rgb(139,69,19),psym=sym(7)
;Cayrel
OPLOT,feh[415:436],mgfe[415:436],color=rgb(0,0,1.D0),psym=sym(8)
endif

endif
;;;

;;;;
;EU;
;;;;
if ggg eq 'eu' then begin

if freb eq 1 then begin
OPLOT,feh,mgfe,color=rgb(0,0,1.0),psym=sym(12)
LEGENDE,/RIGHT,/TOP,/DRAWFRAME,FRAMETOP=10,FRAMEBOTTOM=8,FRAMERIGHT=6,FRAMELEFT=4.8,CHARSIZE=1.2,XPOS=1.25,YPOS=1.02
legende_symbol,sym(12),'Frebel',color=rgb(0,0,1.d0)
endif

;oplot,fe_h_model,eu_fe_model_s,color=rgb(0,1d0,1d0),THICK=4
;oplot,fe_h_model,eu_fe_model_r,color=rgb(100,50,50),THICK=4
;c=isocolors(3,saturation=.05)
;c[0]=rgb(0,1d0,0)
;c[1]=rgb(100,50,50)
;c[2]=rgb(0,1d0,1d0)
;lab_proc=['R+S','R-Process','S-Process']
;;LEGENDE,/LEFT,/TOP
;LEGENDE,/RIGHT,/BOTTOM
;FOR i=0,N_ELEMENTS(lab_proc)-1 DO BEGIN
;   LEGENDE_COL,c[i],lab_proc[i]
;ENDFOR

if freb eq 0 then begin
;Mashonkina
OPLOT,feh[0:13],mgfe[0:13],color=rgb(139,69,19),psym=sym(2)
;Del Peloso
;OPLOT,feh[14:33],mgfe[14:33],color=rgb(0,0,1.D0),psym=sym(8)
;Koch
OPLOT,feh[14:87],mgfe[14:87],color=rgb(1.D0,0,0),psym=sym(1)
;Gonzalez
OPLOT,feh[88:115],mgfe[88:115],color=rgb(0,1.D0,1.D0),psym=sym(10)
;Fulbright
OPLOT,feh[116:200],mgfe[116:200],color=rgb(255,255,0),psym=sym(9)
;Barklem
OPLOT,feh[201:268],mgfe[201:268],color=rgb(0,100,0),psym=sym(14)
endif

endif
;;;

;;;;
;Mn;
;;;;
if ggg eq 'mn' then begin
;Zhang
OPLOT,feh[0:29],mgfe[0:29],color=rgb(1.D0,0,0),psym=sym(1)
;Cayrel
OPLOT,feh[30:64],mgfe[30:64],color=rgb(0,0,1.D0),psym=sym(8)
;Barklem
OPLOT,feh[65:301],mgfe[65:301],color=rgb(0,100,0),psym=sym(14)
endif
;;;

full=where(fe_h_model ge -4.0)
hypen=where(fe_h_model lt -4.0)
oplot,fe_h_model[full],y_model[full],color=rgb(0,255,0),THICK=8
oplot,fe_h_model[hypen],y_model[hypen],color=rgb(0,255,0),THICK=8,line=2

ii=WHERE(w1 GT 0.1)
;ii=WHERE(w1 GT 0.01)

c=rgb(1.,0,1.)

;oplot,x[ii],av[ii],line=0,color=c
;oplot,x[ii],av[ii]+sig[ii],line=2,color=c
;oplot,x[ii],av[ii]-sig[ii],line=2,color=c

jj=(max(ii)-ii)+min(ii)

xx=[x[ii],x[jj]]
yy=[av[ii]+sig[ii],av[jj]-sig[jj]]

;uncomment for hatches
;POLYFILL,xx,yy,/line_fill,COLOR=c,orientation=30

;the following plots for too low metallicities for frebel data set
if freb eq 0 then begin
plots,x[ii],av[ii]+sig[ii],color=c,thick=4
plots,x[ii],av[ii]-sig[ii],color=c,thick=4
plots,x[ii],av[ii],color=c,thick=8
endif

if ggg eq 'pb' then skip=1
if ggg ne 'pb' then skip=0

if freb eq 1 and skip eq 0 then begin
xxx=x[ii]
aaa=av[ii]
sss=sig[ii]
kk=where(xxx gt -4.0)
plots,xxx[kk],aaa[kk]+sss[kk],color=c,thick=10
plots,xxx[kk],aaa[kk]-sss[kk],color=c,thick=10
plots,xxx[kk],aaa[kk],color=c,thick=10
endif

barr=dblarr(2,n_elements(aaa[kk]))
barr[0,*]=xxx[kk]
barr[1,*]=aaa[kk]

plot,barr[0,*],barr[1,*],psym=2
stop

binwrite,array=barr,filename='~/kepler/zdep/frebel_av/'+ggg+'fe'

stop

PLOT,x,y,/noerase,/nodata,xs=1,ys=1,XR=xr,YR=yr,$
      XMARGIN=xmargin, $
      YMARGIN=ymargin, $
      XTITLE='[Fe/H]', $
      YTITLE=y_title

if skip eq 0 then begin
bartitle='Log(Relative Contribution)'
;; stop
ch=CONVERT_COORD(!D.X_CH_SIZE,!D.Y_CH_SIZE,/DEVICE,/TO_NORMAL)

colorbar, $
   1-(colbarwidth+colbarmargin)*ch[0], $
   (ymargin[0])*ch[1], $
   colbarwidth*ch[0], $
   1-(ymargin[0]+ymargin[1])*ch[1],$
   mrange[[1,0]], $
   colorscheme, $
   anno=bartitle
endif

if psmake eq 1 then begin
closeps
endif

stop
openpsfl,'~/Desktop/gce/plots/mgfe_10_18_2011.ps'
filename='~/Desktop/gce/plots/mgfe_10_18_2011_' + cstr + '.jpeg'
write_jpeg,filename,tvrd()

stop
END

;-----------------------------------------------------------------------
;ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
;-----------------------------------------------------------------------

PRO test26,x,o1,o2,o3, NAME=name

NAME='gray gamma color function'

IF N_ELEMENTS(x) EQ 0 THEN RETURN


col_min=0
col_max=1
y=(x>col_min)<col_max

gamma=2.D0
y^=gamma

MAP=[$
    [0,1,1,1,1], $
    [1,0,0,0,1]]


;backup
;MAP=[$
;    [0,1,1,1,1], $
;    [1,0,0,0,1]]

MAP=DOUBLE(map)
MAP[1:3,*]*=255

colormap, y, o1, o2, o3, $
  MAP_MODE=0, $
  MAP_DATA=map

END

;-----------------------------------------------------------------------
;ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
;-----------------------------------------------------------------------

function CONGRID, arr, x, y, z, $
    CENTER=center, $
    CUBIC = cubicIn, $
    INTERP=interp, $
    MINUS_ONE=minus_one

    COMPILE_OPT idl2
    ON_ERROR, 2     ;Return to caller if error

    ndim = SIZE(arr, /N_DIMENSIONS)
    dims = SIZE(arr, /DIMENSIONS)

    if ((ndim lt 1) or (ndim gt 3)) then $
      Message, 'Array must have 1, 2, or 3 dimensions.'

    ;;  Supply defaults = no interpolate, and no minus_one.
    int = KEYWORD_SET(interp)
    m1 = KEYWORD_SET(minus_one)
    cubic = (N_ELEMENTS(cubicIn) gt 0) ? cubicIn : 0
    if (cubic ne 0) then int = 1    ;Cubic implies interpolate
    offset = KEYWORD_SET(center) ? 0.5 : 0.0

    ; Construct new interpolate coordinates.
    ; Skip this for 2D nearest-neighbor since we use POLY_2D instead.
    if ((ndim ne 2) or ((ndim eq 2) and int)) then begin
        ; Note that we need to use "offset" twice: Once to shift the new
        ; coordinates to the midpoint, and again to shift the location of
        ; the original coordinates to their midpoint.
        switch ndim of  ; Fall through for ndim>1.
            3: srz = float(dims[2]-m1)/(z-m1)*(findgen(z) + offset) - offset
            2: sry = float(dims[1]-m1)/(y-m1)*(findgen(y) + offset) - offset
            1: srx = float(dims[0]-m1)/(x-m1)*(findgen(x) + offset) - offset
        endswitch
    endif

    case ndim of
        1: begin                ; *** ONE DIMENSIONAL ARRAY
            arr_r = (int) ? INTERPOLATE(arr, srx, CUBIC = cubic) : $
                arr[ROUND(srx)]
           end
        2: begin                ; *** TWO DIMENSIONAL ARRAY
            if (int) then begin  ; bilinear or cubic
                arr_r = INTERPOLATE(arr, srx, sry, /GRID, CUBIC=cubic)
            endif else begin  ; nearest neighbor
                ; Note: For expansion, divide by (x-1) so that CONGRID
                ; will agree with REBIN.
                expand = (x gt dims[0])
                xm1 = (m1 or expand) ? x-1 : x
                arr_r = POLY_2D(arr, $
                    [[0,0],[(dims[0]-m1)/float(xm1),0]], $ ;Use poly_2d
                    [[0,(dims[1]-m1)/float(y-m1)],[0,0]],int,x,y)
            endelse
           end
        3: begin                ; *** THREE DIMENSIONAL ARRAY
            ; Only supports linear interpolation.
            arr_r = INTERPOLATE(arr, srx, sry, srz, /GRID)
           end
    endcase

    return, arr_r
end
