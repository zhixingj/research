PRO ws

;CALCULATES FULL S PROCESS NUCLEOSYNTHESIS VIA MATRIX METHOD
;USES b+/b- RATIA FROM KEPLER DECAY.DAT 
;USES N CAPTURE CROSS SECTIONS FROM KEPLER
;USES n/b- RATIA CALCULATED FROM HALF LIVES 

erase
!P.BACKGROUND = 'FFFFFF'XL   ;Set background color to white
!P.COLOR = 'OOOOOO'XL  
erase

;000000000000000000000
;DO WEAK S PATH FIRST
;000000000000000000000

numiso=287

;GET SOLAR ABUN
OPENR,lun,'~/Desktop/gce/sollo09.dat', /GET_LUN
header = strarr(28)
READF, lun, header
abu_array = strarr(287)
READF, lun, abu_array
close,lun
free_lun,lun

abu_array2 = strarr(287*2)

for i=0,numiso-1 do begin
abu_array(i) = STRMID(abu_array(i),5,15)
abu_array(i) = DOUBLE(abu_array(i))
endfor

;MAKE ION ARRAY
OPENR,lun,'~/Desktop/gce/isotope_As.dat',/GET_LUN
isomass=DBLARR(numiso)
READF,lun,isomass
close,lun
free_lun,lun

OPENR,lun,'~/Desktop/gce/isotope_names.dat',/GET_LUN
isoname=STRARR(287)
READF,lun,isoname
close,lun
free_lun,lun

ion=string(uint(isomass))
for i=0,numiso-1 do begin
ion[i]=isoname[i]+strsplit(ion[i],/extract)
endfor

sr86=where(ion eq 'Sr86')
sr87=where(ion eq 'Sr87')

;GET MAIN S PROCESS ABUNDANCES
;openr,lun,'~/Desktop/gce/sfrac_arlandini.dat'
;openr,lun,'~/Desktop/gce/sfrac_gallino.dat'
openr,lun,'~/Desktop/gce/sfrac_sollo09_gallino.dat'
header=strarr(1)
readf,lun,header
sfrac=dblarr(numiso)
readf,lun,sfrac
close,lun
free_lun,lun
mains=sfrac*abu_array

;get weak s ion targets in an array
sion=weak_s_target(ion)

;get weak s ion daughters in an array
dion=weak_s_daughter(ion)

;get unique iso names for defining matrix rows and columns
array=sion[uniq(sion)]
uniq=[array,dion[n_elements(dion)-1]]

;get indicies for weak s path relative to full ion array
;these are used for:
;reconstituting final weak s into ion format
;parsing solar abun and isomass to compare weak s with them
;note s_index will NOT contain unstable isos 
s_index=weak_s_index(uniq,ion)

;TAKE SOLAR TO BE SOLAR - MAIN S-PROCESS, FIT THE WEAK S-PROCESS TO THIS!
;TO UNDO, JUST REMOVE -MAINS(XXX) TERMS FROM BELOW 

  ge70=where(ion eq 'Ge70')
  se76=where(ion eq 'Se76')
  kr80=where(ion eq 'Kr80')
  kr82=where(ion eq 'Kr82')
  sr86=where(ion eq 'Sr86')
  sr87=where(ion eq 'Sr87')

  Ge70=where(uniq eq 'Ge70')
  Se76=where(uniq eq 'Se76')
  Kr80=where(uniq eq 'Kr80')
  Kr82=where(uniq eq 'Kr82')
  Sr86=where(uniq eq 'Sr86')
  Sr87=where(uniq eq 'Sr87')

;get solar/mass arrays (careful: must use abundances after decays)
; Se79 -> Br89
; Zr93 -> Nb93
; Tc99 -> Ru99
solar=abu_array(s_index);-mains(s_index)
mass=isomass(s_index)
;Se79=where(uniq eq 'Se79')
Zr93=where(uniq eq 'Zr93')
Tc99=where(uniq eq 'Tc99')
;Br79=where(ion eq 'Br79')
Nb93=where(ion eq 'Nb93')
Ru99=where(ion eq 'Ru99')
;solar[Se79]=abu_array[Br79];-mains(Br79)
solar[Zr93]=abu_array[Nb93];-mains(Nb93)
solar[Tc99]=abu_array[Ru99];-mains(Ru99)
;mass[Se79]=79d0
mass[Zr93]=93d0
mass[Tc99]=99d0

solarM=abu_array(s_index)-mains(s_index)
mass=isomass(s_index)
;Se79=where(uniq eq 'Se79')
Zr93=where(uniq eq 'Zr93')
Tc99=where(uniq eq 'Tc99')
;Br79=where(ion eq 'Br79')
Nb93=where(ion eq 'Nb93')
Ru99=where(ion eq 'Ru99')
;solarM[Se79]=abu_array[Br79]-mains(Br79)
solarM[Zr93]=abu_array[Nb93]-mains(Nb93)
solarM[Tc99]=abu_array[Ru99]-mains(Ru99)
;mass[Se79]=79d0
mass[Zr93]=93d0
mass[Tc99]=99d0

;15% Kr80 and 3% Kr82 made by p-process
solarM[Kr80] = solarM[Kr80]-0.15d0*solar[Kr80]
solarM[Kr82] = solarM[Kr82]-0.03d0*solar[Kr82]

submain=1
zero=where(solar eq 0d0)


;get weak s cross sections: they match sion array
sig=w_cs(ion,sion)

;get branching ratios
b=branch(sion)

;get coefficient matrix
a=matrix(uniq,sion,dion,sig,b)

;q=strarr(1,n_elements(sion))
;q[0,*]=sion[*]
;print,q

;USING dN/dt = A*N, so N is column matrix

;0000000000000000000000000000
;BEGIN LOOP FOR FINDING TAU
;0000000000000000000000000000

;tau=dindgen(1000)/1000+0.01
;tau[0]=tau[1]/2d0
;i=0

tt=dindgen(100)/100+.01
ttgood=dblarr(n_elements(tt))
;FOR x=0,n_elements(tt)-1 do begin
;tau=[.4,.3,tt[x]]
;tau=[.1,.2,.25,.3,.35,.4]
;tau=[.1,.15,.18,.22,.25,.3,.33,.35,.37,.4]
tau=[.2,.25,.3,.35,.4]
;tau=[.1,.2,.3,.4]
;tau=[.2,.3]
sr=dblarr(2,n_elements(tau))
sr=string(sr)
;a=transpose(a)
piece=0
if piece eq 1 then begin
upp=4
a=a[0:upp,0:upp]
uniq=uniq[0:upp]
solar=solar[0:upp]
abu_array=abu_array[0:upp]
sion=sion[0:upp]
dion=dion[0:upp]
endif
aa=dblarr(n_elements(tau),6)
bb=[[1d0],[1d0],[1d0],[1d0],[1d0],[1d0]]
sfits=dblarr(n_elements(tau),n_elements(solar))

FOR i=0,n_elements(tau)-1 DO BEGIN
  m=a*tau[i]
  eval=double(la_eigenproblem(m,EIGENVECTORS=evec))
  evec=double(evec)
    ; Check the results using the eigenvalue equation:  
    maxErr = 0d  
    for j = 0, n_elements(uniq)-1 do begin 
    ; A*z = lambda*z  
    alhs = m ## evec[*,j]  
    arhs = eval[j]*evec[*,j]  
    maxErr = maxErr > MAX(ABS(alhs - arhs))  
    endfor 
  expA = evec#diag_matrix(exp(eval))#invert(evec)
  ;n=expA[0,*]
  n=(expA[*,0])
  
  if submain eq 1 then begin
    if zero ne -1 then  n(zero) = 0d0
  endif
  Ge70=where(uniq eq 'Ge70')
  Se76=where(uniq eq 'Se76')
  Kr80=where(uniq eq 'Kr80')
  Kr82=where(uniq eq 'Kr82')
  Sr86=where(uniq eq 'Sr86')
  Sr87=where(uniq eq 'Sr87')
if piece eq 0 then begin 
  ;COMPARE WITH SOLAR
  ;PICK OUT S-ONLY ISOTOPES TO MATCH
  if max(n/solar) eq n[Ge70]/solar[Ge70] then begin
    print,tau[i],'Ge70'
    sr[0,i]='Ge70'
    sr[1,i]=string(tau[i])    
    
  endif  
  if max(n/solar) eq n[Se76]/solar[Se76] then begin
    print,tau[i],'Se76'
    sr[0,i]='Se76'
    sr[1,i]=string(tau[i])   
    
  endif 
  if max(n/solar) eq n[Kr80]/solar[Kr80] then begin
    print,tau[i],'Kr80'
    sr[0,i]='Kr80'
    sr[1,i]=string(tau[i]) 
      
  endif 
  if max(n/solar) eq n[Kr82]/solar[Kr82] then begin
    print,tau[i],'Kr82'
    sr[0,i]='Kr82'
    sr[1,i]=string(tau[i]) 
      
  endif 
  if max(n/solar) eq n[Sr86]/solar[Sr86] then begin
    print,tau[i],'Sr86'
    sr[0,i]='Sr86'
    sr[1,i]=string(tau[i])
 
  endif
  if max(n/solar) eq n[Sr87]/solar[Sr87] then begin
    print,tau[i],'Sr87'
    sr[0,i]='Sr87'
    sr[1,i]=string(tau[i])
 
  endif
endif

  ;fac=n[Sr86]/solar[Sr86]
  
  ;convert n to mass fractions
  ;n=atoms2massfrac(n,solar,abu_array)
  

  if submain eq 0 then begin
    fac=max(n/solar)
    maxx=where(n/solar eq fac)
    for k=0,n_elements(uniq)-1 do begin
      n[k]=n[k]/fac
    endfor 
    ratio=n/solar
    plot,mass,alog10(n/solar),psym=2,xrange=[mass[0],mass[n_elements(sion)-1]]
  endif

  if submain eq 1 then begin
    nonzero=where(n ne 0d0 and solarM ne 0d0)
    fac=max(n[nonzero]/solarM[nonzero])
    maxx=where(n/solarM eq fac)
    for k=0,n_elements(uniq)-1 do begin
      n[k]=n[k]/fac
    endfor 
    ratio=n/solarM
    plot,mass,alog10(n/solarM),psym=2,xrange=[mass[0],mass[n_elements(uniq)-1]]  
    ;plot,mass,alog10(n),psym=2,xrange=[mass[0],mass[n_elements(uniq)-1]] 
    vline,70
    vline,76
    vline,80
    vline,82
    vline,86
    vline,87
    ;print,tau[i]
    ;stop
  endif

;if tau[i] eq .277d0 then stop
   ;plot,mass,alog10(n),psym=2,xrange=[mass[0],mass[n_elements(sion)-1]]
  ;print,'----'
  ;print,i,tau[i]
  ;print,ratio[Ge70]
  ;print,ratio[Se76]
  ;print,ratio[Kr82]
  ;print,ratio[Sr86]
  ;print,'----'
  
  aa[i,*]=[ratio[Ge70],ratio[Se76],ratio[Kr80],ratio[Kr82],ratio[Sr86],ratio[Sr87]]
  sfits[i,*]=ratio

ENDFOR	

;a2=0.96959363d0 ;.2
;a3=0.0d0 ;.25
;a4=1.3878442d0 ;.3

ttee=1
if ttee eq 1 then begin
;a1=0.0d0 ;.1
;a2=1.5d0 ;.2
;a3=0.2d0 ;.25
;a4=0.4d0 ;.3
;a5=0.6d0 ;.35
;a6=0.7d0 ;.4
a1=1.5d0 ;.2
a2=0.2d0 ;.25
a3=0.4d0 ;.3
a4=0.6d0 ;.35
a5=0.7d0 ;.4
t1=sfits[0,*]*a1
t2=sfits[1,*]*a2
t3=sfits[2,*]*a3
t4=sfits[3,*]*a4
t5=sfits[4,*]*a5
;plot,mass,t1+t2+t3+t4+t5+t6,psym=2,xrange=[mass[0],mass[n_elements(uniq)-1]],yrange=[0,2]  
plot,mass,alog10(t1+t2+t3+t4+t5),psym=2,xrange=[mass[0],mass[n_elements(uniq)-1]],yrange=[-2,0.5] 
hline,0,line=2
vline,70
vline,76
vline,80
vline,82
vline,86
vline,87
array=t1+t2+t3+t4+t5
stop
endif

tee=0
if tee eq 1 then begin
cc=la_least_squares(aa,bb)
print,la_least_squares(aa,bb)
array=dblarr(n_elements(solar))
for i=0,n_elements(tau)-1 do begin
for j=0,n_elements(solar)-1 do begin
  array[j]+=cc[i]*sfits[i,j]
endfor
endfor
plot,mass,alog10(array),psym=2
hline,0
vline,70
vline,76
vline,80
vline,82
vline,86
vline,87
endif
stop
openpsfl,'~/Desktop/gce/plots/weak_s/GCE_PUBL_weaks_solar_minus_main_3_14_2012.ps'
xr=[55,102]
yr=[-3.5,0.5]
s_only_shape=6
plot,xr,yr,/nodata,XTITLE='Mass Number',YTITLE='Log[weak/(solar-main)]',/YS,/XS
oplot,mass[0:Ge70-1],alog10(array[0:Ge70-1]),psym=sym(12)
oplot,mass[Ge70+1:Se76-1],alog10(array[Ge70+1:Se76-1]),psym=sym(12)
oplot,mass[Se76+1:Kr80-1],alog10(array[Se76+1:Kr80-1]),psym=sym(12)
oplot,mass[Kr80+1:Kr82-1],alog10(array[Kr80+1:Kr82-1]),psym=sym(12)
oplot,mass[Kr82+1:Sr86-1],alog10(array[Kr82+1:Sr86-1]),psym=sym(12)
oplot,mass[Sr86+1:Sr87-1],alog10(array[Sr86+1:Sr87-1]),psym=sym(12)
oplot,mass[Sr87+1:n_elements(mass)-1],alog10(array[Sr87+1:n_elements(mass)-1]),psym=sym(12)
oplot,mass[Ge70],alog10(array[Ge70]),psym=s_only_shape,color=rgb(1.0,0,0)
oplot,mass[Se76],alog10(array[Se76]),psym=s_only_shape,color=rgb(1.0,0,0)
oplot,mass[Kr80],alog10(array[Kr80]),psym=s_only_shape,color=rgb(1.0,0,0)
oplot,mass[Kr82],alog10(array[Kr82]),psym=s_only_shape,color=rgb(1.0,0,0)
oplot,mass[Sr86],alog10(array[Sr86]),psym=s_only_shape,color=rgb(1.0,0,0)
oplot,mass[Sr87],alog10(array[Sr87]),psym=s_only_shape,color=rgb(1.0,0,0)
hline,0,line=2
LEGENDE,/BOTTOM,/RIGHT,/DRAWFRAME,FRAMETOP=2,FRAMERIGHT=7,FRAMELEFT=5,CHARSIZE=1.2,XPOS=1.1,YPOS=0.0
legende_symbol,6,'S-only',color=rgb(1.0,0,0)
legende_symbol,7,'S & R'
closeps
stop

togg=0
for p=0,n_elements(cc)-1 do begin
  if cc[p] gt 0.0 then togg++
endfor
if togg eq 3 then begin
ttgood[x]=tt[x]
endif

;ENDFOR ;TT LOOP


toobig=where(array gt 1d0)
array[toobig]=1d0

;MUST RENAME SION SO THEY MATCH THEIR DECAY PRODUCTS
; Se79 -> Br89
; Zr93 -> Nb93
; Tc99 -> Ru99
;uniq[Se79]='Br79'
uniq[Zr93]='Nb93'
uniq[Tc99]='Ru99'

array=array*solarM/solar

s_index=weak_s_index(uniq,ion)
wsfrac=DBLARR(numiso)
for i=0,n_elements(s_index)-1 do begin
  wsfrac[s_index[i]]=array[i]
;  wsfrac[s_index[i]]=wsabun[i]
endfor

ga69=where(ion eq 'Ga69')
wsabun=wsfrac
p_in=get_p(ion)
rfrac=dblarr(numiso)
for i=ga69[0],numiso-1 do begin
  rfrac[i]=1d0-(sfrac[i]+wsfrac[i])
endfor
rfrac[p_in]=0d0

;openpsfl,'~/Desktop/gce/sprocess/MAIN_WEAK_COLOR.ps'
xr=[0,240]
yr=[-5,1]
plot,xr,yr,/nodata,TITLE='Red: Weak   Blue: Main',XTITLE='Mass Number',YTITLE='Log(abun/solar)',/YS,/XS
hline,0
oplot,mass,alog10(array),color=rgb(1d0,0d0,0d0)
nozero=where(sfrac gt 0.0)
oplot,isomass[nozero],alog10(sfrac[nozero]),color=rgb(0,0,1d0)
nozero=where(rfrac gt 0.0)
oplot,isomass[nozero],alog10(rfrac[nozero]),color=rgb(0,1d0,0)

;closeps
stop

plot,isomass,alog10(wsfrac),psym=2,yrange=[-8,1]
oplot,isomass,alog10(sfrac),psym=6



stot=wsfrac+sfrac
big=where(stot gt 1.0)
if big ne -1 then begin
  print,'[WS.PRO] overproductions!'
  print, ion[big]
  stop
  overprod=stot[big]-1d0
  wsfrac[big]-=overprod
endif


plot,isomass,alog10(wsfrac+sfrac),psym=2,yrange=[-8,1]


wsfrac=string(wsfrac)
for i=0,numiso-1 do begin
  wsfrac[i]=strsplit(wsfrac[i],/extract)
endfor

stop
openw,lun,'~/Desktop/gce/wsfrac_gallino.dat',/get_lun
;printf,lun,';WS abundance fractions relative to solar, average of tau=[.3,.4], exp(matrix) method employed. 6/28/2011' 
printf,lun,';WS abundance fractions relative to solar, fit with tau=[.2,.25,.3,.35,.4] & cc=[1.5,0.2,0.4,0.6,0.7], exp(matrix) method employed. 7/5/2011' 
for i=0,numiso-1 do begin
  printf,lun,wsfrac[i]
endfor
close,lun
free_lun,lun

stop

good=where(ttgood ne 0d0)
print,'--------------'
print,ttgood[good]
print,'--------------'

kk=double(sr[1,*])
pgood=where(kk ne 0d0)
srgood=strarr(2,n_elements(pgood))
srgood[0,*]=sr[0,pgood]
srgood[1,*]=sr[1,pgood]
;print,'-----------------'
;print,srgood
;print,'-----------------'

;0000000000000000000000000000
;END LOOP FOR FINDING TAU
;0000000000000000000000000000

STOP

;EXPORT TO DATA FILE
;MUST RENAME SION SO THEY MATCH THEIR DECAY PRODUCTS
; Se79 -> Br89
; Zr93 -> Nb93
; Tc99 -> Ru99
uniq[Se79]='Br79'
uniq[Zr93]='Nb93'
uniq[Tc99]='Ru99'

;FIX s_index TO REFLECT UNSTABLE -> STABLE

;----finish export here
openw,lun,'~/Desktop/gce/weaksabun.dat',/get_lun
for i=0,n_elements(n)-1 do begin
  printf,lun,uniq[i],n[i]
endfor
close,lun
free_lun,lun




STOP
END