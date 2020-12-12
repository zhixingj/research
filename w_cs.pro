FUNCTION w_cs,ion,sion

;returns array of weak s ions cross sections from various sources:
kadonis=1
kepler=0


;00000000000000000000000000000000000000000000000000000000000000
;-KADONIS MACS-

IF kadonis EQ 1 THEN BEGIN

str=''
OPENR,lun,'~/Desktop/gce/sprocess/kadonis_macs.dat',/GET_LUN
while not eof(lun) do begin
    kep=''
    readf,lun,kep,FORMAT="(A)"
    str=[str,kep]
endwhile
d=strarr(n_elements(str)-1)
d=str[1:n_elements(str)-2]
str=d
close,lun
free_lun,lun

kad_n=strarr(n_elements(str))
kad_s=dblarr(n_elements(str))
for i=0,n_elements(str)-1 do begin
  kad_n[i]=strsplit(strmid(str[i],0,5),/extract)
  kad_s[i]=double(strmid(str[i],5,9))
endfor

sig=dblarr(n_elements(sion))
for i=0,n_elements(sion)-1 do begin
  sig[i]=kad_s(where(sion[i] eq kad_n))
endfor

ENDIF

;00000000000000000000000000000000000000000000000000000000000000





;00000000000000000000000000000000000000000000000000000000000000
;-KEPLER RATES-
;SCALE KEPLER RATES TO GET CROSS SECTIONS
;SCALE CHOSEN TO MATCH CROSS SECTION FOR FE56 FROM:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Stellar Nucleosynthesis Data from the Tables of Reaction Rates        ;
; for Nucleosynthsis Charged Particle, Weak, and Neutrino Interactions, ;
; Version 92.1, by R.D. Hoffman and S.E. Woosley (1992).                ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; The 30 keV Cross section (mb) was calculated from the (n,g) reaction  ;
; rates at T9=0.3 (rngp3) and T9=0.4 (rngp4) using the prescription of  ;
; Woosley et al. (OAP-422, 1978), eq. 41                                ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CROSS SECTION FE56(mb) = 13
;NORMALIZE KEPLER RATES::

IF kepler EQ 1 THEN BEGIN

str=''
OPENR,lun,'~/Desktop/gce/sprocess/rates.dat',/GET_LUN
while not eof(lun) do begin
    kep=''
    readf,lun,kep,FORMAT="(A)"
    str=[str,kep]
endwhile
d=strarr(n_elements(str)-1)
d=str[1:n_elements(str)-2]
str=d
close,lun
free_lun,lun

kep_n=strarr(n_elements(str))
kep_s=dblarr(n_elements(str))
for i=0,n_elements(str)-1 do begin
  kep_n[i]=strsplit(strmid(str[i],0,5),/extract)
  kep_s[i]=double(strmid(str[i],5,8))
endfor

sig=dblarr(n_elements(sion))
for i=0,n_elements(sion)-1 do begin
  sig[i]=kep_s(where(sion[i] eq kep_n))
endfor

;now normalize using the prescription above
f1=13d0/sig[0]
for i=0,n_elements(sig)-1 do begin 
  sig[i]*=f1
endfor

ENDIF

;0000000000000000000000000000000000000000000000000000000000

RETURN,sig

END
