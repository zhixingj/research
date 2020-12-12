function snii, ii, z, z_max, snii_z, snii_a, numiso, massive, ion, fp

isoii = DBLARR(numiso,z_max)
isol = DBLARR(numiso,z_max)

;Get index for where metallicity array is certain values
zero = where(abs(z - 1.0) eq min(abs(z - 1.0)))
zlow = where(abs(alog10(z)+2.5) eq min(abs(alog10(z)+2.5)))
;zlow=-2.184d0

;FIND FE56 ABUNDANCE VALUE AT Z = -3 USING SNII ABUNDANCE FROM LODDERS
fe56low = where(snii_z eq 26 and snii_a eq 56)
felowabun = ii(fe56low)*10.0^(-2.5)

;
;THE IDEA IS TO NORMALIZE ALL FE PEAK ELEMENTS TO THIS ABUNDANCE
;

;GET NORMALIZATION FACTOR
NormFactor = felowabun/massive(fe56low)

;NOW NORMALIZE ALL FE PEAK HEGER DATA BY THIS VALUE 
for i=0,numiso-1 do begin
  massive[i]*=NormFactor
endfor

openw,lun,'~/gce/massivedatascaledforminus3.dat',/get_lun
;openw,lun,'~/gce/massivedatascaledforminus3BU_11_29_2012.dat',/get_lun ;(cayrelfit linear)
for i=0,n_elements(massive)-1 do begin
printf,lun,massive[i]
endfor
close,lun
free_lun,lun

;
;NOW WE WANT TO LINEARIZE THESE DATA POINTS AT Z = -2.5 TO THE LODDERS ABUN FOR FE PEAK ISOTOPES
;

;FIRST GET THE MATRIX OF SLOPES
Slopes = DBLARR(1,numiso)
Slopesl=DBLARR(1,numiso)
for i = 0,numiso-1 do begin
   Slopesl(i) = (ii(i) - massive(i))/(z(zero) - z(zlow))
   Slopes(i) = (alog10(ii(i)) - alog10(massive(i)))/(alog10(z(zero)) - alog10(z(zlow)))
endfor

Slopes[0:8]=0.d0
Slopes[76:286]=0.d0

Slopesl[0:8]=0.d0
Slopesl[76:286]=0.d0

;openw,lun,'~/gce/slopeslog.dat',/get_lun
;for i=0,n_elements(Slopes)-1 do begin
;printf,lun,Slopes[i]
;endfor
;close,lun
;free_lun,lun


;NOW MAKE MATRIX OF ABUNDANCES
for i = 0l,z_max-1 do begin
   for j=0,numiso-1 do begin
      isoii(j,i) = 10^(Slopes(j)*(alog10(z(i)) - alog10(z(zlow))) + alog10(massive(j)))

    ;LINEAR SCALING IN LINEAR SPACE -- FOR COMPARISON;;;;;;;;;;;
      if alog10(z(i)) ge -2.5 then begin                       ;
          isol(j,i) = Slopesl(j)*(z(i) - z(zlow)) + massive(j) ;
      endif                                                    ;
      if alog10(z(i)) lt -2.5 then begin                       ;
          isol(j,i) =  massive(j)*(z(i)/z(zlow))               ;
      endif                                                    ;
    ;LINEAR SCALING IN LINEAR SPACE -- FOR COMPARISON;;;;;;;;;;;          
   endfor
endfor

;openpsfl,'logVlin.ps'
;plot,alog10(z),alog10(isol(13,*)),xtitle='log(xi)',ytitle='Log(16O) (log:red lin:black)',charsize=1.5,title='Log vs Linear'
;oplot,alog10(z),alog10(isoii(13,*)),color=rgb(1.,0,0)
;closeps

;stop
RETURN, isoii

END
