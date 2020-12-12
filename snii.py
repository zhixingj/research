import numpy as np
import os

def snii(ii, z, z_max, snii_z, snii_a, numiso, massive, ion, fp
):
    isoii = np.zeros((z_max,numiso))
    isol = np.zeros((z_max, numiso))

    #Get index for where metallicity array is certain values
    zero = np.where(abs(z - 1.0) == min(abs(z - 1.0)))
    zlow = np.where(abs(np.log10(z)+2.5) == min(abs(np.log10(z)+2.5)))
    #zlow=-2.184

    #FIND FE56 ABUNDANCE VALUE AT Z = -3 USING SNII ABUNDANCE FROM LODDERS
    fe56low = np.asarray(np.where((snii_z == 26) & (snii_a ==56)))
    felowabun = ii[fe56low]*10.0**(-2.5)

    #THE IDEA IS TO NORMALIZE ALL FE PEAK ELEMENTS TO THIS ABUNDANCE

    #GET NORMALIZATION FACTOR
    NormFactor = felowabun/massive[fe56low]

    #NOW NORMALIZE ALL FE PEAK HEGER DATA BY THIS VALUE 
    for i in range(0,numiso):
        massive[i]*=NormFactor

    dirname = os.path.dirname(__file__)
    massiveFile = os.path.join(dirname, 'massivedatascaledforminus3.dat')
    np.savetxt(massiveFile, massive)
        

    # NOW WE WANT TO LINEARIZE THESE DATA POINTS AT Z = -2.5 TO THE LODDERS ABUN FOR FE PEAK ISOTOPES

    #FIRST GET THE MATRIX OF SLOPES
    Slopes = np.zeros((numiso, 1))
    Slopesl = np.zeros((numiso, 1))
    for i in range(0,numiso):
        Slopesl[i] = (ii[i] - massive[i])/(z[zero] - z[zlow])
        Slopes[i] = (np.log10(ii[i]) - np.log10(massive[i]))/(np.log10(z[zero]) - np.log10(z[zlow]))

    Slopes[0:9]=0.0
    Slopes[76:287]=0.0

    Slopesl[0:9]=0.0
    Slopesl[76:287]=0.0

    ''';openw,lun,'~/gce/slopeslog.dat',/get_lun
    ;for i=0,n_elements(Slopes)-1 do begin
    ;printf,lun,Slopes[i]
    ;endfor
    ;close,lun
    ;free_lun,lun'''

    #NOW MAKE MATRIX OF ABUNDANCES
    for i in range(0, z_max):
        for j in range(0, numiso):
            isoii[i, j] = 10**(Slopes[j]*(np.log10(z[i]) - np.log10(z[zlow])) + np.log10(massive[j]))

        #LINEAR SCALING IN LINEAR SPACE -- FOR COMPARISON;;;;;;;;;;;
            if np.log10(z[i]) >= -2.5:
                isol[i, j] = Slopesl[j]*(z[i] - z[zlow]) + massive[j]
            if np.log10(z[i]) <= -2.5:
                isol[i, j] =  massive[j]*(z[i]/z[zlow])
        #LINEAR SCALING IN LINEAR SPACE -- FOR COMPARISON;;;;;;;;;;;          


    ''';openpsfl,'logVlin.ps'
    ;plot,alog10(z),alog10(isol(13,*)),xtitle='log(xi)',ytitle='Log(16O) (log:red lin:black)',charsize=1.5,title='Log vs Linear'
    ;oplot,alog10(z),alog10(isoii(13,*)),color=rgb(1.,0,0)
    ;closeps'''
    return isoii
    