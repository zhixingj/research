import numpy as np
import os
import math
def gce_data(freb, ggg = None):
    dirname = os.path.dirname(__file__)+'/frebel'
    if freb == 0:
        mgfeFile = os.path.join(dirname, 'mgfe.dat')
        mgfe = np.loadtxt(mgfeFile, dtype = 'float')

        mgfeErrorFile = os.path.join(dirname, 'mgfe_error.dat')
        mgfe_sig = np.loadtxt(mgfeErrorFile, dtype = 'float')

        fehFile = os.path.join(dirname, 'feh_mg.dat')
        feh = np.loadtxt(fehFile, dtype = 'float')

        fehmgErrorFile = os.path.join(dirname, 'feh_mg_error.dat')
        feh_sig = np.loadtxt(fehmgErrorFile, dtype = 'float')

    if freb == 1:
        frebmgfeFile = os.path.join(dirname, 'frebelmg_mgfe.dat')
        frebmgfe = np.loadtxt(frebmgfeFile, skiprows = 2, dtype = 'object')

        frebfehFile = os.path.join(dirname, 'frebelmg_feh.dat')
        frebfeh = np.loadtxt(frebfehFile, skiprows = 2, dtype = 'object')

        mgfe = frebmgfe[0:frebmgfe.size-1]
        feh = frebfeh[0:frebfeh.size-1]

        #remove dwarfs
        dummymg = mgfe
        dummyfe = feh
        dummyfe = dummyfe.astype(np.float)
        # ;following omit boo,car,comBer,draco,fnx,g4,sextans,leo
        # ;dummymg[725:732]='*'
        # ;dummymg[751:778]='*'
        # ;dummymg[1161:1163]='*'
        # ;dummymg[1187:1209]='*'
        # ;following omit extra to be safe
        dummymg[725:854] = float('NaN')
        dummymg[1161:1223] = float('NaN')
        star = np.where(dummymg == '*')[0]
        dummymg[star] = float('NaN')
        dummymg = dummymg.astype(np.float)

        for i in range (0, mgfe.size):
            dummymg[i] = dummymg[i]
        
        good = np.where(~np.isnan(dummymg))[0]
        mgfe = np.zeros(good.size)
        feh = np.zeros(good.size)
        k = 0

        # ;omit data that is [Fe/H] <-4
        for i in range(0, good.size):
            if dummyfe[good[i]] >= -4.0:
                mgfe[k] = dummymg[good[i]]
                feh[k] = dummyfe[good[i]]
                k+=1

        # ;get soubiran data
        soubiranmgfeFile = os.path.join(dirname, 'soubiranmg_mgfe.dat')
        mgfe_s= np.loadtxt(soubiranmgfeFile, dtype = 'float')
        

        soubiranfehFile = os.path.join(dirname, 'soubiranmg_feh.dat')
        feh_s = np.loadtxt(soubiranfehFile, dtype = 'float')

        #mgfe.size = 724

        # ;convert from solag89 to sollo09
        # ;[Mg/Fe]_l = [Mg/Fe]_ag + log(Mg/Fe)_solag - log(Mg/Fe)_sollo09
        # ;          = [Mg/Fe]_ag + log(Mg_ag/Mg_lod) - log(Fe_ag/Fe_lod) <-- format used
        mgfe += np.log10(0.0006539/0.000678037) - np.log10(.001252/.001288)
        feh += np.log10(.001252/.001288)-np.log10(.7070/0.71126)

        # ;soubiran data is already converted from solag89
        # ;add frebel and soubiran into one array
        mgfe=np.hstack((mgfe,mgfe_s))
        feh=np.hstack((feh,feh_s))

        mgfe_sig = np.zeros(mgfe.size)
        feh_sig = np.zeros(mgfe.size)
        mgfe_sig[0 : mgfe_sig.size] = 0.1
        feh_sig[0 : feh_sig.size] = 0.1

        yr=[-1.0,2.0]

        data=np.zeros((mgfe.size,4))
        data[:, 0] = feh
        data[:, 1] = mgfe
        data[:, 2] = feh_sig
        data[:, 3] = mgfe_sig

        # END MAGNESIUM
        dirname = os.path.dirname(__file__)
        resultFile = os.path.join(dirname, 'gce_data.csv')
        header = '[Fe/H], [Mg/Fe], feh_sig, mgfe_sig'
        np.savetxt(resultFile, data, header = header, delimiter = ',')
        return data



