import os
import numpy as np

def w_cs(ion, sion):
    # ;returns array of weak s ions cross sections from various sources:
    kadonis=1
    kepler=0

    # ;-KADONIS MACS-

    if kadonis == 1:
        dirname = os.path.dirname(__file__)
        kadonisFile = os.path.join(dirname, 'kadonis_macs.dat')
        kadonis_macs = np.loadtxt(kadonisFile, usecols = 1)
        # kadonis_macs = kadonis_macs[1:-2]

        kad_s = kadonis_macs
        kad_n = np.loadtxt(kadonisFile, usecols = 0, dtype = "str")
        kad_n = np.char.lower(kad_n)
        # kad_n = np.empty(kadonis_macs.size, dtype = 'object')
        # kad_s = np.empty(kadonis_macs.size, dtype = 'object')
        # for i in range(0, kadonis_macs.size):
        #     kad_n[i] = kadonis_macs[i][0:5]
        #     kad_s[i] = kadonis_macs[i][5:14]


        sig = np.zeros(sion.size)
        for i in range(0, sion.size):
            sig[i] = kad_s[np.where(sion[i]==kad_n)]

    # ;-KEPLER RATES-
    # ;SCALE KEPLER RATES TO GET CROSS SECTIONS
    # ;SCALE CHOSEN TO MATCH CROSS SECTION FOR FE56 FROM:
    # ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    # ; Stellar Nucleosynthesis Data from the Tables of Reaction Rates        ;
    # ; for Nucleosynthsis Charged Particle, Weak, and Neutrino Interactions, ;
    # ; Version 92.1, by R.D. Hoffman and S.E. Woosley (1992).                ;
    # ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    # ; The 30 keV Cross section (mb) was calculated from the (n,g) reaction  ;
    # ; rates at T9=0.3 (rngp3) and T9=0.4 (rngp4) using the prescription of  ;
    # ; Woosley et al. (OAP-422, 1978), eq. 41                                ;
    # ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    # ;CROSS SECTION FE56(mb) = 13
    # ;NORMALIZE KEPLER RATES::

    if kepler == 0:
        dirname = os.path.dirname(__file__)
        ratesFile = os.path.join(dirname, 'rates.dat')
        rates = np.loadtxt(ratesFile, usecols = 1)
        kep_s = rates
        kep_n = np.loadtxt(ratesFile, usecols = 0, dtype = "str")
        kep_n = np.char.lower(kep_n)


        sig = np.zeros(sion.size)
        for i in range(0, sion.size):
            sig[i] = kep_s[np.where(sion[i] == kep_n)]
        # ;now normalize using the prescription above
        f1 = 13/sig[0]
        for i in range(0, sig.size):
            sig[i] *= f1


    return sig


'''
    kep_n = np.empty(rates.size, dtype = 'object')
    kep_s = np.zeros(rates.size)
    for i in range(0, rates.size):
        kep_n[i] = rates[i][0: 5]
        kep_s[i] = rates[i][5: 13]

    sig = np.zeros(sion.size)
    for i in range(0, sion.size):
        sig[i] = kep_s[np.where(sion[i] == kep_n)]

    # ;now normalize using the prescription above
    f1 = 12/sig[0]
    for i in range(0, sig.size):
        sig[i] *= f1
'''
        # return sig
