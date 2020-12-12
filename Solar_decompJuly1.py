import numpy as np
import os
def solar_decomp(abu_array, amu, iso_zs, fp, ion, massive):
    numiso = 287
    h1 = 0.75128 / 1.007825
    h2 = 4.31998e-5 / 2.014102
    he3 = 2.13003e-5 / 3.016029
    he4 = 0.248659 / 4.002603
    li7 = 1.912209E-09 / 7.016004


    bbn = np.zeros(numiso)
    h_burn = np.zeros(numiso)
    gcr = np.zeros(numiso)
    novae_priGCR_nu = np.zeros(numiso)


    bbn[0:4] =  np.array([h1, h2, he3, he4])
    bbn[5] = li7

    h_burn[2:4] = [abu_array[2]-bbn[2], abu_array[3]-bbn[3]]

    gcr[4] = 0.3 * abu_array[4] #Li6       
    gcr[5] = 0.3 * abu_array[5] #Li7
    gcr[6] = 0.25 * abu_array[6] #Be9
    gcr[7] = 0.25 * abu_array[7] #B10
    gcr[8] = 0.25 * abu_array[8] #B11

    novae_priGCR_nu[4] = abu_array[4] - gcr[4]
    novae_priGCR_nu[5] = abu_array[5] - gcr[5] - bbn[5]
    novae_priGCR_nu[6] = abu_array[6] - gcr[6]
    novae_priGCR_nu[7] = abu_array[7] - gcr[7]
    novae_priGCR_nu[8] = abu_array[8] - gcr[8]

    light = np.zeros(numiso)
    light[0:9] = abu_array[0:9]

    #END LIGHTEST, BEGIN P_IOTOPES

    p_index = np.array([85,93,101,112,113,119,120,127,135,136,143,145,146,147,157,166,167,176,
            177,183,185,186,197,214,215,222,223,229,238,244,246,253,262,269])

    ru98 = np.where(ion == 'ru98')
    pd102 = np.where(ion == 'pd102')
    er164 = np.where(ion == 'er164')

    #DEFINE ARRAYS


    nproc = np.zeros(numiso)
    gproc = np.zeros(numiso)

    i_index = np.arange(287)
  
    cut =  np.where(p_index == ru98[0])
    mix = np.where(p_index == pd102[0])
    low = np.where(p_index <= ru98[0])
    high = np.where(p_index > pd102[0])
    nproc[p_index[low]] = abu_array[p_index[low]]
    nproc[p_index[mix]] = 0.25 * abu_array[p_index[mix]]
    gproc[p_index[mix]] = 0.75 * abu_array[p_index[mix]]
    gproc[p_index[high]] = abu_array[p_index[high]]

    #15% Kr80 and 3% Kr82 made by p-process
    Kr80 = np.where(ion == 'kr80')
    Kr82 = np.where(ion == 'kr82')
    nproc[Kr80] = 0.15 * abu_array[Kr80]
    nproc[Kr82] = 0.03 * abu_array[Kr82]

    #50% of Se74, Kr78, Sr84 made by photodisintegration on weak s-process isotopes, assign as weak s-process done below
    Kr78 = np.where(ion == 'kr78')
    Se74 = np.where(ion == 'se74')
    Sr84 = np.where(ion == 'sr84')
    nproc[Kr78] = 0.5 * abu_array[Kr78]
    nproc[Se74] = 0.5 * abu_array[Se74]
    nproc[Sr84] = 0.5 * abu_array[Sr84]

    #END P-ISOTOPES

    #BEGIN WEAK S, MAIN S, AND R

    sproc = np.zeros(numiso)
    wsproc = np.zeros(numiso)
    rproc = np.zeros(numiso)

    #GET MAIN S-PROCESS YIELDS
    dirname = os.path.dirname(__file__)
    sfracFile = os.path.join(dirname, 'sfrac_sollo09_gallino.dat')
    sfrac = np.loadtxt(sfracFile, skiprows = 1)

    sproc = sfrac * abu_array

    dirname = os.path.dirname(__file__)
    wsfracFile = os.path.join(dirname, 'wsfrac_gallino.dat')
    wsfrac = np.loadtxt(wsfracFile, skiprows = 1)

    # fix all NaNs/negative values
    fixws = np.where(wsfrac != np.absolute(wsfrac))
    wsfrac[fixws] = 0

    wsproc = wsfrac * abu_array

    # 50% of Se74, Kr78, Sr84 made by photodisintegration on weak s-process isotopes, assign as weak s-process
    wsproc[Kr78] = 0.5 * abu_array[Kr78]
    wsproc[Se74] = 0.5 * abu_array[Se74]
    wsproc[Sr84] = 0.5 * abu_array[Sr84]

    #fix Kr80 overabun & rescale Kr82 underabun
    wsproc[Kr80] = abu_array[Kr80] - sproc[Kr80] - nproc[Kr80]
    rat = sproc[Kr82] / wsproc[Kr82]
    left = abu_array[Kr82] - wsproc[Kr82] - sproc[Kr82] - nproc[Kr82]
    sproc[Kr82] += rat * left
    wsproc[Kr82] += (1 - rat) * left

    #mesh weak and main together, take r-process as the residual for ion > Ga69
    #For ion < Ga69, no rproc, only massive & Ia
    #scale so that isotopic ratios preserved across sproc and wsproc

    toobig = np.where(wsproc > abu_array)
    if toobig != -1:
        wsproc[toobig] = abu_array[toobig]


    ratio = (wsproc + sproc + nproc) / abu_array
    overabun = np.where(ratio > 1)
    underabun = np.where(ratio <= 1)

    ga69 = np.where(ion == 'ga69')

    #These s-only isotopes are underproduced by Gallino.
    gd152 = np.where(ion == 'gd152')
    os187 = np.where(ion == 'os187')
    pt192 = np.where(ion == 'pt192')

    # set gamma process for these equal to solar - main s-process
    gproc[gd152] = abu_array[gd152] - sproc[gd152]
    gproc[os187] = abu_array[os187] - sproc[os187]
    gproc[pt192] = abu_array[pt192] - sproc[pt192]

    # Er164 is made by Gallino, which is proton rich
    gproc[er164] = abu_array[er164] - sproc[er164]

    # FIX TA180, W180, SN115, CD108, MO94, which are also made by Gallino, which are proton rich
    ta180 = np.where(ion == 'ta180')
    w180 = np.where(ion == 'w180')
    sn115 = np.where(ion == 'sn115')
    cd108 = np.where(ion == 'cd108')
    mo94 = np.where(ion == 'mo94')
    sproc[ta180] =0
    sproc[w180] = 0
    gproc[sn115] = abu_array[sn115] - sproc[sn115]
    gproc[cd108] = abu_array[cd108] - sproc[cd108]
    nproc[mo94] = abu_array[mo94] - sproc[mo94]

    for i in range(76, numiso):
        rproc[i] = abu_array[i] - wsproc[i] - sproc[i] - nproc[i] - gproc[i]

    #set wsproc for fe56 to zero
    #in reality there is a small (negligible) contribution
    wsproc[61]=0

    # END WEAK S, MAIN S, AND R


    #FIT MASSIVE, WEAK S, MAIN S, AND TYPE IA

    #GET TYPE IA YIELDS
    dirname = os.path.dirname(__file__)
    w7File = os.path.join(dirname, 'w7.dat')
    w7 = np.loadtxt(w7File, skiprows = 1)


    #convert w7 to mol frac
    w7 = w7 / amu

    #fe57 =< A <= zn68: take massive contributions as residuals of sproc and wsproc
    masshigh = np.zeros(numiso)
    fe56 = np.where(ion == 'fe56')
    ga69 = np.where(ion == 'ga69')
    for i in range(fe56[0][0] + 1,76):
        masshigh[i]=abu_array[i]-(wsproc[i]+sproc[i])


    #A<=fe56: PERFORM 2 SCALINGS
    #FIRST SCALING: SCALE ACCORDING TO TYPE IA FE-FRACTION PARAMETER
    masslow = np.zeros(numiso, dtype="float")
    m_fac = (1 - fp) * abu_array[fe56] / massive[fe56]
    i_fac = fp * abu_array[fe56] / w7[fe56]


    for i in range(0,fe56[0][0]+1):
        masslow[i] = m_fac * massive[i]
        w7[i] = i_fac * w7[i]
        
     

    #binwrite, array=masslow, filename='massivePRE'
    #binwrite, array=w7, filename='onePRE'

    #SECOND SCALING: PRESERVE ISOTOPIC RATIOS BETWEEN MASSIVE AND TYPE IA BETWEEN c12 and fe56
    #ZERO ALL MASSIVE CONTRIBUTIONS PAST ru102
    c12 = np.where(ion == 'c12')
    fe56 = np.where(ion == 'fe56')
   


    w_new = np.zeros(numiso)
    m_new = np.zeros(numiso)

    for i in range(c12[0][0], fe56[0][0] + 1):
        w_new[i] = abu_array[i] / (masslow[i] / w7[i] + 1)
        m_new[i] = abu_array[i] / (w7[i] / masslow[i] + 1)
       

    #binwrite,array=m_new,filename='massivePOST'
    #binwrite,array=w_new,filename='onePOST'
    #stop

    '''for i=c12[0],fe56[0] do begin
    masslow(i) = masslow(i)*abu_array(i)/(masslow(i)+w7(i))
    w72(i) = w7(i)*abu_array(i)/(masslow(i)+w7(i))
    masslow(i)*=abu_array(i)/(masslow(i)+w7(i))
    w7(i)*=abu_array(i)/(masslow(i)+w7(i))
    endfor'''

    massive2 = np.zeros(numiso)
    massive2 = m_new + masshigh
    w7 = np.zeros(numiso)
    w7 = w_new
   

    #binwrite,array=massive2,filename='massivePOST' ????? include or not?
    #binwrite,array=w7,filename='onePOST'

    #for testing purposes

    s_ws_r_n_g_mg = (sproc + wsproc + rproc + nproc + gproc + masshigh) / abu_array
    s_ws_r_n_g_mh = (sproc + wsproc + rproc + nproc + gproc + masshigh) / abu_array
    s_ws = (sproc + wsproc) / abu_array
    all = (sproc + wsproc + rproc + nproc + gproc + w7 + massive2 + novae_priGCR_nu + bbn + h_burn + gcr) / abu_array
    sn = (w7 + massive2) / abu_array
    #stop
    #END FIT MASSIVE AND TYPE IA

    #MAKE FINAL ARRAY TO SEND
    isotopes = np.stack((iso_zs, amu, sproc, wsproc, rproc, nproc, gproc, w7, massive2, bbn, gcr, novae_priGCR_nu, h_burn), axis = 1)
    resultFile = os.path.join(dirname, 'result_isotopes.csv')
    header = 'iso_Zs, amu, abu_array, sproc, wsproc, rproc, nproc, gproc, w7, massive2, bbn, gcr, novae_priGCR_nu, h_burn'
    np.savetxt(resultFile, isotopes, header=header, delimiter=",", fmt = '%s')

    return isotopes

def test():
    fp = 0.7
    dirname = os.path.dirname(__file__)
    amuFile = os.path.join(dirname, 'amu.dat')
    amu = np.loadtxt(amuFile, skiprows = 5)

    dirname = os.path.dirname(__file__)
    iso_zsFile = os.path.join(dirname, 'isotope_Zs.dat')
    iso_zs = np.loadtxt(iso_zsFile)

    dirname = os.path.dirname(__file__)
    ionFile = os.path.join(dirname, 'ions.dat')
    ion = np.loadtxt(ionFile, dtype = 'str')

    dirname = os.path.dirname(__file__)
    massiveFile = os.path.join(dirname, 'starfit_274823_6_13_2013.dat')
    massive = np.loadtxt(massiveFile)
    massive[0:9] = 0
    massive[62:287] = 0

    dirname = os.path.dirname(__file__)
    abuFile = os.path.join(dirname, 'sollo20.dat')
    abu_array = np.loadtxt(abuFile, skiprows = 2, usecols = (1))
    abu_array =  abu_array / amu

    # dirname = os.path.dirname(__file__)
    # abuFile = os.path.join(dirname, 'sollo20.dat')
    # abu_array = np.loadtxt(abuFile, skiprows = 2, usecols = (0))
    # abu_array =  abu_array / amu



    result = solar_decomp(abu_array, amu, iso_zs, fp, ion, massive)
    dirname = os.path.dirname(__file__)
    resultFile = os.path.join(dirname, 'result.csv')
    header = 'iso_Zs, amu, abu_array, sproc, wsproc, rproc, nproc, gproc, w7, massive2, bbn, gcr, novae_priGCR_nu, h_burn'
    np.savetxt(resultFile, result, header=header, delimiter=",")

if __name__=="__main__":
    test()
