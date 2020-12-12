import numpy as np

def light(lightest, z, z_max, numiso, sp, rp):
    # ;BBN values
    # ;Mass frac
    iso_lt = np.zeros((z_max, numiso))
    h1=0.75128
    h2=4.31998e-5
    he3=2.13003e-5
    he4=.248659
    li6=5.96248E-14  
    li7=1.912209E-09 

    # mol/g
    h1=0.75128/1.007825
    h2=4.31998e-5/2.014102
    he3=2.13003e-5/3.016029
    he4=.248659/4.002603
    li7 = 1.912209E-09/7.016004

    hydrogen = np.zeros((z_max,1))

    # ;HELIUM -> ASSUMPTION IS THAT He ~ primary
    for i in range(0, z_max):
        iso_lt[i, 2] = (lightest[2]-he3)/(1.0-0.0)*z[i]**(rp)+he3
        iso_lt[i, 3] = (lightest[3]-he4)/(1.0-0.0)*z[i]**(rp)+he4

        # ;DEUTERIUM -> ASSUMPTION IS THAT D ~ Z
        iso_lt[i, 1] = (lightest[1]-h2)/(1.0-0.0)*z[i]**(rp)+h2

        # ;HYDROGEN -> Use H1 = 1.0 - Y - Z - D
        iso_lt[i, 0] = 1.0 - z[i]*0.01530 - (iso_lt[i, 2]*3.016029 + iso_lt[i, 3]*4.002603) - iso_lt[i, 1]*2.014102

        ''';;HYDROGEN -> H = 1.057 - Z - He (in mass fractions)(OUTDATED)
        ;;H_1/H_2 = 51,550 AT SOLAR, PERSERVE THIS RATIO FOR THE TIME BEING (OUTDATED)
        ;hydrogen(i) = 1.d0 - z(i)*0.01530 - (iso_lt(2,i)*3.016029 + iso_lt(3,i)*4.002603)
        ;iso_lt(0,i) = hydrogen(i)*(0.999986)
        ;iso_lt(1,i) = hydrogen(i)*(1.0 - 0.999986)'''

        # Li6
        iso_lt[i, 4] = lightest[4]*z[i]**sp

        # Li7
        iso_lt[i, 5] = (lightest[5]-li7)*z[i]**(rp)+li7

        # Be9
        iso_lt[i, 6] = lightest[6]*z[i]**(rp)

        # B10 and B11
        iso_lt[i, 7] = lightest[7]*z[i]**sp
        iso_lt[i, 8] = 0.6 *lightest[8]*z[i]**sp+0.4 *lightest[8]*z[i]**rp


        # zero out all other columns
        iso_lt[i, 9:286] = 0.0
    
    return iso_lt
