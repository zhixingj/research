import numpy as np
import os

def weak_s_target(ion):
    # ;returns array of weak s ion target nuclei
    #
    # ;sion=strarr(nums+1)
    # ;sion[0:7]=s_name1[0:7]
    # ;sion[8]=s_name2[8]
    # ;sion[9]=s_name1[8]
    # ;sion[10:nums]=s_name1[9:nums-1]

    # ;a=strarr(n_elements(s_name1)+n_elements(s_name2))
    # ;a[0:n_elements(s_name1)-1]=s_name1
    # ;a[n_elements(s_name1):n_elements(s_name1)+n_elements(s_name2)-1]=s_name2
    # ;unique_index=uniq(a)

    sion=['Fe56',\
    'Fe57',\
    'Fe58',\
    'Co59',\
    'Ni60',\
    'Ni61',\
    'Ni62',\
    'Cu63',\
    'Cu63',\
    'Ni64',\
    'Zn64',\
    'Cu65',\
    'Zn66',\
    'Zn67',\
    'Zn68',\
    'Ga69',\
    'Ge70',\
    'Ga71',\
    'Ge72',\
    'Ge73',\
    'Ge74',\
    'As75',\
    'Se76',\
    'Se77',\
    'Se78',\
    'Br79',\
    'Br79',\
    'Se80',\
    'Kr80',\
    'Br81',\
    'Kr82',\
    'Kr83',\
    'Kr84',\
    'Rb85',\
    'Sr86',\
    'Sr87',\
    'Sr88',\
    'Y89',\
    'Zr90',\
    'Zr91',\
    'Zr92',\
    'Zr93',\
    'Zr94',\
    'Mo95',\
    'Mo96',\
    'Mo97',\
    'Mo98',\
    'Tc99',\
    'Ru100',\
    'Ru101']

    # ;THIS CODE TAKES INTO ACCOUNT n\b- BRANCHINGS
    # ;nums=51
    # ;sion=strarr(nums)
    # ;sion[0:3]=ion[61:64]
    # ;sion[4:20]=ion[66:82]
    # ;sion[21]=ion[84]
    # ;sion[22:25]=ion[86:89]
    # ;sion[26:27]=ion[91:92]
    # ;sion[28:31]=ion[94:97]
    # ;sion[32]=ion[99]
    # ;sion[33:40]=ion[102:109]
    # ;sion[41]=ion[111]
    # ;sion[42:46]=ion[113:117]
    # ;sion[47:50]=ion[121:124]

    return sion
