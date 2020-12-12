import numpy as np
import os

def branch(sion):
    # ;returns array of barnching ratia for weak s path
    # ;array element is defined to be 1d0 unless there is a branching
    # ;branching ratios calculated from: Takahashi, K., Yokoi, K., At. Data Nucl. Data Tables 36, 375 (1987)
    #
    # ;branchings (defined by b+)
    #

    b= np.zeros(sion.size)
    b[:]=1

    Cu63= np.where(sion == 'cu63')

    b[Cu63[0][0]]=0.3311
    b[Cu63[0][1]]=1-b[Cu63[0][0]]


    Br79=np.where(sion == 'br79')
    b[Br79[0][0]]=0.0325
    b[Br79[0][1]]=1-b[Br79[0][0]]


    return b



    # return b
