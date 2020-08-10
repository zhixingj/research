
import numpy as np
import os

def weak_s_index(uniq,ion):

#returns array of indicies of ion array that correspond to weak s path isotopes

    nums = uniq.size
    s_index=np.zeros(nums)
    k=0
    for i in range(0,nums):
        index = np.asarray(np.where(uniq[i] == ion))[0]
        if index.size == 0:
            s_index[k] = 0
        if index.size > 0:
            s_index[k] = index[0]
        k = k +1
    # s_index = s_index[0:nums-2]
    s_index = s_index.astype(int)
    return s_index
