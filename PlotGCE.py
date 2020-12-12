import numpy as np
import os
import pandas as pd
from gce_data import gce_data
import matplotlib.pyplot as plt

fig = plt.figure()
data = gce_data(1)
df = pd.DataFrame(data, columns = ['[Fe/H]', '[Mg/Fe]', 'feh_sig', 'mgfe_sig'])
fig, ax = plt.subplots(2,1)
df.plot.scatter(x = '[Fe/H]', y = '[Mg/Fe]', ax = ax )
plt.show(block=True)