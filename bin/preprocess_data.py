import numpy as np
import os
from gwdama.io import GwDataManager
import matplotlib.pyplot as plt
#import tkinter
import matplotlib
matplotlib.use('TkAgg')

rawdata_dama = GwDataManager(os.path.join('GWsimdata', 'BHBH-CE-gwskysim-strain_H1_1241833603-1241843628GPS.hdf5'), mode='r')

rawdata_strain = rawdata_dama['gwskysim_simulation']

rawdata_strain.whiten(4)

#plt.figure()
data_timeseries = rawdata_dama['gwskysim_simulation_whiten'].to_TimeSeries()
data_timeseries.plot()
print(rawdata_dama)

plt.show()
