# coding: utf-8
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob
from baecc.instruments import pip_psd

t0=pd.datetime(2014,3,3,5,0)
t1=pd.datetime(2014,3,3,9,0)
fnames=glob('/home/jussitii/DATA/PIP/a_DSD_Tables/00420140303*.dat')
psd = pip_psd.PipPSD(filenames=fnames, dt_start=t0, dt_end=t1)
psd.store_good_data()
data = psd.good_data()


# plotting
#plt.figure()
#psd.plot()

