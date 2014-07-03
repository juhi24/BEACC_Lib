# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import read
from glob import glob
from os import path
import matplotlib.pyplot as plt
import pandas as pd

def batch_hdf(datadir='../DATA', outname='baecc.h5', dtstr='20140[2-3]??'):
    pipv_files = glob(path.join(datadir, 'PIP/a_Velocity_Tables/004%s/*2.dat' % dtstr))
    dsd_files = glob(path.join(datadir, 'PIP/a_DSD_Tables/004%s_a_d.dat' % dtstr))
    pluvio200_files = glob(path.join(datadir, 'Pluvio200/pluvio200_0?_%s*' % dtstr))
    pluvio400_files = glob(path.join(datadir, 'Pluvio400/pluvio400_0?_%s*' % dtstr))

    pluvio200 = read.Pluvio(pluvio200_files)
    pluvio400 = read.Pluvio(pluvio400_files)
    pipv = read.PipV(pipv_files)
    dsd = read.PipDSD(dsd_files)

    hdf_file = glob(path.join(datadir, outname))

    for instr in [pluvio200, pluvio400, pipv, dsd]:
        instr.to_hdf(hdf_file)

dt_start = pd.datetime(2014, 2, 1, 0, 0, 1)
dt_end = pd.datetime(2014, 2, 28, 23, 45, 0)

m200, m400 = Method1.from_hdf(dt_start, dt_end, autoshift=True, rule='2min')

m200.dsd.data.drop(['26.0'], 1, inplace=True)

case_start = pd.datetime(2014, 2, 2, 16, 0, 1)
case_end = pd.datetime(2014, 2, 2, 18, 0, 0)

#case_start = pd.datetime(2014, 2, 8, 0, 30, 1)
#case_end = pd.datetime(2014, 2, 8, 10, 30, 1)

#case_start = pd.datetime(2014, 2, 12, 0, 0, 1)
#case_end = pd.datetime(2014, 2, 12, 23, 30, 0)

#case_start = pd.datetime(2014, 2, 15, 21, 0, 1)
#case_end = pd.datetime(2014, 2, 16, 1, 0, 0)

#case_start = pd.datetime(2014, 2, 21, 22, 0, 1)
#case_end = pd.datetime(2014, 2, 21, 23, 0, 0)

#case_start = pd.datetime(2014, 2, 23, 0, 0, 1)
#case_end = pd.datetime(2014, 2, 23, 23, 0, 0)

#case_start = pd.datetime(2014, 3, 2, 10, 0, 1)
#case_end = pd.datetime(2014, 3, 2, 11, 0, 0)

#case_start = pd.datetime(2014, 3, 20, 19, 0, 1)
#case_end = pd.datetime(2014, 3, 20, 21, 0, 0)

m = m200.between_datetime(case_start, case_end)
m.noprecip_bias()

#m200.plot()
#m400.plot()
#plt.show()