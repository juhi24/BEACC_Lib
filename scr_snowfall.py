# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import read
from glob import glob
from os import path
import matplotlib.pyplot as plt

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

dt_start = '20140201T00:00:01'
dt_end = '20140228T23:45:00'

#dt_start = '20140221T16:00:01'
#dt_end = '20140221T23:45:00'

#dt_start = '20140223T00:00:01'
#dt_end = '20140223T23:00:00'

#dt_start = '20140208T0:30:01'
#dt_end = '20140208T10:30:00'

m200, m400 = Method1.from_hdf(dt_start, dt_end, unbias=True, rule='2min')

m200.dsd.data.drop(['26.0'], 1, inplace=True)

#m200.plot()
#m400.plot()
#plt.show()