# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import read
from glob import glob
import copy
import matplotlib.pyplot as plt

#dtstr = '20140[2-3]??' #    all feb-mar
#dtstr = '20140212' #        feb 12
dtstr = '20140221' #        feb 21

# feb 21
pipv_files = glob('/home/jussitii/DATA/PIP/a_Velocity_Tables/004%s/*2.dat' % dtstr)
dsd_files = glob('/home/jussitii/DATA/PIP/a_DSD_Tables/004%s_a_d.dat' % dtstr)
pluvio200_files = glob('/home/jussitii/DATA/Pluvio200/pluvio200_0?_%s*' % dtstr)
pluvio400_files = glob('/home/jussitii/DATA/Pluvio400/pluvio400_0?_%s*' % dtstr)

pluvio200 = read.Pluvio(pluvio200_files)
pluvio400 = read.Pluvio(pluvio400_files)
pipv = read.PipV(pipv_files)
dsd = read.PipDSD(dsd_files)

pluvio200_lim = copy.deepcopy(pluvio200)
pluvio400_lim = copy.deepcopy(pluvio400)

dt_start = '20140221T16:00:01'
dt_end = '20140221T23:45:00'

#dt_start = '20140215T21:20:01'
#dt_end = '20140216T1:00:00'

#dt_start = '20140208T0:30:01'
#dt_end = '20140208T10:30:00'

dsd_lim = dsd
for instr_lim in [pluvio200_lim,pluvio400_lim,dsd_lim,pipv]:
    instr_lim.data = instr_lim.data[dt_start:dt_end]

m200 = Method1(dsd_lim,pipv,pluvio200_lim,rule='5min')
m400 = Method1(dsd_lim,pipv,pluvio400_lim,rule='5min')

m200.plot()
#m400.plot()
#m200_s.plot()
#m400_s.plot()
plt.show()