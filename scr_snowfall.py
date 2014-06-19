# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import read
from glob import glob
import copy
import matplotlib.pyplot as plt

# All from feb-march
#pipv_files = glob('/home/jussitii/DATA/PIP/a_Velocity_Tables/00420140[2-3]??/*2.dat')
#dsd_files = glob('/home/jussitii/DATA/PIP/a_DSD_Tables/00420140[2-3]??_a_d.dat')
#pluvio200_files = glob('/home/jussitii/DATA/Pluvio200/pluvio200_0?_20140[2-3]*')
#pluvio400_files = glob('/home/jussitii/DATA/Pluvio400/pluvio400_0?_20140[2-3]*')
        
#pipv_files = glob('/home/jussitii/DATA/PIP/a_Velocity_Tables/00420140212/*2.dat')
#dsd_files = glob('/home/jussitii/DATA/PIP/a_DSD_Tables/00420140212_a_d.dat')
#pluvio200_files = glob('/home/jussitii/DATA/Pluvio200/pluvio200_0?_20140212*')
#pluvio400_files = glob('/home/jussitii/DATA/Pluvio400/pluvio400_0?_20140212*')

# feb 21
pipv_files = glob('/home/jussitii/DATA/PIP/a_Velocity_Tables/00420140221/*2.dat')
dsd_files = glob('/home/jussitii/DATA/PIP/a_DSD_Tables/00420140221_a_d.dat')
pluvio200_files = glob('/home/jussitii/DATA/Pluvio200/pluvio200_0?_20140221*')
pluvio400_files = glob('/home/jussitii/DATA/Pluvio400/pluvio400_0?_20140221*')

# feb 15-16
#dsd_files = glob('/home/jussitii/DATA/PIP/a_DSD_Tables/0042014021[5-6]_a_d.dat')
#pipv_files = glob('/home/jussitii/DATA/PIP/a_Velocity_Tables/0042014021[5-6]/*2.dat')
#pluvio400_files = glob('/home/jussitii/DATA/Pluvio400/pluvio400_0?_2014021[5-6]*')
#pluvio200_files = glob('/home/jussitii/DATA/Pluvio200/pluvio200_0?_2014021[5-6]*')

#pipv_files = glob('/home/jussitii/DATA/PIP/a_Velocity_Tables/00420140208/*2.dat')
#dsd_files = glob('/home/jussitii/DATA/PIP/a_DSD_Tables/00420140208_a_d.dat')
#pluvio200_files = glob('/home/jussitii/DATA/Pluvio200/pluvio200_0?_20140208*')
#pluvio400_files = glob('/home/jussitii/DATA/Pluvio400/pluvio400_0?_20140208*')

pluvio200 = read.Pluvio(pluvio200_files)
pluvio400 = read.Pluvio(pluvio400_files)
pipv = read.PipV(pipv_files)
dsd = read.PipDSD(dsd_files)

pluvio200_lim = copy.deepcopy(pluvio200)
pluvio400_lim = copy.deepcopy(pluvio400)

#shift
pluvio200.data = pluvio200.data.shift(-5).fillna(0)
pluvio400.data = pluvio400.data.shift(-5).fillna(0)

dt_start = '20140221T16:00:01'
dt_end = '20140221T23:45:00'

#dt_start = '20140215T21:20:01'
#dt_end = '20140216T1:00:00'

#dt_start = '20140208T0:30:01'
#dt_end = '20140208T10:30:00'

pluvio200_lim_shift = copy.deepcopy(pluvio200)
pluvio400_lim_shift = copy.deepcopy(pluvio400)
dsd_lim = dsd
for instr_lim in [pluvio200_lim,pluvio400_lim,dsd_lim,pipv]:
    instr_lim.data = instr_lim.data[dt_start:dt_end]

m200 = Method1(dsd_lim,pipv,pluvio200_lim,rule='5min')
m400 = Method1(dsd_lim,pipv,pluvio400_lim,rule='5min')

m200_s = Method1(dsd_lim,pipv,pluvio200_lim_shift,rule='5min')
m400_s = Method1(dsd_lim,pipv,pluvio400_lim_shift,rule='5min')

m200.plot()
m400.plot()
m200_s.plot()
m400_s.plot()
plt.show()