# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import read
from glob import glob
import copy

# All from feb-march
#pipv_files = glob('/home/jussitii/DATA/PIP/a_Velocity_Tables/00420140[2-3]??/*2.dat')
#dsd_files = glob('/home/jussitii/DATA/PIP/a_DSD_Tables/00420140[2-3]??_a_d.dat')
#pluvio200_files = glob('/home/jussitii/DATA/Pluvio200/pluvio200_0?_20140[2-3]*')
#pluvio400_files = glob('/home/jussitii/DATA/Pluvio400/pluvio400_0?_20140[2-3]*')
        
pipv_files = glob('/home/jussitii/DATA/PIP/a_Velocity_Tables/00420140212/*2.dat')
dsd_files = glob('/home/jussitii/DATA/PIP/a_DSD_Tables/00420140212_a_d.dat')
pluvio200_files = glob('/home/jussitii/DATA/Pluvio200/pluvio200_0?_20140212*')
pluvio400_files = glob('/home/jussitii/DATA/Pluvio400/pluvio400_0?_20140212*')

#dsd_files = glob('/home/jussitii/DATA/PIP/a_DSD_Tables/0042014021[5-6]_a_d.dat')
#pipv_files = glob('/home/jussitii/DATA/PIP/a_Velocity_Tables/0042014021[5-6]/*2.dat')
#pluvio400_files = glob('/home/jussitii/DATA/Pluvio400/pluvio400_0?_2014021[5-6]*')
#pluvio200_files = glob('/home/jussitii/DATA/Pluvio200/pluvio200_0?_2014021[5-6]*')

pluvio200 = read.Pluvio(pluvio200_files)
pluvio400 = read.Pluvio(pluvio400_files)
pipv = read.PipV(pipv_files)
dsd = read.PipDSD(dsd_files)

dt_start = '20140212T4:30:01'
dt_end = '20140212T8:30:00'
pluvio200_lim = copy.deepcopy(pluvio200)
pluvio400_lim = copy.deepcopy(pluvio400)
dsd_lim = dsd
for instr_lim in [pluvio200_lim,pluvio400_lim,dsd_lim,pipv]:
    instr_lim.data = instr_lim.data[dt_start:dt_end]

m200 = Method1(dsd_lim,pipv,pluvio200_lim,rule='15min')
m400 = Method1(dsd_lim,pipv,pluvio400_lim,rule='15min')