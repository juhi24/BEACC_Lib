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

def batch_import(dtstr, datadir='../DATA'):
    pipv_files = glob(path.join(datadir, 'PIP/a_Velocity_Tables/004%s/*2.dat' % dtstr))
    dsd_files = glob(path.join(datadir, 'PIP/a_DSD_Tables/004%s_a_d.dat' % dtstr))
    pluvio200_files = glob(path.join(datadir, 'Pluvio200/pluvio200_??_%s*.txt' % dtstr))
    pluvio400_files = glob(path.join(datadir, 'Pluvio400/pluvio400_??_%s*.txt' % dtstr))
    
    pluvio200 = read.Pluvio(pluvio200_files)
    pluvio400 = read.Pluvio(pluvio400_files)
    pipv = read.PipV(pipv_files)
    dsd = read.PipDSD(dsd_files)
    return {'vel': pipv, 'dsd': dsd, 
            'pluvio200': pluvio200, 'pluvio400': pluvio400}

def batch_hdf(datadir='../DATA', outname='baecc.h5', dtstr='20140[2-3]??'):
    instrdict = batch_import(dtstr, datadir)

    hdf_file = glob(path.join(datadir, outname))

    for instr in [instrdict.values()]:
        instr.to_hdf(hdf_file)

dt_start = pd.datetime(2014, 2, 1, 0, 0, 1)
dt_end = pd.datetime(2014, 3, 1, 23, 45, 0)

#m200, m400 = Method1.from_hdf(dt_start, dt_end, autoshift=False, rule='5min')
instr = batch_import(dtstr='201402??', datadir='../DATA')
m200 = Method1(instr['dsd'], instr['vel'], instr['pluvio200'], rule='5min')

m200.dsd.data.drop([26.0], 1, inplace=True)

case_start = pd.datetime(2014, 2, 2, 16, 0, 1)
case_end = pd.datetime(2014, 2, 2, 18, 0, 0)
case2 = [case_start, case_end]

case_start = pd.datetime(2014, 5, 9, 16, 0, 1)
case_end = pd.datetime(2014, 5, 9, 21, 0, 0)
raincase9 = [case_start, case_end]

case_start = pd.datetime(2014, 5, 26, 17, 30, 1)
case_end = pd.datetime(2014, 5, 26, 19, 30, 0)
raincase26 = [case_start, case_end]

case7_start = pd.datetime(2014, 2, 7, 22, 45, 1)
case7_end = pd.datetime(2014, 2, 8, 0, 0, 0)
case7 = [case7_start, case7_end]

case_start = pd.datetime(2014, 2, 8, 0, 30, 1)
case_end = pd.datetime(2014, 2, 8, 10, 30, 0)

case_start = pd.datetime(2014, 2, 12, 0, 0, 1)
case_end = pd.datetime(2014, 2, 12, 23, 30, 0)
case12 = [case_start, case_end]

case16_start = pd.datetime(2014, 2, 16, 0, 0, 1)
case16_end = pd.datetime(2014, 2, 16, 1, 0, 0)
case16 = [case16_start, case16_end]

case21_start = pd.datetime(2014, 2, 21, 22, 50, 1)
case21_end = pd.datetime(2014, 2, 21, 23, 20, 0)
case21 = [case21_start, case21_end]

case_start = pd.datetime(2014, 2, 23, 0, 0, 1)
case_end = pd.datetime(2014, 2, 23, 23, 0, 0)
case23 = [case_start, case_end]

case_start = pd.datetime(2014, 3, 2, 10, 0, 1)
case_end = pd.datetime(2014, 3, 2, 11, 0, 0)
mar2 = [case_start, case_end]

case_start = pd.datetime(2014, 3, 18, 9, 0, 1)
case_end = pd.datetime(2014, 3, 18, 10, 0, 0)
mar18 = [case_start, case_end]

case_start = pd.datetime(2014, 3, 20, 19, 0, 1)
case_end = pd.datetime(2014, 3, 20, 21, 0, 0)
mar20 = [case_start, case_end]

case = [] # initialize
for case_span in [case7, case16, case21]:
    m = m200.between_datetime(*case_span)
    m.autoshift(inplace=True)
    m.noprecip_bias(inplace=True)
    #m.plot()
    #plt.figure()
    #m.plot_v_stuff()
    case.append(m)

#plt.figure()
#case[0].plot_v_stuff()

#m400.plot()
#plt.show()