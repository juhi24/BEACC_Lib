# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import read
from glob import glob
from datetime import datetime
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
    hdf_file = path.join(datadir, outname)
    for key in instrdict:
        instrdict[key].to_hdf(filename=hdf_file)
    
dt_start = pd.datetime(2014, 2, 1, 0, 0, 1)
dt_end = pd.datetime(2014, 7, 31, 23, 40, 0)

m200, m400 = Case.from_hdf(dt_start, dt_end, autoshift=False, rule='6min')
#instr = batch_import(dtstr='2014022[1-2]', datadir='../DATA')
#m200 = Case(instr['dsd'], instr['vel'], instr['pluvio200'], rule='5min',
#               liquid=False)

#m200.dsd.data.drop([26.0], 1, inplace=True)

#c200 = []
#c400 = []
#casetimes = [case7, case16, case21]
#for case in casetimes:
#    c200.append(m200.between_datetime(*case, autoshift=True, autobias=True))
#    c400.append(m400.between_datetime(*case, autoshift=True, autobias=True))
