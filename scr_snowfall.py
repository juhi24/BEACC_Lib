# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *

e = EventsCollection('cases_of_interest.csv', '%d.%m. %H:%M')
e.autoimport_data(autoshift=False)
    
#dt_start = pd.datetime(2014, 2, 1, 0, 0, 1)
#dt_end = pd.datetime(2014, 7, 31, 23, 40, 0)

#m200, m400 = Case.from_hdf(dt_start, dt_end, autoshift=False, rule='6min')

#instr = batch_import(dtstr='2014022[1-2]', datadir='../DATA')
#m200 = Case(instr['dsd'], instr['vel'], instr['pluvio200'], rule='5min',
#               liquid=False)

#m200.dsd.data.drop([26.0], 1, inplace=True)
