# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np

e = EventsCollection('cases/evap.csv', '%d.%m. %H:%M')
e.autoimport_data(autoshift=False, rule='5min')

p4 = e.events.pluvio400[0].pluvio
p4.shift_periods = -6
p4.acc('1min').tshift(periods=6).plot()
raw = p4.data.bucket_nrt - p4.data.bucket_nrt[0]
raw.plot()

#for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    #c.pluvio.acc('1min').plot()    
    #c.pluvio.shift_periods = -6
    #c.plot(pip=False)
    #c.pipv.plots(save=True)

#for c in e.events.pluvio200:
#    c.minimize_lsq()
    
#dt_start = pd.datetime(2014, 2, 1, 0, 0, 1)
#dt_end = pd.datetime(2014, 7, 31, 23, 40, 0)

#m200, m400 = Case.from_hdf(dt_start, dt_end, autoshift=False, rule='6min')

#instr = batch_import(dtstr='2014022[1-2]', datadir='../DATA')
#m200 = Case(instr['dsd'], instr['vel'], instr['pluvio200'], rule='5min',
#               liquid=False)

#m200.dsd.data.drop([26.0], 1, inplace=True)
