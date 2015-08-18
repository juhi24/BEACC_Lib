# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
from datetime import datetime
from os import path

dtformat_default = '%d.%m. %H:%M'
dtformat_snex = '%Y %d %B %H UTC'

e = EventsCollection('cases/pip2015.csv', dtformat_snex)
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.instr['pluvio'].shift_periods = -6
    c.instr['pluvio'].n_combined_intervals = 2
    c.instr['pipv'].use_flip = False
    basename = read.ensure_dir(path.join('../results/pip2015/rho',  c.instr['pluvio'].name))
    c.density(pip_filter=False).to_csv(path.join(basename, c.dtstr('%Y%m%d') + '.csv'))
    c.instr['pluvio'].tdelta().to_csv(basename + 'timedelta_' + c.dtstr('%Y%m%d') + '.csv')
