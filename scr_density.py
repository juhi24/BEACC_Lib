# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

dtformat_default = '%d.%m. %H:%M'
dtformat_default_year = '%d.%m.%Y %H:%M'
dtformat_snex = '%d %B %H UTC'

e = EventsCollection('cases/test.csv', dtformat_default_year)
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True, datafile=['../DATA/test03.h5'])

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.pluvio.shift_periods = -6
    basename = '/home/dori/SnowCases_BAEC/DensityJussi/to_davide/' + datetime.strftime(c.pluvio.dt_start().to_datetime(),'%Y%m%d%H')
    c.density(pip_filter=False).to_csv(basename + 'density_' + c.pluvio.name + '.csv')
    c.pluvio.tdelta().to_csv(basename + 'timedelta_' + c.pluvio.name + '.csv')

#c.plot_velfitcoefs(rhomax=600, countmin=2000)
