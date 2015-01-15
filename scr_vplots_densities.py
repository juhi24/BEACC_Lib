# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
import os

dtformat_default = '%d.%m. %H:%M'
dtformat_snex = '%d %B %H UTC'

e = EventsCollection('cases/pip2015.csv', dtformat_snex)
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

basedir = '/home/jussitii/results/pip2015'

plt.ioff()

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.pluvio.shift_periods = -6
    dtstr = datetime.strftime(c.pluvio.dt_start().to_datetime(),'%Y%m%d%H')
    density_filepath = os.path.join(basedir, dtstr + 'density_' + c.pluvio.name + '.csv')
    delta_filepath = os.path.join(basedir, dtstr + 'timedelta_' + c.pluvio.name + '.csv')
    c.density(pip_filter=False).to_csv(density_filepath)
    c.pluvio.tdelta().to_csv(delta_filepath)
    savedir = os.path.join(basedir, 'velplots', dtstr[0:8], c.pluvio.name)
    ensure_dir(savedir)
    c.pipv.plots(save=True, savedir=savedir, suffix='.eps', rule=c.rule, ymax=3)
    plt.close('all')

#c.plot_velfitcoefs(rhomax=600, countmin=2000)
