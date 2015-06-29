# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
from read import *
import numpy as np
import matplotlib.pyplot as plt
from os import path

dtformat_default = '%d.%m. %H:%M'
dtformat_snex = '%Y %d %B %H UTC'

e = EventsCollection('cases/pip2015.csv', dtformat_snex)
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

#plt.close('all')
#plt.ioff()

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.pluvio.shift_periods = -6
    c.pluvio.n_combined_intervals = 2
    ax = c.plot_vfits_in_density_ranges()
    savepath = read.ensure_dir(path.join('../results/pip2015/vfits_density_ranges', c.pluvio.name))
    plt.savefig(path.join(savepath, c.dtstr('%Y%m%d%H%M.eps')))