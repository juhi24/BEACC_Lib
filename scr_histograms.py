# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
from os import path
import seaborn as sns

dtformat_default = '%d.%m. %H:%M'
dtformat_snex = '%Y %d %B %H UTC'

e = EventsCollection('cases/pip2015.csv', dtformat_snex)
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

plt.close('all')
#plt.ioff()
plt.ion()
kwargs = {'kde': True, 'rug': True, 'kde_kws': {'label': 'KDE'}}

def subplots():
    return plt.subplots(e.events.pluvio400.count(), sharex=True, sharey=True,
                        tight_layout=False)

def remove_gaps(f):
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

fd, axarrd = subplots()
fm, axarrm = subplots()
fn, axarrn = subplots()

for i, c in enumerate(e.events.pluvio400.values):
    c.instr['pluvio'].shift_periods = -6
    c.instr['pluvio'].n_combined_intervals = 2
    rng = (0,6)
    sns.distplot(c.d_0().dropna(), ax=axarrd[i], label=c.dtstr(), bins=17, 
                 hist_kws={'range':rng}, **kwargs)
    axarrd[i].set_xlim(rng)
    rng = (-2, 8)
    axm = sns.distplot(c.mu().dropna(), ax=axarrm[i], label=c.dtstr(), bins=20,
                       hist_kws={'range':rng}, **kwargs)
    axarrm[i].set_xlim(rng)
    rng = (0, 200000)
    axn = sns.distplot(c.n_w().dropna(), ax=axarrn[i], label=c.dtstr(), bins=30,
                       hist_kws={'range':rng}, **kwargs)
    axarrn[i].set_xlim(rng)

for f in (fd, fm, fn):
    remove_gaps(f)