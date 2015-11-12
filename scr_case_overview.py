# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import matplotlib.pyplot as plt
from os import path
import read

from scr_snowfall import pip2015events, test_events

debug = False

if debug:
    e = test_events()
else:
    e = pip2015events()

cdf_path = '/home/jussitii/DATA/arm/tmpmwrlosM1.b1.20140201.000113.cdf'
vap = read.cdf_to_series(cdf_path, 'vap')

case = e.events.paper[0]
data = case.summary()
vap = vap[data.index[0]:data.index[-1]].copy()

params = ['intensity', 'density', 'D_0', 'N_0', 'N_w', '']
axarr = data.loc[:, params].plot(subplots=True, drawstyle='steps')
axdict = dict(zip(params, axarr))
fig = axarr[0].get_figure()
vap.plot(ax=axarr[-1], label='LWP')


for param in ['N_0', 'N_w']:
    axdict[param].set_yscale('log')

plt.legend()
plt.tight_layout()
