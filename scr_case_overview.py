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
include_vap = False
savedir = '../results/pip2015/case_overview'

read.ensure_dif(savedir)

if debug:
    e = test_events()
else:
    e = pip2015events()

ïf include_vap:
    df_path = '/home/jussitii/DATA/arm/tmpmwrlosM1.b1.20140201.000113.cdf'
    vap = read.cdf_to_series(cdf_path, 'vap')

for case in e.event.paper.values:
    data = case.summary()
    ïf include_vap:
        vap = vap[data.index[0]:data.index[-1]].copy()
    
    params = ['intensity', 'density', 'D_0', 'N_w']
    ïf include_vap:
        params.append('')
    fig = plt.figure(dpi=150, figsize=(5,10))
    axarr = data.loc[:, params].plot(subplots=True, drawstyle='steps')
    axdict = dict(zip(params, axarr))
    fig = axarr[0].get_figure()
    ïf include_vap:
        vap.plot(ax=axarr[-1], label='LWP')
    
    
    for param in ['N_w']:
        axdict[param].set_yscale('log')
    
    plt.legend()
    plt.tight_layout()
    fig.savefig(path.join(savedir, case.dtstr('%Y%m%d') + '.eps'))