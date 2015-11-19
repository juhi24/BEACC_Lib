# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import matplotlib.pyplot as plt
from os import path
import read

from scr_snowfall import pip2015events, test_events

plt.close('all')
plt.ion()
debug = True
include_vap = False
savedir = '../results/pip2015/case_overview'

if debug:
    e = test_events()
    savedir += '/test'
else:
    e = pip2015events()

read.ensure_dir(savedir)

if include_vap:
    df_path = '/home/jussitii/DATA/arm/tmpmwrlosM1.b1.20140201.000113.cdf'
    vap = read.cdf_to_series(cdf_path, 'vap')

for case in e.events.paper.values:
    data = case.summary()
    if include_vap:
        vap = vap[data.index[0]:data.index[-1]].copy()
    
    params = ['intensity', 'density', 'D_0', 'N_w']
    data.density[data.density>800] = np.nan
    if include_vap:
        params.append('')
    axarr = data.loc[:, params].plot(figsize=(5,7), subplots=True,
                                     drawstyle='steps')
    axdict = dict(zip(params, axarr))
    fig = axarr[0].get_figure()
    if include_vap:
        vap.plot(ax=axarr[-1], label='LWP')
    
    #axdict['density'].axis((None, None, 0, 600))
    
    for param in ['N_w']:
        axdict[param].set_yscale('log')
    
    plt.legend()
    plt.tight_layout()
    fig.savefig(path.join(savedir, case.dtstr('%Y%m%d') + '.eps'), dpi=150)