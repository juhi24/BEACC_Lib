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
debug = False
include_vap = False
savedir = '../results/pip2015/case_overview'

if debug:
    e = test_events()
    savedir += '/test'
else:
    e = pip2015events()

read.ensure_dir(savedir)

def gray_out(data, axdict, fltr_label='d0_fltr', labels=['D_0', 'density']):
    for t_end, row in data[data[fltr_label]].iterrows():
        for label in labels:
            ax = axdict[label]
            ax.axvspan(row.start, t_end, edgecolor='none', facecolor='0.8',
                       alpha=0.8)

if include_vap:
    df_path = '/home/jussitii/DATA/arm/tmpmwrlosM1.b1.20140201.000113.cdf'
    vap = read.cdf_to_series(cdf_path, 'vap')

for case in e.events.paper.values:
    data = case.summary()
    if include_vap:
        vap = vap[data.index[0]:data.index[-1]].copy()
    
    params = ['intensity', 'density', 'D_0', 'N_w']
    start_time = case.instr['pluvio'].start_time()
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
    data = read.merge_series(data, case.instr['pluvio'].start_time())
    data['d0_fltr'] = data.D_0 < 0.63
    gray_out(data, axdict)
    plt.legend()
    plt.tight_layout()
    fig.savefig(path.join(savedir, case.dtstr('%Y%m%d') + '.eps'), dpi=150)