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

def gray_out(data, axdict, fltr_label='d0_fltr', labels=['D_0', 'density']):
    for t_end, row in data[data[fltr_label]].iterrows():
        for label in labels:
            ax = axdict[label]
            ax.axvspan(row.start, t_end, edgecolor='none', facecolor='0.8',
                       alpha=0.8)

def plot_overview(data, params=['intensity', 'density', 'D_0', 'N_w'],
                  axlist=None):
    data.density[data.density>800] = np.nan
    if axlist is None:
        axlist = data.loc[:, params].plot(figsize=(5,7), subplots=True,
                                         drawstyle='steps')
    else:
        for i, param in enumerate(params):
            ax = axlist[i]
            data[param].plot(ax=ax, drawstyle='steps')
    axdict = dict(zip(params, axlist))
    fig = axlist[0].get_figure()
    #axdict['density'].axis((None, None, 0, 600))
    for param in ['N_w']:
        axdict[param].set_yscale('log')
    data = read.merge_series(data, case.instr['pluvio'].start_time())
    data['d0_fltr'] = data.D_0 < 0.63
    gray_out(data, axdict)
    return fig, axlist

for case in e.events.paper.values:
    data = case.summary()
    fig, axlist = plot_overview(data)
    plt.legend()
    plt.tight_layout()
    fig.savefig(path.join(savedir, case.dtstr('%Y%m%d') + '.eps'), dpi=150)
