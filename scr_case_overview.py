# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.dates import DateFormatter
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

def d0fltr(data, case):
    data = read.merge_series(data, case.instr['pluvio'].start_time())
    data['d0_fltr'] = data.D_0 < 0.63
    return data

def plot_overview(data, params=['intensity', 'density', 'D_0', 'N_w'],
                  axlist=None):
    data.density[data.density>800] = np.nan
    if axlist is None:
        axlist = data.loc[:, params].plot(figsize=(5, 7), subplots=True,
                                         drawstyle='steps')
    else:
        for i, param in enumerate(params):
            ax = axlist[i]
            data[param].plot(ax=ax, drawstyle='steps', label='$'+param+'$')
    axdict = dict(zip(params, axlist))
    fig = axlist[0].get_figure()
    #axdict['density'].axis((None, None, 0, 600))
    for param in ['N_w']:
        axdict[param].set_yscale('log')
    data = d0fltr(data, case)
    gray_out(data, axdict)
    return fig, axdict

params=['intensity', 'density', 'D_0', 'N_w']

for case in e.events.paper.values:
    data = case.summary()
    fig, axlist = plot_overview(data, params=params)
    plt.legend()
    plt.tight_layout()
    fig.savefig(path.join(savedir, case.dtstr('%Y%m%d') + '.eps'), dpi=150)
fig = plt.figure(figsize=(6, 9))
extent = (0.375, 5, 0.5, 2.5)
dtlist = ['2015-01-14 02:53:00', '2015-01-14 03:30:00', '2015-01-14 03:40:00']
series_ax = []
fit_ax = []
gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
gs_series = gridspec.GridSpecFromSubplotSpec(len(params), 1, subplot_spec=gs[0],
                                             hspace=0.15)
gs_fit = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1],
                                          wspace=0.05)
for i in range(len(params)):
    series_ax.append(plt.subplot(gs_series[i]))
for i, dt in enumerate(dtlist):
    vfit = case.instr['pipv'].fits.polfit[dt]
    ax = plt.subplot(gs_fit[i])
    vfit.plot(source_style='hex', unfiltered=True, ax=ax,
              source_kws={'gridsize': 40, 'extent': extent})
    ax.axis(extent)
    fit_ax.append(ax)
fig, axdict = plot_overview(data, axlist=series_ax, params=params)
data['D_max'].plot(ax=axdict['D_0'], drawstyle='steps', label='$D_{max}$')
axdict['D_0'].legend()
for ax in fit_ax:
    ax.legend()
for ax in series_ax + fit_ax:
    ax.set_xlabel('')
fit_ax[0].set_ylabel('Fall velocity (m/s)')
fit_ax[1].set_xlabel('Equivalent diameter (mm)')
labels = [a.get_xticklabels() for a in series_ax[:-1]]
labels.extend([a.get_yticklabels() for a in fit_ax[1:]])
plt.setp(labels, visible=False)
axdict['intensity'].set_ylabel('$R$ (mm/h)')
axdict['density'].set_ylabel('$\\rho$ (kg/m$^3$)')
axdict['D_0'].set_ylabel('mm')
axdict['N_w'].set_ylabel('$N_w$')
tfmt = DateFormatter('%H:%M')
series_ax[-1].xaxis.set_major_formatter(tfmt)
