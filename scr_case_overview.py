# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
import snowfall as sf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.dates import DateFormatter
from os import path
import read
import warnings

from scr_snowfall import pip2015events, test_events

plt.close('all')
plt.ioff()
debug = False
resultsdir = '../results/pip2015'
paperdir = read.ensure_dir(path.join(resultsdir, 'paper'))
savedir = path.join(resultsdir, 'case_overview')
d0_col = 'D_0_gamma'

if debug:
    e = test_events()
    savedir += '/test'
else:
    e = pip2015events()

read.ensure_dir(savedir)


def select_rows(data, dtlist):
    return data.loc[tuple(map(pd.to_datetime, dtlist)),:]


def gray_out(data, axdict, fltr_label='d0_fltr', labels=['D_0', 'density']):
    for t_end, row in data[data[fltr_label]].iterrows():
        for label in labels:
            ax = axdict[label]
            ax.axvspan(row.start, t_end, edgecolor='none', facecolor='0.8',
                       alpha=0.8)


def markers(data, ax, ycol='density', labelcol='label'):
    for t_end, row in data.iterrows():
        ax.text(row.middle, row[ycol], row[labelcol], ha='center', va='bottom',
                weight='heavy')


def plot_overview(data, params=['intensity', 'density', 'D_0', 'N_w'],
                  axlist=None):
    data.density[data.density>800] = np.nan
    if axlist is None:
        axlist = data.loc[:, params].plot(figsize=(5, 7), subplots=True,
                                         drawstyle='steps')
    else:
        for i, param in enumerate(params):
            ax = axlist[i]
            label = '$'+param+'$'
            if param==d0_col:
                label = '$D_0$'
            data[param].plot(ax=ax, drawstyle='steps', label=label)
    axdict = dict(zip(params, axlist))
    fig = axlist[0].get_figure()
    #axdict['density'].axis((None, None, 0, 600))
    for param in ['N_w']:
        axdict[param].set_yscale('log')
    data = sf.d0fltr(data)
    gray_out(data, axdict, labels=[d0_col, 'density'])
    return fig, axdict


def all_cases_simple_overview(e):
    for case in e.events.paper.values:
        data = case.summary()
        fig, axlist = plot_overview(data, params=params)
        plt.legend()
        plt.tight_layout()
        sdir = read.ensure_dir(path.join(savedir, 'simple'))
        fig.savefig(path.join(sdir, case.dtstr('%Y%m%d') + '.eps'), dpi=150)

params=['intensity', 'density', d0_col, 'N_w']
extent = (0.375, 4, 0.5, 1.5)
xtick_pos = (1, 2, 3, 4)

#all_cases_simple_overview(e)

for ievent, event in e.events.iterrows():
    case = event.paper
    data = case.summary()
    fig = plt.figure(figsize=(6, 9))
    if event.a is np.nan:
        last_ts = case.instr['pipv'].fits.polfit.index[-1]
        dtlist = (last_ts, last_ts, last_ts)
    else:
        dtlist = (event.a, event.b, event.c)
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
                  source_kws={'gridsize': 20, 'extent': extent})
        ax.axis(extent)
        fit_ax.append(ax)
    fig, axdict = plot_overview(data, axlist=series_ax, params=params)
    data['D_max'].plot(ax=axdict[d0_col], drawstyle='steps', label='$D_{max}$')
    # Reverse label order to match line order
    handles, labels = axdict[d0_col].get_legend_handles_labels()
    axdict[d0_col].legend(handles[::-1], labels[::-1], loc='upper left', frameon=False)
    for ax in fit_ax:
        ax.set_xticks(xtick_pos)
        ax.tick_params(axis='both', direction='out', length=4)
        ax.legend()
    for ax in series_ax + fit_ax:
        ax.set_xlabel('')
    sample=data.loc[tuple(map(pd.to_datetime, dtlist)),:]
    sample['label'] = ['a', 'b', 'c']
    sample['ax'] = fit_ax
    for i, row in sample.iterrows():
        row.ax.text(extent[1]-0.2, extent[2], row.label, ha='right', va='bottom',
                    weight='heavy')
    markers(sample, ax=axdict['density'])
    series_ax[0].set_title(case.dtstr())
    fit_ax[0].set_ylabel('$v$, m$\,$s$^{-1}$')
    fit_ax[1].set_xlabel('$D$, mm')
    labels = [a.get_xticklabels() for a in series_ax[:-1]]
    labels.extend([a.get_yticklabels() for a in fit_ax[1:]])
    plt.setp(labels, visible=False)
    axdict['intensity'].set_ylabel('$LWE$, mm$\,$h$^{-1}$')
    axdict['density'].set_ylabel('$\\rho$, ' + read.RHO_UNITS)
    read.rho_scale(axdict['density'].yaxis)
    axdict[d0_col].set_ylabel('mm')
    axdict['N_w'].set_ylabel('$N_w$, mm$^{-1}\,$m$^{-3}$')
    tfmt = DateFormatter('%H:%M')
    series_ax[-1].xaxis.set_major_formatter(tfmt)
    savekws = {'dpi': 150, 'bbox_inches': 'tight'}
    savename = case.dtstr('%Y%m%d') + '.eps'
    try:
        fig.savefig(path.join(savedir, savename), **savekws)
    except ValueError as err:
        warnings.warn("ValueError: {0} Skipping plot for the {1} case.".format(err, case.dtstr()))
    if not debug and event.case_study==True: # yes, necessary to check for True
        fig.savefig(path.join(paperdir, savename), **savekws)
