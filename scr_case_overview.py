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


def isnan(var):
    return var != var


def select_rows(data, dtlist):
    return data.loc[tuple(map(pd.to_datetime, dtlist)),:]


def gray_out(data, axdict, fltr_label='d0_fltr', labels=['D_0', 'density']):
    for t_end, row in data[data[fltr_label]].iterrows():
        for label in labels:
            ax = axdict[label]
            ax.axvspan(row.start, t_end, edgecolor='none', facecolor='0.8',
                       alpha=0.8)


def markers(data, ax, ycol='density', labelcol='label', mark_below=True, **kws):
    for t_end, row in data.iterrows():
        label = row[labelcol]
        va = row['va']
        if isnan(va):
            va = 'bottom'
        if mark_below:
            y = 0
            va = 'top'
        else:
            y = row[ycol]
        ax.text(row.middle, y, label, ha='center', va=va,
                weight='heavy', **kws)


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
    for param in ['N_w']:
        axdict[param].set_yscale('log')
    data = sf.d0fltr(data)
    gray_out(data, axdict, labels=[d0_col, 'density'])
    return fig, axdict


def plot_vfit(case, dt, ax=None, extent=(0.375, 4, 0.5, 1.5),
              xtick_pos=(1, 2, 3, 4)):
    if ax is None:
        ax = plt.gca()
    vfit = case.instr['pipv'].fits.polfit[dt]
    vfit.plot(source_style='hex', unfiltered=True, ax=ax,
              source_kws={'gridsize': 20, 'extent': extent})
    ax.axis(extent)
    dt_start = case.instr['pluvio'].start_time()[dt]
    dt_end = pd.to_datetime(dt)
    ax.set_title('{0}–{1}'.format(dt_start.strftime(tformat),
                                  dt_end.strftime(tformat)))
    ax.set_xticks(xtick_pos)
    ax.tick_params(axis='both', direction='out', length=4)
    ax.legend()


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
tformat = '%H:%M'
#all_cases_simple_overview(e)

for ievent, event in e.events.iterrows():
    case = event.paper
    data = case.summary()
    selection = data.tdelta>pd.datetools.Minute(30).delta
    data.loc[selection,['density','D_max','D_0_gamma', 'N_w']]=np.nan
    if data.density.dropna().empty:
        continue
    fig = plt.figure(figsize=(6, 9))
    if event.a is np.nan:
        last_ts = case.instr['pipv'].fits.polfit.index[-1]
        dtlist = (last_ts, last_ts, last_ts)
        valist = (np.nan, np.nan, np.nan)
    else:
        dtlist = (event.a, event.b, event.c)
        valist = (event.va_a, event.va_b, event.va_c)
    series_ax = []
    fit_ax = []
    gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
    gs_series = gridspec.GridSpecFromSubplotSpec(len(params), 1, subplot_spec=gs[0],
                                                 hspace=0.15)
    gs_fit = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1],
                                              wspace=0.05)
    gs.update(hspace=0.25)
    for i in range(len(params)):
        series_ax.append(plt.subplot(gs_series[i]))
    for i, dt in enumerate(dtlist):
        ax = plt.subplot(gs_fit[i])
        plot_vfit(case, dt, extent=extent, xtick_pos=xtick_pos, ax=ax)
        fit_ax.append(ax)
    fig, axdict = plot_overview(data, axlist=series_ax, params=params)
    data['D_max'].plot(ax=axdict[d0_col], drawstyle='steps', label='$D_{max}$')
    # Reverse label order to match line order
    handles, labels = axdict[d0_col].get_legend_handles_labels()
    axdict[d0_col].legend(handles[::-1], labels[::-1], loc='upper center',
                          frameon=False)
    for ax in series_ax + fit_ax:
        ax.set_xlabel('')
    sample=data.loc[tuple(map(pd.to_datetime, dtlist)),:]
    sample['label'] = ['a', 'b', 'c']
    sample['va'] = valist
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
    axdict['density'].set_ylim(0, 250)
    axdict['intensity'].set_ylim(0, 3.5)
    axdict[d0_col].set_ylim(0, 16)
    axdict['N_w'].set_ylim(1e2, 1e6)
    read.rho_scale(axdict['density'].yaxis)
    axdict[d0_col].set_ylabel('mm')
    axdict['N_w'].set_ylabel('$N_w$, mm$^{-1}\,$m$^{-3}$')
    tfmt = DateFormatter('%H:%M')
    series_ax[-1].xaxis.set_major_formatter(tfmt)
    # saving
    savekws = {'dpi': 150, 'bbox_inches': 'tight'}
    savename = case.dtstr('%Y%m%d') + '.eps'
    try:
        fig.savefig(path.join(savedir, savename), **savekws)
    except ValueError as err:
        warnings.warn("ValueError: {0} Skipping plot for the {1} case.".format(err, case.dtstr()))
    if not debug and event.case_study==True: # yes, necessary to check for True
        fig.savefig(path.join(paperdir, savename), **savekws)
