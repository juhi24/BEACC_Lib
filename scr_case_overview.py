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
plt.ioff()
debug = False
savedir = '../results/pip2015/case_overview'

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
            data[param].plot(ax=ax, drawstyle='steps', label='$'+param+'$')
    axdict = dict(zip(params, axlist))
    fig = axlist[0].get_figure()
    #axdict['density'].axis((None, None, 0, 600))
    for param in ['N_w']:
        axdict[param].set_yscale('log')
    data = d0fltr(data)
    gray_out(data, axdict)
    return fig, axdict


def all_cases_simple_overview(e):
    for case in e.events.paper.values:
        data = case.summary()
        fig, axlist = plot_overview(data, params=params)
        plt.legend()
        plt.tight_layout()
        sdir = read.ensure_dir(path.join(savedir, 'simple'))
        fig.savefig(path.join(sdir, case.dtstr('%Y%m%d') + '.eps'), dpi=150)

params=['intensity', 'density', 'D_0', 'N_w']
extent = (0.375, 4, 0.5, 1.5)
xtick_pos = (1, 2, 3, 4)

#all_cases_simple_overview(e)

for icase, case in e.events.paper.iteritems():
    data = case.summary()
    fig = plt.figure(figsize=(6, 9))
    if debug:
        dtlistlist = (('2014-03-20 17:27:00', '2014-03-20 18:31:00', '2014-03-20 19:11:00'),
                      ('2014-12-18 15:55:00', '2014-12-18 16:14:00', '2014-12-18 16:43:00'),
                      ('2015-01-14 02:53:00', '2015-01-14 03:30:00', '2015-01-14 03:40:00'))
    else:
        dtlistlist = (('2014-02-01 00:21:00', '2014-02-01 01:55:00', '2014-02-01 02:19:00'),
                      ('2014-02-12 06:21:00', '2014-02-12 07:28:00', '2014-02-12 07:42:00'),
                      ('2014-02-15 22:50:00', '2014-02-15 23:33:00', '2014-02-15 23:48:00'),
                      ('2014-02-21 19:04:00', '2014-02-21 23:35:00', '2014-02-22 00:31:00'),
                      ('2014-03-18 09:37:00', '2014-03-18 10:19:00', '2014-03-18 16:39:00'),
                      ('2014-03-20 17:27:00', '2014-03-20 18:31:00', '2014-03-20 19:11:00'),
                      ('2014-12-18 15:55:00', '2014-12-18 16:14:00', '2014-12-18 16:43:00'),
                      ('2014-12-30 04:41:00', '2014-12-30 07:24:00', '2014-12-30 08:26:00'),
                      ('2015-01-03 14:21:00', '2015-01-03 15:46:00', '2015-01-03 19:46:00'),
                      ('2015-01-08 11:56:00', '2015-01-08 12:18:00', '2015-01-08 13:16:00'),
                      ('2015-01-09 20:24:00', '2015-01-09 21:02:00', '2015-01-09 23:48:00'),
                      ('2015-01-11 05:23:00', '2015-01-11 05:23:00', '2015-01-11 05:23:00'),
                      ('2015-01-14 03:07:00', '2015-01-14 03:34:00', '2015-01-14 04:03:00'),
                      ('2015-01-18 17:08:00', '2015-01-18 18:11:00', '2015-01-18 19:17:00'),
                      ('2015-01-22 22:30:00', '2015-01-23 01:16:00', '2015-01-23 02:14:00'),
                      ('2015-01-23 16:07:00', '2015-01-23 18:49:00', '2015-01-23 20:16:00'),
                      ('2015-01-25 10:12:00', '2015-01-25 12:33:00', '2015-01-25 13:36:00'))
    dtlist = dtlistlist[icase[1]]
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
    data['D_max'].plot(ax=axdict['D_0'], drawstyle='steps', label='$D_{max}$')
    # Reverse label order to match line order
    handles, labels = axdict['D_0'].get_legend_handles_labels()
    axdict['D_0'].legend(handles[::-1], labels[::-1], loc='upper left', frameon=False)
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
    axdict['density'].set_ylabel('$\\rho$, kg$\,$m$^{-3}$')
    axdict['D_0'].set_ylabel('mm')
    axdict['N_w'].set_ylabel('$N_w$')
    tfmt = DateFormatter('%H:%M')
    series_ax[-1].xaxis.set_major_formatter(tfmt)
    fig.savefig(path.join(savedir, case.dtstr('%Y%m%d') + '.eps'), dpi=150,
                bbox_inches='tight')
