# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
from os import path

dtformat_default = '%d.%m. %H:%M'
dtformat_snex = '%Y %d %B %H UTC'

e = EventsCollection('cases/pip2015.csv', dtformat_snex)
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

#plt.close('all')
#plt.ioff()
plt.ion()

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.instr['pluvio'].shift_periods = -6
    c.instr['pluvio'].n_combined_intervals = 2

comb200 = e.events.pluvio200.sum()
comb400 = e.events.pluvio400.sum()
del(e)

limslist = limitslist((0, 150, 300, 800))
n_ranges = len(limslist)
separate = False

basepath = '../results/pip2015/paper/psd'

for c in (comb200, comb400):
    d0 = c.d_0()
    nt = c.n_t()
    n = c.intervalled(c.instr['dsd'].n_all)
    y = n*d0/nt
    d = y*0+y.columns.values
    x = d#/d0
    if not separate:
        fig, axarr = plt.subplots(1, n_ranges, dpi=100, sharex=True, sharey=True,
                                  tight_layout=True, figsize=(n_ranges*6, 6))
    for i, rhorange in enumerate(limslist):
        if separate:
            fig, ax = plt.subplots(dpi=100)
        else:
            ax = axarr[i]
        rhomin = rhorange[0]
        rhomax = rhorange[1]
        xpart = c.data_in_density_range(x, rhomin, rhomax)
        ypart = c.data_in_density_range(y, rhomin, rhomax)
        if xpart.empty or ypart.empty:
            ax.semilogy()
            continue
        ax.semilogy(xpart, ypart, marker='o', linestyle='None', color='black')
        ax.set_ylabel('$N_D D_0 N_T^{-1}$')
        ax.set_xlabel('$D D_0^{-1}$')
        ax.set_title('$%s < \\rho < %s$' % (rhomin, rhomax))
        ax.grid(True)
    plt.axis([0, 5, 10e-5, 10])
    ax.set_title('$\\rho > ' + str(rhomin) + '$')
