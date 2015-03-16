# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
from os import path

plt.close('all')
plt.ioff()

def init_dataset(e, basedir='..', shift=-6):
    for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
        c.dsd.store_good_data() # improve performance by storing filtered dsd tables in memory
        c.pluvio.shift_periods = shift
        c.reset() # reset memory cache after changing pluvio timeshift
        dtstr = datetime.strftime(c.pluvio.dt_start().to_datetime(),'%Y%m%d%H')
        savedir = path.join(basedir, dtstr[0:8], c.pluvio.name)
        read.ensure_dir(savedir)
        fig1 = plt.figure(dpi=100)
        ax1 = plt.gca()
        c.density().plot(ax=ax1, ylim=[0,600], style='.', title=c.dtstr())
        ax1.set_ylabel('bulk density [kg/m^3]')
        fig1.savefig(path.join(savedir, 'rho.eps'))
        fig2 = plt.figure(dpi=100)
        ax2 = plt.gca()
        c.pipv.fit_params().plot(ax=ax2)
        plt.legend()
        fig2.savefig(path.join(savedir, 'ab.eps'))
        c.summary().drop('case',axis=1).to_csv(path.join(savedir,'summary.csv'))
        c.pluvio.tdelta().to_csv(path.join(savedir, 'timedelta.csv'))
        plt.close('all')

dtformat_snex = '%Y %d %B %H UTC'
dtformat_davide = '%d.%m.%y %H:%M'

w1415 = EventsCollection('cases/2015.csv', dtformat_davide)
w1415.autoimport_data(autoshift=False, autobias=False, rule='5min',
                      varinterval=True, datafile=['../DATA/winter1415.h5'])

w14 = EventsCollection('cases/pip2015.csv', dtformat_snex)
w14.autoimport_data(autoshift=False, autobias=False, rule='5min',
                    varinterval=True, datafile=['../DATA/baecc.h5'])

basedir = '/home/jussitii/results/pip2015'

init_dataset(w14, basedir=basedir, shift=-6)
init_dataset(w1415, basedir=basedir, shift=-5)

#e = combine_datasets(w14, w1415)
