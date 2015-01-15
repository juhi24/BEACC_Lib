# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt

dtformat_default = '%d.%m. %H:%M'
dtformat_snex = '%d %B %H UTC'

e = EventsCollection('cases/pip2015.csv', dtformat_snex)
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

basedir = '/home/jussitii/results/pip2015'

plt.ion()
#ax = plt.gca()
#plt.figure(dpi=120)

for c in e.events.pluvio200.values:
    c.pluvio.shift_periods = -6
    #c.plot_velfitcoefs(ax=ax, rhomax=600, countmin=2000)

#c.plot_velfitcoefs(rhomax=600, countmin=2000)

def plot_velfitcoefs(c, fig=None, ax=None, rhomin=None, rhomax=None, countmin=1, **kwargs):
    rho = c.density().replace(np.inf, np.nan)
    params = c.pipv.fits.polfit.apply(lambda fit: fit.params)
    selection = pd.DataFrame([rho.notnull(), c.partcount() > countmin]).all()
    rho = rho[selection]
    params = params[selection]
    a = params.apply(lambda p: p[0]).values
    b = params.apply(lambda p: p[1]).values
    if fig is None:
        fig = plt.figure(dpi=120)
    if ax is None:
        ax = plt.gca()
    if rhomin is None:
        vmin = rho.min()
    if rhomax is None:
        vmax = rho.max()
    choppa = ax.scatter(b,rho.values,c=a, vmin=rhomin, vmax=rhomax,
                        **kwargs)
    cb = fig.colorbar(choppa, label='a')
    ax.set_xlabel('$b_u$', fontsize=15)
    ax.set_ylabel('$\rho_u$', fontsize=15)
    return ax