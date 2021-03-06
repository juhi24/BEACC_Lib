# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
from os import path

dtformat_snex = '%Y %d %B %H UTC'
dtformat_davide = '%d.%m.%y %H:%M'

e = EventsCollection('cases/2015.csv', dtformat_davide)
e.autoimport_data(autoshift=False, autobias=False, rule='5min', 
                  varinterval=True, datafile=['../DATA/test_winter1415.h5'])

basedir = '/home/jussitii/results/pip2015'

plt.close('all')
plt.ioff()

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.dsd.store_good_data() # improve performance by storing filtered dsd tables in memory
    c.pluvio.shift_periods = -5
    c.reset() # reset memory cache after changing pluvio timeshift

for c in e.events.pluvio200.values:
    fig = plt.figure(dpi=100)
    ax = c.density().plot(ylim=[0,600], style='.', title=c.dtstr())
    ax.set_ylabel('bulk density [kg/m^3]')
    plt.savefig(path.join(basedir,'rho_'+c.dtstr().replace(' ','')+'.eps'))
    dtstr = datetime.strftime(c.pluvio.dt_start().to_datetime(),'%Y%m%d%H')
    savedir = os.path.join(basedir, 'velplots', dtstr[0:8], c.pluvio.name)
    ensure_dir(savedir)
    c.pipv.plots(save=True, savedir=savedir, suffix='.eps', rule=c.rule, ymax=3)

#fig1 = plt.figure(dpi=120)
#ax1 = e.plot_pairs(c='density', sizecol='count', vmin=0, vmax=300,
#             query='density<600 & count>1000 & b>0', colorbar=True, 
#             xlim=[0.5,2.5])
#plt.tight_layout()

#fig2 = plt.figure(dpi=120)
#ax2 = e.plot_pairs(x='D_0', y='density', sizecol='count', vmin=0, vmax=800,
#             query='density<600 & count>1000 & b>0', colorbar=False, xlim=[0,8], ylim=[0,600])
#plt.tight_layout()

fig4 = plt.figure(dpi=120)
ax4 = e.plot_pairs(c='k', x='D_0_gamma', y='density', sizecol='count', scale=0.5,
             query='density<600 & count>1000 & b>0', colorbar=False, xlim=[0,6], ylim=[0,500])
plt.tight_layout()


#fig3 = plt.figure(dpi=120)
#ax3 = e.plot_pairs(ax=ax3, x='D_max', y='b',c='density', sizecol='count', vmin=0, vmax=300,
#             query='density<600 & count>1000 & b>0', colorbar=True)
#plt.tight_layout()

brandes = read.PolFit(params=[178, -0.922])
brandes.plot(ax=ax4, label='Brandes et al.')

s = e.summary()
rho_d0 = read.PolFit(x=s.D_0_gamma, y=s.density, sigma=1/s['count'], xname='D_0')
rho_d0.find_fit()
rho_d0.plot(ax=ax4)

plt.legend()
