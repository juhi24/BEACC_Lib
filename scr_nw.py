# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
import fit
from scr_snowfall import param_table

def nw_dmitri(rho, d0):
    return 10**(5.5-(3.57*(rho/1000)**(1/3)-1)*d0)

def nw(rho, d0):
    return 10**(5.6-(0.403*rho**(1/3)-1.2)*d0)

def alpha_fit():
    alpha = np.array([.465, .851, 1.301])
    rho = np.array([72.7, 134, 243.8])
    afit = fit.LinFit(x=rho**(1/3), y=alpha)
    afit.find_fit()
    return afit

plt.close('all')

#data = param_table()
#data['D_0_rho'] = data.D_0_gamma*data.density
#ax2 = data.plot.scatter(x='D_0_rho', y='N_w', c='density', logy=True)
data['D_0_rho1'] = data.D_0_gamma*data.density**(1/3)
ax3 = data.plot.scatter(x='D_0_rho1', y='N_w', c='density', logy=True)
#data['D_0_rho2'] = data.D_0_gamma*(data.density**(1/3)*0.403-1.2)
#ax4 = data.plot.scatter(x='D_0_rho2', y='N_w', c='density', logy=True)
data['nw_rho'] = nw(data.density, data.D_0_gamma)
ax = data.plot.scatter(x='N_w', y='nw_rho', c='D_0_gamma', logx=True, logy=True)
ax5 = data.plot.scatter(x='D_0_gamma', y='N_w', c='density', logy=True)
#ax.set_xlim([data.N_w.min(), data.N_w.max()])
#ax.set_ylim([nw1.min(), nw1.max()])
ax.set_xlabel('$N_w$')
ax.set_ylabel('$N_w(\\rho)$')
ax.plot([0,1e8],[0, 1e8])
ax.axis([1e1,1e6,1e1,1e6])

d0_nw = fit.LinFit(x=data.D_0_gamma, y=np.log10(data.N_w), xname='D_0')
d0_nw.find_fit()
plt.figure()
d0_nw.plot(source_style='raw')

d0rho_nw = fit.LinFit(x=data.D_0_rho1, y=np.log10(data.N_w), xname='D_0*rho^(1/3)')
d0rho_nw.find_fit()
plt.figure()
d0rho_nw.plot(source_style='raw')

#plt.figure()
#afit = alpha_fit()
#afit.plot(source_style='raw')