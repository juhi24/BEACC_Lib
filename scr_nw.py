# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
import numpy as np
import matplotlib.pyplot as plt
import fit
from scr_snowfall import param_table, paths
from os import path

savepath = path.join(paths['results'], 'psd')

plt.close('all')
plt.ioff()


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

data = param_table()

#data['D_0_rho'] = data.D_0_gamma*data.density
#ax2 = data.plot.scatter(x='D_0_rho', y='N_w', c='density', logy=True)
data['D_0_rho1'] = data.D_0_gamma*data.density**(1/3)
#ax3 = data.plot.scatter(x='D_0_rho1', y='N_w', c='density', logy=True)
#data['D_0_rho2'] = data.D_0_gamma*(data.density**(1/3)*0.403-1.2)
#ax4 = data.plot.scatter(x='D_0_rho2', y='N_w', c='density', logy=True)
data['nw_rho'] = nw(data.density, data.D_0_gamma)
#ax = data.plot.scatter(x='N_w', y='nw_rho', c='D_0_gamma', logx=True, logy=True)
#ax5 = data.plot.scatter(x='D_0_gamma', y='N_w', c='density', logy=True)
#ax.set_xlim([data.N_w.min(), data.N_w.max()])
#ax.set_ylim([nw1.min(), nw1.max()])
#ax.set_xlabel('$N_w$')
#ax.set_ylabel('$N_w(\\rho)$')
#ax.plot([0,1e8],[0, 1e8])
#ax.axis([1e1,1e6,1e1,1e6])
legendkws = {'frameon':True}

txtkws = {'x':0.15, 'y':0.1, 'ha':'right', 'va':'bottom'}
skws = {'c':'', 'label':'data points'}
legendkws = {'frameon':True}
fig, (axd0, axd0rho) = plt.subplots(ncols=2, sharey=True, dpi=150,
                                    figsize=(7,4), tight_layout=True)
d0_nw = fit.ExponentialFit(x=data.D_0_gamma, y=data.N_w, xname='D_0',
                           use_latex_fmt=True)
d0_nw.find_fit()
fitlabel = ''
d0_nw.plot(ax=axd0, source_style='raw', source_kws=skws)
axd0.set_yscale('log')
axd0.legend(**legendkws)
axd0.text(s='(a)', transform=axd0.transAxes, **txtkws)

d0rho = 0.001**(1/3)*data.D_0_rho1
d0rho_nw = fit.ExponentialFit(x=d0rho, y=data.N_w, xname='D_0\\rho^{^1\!/_3}',
                              use_latex_fmt=True)
d0rho_nw.find_fit()
axd0rho = d0rho_nw.plot(ax=axd0rho, source_style='raw', source_kws=skws)
axd0rho.set_yscale('log')
axd0rho.legend(**legendkws)
axd0rho.text(s='(b)', transform=axd0rho.transAxes, **txtkws)

axd0rho.axis([None, 2.5, 1e1, 1e6])
axd0.set_ylabel('$N_w$, mm$^{-1}$m$^{-3}$')
axd0.set_xlabel('$D_0$, mm')
axd0rho.set_xlabel('$D_0 \\rho^{^1\!/_3}$, g$^{^1\!/_3}$')
fig.savefig(path.join(savepath, 'nw_d0.eps'))
fig.savefig(path.join(paths['paper'], 'nw_d0.eps'))
#plt.figure()
#afit = alpha_fit()
#afit.plot(source_style='raw')