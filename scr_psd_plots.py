# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
import snowfall as sf
import read
import numpy as np
import matplotlib.pyplot as plt
from os import path
import fit

from scr_snowfall import param_table

plt.close('all')
plt.ioff()
debug = False

rholims = (0, 100, 200, 800)
#rholims = (0, 150, 300, 800)
limslist = sf.limitslist(rholims)
n_ranges = len(limslist)
separate = False
d0_col = 'D_0_gamma'

resultspath = '../results/pip2015'
paperpath = read.ensure_dir(path.join(resultspath, 'paper'))
savepath = path.join(resultspath, 'psd')

rhorangestr = '$%s < \\rho \leq %s$'


def d0_nw_paper(data, rholimits):
    rhopairs = sf.limitslist(rholimits)
    fig, axarr = plt.subplots(nrows=1, ncols=3, dpi=150, figsize=(11,4),
                              sharex=True, sharey=True, tight_layout=True)
    for i, (rhomin, rhomax) in enumerate(rhopairs):
        limitstr = rhorangestr % (rhomin*read.RHO_SCALE, rhomax*read.RHO_SCALE)
        ax = axarr[i]
        datarho = data[data.rhomin==rhomin]#.loc['second']
        #datarho.plot(kind='scatter', x='D_0', y='N_w', ax=axarr[i], logy=True)
        datarho.plot(kind='scatter', x=d0_col, y='log_nw', c='', ax=axarr[i])
        lfit = fit.LinFit(x=datarho[d0_col], y=datarho.log_nw, xname='D_0')
        lfit.find_fit()
        lfit.plot(source_style=None, ax=ax)
        #efit = fit.ExponentialFit(x=datarho.D_0, y=datarho.N_w)
        #efit.find_fit()
        #efit.plot(source_style='raw', ax=ax)
        ax.legend()
        ax.set_title(limitstr)
        ax.set_xlabel('')
    axarr[0].set_ylabel('$log(N_w)$')
    axarr[1].set_xlabel('$D_0$, mm')
    axarr[-1].set_title('$\\rho > %s$' % (rhomin*read.RHO_SCALE))
    plt.axis([0, 5, 1, 6])
    #remove_subplot_gaps(fig, axis='row')
    return fig, axarr


data = param_table(debug=debug)
data['log_nw'] = np.log10(data['N_w'])
data = sf.apply_rho_intervals(data, limits=rholims)
fig, axarr = d0_nw_paper(data, rholimits=rholims)
if debug:
    savepath += '/test'
read.ensure_dir(savepath)
fig.savefig(path.join(savepath, 'nw_d0.eps'))
fig.savefig(path.join(paperpath, 'nw_d0.eps'))
