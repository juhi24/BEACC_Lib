# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
import snowfall as sf
import numpy as np
import read
import matplotlib.pyplot as plt
from os import path
import fit
import pandas as pd

from scr_snowfall import pip2015events, test_events, rholimits, param_table, paths

debug = False
unfiltered = False
legacy = False
tld = '.eps'
bnds = (0.25,3.0,0.5,1.8)
savepath = path.join(paths['results'], 'vfits_density_ranges')

def dict2tuple(di):
    return tuple(di[key] for key in sorted(di))

if debug:
    e = test_events()
else:
    e = pip2015events()

#plt.close('all')
plt.ioff()

hextent = np.array(bnds)+(-0.1,0.1,-0.1,0.1)
kws = {'separate': True,
       'source_style': 'hex',
       'source_kws': {'gridsize': 26, 'extent': hextent},
       'unfiltered': True}
fitargs = {'force_flip': False,
           'try_flip': False,
           'fitclass': fit.PolFit}

limslist=sf.limitslist(rholimits)
data = param_table(e=e, debug=debug)
data.index = data.index.droplevel()
vdata = read.merge_series(e.pluv_grouper(), e.vel_data())
vtable = pd.merge(vdata, pd.DataFrame(data.rhomin), left_on='group', right_index=True)
# TODO: check duplicates
vfits = dict()
for rhomin in rholimits[:-1]:
    pipv = read.PipV(data=vtable.query('rhomin=={0}'.format(rhomin)))
    vfits[rhomin], std, etc = pipv.pickler('v-d_rho_range', pipv.find_fit, **fitargs)
fig, axarr = sf.plot_vfits_rho_intervals(dict2tuple(vfits), limslist, **kws)
axarr[0].axis(bnds)
for ax in axarr:
    ax.tick_params(direction='out', top=False)

pc = pd.DataFrame(index=rholimits[:-1], columns=['all', 'hwfm'])
for rhomin in vfits:
    pc.loc[rhomin] = vfits[rhomin].x_unfiltered.size, vfits[rhomin].x.size
pc.loc['total'] = pc.sum()

if debug:
    savepath += '/test'
read.ensure_dir(savepath)
fname = 'combined'
if unfiltered:
    fname += '_unfiltered'
fig.savefig(path.join(savepath, fname + tld))
#fig_fltr.savefig(path.join(savepath, fname + '_d0fltr' + tld))
fig.savefig(path.join(paths['paper'], 'vfits_rho_ranges' + tld))
pc.to_csv(path.join(paths['tables'], 'pc.csv'), sep='\t')
