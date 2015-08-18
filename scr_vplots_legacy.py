# -*- coding: utf-8 -*-
"""
@author: jussitii
"""
from snowfall import *
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.rc('font', size=15)
plt.rc('grid', linewidth=1)
plt.rc('xtick.major', width=1.5)
plt.rc('ytick.major', width=1.5, size=5)
plt.rc('legend', fontsize='small')
plt.rc('axes', titlesize='medium')

plt.close('all')

e = EventsCollection('cases/cases_of_interest.csv', '%d.%m. %H:%M')
e.autoimport_data(autoshift=False, rule='6min')

for c in e.events.pluvio200.values:
    c.pipv.find_fits(rule='6min')
    c.pipv.plots(save=True, suffix='.eps', grid=False, xmax=4, ymax=3,
                 xticks=[0, 1, 2, 3, 4], yticks=[0, 1, 2, 3],
                 colorbar=False, hexsize=8)
          
f = plt.figure(dpi=175, figsize=(1, 3))
ax = f.add_axes([0.05, 0.05, 0.3, 0.9])
norm = mpl.colors.Normalize(vmin=0, vmax=100)
cb = mpl.colorbar.ColorbarBase(ax, cmap='gray_r', norm=norm, 
                               orientation='vertical')
cb.set_label('%')