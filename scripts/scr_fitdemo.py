# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.optimize import curve_fit

e = EventsCollection('cases/cases_of_interest.csv', '%d.%m. %H:%M')
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.pluvio.shift_periods = -6

c=e.events.pluvio400[2]

plt.close('all')
short = c.between_datetime('2014-02-21 23:22:01', '2014-02-21 23:25')
axarr = short.pipv.plots(rule=short.rule, hexbin=True, plotfit=False)
fit = read.PolFit()
vdata = short.pipv.good_data().between_time('23:24:01', '23:25')
params, pcov = curve_fit(fit.func, vdata.Wad_Dia, vdata.vel_v)
fit.params = params
fit.x=vdata.Wad_Dia
fit.y=vdata.vel_v
fit.plot(ax=axarr[2])
plt.legend()
#plt.pcolor(D,V,Z,cmap=plt.cm.gist_earth_r)
plt.figure(dpi=120)
plt.plot(V[:,9],Z[:,9])
plt.title('KDE, D = ' + str(D[:,9][0]) + 'mm')
plt.xlabel('fall velocity (m/s)')
plt.ylabel('KDE value')