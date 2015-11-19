# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
from datetime import datetime
from os import path

from scr_snowfall import pip2015events

#plt.close('all')
plt.ion()

basepath = read.ensure_dir('../results/pip2015/paper/density')
dtfmt = '%Y%m%d'

e = pip2015events()

rho_label = 'bulk density (kg m$^{-3}$)'
t_label = 'time'

for c in e.events.paper.values:
    savepath = basepath
    rho = c.density()
    rho.to_csv(path.join(savepath, c.dtstr(dtfmt) + '.csv'))
    c.instr['pluvio'].tdelta().to_csv(path.join(savepath, 'timedelta_' + c.dtstr(dtfmt) + '.csv'))
    plt.figure(dpi=120)
    rho.plot(drawstyle='steps')
    plt.title(c.dtstr())
    plt.xlabel(t_label)
    plt.ylabel(rho_label)
    plt.ylim((0,500))
    plt.savefig(path.join(savepath, c.dtstr(dtfmt) + '_single.eps'))

for i, ev in e.events.iterrows():
    plt.figure(dpi=120)
    for c in (ev.pluvio200, ev.pluvio400):
        if c is None:
            continue
        rho = c.density()
        rho.plot(drawstyle='steps', label=c.instr['pluvio'].name)
    plt.title(c.dtstr())
    plt.xlabel(t_label)
    plt.ylabel(rho_label)
    plt.ylim((0,500))
    plt.legend(loc='lower right')
    plt.savefig(path.join(basepath, c.dtstr(dtfmt) + '.eps'))
