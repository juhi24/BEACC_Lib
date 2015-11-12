# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import matplotlib.pyplot as plt
from os import path
import read

from scr_snowfall import pip2015events, test_events

debug = False

if debug:
    e = test_events()
else:
    e = pip2015events()

case = e.events.paper[0]
data = case.summary()

params = ['intensity', 'density', 'D_0', 'N_0', 'N_w']
axarr = data.loc[:, params].plot(subplots=True, drawstyle='steps')

axdict = dict(zip(params, axarr))

for param in ['N_0', 'N_w']:
    axdict[param].set_yscale('log')
