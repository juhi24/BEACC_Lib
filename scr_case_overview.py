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

data.loc[:,['density', 'D_0', 'N_0']].plot(subplots=True, drawstyle='steps')