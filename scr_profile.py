# -*- coding: utf-8 -*-
"""
A script for profiling purposes.

@author: Jussi Tiira
"""
from snowfall import *
import matplotlib.pyplot as plt
from os import path
import read

from scr_snowfall import pip2015events, test_events

debug = True

if debug:
    e = test_events()
else:
    e = pip2015events()

data = e.summary()