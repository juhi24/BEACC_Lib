# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""

import read
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from os import path

datapath = '../DATA/tartu2015'

dT = 60.0 # disdrometer sampling time, seconds

# Taken from http://www.helsinki.fi/~mleskine/sateet/RD69.html
dsd_D = np.array([0.3671, 0.4643, 0.5611, 0.6667, 0.7822, 0.9245, 1.1275, 
                  1.3424, 1.5171, 1.6759, 1.9227, 2.2695, 2.5945, 2.8797, 
                  3.2091, 3.5557, 3.9284, 4.3634, 4.8738, 5.3892])
dsd_dD = np.array([0.0932, 0.1010, 0.0917, 0.1196, 0.1123, 0.1722, 0.2329,
                   0.1968, 0.1528, 0.1658, 0.3287, 0.3639, 0.2861, 0.2843,
                   0.3745, 0.3196, 0.4239, 0.4471, 0.5735, 0.4563])
dsd_v = np.array([1.4702, 1.9005, 2.3074, 2.7331, 3.1948, 3.7564, 4.4181,
                  5.0182, 5.4520, 5.8193, 6.3375, 7.0271, 7.5608, 7.9154,
                  8.2681, 8.5640, 8.7901, 8.9693, 9.0791, 9.1400])

dsd_data = pd.read_csv(path.join(datapath, 'dsd.csv'), index_col='datetime',
                  parse_dates=True)
dsd_data.columns = dsd_D

dsd = read.DSD(data=dsd_data, binwidth=dsd_dD, bin_v=dsd_v)

sep1 = dsd.between_datetime('2012-09-01 08:00', '2012-09-2 00:00')

#sep1.plot(img=False)

