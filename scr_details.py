# -*- coding: utf-8 -*-
"""
Details and statistics script.
@author: Jussi Tiira
"""

import snowfall as sf
from scr_snowfall import param_table, QSTR
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

plt.ion()
#QSTR_FILTERED = 'density>599 | count<801 | b<0 | density!=density'
QSTR_FILTERED = 'D_0_gamma<0.6 | intensity<0.2 | count<801'

def bins(xmin=0, xmax=6, binsize=0.5):
    bins = np.arange(xmin, xmax+binsize, binsize)
    return bins

def param_table_unfiltered():
    return param_table(query_str='')

def param_table_filtered_rows():
    qstr = QSTR_FILTERED
    return param_table(query_str=qstr)

def plot_intensity_hist(data, **kws):
    ax = data.intensity.hist(**kws)
    ax.set_ylabel('Frequency')
    ax.set_xlabel('LWE')
    return ax

def plot_intensity_binned_cdf(data, **kws):
    ax = data.intensity.hist(cumulative=True, normed=True, **kws)
    ax.set_ylabel('CDF')
    ax.set_xlabel('LWE')
    ax.axis([None, None, 0, 1])
    return ax

def plot_intensity_cdf(data, xmin=0, xmax=6):
    cdf = sf.series_cdf(data.intensity)
    ax = cdf.plot(drawstyle='steps')
    ax.set_ylabel('CDF')
    ax.set_xlabel('LWE')    
    ax.axis([xmin, xmax, 0, 1])
    return ax
