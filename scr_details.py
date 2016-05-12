# -*- coding: utf-8 -*-
"""
Details and statistics script.
@author: Jussi Tiira
"""

from scr_snowfall import param_table
import matplotlib.pyplot as plt
import seaborn as sns

plt.ion()
QSTR_FILTERED = 'density>599 | count<801 | b<0 | density!=density'


def param_table_unfiltered():
    return param_table(query_str='')

def param_table_filtered_rows():
    qstr = QSTR_FILTERED
    return param_table(query_str=qstr)

def plot_intensity_hist(data, bins=18):
    ax = data.intensity.hist(bins=bins)
    ax.set_ylabel('frequency')
    ax.set_xlabel('LWE')