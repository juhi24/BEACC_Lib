# -*- coding: utf-8 -*-
"""
Details and statistics script.
@author: Jussi Tiira
"""

from scr_snowfall import param_table

def param_table_unfiltered():
    return param_table(query_str='')

def param_table_filtered_rows():
    qstr = 'density>599 | count<801 | b<0 | density!=density'
    return param_table(query_str=qstr)
