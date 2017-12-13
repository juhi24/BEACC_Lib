# coding: utf-8
from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import matplotlib.pyplot as plt
import pandas as pd
from glob import glob
from baecc.instruments import pluvio

t0=pd.datetime(2014,3,3,5,0)
t1=pd.datetime(2014,3,3,9,0)
fnames=glob('/home/jussitii/DATA/Pluvio200/pluvio200_*_20140303*.txt')
pluv = pluvio.Pluvio(filenames=fnames, dt_start=t0, dt_end=t1, name='pluvio200')



