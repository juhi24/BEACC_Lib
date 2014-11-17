# -*- coding: utf-8 -*-
"""
@author: jussitii
"""
from snowfall import *
import numpy as np

instr = batch_import(dtstr='2014022[1-2]', datadir='../DATA')
m200 = Case(instr['dsd'], instr['vel'], instr['pluvio200'], rule='5min',
               liquid=False)

m200.dsd.data.drop([26.0], 1, inplace=True)