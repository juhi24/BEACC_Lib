# coding: utf-8
import locale
from os import path
from j24 import home
from baecc.tools import merge_series, merge_multiseries

locale.setlocale(locale.LC_ALL, 'C')

# general configuration
DEBUG = False
CGS_UNITS = True # display units in cgs instead of SI

# CONFIG default paths
HOME = home()
DATA_DIR = path.join(HOME, 'DATA')
USER_DIR = path.join(HOME, '.baecc')
RESULTS_DIR = path.join(HOME, 'results')
H5_PATH = path.join(DATA_DIR, 'baecc.h5')

# constants
if CGS_UNITS:
    RHO_SCALE = 1e-3
    RHO_UNITS = 'g$\,$cm$^{-3}$'
else:
    RHO_SCALE = 1
    RHO_UNITS = 'kg$\,$m$^{-3}$'
