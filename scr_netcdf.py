# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
from os import path
from glob import glob
from netCDF4 import Dataset
from datetime import datetime, date

resultspath = '../results/netcdf'

title = 'BAECC 2014 PIP Particle Imaging Package'
subject = 'particle size distribution, fall velocity'
desc = 'Data set includes daily measurements of particle size distribution and fall velocity observed by video camera with high frame rate.'
lang = 'en'
place = 'projection=wgs84 north=61.8436892 east=24.28776892 name=Hyytiälä elevation=150 zunits=m'
dtype = 'netCDF4'
creator = 'Matti Leskinen'
owner = 'University of Helsinki, Department of Physics, Division of Atmospheric Sciences'
contributor = 'Annakaisa von Lerber'
cemail = 'dmitri.moisseev@helsinki.fi'
cphone = '+358294150866'
comment = 'Instrument located on the BAECC measurement field.'
version = 'PIP_rev_1308a'

date_start = date(2014, 2, 21)
#date_end = date(2014, 9, 12)
date_end = date(2014, 2, 22)

for day in daterange(date_start, date_end):
    dtstr = datetime.strftime(day, '%Y%m%d')
    print(dtstr)
    pipv_files = datafilelist(pipv_subpath % dtstr)
    dsd_files = datafilelist(dsd_subpath % dtstr)
    dsd = read.PipDSD(dsd_files)
    pipv = read.PipV(pipv_files)