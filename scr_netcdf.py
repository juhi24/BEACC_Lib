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
instr_prefix = 'PIP_004_'
ext = '.cdf'

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

date_start = date(2014, 2, 4)
#date_end = date(2014, 9, 12)
date_end = date(2014, 2, 5)

for day in daterange(date_start, date_end):
    dtstr = datetime.strftime(day, '%Y%m%d')
    print(dtstr)
    out_filepath = path.join(resultspath, instr_prefix + dtstr + ext)
    pipv_files = datafilelist(pipv_subpath % dtstr)
    dsd_files = datafilelist(dsd_subpath % dtstr)
    dsd = read.PipDSD(dsd_files)
    pipv = read.PipV(pipv_files)
    
    nc = Dataset(out_filepath, 'w', format='NETCDF4')
    
    dbinsize = nc.createDimension('bin size', 105)
    vsize = nc.createVariable('Particle_size', 'f8', 'bin size')
    vsize.description = 'Particle bin center of an area-equivalent diameter'
    vsize.units = 'mm'
    
    dtime = nc.createDimension('time', 18)
    vtime = nc.createVariable('Time', str, 'time')
    vtime.description = 'Timestamp of the PIP in a minute interval'
    vtime.units = 'UTC'
    
    vtime_v = nc.createVariable('Time_velocity', str, 'time')
    vtime_v.description = 'Time stamp of the particles observed falling in a minute interval of the fall velocity data'
    vtime_v.units = 'yyyy-mm-dd HH:MM:00'
    
    
    dvel = nc.createDimension('time_velocity')
    vvel = nc.createVariable('Velocity', 'f4', 'time_velocity')
    vvel.description = 'Fall velocity of the particle. The observation has more than two video frames.'
    vvel.units = 'm/s'
    vvel[:] = pipv.data.vel_v.values
    
    vsize_v = nc.createVariable('Particle_size_velocity', 'f4', 'time_velocity')
    vsize_v.description = 'The area-equivalent diameter of particle of fall velocity data.'
    vsize_v.units = 'mm'
    vsize_v[:] = pipv.data.Wad_Dia.values
    nc.close()