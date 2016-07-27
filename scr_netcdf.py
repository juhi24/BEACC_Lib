# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
import snowfall as sf
import read
import numpy as np
from os import path
import netCDF4 as cdf
from datetime import datetime, date, timedelta

resultspath = '../results/netcdf'
identifier = 'PIP_004'
instr_prefix = identifier + '_'
ext = '.cdf'

dtformat = '%Y-%m-%dT%H:%M:%S'
timeunits = 'seconds since 1970-01-01 00:00:00 UTC'

date_start = date(2014, 1, 31)
#date_end = date(2014, 9, 12)
date_end = date(2014, 9, 13)

def handle_empty(data, column):
    if data.empty:
        return np.array([])
    return data[column].values

for day in sf.daterange(date_start, date_end):
    dtstr = day.strftime('%Y%m%d')
    print(dtstr)
    out_filepath = path.join(resultspath, instr_prefix + dtstr + ext)
    pipv_files = read.datafilelist(read.PIPV_SUBPATH % dtstr)
    dsd_files = read.datafilelist(read.DSD_SUBPATH % dtstr)
    pipv_is_empty = len(pipv_files) < 1
    dsd_is_empty = len(dsd_files) < 1
    if not dsd_is_empty:
        dsd = read.PipDSD(dsd_files)
        dsd.data.drop_duplicates(inplace=True)
    if not pipv_is_empty:
        pipv = read.PipV(pipv_files)
    
    nc = cdf.Dataset(out_filepath, 'w', format='NETCDF4')
    nc.Title = 'BAECC 2014 PIP Particle Imaging Package'
    nc.setncattr('Title.lan', 'eng')
    nc.Subject = 'particle size distribution, fall velocity'
    nc.Description = 'Data set includes daily measurements of particle size distribution and fall velocity observed by video camera with high frame rate.'
    nc.Language = 'eng'
    nc.Identifier = identifier
    nc.Modified = datetime.utcnow().strftime(dtformat + 'Z')
    nc.TemporalCoverage = 'start=%s; stop=%s UTC' % (day.strftime(dtformat), (day + timedelta(days=1)).strftime(dtformat))
    nc.SpatialCoverage = 'projection=wgs84 north=61.8436892 east=24.28776892 name=Hyytiälä elevation=150 zunits=m'
    nc.Type = 'netCDF4'
    nc.Creator = 'Matti Leskinen'
    nc.Owner = 'University of Helsinki, Department of Physics, Division of Atmospheric Sciences'
    nc.Contributor = 'Annakaisa von Lerber'
    nc.setncattr('Contact.email', 'dmitri.moisseev@helsinki.fi')
    nc.setncattr('Contact.phone', '+358294150866')
    nc.setncattr('Contact.type', 'person')
    nc.setncattr('Project.funder', 'DoE ARM/NASA')
    nc.setncattr('Project.name', 'BAECC SNEX')
    nc.setncattr('Project.homepage', 'http://www.arm.gov/campaigns/amf2014baecc')
    nc.setncattr('Rights.category', 'licensed')
    nc.setncattr('Rights.declaration', 'https://creativecommons.org/licenses/by/4.0/')
    nc.setncattr('Measurement.comment', 'Instrument located on the BAECC measurement field.')
    nc.SoftwareVersion = 'PIP_rev_1308a'
    
    if not dsd_is_empty:
        dtime = nc.createDimension('time', dsd.data.index.values.size)
        vtime = nc.createVariable('time', 'i4', 'time')
        vtime.description = 'POSIX timestamp of the PIP in a minute interval'
        vtime.units = timeunits
        vtime[:] = dsd.data.index.astype(datetime).map(datetime.timestamp).astype(int)
        
        dbinsize = nc.createDimension('bin size', dsd.data.columns.values.size)
        vbins = nc.createVariable('particle_size', 'f4', 'bin size')
        vbins.description = 'Particle bin center of an area-equivalent diameter'
        vbins.units = 'mm'
        vbins[:] = dsd.data.columns.values 
        
        vdsd = nc.createVariable('psd', 'f4', ('time', 'bin size'))
        vdsd.description = 'Particle size distribution of the minutes with more than 10 particles observed'
        vdsd.units = 'm^-3mm^-1'
        vdsd[:] = 2*dsd.data.values
    
    if not pipv_is_empty:
        dvtime = nc.createDimension('time velocity', pipv.data.index.values.size)
        vtime_v = nc.createVariable('time_velocity', 'i4', 'time velocity')
        vtime_v.description = 'POSIX timestamp of the particles observed falling in a minute interval of the fall velocity data'
        vtime_v.units = timeunits
        vtime_v[:] = pipv.data.index.astype(datetime).map(datetime.timestamp).astype(int)
        
        vvel = nc.createVariable('velocity', 'f4', 'time velocity')
        vvel.description = 'Fall velocity of the particle. The observation has more than two video frames.'
        vvel.units = 'm/s'
        vvel[:] = handle_empty(pipv.data, 'vel_v')
        
        vsize_v = nc.createVariable('particle_size_velocity', 'f4', 'time velocity')
        vsize_v.description = 'The area-equivalent diameter of particle of fall velocity data.'
        vsize_v.units = 'mm'
        vsize_v[:] = handle_empty(pipv.data, 'Wad_Dia')
        nc.close()
    