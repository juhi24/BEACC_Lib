# coding: utf-8
import numpy as np
import pandas as pd
import datetime
from scipy import io
from os import path
from baecc import read


class Radar(read.InstrumentData):
    """Radar reflectivity at lowest level"""
    def __init__(self, filenames=None, dt_start=None, dt_end=None, **kwargs):
        """Create vertical pointing Radar object using data from various radar
        modes"""
        self._time_lag = pd.to_timedelta(0.0, unit='s')
        read.InstrumentData.__init__(self, filenames, **kwargs)
        if self.data.empty and filenames:
            print('Reading Radar data...')
            self.name = (path.basename(path.dirname(self.filenames[0])))
            for filename in filenames:
                print(filename)
                radardata = io.netcdf.netcdf_file(filename)
                radarvariables = radardata.variables
                if filename.endswith('.nc'):
                    if 'XSACR' in radardata.title.decode():
                        range_idx = 1
                    if 'KaSACR' in radardata.title.decode():
                        range_idx = 0
                    refl = radarvariables['reflectivity']
                    reflectivity = 10.0**(0.1*(refl.data[:, range_idx]*refl.scale_factor + refl.add_offset))
                    basetime = datetime.datetime.strptime(radarvariables['time'].units.decode(),
                                                          'seconds since %Y-%m-%dT%H:%M:%SZ')
                    delta = radarvariables['time'].data
                    deltatime = pd.to_timedelta(np.round(delta), unit='s')
                    time = basetime + deltatime
                    elevation = radarvariables['elevation'].data
                    #rng = radarvariables['range'].data
                    VP = np.abs(elevation-90.0) < 0.5
                    tmpDF = pd.DataFrame(reflectivity[VP], index=time[VP],
                                         columns=['reflectivity'],
                                         dtype=np.float64)
                    self.data = self.data.append(tmpDF)
                elif filename.endswith('.cdf'):
                    if 'reflectivity_copol' in radarvariables.keys():
                        range_idx = 10
                        reflectivity = 10.0**(0.1*radarvariables['reflectivity_copol'].data[:, range_idx].byteswap().newbyteorder())
                    elif 'reflectivity' in radarvariables.keys():
                        range_idx = 6
                        ref1 = 10.0**(0.1*radarvariables['reflectivity'].data[:, range_idx])
                        ref2 = 10.0**(0.1*radarvariables['reflectivity'].data[:, range_idx+1])
                        reflectivity = 0.5*(ref1+ref2)
                    basetime = datetime.datetime.strptime(radarvariables['time'].units.decode(),
                                                          'seconds since %Y-%m-%d %H:%M:%S 0:00')
                    delta = radarvariables['time'].data
                    deltatime = pd.to_timedelta(np.round(delta), unit='s')
                    time = basetime + deltatime
                    self.data = pd.DataFrame(reflectivity, index=time,
                                             columns=['reflectivity'],
                                             dtype=np.float64)
        self.finish_init(dt_start, dt_end)

    def good_data(self):
        """Return useful data with filters and corrections applied."""
        if self.stored_good_data is not None:
            return self.stored_good_data
        data = self.data.copy()
        data.index = self.data.index + self.time_lag
        ns1min = 1*60*1000000000
        time = pd.DatetimeIndex((np.round(data.index.astype(np.int64)/ns1min))*ns1min)
        data.index = time
        return data

    def z(self, rule=None, varinterval=True):
        """Reflectivity time series"""
        grp = self.grouped(rule=rule, varinterval=varinterval)
        z = grp.mean()
        zs = z[z.columns[0]]
        zs.name = self. name + ' reflectivity'
        zs.index.name = 'datetime'
        return zs

    @property
    def time_lag(self):
        return self._time_lag

    @time_lag.setter
    def time_lag(self, lag):
        self._time_lag = lag