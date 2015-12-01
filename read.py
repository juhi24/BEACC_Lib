# -*- coding: utf-8 -*-
"""
tools for reading and working with baecc data
@author: Jussi Tiira
"""
import time
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import datetime
import pandas as pd
import numpy as np
import linecache
import copy
import fit
from scipy import stats, io
from scipy.optimize import fmin, minimize
import pickle
import warnings
from glob import glob
import locale
import netCDF4 as nc
import hashlib

locale.setlocale(locale.LC_ALL, 'en_GB.UTF-8')

# general configuration
DEBUG = False

# CONFIG default paths
DATA_DIR = '../DATA'
H5_FILE = 'baecc.h5'
H5_PATH = os.path.join(DATA_DIR, H5_FILE)
PIPV_SUBPATH = 'PIP/a_Velocity_Tables/004%s/*2.dat'
DSD_SUBPATH = 'PIP/a_DSD_Tables/004%s*.dat'
P200_SUBPATH = 'Pluvio200/pluvio200_??_%s*.txt'
P400_SUBPATH = 'Pluvio400/pluvio400_??_%s*.txt'
RADAR_SUBPATH = 'Radar/%s/tmp%s*M1.a1.%s.*'
RESULTS_DIR = '../results'
CACHE_DIR = 'cache'
MSGTLD = '.msg'
PICKLETLD = '.pkl'

ns1min = 1.0*60.0*1000000000.0


def hash_dict(d):
    return fingerprint(str(sorted(d.items())))


def fingerprint(string):
    return hashlib.sha256(string.encode('utf-8')).hexdigest()[-12:]


def cdf_to_series(path, vname, tname='time'):
    cdf = nc.Dataset(path)
    h = cdf.variables[vname]
    t = cdf.variables[tname]
    dt = nc.num2date(t[:], t.units)
    hs = pd.Series(h[:], index=dt)
    hs.name = vname
    return hs


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def file_shorter_than(fname, limit):
    with open(fname) as f:
        for i, l in enumerate(f):
            if i+2>limit:
                return False
    return True


def datenum2datetime(matlab_datenum):
    """Convert MATLAB datenum to datetime."""
    return datetime.datetime.fromordinal(int(matlab_datenum)) + datetime.timedelta(days=matlab_datenum % 1) - datetime.timedelta(days=366)


def ensure_dir(directory):
    """Make sure the directory exists. If not, create it."""
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory


def merge_series(s1, s2, **kwargs):
    """Merge pandas Series and/or DataFrames on index"""
    return pd.merge(pd.DataFrame(s1), pd.DataFrame(s2),
                    left_index=True, right_index=True, **kwargs)


def merge_multiseries(s1, s2, *series, **kwargs):
    """Merge multiple pandas Series and/or DataFrames on index"""
    s_all = merge_series(s1, s2, **kwargs)
    if len(series) > 0:
        for s in series:
            s_all = merge_series(s_all, s, **kwargs)
    return s_all


def kde(x, y):
    values = np.vstack((x, y))
    return stats.gaussian_kde(values)


def filter_outlier(X, Y, Z, data, xname='x', yname='y', frac=0.5,
                   bin_limit_multiplier=0.05):
    filtered = pd.DataFrame()
    std = []
    x = X[0, :]
    y = Y[:, 0]
    xwidth = (x[-1]-x[0])/len(x)    # TODO: check if correct
    #print('xw: %s' % xwidth)
    xlims = np.append(x-0.5*xwidth, x[-1]+0.5*xwidth)
    xmin = xlims[:-1]
    xmax = xlims[1:]
    ymin = []
    ymax = []
    data_count = data.count()[0]
    for i in range(0, Z.shape[1]):  # loop through bins
        z = Z[:, i]
        z_lim = z.max()*frac
        y_fltr = y[z > z_lim]       # FWHM when frac=0.5
        if y_fltr.size == 0:
            ymin.append(None)
            ymax.append(None)
            continue
        ymin.append(y_fltr[0])
        ymax.append(y_fltr[-1])
        bin_data = bindata(xmin[i], xmax[i], data, ymin=ymin[-1], ymax=ymax[-1],
                           xname=xname, yname=yname)
        bin_count = bin_data.count()[0]
        count_limit = bin_limit_multiplier*xwidth*data_count
        if bin_count < count_limit:
            ymin[-1] = None
            ymax[-1] = None
            continue
        std.append(bin_data[yname].std())
        filtered = filtered.append(bin_data)
    # return filtered data, stds and half width at frac*kde_max
    return filtered, np.array(std), xlims, ymin, ymax


def bindata(xmin, xmax, data, ymin=None, ymax=None, xname='x', yname='y'):
    """Return data that falls into given x bin."""
    cond = '%s > %s and %s < %s' % (xname, xmin, xname, xmax)
    if ymin is not None and ymax is not None:
        cond += ' and %s > %s and %s < %s' % (yname, ymin, yname, ymax)
    return data.query(cond)


def datafilelist(subpath, datadir=DATA_DIR):
    return glob(os.path.join(datadir, subpath))


def datafilelistloop(subpath, dtstrlist, datadir=DATA_DIR):
    listout = []
    for dtstr in dtstrlist:
        listout.extend(datafilelist(subpath % dtstr, datadir=datadir))
    return listout


def batch_import(dtstrlist, datadir=DATA_DIR, radar=False):
    """Read ASCII data according to a datestring pattern."""
    pipv_files = datafilelistloop(PIPV_SUBPATH, dtstrlist, datadir=datadir)
    dsd_files = datafilelistloop(DSD_SUBPATH, dtstrlist, datadir=datadir)
    pluvio200_files = datafilelistloop(P200_SUBPATH, dtstrlist, datadir=datadir)
    pluvio400_files = datafilelistloop(P400_SUBPATH, dtstrlist, datadir=datadir)
    if radar:
        xsacr_files = datafilelistloop(RADAR_SUBPATH, [('XSACR', 'xsacr', dtstr) for dtstr in dtstrlist],
                                   datadir=datadir)
        kasacr_files = datafilelistloop(RADAR_SUBPATH, [('KASACR', 'kasacr', dtstr) for dtstr in dtstrlist],
                                    datadir=datadir)
        kazr_files = datafilelistloop(RADAR_SUBPATH, [('KAZR', 'kazrge', dtstr) for dtstr in dtstrlist],
                                  datadir=datadir)
        mwacr_files = datafilelistloop(RADAR_SUBPATH, [('MWACR', 'mwacr', dtstr) for dtstr in dtstrlist],
                                   datadir=datadir)
    pluvio200 = Pluvio(pluvio200_files)
    pluvio400 = Pluvio(pluvio400_files)
    pipv = PipV(pipv_files)
    dsd = PipDSD(dsd_files)
    if radar:
        xsacr = Radar(xsacr_files)
        kasacr = Radar(kasacr_files)
        kazr = Radar(kazr_files)
        mwacr = Radar(mwacr_files)
        return {'vel': pipv, 'dsd': dsd, 'pluvio200': pluvio200,
                'pluvio400': pluvio400, 'xsacr': xsacr, 'kasacr': kasacr,
                'kazr': kazr, 'mwacr': mwacr}
    return {'vel': pipv, 'dsd': dsd, 'pluvio200': pluvio200,
            'pluvio400': pluvio400}


def batch_create_hdf(instrdict=None, datadir=DATA_DIR, outname=H5_FILE,
                     dtstrlist=('20140[2-3]??')):
    """Read ASCII data and export to hdf."""
    if instrdict is None:
        instrdict = {PIPV_SUBPATH: PipV,
                     DSD_SUBPATH: PipDSD,
                     P200_SUBPATH: Pluvio,
                     P400_SUBPATH: Pluvio}
    hdf_file = os.path.join(datadir, outname)
    for key in instrdict:
        instr = instrdict[key].from_raw(dtstrlist, subpath=key, datadir=datadir)
        instr.to_hdf(filename=hdf_file)


class Cacher:
    """common methods to use msg cache"""
    def __init__(self, use_cache=True, storefilename='store.h5',
                 parent=None):
        self.use_cache = use_cache
        self.storefilename = storefilename
        self.parent = parent

    def msger(self, name, func, **kwargs):
        """Read from msgpack if caching is in use."""
        if self.use_cache:
            return self.msg_io(name, func, **kwargs)
        return func(**kwargs)

    def pickler(self, name, func, **kwargs):
        if self.use_cache:
            return self.pkl_io(name, func, **kwargs)
        return func(**kwargs)

    def cache_dir(self):
        """Return full path to cache directory."""
        if self.parent is None:
            return os.path.join(CACHE_DIR, self.fingerprint())
        return os.path.join(CACHE_DIR, self.parent.fingerprint(),
                            self.fingerprint())

    def store_path(self):
        """Return full path to hdf store file."""
        return os.path.join(self.cache_dir(), self.storefilename)

    def store_read(self, tablename, default_value=None, nocache_value=None):
        """Read from hdf store if using caching."""
        if self.use_cache:
            ensure_dir(self.cache_dir())
            try:
                with pd.HDFStore(self.store_path()) as store:
                    return store.get(tablename)
            except KeyError as err:
                warnings.warn("KeyError: {0} Using default value.".format(err))
                return default_value
        return nocache_value

    def store_write(self, tablename, data):
        """Write data to hdf store."""
        ensure_dir(self.cache_dir())
        with pd.HDFStore(self.store_path()) as store:
                store[tablename] = data

    def msg_io(self, name, func, **kwargs):
        """Read data from msgpack. If not available, calculate and store."""
        cd = self.cache_dir()
        msgpath = os.path.join(cd, name + MSGTLD)
        if os.path.exists(msgpath):
            data = pd.read_msgpack(msgpath)
        else:
            ensure_dir(cd)
            data = func(**kwargs)
            data.to_msgpack(msgpath)
        return data

    def pkl_io(self, name, func, **kwargs):
        cd = self.cache_dir()
        pklpath = os.path.join(cd, name + PICKLETLD)
        if os.path.exists(pklpath):
            with open(pklpath, 'rb') as cachefile:
                data = pickle.load(cachefile)
        else:
            ensure_dir(cd)
            data = func(**kwargs)
            with open(pklpath, 'wb') as cachefile:
                pickle.dump(data, cachefile, pickle.HIGHEST_PROTOCOL)
        return data

    def clear_cache(self, extra_files=None):
        """Remove cache files used by the Cacher object."""
        store = self.store_path()
        filelist = []
        if os.path.exists(store):
            filelist.append(store)
        if extra_files is not None:
            filelist.extend(extra_files)
        for f in filelist:
            os.remove(f)

    def fingerprint(self):
        """state-aware object identifier, immutable between sessions"""
        pass


class PrecipMeasurer:
    """parent for classes with precipitation measurements
    Either amount or acc (or both) methods should be overridden."""
    def __init__(self):
        pass

    def amount(self, **kwargs):
        """timestep precipitation in mm"""
        am = self.acc(**kwargs).diff()
        am.name = 'amount'
        return am

    def acc(self, **kwargs):
        """precipitation accumulation in mm"""
        acc = self.amount(**kwargs).cumsum()
        acc.name = 'accumulation'
        return acc

    def intensity(self, tdelta=None, **kwargs):
        """precipitation intensity in mm/h"""
        r = self.amount(**kwargs)
        if tdelta is None:
            tdelta = r.index.freq.delta
        frac = tdelta.apply(lambda t: np.timedelta64(1,'h')/t)
        intensity = frac * r
        intensity.name = 'intensity'
        return intensity


class InstrumentData(Cacher):
    """Parent for instrument data classes."""
    # TODO: Separate read_csv and __init__
    def __init__(self, filenames, hdf_table=None, use_cache=True):
        """Read from either ASCII data file or hdf5."""
        self.filenames = filenames
        self.data = pd.DataFrame()
        # if filtered data needed often, keep in memory
        self.stored_good_data = None    # set to None to disable
        if hdf_table is not None:
            self.name = hdf_table
            self.data = self.data.append(pd.read_hdf(filenames[0], hdf_table))
        Cacher.__init__(self, use_cache=use_cache)

    def __add__(self, other):
        combined = copy.deepcopy(self)
        combined.data = pd.concat([self.data, other.data])
        combined.clear_cache()
        return combined

    @classmethod
    def from_raw(cls, dtstrlist, subpath='', datadir=DATA_DIR):
        filelist = datafilelistloop(subpath, dtstrlist, datadir=datadir)
        return cls(filelist)

    def fingerprint(self):
        return fingerprint(str(self.data))

    def finish_init(self, dt_start, dt_end):
        """Sort and name index, cut time span."""
        self.data.sort_index(inplace=True)
        self.data.index.names = ['datetime']
        self.set_span(dt_start, dt_end)
        Cacher.__init__(self, storefilename=self.name + '.h5')

    def store_good_data(self, **kwargs):
        """Store good data to memory (to bypass recalculation of filters)."""
        self.stored_good_data = self.good_data(**kwargs)

    def parse_datetime(self):
        """Parse timestamps in data files. Used by class constructor."""
        pass

    def good_data(self):
        """Return useful data with filters and corrections applied."""
        if self.stored_good_data is not None:
            return self.stored_good_data
        return self.data

    def to_hdf(self, filename='../DATA/baecc.h5'):
        """Save object in hdf5 format."""
        self.data.to_hdf(filename, self.name, format='table', append=True)

    def between_datetime(self, date_start, date_end, inplace=False):
        """Limit the time span of data."""
        if inplace:
            instr = self
        else:
            instr = copy.deepcopy(self)
        instr.set_span(date_start, date_end)
        return instr

    def set_span(self, dt_start, dt_end):
        """Limit data to given time interval."""
        for dt in [dt_start, dt_end]:
            dt = pd.datetools.to_datetime(dt)
        self.data = self.data[dt_start:dt_end].copy()

    def grouped(self, varinterval=False, rule=None, col=None):
        if rule is None:
            rule = self.rule
        data = self.good_data()
        if col is not None:
            data = pd.DataFrame(data[col])
        if varinterval:
            grpd_data = pd.merge(data, rule, left_index=True, right_index=True)
            return grpd_data.groupby('group')
        return data.groupby(pd.Grouper(freq=rule, closed='right',
                                                   label='right'))


class Radar(InstrumentData):
    """Radar reflectivity at lowest level"""
    def __init__(self, filenames, dt_start=None, dt_end=None, **kwargs):
        """Create vertical pointing Radar object using data from various radar
        modes"""
        print('Reading Radar data...')
        self._time_lag = pd.to_timedelta(0.0, unit='s')
        InstrumentData.__init__(self, filenames, **kwargs)
        if self.data.empty and filenames:
            self.name = (os.path.basename(os.path.dirname(self.filenames[0])))
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
        data = copy.deepcopy(self.data)
        data.index = self.data.index + self.time_lag
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


class Pluvio(InstrumentData, PrecipMeasurer):
    """Pluviometer data handling"""
    def __init__(self, filenames, dt_start=None, dt_end=None, **kwargs):
        """Create a Pluvio object using data from a list of files."""
        print('Reading pluviometer data...')
        InstrumentData.__init__(self, filenames, **kwargs)
        self.bias = 0
        self._shift_periods = 0
        self._shift_freq = '1min'
        self.lwc = None
        self.use_bucket = False
        self._varinterval = True
        self.n_combined_intervals = 1
        col_suffix = 'nrt'
        self.amount_col = 'acc_' + col_suffix
        self.bucket_col = 'bucket_' + col_suffix
        if self.data.empty:
            self.name = os.path.basename(os.path.dirname(self.filenames[0])).lower()
            self.col_description = ['date string',
                                    'intensity RT [mm h]',
                                    'accumulated RT/NRT [mm]',
                                    'accumulated NRT [mm]',
                                    'accumulated total NRT [mm]',
                                    'bucket RT [mm]',
                                    'bucket NRT [mm]',
                                    'temperature load cell [degC]',
                                    'heating status',
                                    'status',
                                    'temperature electronics unit',
                                    'supply voltage',
                                    'ice rim temperature']
            col_abbr = ['datestr',
                        'group',#'i_rt', #don't know why, but needs this column to be exsisting before
                        'acc_rt',
                        'acc_nrt',
                        'acc_tot_nrt',
                        'bucket_rt',
                        'bucket_nrt',
                        't_load',
                        'heating',
                        'status',
                        't_elec',
                        'volt',
                        't_rim']
            for filename in filenames:
                print('.', end='')
                #num_lines = file_len(filename)
                self.current_file = filename
                try:
                    self.data = self.data.append(pd.read_csv(filename, sep=';',
                                     names=col_abbr,
                                     skip_blank_lines=True,
                                     error_bad_lines=False,
                                     warn_bad_lines=True,
                                     parse_dates={'datetime':['datestr']},
                                     date_parser=self.parse_datetime,
                                     index_col='datetime',
                                     verbose=DEBUG))
                except NotImplementedError as err:
                    print('\n%s: %s' % (filename, format(err)))
            print()
            #self.data.drop(['i_rt'], 1, inplace=True) # crap format
        self.buffer = pd.datetools.timedelta(0)
        self.finish_init(dt_start, dt_end)
        self.data['group'] = self.data.acc_nrt.astype(bool).astype(int).cumsum().shift(1).fillna(0)
        self.data['heating'] = self.data.astype(float)  # sometimes this are interpreted as int
        self.data['status'] = self.data.astype(float)   # sometimes this are interpreted as int

    @property
    def varinterval(self):
        return self._varinterval

    @varinterval.setter
    def varinterval(self, varinterval):
        self._varinterval = varinterval
        self.use_bucket = not varinterval

    @property
    def shift_periods(self):
        return self._shift_periods

    @shift_periods.setter
    def shift_periods(self, shift_periods):
        self._shift_periods = shift_periods
        if self.use_bucket:
            self.noprecip_bias(self.lwc, inplace=True)

    @property
    def shift_freq(self):
        return self._shift_freq

    @shift_freq.setter
    def shift_freq(self, shift_freq):
        self._shift_freq = shift_freq
        if self.use_bucket:
            self.noprecip_bias(self.lwc, inplace=True)

    @classmethod
    def from_raw(cls, *args, subpath=P200_SUBPATH, **kwargs):
        return super().from_raw(*args, subpath=subpath, **kwargs)

    def fingerprint(self):
        identifiers = [self.name, self.shift_periods, self.shift_freq,
                       self.varinterval]
        if self.varinterval:
            identifiers.extend([self.n_combined_intervals])
        idstr = ''.join(tuple(map(str, identifiers)))
        return fingerprint(super().fingerprint() + idstr)

    def parse_datetime(self, datestr, include_sec=False):
        datestr = str(int(datestr))
        t = time.strptime(datestr, '%Y%m%d%H%M%S')
        if include_sec:
            t_end = 6
        else:
            t_end = 5
        #return datetime.datetime(*t[:t_end], tzinfo=datetime.timezone.utc)
        return datetime.datetime(*t[:t_end])

    def good_data(self):
        if self.stored_good_data is not None:
            return self.stored_good_data
        data = copy.deepcopy(self.data)
        swap_date = pd.datetime(2014, 5, 16, 8, 0, 0)#, tzinfo=datetime.timezone.utc)
        swap_date2 = pd.datetime(2014, 8, 31, 8, 0, 0)#, tzinfo=datetime.timezone.utc ) # TODO put correct switch date
        if self.data.index[-1] > swap_date and self.data.index[-1] < swap_date2:
            precip_cols = ['acc_rt', 'acc_nrt', 'acc_tot_nrt', 'bucket_rt',
                           'bucket_nrt']
            if self.name == 'pluvio200':
                correction = 2
            elif self.name == 'pluvio400':
                correction = 0.5
            for col in precip_cols:
                data[col] = self.data[col]*correction
        return data

    def set_span(self, dt_start, dt_end):
        """Set time span with a buffer for timeshift."""
        if dt_start is None or dt_end is None:
            super().set_span(dt_start, dt_end)
            return
        for dt in [dt_start, dt_end]:
            dt = pd.datetools.to_datetime(dt)
        self.buffer = pd.datetools.timedelta(hours=2)
        if dt_start is None or dt_end is None:
            self.buffer = pd.datetools.timedelta(0)
        elif dt_start-self.buffer < self.data.index[0] or dt_end+self.buffer > self.data.index[-1]:
            self.buffer = pd.datetools.timedelta(0)
        self.data = self.data[dt_start-self.buffer:dt_end+self.buffer]

    def timeshift(self):
        """Return timeshift as timedelta."""
        if self.shift_periods == 0:
            return pd.datetools.timedelta(0)
        return self.shift_periods*pd.datetools.to_offset(self.shift_freq)

    def dt_start(self):
        return self.data.index[0] + self.buffer

    def dt_end(self):
        return self.data.index[-1] - self.buffer

    def shift_reset(self):
        """Reset time shift."""
        self.shift_periods = 0
        self.shift_freq = '1min'

    def constinterval_amount(self, rule='1H', upsample=True, **kwargs):
        """Calculate precipitation amount"""
        if upsample:
            acc_1min = self.constinterval_acc('1min', **kwargs)
        else:
            acc_1min = self.acc_raw()
        r = acc_1min.diff().resample(rule, how=np.sum, closed='right',
                                     label='right')
        if not upsample:
            return r.fillna(0)
        t_r0 = r.index[0]
        r[0] = acc_1min[t_r0]-acc_1min[0]
        return r

    def amount(self, crop=True, shift=True, **bucketkwargs):
        if not self.varinterval:
            return self.constinterval_amount(shift=shift, **bucketkwargs)
        am = self.good_data()[self.amount_col]
        n = self.n_combined_intervals
        am = pd.stats.moments.rolling_sum(am[am > 0], window=n).iloc[n-1::n]
        if shift:
            am = am.tshift(periods=self.shift_periods, freq=self.shift_freq)
        if crop:
            am = am[self.dt_start():self.dt_end()]
        return am

    def intensity(self, **kwargs):
        if self.varinterval:
            return super().intensity(tdelta=self.tdelta(), **kwargs)
        return super().intensity(tdelta=None, **kwargs)

    def constinterval_acc(self, rule='1H', interpolate=True, unbias=True,
                          shift=True, filter_evap=True):
        """Resample unbiased accumulated precipitation in mm."""
        accum = self.acc_raw().asfreq('1min')
        if interpolate:
            accum.interpolate(method='time', inplace=True)
        else:
            accum.fillna(method='bfill', inplace=True)
        if shift:
            accum = accum.tshift(periods=self.shift_periods,
                                 freq=self.shift_freq)
        if unbias:
            accum -= self.bias
        accum = accum[self.dt_start():self.dt_end()]
        if filter_evap:
            amount = accum.diff()
            evap = amount[amount < 0]
            evap_accum = evap.reindex(accum.index).fillna(0).cumsum()
            #evap_accum.plot()
            #accum.plot()
            accum -= evap_accum     # accumulation should never drop
        return accum.resample(rule, how='last', closed='right', label='right')

    def acc_raw(self):
        """accumulation from raw data"""
        return self.good_data()[self.bucket_col]-self.good_data()[self.bucket_col][0]

    def noprecip_bias(self, lwc, inplace=False):
        """Calculate accumulated bias using LWC."""
        self.lwc = lwc
        accum = self.acc(rule='1min', shift=True, unbias=False,
                         filter_evap=False)
        lwc_filled = lwc.reindex(accum.index).fillna(0)
        bias_amount = accum.diff().fillna(0)[lwc_filled == 0]
        #bias_amount[bias_amount > 0] = 0
        bias_acc = bias_amount.cumsum()
        if bias_acc.empty:
            if inplace:
                self.bias = 0
            return 0
        bias_acc_filled = bias_acc.reindex(accum.index).asfreq('1min').fillna(method='bfill').fillna(method='ffill')
        if inplace:
            self.bias = bias_acc_filled
        return bias_acc_filled

    def tdelta(self):
        """Return lengths of timesteps as Series of timedeltas."""
        a = self.amount(crop=False)
        delta = pd.Series(a.index.to_pydatetime(), index=a.index).diff()
        longest_delta = pd.datetools.timedelta(hours=1)
        delta[delta > longest_delta] = longest_delta
        delta = pd.to_timedelta(delta)
        delta.name = 'tdelta'
        return delta[self.dt_start():self.dt_end()].fillna(longest_delta)

    def start_time(self):
        """timestep starting times"""
        tdelta = self.tdelta()
        td = pd.DataFrame(tdelta)
        td['end'] = td.index
        start = td.end.sub(tdelta)
        start.name = 'start'
        return start

    def grouper(self, shift=True):
        """data group names (timestamp) for each data timestamp"""
        ticks = self.good_data()[self.amount_col].astype(bool)
        if shift:
            ticks = ticks.tshift(periods=self.shift_periods,
                                 freq=self.shift_freq)
        ticktime = self.amount().index
        dtgroups = pd.Series(ticktime, index=ticktime).reindex(ticks.index).bfill()[self.dt_start():self.dt_end()]
        dtgroups.name = 'group'
        last_index = self.tdelta().index[-1]
        return pd.DataFrame(dtgroups[dtgroups.notnull()])[:last_index]

    def groupby_interval(self, data):
        """Group data by integration time intervals."""
        return merge_series(data, self.grouper()).groupby('group')

class PipDSD(InstrumentData):
    """PIP particle size distribution data handling"""
    def __init__(self, filenames, dt_start=None, dt_end=None, d_bin=.25, **kwargs):
        """Create a PipDSD object using data from a list of PIP DSD table files."""
        print('Reading PIP PSD data...')
        InstrumentData.__init__(self, filenames, **kwargs)
        self.d_bin = d_bin
        self.name = 'pip_dsd'
        common_csv_kws = {'skiprows': 8,
                          'header': 3,
                          'parse_dates': {'datetime':['hr_d', 'min_d']},
                          'date_parser': self.parse_datetime,
                          'index_col': 'datetime',
                          'verbose': DEBUG}
        if self.data.empty:
            for filename in filenames:
                if DEBUG:
                    print(filename)
                else:
                    print('.', end='')
                if file_shorter_than(filename, 14):
                    # file has no data
                    continue
                self.current_file = filename
                # File format changed
                if int(os.path.split(filename)[1][3:11]) > 20140831: # TODO fixme
                    self.data = self.data.append(pd.read_csv(filename,
                                                             engine='python',
                                                             sep='\t',
                                                             skipfooter=1,
                                                             **common_csv_kws))
                else:
                    self.data = self.data.append(pd.read_csv(filename,
                                                             delim_whitespace=True,
                                                             **common_csv_kws))
            print()
            avg = pd.read_csv(self.current_file, delim_whitespace=True, skiprows=8,
                                   nrows=1, header=None).drop([0, 1, 2, 3, 4], axis=1)
            #self.num_d = self.data[['Num_d']]
            # 1st size bin is crap data, last sometimes nans
            self.d_bin = float(linecache.getline(self.current_file,
                                                 11).split()[5])
            self.data.drop(['day_time', 'Num_d', 'Bin_cen'], 1,
                           inplace=True)
            self.data.columns = pd.Index([float(i) for i in self.data.columns])
            self.data.sort_index(axis=1)
            avg.columns = self.data.columns
            self.avg = avg.T[0]
            self.avg.name = 'dsd_avg'
            self.data = self.data.astype(float)
        self.data.drop_duplicates(inplace=True)
        self.data = self.data.resample('1min').fillna(0)
        self.finish_init(dt_start, dt_end)

    @classmethod
    def from_raw(cls, *args, subpath=DSD_SUBPATH, **kwargs):
        return super().from_raw(*args, subpath=subpath, **kwargs)

    def parse_datetime(self, hh, mm):
        dateline = linecache.getline(self.current_file, 6)
        datearr = [int(x) for x in dateline.split()]
        d = datetime.date(*datearr)
        t = datetime.time(int(hh), int(mm))
        return datetime.datetime.combine(d, t)#.replace(tzinfo=datetime.timezone.utc)

    def bin_cen(self):
        """Return array of bin centers."""
        return self.good_data().columns.values

    def n(self, d, **kwargs):
        """number concentrations for given diameter"""
        n = self.psd(col=d, **kwargs)
        ns = n[n.columns[0]] # convert to Series
        ns.name = 'N_' + str(d)
        return ns

    def psd(self, rule='1min', varinterval=False, **kwargs):
        grp = self.grouped(rule=rule, varinterval=varinterval, **kwargs)
        n = grp.mean()
        return n

    def plot(self, img=True, **kwargs):
        """Plot particle size distribution over time."""
        if img:
            plt.matshow(self.good_data(**kwargs).transpose(), norm=LogNorm(),
                        origin='lower')
        else:
            plt.pcolor(self.good_data(**kwargs).transpose(), norm=LogNorm())
        plt.colorbar()
        plt.title('PIP DSD')
        plt.xlabel('time (UTC) BROKEN')
        plt.ylabel('D (mm) BROKEN')

    def filter_cats_and_dogs(self, data=None, window=5):
        """a rolling window filter for isolated data points"""
        if data is None:
            data = self.data
        # Any datapoint after <window-1> bins of zeros will be flagged
        is_dog = pd.rolling_count(data.mask(data == 0).T, window).T == 1
        is_dog.ix[:, :window] = False # unflag <window> first columns
        is_dog[is_dog == False] = np.nan # False --> NaN, True --> 1
        # In a time interval flag anything that's bigger than any previously
        # flagged bin (forward fill).
        is_dog = is_dog.ffill(axis=1).fillna(False)
        is_dog = is_dog.astype(np.bool)     # convert back to boolean
        filtered = copy.deepcopy(data)
        filtered[is_dog] = 0                # apply filter
        return filtered

    def good_data(self, filter_large=True, **kwargs):
        if self.stored_good_data is not None:
            return self.stored_good_data
        gain_correction = 2
        update_date = pd.datetime(2014, 11, 25, 8, 0, 0)
        if self.data.index[-1] > update_date:
            gain_correction = 1
        data = self.data
        if filter_large:
            data = self.filter_cats_and_dogs(data=data, **kwargs)
        bin_cen = self.data.columns.values
        too_small = bin_cen[bin_cen < 0.3]
        return gain_correction*data.drop(too_small, 1)


class PipV(InstrumentData):
    """PIP particle velocity and diameter data handling"""
    def __init__(self, filenames, dt_start=None, dt_end=None, **kwargs):
        """Create a PipV object using data from a list of PIP velocity table files."""
        print('Reading PIP particle velocity data...')
        InstrumentData.__init__(self, filenames, **kwargs)
        self.name = 'pip_vel'
        self.dmin = 0.375   # shortest diameter where data is good
        self._fits = pd.DataFrame()
        # num=511 --> binwidth 0.05
        if DEBUG:
            num = 103
        else:
            num = 409
        self.dbins = np.linspace(self.dmin, 25.875, num=num)
        self._std = pd.DataFrame(columns=self.dbins)
        # half width at fraction of maximum
        self._hwfm = pd.DataFrame(columns=self.dbins)
        self.default_fit = fit.PolFit
        self.flip = False
        if self.data.empty:
            for filename in filenames:
                print('.', end='')
                self.current_file = filename
                if int(filename[-23:-15]) > 20140831: # TODO fixme
                    newdata = pd.read_csv(filename,
                                          engine='python', sep='\t',
                                          skipinitialspace=True, skiprows=8,
                                          skip_footer=1,
                                          parse_dates={'datetime':['minute_p']},
                                          date_parser=self.parse_datetime,
                                          verbose=DEBUG)
                    newdata = newdata[newdata['RecNum'] > -99]
                else:
                    newdata = pd.read_csv(filename,
                                          delim_whitespace=True,
                                          skiprows=8,
                                          parse_dates={'datetime':['minute_p']},
                                          date_parser=self.parse_datetime,
                                          verbose=DEBUG)
                    newdata = newdata[newdata['RecNum']>-99]
                if not newdata.empty:
                    newdata.rename_axis({'vel_v_1': 'vel_v',
                                         'vel_h_1': 'vel_h',
                                         'vel_v_2': 'vel_v',
                                         'vel_h_2': 'vel_h'},
                                        axis=1, inplace=True)
                    newdata.set_index(['datetime', 'Part_ID', 'RecNum'],
                                      inplace=True)
                    g = newdata.groupby(level=['datetime', 'Part_ID'])
                    newdata = g.mean()
                    self.data = self.data.append(newdata)
            print()
            if len(self.data.index):
                self.data = self.data[self.data.vel_v.notnull()]
            self.data.reset_index(level=1, inplace=True)
            self.data = self.data.astype(float)
        self.finish_init(dt_start, dt_end)

    @property
    def rule(self):
        if self.fits.empty:
            return None
        #return self.fits.index.freq.freqstr
        return self.fits.index.freqstr

    @property
    def binwidth(self):
        d = self.dbins
        return (d[-1]-d[0])/(len(d)-1)

    @property
    def fits(self):
        return self.store_read('fits', default_value=pd.DataFrame(),
                               nocache_value=self._fits)

    @fits.setter
    def fits(self, fits):
        if self.use_cache:
            self.store_write('fits', fits)
        else:
            self._fits = fits

    @property
    def std(self):
        return self.store_read('std',
                               default_value=pd.DataFrame(columns=self.dbins),
                               nocache_value=self._std)

    @std.setter
    def std(self, std):
        if self.use_cache:
            self.store_write('std', std)
        else:
            self._std = std

    @property
    def hwfm(self):
        return self.store_read('hwfm',
                               default_value=pd.DataFrame(columns=self.dbins),
                               nocache_value=self._hwfm)

    @hwfm.setter
    def hwfm(self, hwfm):
        if self.use_cache:
            self.store_write('hwfm', hwfm)
        else:
            self._hwfm = hwfm

    @classmethod
    def from_raw(cls, *args, subpath=PIPV_SUBPATH, **kwargs):
        return super().from_raw(*args, subpath=subpath, **kwargs)

    def fingerprint(self):
        identifiers = (self.flip, self.dbins)
        idstr = ''.join(tuple(map(str, identifiers)))
        return fingerprint(super().fingerprint() + idstr)

    def v(self, d, fitclass=None, varinterval=True, rule=None):
        """velocities according to fits for given diameter"""
        if fitclass is None:
            fitclass = self.default_fit
        if rule is None:
            rule = self.rule
        if self.fits.empty:
            self.find_fits(rule, fitclass=fitclass, varinterval=varinterval)
        elif not varinterval:
            if pd.datetools.to_offset(rule) != self.fits.index.freq:
                print('different sampling freq')
                self.find_fits(rule, fitclass=fitclass,
                               varinterval=varinterval)
        v = []
        for vfit in self.fits[fitclass.name].values:
            v.append(vfit.func(d))
        return pd.Series(v, index=self.fits.index, name='v')

    def lwc(self, rule='1min'):
        """liquid water content"""
        d3 = self.good_data().Wad_Dia**3
        return d3.resample(rule, how=np.sum, closed='right', label='right')

    def parse_datetime(self, mm):
        datestr = self.current_file.split('/')[-1].split('_')[0]
        yr = int(datestr[3:7])
        mo = int(datestr[7:9])
        dd = int(datestr[9:11])
        hh = int(datestr[11:13])
        return datetime.datetime(yr, mo, dd, hh, int(mm))#,
                                 #tzinfo=datetime.timezone.utc)

    def good_data(self):
        if self.stored_good_data is not None:
            return self.stored_good_data
        return self.data[self.data.Wad_Dia > self.dmin]

    def filter_outlier(self, data, frac=0.5, flip=False):
        """Filter outliers using KDE"""
        if data is None:
            data = self.good_data()
        if flip:
            X, Y, Z = self.kde_grid(data)
            return filter_outlier(Y.T, X.T, Z.T, data, xname='vel_v',
                                  frac=frac, yname='Wad_Dia')
        return filter_outlier(*self.kde_grid(data), data=data,
                              xname='Wad_Dia', yname='vel_v', frac=frac)

    def frac_larger(self, d):
        """Return fraction of particles that are larger than d."""
        vdata = self.good_data()
        return vdata[vdata.Wad_Dia > d].vel_v.count()/vdata[vdata.Wad_Dia < d].vel_v.count()

    def d_cut(self, frac=0.05, d_guess=2):
        """Return d for which given fraction of particles are larger."""
        dcost = lambda d: abs(self.frac_larger(d[0])-frac)
        return fmin(dcost, d_guess)[0]

    def find_fit(self, fitclass=None, data=None, use_kde_peak=False,
                 cut_d=False, frac=0.5, use_curve_fit=True, bin_num_min=5,
                 filter_outliers=True, name=None, try_flip=True,
                 plot_flip=False, force_flip=False, cut_kws={}, **kwargs):
        """Find and store a fit for either raw data or kde."""
        # TODO: clean this mess
        def too_few_particles(use_curve_fit, use_kde_peak):
            print('Too few particles.')
            return False, False
        std = pd.DataFrame()
        hwfm = pd.DataFrame()
        if force_flip:
            try_flip = True
        if data is None:
            data = self.good_data()
        origdata = copy.deepcopy(data)  # do not rewrite
        if fitclass is None:
            fitclass = self.default_fit
        vfit = fitclass()
        partcount = data.count()[0]
        # increased from 5 to 10
        if partcount < 10 and (use_curve_fit or use_kde_peak):
            use_curve_fit, use_kde_peak = too_few_particles(use_curve_fit,
                                                            use_kde_peak)
        elif filter_outliers:
            data, stdarr, xlims, ymin, ymax = self.filter_outlier(data=data,
                                                             frac=frac)
            vfit.fltr_upper_x = xlims
            vfit.fltr_lower_x = xlims
            vfit.fltr_upper_y = ymax + [ymax[-1]]
            vfit.fltr_lower_y = ymin + [ymin[-1]]
            # TODO: Make a proper conditional to include this when needen, if needed
            if False:
                datao = self.filter_outlier(data=data, frac=frac, flip=True)[0]
                fltrcount = datao.count()[0]
                if fltrcount < 2 and (use_curve_fit or use_kde_peak):
                    use_curve_fit, use_kde_peak = too_few_particles(use_curve_fit,
                                                                    use_kde_peak)
            if name is not None:
                stdarr.resize(self.dbins.size, refcheck=False)
                std = pd.DataFrame(pd.Series(stdarr,
                                             index=self.dbins, name=name)).T
                # TODO: empty hwfm
                hwfm = pd.DataFrame(pd.Series(index=self.dbins, name=name)).T
                for df in [std, hwfm]:
                    df.index.name = 'datetime'
        else:
            print('Did not apply filter.')
        if use_kde_peak:
            d, v = self.kde_peak(data=data)
        else:
            d = data.Wad_Dia.values
            v = data.vel_v.values
        if cut_d:
            dcut = self.d_cut(**cut_kws)
            d = d[d < dcut]
            v = v[d < dcut]
        if use_kde_peak:
            num = np.array([bindata(diam, self.binwidth, data=data,
                                    xname='Wad_Diam',
                                    yname='vel_v').vel_v.count() for diam in d])
            d = d[num > bin_num_min]
            v = v[num > bin_num_min]
            sig = [self.dbin(diam, self.binwidth, data=data, xname='Wad_Diam',
                             yname='vel_v').vel_v.sem() for diam in d]
        else:
            sig = np.ones(d.size)
        vfit.x = d
        vfit.y = v
        vfit.x_unfiltered = origdata.Wad_Dia.values
        vfit.y_unfiltered = origdata.vel_v.values
        if use_curve_fit:
            params, pcov = vfit.find_fit(**kwargs)
            perr = vfit.perr()   # standard errors of d, v
            if try_flip and not use_kde_peak:
                fiti = fitclass(flipped=True)
                datai, stdarri, xlimsi, ymini, ymaxi = self.filter_outlier(data=data,
                                                               frac=frac,
                                                               flip=True)
                fiti.x = datai.Wad_Dia.values
                fiti.y = datai.vel_v.values
                fiti.x_unfiltered = origdata.Wad_Dia.values
                fiti.y_unfiltered = origdata.vel_v.values
                paramsi, pcovi = fiti.find_fit(**kwargs)
                perri = fiti.perr()
                if plot_flip:
                    f, axarr = plt.subplots(1, 3, sharey=True, sharex=True, figsize=(12, 6))
                    for ax in axarr:
                        fiti.plot(ax=ax, label='flipped %.4f' % perri[1])
        else:
            result = minimize(vfit.cost, vfit.quess, method='Nelder-Mead',
                              args=(d, v, sig))
            vfit.params = result.x
        if plot_flip and use_curve_fit:
            self.plot_flip(axarr, f, vfit, data, datai, origdata, perr)
        desfmt = '{0:.4f}'
        fitstr = 'standard'
        errstr = ''
        fitout = vfit
        if use_curve_fit:
            errstr += 'err: std ' + desfmt.format(perr[1])
            if try_flip:
                errstr += ', inv ' + desfmt.format(perri[1])
                if (perr[1] > perri[1] and fiti.is_good()) or force_flip:
                    fitout = fiti
                    fitstr = 'flipped'
        print(fitstr + ' fit: ' + str(fitout) + '; ' + str(partcount) + ' particles; ' + errstr)
        return fitout, std, hwfm    # TODO: wrong std, hwfm when flipped

    def plot_flip(self, axarr, f, vfit, data, datai, origdata, perr):
        filterstr = ['D-binned filter', 'v-binned filter', 'unfiltered']
        for ax in axarr:
            vfit.plot(ax=ax, label='original %.4f' % perr[1])
        self.plot(ax=axarr[0], data=data, ymax=3)
        self.plot(ax=axarr[1], data=datai, ymax=3)
        self.plot(ax=axarr[2], data=origdata, ymax=3)
        for i, ax in enumerate(axarr):
            ax.set_title(ax.get_title() + ', ' + filterstr[i])
        plt.legend()
        f.tight_layout()
        fname = data.index[-1].strftime('%H%M.eps')
        datedir = data.index[-1].strftime('%Y%m%d')
        f.savefig(os.path.join(ensure_dir(os.path.join(RESULTS_DIR, 'pip2015',
                                                       'fitcomparison',
                                                       datedir)), fname))
        return axarr

    def find_fits(self, rule, varinterval=True, empty_on_fail=True, **kwargs):
        print('Calculating velocity fits for given sampling frequency...')
        names = []
        fits = []
        stds = []
        hwfms = []
        for name, group in self.grouped(rule=rule, varinterval=varinterval):
            try:
                newfit, std, hwfm = self.find_fit(data=group, name=name,
                                                  try_flip=self.flip,
                                                  **kwargs)
            except RuntimeError as err:
                print('%s: %s' % (name, err))
                print('Particle count: %s' % group.vel_v.count())
                if len(fits) == 0 or empty_on_fail:
                    print('Using an empty fit')
                    newfit = self.default_fit()
                    std = np.nan
                    hwfm = np.nan
                else:
                    print('Using fit from previous time step.')
                    newfit = fits[-1]
                    std = stds[-1]
                    hwfm = hwfms[-1]
            fits.append(newfit)
            names.append(name)
            stds.append(std)
            hwfms.append(hwfm)
        self.std = pd.concat(stds)
        self.hwfm = pd.concat(hwfms)
        if varinterval:
            timestamps = names
        else:
            timestamps = pd.DatetimeIndex(names, freq=rule)
        if self.fits.empty:
            self.fits = pd.DataFrame(fits, index=timestamps,
                                     columns=[newfit.name])
        elif self.fits.index.equals(timestamps):
            self.fits[newfit.name] = fits
        else:
            self.fits = pd.DataFrame(fits, index=timestamps,
                                     columns=[newfit.name])
        return self.fits

    def fit_params(self, fit_type=None):
        """Return DataFrame of fit parameters."""
        if fit_type is None:
            fit_type = self.default_fit.name
        letter = 'abcdef'
        params = self.fits[fit_type].apply(lambda vfit: vfit.params)
        paramlist = []
        for i in range(params.values[0].size):
            param = params.apply(lambda p: p[i])
            param.name = letter[i]
            paramlist.append(param)
        return merge_multiseries(*paramlist)

    def partcount(self, rule, varinterval):
        return self.grouped(rule=rule, varinterval=varinterval).Part_ID.count()

    def kde(self, data=None):
        """kernel-density estimate of d,v data using gaussian kernels"""
        if data is None:
            data = self.good_data()
        d = data.Wad_Dia.values
        v = data.vel_v.values
        return kde(d, v)

    def grids(self, data=None):
        if data is None:
            data = self.good_data()
        d = data.Wad_Dia.values
        v = data.vel_v.values
        dmax = d.max()+20*self.binwidth
        dbins = self.dbins[self.dbins < dmax]
        num_vbins = round(len(self.dbins)/5)
        return np.meshgrid(dbins, np.linspace(v.min(), v.max(), num_vbins))

    def kde_grid(self, data=None):
        """Calculate kernel-density estimate with given resolution."""
        X, Y = self.grids(data)
        points = np.vstack([X.ravel(), Y.ravel()])
        kernel = self.kde(data)
        Z = np.reshape(kernel(points).T, X.shape)
        return X, Y, Z

    def kde_peak(self, **kwargs):
        """the most propable velocities for each diameter in grid"""
        D, V, Z = self.kde_grid(**kwargs)
        x = D[0, :]
        y = V[:, 0][Z.argmax(axis=0)]
        return x, y

    def plot_kde(self, ax=None, **kwargs):
        """Plot kde grid."""
        if ax is None:
            ax = plt.gca()
        D, V, Z = self.kde_grid(**kwargs)
        pc = ax.pcolor(D, V, Z, cmap=plt.cm.gist_earth_r)
        return pc

    def plot_fits(self, fit_type=None, savefigs=False, path='.',
                  fname='%Y%m%d_%H%M.eps', **kwargs):
        if fit_type is None:
            fit_type = self.default_fit.name    # as a string
        fits = self.fits[fit_type]
        axlist = []
        flist = []
        for i, vfit in fits.iteritems():
            flist.append(plt.figure())
            axlist.append(vfit.plot(**kwargs))
            plt.savefig(os.path.join(path, i.strftime(fname)))
        return pd.DataFrame(index=fits.index, data={'fit':fits, 'fig':flist,
                                                    'ax':axlist})

    def plot_fit(self, tstep=None, **kwargs):
        if tstep is None:
            fits = [self.find_fit()[0]]
        else:
            fits = self.fits.loc[tstep].values
        for vfit in fits:
            vfit.plot(**kwargs)

    def plot(self, data=None, hexbin=True, ax=None, xmax=None, ymax=None,
             show_particle_count=False, colormap='gray_r', ygrid=True,
             hexsize=12, **kwargs):
        """Plot velocity data."""
        if ax is None:
            ax = plt.gca()
        if data is None:
            data = self.good_data()
        margin = 0.1
        if xmax is None:
            xmax = np.ceil(self.d_cut(frac=0.05))
        right = xmax-2+margin
        partcount = data.Part_ID.count()
        if partcount < 1:
            return ax
        if hexbin:
            data.plot(x='Wad_Dia', y='vel_v', kind='hexbin', label='hexbinned',
                      ax=ax, gridsize=int(hexsize*data.Wad_Dia.max()**0.5),
                      colormap=colormap, **kwargs)
        else:
            data.plot(x='Wad_Dia', y='vel_v', style=',', ax=ax,
                      alpha=0.2, color='black', label='pip raw', **kwargs)
        #fit.gunn_kinzer.plot(dmax=20, label='Gunn&Kinzer', ax=ax, zorder=5, ls='--')
        if ymax is None:
            ymax = data.vel_v.max() + margin
        ax.axis([0, xmax, 0, ymax])
        ax.yaxis.grid(ygrid)
        t_start = data.index[0]-datetime.timedelta(minutes=1)
        t_end = data.index[-1]
        label_format = '%H:%M'
        ax.set_title('%s-%s UTC' % (t_start.strftime(label_format),
                                    t_end.strftime(label_format)))
        if show_particle_count:
            ax.text(right, margin, 'particle count: %s' % str(partcount))
        ax.set_ylabel('Fall velocity (m/s)')
        ax.set_xlabel('D (mm)')
        ax.legend(loc='upper right')
        return ax

    def plots(self, rule=None, separate=True, peak=False, save=False, ncols=1,
              prefix='', suffix='.png', ymax=None, plotfit=True, savedir=None,
              **kwargs):
        """Plot datapoints and fit for each timestep."""
        ngroups = self.grouped(rule=rule).ngroups
        #nrows = int(np.ceil(ngroups/ncols))
        if save:
            home = os.curdir
            if 'HOME' in os.environ:
                home = os.environ['HOME']
            if savedir is None:
                savedir = os.path.join(home, 'Pictures', 'vel_plots')
            #suffix = '_' + self.rule
        if not separate:
            f, axarr = plt.subplots(1, ngroups, sharex='col', sharey='row',
                                    figsize=(ngroups*8, 7), tight_layout=True)
        else:
            axarr = []
            farr = []
        if ymax is None:
            self.good_data().vel_v.max()
        for i, (name, group) in enumerate(self.grouped(rule=rule)):
            if separate:
                f = plt.figure(dpi=175, figsize=(3.5, 3))
                ax = plt.gca()
                farr.append(f)
                axarr.append(ax)
            if group.Part_ID.count() < 1:
                continue
            if plotfit:
                self.plot_fit(tstep=name, zorder=6, ax=axarr[i], marker=',',
                              alpha=0.3)
            self.plot(data=group, ax=axarr[i],
                      ymax=ymax, **kwargs)
            f.tight_layout()
            if save and separate:
                t = group.index[-1]
                fname = t.strftime(prefix + '%Y%m%d%H%M' + suffix)
                f.savefig(os.path.join(savedir, fname))
            if peak:
                axarr[i].scatter(*self.kde_peak(data=group), label='kde peak')
        if save and not separate:
            fname = t.strftime('%Y%m%d' + suffix)
            f.savefig(os.path.join(savedir, fname))
        return axarr


class PipPart(InstrumentData):
    """PIP particle tables"""
    def __init__(self, filenames, dt_start=None, dt_end=None, **kwargs):
        print('Reading PIP particle data...')
        InstrumentData.__init__(self, filenames, **kwargs)
        self.name = 'pip_part'
        dtype = {'Year': np.int32, 'Month': np.int32, 'Day': np.int32,
                 'Hr': np.int32, 'Min': np.int32, 'Sec': np.int32}
        if self.data.empty:
            for filename in filenames:
                newdata = pd.read_csv(filename, delim_whitespace=True,
                                      skiprows=[0, 1, 2, 3, 4, 5, 6, 7, 9],
                                      index_col='datetime',
                                      parse_dates={'datetime': ['Year',
                                                                'Month',
                                                                'Day', 'Hr',
                                                                'Min', 'Sec']},
                                      date_parser=datetime.datetime,
                                      dtype=dtype)
                self.data = self.data.append(newdata)
            self.finish_init(dt_start, dt_end)


RHO_0 = 1.225 # sea level air density, from http://www.aviation.ch/tools-atmosphere.asp
#RHO_H = 1.218 # density at 56 m, same source
RHO_W = 1000.0

## Disdrometer sampling area, from 
## "RAINDROP SIZE DISTRIBUTION OVER NORTHEASTERN COAST OF BRAZIL"
#DSD_A = 0.0050 

TAU = 2*np.pi

def ar(d):
    return 1/aux.dsr_thurai_2007(d)

def drop_zeros_rows(df):
    return df[(df.T != 0).any()]

class PSD(InstrumentData):
    """General PSD"""
    # TODO: fix overlap with PipDSD
    # TODO: make PipPSD inherit this
    def __init__(self, data=None, binwidth=None, bin_v=None):
        if data is not None:
            self.data = data.resample('1min').fillna(0)
        self.binwidth = binwidth
        self.bin_v = pd.Series(bin_v, index=data.columns.values)
        self.drop_empty = True

    def plot(self, img=True, **kwargs):
        """Plot particle size distribution over time."""
        if img:
            plt.matshow(self.good_data(**kwargs).transpose(), norm=LogNorm(),
                        origin='lower', vmin=0.00001)
        else:
            plt.pcolor(self.good_data(**kwargs).transpose(), norm=LogNorm(),
                       vmin=0.00001)
        plt.colorbar()
        plt.title('DSD')
        plt.xlabel('time (UTC) BROKEN')
        plt.ylabel('D (mm) BROKEN')

    def binwidth_df(self):
        return pd.DataFrame(data=self.binwidth, index=self.bin_cen()).T

    def good_data(self, drop_empty=True):
        data = self.data
        if self.drop_empty and drop_empty:
            data = drop_zeros_rows(data)
        return data

    def intensity(self):
        return self.sum_over_d(self.rr)

    def v(self, d):
        return self.bin_v[self.bin_select(d)]

    def n(self, d):
        return self.good_data()[d]

    def bin_cen(self):
        return self.good_data().columns.values

    def bin_lower(self):
        return self.bin_cen() - 0.5*self.binwidth

    def bin_upper(self):
        return self.bin_lower() + self.binwidth

    def bin_edge(self):
        return np.append(self.bin_lower(), self.bin_upper()[-1])

    def bin_edge_df(self):
        edge = pd.DataFrame(data=[self.bin_lower(), self.bin_upper()]).T
        edge.columns = ['lower', 'upper']
        edge.index = self.bin_cen()
        return edge

    def bin_select(self, d):
        for i, edge in self.bin_edge_df().iterrows():
            if d > edge.lower and d < edge.upper:
                return edge.name
        return

    def series_zeros(self):
        s = self.good_data()[self.good_data().columns[0]]*0
        s.name = 'series'
        return s

    def sum_over_d(self, func, **kwargs):
        dD = self.binwidth
        result = self.series_zeros()
        for d in self.bin_cen():
            result = result.add(func(d, **kwargs)*dD[d], fill_value=0)
        return result

    def rr(self, d):
        return 3.6e-3*TAU/12*(ar(d)*d)**2*1/ar(d)*d*self.v(d)*self.n(d)

    def to_tm(self, data=None):
        if data is None:
            data = self.good_data().mean()
        return psd.BinnedPSD(bin_edges=self.bin_edge(),
                             bin_psd=data.values)

    def to_tm_series(self, resample=None):
        if resample is None:
            data = self.good_data()
        else:
            data = self.good_data(drop_empty=False).resample(resample, how=np.mean,
                                                             closed='right',
                                                             label='right')
        return data.apply(self.to_tm, axis=1)
