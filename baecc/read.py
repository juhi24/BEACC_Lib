# -*- coding: utf-8 -*-
"""
tools for reading and working with baecc data
@author: Jussi Tiira
"""
import os
import datetime
import pandas as pd
import numpy as np
import copy
import pickle
import warnings
import netCDF4 as nc
import hashlib
import gc
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import ticker
from glob import glob
from j24 import ensure_dir, home
from baecc.instruments import Pluvio, PipPSD, PipV, Radar

# general configuration
DEBUG = False
CGS_UNITS = True # display units in cgs instead of SI

# CONFIG default paths
HOME = home()
DATA_DIR = os.path.join(HOME, 'DATA')
USER_DIR = os.path.join(HOME, '.baecc')
RESULTS_DIR = os.path.join(HOME, 'results')
CACHE_DIR = os.path.join(HOME, '.cache', 'baecc')
H5_FILE = 'baecc.h5'
H5_PATH = os.path.join(DATA_DIR, H5_FILE)
PIPV_SUBPATH = 'PIP/a_Velocity_Tables/004%s/*2.dat'
DSD_SUBPATH = 'PIP/a_DSD_Tables/004%s*.dat'
P200_SUBPATH = 'Pluvio200/pluvio200_??_%s*.txt'
P400_SUBPATH = 'Pluvio400/pluvio400_??_%s*.txt'
RADAR_SUBPATH = 'Radar/%s/tmp%s*M1.a1.%s.*'
MSGTLD = '.msg'
PICKLETLD = '.pkl'

# PIP-observed particle size to the volume equivalent diameter
PHI = 0.92

# constants
if CGS_UNITS:
    RHO_SCALE = 1e-3
    RHO_UNITS = 'g$\,$cm$^{-3}$'
else:
    RHO_SCALE = 1
    RHO_UNITS = 'kg$\,$m$^{-3}$'


def rho_scale(axis):
    axis_scale(axis, RHO_SCALE)


def axis_scale(axis, scale):
    formatter = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*scale))
    axis.set_major_formatter(formatter)


def combine2str(*identifiers):
    return ''.join(tuple(map(str, identifiers)))


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


def datenum2datetime(matlab_datenum):
    """Convert MATLAB datenum to datetime."""
    return datetime.datetime.fromordinal(int(matlab_datenum)) + datetime.timedelta(days=matlab_datenum % 1) - datetime.timedelta(days=366)


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
    dsd = PipPSD(dsd_files)
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
    """Read ASCII data and export to hdf5."""
    if instrdict is None:
        instrdict = {PIPV_SUBPATH: PipV,
                     DSD_SUBPATH: PipPSD,
                     P200_SUBPATH: Pluvio,
                     P400_SUBPATH: Pluvio}
    hdf_file = os.path.join(datadir, outname)
    for key in instrdict:
        instr = instrdict[key].from_raw(dtstrlist, subpath=key, datadir=datadir)
        instr.to_hdf(filename=hdf_file)
        del(instr)
        gc.collect()


def plot_psd(data, ax=None):
    """Plot particle size distribution over time."""
    if ax is None:
        ax = plt.gca()
    qmesh = ax.pcolormesh(data.index.values, data.columns.values, data.transpose(),
                  norm=LogNorm())
    plt.colorbar(qmesh, ax=ax)
    ax.set_title('PIP PSD')
    ax.set_xlabel('time (UTC)')
    ax.set_ylabel('D (mm)')
    return ax


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
        if os.path.isfile(msgpath):
            data = pd.read_msgpack(msgpath)
        else:
            ensure_dir(cd)
            data = func(**kwargs)
            data.to_msgpack(msgpath)
        return data

    def pkl_io(self, name, func, **kwargs):
        cd = self.cache_dir()
        pklpath = os.path.join(cd, name + PICKLETLD)
        if os.path.isfile(pklpath):
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
    def __init__(self, filenames=None, data=None, hdf_table=None, use_cache=True):
        """Read from either ASCII data file or hdf5."""
        self.filenames = filenames
        if data is None:
            self.data = pd.DataFrame()
        else:
            self.data = data
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


class PipPart(InstrumentData):
    """PIP particle tables"""
    def __init__(self, filenames=None, dt_start=None, dt_end=None, **kwargs):
        InstrumentData.__init__(self, filenames, **kwargs)
        self.name = 'pip_part'
        dtype = {'Year': np.int32, 'Month': np.int32, 'Day': np.int32,
                 'Hr': np.int32, 'Min': np.int32, 'Sec': np.int32}
        if self.data.empty:
            print('Reading PIP particle data...')
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

    def plot(self, **kwargs):
        """wrapper for plot_psd"""
        return plot_psd(self.good_data(), **kwargs)

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
