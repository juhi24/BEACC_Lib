# coding: utf-8
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
