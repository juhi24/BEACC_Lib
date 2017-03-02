# coding: utf-8
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import linecache
import copy
import datetime
from os import path
import baecc
from baecc import instruments, caching


def file_shorter_than(fname, limit):
    with open(fname) as f:
        for i, l in enumerate(f):
            if i+2>limit:
                return False
    return True


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


class PipPSD(instruments.InstrumentData):
    """PIP particle size distribution data handling"""
    def __init__(self, filenames=None, dt_start=None, dt_end=None, **kwargs):
        """Create a PipDSD object using data from a list of PIP DSD table
        files."""
        instruments.InstrumentData.__init__(self, filenames, **kwargs)
        self.name = 'pip_dsd'
        self.use_voleq_d = True
        common_csv_kws = {'skiprows': 8,
                          'header': 3,
                          'parse_dates': {'datetime':['hr_d', 'min_d']},
                          'date_parser': self.parse_datetime,
                          'index_col': 'datetime',
                          'verbose': baecc.DEBUG}
        if self.data.empty:
            print('Reading PIP PSD data...')
            for filename in filenames:
                if baecc.DEBUG:
                    print(filename)
                else:
                    print('.', end='')
                if file_shorter_than(filename, 14):
                    # file has no data
                    continue
                self.current_file = filename
                # File format changed
                if int(path.split(filename)[1][3:11]) > 20141124:
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
    def from_raw(cls, *args, subpath=instruments.PSD_SUBPATH, **kwargs):
        return super().from_raw(*args, subpath=subpath, **kwargs)

    def parse_datetime(self, hh, mm):
        dateline = linecache.getline(self.current_file, 6)
        datearr = [int(x) for x in dateline.split()]
        d = datetime.date(*datearr)
        t = datetime.time(int(hh), int(mm))
        return datetime.datetime.combine(d, t)#.replace(tzinfo=datetime.timezone.utc)

    def bin_cen(self):
        """Return array of bin centers as area equivalent diameter."""
        return self.good_data().columns.values

    def bin_width(self):
        """as a Series of area equivalent diameter"""
        # TODO: should be calculated using actual bin edges
        return pd.Series(self.bin_cen(), index=self.bin_cen()).diff().bfill()

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

    def plot(self, data_kws={}, **kws):
        """wrapper for plot_psd"""
        return plot_psd(self.good_data(**data_kws), **kws)

    def fingerprint(self):
        identifiers = (super().fingerprint(), self.use_voleq_d)
        idstr = caching.combine2str(*identifiers)
        return caching.fingerprint(idstr)

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
        too_large = bin_cen[bin_cen > 25]
        data_fltr = gain_correction*data.drop(too_small, 1).drop(too_large, 1)
        if self.use_voleq_d:
            data_fltr.columns = data_fltr.columns/baecc.PHI
            data_fltr = data_fltr*baecc.PHI
        return data_fltr
