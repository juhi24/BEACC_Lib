# coding: utf-8
import copy
import pandas as pd
from baecc.caching import Cacher, fingerprint


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