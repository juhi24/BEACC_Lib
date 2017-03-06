# coding: utf-8
import pandas as pd
import baecc
from baecc import caching, case
from datetime import datetime, timedelta

class EventsCollection(caching.Cacher):
    """Manage a table of precipitation events."""
    def __init__(self, csv, dtformat='%d %B %H UTC', default_col='paper',
                 **cacher_kws):
        """Read event metadata from a csv file."""
        self.dtformat = dtformat
        self.default_col = default_col
        self.events = pd.read_csv(csv, parse_dates=['start', 'end'],
                                  date_parser=self.parse_datetime)
        self.events.sort_values(by=['start', 'end'], inplace=True)
        self.events.start += timedelta(seconds=1)
        caching.Cacher.__init__(self, **cacher_kws)

    def parse_datetime(self, dtstr):
        #date = datetime.strptime(dtstr+'+0000', self.dtformat+'%z')
        date = datetime.strptime(dtstr, self.dtformat)
        return date

    def fingerprint(self):
        identifiers = [self.events]
        for c in self.events[self.default_col]:
            identifiers.extend(c.fingerprint())
        idstr = caching.combine2str(*identifiers)
        return caching.fingerprint(idstr)

    def total_duration(self, events_col=None, winter=None):
        if events_col is None:
            events_col = self.default_col
        t_tot = pd.timedelta_range(0,0)[0]
        cases = self.events[events_col]
        if winter is not None:
            cases = cases.loc[winter]
        for c in cases:
            t_tot += c.duration()
        return t_tot

    def duration_weights(self, events_col=None):
        if events_col is None:
            events_col = self.default_col
        weights = []
        t_tot = self.total_duration()
        for c in self.events[events_col]:
            weights.append(c.duration()/t_tot)
        return pd.Series(weights, index=self.events.index)


    def vel_data(self, events_col=None, winter=None):
        if events_col is None:
            events_col = self.default_col
        cases = self.events[events_col]
        if winter is not None:
            cases = cases.loc[winter]
        vdatalist = []
        for c in cases:
            data = c.instr['pipv'].good_data()
            vdatalist.append(data)
        vdata = pd.concat(vdatalist)
        return vdata


    def add_data(self, data, autoshift=True, autobias=True):
        """Add data from a Case object."""
        cases = []
        for (i, e) in self.events.iterrows():
            cases.append(data.between_datetime(e.start, e.end,
                                               autoshift=autoshift,
                                               autobias=autobias))
        self.events[data.instr['pluvio'].name] = cases

    def autoimport_data(self, datafile=baecc.H5_PATH, autoshift=False,
                        autobias=False, radar=False, **casekwargs):
        """Import data from a hdf file."""
        timemargin = timedelta(hours=3)
        dt_start = self.events.iloc[0].start - timemargin
        dt_end = self.events.iloc[-1].end + timemargin
        for pluvio_name in ('pluvio200', 'pluvio400'):
            data = case.Case.from_hdf(dt_start, dt_end, autoshift=False,
                                 filenames=[datafile], radar=radar,
                                 pluvio_name=pluvio_name, **casekwargs)
            if data is not None:
                self.add_data(data, autoshift=autoshift, autobias=autobias)

    def summary(self, col=None, dtformat='%Y %b %d', concatkws={},
                **kwargs):
        if col is None:
            col = self.default_col
        sumlist = []
        for c in self.events[col]:
            sumlist.append(c.summary(dtformat=dtformat, **kwargs))
        return pd.concat(sumlist, **concatkws)

    def pluv_grouper(self, events_col=None, winter=None):
        if events_col is None:
            events_col = self.default_col
        cases = self.events[events_col]
        if winter is not None:
            cases = cases.loc[winter]
        grouperlist = []
        for c in self.events[events_col]:
            grouperlist.append(c.instr['pluvio'].grouper())
        return pd.concat(grouperlist)

    def split_index(self, date=pd.datetime(2014,7,1),
                    names=('first', 'second')):
        isfirst = self.events.start < date
        idf = isfirst.copy()
        idf[isfirst] = names[0]
        idf[-isfirst] = names[1]
        tuples = list(zip(*(idf.values, idf.index.values)))
        index = pd.MultiIndex.from_tuples(tuples, names=('winter', 'i'))
        self.events.index = index
