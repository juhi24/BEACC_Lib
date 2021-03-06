# coding: utf-8
from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import time
import os
import numpy as np
import pandas as pd
from warnings import warn
from datetime import datetime, timedelta
import baecc
from baecc import instruments, caching


P200_SUBPATH = 'Pluvio200/pluvio200_??_%s*.txt'
P400_SUBPATH = 'Pluvio400/pluvio400_??_%s*.txt'


def parse_datetime(datestr, include_sec=False):
    datestr = str(int(datestr))
    t = time.strptime(datestr, '%Y%m%d%H%M%S')
    if include_sec:
        t_end = 6
    else:
        t_end = 5
    #return datetime.datetime(*t[:t_end], tzinfo=datetime.timezone.utc)
    return datetime(*t[:t_end])


class Pluvio(instruments.InstrumentData, instruments.PrecipMeasurer):
    """Pluviometer data handling"""
    def __init__(self, filenames=[], dt_start=None, dt_end=None, name=None,
                 **kwargs):
        """Create a Pluvio object using data from a list of files."""
        instruments.InstrumentData.__init__(self, filenames, **kwargs)
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
        self.maxdelta = timedelta(hours=1)
        self.name = name.lower()
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
        if self.data.empty:
            if len(filenames)<1:
                warn('No input files or data given.')
                return
            print('Reading pluviometer data...')
            if self.name is None:
                self.name = os.path.basename(os.path.dirname(self.filenames[0])).lower()
            col_abbr = ['datestr',
                        'i_rt',
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
                    newdata = pd.read_csv(filename, sep=';',
                                          names=col_abbr,
                                          parse_dates={'datetime':['datestr']},
                                          date_parser=parse_datetime,
                                          index_col='datetime',
                                          skip_blank_lines=True,
                                          error_bad_lines=False,
                                          warn_bad_lines=True,
                                          verbose=baecc.DEBUG)
                    self.data = self.data.append(newdata)
                except NotImplementedError as err:
                    print('\n%s: %s' % (filename, format(err)))
            print()
            #self.data.drop(['i_rt'], 1, inplace=True) # crap format
        self.buffer = timedelta(0)
        self.finish_init(dt_start, dt_end)
        self.data['group'] = self.data.acc_nrt.astype(bool).astype(int).cumsum().shift(1).fillna(0)
        self.data['heating'] = self.data.heating.astype(float)  # sometimes this are interpreted as int
        self.data['status'] = self.data.status.astype(float)   # sometimes this are interpreted as int

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
    def from_raw(cls, subpath=P200_SUBPATH, *args, **kwargs):
        return super().from_raw(*args, subpath=subpath, **kwargs)

    def fingerprint(self):
        identifiers = [super().fingerprint(), self.name, self.shift_periods,
                       self.shift_freq, self.varinterval]
        if self.varinterval:
            identifiers.extend([self.n_combined_intervals])
        idstr = caching.combine2str(*identifiers)
        return caching.fingerprint(idstr)

    def good_data(self):
        if self.stored_good_data is not None:
            return self.stored_good_data
        data = self.data.copy()
        data = data[data.status==0].copy()
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
            dt = pd.to_datetime(dt)
        self.buffer = timedelta(hours=2)
        if dt_start is None or dt_end is None:
            self.buffer = timedelta(0)
        elif dt_start-self.buffer < self.data.index[0] or dt_end+self.buffer > self.data.index[-1]:
            self.buffer = timedelta(0)
        self.data = self.data[dt_start-self.buffer:dt_end+self.buffer]

    def timeshift(self):
        """timeshift as timedelta"""
        if self.shift_periods == 0:
            return timedelta(0)
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

    def amount(self, crop=True, shift=True, upsample_noprecip=True, **bucketkwargs):
        if not self.varinterval:
            return self.constinterval_amount(shift=shift, **bucketkwargs)
        am0 = self.good_data()[self.amount_col]
        if shift and self.shift_periods>0:
            am0 = am0.tshift(periods=self.shift_periods, freq=self.shift_freq)
        n = self.n_combined_intervals
        am = am0[am0 > 0].rolling(center=False, window=n).sum().iloc[n-1::n]
        am = pd.concat([am0[:1], am, am0[-1:]])
        am = am.groupby(am.index).first() # drop index duplicates
        if crop:
            am = am[self.dt_start():self.dt_end()]
        if upsample_noprecip:
            dt = self.tdelta(limit_maxdelta=False, upsample_noprecip=False)
            for t, ddt in dt[dt>self.maxdelta].iteritems():
                h = int(ddt.seconds/3600)
                for hh in range(h):
                    tt = t-timedelta(hours=hh+1)
                    am[tt] = 0
            am.sort_index(inplace=True)
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

    def tdelta(self, limit_maxdelta=True, **kws):
        """lengths of timesteps as Series of timedeltas"""
        a = self.amount(crop=False, **kws)
        delta = pd.Series(a.index.to_pydatetime(), index=a.index).diff()
        if limit_maxdelta:
            delta[delta > self.maxdelta] = self.maxdelta
        delta = pd.to_timedelta(delta)
        delta.name = 'tdelta'
        out = delta[self.dt_start():self.dt_end()].fillna(self.maxdelta)
        out.iloc[0] = timedelta(0) # First delta is unknown
        return out

    def start_time(self):
        """timestep starting timestamps"""
        tdelta = self.tdelta()
        td = pd.DataFrame(tdelta)
        td['end'] = td.index
        start = td.end.sub(tdelta)
        start.name = 'start'
        return start

    def half_time(self):
        """timestep middle timestamps"""
        t_half = self.start_time()+self.tdelta().div(2)
        t_half.name = 'middle'
        return t_half

    def grouper(self, shift=True):
        """data group names (timestamp) for each data timestamp"""
        ticks = self.good_data()[self.amount_col].astype(bool)
        if shift:
            ticks = ticks.tshift(periods=self.shift_periods,
                                 freq=self.shift_freq)
        ticktime = self.amount(shift=shift).index
        dtgroups = pd.Series(ticktime, index=ticktime).reindex(ticks.index).bfill()[self.dt_start():self.dt_end()]
        dtgroups.name = 'group'
        last_index = self.tdelta().index[-1]
        return pd.DataFrame(dtgroups[dtgroups.notnull()])[:last_index]

    def groupby_interval(self, data):
        """Group data by integration time intervals."""
        return baecc.merge_series(data, self.grouper()).groupby('group')

