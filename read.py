"""tools for reading and working with baecc data"""
import time
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import datetime
import pandas as pd
import numpy as np
import linecache
import copy
from scipy import stats
from scipy.optimize import curve_fit, fmin, minimize

GUNN_KINZER = (9.65, 10.30/9.65, 0.6)
CACHE_DIR = 'cache'
MSGTLD = '.msg'

def datenum2datetime(matlab_datenum):
    """Convert MATLAB datenum to datetime."""
    return datetime.datetime.fromordinal(int(matlab_datenum)) + datetime.timedelta(days=matlab_datenum%1) - datetime.timedelta(days=366)

def msg_io(msgpath, func, **kwargs):
    if os.path.exists(msgpath):
        data = pd.read_msgpack(msgpath)
    else:
        data = func(**kwargs)
        data.to_msgpack(msgpath)
    return data

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

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

class Cacher:
    """common methods to use msg cache"""
    def __init__(self, use_cache=True, storefilename='store.h5'):
        self.use_cache = use_cache
        self.storefilename = storefilename

    def msger(self, name, func, *cache_dir_args, **kwargs):
        if self.use_cache:
            msgpath = os.path.join(self.cache_dir(*cache_dir_args), name + MSGTLD)
            data = msg_io(msgpath, func, **kwargs)
        else:
            data = func(**kwargs)
        return data

    def cache_dir(self, dt_start, dt_end, *extra_dir_names):
        dtstrformat = '%Y%m%d%H%M'
        dtstr = dt_start.strftime(dtstrformat) + '-' + dt_end.strftime(dtstrformat)
        cache_dir = os.path.join(CACHE_DIR, dtstr, *extra_dir_names)
        ensure_dir(cache_dir)
        return cache_dir

    def msgpath(self, name, *cache_dir_args):
        """wrapper for full msgpack file path"""
        return os.path.join(self.cache_dir(*cache_dir_args), name + MSGTLD)

    def store_path(self, *cache_dir_args):
        return os.path.join(self.cache_dir(*cache_dir_args), self.storefilename)

class Fit:
    """parent for different fit types"""
    def __init__(self, x=None, y=None, sigma=None, params=None, name='fit'):
        self.params = params
        self.name = name
        self.x = x
        self.y = y
        self.sigma = sigma

    def func(self, x, a=None):
        """Fit function. If no coefficients are given use stored ones."""
        if a is None:
            return self.func(x, *self.params)
        pass

    def penalty(self, params):
        """penalty function used by the cost function"""
        return 0

    def plot(self, xmax=20, samples=1000, ax=None, label=None,
             source_data=False, marker='ro', **kwargs):
        """Plot fit curve and fitted data."""
        if ax is None:
            ax = plt.gca()
        if self.params is None:
            return ax
        x = np.linspace(0, xmax, samples)
        y = [self.func(xi, *self.params) for xi in x]
        if label is None:
            label = r'$' + str(self) + r'$'
        ax.plot(x, y, label=label, linewidth=2)
        if source_data:
            ax.plot(self.x, self.y, marker, **kwargs)
        return ax

    def cost(self, params, xarr, yarr, sigarr):
        """Cost function that can be used to find fit coefs by minimization."""
        cost = 0
        for i, x in enumerate(xarr):
            y = yarr[i]
            sig = sigarr[i]
            cost += 1/sig**2*(y - self.func(x, *params))**2 + self.penalty(params)
        return cost

    def find_fit(self, store_params=True, **kwargs):
        selection = pd.DataFrame([self.x.notnull(), self.y.notnull()]).all()
        x = self.x[selection]
        y = self.y[selection]
        if self.sigma is not None:
            kwargs['sigma'] = self.sigma[selection]
        params, cov = curve_fit(self.func, x, y, **kwargs)
        if store_params:
            self.params = params
        return params, cov

class ExpFit(Fit):
    """exponential fit of form a*(1-b*exp(-c*D))"""
    def __init__(self, params=None, **kwargs):
        super().__init__(params=params, name='expfit', **kwargs)
        self.quess = (1., 1., 1.)

    def __repr__(self):
        if self.params is None:
            paramstr = 'abc'
        else:
            paramstr = ['{0:.3f}'.format(p) for p in self.params]
        s = '%s(1-%s\exp(-%sD))' % (paramstr[0], paramstr[1], paramstr[2])
        return s.replace('--', '+')

    def func(self, x, a=None, b=None, c=None):
        if a is None:
            return self.func(x, *self.params)
        return a*(1-b*np.exp(-c*x))

    def penalty(self, params):
        return 0
        return 1000*(max(0, 0.1-params[1]) + max(0, 0.4-params[2]))

class PolFit(Fit):
    """polynomial fit of form a*D**b"""
    def __init__(self, params=None, **kwargs):
        super().__init__(params=params, name='polfit', **kwargs)
        self.quess = (1., 1.)

    def __repr__(self):
        if self.params is None:
            paramstr = 'ab'
        else:
            paramstr = ['{0:.3f}'.format(p) for p in self.params]
        return '%sD^{%s}' % (paramstr[0], paramstr[1])

    def func(self, x, a=None, b=None):
        if a is None:
            return self.func(x, *self.params)
        return a*x**b

    def penalty(self, params):
        #return 0
        return 1000*max(0, 0.2-params[1])

gunn_kinzer = ExpFit(params=GUNN_KINZER)

class PrecipMeasurer:
    """parent for classes with precipitation measurements
    Either amount or acc (or both) methods should be overridden."""
    def __init__(self):
        pass

    def amount(self, **kwargs):
        """timestep precipitation in mm"""
        return self.acc(**kwargs).diff()

    def acc(self, **kwargs):
        """precipitation accumulation in mm"""
        return self.amount(**kwargs).cumsum()

    def intensity(self, tdelta=None, **kwargs):
        """precipitation intensity in mm/h"""
        r = self.amount(**kwargs)
        if tdelta is None:
            tdelta = r.index.freq.delta
        frac = pd.datetools.timedelta(hours=1)/tdelta
        return frac * r

class InstrumentData(Cacher):
    """Parent for instrument data classes."""
    def __init__(self, filenames, hdf_table=None, use_cache=True, case=None):
        """Read from either ASCII data file or hdf5."""
        self.filenames = filenames
        self.data = pd.DataFrame()
        self.use_cache = use_cache
        self.case = case # to be used only for msg cache!
        # if filtered data needed often, keep in memory
        self.stored_good_data = None # set to None to disable
        if hdf_table is not None:
            self.name = hdf_table
            self.data = self.data.append(pd.read_hdf(filenames[0], hdf_table))

    def finish_init(self, dt_start, dt_end):
        """Sort and name index, cut time span."""
        self.data.sort_index(inplace=True)
        self.data.index.names = ['datetime']
        self.set_span(dt_start, dt_end)
        self.storefilename = self.name + '.h5'

    def append_data(self, appended):
        #print(len(appended.data.index))
        #print(self.data.dtypes)
        #print(appended.data.dtypes)
        if len(appended.data.index):
            self.data = self.data.append(appended.data)

    def store_good_data(self, **kwargs):
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
        for dt in [dt_start, dt_end]:
            dt = pd.datetools.to_datetime(dt)
        self.data = self.data[dt_start:dt_end]

    def grouped(self, varinterval=True, rule=None, col=None):
        if rule is None:
            rule = self.rule
        if varinterval:
            data = self.good_data()
            if col is not None:
                data = pd.DataFrame(data[col])
            grpd_data = pd.merge(data, rule, left_index=True, right_index=True)
            return grpd_data.groupby('group')
        return self.good_data().groupby(pd.Grouper(freq=rule, closed='right', label='right'))

    def cache_dir(self):
        if self.case is None:
            return
        dt_start, dt_end = self.case.dt_start_end()
        return super().cache_dir(dt_start, dt_end, self.case.pluvio.name)

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
        self.col_suffix = 'nrt'
        self.use_bucket = False
        self._varinterval = True
        self.amount_col = 'acc_' + self.col_suffix
        self.bucket_col = 'bucket_' + self.col_suffix
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
                num_lines = sum(1 for line in open(filename))
                self.current_file = filename
                self.data = self.data.append(pd.read_csv(filename, sep=';',
                            names=col_abbr,
                            skiprows=list(range(1, num_lines, 2)),
                            parse_dates={'datetime':['datestr']},
                            date_parser=self.parse_datetime,
                            index_col='datetime'))
            self.data.drop(['i_rt'], 1, inplace=True) # crap format
        self.buffer = pd.datetools.timedelta(0)
        self.finish_init(dt_start, dt_end)
        self.data['group'] = self.data.acc_nrt.astype(bool).astype(int).cumsum().shift(1).fillna(0)

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

    def parse_datetime(self, datestr, include_sec=False):
        datestr = str(int(datestr))
        t = time.strptime(datestr, '%Y%m%d%H%M%S')
        if include_sec:
            t_end = 6
        else:
            t_end = 5
        return datetime.datetime(*t[:t_end])

    def good_data(self):
        if self.stored_good_data is not None:
            return self.stored_good_data
        data = copy.deepcopy(self.data)
        swap_date = pd.datetime(2014, 5, 16, 8, 0, 0)
        swap_date2 = pd.datetime(2014, 8, 31, 8, 0, 0) # TODO put correct switch date
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
        am = am[am > 0]
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
            accum = accum.tshift(periods=self.shift_periods, freq=self.shift_freq)
        if unbias:
            accum -= self.bias
        accum = accum[self.dt_start():self.dt_end()]
        if filter_evap:
            amount = accum.diff()
            evap = amount[amount < 0]
            evap_accum = evap.reindex(accum.index).fillna(0).cumsum()
            #evap_accum.plot()
            #accum.plot()
            accum -= evap_accum # accumulation should never drop
        return accum.resample(rule, how='last', closed='right', label='right')

    def acc_raw(self):
        """accumulation from raw data"""
        return self.good_data()[self.bucket_col]-self.good_data()[self.bucket_col][0]

    def noprecip_bias(self, lwc, inplace=False):
        """Calculate accumulated bias using LWC."""
        self.lwc = lwc
        accum = self.acc(rule='1min', shift=True, unbias=False, filter_evap=False)
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
        a = self.amount(crop=False)
        delta = pd.Series(a.index.to_pydatetime(), index=a.index).diff()
        longest_delta = pd.datetools.timedelta(hours=1)
        delta[delta > longest_delta] = longest_delta
        delta = pd.to_timedelta(delta)
        return delta[self.dt_start():self.dt_end()].fillna(longest_delta)

    def grouper(self, shift=True):
        ticks = self.good_data()[self.amount_col].astype(bool)
        if shift:
            ticks = ticks.tshift(periods=self.shift_periods, freq=self.shift_freq)
        dtgroups = pd.Series(ticks.index[ticks], index=ticks.index[ticks]).reindex(ticks.index).bfill()[self.dt_start():self.dt_end()]
        #numgroups = ticks.astype(int).cumsum().shift(1).fillna(0)[self.dt_start():self.dt_end()]
        dtgroups.name = 'group'
        last_index = self.tdelta().index[-1]
        return pd.DataFrame(dtgroups[dtgroups.notnull()])[:last_index]

class PipDSD(InstrumentData):
    """PIP particle size distribution data handling"""
    def __init__(self, filenames, dt_start=None, dt_end=None, **kwargs):
        """Create a PipDSD object using data from a list of PIP DSD table files."""
        print('Reading PIP DSD data...')
        InstrumentData.__init__(self, filenames, **kwargs)
        self.d_bin = 0.25
        self.name = 'pip_dsd'
        if self.data.empty:
            for filename in filenames:
                self.current_file = filename
                if int(filename[-16:-8]) > 0: # TODO fixme
                    self.data = self.data.append(pd.read_csv(filename,
                                    engine='python', sep='\t', skiprows=8, skip_footer=1, header=3,
                                    parse_dates={'datetime':['hr_d', 'min_d']},
                                    date_parser=self.parse_datetime,
                                    index_col='datetime'))
                else:
                    self.data = self.data.append(pd.read_csv(filename,
                                    delim_whitespace=True, skiprows=8, header=3,
                                    parse_dates={'datetime':['hr_d', 'min_d']},
                                    date_parser=self.parse_datetime,
                                    index_col='datetime'))
            #self.num_d = self.data[['Num_d']]
            # 1st size bin is crap data, last sometimes nans
            self.d_bin = float(linecache.getline(self.current_file, 11).split()[5])
            self.data.drop(['day_time', 'Num_d', 'Bin_cen', self.data.columns.values[3]], 1,
                           inplace=True)
            self.data.columns = pd.Index([float(i) for i in self.data.columns])
            self.data.sort_index(axis=1)
        self.data.drop_duplicates(inplace=True)
        self.data = self.data.resample('1min').fillna(0)
        self.finish_init(dt_start, dt_end)

    def parse_datetime(self, hh, mm):
        dateline = linecache.getline(self.current_file, 6)
        datearr = [int(x) for x in dateline.split()]
        d = datetime.date(*datearr)
        t = datetime.time(int(hh), int(mm))
        return datetime.datetime.combine(d, t)

    def n(self, d, rule=None, varinterval=True):
        if varinterval:
            grp = self.grouped(rule=rule, varinterval=varinterval, col=d)
            n = grp.mean()
            ns = n[n.columns[0]]
            ns.name = 'N'
            ns.index.name = 'datetime'
            return ns
        return self.good_data()[d].resample(rule, how=np.mean, closed='right',
                                            label='right')

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

    def filter_cat_and_dog(self, data=None, window=5):
        """a rolling window filter for isolated data points"""
        if data is None:
            data = self.data
        is_dog = pd.rolling_count(data.mask(data == 0).T, window).T == 1
        is_dog.ix[:, :window] = False # ignore first columns
        is_dog[is_dog == False] = np.nan
        is_dog = is_dog.ffill(axis=1).fillna(False)
        is_dog = is_dog.astype(np.bool)
        filtered = copy.deepcopy(data)
        filtered[is_dog] = 0
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
            data = self.filter_cat_and_dog(data=data, **kwargs)
        return gain_correction*data

class PipV(InstrumentData):
    """PIP particle velocity and diameter data handling"""
    def __init__(self, filenames, dt_start=None, dt_end=None, **kwargs):
        """Create a PipV object using data from a list of PIP velocity table files."""
        print('Reading PIP particle velocity data...')
        InstrumentData.__init__(self, filenames, **kwargs)
        self.name = 'pip_vel'
        self.dmin = 0.375 # smallest diameter where data is good
        self._fits = pd.DataFrame()
        self.dbins = np.linspace(0.375, 25.875, num=206)
        self.std = pd.DataFrame(columns=self.dbins)
        # half width at fraction of maximum
        self.hwfm = pd.DataFrame(columns=self.dbins)
        self._default_fit = PolFit()
        if self.data.empty:
            for filename in filenames:
                self.current_file = filename
                if int(filename[-23:-15]) > 0: # TODO fixme
                    newdata = pd.read_csv(filename,
                                          engine='python', sep='\t',
                                          skipinitialspace=True, skiprows=8, skip_footer=1,
                                          parse_dates={'datetime':['minute_p']},
                                          date_parser=self.parse_datetime)
                else:
                    newdata = pd.read_csv(filename,
                                          delim_whitespace=True, skiprows=8,
                                          parse_dates={'datetime':['minute_p']},
                                          date_parser=self.parse_datetime)
                newdata.rename_axis(
                    {'vel_v_1':'vel_v', 'vel_h_1':'vel_h',
                     'vel_v_2':'vel_v', 'vel_h_2':'vel_h'}, axis=1, inplace=True)
                if len(newdata.index):
                    self.data = self.data.append(newdata)
                #print(filename)
            #print(len(self.data.index))
            if len(self.data.index):
                self.data.set_index(['datetime', 'Part_ID', 'RecNum'], inplace=True)
                self.data = self.data.groupby(level=['datetime', 'Part_ID']).mean()
                self.data = self.data[self.data.vel_v.notnull()]
            self.data.reset_index(level=1, inplace=True)
        self.finish_init(dt_start, dt_end)

    @property
    def rule(self):
        if self.fits.empty:
            return None
        return self.fits.index.freq.freqstr

    @property
    def binwidth(self):
        d = self.dbins
        return (d[-1]-d[0])/(len(d)-1)

    @property
    def fits(self):
        if self.use_cache:
            try:
                with pd.HDFStore(self.store_path()) as store:
                    return store.get('fits')
            except KeyError:
                return pd.DataFrame()
        return self._fits

    @fits.setter
    def fits(self, fits):
        if self.use_cache:
            with pd.HDFStore(self.store_path()) as store:
                store['fits'] = fits
        else:
            self._fits = fits

    @property
    def default_fit(self):
        return copy.deepcopy(self._default_fit)

    @default_fit.setter
    def default_fit(self, emptyfit):
        self._default_fit = emptyfit

    def grids(self, data=None):
        if data is None:
            data = self.good_data()
        d = data.Wad_Dia.values
        v = data.vel_v.values
        dmax = d.max()+20*self.binwidth
        dbins = self.dbins[self.dbins < dmax]
        num_vbins = round(len(self.dbins)/5)
        return np.meshgrid(dbins, np.linspace(v.min(), v.max(), num_vbins))

    def v(self, d, emptyfit=None, varinterval=True, rule=None):
        if emptyfit is None:
            emptyfit = self.default_fit
        if rule is None:
            rule = self.rule
        if self.fits.empty:
            self.find_fits(rule, fit=emptyfit, varinterval=varinterval)
        elif not varinterval:
            if pd.datetools.to_offset(rule) != self.fits.index.freq:
                # different sampling freq
                self.find_fits(rule, fit=emptyfit, varinterval=varinterval)
        v = []
        for fit in self.fits[emptyfit.name].values:
            v.append(fit.func(d))
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
        return datetime.datetime(yr, mo, dd, hh, int(mm))

    def good_data(self):
        if self.stored_good_data is not None:
            return self.stored_good_data
        return self.data[self.data.Wad_Dia > self.dmin]

    def dbin(self, d, binwidth=None, data=None, vmin=None, vmax=None):
        """Return data that falls into given size bin."""
        if binwidth is None:
            binwidth = self.binwidth
        if data is None:
            data = self.good_data()
        dmin = d-0.5*binwidth
        dmax = d+0.5*binwidth
        vcond = 'Wad_Dia > %s and Wad_Dia < %s' % (dmin, dmax)
        if vmin is not None and vmax is not None:
            vcond += ' and vel_v > %s and vel_v < %s' % (vmin, vmax)
        return data.query(vcond)

    def filter_outlier(self, data=None, frac=0.5):
        filtered = pd.DataFrame()
        HWfracM = []
        std = []
        X, Y, Z = self.kde_grid(data)
        y = Y[:, 0]
        for i in range(0, Z.shape[1]):
            z = Z[:, i]
            z_lim = z.max()*frac
            y_fltr = y[z > z_lim] # FWHM when frac=0.5
            if y_fltr.size == 0:
                continue
            vmin = y_fltr[0]
            vmax = y_fltr[-1]
            HWfracM.append(0.5*(vmax-vmin))
            d = X[0, i]
            bin_data = self.dbin(d=d, data=data, vmin=vmin, vmax=vmax)
            std.append(bin_data.vel_v.std())
            filtered = filtered.append(bin_data)
        # return filtered data, stds and half width at frac*kde_max
        return filtered, np.array(std), np.array(HWfracM)

    def frac_larger(self, d):
        """Return fraction of particles that are larger than d."""
        vdata = self.good_data()
        return vdata[vdata.Wad_Dia > d].vel_v.count()/vdata[vdata.Wad_Dia < d].vel_v.count()

    def d_cut(self, frac=0.05, d0=2):
        """Return d for which given fraction of particles are larger."""
        dcost = lambda d: abs(self.frac_larger(d[0])-frac)
        return fmin(dcost, d0)[0]

    def find_fit(self, fit=None, data=None, kde=False, cut_d=False, frac=0.5,
                 use_curve_fit=True, bin_num_min=5, filter_outliers=True,
                 name=None, **kwargs):
        """Find and store a fit for either raw data or kde."""
        if data is None:
            data = self.good_data()
        if fit is None:
            fit = self.default_fit
        print(data.count()[0])
        if data.count()[0] < 5 and (use_curve_fit or kde):
            print('Too few particles.')
            kde = False
            use_curve_fit = False
        elif filter_outliers:
            data, std, HWfracM = self.filter_outlier(data=data, frac=frac)
            if name is not None:
                std.resize(self.dbins.size, refcheck=False)
                HWfracM.resize(self.dbins.size, refcheck=False)
                self.std = self.std.append(pd.DataFrame(pd.Series(std, index=self.dbins, name=name)).T)
                self.hwfm = self.hwfm.append(pd.DataFrame(pd.Series(HWfracM, index=self.dbins, name=name)).T)
        else:
            print('Could not apply filter.')
        if kde:
            d, v = self.kde_peak(data=data)
        else:
            d = data.Wad_Dia.values
            v = data.vel_v.values
        if cut_d:
            dcut = self.d_cut(**kwargs)
            d = d[d < dcut]
            v = v[d < dcut]
        if kde:
            num = np.array([self.dbin(diam, self.binwidth,
                                      data=data).vel_v.count() for diam in d])
            d = d[num > bin_num_min]
            v = v[num > bin_num_min]
            sig = [self.dbin(diam, self.binwidth, data=data).vel_v.sem() for diam in d]
            sigargs = {'sigma': sig, 'absolute_sigma': True}
        else:
            sig = np.ones(d.size)
            sigargs = {}
        if use_curve_fit:
            params, pcov = curve_fit(fit.func, d, v, **sigargs)
        else:
            result = minimize(fit.cost, fit.quess, method='Nelder-Mead', args=(d, v, sig))
            params = result.x
        fit.params = params
        fit.x = d
        fit.y = v
        return fit

    def find_fits(self, rule, varinterval=True, draw_plots=False, empty_on_fail=True, **kwargs):
        print('Calculating velocity fits for given sampling frequency...')
        names = []
        fits = []
        for name, group in self.grouped(rule=rule, varinterval=varinterval):
            try:
                newfit = copy.deepcopy(self.find_fit(data=group, name=name,
                                                     **kwargs))
            except RuntimeError as err:
                print('%s: %s' % (name, err))
                print('Particle count: %s' % group.vel_v.count())
                if len(fits) == 0 or empty_on_fail:
                    print('Using an empty fit')
                    newfit = self.default_fit
                else:
                    print('Using fit from previous time step.')
                    newfit = fits[-1]
            fits.append(newfit)
            names.append(name)
            if draw_plots:
                f, ax = plt.subplots()
                self.plot_kde(data=group, ax=ax)
                self.plot(data=group, hexbin=False, ax=ax)
                plt.show()
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
        if fit_type is None:
            fit_type = self.default_fit.name
        letter = 'abcdef'
        params = self.fits[fit_type].apply(lambda fit: fit.params)
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
        values = np.vstack([d, v])
        return stats.gaussian_kde(values)

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

    def plot_fit(self, tstep=None, **kwargs):
        if tstep is None:
            fits = [self.find_fit()]
        else:
            fits = self.fits.loc[tstep].values
        for fit in fits:
            fit.plot(**kwargs)

    def plot(self, data=None, hexbin=True, ax=None, xmax=None, ymax=None,
             show_particle_count=False, colormap='gray_r', ygrid=True,
             hexsize=12, **kwargs):
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
        #gunn_kinzer.plot(dmax=20, label='Gunn&Kinzer', ax=ax, zorder=5, ls='--')
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

class PIPpart(InstrumentData):
    """PIP particle tables"""
    def __init__(self, filenames, dt_start=None, dt_end=None, **kwargs):
        print('Reading PIP particle data...')
        InstrumentData.__init__(self, filenames, **kwargs)
        self.name = 'pip_part'
        dtype = {'Year':np.int32, 'Month':np.int32, 'Day':np.int32,
                 'Hr':np.int32, 'Min':np.int32, 'Sec':np.int32}
        if self.data.empty:
            for filename in filenames:
                newdata = pd.read_csv(filename, delim_whitespace=True,
                                      skiprows=[0, 1, 2, 3, 4, 5, 6, 7, 9],
                                      index_col='datetime',
                                      parse_dates={'datetime':['Year', 'Month', 'Day', 'Hr', 'Min', 'Sec']},
                                      date_parser=datetime.datetime,
                                      dtype=dtype)
                self.data = self.data.append(newdata)
            self.finish_init(dt_start, dt_end)
