"""tools for reading and working with baecc data"""
import time
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import datetime
import pandas as pd
import numpy as np
import linecache
from os import path
import copy
from scipy import stats
from scipy.optimize import curve_fit, fmin

GUNN_KINZER = (9.65, 10.30/9.65, 0.6)

def v_fit(d, a, b, c):
    return a*(1 - b*np.exp(-c*d))
    
def v_polfit(d, a, b):
    return a*d**b

def v_gunn_kinzer(d):
    """d in mm --> return v in m/s"""
    return v_fit(d, *GUNN_KINZER)
    
def plot_gunn_kinzer(dmax, samples=100, ax=None, **kwargs):
    """Plot Gunn&Kinzer v(d) relation."""
    return plot_v_fit(*GUNN_KINZER, func=v_fit, dmax=dmax, samples=samples, ax=ax, **kwargs)
    
def plot_v_fit(*args, func=v_fit, dmax=20, samples=100, ax=None, **kwargs):
    """Plot Gunn&Kinzer shape v(d) relation with custom parameters."""
    if ax is None:
        ax = plt.gca()
    diam = np.linspace(0,dmax,samples)
    vel = [func(d, *args) for d in diam]
    ax.plot(diam, vel, **kwargs)
    return ax

class PrecipMeasurer:
    """parent for classes with precipitation measurements
    Either rainrate or acc (or both) methods should be overridden."""
    def __init__(self):
        pass
    
    def rainrate(self, **kwargs):
        return self.acc(**kwargs).diff()
    
    def acc(self, **kwargs):
        return self.rainrate(**kwargs).cumsum()
        
    def intensity(self, **kwargs):
        """precipitation intensity in mm/h"""
        r = self.rainrate(**kwargs)
        frac = pd.datetools.timedelta(hours=1)/r.index.freq.delta
        return frac * r

class InstrumentData:
    """Parent for instrument data classes."""
    def __init__(self, filenames, hdf_table=None):
        """Read from either instrument output data file or hdf5."""
        self.filenames = filenames
        self.data = pd.DataFrame()
        if hdf_table is not None:
            self.name = hdf_table
            self.data = self.data.append(pd.read_hdf(filenames[0], hdf_table))
            
    def finish_init(self, dt_start, dt_end):
        """Sort and name index, cut time span."""
        self.data.sort_index(inplace=True)
        self.data.index.names = ['datetime']
        self.set_span(dt_start, dt_end)
        
    def parse_datetime(self):
        """Parse timestamps in data files. Used by class constructor."""
        pass
    
    def good_data(self):
        """Return useful data with filters and corrections applied."""
        return self.data
        
    @staticmethod    
    def _sum(x):
        if len(x) < 1: return 0
        else: return sum(x)
        
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
        
class Pluvio(InstrumentData, PrecipMeasurer):
    """Pluviometer data handling"""
    def __init__(self, filenames, dt_start=None, dt_end=None, **kwargs):
        """Create a Pluvio object using data from a list of files."""
        print('Reading pluviometer data...')
        InstrumentData.__init__(self, filenames, **kwargs)
        self.bias = 0
        self.shift_periods = 0
        self.shift_freq = None
        if self.data.empty:
            self.name = path.basename(path.dirname(self.filenames[0])).lower()
            self.col_description = ['date string',
                    'intensity RT [mm h]',
                    'accumulated RT NRT [mm]',
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
                            skiprows=list(range(1,num_lines,2)),
                            parse_dates={'datetime':['datestr']},
                            date_parser=self.parse_datetime,
                            index_col='datetime'))
            self.data.drop(['i_rt'], 1, inplace=True) # crap format
        self.buffer = pd.datetools.timedelta(0)
        self.finish_init(dt_start, dt_end)
        
    def parse_datetime(self, datestr, include_sec=False):
        datestr = str(int(datestr))
        t = time.strptime(datestr, '%Y%m%d%H%M%S')
        if include_sec:
            t_end = 6
        else:
            t_end = 5
        return datetime.datetime(*t[:t_end])
    
    def set_span(self, dt_start, dt_end):
        if dt_start is None or dt_end is None:
            super().set_span(dt_start, dt_end)
            return
        for dt in [dt_start, dt_end]:
            dt = pd.datetools.to_datetime(dt)
        self.buffer = pd.datetools.timedelta(minutes=15)
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
        return self.data.index[0] + self.buffer #+ self.timeshift()
        
    def dt_end(self):
        return self.data.index[-1] - self.buffer #+ self.timeshift()
        
    def shift_reset(self):
        """Reset time shift."""
        self.shift_periods = 0
        self.shift_freq = None
        
    def rainrate(self, rule='1H', upsample=True, **kwargs):
        """Calculate rainrate for given time rule"""
        if upsample:
            acc_1min = self.acc('1min', **kwargs)
        else:
            acc_1min = self.acc_raw()
        r = acc_1min.diff().resample(rule, how=self._sum, closed='right', label='right')
        if not upsample:
            return r.fillna(0)
        t_r0 = r.index[0]
        r[0] = acc_1min[t_r0]-acc_1min[0]
        return r
        
    def acc(self, rule='1H', interpolate=True, unbias=True, shift=True):
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
        return accum.resample(rule, how='last', closed='right', label='right')
                
    def acc_raw(self):
        """accumulation from raw data"""
        return self.data.bucket_nrt-self.data.bucket_nrt[0]
        
    def noprecip_bias(self, lwc, inplace=False):
        """Calculate accumulated bias using LWC."""
        accum = self.acc(rule='1min', unbias=False)
        lwc_filled = lwc.reindex(accum.index).fillna(0)
        bias_acc = accum.diff()[lwc_filled==0].cumsum()
        if bias_acc.empty:
            if inplace:
                self.bias = 0
            return 0
        bias_acc_filled = bias_acc.reindex(accum.index).asfreq('1min').fillna(method='bfill').fillna(method='ffill')
        if inplace:
            self.bias = bias_acc_filled
        return bias_acc_filled

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
                self.data = self.data.append(pd.read_csv(filename, 
                                delim_whitespace=True, skiprows=8, header=3,
                                parse_dates={'datetime':['hr_d','min_d']},
                                date_parser=self.parse_datetime,
                                index_col='datetime'))
            #self.num_d = self.data[['Num_d']]
            # 1st size bin is crap data, last sometimes nans
            self.data.drop(['day_time', 'Num_d', 'Bin_cen', '0.125'], 1, inplace=True)
            self.data.columns = pd.Index([float(i) for i in self.data.columns])
            self.data.sort_index(axis=1)
        self.finish_init(dt_start, dt_end)

    def parse_datetime(self, hh, mm):
        dateline = linecache.getline(self.current_file, 6)
        datearr = [int(x) for x in dateline.split()]
        d = datetime.date(*datearr)
        t = datetime.time(int(hh), int(mm))
        return datetime.datetime.combine(d, t)
        
    def plot(self, img=True):
        """Plot particle size distribution over time."""
        if img:
            plt.matshow(self.good_data().transpose(), norm=LogNorm(), origin='lower')
        else:
            plt.pcolor(self.good_data().transpose(), norm=LogNorm())
        plt.colorbar()
        plt.title('PIP DSD')
        plt.xlabel('time (UTC) BROKEN')
        plt.ylabel('D (mm) BROKEN')
        
    def filter_cat_and_dog(self, window=5):
        """a rolling window filter for isolated data points"""
        is_dog = pd.rolling_count(self.data.mask(self.data==0).T, window).T == 1
        is_dog.ix[:,:window] = False # ignore first columns
        filtered = copy.deepcopy(self.data)
        filtered[is_dog] = 0
        return filtered
        
    def good_data(self, **kwargs):
        return 2*self.filter_cat_and_dog(**kwargs)
        
class PipV(InstrumentData):
    """PIP particle velocity and diameter data handling"""
    def __init__(self, filenames, dt_start=None, dt_end=None, **kwargs):
        """Create a PipV object using data from a list of PIP velocity table files."""
        print('Reading PIP particle velocity data...')
        InstrumentData.__init__(self, filenames, **kwargs)
        self.name = 'pip_vel'
        self.dmin = 0.375 # smallest diameter where data is good
        self._fit_func = v_fit
        self._fit_params = None
        self.D = None
        self.V = None
        self.Z = None
        if self.data.empty:
            for filename in filenames:
                self.current_file = filename
                newdata = pd.read_csv(filename,
                                        delim_whitespace=True, skiprows=8,
                                        parse_dates={'datetime':['minute_p']},
                                        date_parser=self.parse_datetime)
                newdata.rename_axis(
                    {'vel_v_1':'vel_v', 'vel_h_1':'vel_h', 
                     'vel_v_2':'vel_v', 'vel_h_2':'vel_h'}, axis=1, inplace=True)
                self.data = self.data.append(newdata)
            self.data.set_index(['datetime', 'Part_ID', 'RecNum'], inplace=True)
            self.data = self.data.groupby(level=['datetime','Part_ID']).mean()
            self.data.reset_index(level=1, inplace=True)
        self.finish_init(dt_start, dt_end)
        
    @property
    def fit_func(self):
        return self._fit_func
        
    @fit_func.setter
    def fit_func(self, func):
        self._fit_func = func
        self.find_fit()
        
    @property   
    def fit_params(self):
        if self._fit_params is None:
            self.find_fit()
        return self._fit_params
    
    @fit_params.setter
    def fit_params(self, params):
        self._fit_params = params
    
    @fit_params.deleter
    def fit_params(self):
        self._fit_params = None
        
    def v(self, d):
        return self.fit_func(d, *self.fit_params)
        
    def lwc(self, rule='1min'):
        """liquid water content"""
        d3 = self.good_data().Wad_Dia**3
        return d3.resample(rule, how=self._sum, closed='right', label='right')
    
    def parse_datetime(self, mm):
        datestr = self.current_file.split('/')[-1].split('_')[0]
        yr = int(datestr[3:7])
        mo = int(datestr[7:9])
        dd = int(datestr[9:11])
        hh = int(datestr[11:13])
        return datetime.datetime(yr, mo, dd, hh, int(mm))
        
    def good_data(self):
        return self.data[self.data.Wad_Dia>self.dmin]
        
    def dbin(self, d, binwidth):
        """Return data that falls into given size bin."""
        vcond = 'Wad_Dia > %s and Wad_Dia < %s' % (d-0.5*binwidth, d+0.5*binwidth)
        return self.good_data().query(vcond)
        
    def frac_larger(self, d):
        """Return fraction of particles that are larger than d."""
        vdata = self.good_data()
        return vdata[vdata.Wad_Dia>d].vel_v.count()/vdata[vdata.Wad_Dia<d].vel_v.count()
        
    def d_cut(self, frac=0.05, d0=2):
        """Return d for which given fraction of particles are larger."""
        dcost = lambda d: abs(self.frac_larger(d[0])-frac)
        return fmin(dcost, d0)[0]
        
    def find_fit(self, kde=True, cut_d=True, **kwargs):
        """Find and store a Gunn&Kinzer shape fit for either raw data or kde.
        """
        if kde:
            d, v = self.kde_peak()
        else:
            d = self.good_data().Wad_Dia.values
            v = self.good_data().vel_v.values
        if cut_d:
            dcut = self.d_cut(**kwargs)
            d = d[d<dcut]
            v = v[d<dcut]
        self.fit_params, pcov = curve_fit(self.fit_func, d, v)
        return self.fit_params, pcov
        
    def kde(self):
        """kernel-density estimate of d,v data using gaussian kernels"""
        d = self.good_data().Wad_Dia.values
        v = self.good_data().vel_v.values
        values = np.vstack([d, v])
        return stats.gaussian_kde(values)
        
    def set_kde_grid(self, resolution=100):
        """Calculate and store kernel-density estimate with given resolution.
        """
        d = self.good_data().Wad_Dia.values
        v = self.good_data().vel_v.values
        X, Y = np.meshgrid(np.linspace(d.min(),d.max(),resolution), 
                           np.linspace(v.min(),v.max(),resolution))
        points = np.vstack([X.ravel(), Y.ravel()])
        kernel = self.kde()       
        Z = np.reshape(kernel(points).T, X.shape)
        self.D = X
        self.V = Y
        self.Z = Z
        return X, Y, Z
        
    def kde_peak(self):
        """the most propable velocities for each diameter in grid"""
        if self.Z is None:
            self.set_kde_grid()
        x = self.D[0,:]
        y = self.V[:,0][self.Z.argmax(axis=0)]
        return x, y
        
    def plot_kde(self, ax=None):
        """Plot kde grid."""
        if ax is None:
            ax = plt.gca()
        if self.Z is None:
            self.set_kde_grid()
        pc = ax.pcolor(self.D,self.V,self.Z, 
                  cmap=plt.cm.gist_earth_r)
        return pc
        
    def plot_fit(self, **kwargs):
        plot_v_fit(*self.fit_params, func=self.fit_func, **kwargs)
        
    def plot(self, ax=None, style=',', label='pip raw', **kwargs):
        """Plot datapoints and kde."""
        if ax is None:
            ax = plt.gca()
        self.plot_kde(ax=ax)
        self.good_data().plot(x='Wad_Dia', y='vel_v', ax=ax, style=style, 
                              label=label, alpha=0.2, color='black', 
                              **kwargs)
        partcount = self.good_data().Part_ID.count()
        margin = 0.1
        xmax = self.good_data().Wad_Dia.max() + margin
        ymax = self.good_data().vel_v.max() + margin
        ax.axis([0, xmax, 0, ymax])
        ax.set_title('PIP velocity data %s' % self.good_data().index[0].date().isoformat())
        ax.text(0.1, round(ymax)-1, 'particle count: %s' % str(partcount))
        ax.set_xlabel('D (mm)')
        ax.set_ylabel('Vertical velocity (m/s)')
        return ax