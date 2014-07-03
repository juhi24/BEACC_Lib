import time
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import datetime
import pandas as pd
import linecache
from os import path
import copy

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
        
class Pluvio(InstrumentData):
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
        self.buffer = pd.datetools.timedelta(minutes=10)
        if dt_start is None or dt_end is None:
            self.buffer = pd.datetools.timedelta(0)
        elif dt_start-self.buffer < self.data.index[0] or dt_end+self.buffer > self.data.index[-1]:
            self.buffer = pd.datetools.timedelta(0)
        self.data = self.data[dt_start-self.buffer:dt_end+self.buffer]
        
    def timeshift(self):
        if self.shift_periods == 0:
            return pd.datetools.timedelta(0)
        return self.shift_periods*pd.datetools.to_offset(self.shift_freq)
        
    def dt_start(self):
        return self.data.index[0] + self.buffer #+ self.timeshift()
        
    def dt_end(self):
        return self.data.index[-1] - self.buffer #+ self.timeshift()
        
    def shift_reset(self):
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
        return self.data.bucket_nrt-self.data.bucket_nrt[0]
        
    def noprecip_bias(self, lwc, inplace=False):
        """Calculate accumulated bias using LWC."""
        accum = self.acc(rule='1min', unbias=False)
        lwc_filled = lwc.reindex(accum.index).fillna(0)
        bias_acc = accum.diff()[lwc_filled==0].cumsum()
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
            sorted_column_index = list(map(str,(sorted(self.bin_cen())))) # trust me
            self.data = self.data.reindex_axis(sorted_column_index, axis=1)
        self.finish_init(dt_start, dt_end)

    def parse_datetime(self, hh, mm):
        dateline = linecache.getline(self.current_file, 6)
        datearr = [int(x) for x in dateline.split()]
        d = datetime.date(*datearr)
        t = datetime.time(int(hh), int(mm))
        return datetime.datetime.combine(d, t)
        
    def bin_cen(self):
        """Return size bin centers in numeric format."""
        return list(map(float, self.data.columns))
        
    def plot(self, img=True):
        """Plot particle size distribution over time."""
        if img:
            plt.matshow(self.data.transpose(), norm=LogNorm(), origin='lower')
        else:
            plt.pcolor(self.data.transpose(), norm=LogNorm())
        plt.colorbar()
        plt.title('PIP DSD')
        plt.xlabel('time (UTC)')
        plt.ylabel('D (mm)')
        
    def corrected_data(self):
        return 2*self.data
        
class PipV(InstrumentData):
    """PIP particle velocity and diameter data handling"""
    def __init__(self, filenames, dt_start=None, dt_end=None, **kwargs):
        """Create a PipV object using data from a list of PIP velocity table files."""
        print('Reading PIP particle velocity data...')
        InstrumentData.__init__(self, filenames, **kwargs)
        self.name = 'pip_vel'
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
        
    def lwc(self, rule='1min'):
        d3 = self.data.Wad_Dia**3
        return d3.resample(rule, how=self._sum, closed='right', label='right')
    
    def parse_datetime(self, mm):
        datestr = self.current_file.split('/')[-1].split('_')[0]
        yr = int(datestr[3:7])
        mo = int(datestr[7:9])
        dd = int(datestr[9:11])
        hh = int(datestr[11:13])
        return datetime.datetime(yr, mo, dd, hh, int(mm))