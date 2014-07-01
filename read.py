import time
import glob
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import datetime
import pandas as pd
from collections import defaultdict
import linecache
from os import path
    
def hotplate(date):
    #file_ = open(filename)
    date_str = time.strftime("%Y%m%d",date)
    files = glob.glob("../data/Hotplate/hot_plate_100901_"+date_str+"*")
    
    data=defaultdict(list)

    #descriptin of file
    file_format = {
    1: 'Output format',
    2: 'Fault indicator',
    3: 'Unix timestamp',
    4: 'Voltage, sensor, instantaneous (V)',
    5: 'Voltage, reference, instantaneous (V)',
    6: 'Current, sensor, instantaneous (A)',
    7: 'Current, reference, instantaneous (A)',
    8: 'Resistance, sensor, ratio of previous fields (ohm)',
    9: 'Resistance, reference, ratio of previous fields (ohm)',
    10: 'Power, sensor, 1-minute running average, 1 Hz samples (W)',
    11: 'Power, reference, 1-minute running averag    file_ = open(filename)e, 1 Hz samples (W)',
    12: 'Control effort (PWM), sensor, instantaneous (%)',
    13: 'Control effort (PWM), reference, instantaneous (%)',
    14:'Ambient temperature, 1-minute running ave. 1Hz samples (degC)',
    15: 'Enclosure temperature, 1-minute running ave. 1Hz samples (degC)',
    16: 'Solar IR sensor temperature, 1-minute running ave. 1Hz samples (degC)',
    17: 'Solar radiation, 1-minute ave., 1Hz samples (Wm-2)',
    18: 'Net IR radiation ground to sky, 1-minute running ave., 1 Hz samples (Wm-2)',
    19: 'Barometric pressure, referenced to sea level (mbar)',
    20: 'Temperature of humidity sensor (degC)',
    21: 'Relative humidity (%)',
    22: 'Wind speed, 1- minute running ave, 1 Hz samples (m s)',
    23: 'Collection efficiency, 1-minute running ave, 1Hz samples (W)',
    24: 'Power offset, 1-minute running ave., 1Hz samples (W)',
    25: 'Power offset due to radiation effects, 1Hz samples (W)',
    26: 'Raw precip. rate, 1-minute running ave, 1import os.path Hz samples (mm hr)',
    27: 'Power, sensor, 5-minute running average, 1 minute samples (W)',
    28: 'Power, reference, 5-minute running average, 1 minute samples (W)',
    29: 'delta Power, 5-minute running ave, 1 minute samples (W)',
    30: 'Ambient temperature, 5-minute running ave., 1 minute samples (degC)',
    31: 'Power offset, 5-minute ave., 1 minute samples (W)',
    32: 'Raw precip.rate, 5-minute running ave, 1 min. samples (mm hr)',
    33: 'Current precipitation rate (mm hr)',
    34: 'Total accumulated liquid precipitation (mm)'
     }

    for filename in files:
        file_ = open(filename)
        lines = file_.readlines()
        file_.close
        time_ = []
        for i in lines:
            var = i.split(',')
            if len(var) > 36:
                time_tmp = time.strptime(var[0],'%Y%m%d%H%M%S')
                time_tmp = time.mktime(time_tmp)
                firmware = var[2]
                dev_num = var[1]
                
                data['hotplate_time'].append(time_tmp)
                data['hotplate_firmware'].append(firmware)
                data['hotpalte_devnum'].append(dev_num)
                
                for key,value in file_format.iteritems():
                    data["hotplate_"+value].append(var[key+2])
            
    return pd.DataFrame(data,columns=data.keys())

def jenoptik(date):
    date_str = time.strftime("%Y%m%d",date)
    files = glob.glob("../data/Jenoptik/"+date_str+"*")
    snow = []
    time_ = []
    signal=[]
    temp = []
    for filename in files:
        file_ = open(filename)
        lines = file_.readlines()
        file_.close
        for i in lines:
            var = i.split(',')
            if len(var) > 0:
                time_tmp = time.strptime(var[0],'%Y-%m-%d %H:%M:%S')
                time_tmp = time.mktime(time_tmp)
                time_.append(time_tmp)
                snow.append(float(var[1])-0.034)
                signal.append(float(var[2]))
                temp.append(float(var[3]))
    #print acc
    d = {'jenoptik_time' : time_, 'jenoptik_snow_depth': snow,'jenoptik_signal_strength':signal,'jenoptik_temperature':temp}
    return pd.DataFrame(d)

def pluvio(date):
    """DEPRECATED"""
    date_str = time.strftime("%Y%m%d",date)
    files = glob.glob("../data/Pluvio200/pluvio200_02_"+date_str+"*")
    data=defaultdict(list)


    file_format = {
    2: 'Intensity RT  [mm h]',
    3: 'Accumulated RT NRT [mm]',
    4: 'Accumulated NRT [mm]',
    5: 'Accumulated total NRT [mm]',
    6: 'Bucket RT [mm]',
    7: 'Bucket NRT [mm]',
    8: 'Temperature load cell [degC]',
    9: 'Heating status',
    10: 'Status',
    11: 'Temperature electronics unit',
    12: 'Supply Voltage',
    13: 'Temperature orfice ring rim'
    }

    for filename in files:
        file_ = open(filename)
        lines = file_.readlines()
        file_.close
        time_ = []
        for i in lines:
            var = i.split(';')
            if len(var) > 3:
                time_tmp = time.strptime(var[0],'%Y%m%d%H%M%S')
                time_tmp = time.mktime(time_tmp)
                data['pluvio_time'].append(time_tmp)
                print(len(var))
                for key,value in file_format.iteritems():
                    data['pluvio '+value].append(var[key-1])
            
    #data["pluvio_time"] = time_
    #print len(data['status'])
    return pd.DataFrame(data)

class InstrumentData:
    """Parent for instrument data classes."""
    def __init__(self, filenames, hdf_table=None):
        self.filenames = filenames
        self.data = pd.DataFrame()
        if hdf_table is not None:
            self.name = hdf_table
            self.data = self.data.append(pd.read_hdf(filenames[0], hdf_table))
            
    def finish_init(self, dt_start, dt_end):
        self.data.sort_index(inplace=True)
        self.data.index.names = ['datetime']
        self.data = self.data[dt_start:dt_end]
        
    def parse_datetime(self):
        """Parse timestamps in data files. Used by class constructor."""
        pass
        
    @staticmethod    
    def _sum(x):
        if len(x) < 1: return 0
        else: return sum(x)
        
    def to_hdf(self, filename='../DATA/baecc.h5'):
        self.data.to_hdf(filename, self.name, format='table', append=True)
        
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
                    'intensity RT  [mm h]',
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
        self.finish_init(dt_start, dt_end)
        
    def parse_datetime(self, datestr, include_sec=False):
        datestr = str(int(datestr))
        t = time.strptime(datestr, '%Y%m%d%H%M%S')
        if include_sec:
            t_end = 6
        else:
            t_end = 5
        return datetime.datetime(*t[:t_end])
        
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
        if unbias:
            accum -= self.bias
        if shift:
            accum = accum.tshift(periods=self.shift_periods, freq=self.shift_freq)
        return accum.resample(rule, how='last', closed='right', label='right')
                
    def acc_raw(self):
        return self.data.bucket_nrt-self.data.bucket_nrt[0]
        
    def noprecip_bias(self, lwc, inplace=False):
        """Calculate accumulated bias using LWC."""
        accum = self.acc(rule='1min', unbias=False)
        lwc_filled = lwc.reindex(accum.index).fillna(0)
        bias_acc = accum.diff()[lwc_filled==0].cumsum()
        bias_acc_filled = bias_acc.reindex(accum.index).asfreq('1min').fillna(method='bfill')
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
            # 1st size bin is crap data
            self.data.drop(['day_time', 'Num_d', 'Bin_cen', '0.125'], 1, inplace=True)
            sorted_column_index = list(map(str,(sorted(self.bin_cen())))) # trust me
            self.data = self.data.reindex_axis(sorted_column_index, axis=1)
        self.finish_init(dt_start, dt_end)

    def parse_datetime(self, hh, mm):
        dateline = linecache.getline(self.current_file, 6)
        datearr = [int(x) for x in dateline.split()]
        date = datetime.date(*datearr)
        time = datetime.time(int(hh), int(mm))
        return datetime.datetime.combine(date, time)
        
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