from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from datetime import timedelta
import os
import pandas as pd
import scipy.io as matio
import scipy.optimize as opt
from scipy.special import gamma
from scipy import sqrt
import matplotlib.pyplot as plt
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20141217')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20141218')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20141220')
#####batch_create_hdf(datadir='../DATA', outname='test_winter1415.h5',dtstr='20141224')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20141227')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20141228')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20141230')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20150103')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20150108')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20150109')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20150110')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20150111')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20150112')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20150113')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20150114')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20150116')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20150118')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20150122')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20150123')
#batch_create_hdf(datadir='../DATA', outname='winter.h5',dtstr='20150125')

###batch_create_hdf(datadir='../DATA', outname='new_winter.h5',dtstr='20150127') #NODATA
###batch_create_hdf(datadir='../DATA', outname='new_winter.h5',dtstr='20150131') #NODATA
###batch_create_hdf(datadir='../DATA', outname='new_winter.h5',dtstr='20150201') #NODATA
#batch_create_hdf(datadir='../DATA', outname='baecc2.h5',dtstr='20140131')
#batch_create_hdf(datadir='../DATA', outname='baecc2.h5',dtstr='20140201')
#batch_create_hdf(datadir='../DATA', outname='baecc2.h5',dtstr='20140212')
#exit()

radarfolder='../../Radars/'
RadarVP = pd.DataFrame()
for root, dirs, files in os.walk(radarfolder):
    for radarfile in files:
        if radarfile.startswith('tmpxsacr') and radarfile.endswith('.nc'):
            radardata = matio.netcdf.netcdf_file(root + '/' + radarfile)
            radarvariables = radardata.variables
            if 'RHI' in radardata.scan_name.decode():
                print('RHI')
                range_idx = 1
            if 'v' in radardata.scan_name.decode():
                print('VPT')
                range_idx = 0
            print(radardata.scan_name.decode())
            if 'reflectivity' in radarvariables.keys():
                reflectivity = radarvariables['reflectivity'].data[:,range_idx]*radarvariables['reflectivity'].scale_factor + radarvariables['reflectivity'].add_offset
            if 'time' in radarvariables.keys():
                basetime = datetime.strptime(radarvariables['time'].units.decode(),'seconds since %Y-%m-%dT%H:%M:%SZ')
                time_lag = timedelta(minutes=4.4)
                if basetime > datetime(2014,2,15):
                    time_lag = timedelta(minutes=3.8)
                if basetime > datetime(2014,2,21):
                    time_lag = timedelta(minutes=1.0)
                deltatime = pd.to_timedelta(radarvariables['time'].data,unit='s')
                time = basetime + deltatime + time_lag
            if 'elevation' in radarvariables.keys():
                elevation = radarvariables['elevation'].data
            if 'range' in radarvariables.keys():
                rng = radarvariables['range'].data
            print(radarfile,rng[0])
            VP = np.abs(elevation-90) < 0.5
            tmpDF = pd.DataFrame(reflectivity[VP],index=time[VP])
            RadarVP = RadarVP.append(tmpDF)
RadarVP.sort_index(inplace=True)

dtformat_default = '%d.%m. %H:%M'
dtformat_default_year = '%d.%m.%y %H:%M'
dtformat_snex = '%y %d %B %H UTC'
dtformat_print = '%y%m%d%H%M'

folder = '/home/dori/SnowCases_BAEC/DensityJussi/beta/'

e = EventsCollection('cases/cases_of_interest_radar.csv', dtformat_default_year)
e.autoimport_data(autoshift=False, autobias=False, rule='5min', varinterval=True, datafile=['../DATA/baecc.h5'])

#for c in np.append(e.events.pluvio200.values,e.events.pluvio400.values):

def func_all(mu,delta,deltaZ):
    def func_beta(beta):
        m=mu
        #if mu < 0.0:
        #    m= 0.0
        d=delta
        dZ=10**(deltaZ/20.0)
        return sqrt(gamma(2.0*beta+m+1.0)/gamma(m+7.0))*(gamma(d+m+4.0)/gamma(beta+d+m+1.0))-dZ
    return func_beta

for c in e.events.pluvio200.values:
    c.pluvio.shift_periods = -6
    basename = folder + datetime.strftime(c.pluvio.dt_start().to_datetime(),'%Y%m%d%H')
    print(datetime.strftime(c.pluvio.dt_start().to_datetime(),'%Y%m%d%H'))
    #depth = c.amount(params=[100],simple=True)
    #depth.to_csv(basename + 'depth_' + c.pluvio.name + '.csv')
    #print(depth.sum())
    #print(depth)
    #print(c.pluvio.amount(rule=c.rule))
    c.pluvio.n_combined_intervals = 2
#    c.density(pluvio_filter=True,pip_filter=False).to_csv(basename + 'density_' + c.pluvio.name + '.csv')
    c.pipv.plots(save=True, suffix='.eps', grid=False, xmax=4, ymax=3, xticks=[0,1,2,3,4], yticks=[0,1,2,3],colorbar=False, hexsize=8)
#    c.summary().to_csv(basename + 'summary_' + c.pluvio.name + '.csv')
#    c.pluvio.tdelta().to_csv(basename + 'timedelta_' + c.pluvio.name + '.csv')
#    c.density(pluvio_filter=True,pip_filter=False).plot()
#    axes=plt.gca()
#    axes.set_ylim([0, 1000])
#    plt.savefig(basename + 'density_' + c.pluvio.name + '.png')
#    plt.close("all")
    
    # Faccio le medie
    delta = pd.Series(c.pluvio.amount(crop=True).index.to_datetime(),index=c.pluvio.amount(crop=True).index).diff()
    time_Zavg = pd.DataFrame()
    for idx,val in delta.iteritems():
        if pd.isnull(val):
            tmpTimeZ = pd.DataFrame([np.nan],index=[idx.to_datetime()],columns=['Zmea'])
        else:
            RR=RadarVP[idx.to_datetime()-val:idx.to_datetime()]
            #print(RR.shape)
            Zmean = 10*np.log10((10.0**(RR*0.1)).mean()).values
            tmpTimeZ = pd.DataFrame(Zmean,index=[idx.to_datetime()],columns=['Zmea'])
        time_Zavg = time_Zavg.append(tmpTimeZ)
            
    Zfile = pd.read_csv(folder + 'Z_' + c.pluvio.name + '_' + datetime.strftime(c.pluvio.dt_start().to_datetime(),'%Y%m%d%H') + '30.89598.csv',header=None,names=['rho','n','Z'],parse_dates=True)
    
    Zestimate = pd.DataFrame(Zfile['Z'].values,index=time_Zavg.index,columns=['Zest'])
    dataZ = pd.concat([time_Zavg, c.mu(), c.lam(),c.pipv.fit_params()['b'],c.density(),Zestimate], join='outer', axis = 1)
    FinalData = pd.DataFrame()
    for index, row in dataZ.iterrows():
        dZ = row['Zmea']-row['Zest']
        f = func_all(mu=row['mu'],delta=row['b'],deltaZ=dZ)
        if np.isnan(dZ) or f(0.5)*f(4.5) > 0.0:
            betaopt=np.nan
        else:
            betaopt = opt.brentq(f,a=0.5,b=4.5,xtol=1.0e-04)
        print(dZ,betaopt,row['b'],row['mu'])
        tmpFinalData = pd.DataFrame([betaopt],index=[index.to_datetime()],columns=['beta'])
        #,row['mu'],row['b'],row['Zest'],row['Zmea']
        FinalData = FinalData.append(tmpFinalData)
    
    FinalData.plot()
    plt.title(datetime.strftime(c.pluvio.dt_start().to_datetime(),'%Y %m %d'))
    plt.savefig(basename + 'beta_' + c.pluvio.name + '.png')
    plt.close("all")
    
    

#c.plot_velfitcoefs(rhomax=600, countmin=2000)

#03.01.15 00:00,03.01.15 23:50
#12.01.15 21:00,13.01.15 08:00
#31.01.14 21:00,01.02.14 04:00


#14.01.15 01:00,14.01.15 05:00
#16.01.15 01:00,16.01.15 07:00
