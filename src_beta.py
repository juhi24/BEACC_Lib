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

#from pytmatrix.test import test_tmatrix
#test_tmatrix.run_tests()

#batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20141217')
#batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20141218')
#batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20141220')
######batch_create_hdf(datadir='../DATA', outname='test_winter1415.h5',dtstr='20141224')
#batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20141227')
#batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20141228')
###batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20141230')
###batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20150103')
###batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20150108')
###batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20150109')
###batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20150110')
###batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20150111')
###batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20150112')
###batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20150113')
###batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20150114')
###batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20150116')
###batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20150118')
###batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20150122')
###batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20150123')
###batch_create_hdf(datadir='../DATA', outname='winter1.h5',dtstr='20150125')
##
####batch_create_hdf(datadir='../DATA', outname='new_winter.h5',dtstr='20150127') #NODATA
####batch_create_hdf(datadir='../DATA', outname='new_winter.h5',dtstr='20150131') #NODATA
####batch_create_hdf(datadir='../DATA', outname='new_winter.h5',dtstr='20150201') #NODATA
#batch_create_hdf(datadir='../DATA', outname='baecc3.h5',dtstr='20140131')
#batch_create_hdf(datadir='../DATA', outname='baecc3.h5',dtstr='20140201')
batch_create_hdf(datadir='../DATA', outname='baecc3.h5',dtstr='20140212')
#batch_create_hdf(datadir='../DATA', outname='baecc3.h5',dtstr='20140215')
#batch_create_hdf(datadir='../DATA', outname='baecc3.h5',dtstr='20140216')
#batch_create_hdf(datadir='../DATA', outname='baecc3.h5',dtstr='20140221')
#batch_create_hdf(datadir='../DATA', outname='baecc3.h5',dtstr='20140222')
print('FINE FINE FINE')

#%%

dtformat_default = '%d.%m. %H:%M'
dtformat_default_year = '%d.%m.%y %H:%M'
dtformat_snex = '%y %d %B %H UTC'
dtformat_print = '%y%m%d%H%M'

folder = '/home/dori/SnowCases_BAEC/DensityJussi/test/'
#RadarVP.to_csv(folder + 'radar_data.csv')

e = EventsCollection('cases/cases_of_interest_radar.csv', dtformat_default_year)
e.autoimport_data(autoshift=False, autobias=False, rule='5min',
                  varinterval=True, radar=True, datafile=['../DATA/baecc3.h5'])

def func_all_beta(mu,delta,deltaZ):
    def func_beta(beta):
        m=mu
        #if mu < 0.0:
        #    m= 0.0
        d=delta
        dZ=10**(deltaZ/20.0)
        return sqrt(gamma(2.0*beta+m+1.0)/gamma(m+7.0))*(gamma(d+m+4.0)/gamma(beta+d+m+1.0))-dZ
    return func_beta

def func_all_b(mu,delta,deltaZ):
    def func_b(b):
        m=mu
        #if mu < 0.0:
        #    m= 0.0
        d=delta
        dZ=10**(deltaZ/20.0)
        return sqrt(gamma(2.0*b+m+7.0)/gamma(m+7.0))*(gamma(d+m+4.0)/gamma(b+d+m+4.0))-dZ
    return func_b
    
huang_dir = '../huang/'
b_huang = pd.DataFrame()
for root, dirs, files in os.walk(huang_dir):
    for huangfile in files:
        if huangfile.endswith('.mat'):
            mat = matio.loadmat(huang_dir + huangfile)
            #aa = mat['aa']
            ba = mat['ba']
            t_cnt = mat['t_cnt']
            day = datetime.strptime(huangfile[1:6],'%y%j')
            td = pd.to_timedelta(t_cnt[0]*3600+120+30,unit='s')
            time = day + td
            tmp_b = pd.DataFrame(np.concatenate((3+ba.T, ba.T),axis=1),index=time,columns=['beta_huang','b_huang'])
            b_huang = b_huang.append(tmp_b)
b_huang.sort_index(inplace=True)            


#%%
#for c in np.append(e.events.pluvio200.values,e.events.pluvio400.values):
Zray = pd.Series()
for c in e.events.pluvio200.values:
    c.pluvio.shift_periods = -6
    basename = folder + datetime.strftime(c.pluvio.dt_start().to_datetime(),'%Y%m%d%H')
    print(datetime.strftime(c.pluvio.dt_start().to_datetime(),'%Y%m%d%H'))
    c.pluvio.n_combined_intervals = 1
    zx = 10.0*np.log10(c.z('XSACR'))
    zk = 10.0*np.log10(c.z('KASACR'))
    zkz = 10.0*np.log10(c.z('KAZR'))
    zmw = 10.0*np.log10(c.z('MWACR'))
    plt.figure()
    zx.plot(label='xsacr')
    #zk.plot(label='kasacr')
    #zkz.plot(label='kazr')
    #zmw.plot(label='mwacr')
    zxtm = c.tmatrix(wl=30.87)
    zxray = c.Z_rayleigh_Xband()
    zxtm.plot(label='Xtm')
    zxray.plot(label='Xray')
    plt.legend(loc=0)
    plt.savefig(basename + 'radar_avg.png')
    plt.close()
    #depth = c.amount(params=[100],simple=True)
    #depth.to_csv(basename + 'depth_' + c.pluvio.name + '.csv')
    #print(depth.sum())
    #print(depth)
    #print(c.pluvio.amount(rule=c.rule))

####
#    c.pluvio.n_combined_intervals = 2
    #c.tmatrix(wl=30.89598)
#    Zray = Zray.append(c.Z_rayleigh_Xband())
    #den_dtfr = c.density(pluvio_filter=True,pip_filter=False)
    #den_dtfr.to_csv(basename + 'density_' + c.pluvio.name + '.csv')
    #c.pipv.plots(save=True, suffix='.eps', grid=False, xmax=4, ymax=3, xticks=[0,1,2,3,4], yticks=[0,1,2,3],colorbar=False, hexsize=8)
    #c.summary().to_csv(basename + 'summary_' + c.pluvio.name + '.csv')
#    c.pluvio.tdelta().to_csv(basename + 'timedelta_' + c.pluvio.name + '.csv')
#    c.density(pluvio_filter=True,pip_filter=False).plot()
#    axes=plt.gca()
#    axes.set_ylim([0, 1000])
#    plt.savefig(basename + 'density_' + c.pluvio.name + '.png')
#    plt.close("all")

    
    # Faccio le medie
#    delta = pd.Series(c.pluvio.amount(crop=True).index.to_datetime(),index=c.pluvio.amount(crop=True).index).diff()
#    time_Zavg = pd.DataFrame()
#    for idx,val in delta.iteritems():
#        if pd.isnull(val):
#            tmpTimeZ = pd.DataFrame([np.nan],index=[idx.to_datetime()],columns=['Zmea'])
#        else:
#            RR=RadarVP[idx.to_datetime()-val:idx.to_datetime()]
#            #print(RR.shape)
#            Zmean = 10*np.log10((10.0**(RR*0.1)).mean()).values
#            tmpTimeZ = pd.DataFrame(Zmean,index=[idx.to_datetime()],columns=['Zmea'])
#        time_Zavg = time_Zavg.append(tmpTimeZ)
#            
#    Zfile = pd.read_csv(folder + 'Z_' + c.pluvio.name + '_' + datetime.strftime(c.pluvio.dt_start().to_datetime(),'%Y%m%d%H') + '30.89598.csv',header=None,names=['rho','n','Z'],parse_dates=True)
#    
#    Zestimate = pd.DataFrame(Zfile['Z'].values,index=time_Zavg.index,columns=['Zest'])
#    dataZ = pd.concat([time_Zavg, c.mu(), c.lam(),c.pipv.fit_params()['b'],c.density(),Zestimate], join='outer', axis = 1)
#    FinalData = pd.DataFrame()
#    for index, row in dataZ.iterrows():
#        dZ = row['Zmea']-row['Zest']
#        f = func_all_beta(mu=row['mu'],delta=row['b'],deltaZ=dZ)
#        g = func_all_b(mu=row['mu'],delta=row['b'],deltaZ=dZ)
#        if np.isnan(dZ) or f(0.5)*f(4.5) > 0.0:
#            betaopt = np.nan
#            bopt = np.nan
#        else:
#            betaopt = opt.brentq(f,a=0.5,b=4.5,xtol=1.0e-04)
#            bopt = opt.brentq(g,a=-2.5,b=1.5,xtol=1.0e-04)
#        print(dZ,betaopt,bopt,row['b'],row['mu'])
#        tmpFinalData = pd.DataFrame(np.array([[betaopt,bopt]]),index=[index.to_datetime()],columns=['beta','b'])
#        FinalData = FinalData.append(tmpFinalData)
#    
#    axe = FinalData.plot()
#    b_huang[c.pluvio.dt_start().to_datetime():c.pluvio.dt_end().to_datetime()].plot(ax=axe)
#    plt.title(datetime.strftime(c.pluvio.dt_start().to_datetime(),'%Y %m %d'))
#    
#    plt.savefig(basename + 'beta_b_' + c.pluvio.name + '.png')
#    plt.close("all")
    

#c.plot_velfitcoefs(rhomax=600, countmin=2000)

#03.01.15 00:00,03.01.15 23:50
#12.01.15 21:00,13.01.15 08:00
#31.01.14 21:00,01.02.14 04:00


#14.01.15 01:00,14.01.15 05:00
#16.01.15 01:00,16.01.15 07:00
