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
#batch_create_hdf(datadir='../DATA', outname='baecc3.h5',dtstr='20140212')
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

folder = '/home/dori/SnowCases_BAEC/DensityJussi/corr082_nomu/'
#RadarVP.to_csv(folder + 'radar_data.csv')
huang = False
e = EventsCollection('cases/cases_of_interest_radar.csv', dtformat_default_year)
#e = EventsCollection('cases/test2.csv', dtformat_default_year)
if huang:
    e.autoimport_data(autoshift=False, autobias=False, rule='5min', unbias=True,
                  varinterval=False, radar=True, datafile=['../DATA/baecc3.h5'])
else:
    e.autoimport_data(autoshift=False, autobias=False, rule='5min',
                  varinterval=True, radar=True, datafile=['../DATA/baecc3.h5'])

def func_all_beta(mu,delta,deltaZ):
    def func_beta(beta):
        m=mu
        #if mu < 0.0:
        m= 0.0
        d=delta
        dZ=10**(deltaZ/20.0)
        return sqrt(gamma(2.0*beta+m+1.0)/gamma(m+7.0))*(gamma(d+m+4.0)/gamma(beta+d+m+1.0))-dZ
    return func_beta

def func_all_b(mu,delta,deltaZ):
    def func_b(b):
        m=mu
        #if mu < 0.0:
        m= 0.0
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
            bulk = mat['bulk']
            day = datetime.strptime(huangfile[1:6],'%y%j')
            td = pd.to_timedelta(t_cnt[0]*3600+120+30,unit='s')
            time = day + td
            tmp_b = pd.DataFrame(np.concatenate((3+ba.T, ba.T, 1000.0*bulk.T),axis=1),index=time,columns=['beta_huang','b_huang','density'])
            b_huang = b_huang.append(tmp_b)
b_huang.sort_index(inplace=True)
b_huang = b_huang.drop_duplicates()

           


#%%
Zray = pd.Series()
for c in np.append(e.events.pluvio200.values,e.events.pluvio400.values):
#for c in e.events.pluvio200.values:
    #c.pluvio.lwc = c.pipv.lwc(rule=c.rule)
    #c.pluvio.shift_periods = -6
    basename = folder + datetime.strftime(c.pluvio.dt_start().to_datetime(),'%Y%m%d%H')
    print(datetime.strftime(c.pluvio.dt_start().to_datetime(),'%Y%m%d%H'))
    start = c.pluvio.dt_start()
    if c.varinterval:
        if start.month == 2:
            if start.day == 12:
                c.pluvio.shift_periods = -5
                c.xsacr.time_lag = pd.to_timedelta(60.0*2.0,unit='s')
                c.pluvio.n_combined_intervals = 1
#                c.xsacr.time_lag = pd.to_timedelta(60.0*2.0,unit='s') #fixme
            elif start.day == 15 or start.day == 16:
                c.pluvio.shift_periods = -5
                c.xsacr.time_lag = pd.to_timedelta(60.0*4.0,unit='s')
#                c.pluvio.n_combined_intervals = 1
                c.xsacr.time_lag = pd.to_timedelta(60.0*3.8,unit='s')
            elif start.day == 21 or start.day ==22:
                c.pluvio.shift_periods = -6
                c.xsacr.time_lag = pd.to_timedelta(60.0*1.0,unit='s')
                c.pluvio.n_combined_intervals = 1

    print('estraggo Z')
    zx = 10.0*np.log10(c.z('XSACR'))
    zk = 10.0*np.log10(c.z('KASACR'))
    zkz = 10.0*np.log10(c.z('KAZR'))
    zmw = 10.0*np.log10(c.z('MWACR'))
    plt.figure()
    zx.plot(label='xsacr')
    print('densitá')
    print(c.density(pluvio_filter=False))
    if huang:
        print('Estimate huang reflectivities')
        rho = b_huang['density']
        rhocut = rho.loc[c.pluvio.dt_start():c.pluvio.dt_end()]
        zxtm = c.tmatrix(wl=30.8,pluvio_filter=False,density=rhocut)
        #zktm = c.tmatrix(wl=8.49,pluvio_filter=False,density=rhocut)
        #zwtm = c.tmatrix(wl=3.15,pluvio_filter=False,density=rhocut)
        zxray = c.Z_rayleigh_Xband(density=rhocut)
    else:
        print('TM')
        zxtm = c.tmatrix(wl=30.8,pluvio_filter=False)
        #zktm = c.tmatrix(wl=8.49,pluvio_filter=False)
        #zwtm = c.tmatrix(wl=3.15,pluvio_filter=False)
        print('RAY')
        zxray = c.Z_rayleigh_Xband()
    plt.figure()
    zx.plot(label='xsacr')
    zxtm.plot(label='Xtm')
    zxray.plot(label='Xray')
    plt.legend(loc=0)
    plt.savefig(basename + c.pluvio.name + 'radar_avgX.png')
    plt.close('all')
#    plt.figure()
#    zk.plot(label='kasacr')
#    zkz.plot(label='kazr')
#    zktm.plot(label='Ktm')
#    plt.legend(loc=0)
#    plt.savefig(basename + c.pluvio.name + 'radar_avgK.png')
#    plt.close('all')
#    plt.figure()
#    zmw.plot(label='mwacr')
#    zwtm.plot(label='Wtm')
#    plt.legend(loc=0)
#    plt.savefig(basename + c.pluvio.name + 'radar_avgW.png')
#    plt.close('all')
    
#    depth = c.amount(params=[100],simple=True)
#    accum = c.pluvio.amount(rule=c.rule)
#    xsacr = c.z(radarname = 'XSACR')
#    xsacr = xsacr.fillna(xsacr.interpolate(method='cubic'))
#    xsacr = 10.0*np.log10(xsacr)
#    kasacr = 10.0*np.log10(c.z(radarname = 'KASACR'))
#    kazr = 10.0*np.log10(c.z(radarname = 'KAZR'))
#    mwacr = 10.0*np.log10(c.z(radarname = 'MWACR'))
#    
#    cor_plu_pip = np.correlate(depth,accum,'same')
#    cor_plu_xsa = np.correlate(accum,xsacr,'same')
#    cor_xsa_pip = np.correlate(xsacr,depth,'same')
#    cor_xsa_ksa = np.correlate(xsacr,kasacr,'same')
#    cor_xsa_mwa = np.correlate(xsacr,mwacr,'same')
#    if depth.size % 2:
#        print('odd  ',depth.size,cor_plu_pip.argmax(),cor_plu_pip.argmax()-(depth.size-1)//2)
#        print('odd  ',depth.size,cor_xsa_pip.argmax(),cor_xsa_pip.argmax()-(depth.size-1)//2)
#    else:
#        print('even ',depth.size,cor_plu_pip.argmax(),cor_plu_pip.argmax()-depth.size//2)
#        print('even ',depth.size,cor_xsa_pip.argmax(),cor_xsa_pip.argmax()-depth.size//2)
#        
#    if accum.size % 2:
#        print('odd  ',accum.size,cor_plu_xsa.argmax(),cor_plu_xsa.argmax()-(accum.size-1)//2)
#    else:
#        print('even ',accum.size,cor_plu_xsa.argmax(),cor_plu_xsa.argmax()-accum.size//2)
    #depth.to_csv(basename + 'depth_' + c.pluvio.name + '.csv')
    #print(depth.sum())
    #print(depth)
    #print(c.pluvio.amount(rule=c.rule))

####
    #c.tmatrix(wl=30.89598)
#    Zray = Zray.append(c.Z_rayleigh_Xband())
    #den_dtfr = c.density(pluvio_filter=True,pip_filter=False)
    #den_dtfr.to_csv(basename + 'density_' + c.pluvio.name + '.csv')
    #c.pipv.plots(save=True, suffix='.eps', grid=False, xmax=4, ymax=3, xticks=[0,1,2,3,4], yticks=[0,1,2,3],colorbar=False, hexsize=8)
#    c.summary().to_csv(basename + 'summary_' + c.pluvio.name + '.csv')
#    c.pluvio.tdelta().to_csv(basename + 'timedelta_' + c.pluvio.name + '.csv')
    c.density(pluvio_filter=True,pip_filter=False).plot()
    b_huang.density[c.pluvio.dt_start():c.pluvio.dt_end()].plot()
    axes=plt.gca()
    axes.set_ylim([0, 800])
    plt.savefig(basename + 'density_' + c.pluvio.name + '.png')
    plt.close("all")

    
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
#    Zestimate = pd.DataFrame(Zfile['Z'].values,index=time_Zavg.index,columns=['Zest'])

    #Zestimate = pd.DataFrame(zxtm)
    Zestimate = pd.DataFrame(zxtm)
    time_Zavg = pd.DataFrame(zx)
    dataZ = pd.concat([time_Zavg, c.mu(), c.lam(),c.pipv.fit_params()['b'],c.density(),Zestimate], join='outer', axis = 1)
    FinalData = pd.DataFrame()
    for index, row in dataZ.iterrows():
        dZ = row['XSACR reflectivity']-row[Zestimate.columns[0]]
        f = func_all_beta(mu=row['mu'],delta=row['b'],deltaZ=dZ)
        g = func_all_b(mu=row['mu'],delta=row['b'],deltaZ=dZ)
        if huang:
            f = func_all_beta(mu=row['mu'],delta=0.0,deltaZ=dZ)
            g = func_all_b(mu=row['mu'],delta=0.0,deltaZ=dZ)
        if np.isnan(dZ) or f(0.5)*f(4.5) > 0.0:
            betaopt = np.nan
            bopt = np.nan
        else:
            betaopt = opt.brentq(f,a=0.5,b=4.5,xtol=1.0e-04)
            bopt = opt.brentq(g,a=-2.5,b=1.5,xtol=1.0e-04)
        print(dZ,betaopt,bopt,row['b'],row['mu'])
        tmpFinalData = pd.DataFrame(np.array([[betaopt,bopt]]),index=[index.to_datetime()],columns=['beta','b'])
        FinalData = FinalData.append(tmpFinalData)
    plt.figure()
    axe = FinalData['beta'].plot(marker='*')
    b_huang['beta_huang'][c.pluvio.dt_start().to_datetime():c.pluvio.dt_end().to_datetime()].plot(ax=axe)
    plt.title(datetime.strftime(c.pluvio.dt_start().to_datetime(),'%Y %m %d'))
    
    plt.savefig(basename + 'beta_b_' + c.pluvio.name + '.png')
    plt.close("all")
    

#c.plot_velfitcoefs(rhomax=600, countmin=2000)