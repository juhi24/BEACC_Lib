"""Tools for estimating density and other properties of falling snow"""
import numpy as np
import pandas as pd
import read
from datetime import datetime, timedelta
from scipy.optimize import minimize
from scipy.special import gamma
from glob import glob
from itertools import cycle
from multiprocessing import Pool
import matplotlib.pyplot as plt
import copy
import locale
import os
import warnings

from pytmatrix import tmatrix, psd, refractive, radar
from pytmatrix import tmatrix_aux as tm_aux

# general configuration
DEBUG = True

# CONFIG default paths
DATA_DIR = '../DATA'
H5_FILE = 'baecc.h5'
H5_PATH = os.path.join(DATA_DIR, H5_FILE)
PIPV_SUBPATH = 'PIP/a_Velocity_Tables/004%s/*2.dat'
DSD_SUBPATH = 'PIP/a_DSD_Tables/004%s_a_d.dat'
P200_SUBPATH = 'Pluvio200/pluvio200_??_%s*.txt'
P400_SUBPATH = 'Pluvio400/pluvio400_??_%s*.txt'
RADAR_SUBPATH = 'Radar/%s/tmp%s*M1.a1.%s.*'

locale.setlocale(locale.LC_ALL, 'en_GB.UTF-8')

if DEBUG:
    warnings.simplefilter('default')
else:
    warnings.simplefilter('ignore')

TAU = 2*np.pi
RHO_W = 1000

# Wood et al. 2013 Characterization of video disdrometer uncertainties ...
#PIP_CORR = 1.0/0.82
PIP_CORR = 1/0.9    # Davide


def deprecation(message, stacklevel=2):
    """Issue DeprecationWarning"""
    warnings.warn(message, DeprecationWarning, stacklevel=stacklevel)


def switch_wl(x):
    return {tm_aux.wl_C: "C", tm_aux.wl_X: "X", tm_aux.wl_Ku: "Ku",
            tm_aux.wl_Ka: "Ka", tm_aux.wl_W: "W"}.get(x, str(x))


def daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days)):
        yield start_date + timedelta(n)


def datafilelist(subpath, datadir=DATA_DIR):
    return glob(os.path.join(datadir, subpath))


def batch_import(dtstr, datadir=DATA_DIR):
    """Read ASCII data according to a datestring pattern."""
    pipv_files = datafilelist(PIPV_SUBPATH % dtstr, datadir=datadir)
    dsd_files = datafilelist(DSD_SUBPATH % dtstr, datadir=datadir)
    pluvio200_files = datafilelist(P200_SUBPATH % dtstr, datadir=datadir)
    pluvio400_files = datafilelist(P400_SUBPATH % dtstr, datadir=datadir)
    xsacr_files = datafilelist(RADAR_SUBPATH % ('XSACR', 'xsacr', dtstr),
                               datadir=datadir)
    kasacr_files = datafilelist(RADAR_SUBPATH % ('KASACR', 'kasacr', dtstr),
                                datadir=datadir)
    kazr_files = datafilelist(RADAR_SUBPATH % ('KAZR', 'kazrge', dtstr),
                              datadir=datadir)
    mwacr_files = datafilelist(RADAR_SUBPATH % ('MWACR', 'mwacr', dtstr),
                               datadir=datadir)
    pluvio200 = read.Pluvio(pluvio200_files)
    pluvio400 = read.Pluvio(pluvio400_files)
    pipv = read.PipV(pipv_files)
    dsd = read.PipDSD(dsd_files)
    xsacr = read.Radar(xsacr_files)
    kasacr = read.Radar(kasacr_files)
    kazr = read.Radar(kazr_files)
    mwacr = read.Radar(mwacr_files)
    return {'vel': pipv, 'dsd': dsd, 'pluvio200': pluvio200,
            'pluvio400': pluvio400, 'xsacr': xsacr, 'kasacr': kasacr,
            'kazr': kazr, 'mwacr': mwacr}


def batch_create_hdf(datadir=DATA_DIR, outname=H5_FILE,
                     dtstr='20140[2-3]??'):
    """Read ASCII data and export to hdf."""
    instrdict = batch_import(dtstr, datadir)
    hdf_file = os.path.join(datadir, outname)
    for key in instrdict:
        instrdict[key].to_hdf(filename=hdf_file)


def scatterplot(x, y, c=None, kind='scatter', **kwargs):
    """scatter plot of two Series objects"""
    plotdata = read.merge_series(x, y)
    if c is not None:
        kwargs['c'] = c
    return plotdata.plot(kind=kind, x=x.name, y=y.name, **kwargs)


def combine_datasets(*datasets):
    eventslist = []
    for e in datasets:
        eventslist.append(e.events)
    combined = copy.deepcopy(datasets[0])
    combined.events = pd.concat(eventslist)
    combined.events.sort(columns='start', inplace=True)
    combined.events.reset_index(inplace=True)
    return combined


class MultiSeries:
    """Provide calculated parameters as one DataFrame and use it for plotting.
    """
    def __init__(self):
        pass

    def summary():
        pass

    def plot_pairs(self, x='a', y='b', c=None, sizecol=None, scale=1,
                   kind='scatter', grouped=True, pluvio=None, query=None,
                   ax=None, colorbar=False, markers='os^vD*p><',
                   edgecolors='none', dtformat='%Y %b %d', **kwargs):
        """Easily plot parameters against each other."""
        sumkwargs = {}
        if ax is None:
            ax = plt.gca()
        if pluvio is not None:
            sumkwargs['pluvio'] = pluvio
        data = self.summary(dtformat=dtformat, **sumkwargs)
        if query is not None:
            data = data.query(query)
        if c is not None:
            kwargs['c'] = c
        if sizecol is not None:
            kwargs['s'] = scale*np.sqrt(data[sizecol])
        if grouped:
            groups = data.groupby('case')
            for (name, group), marker in zip(groups, cycle(markers)):
                colorbar = groups.case.first().iloc[0] == name and colorbar
                group.plot(ax=ax, x=x, y=y, marker=marker, kind=kind,
                           label=name, colorbar=colorbar,
                           edgecolors=edgecolors, **kwargs)
            return ax
        return data.plot(ax=ax, x=x, y=y, kind=kind, colorbar=colorbar,
                         edgecolors=edgecolors, **kwargs)


class EventsCollection(MultiSeries):
    """Manage multiple snow/rain events."""
    def __init__(self, csv, dtformat='%d %B %H UTC'):
        """Read event metadata from a csv file."""
        self.dtformat = dtformat
        self.events = pd.read_csv(csv, parse_dates=['start', 'end'],
                                  date_parser=self.parse_datetime)
        self.events.sort(columns=['start', 'end'], inplace=True)
        self.events.start += pd.datetools.timedelta(seconds=1)

    def parse_datetime(self, dtstr):
        #date = datetime.strptime(dtstr+'+0000', self.dtformat+'%z')
        date = datetime.strptime(dtstr, self.dtformat)
        return date

    def add_data(self, data, autoshift=True, autobias=True):
        """Add data from a Case object."""
        cases = []
        for (i, e) in self.events.iterrows():
            cases.append(data.between_datetime(e.start, e.end,
                                               autoshift=autoshift,
                                               autobias=autobias))
        self.events[data.instr['pluvio'].name] = cases

    def autoimport_data(self, datafile=[H5_PATH], autoshift=False,
                        autobias=False, radar=False, **casekwargs):
        """Import data from a hdf file."""
        timemargin = pd.datetools.timedelta(hours=3)
        dt_start = self.events.iloc[0].start - timemargin
        dt_end = self.events.iloc[-1].end + timemargin
        data = Case.from_hdf(dt_start, dt_end, autoshift=False,
                             filenames=datafile, radar=radar, **casekwargs)
        for d in data:
            if d is not None:
                self.add_data(d, autoshift=autoshift, autobias=autobias)
        return

    def summary(self, pluvio='pluvio200', dtformat='%Y %b %d', **kwargs):
        sumlist = []
        for c in self.events[pluvio]:
            sumlist.append(c.summary(dtformat=dtformat))
        return pd.concat(sumlist, **kwargs)


class Case(read.PrecipMeasurer, read.Cacher, MultiSeries):
    """Calculate snowfall rate from particle size and velocity data."""
    def __init__(self, dsd, pipv, pluvio, xsacr=None, kasacr=None,
                 kazr=None, mwacr=None, varinterval=True, unbias=False,
                 autoshift=False, liquid=False, quess=(0.01, 2.1),
                 bnd=((0, 0.1), (1, 3)), rule='15min', use_cache=True):
        self._use_cache = use_cache
        if xsacr is None:
            self.instr = {'pluvio': pluvio, 'dsd': dsd, 'pipv': pipv}
        else:
            self.instr = {'pluvio': pluvio, 'dsd': dsd, 'pipv': pipv,
                          'xsacr': xsacr, 'kasacr': kasacr, 'kazr': kazr,
                          'mwacr': mwacr}
        self.instr_depr_args = {'message': 'Please use new syntax: case.instr[instrument_name]',
                                'stacklevel': 3}
        self._varinterval = varinterval
        self.instr['pluvio'].varinterval = varinterval
        self.quess = quess
        self.bnd = bnd
        if varinterval:
            self._rule = None
        else:
            self._rule = rule
        self.liquid = liquid
        self._ab = None         # alpha, beta
        for instr in self.instr.values():
            instr.case = self   # maybe this is kind of silly
        if autoshift:
            self.autoshift()
        if unbias:
            self.noprecip_bias()
        read.Cacher.__init__(self)

    def __repr__(self):
        if self.liquid:
            casetype = 'rain'
        else:
            casetype = 'snow'
        dt_start, dt_end = self.dt_start_end()
        if self.varinterval:
            sampling_label = 'adaptive'
        else:
            sampling_label = self.rule
        return '%s case from %s to %s, %s' % (casetype, dt_start,
                                              dt_end, sampling_label)

    def __add__(self, other):
        combined = copy.deepcopy(self)
        for key in set(list(self.instr.keys()) + list(other.instr.keys())):
            if key in self.instr.keys():
                if key in other.instr.keys():
                    combined.instr[key] = self.instr[key] + other.instr[key]
            elif key in other.instr.keys():
                combined.instr[key] = copy.deepcopy(other.instr[key])
        #combined.clear_cache() # TODO: check if needed
        return combined

    #==========================================================================
    # TODO: remove these getters and setters when no longer needed for
    #       backwards compatibility

    @property
    def dsd(self):
        deprecation(**self.instr_depr_args)
        return self.instr['dsd']

    @dsd.setter
    def dsd(self, dsd):
        deprecation(**self.instr_depr_args)
        self.instr['dsd'] = dsd

    @property
    def pluvio(self):
        deprecation(**self.instr_depr_args)
        return self.instr['pluvio']

    @pluvio.setter
    def pluvio(self, pluvio):
        deprecation(**self.instr_depr_args)
        self.instr['pluvio'] = pluvio

    @property
    def pipv(self):
        deprecation(**self.instr_depr_args)
        return self.instr['pipv']

    @pipv.setter
    def pipv(self, pipv):
        deprecation(**self.instr_depr_args)
        self.instr['pipv'] = pipv

    @property
    def xsacr(self):
        deprecation(**self.instr_depr_args)
        return self.instr['xsacr']

    @xsacr.setter
    def xsacr(self, xsacr):
        deprecation(**self.instr_depr_args)
        self.instr['xsacr'] = xsacr

    @property
    def kasacr(self):
        deprecation(**self.instr_depr_args)
        return self.instr['kasacr']

    @kasacr.setter
    def kasacr(self, kasacr):
        deprecation(**self.instr_depr_args)
        self.instr['kasacr'] = kasacr

    @property
    def kazr(self):
        deprecation(**self.instr_depr_args)
        return self.instr['kazr']

    @kazr.setter
    def kazr(self, kazr):
        deprecation(**self.instr_depr_args)
        self.instr['kazr'] = kazr

    @property
    def mwacr(self):
        deprecation(**self.instr_depr_args)
        return self.instr['mwacr']

    @mwacr.setter
    def mwacr(self, mwacr):
        deprecation(**self.instr_depr_args)
        self.instr['mwacr'] = mwacr

    #==========================================================================

    @property
    def use_cache(self):
        return self._use_cache

    @use_cache.setter
    def use_cache(self, use_cache):
        self._use_cache = use_cache
        for instr in self.instr.values():
            instr.use_cache = use_cache

    @property
    def varinterval(self):
        return self._varinterval

    @varinterval.setter
    def varinterval(self, varinterval):
        self._varinterval = varinterval
        self.instr['pluvio'].varinterval = varinterval
        self.reset()

    @property
    def rule(self):
        if self.varinterval: #and self._rule is None:
            # TODO: needs to be reset on changes for pluvio data
            self._rule = self.instr['pluvio'].grouper()
        return self._rule

    @rule.setter
    def rule(self, rule):
        self._rule = rule

    @property
    def ab(self):
        if self._ab is None:
            print('Parameters not defined. Will now find them via minimization.')
            self.minimize_lsq()
        return self._ab

    @ab.setter
    def ab(self, ab):
        self._ab = ab

    @classmethod
    def from_hdf(cls, dt_start, dt_end, filenames=[H5_PATH], radar=False,
                 **kwargs):
        """Create Case object from a hdf file."""
        for dt in [dt_start, dt_end]:
            dt = pd.datetools.to_datetime(dt)
        pluvio200 = read.Pluvio(filenames, hdf_table='pluvio200')
        pluvio400 = read.Pluvio(filenames, hdf_table='pluvio400')
        dsd = read.PipDSD(filenames, hdf_table='pip_dsd')
        pipv = read.PipV(filenames, hdf_table='pip_vel')
        instr_lst = [pluvio200, pluvio400, dsd, pipv]
        if radar:
            xsacr = read.Radar(filenames, hdf_table='XSACR')
            kasacr = read.Radar(filenames, hdf_table='KASACR')
            kazr = read.Radar(filenames, hdf_table='KAZR')
            mwacr = read.Radar(filenames, hdf_table='MWACR')
            instr_lst = [pluvio200, pluvio400, dsd, pipv, xsacr, kasacr, kazr,
                         mwacr]
        for instr in instr_lst:
            instr.set_span(dt_start, dt_end)
        if radar:
            m200 = cls(dsd, pipv, pluvio200, xsacr, kasacr, kazr, mwacr,
                       **kwargs)
            m400 = cls(dsd, pipv, pluvio400, xsacr, kasacr, kazr, mwacr,
                       **kwargs)
        else:
            m200 = cls(dsd, pipv, pluvio200, **kwargs)
            m400 = cls(dsd, pipv, pluvio400, **kwargs)
        return m200, m400

    def dtstr(self, dtformat='%b %d'):
        """date string in simple format"""
        start, end = self.dt_start_end()
        dtstr = start.strftime(dtformat)
        if start.date() != end.date():
            dtstr += '-' + end.strftime(dtformat)
        return dtstr

    def between_datetime(self, dt_start, dt_end, inplace=False,
                         autoshift=False, autobias=False):
        """Select data only in chosen time frame."""
        dt_start = pd.datetools.to_datetime(dt_start)
        dt_end = pd.datetools.to_datetime(dt_end)
        if inplace:
            m = self
        else:
            m = copy.deepcopy(self)
        for instr in m.instr.values():
            instr.between_datetime(dt_start, dt_end, inplace=True)
            instr.case = m
        m.instr['pluvio'].bias = 0
        if autoshift:
            m.autoshift(inplace=True)
        if autobias:
            m.noprecip_bias(inplace=True)
        m.reset()
        return m

    def reset(self):
        """Reset memory cache."""
        if self.varinterval:
            self.rule = None

    def intensity(self, params=None, simple=False):
        """Calculate precipitation intensity using given or saved parameters.
        """
        if params is None and not self.liquid:
            params = self.ab
        if self.liquid:
            fits = self.series_nans()
            fits.loc[:] = read.gunn_kinzer
            fits.name = read.gunn_kinzer.name
            self.instr['pipv'].fits = pd.DataFrame(fits)
            r = self.sum_over_d(self.r_rho, rho=RHO_W)
        elif simple:
            r = self.sum_over_d(self.r_rho, rho=params[0])
        else:
            r = self.sum_over_d(self.r_ab, alpha=params[0], beta=params[1])
        if self.varinterval:
            return r
        return r.reindex(self.instr['pluvio'].amount(rule=self.rule).index).fillna(0)

    def amount(self, **kwargs):
        """Calculate precipitation in mm using given or saved parameters."""
        i = self.intensity(**kwargs)
        if self.varinterval:
            delta = self.instr['pluvio'].tdelta()
        else:
            delta = i.index.freq.delta
        return i*(delta/pd.datetools.timedelta(hours=1))

    def sum_over_d(self, func, **kwargs):
        """numerical integration over particle diameter"""
        dD = self.instr['dsd'].d_bin
        result = self.series_zeros()
        for d in self.instr['dsd'].bin_cen():#good_data().columns:
            result = result.add(func(d, **kwargs)*dD, fill_value=0)
        return result

    def r_ab(self, d, alpha, beta):
        """(mm/h)/(m/s)*kg/mg / kg/m**3 * mg/mm**beta * mm**beta * m/s * 1/(mm*m**3)
        """
        return 3.6/RHO_W*alpha*(PIP_CORR*d)**beta*self.v(d)*self.n(d)
        #dBin = self.instr['dsd'].d_bin
        #av = self.instr['pipv'].fit_params()['a']
        #bv = self.instr['pipv'].fit_params()['b']
        #return 3.6/RHO_W*alpha*self.n(d)*av/(dBin*(bv+beta+1))*((d+dBin*0.5)**(bv+beta+1)-(d-dBin*0.5)**(bv+beta+1))

    def r_rho(self, d, rho):
        """(mm/h)/(m/s)*m**3/mm**3 * kg/m**3 / (kg/m**3) * mm**3 * m/s * 1/(mm*m**3)
        """
        return 3.6e-3*TAU/12*rho/RHO_W*(PIP_CORR*d)**3*self.v(d)*self.n(d)
        #self.v(d)
        #dBin = self.instr['dsd'].d_bin
        #av = self.instr['pipv'].fit_params()['a']
        #bv = self.instr['pipv'].fit_params()['b']
        #return 3.6e-3*TAU/12*rho/RHO_W*self.n(d)*av/(dBin*(bv+4))*((d+dBin*0.5)**(bv+4)-(d-dBin*0.5)**(bv+4))

    def v(self, d):
        """velocity wrapper"""
        return self.instr['pipv'].v(d, varinterval=self.varinterval,
                                    rule=self.rule)

    def n(self, d):
        """N wrapper"""
        return self.instr['dsd'].n(d, varinterval=self.varinterval,
                                   rule=self.rule)

    def n_t(self):
        """total concentration"""
        name = 'N_t'
        def func():
            nt = self.sum_over_d(self.n)
            nt.name = name
            return nt
        return self.msger(name, func)

    def cache_dir(self):
        dt_start, dt_end = self.dt_start_end()
        return super().cache_dir(dt_start, dt_end, self.instr['pluvio'].name)

    def d_m(self):
        """mass weighted mean diameter"""
        name = 'D_m'
        def func():
            dm = self.n_moment(4)/self.n_moment(3)
            dm.name = name
            return dm
        return self.msger(name, func)

    def d_0(self):
        """median volume diameter"""
        name = 'D_0'
        def func():
            idxd = self.instr['dsd'].good_data().columns
            dd = pd.Series(idxd)
            dD = self.instr['dsd'].d_bin
            d3n = lambda d: (PIP_CORR*d)**3*self.n(d)*dD
            #d3n = lambda d: dD*self.n(d)*((d+dD*0.5)**4.0-(d-dD*0.5)**4.0)/(dD*4.0)
            cumvol = dd.apply(d3n).cumsum().T
            cumvol.columns = idxd
            sumvol = cumvol.iloc[:, -1]
            diff = cumvol-sumvol/2
            dmed = diff.abs().T.idxmin()
            dmed[sumvol < 0.0001] = 0
            dmed.name = name
            return dmed
        return self.msger(name, func)

    def d_max(self):
        """maximum diameter from PSD tables"""
        name = 'D_max'
        def func():
            idxd = self.instr['dsd'].good_data().columns
            dd = pd.Series(idxd)
            nd = dd.apply(self.n).T
            nd.columns = idxd
            dmax = (PIP_CORR*nd[nd > 0.0001].T.apply(pd.Series.last_valid_index).fillna(0))
            dmax.name = name
            return dmax
        return self.msger(name, func)

    def n_moment(self, n):
        moment = lambda d: (PIP_CORR*d)**n*self.n(d)
        #dD = self.instr['dsd'].d_bin
        #moment = lambda d: self.n(d)*((d+dD*0.5)**(n+1.0)-(d-dD*0.5)**(n+1.0))/(dD*(n+1.0))
        nth_mo = self.sum_over_d(moment)
        nth_mo.name = 'M' + str(n)
        return nth_mo

    def eta(self):
        eta = self.n_moment(4)**2/(self.n_moment(6)*self.n_moment(2))
        eta.name = 'eta'
        return eta

    def mu(self):
        eta = self.eta()
        mu = ((7-11*eta)-np.sqrt(eta**2+14*eta+1))/(2*(eta-1))
        mu.name = 'mu'
        return mu

    def lam(self):
        mu = self.mu()
        lam = np.sqrt(self.n_moment(2)*gamma(mu+5)/(self.n_moment(4)*gamma(mu+3)))
        lam.name = 'lambda'
        return lam

    def n_0(self):
        mu = self.mu()
        n0 = self.n_moment(2)*self.lam()**(mu+3)/gamma(mu+3)
        n0.name = 'N_0'
        return n0

    def d_0_gamma(self):
        name = 'D_0_gamma'
        def func():
            d0 = (3.67+self.mu())/self.lam()
            d0.name = name
            return d0
        return self.msger(name, func)

    def partcount(self):
        """particle count"""
        count = self.instr['pipv'].partcount(rule=self.rule,
                                             varinterval=self.varinterval)
        count.name = 'count'
        return count

    def series_zeros(self):
        """Return series of zeros of the shape of timestep averaged data."""
        return self.instr['pluvio'].acc(rule=self.rule)*0.0

    def series_nans(self):
        """Return series of nans of the shape of timestep averaged data."""
        return self.series_zeros()*np.nan

    def noprecip_bias(self, inplace=True):
        """Wrapper to unbias pluvio using LWC calculated from PIP data."""
        return self.instr['pluvio'].noprecip_bias(self.instr['pipv'].lwc(),
                                                  inplace=inplace)

    def pluvargs(self):
        args = {}
        if not self.varinterval:
            args['rule'] = self.rule
        return args

    def cost(self, c, use_accum=True):
        """Cost function for minimization"""
        if use_accum:
            pip_precip = self.acc(params=c)
            cost_method = self.instr['pluvio'].acc
        else:
            pip_precip = self.intesity(params=c)
            cost_method = self.instr['pluvio'].intensity()
        return abs(pip_precip.add(-1*cost_method(**self.pluvargs())).sum())

    def cost_lsq(self, beta):
        """Single variable cost function using lstsq to find linear coef."""
        alpha = self.alpha_lsq(beta)
        return self.cost([alpha, beta])

    def const_lsq(self, c, simple):
        acc_arr = self.acc(params=c, simple=simple).values
        A = np.vstack([acc_arr, np.ones(len(acc_arr))]).T
        y = self.instr['pluvio'].acc(**self.pluvargs()).values
        return np.linalg.lstsq(A, y)[0][0]

    def alpha_lsq(self, beta):
        """Wrapper for const_lsq to calculate alpha"""
        return self.const_lsq(c=[1, beta], simple=False)

    def density_lsq(self):
        """Wrapper for const_lsq to calculate least square particle density"""
        return self.const_lsq(c=[1], simple=True)

    def density(self, pluvio_filter=True, pip_filter=False, rhomax=None):
        """Calculates mean density estimate for each timeframe."""
        name = 'density'
        def func():
            rho_r_pip = self.amount(params=[1], simple=True)
            if pluvio_filter:
                rho_r_pip[self.instr['pluvio'].intensity() < 0.1] = np.nan
            if pip_filter and self.ab is not None:
                rho_r_pip[self.intensity() < 0.1] = np.nan
            rho = self.instr['pluvio'].amount(rule=self.rule)/rho_r_pip
            rho.name = name
            if rhomax is not None:
                rho[rho > rhomax] = np.nan
            return rho.replace(np.inf, np.nan)
        return self.msger(name, func)

    def data_in_density_range(self, data, rhomin, rhomax):
        grouped = read.merge_series(data, self.instr['pluvio'].grouper())
        outdata = pd.merge(grouped, pd.DataFrame(self.density()),
                           left_on='group', right_index=True)
        return outdata.query('%s < density < %s' % (rhomin, rhomax))

    def vfit_density_range(self, lims, **fitargs):
        rhomin = lims[0]
        rhomax = lims[1]
        data = self.data_in_density_range(self.instr['pipv'].good_data(),
                                          rhomin, rhomax)
        fit = self.instr['pipv'].find_fit(data=data, **fitargs)[0]
        fit.x_unfiltered = data.Wad_Dia.values
        fit.y_unfiltered = data.vel_v.values
        return fit

    def vfits_density_range(self, limslist, parallel=True, **fitargs):
        if parallel:
            processes = min(len(limslist), os.cpu_count())
            with Pool(processes) as p:
                fits = p.map(self.vfit_density_range, limslist)
            p.join()
            return fits
        fits = []
        for lims in limslist:
            fits.append(self.vfit_density_range(lims, **fitargs))
        return fits

    def plot_vfits_in_density_ranges(self, rholimits=(0, 150, 300, 800),
                                     separate=False, hide_high_limit=True,
                                     fitargs={}, parallel=True, **kwargs):
        limslist = [(rhomin, rholimits[i+1]) for i, rhomin in enumerate(rholimits[:-1])]
        dlabel = 'Equivalent diameter (mm)'
        vlabel = 'Fall velocity (ms$^{-1}$)'
        fits = self.vfits_density_range(limslist, parallel=parallel, **fitargs)
        n_ranges = len(fits)
        if separate:
            fig, axarr = plt.subplots(1, n_ranges, sharex=True,
                                      sharey=True, dpi=100, tight_layout=True,
                                      figsize=(n_ranges*6, 6))
        else:
            fig, ax = plt.subplots(tight_layout=True)
        for i, fit in enumerate(fits):
            if separate:
                ax = axarr[i]
            lims = limslist[i]
            rhomin = lims[0]
            rhomax = lims[1]
            limitstr = '$%s < \\rho < %s$' % (rhomin, rhomax)
            fitstr = '$' + str(fit) + '$'
            fit.plot(ax=ax, label=fitstr, **kwargs)
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels, loc='lower right')
            ax.set_xlabel(dlabel)
            ax.set_title(limitstr)
        if hide_high_limit:
            limitstr = '$\\rho > ' + str(rhomin) + '$'
            ax.set_title(limitstr)
        if separate:
            axarr[0].set_ylabel(vlabel)
        else:
            ax.set_ylabel(vlabel)
            ax.set_title(self.dtstr())
        if separate:
            return axarr
        return ax

    def z(self, radarname='XSACR'):
        """Radar reflectivity wrapper"""
        if radarname == 'XSACR':
            return self.instr['xsacr'].z(varinterval=self.varinterval,
                                         rule=self.rule)
        elif radarname == 'KASACR':
            return self.instr['kasacr'].z(varinterval=self.varinterval,
                                          rule=self.rule)
        elif radarname == 'KAZR':
            return self.instr['kazr'].z(varinterval=self.varinterval,
                                        rule=self.rule)
        elif radarname == 'MWACR':
            return self.instr['mwacr'].z(varinterval=self.varinterval,
                                         rule=self.rule)
        # TODO: else throw error "unknown radarname"

    def Z_rayleigh_Xband(self, pluvio_filter=True, pip_filter=False,
                         density=None):
        """Use rayleigh formula and maxwell-garnett EMA to compute radar
        reflectivity Z"""
        name = "reflXray"
        constant = 0.2/(0.93*917*917)
        if density is None:
            density = self.density(pluvio_filter=pluvio_filter,
                                   pip_filter=pip_filter)
        Z = 10.0*np.log10(constant*density*density*self.n_moment(6))
        Z.name = name
        return Z

    def volume_avg_density(self, density_size, pluvio_filter=True,
                           pip_filter=False):
        """Calculate volume averaged bulk density for the given density size realation"""
        name = "rho3"
        def density_func(d):
            return density_size((PIP_CORR*d))*self.n(d)    # I am experimenting with precise integration leaved d**3
        def mom3(d):
            return ((PIP_CORR*(d+0.5*self.instr['dsd'].d_bin))**4-(PIP_CORR*(d-0.5*self.instr['dsd'].d_bin))**4)*self.n(d)/(self.instr['dsd'].d_bin*4)
        density = self.sum_over_d(density_func)/self.sum_over_d(mom3)
        density[density.isnull()] = 0
        density.name = name
        return density

    def reflectivity_avg_density(self, density_size, pluvio_filter=True,
                                 pip_filter=False):
        """Calculate volume averaged bulk density for the given density size
        realation"""
        name = "rho6"
        def density_func(d):
            # I am experimenting with precise integration leaved d**6 and squares
            return density_size((PIP_CORR*d))*self.n(d)
        def mom6(d):
            return ((PIP_CORR*(d+0.5*self.instr['dsd'].d_bin))**7-(PIP_CORR*(d-0.5*self.instr['dsd'].d_bin))**7)*self.n(d)/(self.instr['dsd'].d_bin*7)
        density = self.sum_over_d(density_func)/self.sum_over_d(mom6)
        density.name = name
        density[density.isnull()] = 0
        return np.sqrt(density)

    def tmatrix(self, wl, pluvio_filter=True, pip_filter=False, density=None):
        """Calculate radar reflectivity at requested wavelength wl [mm] using
        T-matrix"""
        name = switch_wl(wl) + "reflTM"
        if density is None:
            density = self.density(pluvio_filter=pluvio_filter,
                                   pip_filter=pip_filter)
        Zserie = pd.Series(density)
        dBin = self.instr['dsd'].d_bin
        edges = PIP_CORR*self.instr['dsd'].data.columns.values+0.5*dBin
        grp = self.instr['dsd'].grouped(rule=self.rule,
                                        varinterval=self.varinterval,
                                        col=self.dsd.bin_cen())
        PSDvalues = grp.mean()
        for item in density.iteritems():
            if item[1] > 0.0:
                ref = refractive.mi(wl, 0.001*item[1])
                print(ref)
                flake = tmatrix.Scatterer(wavelength=wl, m=ref,
                                          axis_ratio=1.0/1.0)
                flake.psd_integrator = psd.PSDIntegrator()
                flake.psd_integrator.D_max = 35.0
                flake.psd = psd.BinnedPSD(bin_edges=edges,
                                          bin_psd=PSDvalues.loc[item[0]].values)
                flake.psd_integrator.init_scatter_table(flake)
                Z = 10.0*np.log10(radar.refl(flake))
            else:
                Z = np.nan
            Zserie.loc[item[0]] = Z
        Zserie.name = name
        return Zserie

    def minimize(self, method='SLSQP', **kwargs):
        """Legacy method for determining alpha and beta."""
        print('Optimizing parameters...')
        result = minimize(self.cost, self.quess, method=method, **kwargs)
        self.ab = result.x
        return result

    def minimize_lsq(self):
        """Find beta by minimization and alpha by linear least square."""
        print('Optimizing parameters...')
        result = minimize(self.cost_lsq, self.quess[1], method='Nelder-Mead')
        #self.result = minimize(self.cost_lsq, self.quess[1], method='SLSQP', bounds=self.bnd[1])
        print(result.message)
        beta = result.x[0]
        alpha = self.alpha_lsq(beta)
        self.ab = [alpha, beta]
        return result

    def dt_start_end(self):
        """case start and end time as Timestamp"""
        t = self.time_range()
        return (t[0], t[-1])

    def time_range(self):
        """data time ticks on minute interval"""
        return pd.date_range(self.instr['pluvio'].acc().index[0],
                             self.instr['pluvio'].acc().index[-1], freq='1min')

    def plot(self, axarr=None, kind='line', label_suffix='', pip=True,
             **kwargs):
        """Plot calculated (PIP) and pluvio intensities."""
        if axarr is None:
            f, axarr = plt.subplots(4, sharex=True, dpi=120)
        if pip:
            self.intensity().plot(label='PIP ' + label_suffix, kind=kind,
                                  ax=axarr[0], **kwargs)
        self.instr['pluvio'].intensity(rule=self.rule).plot(label=self.instr['pluvio'].name + ' ' + label_suffix,
                                                   kind=kind, ax=axarr[0],
                                                   **kwargs)
        axarr[0].set_ylabel('mm/h')
        if self.liquid:
            title = 'rain intensity'
        elif not pip:
            title = 'precipitation intensity'
        else:
            title = r'precipitation intensity, $\alpha=%s, \beta=%s$' % (self.ab[0], self.ab[1])
        axarr[0].set_title(title)
        rho = self.density()
        rho.plot(label=label_suffix, ax=axarr[1], **kwargs)
        axarr[1].set_ylabel(r'$\rho_{b}$')
        self.n_t().plot(label=label_suffix, ax=axarr[2], **kwargs)
        axarr[2].set_ylabel(r'$N_{tot} (m^{-3})$')
        self.d_m().plot(label=label_suffix, ax=axarr[3], **kwargs)
        axarr[3].set_ylabel(r'$D_m$ (mm)')
        for ax in axarr:
            ax.legend(loc='upper right')
        for i in [0, 1, 2]:
            axarr[i].set_xlabel('')
        axarr[-1].set_xlabel('time (UTC)')
        plt.show()
        return axarr

    def plot_cost(self, resolution=20, ax=None, cmap='binary', **kwargs):
        """The slowest plot you've drawn"""
        if ax is None:
            ax = plt.gca()
        alpha0 = self.ab[0]
        alpha = np.linspace(0.4*alpha0, 1.4*alpha0, num=resolution)
        beta = np.linspace(self.bnd[1][0], self.bnd[1][1], num=resolution)
        z = np.zeros((alpha.size, beta.size))
        for i, a in enumerate(alpha):
            for j, b in enumerate(beta):
                z[i][j] = self.cost((a, b))
        ax = plt.gca()
        heat = ax.pcolor(beta, alpha, z, cmap=cmap, **kwargs)
        ax.colorbar()
        ax.set_xlabel(r'$\beta$')
        ax.set_ylabel(r'$\alpha$')
        ax.axis('tight')
        ax.set_title('cost function value')
        return z, heat, ax.plot(self.ab[1], self.ab[0], 'ro')

    def plot_cost_lsq(self, resolution, ax=None, *args, **kwargs):
        """Plot cost function value vs. beta."""
        if ax is None:
            ax = plt.gca()
        beta = np.linspace(self.bnd[1][0], self.bnd[1][1], num=resolution)
        cost = np.array([self.cost_lsq(b) for b in beta])
        ax = plt.gca()
        ax.set_xlabel(r'$\beta$')
        ax.set_ylabel('cost')
        ax.set_title('cost function value')
        return ax.plot(beta, cost, *args, **kwargs)

    def plot_velfitcoefs(self, fig=None, ax=None, rhomin=None, rhomax=None,
                         countmin=1, **kwargs):
        rho = self.density()
        params = self.instr['pipv'].fits.polfit.apply(lambda fit: fit.params)
        selection = pd.DataFrame([rho.notnull(),
                                  self.partcount() > countmin]).all()
        rho = rho[selection]
        params = params[selection]
        a = params.apply(lambda p: p[0]).values
        b = params.apply(lambda p: p[1]).values
        if fig is None:
            fig = plt.figure(dpi=120)
        if ax is None:
            ax = plt.gca()
        choppa = ax.scatter(a, b, c=rho.values, vmin=rhomin, vmax=rhomax,
                            **kwargs)
        fig.colorbar(choppa, label='bulk density')
        ax.set_xlabel('$a_u$', fontsize=15)
        ax.set_ylabel('$b_u$', fontsize=15)
        return ax

    def plot_d0_bv(self, rhomin=None, rhomax=None, countmin=1,
                   count_as_size=True, countscale=4, **kwargs):
        rho = self.density()
        params = self.instr['pipv'].fits.polfit.apply(lambda fit: fit.params)
        selection = pd.DataFrame([rho.notnull(),
                                  self.partcount() > countmin]).all()
        count = self.partcount()[selection]
        rho = rho[selection]
        params = params[selection]
        d0 = self.d_0_gamma()[selection]
        b = params.apply(lambda p: p[1])
        b.name = 'b'
        if count_as_size:
            kwargs['s'] = 0.01*countscale*count
        if rhomin is None:
            rhomin = rho.min()
        if rhomax is None:
            rhomax = rho.max()
        return scatterplot(x=d0, y=b, c=rho, vmin=rhomin, vmax=rhomax,
                           **kwargs)

    def summary(self, **kwargs):
        """Return a DataFrame of combined numerical results."""
        casename = self.series_nans().fillna(self.dtstr(**kwargs))
        casename.name = 'case'
        data = read.merge_multiseries(self.partcount(), self.density(),
                                      self.d_0(), self.n_t(), casename,
                                      self.instr['pipv'].fit_params(),
                                      self.d_m(), self.d_max(),
                                      self.d_0_gamma(),
                                      self.amount(params=[100], simple=True),
                                      self.instr['pluvio'].amount(rule=self.rule),
                                      self.eta(), self.mu(), self.lam(),
                                      self.n_0(),
                                      self.n_moment(0), self.n_moment(1),
                                      self.n_moment(2),
                                      self.Z_rayleigh_Xband(),
                                      self.tmatrix(tm_aux.wl_X))
        data.index.name = 'datetime'
        return data.sort_index(axis=1)

    def xcorr(self, rule='1min', ax=None, **kwargs):
        """Plot cross-correlation between lwc estimate and pluvio intensity.
        Extra arguments are passed to pyplot.xcorr.
        """
        if ax is None:
            ax = plt.gca()
        r = self.instr['pluvio'].intensity(rule=rule, unbias=False)
        lwc = self.instr['pipv'].lwc(rule).reindex(r.index).fillna(0)
        return ax.xcorr(lwc, r, **kwargs)

    def autoshift(self, rule='1min', inplace=False):
        """Find and correct pluvio time shift using cross correlation."""
        if self.instr['pluvio'].shift_periods != 0:
            print('Pluvio already timeshifted, resetting.')
            self.instr['pluvio'].shift_reset()
        xc = self.xcorr(rule=rule)
        imaxcorr = xc[1].argmax()
        periods = xc[0][imaxcorr]
        if inplace:
            self.instr['pluvio'].shift_periods = periods
            self.instr['pluvio'].shift_freq = rule
            print('Pluvio timeshift set to %s*%s.'
                  % (str(self.instr['pluvio'].shift_periods),
                     self.instr['pluvio'].shift_freq))
        return periods

    def clear_cache(self):
        xtra = glob(os.path.join(self.cache_dir(), '*' + read.MSGTLD))
        xtra.extend(glob(os.path.join(self.cache_dir(), '*.h5')))
        super().clear_cache(extra_files=xtra)
        self.reset()


class Snow2:
    """UNTESTED.
    Calculate snowfall rate using Szyrmer Zawadski's method from Snow Study II.
    """
    def __init__(self):
        return

    @staticmethod
    def best(re, mh=True):
        if mh:  # MH05
            cl = np.array([3.8233, -1.5211, 0.30065, -0.06104, 0.13074,
                           -0.073429, 0.016006, -0.0012483])
        else:   # KC05
            cl = np.array([3.8816, -1.4579, 0.27749, -0.41521, 0.57683,
                           -0.29220, 0.06467, -0.0053405])
        logx = 0
        for l, c in enumerate(cl):
            logx += c*np.log(re)**l
        return np.exp(logx)

    @staticmethod
    def mass(u, ar, d):
        g = 9.81
        fa = 1
        rho_a = 1.275
        nu_a = 1.544e-5
        re = u*d/nu_a
        return np.pi*rho_a*nu_a**2/(8*g)*Snow2.best(re)*ar*fa
