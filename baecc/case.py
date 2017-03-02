# coding: utf-8
import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
from os import path
import baecc
from baecc import caching, instruments
from scipy.optimize import minimize
from scipy.special import gamma
from pytmatrix import tmatrix, psd, refractive, radar
from pytmatrix import tmatrix_aux as tm_aux

RHO_W = 1000


def scatterplot(x, y, c=None, kind='scatter', **kwargs):
    """scatter plot of two Series objects"""
    plotdata = baecc.merge_series(x, y)
    if c is not None:
        kwargs['c'] = c
    return plotdata.plot(kind=kind, x=x.name, y=y.name, **kwargs)


def split_index(df, date=pd.datetime(2014,7,1), names=('first', 'second')):
    isfirst = df.index < date
    idf = pd.Series(isfirst, index=df.index)
    idf[isfirst] = names[0]
    idf[-isfirst] = names[1]
    tuples = list(zip(*(idf.values, idf.index.values)))
    index = pd.MultiIndex.from_tuples(tuples, names=('winter', 'datetime'))
    df.index = index
    return df


def daterange2str(start, end, dtformat='{day}{month}{year}', delimiter='-',
          hour_fmt='%H', day_fmt='%d.', month_fmt='%m.', year_fmt='%Y'):
    """date range in simple human readable format"""
    same_date = start.date() == end.date()
    start_fmt = dtformat
    end_fmt = dtformat
    for attr in ('minute', 'hour', 'day', 'month', 'year'):
        if getattr(start, attr) == getattr(end, attr) and not same_date:
            start_fmt = start_fmt.replace(('{%s}' % attr),'')
    start_fmt = start_fmt.strip()
    formats = {'hour':hour_fmt, 'day':day_fmt, 'month':month_fmt,
               'year':year_fmt}
    start_str = start.strftime(start_fmt.format(**formats))
    end_str = end.strftime(end_fmt.format(**formats))
    if same_date:
        return start_str
    return start_str + delimiter + end_str


def switch_wl(x):
    return {tm_aux.wl_C: "C", tm_aux.wl_X: "X", tm_aux.wl_Ku: "Ku",
            tm_aux.wl_Ka: "Ka", tm_aux.wl_W: "W"}.get(x, str(x))


def limitslist(limits):
    return [(mini, limits[i+1]) for i, mini in enumerate(limits[:-1])]


def plot_vfits_rho_intervals(fits, limslist, separate=False,
                             hide_high_limit=True, dpi=180, **kwargs):
    dlabel = '$D$, mm'
    vlabel = '$v$, m$\,$s$^{-1}$'
    n_ranges = len(fits)
    if separate:
        fig, axarr = plt.subplots(1, n_ranges, sharex=True,
                                  sharey=True, dpi=dpi, tight_layout=True,
                                  figsize=(n_ranges*3, 3))
    else:
        fig, ax = plt.subplots(tight_layout=True)
    for i, fit in enumerate(fits):
        if separate:
            ax = axarr[i]
        lims = limslist[i]
        rhomin = lims[0]
        rhomax = lims[1]
        limitstr = '$%s < \\rho \leq %s$' % (rhomin*baecc.RHO_SCALE,
                                             rhomax*baecc.RHO_SCALE)
        fitstr = '$' + str(fit) + '$'
        fit.plot(ax=ax, label=fitstr, **kwargs)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='upper right')
        #ax.set_xlabel(dlabel)
        ax.set_title(limitstr)
    middle_ax = axarr[len(axarr)/2]
    middle_ax.set_xlabel(dlabel)
    if hide_high_limit:
        limitstr = '$\\rho > %s$' % (rhomin*baecc.RHO_SCALE)
        ax.set_title(limitstr)
    if separate:
        axarr[0].set_ylabel(vlabel)
    else:
        ax.set_ylabel(vlabel)
    if separate:
        return fig, axarr
    return fig, ax


class Case(instruments.PrecipMeasurer, caching.Cacher):
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
        if autoshift:
            self.autoshift()
        if unbias:
            self.noprecip_bias()
        caching.Cacher.__init__(self)
        for instr in self.instr.values():
            instr.parent = self

    def __repr__(self):
        start, end = self.dt_start_end()
        return '%s case from %s to %s, %s' % (self.casetype(), start, end,
                                                  self.intervalstr())

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
    def from_hdf(cls, dt_start, dt_end, filenames=[baecc.H5_PATH], radar=False,
                 pluvio_name='pluvio200', **kwargs):
        """Create Case object from a hdf file."""
        for dt in [dt_start, dt_end]:
            dt = pd.datetools.to_datetime(dt)
        pluvio = instruments.Pluvio(filenames, hdf_table=pluvio_name)
        dsd = instruments.PipPSD(filenames, hdf_table='pip_dsd')
        pipv = instruments.PipV(filenames, hdf_table='pip_vel')
        if radar:
            xsacr = instruments.Radar(filenames, hdf_table='XSACR')
            kasacr = instruments.Radar(filenames, hdf_table='KASACR')
            kazr = instruments.Radar(filenames, hdf_table='KAZR')
            mwacr = instruments.Radar(filenames, hdf_table='MWACR')
            instr_lst = [dsd, pipv, pluvio, xsacr, kasacr, kazr,
                         mwacr]
        else:
            instr_lst = [dsd, pipv, pluvio]
        for instr in instr_lst:
            instr.set_span(dt_start, dt_end)
        return cls(*instr_lst, **kwargs)

    def casetype(self):
        if self.liquid:
            return 'rain'
        return 'snow'

    def intervalstr(self):
        if self.varinterval:
            return 'adaptive'
        return self.rule

    def fingerprint(self):
        idstr = str(self.dt_start_end()) + self.casetype() + self.intervalstr()
        for key, instr in sorted(self.instr.items()):
            idstr += instr.fingerprint()
        return caching.fingerprint(idstr)

    def dtstr(self, dtformat='{day}{month}{year}', **kws):
        """date string in simple human readable format"""
        return daterange2str(*self.dt_start_end(), dtformat=dtformat, **kws)

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
            instr.case = m # TODO: get rid of this
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
            fits.loc[:] = baecc.fit.gunn_kinzer
            fits.name = baecc.fit.gunn_kinzer.name
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
        dD = self.instr['dsd'].bin_width()
        result = self.series_zeros()
        for d in self.instr['dsd'].bin_cen():
            result = result.add(func(d, **kwargs)*dD[d], fill_value=0)
        return result

    def r_ab(self, d, alpha, beta):
        """(mm/h)/(m/s)*kg/mg / kg/m**3 * mg/mm**beta * mm**beta * m/s * 1/(mm*m**3)
        """
        return 3.6/RHO_W*alpha*d**beta*self.v(d)*self.n(d)
        #dBin = self.instr['dsd'].d_bin
        #av = self.instr['pipv'].fit_params()['a']
        #bv = self.instr['pipv'].fit_params()['b']
        #return 3.6/RHO_W*alpha*self.n(d)*av/(dBin*(bv+beta+1))*((d+dBin*0.5)**(bv+beta+1)-(d-dBin*0.5)**(bv+beta+1))

    def r_rho(self, d, rho):
        """(mm/h)/(m/s)*m**3/mm**3 * kg/m**3 / (kg/m**3) * mm**3 * m/s * 1/(mm*m**3)
        """
        return 3.6e-3*np.pi/6*rho/RHO_W*d**3*self.n(d)*self.v(d)
        #self.v(d)
        #dBin = self.instr['dsd'].d_bin
        #av = self.instr['pipv'].fit_params()['a']
        #bv = self.instr['pipv'].fit_params()['b']
        #return 3.6e-3*TAU/12*rho/RHO_W*self.n(d)*av/(dBin*(bv+4))*((d+dBin*0.5)**(bv+4)-(d-dBin*0.5)**(bv+4))

    def w_slice(self, d, **kwargs):
        rho = self.density(**kwargs)
        return 1e-6*np.pi/6*rho*d**3*self.n(d)

    def w(self, method='gamma', **kwargs):
        """water content in g/m**3"""
        if method == 'integral':
            return self.sum_over_d(self.w_slice)
        if method == 'gamma':
            rho = self.density(**kwargs)
            mu = self.mu()
            return 1e-3*np.pi/6*rho*self.n_0()*gamma(mu+4)*self.d_0()**(mu+4)/((3.67+mu)**(mu+4))
        return

    def intervalled(self, func, *args, **kwargs):
        return func(*args, varinterval=self.varinterval,
                    rule=self.rule, **kwargs)

    def v(self, d):
        """velocity wrapper, m/s"""
        return self.intervalled(self.instr['pipv'].v, d)

    def n(self, d):
        """N wrapper"""
        return self.intervalled(self.instr['dsd'].n, d)

    def f_mu(self):
        mu = self.mu()
        return 6/(3.67)**4*(3.67+mu)**(mu+4)/(gamma(mu+4))

    def n_t(self):
        """total concentration"""
        name = 'N_t'
        def func():
            nt = self.sum_over_d(self.n)
            nt.name = name
            return nt
        return self.msger(name, func)

    def d_m(self):
        """mass weighted mean diameter, mm"""
        name = 'D_m'
        def func():
            dm = self.n_moment(4)/self.n_moment(3)
            dm.name = name
            return dm
        return self.msger(name, func)

    def d_0(self):
        """median volume diameter, mm"""
        name = 'D_0'
        def func():
            idxd = self.instr['dsd'].good_data().columns
            dd = pd.Series(idxd)
            dD = self.instr['dsd'].bin_width()
            d3n = lambda d: d**3*self.n(d)*dD[d]
            #d3n = lambda d: dD[d]*self.n(d)*((d+dD[d]*0.5)**4.0-(d-dD[d]*0.5)**4.0)/(dD[d]*4.0)
            cumvol = dd.apply(d3n).cumsum().T
            cumvol.columns = idxd
            sumvol = cumvol.iloc[:, -1]
            diff = cumvol.sub(sumvol/2, axis=0)
            dmed = diff.abs().T.idxmin()
            dmed[sumvol < 0.0001] = 0
            dmed.name = name
            return dmed
        return self.msger(name, func)

    def d_max(self):
        """maximum diameter from PSD tables, mm"""
        name = 'D_max'
        def func():
            idxd = self.instr['dsd'].good_data().columns
            dd = pd.Series(idxd)
            nd = dd.apply(self.n).T
            nd.columns = idxd
            dmax = nd[nd > 0.0001].T.apply(pd.Series.last_valid_index).fillna(0)
            dmax.name = name
            return dmax
        return self.msger(name, func)

    def n_moment(self, n):
        name = 'M' + str(n)
        def func():
            moment = lambda d: d**n*self.n(d)
            nth_mo = self.sum_over_d(moment)
            nth_mo.name = name
            return nth_mo
        return self.msger(name, func)
    
    #TODO: What is the difference between this and n_moment?
    def mom_n(self, n):
        name = 'mom' + str(n)
        def func():
            dD = self.instr['dsd'].bin_width()
            d_high = lambda d: d+dD[d]*0.5
            d_low =  lambda d: d-dD[d]*0.5
            moment = lambda d: self.n(d)*(d_high(d)**(n+1)-d_low(d)**(n+1))/(dD[d]*(n+1))
            nth_mo = self.sum_over_d(moment)
            nth_mo.name = name
            return nth_mo
        return self.msger(name, func)

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

    def n_w(self):
        name = 'N_w'
        def func():
            integrand = lambda d: d**3*self.n(d)
            nw = 3.67**4/(6*self.d_0()**4)*self.sum_over_d(integrand)
            nw.name = name
            return nw
        return self.msger(name, func)

    def n_w_mu(self, **kwargs):
        mu = self.mu()
        # TODO: check scale
        return self.n_0()/self.f_mu()*self.d_0()**mu

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

    def group_by_density(self, data, rholimits):
        """Add columns rhomin and rhomax."""
        limslist = limitslist(rholimits)
        datalist = []
        for lims in limslist:
            datalist.append(self.data_in_density_range(data, lims[0], lims[1],
                                                       append_limits=True))
        out = pd.concat(datalist)
        return out.sort()

    def group(self, data, merger, drop_grouper=True):
        """Data should have same or higher frequency than merger."""
        grouped = baecc.merge_series(data, self.instr['pluvio'].grouper())
        result = pd.merge(grouped, pd.DataFrame(merger),
                           left_on='group', right_index=True)
        if drop_grouper:
            return result.drop('group', axis=1)
        return result

    def data_in_density_range(self, data, rhomin, rhomax, drop_grouper=True,
                              append_limits=False, **rhokws):
        """Return only data from timesteps where rhomin<rho<rhomax."""
        # TODO: write tests
        rho = self.density(**rhokws)
        outdata = self.group(data, rho, drop_grouper=drop_grouper)
        result = outdata.query('%s < %s < %s' % (rhomin, rho.name, rhomax))
        if append_limits:
            result['rhomin'] = rhomin
            result['rhomax'] = rhomax
        if drop_grouper:
            return result.drop(rho.name, axis=1)
        return result

    def vfit_density_range(self, rhomin, rhomax, data=None, **fitargs):
        pipv = self.instr['pipv']
        if data is None:
            data = pipv.good_data()
        data_in_range = self.data_in_density_range(data, rhomin, rhomax)
        fit = pipv.find_fit(data=data_in_range, **fitargs)[0]
        fit.x_unfiltered = data_in_range[pipv.d_col].values
        fit.y_unfiltered = data_in_range.vel_v.values
        return fit

    def vfits_density_range(self, limslist, **fitargs):
        params_id = str(limslist) + caching.hash_dict(fitargs)
        name = 'vfits_density_range' + caching.fingerprint(params_id)
        def func():
            fits = []
            for lims in limslist:
                fits.append(self.vfit_density_range(*lims, **fitargs))
            return fits
        return self.pickler(name, func)

    def plot_vfits_in_density_ranges(self, rholimits=(0, 100, 200, 1000),
                                     fitargs={}, **kwargs):
        limslist = limitslist(rholimits)
        fits = self.vfits_density_range(limslist, **fitargs)
        return plot_vfits_rho_intervals(fits, limslist, **kwargs)

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
            return density_size(d)*self.n(d)    # I am experimenting with precise integration leaved d**3
        density = self.sum_over_d(density_func)/self.mom_n(3)
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
            return density_size(d)*self.n(d)
        density = self.sum_over_d(density_func)/self.mom_n(6)
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
        z_serie = pd.Series(density)
        dBin = self.instr['dsd'].bin_width()
        edges = self.instr['dsd'].bin_cen()+0.5*dBin
        grp = self.instr['dsd'].grouped(rule=self.rule,
                                        varinterval=self.varinterval,
                                        col=self.dsd.bin_cen())
        psd_values = grp.mean()
        for item in density.iteritems():
            if item[1] > 0.0:
                ref = refractive.mi(wl, 0.001*item[1])
                print(ref)
                flake = tmatrix.Scatterer(wavelength=wl, m=ref,
                                          axis_ratio=1.0/1.0)
                flake.psd_integrator = psd.PSDIntegrator()
                flake.psd_integrator.D_max = 35.0
                flake.psd = psd.BinnedPSD(bin_edges=edges,
                                          bin_psd=psd_values.loc[item[0]].values)
                flake.psd_integrator.init_scatter_table(flake)
                Z = 10.0*np.log10(radar.refl(flake))
            else:
                Z = np.nan
            z_serie.loc[item[0]] = Z
        z_serie.name = name
        return z_serie

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
        if t.size < 1:
            placeholder = self.instr['pluvio'].good_data().index[0]
            return (placeholder, placeholder)
        return (t[0], t[-1])

    def duration(self):
        tstart, tend = self.dt_start_end()
        return tend-tstart

    def time_range(self):
        """data time ticks on minute interval"""
        dt_index = self.instr['pluvio'].acc().index
        if dt_index.size < 1:
            return pd.DatetimeIndex([], freq='T')
        return pd.date_range(dt_index[0], dt_index[-1], freq='1min')

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

    def summary(self, radar=False, include_vfits=False, split_date=None,
                **kwargs):
        """Return a DataFrame of combined numerical results."""
        casename = self.series_nans().fillna(self.dtstr(**kwargs))
        casename.name = 'case'
        pluvio = self.instr['pluvio']
        pipv = self.instr['pipv']
        params = [self.partcount(),
                  self.density(),
                  self.d_0(),
                  self.n_t(),
                  casename,
                  pipv.fit_params(),
                  #self.d(),
                  self.d_m(),
                  self.d_max(),
                  self.d_0_gamma(),
                  #self.amount(params=[100], simple=True), # What is this?
                  pluvio.amount(),
                  pluvio.intensity(),
                  pluvio.start_time(),
                  pluvio.half_time(),
                  pluvio.tdelta(),
                  self.eta(),
                  self.mu(),
                  self.lam(),
                  self.n_0(),
                  self.n_w(),
                  self.n_moment(0),
                  self.n_moment(1),
                  self.n_moment(2)]
        if radar:
            params.extend([self.Z_rayleigh_Xband(), self.tmatrix(tm_aux.wl_X)])
        if include_vfits:
            params.extend(pipv.fits[pipv.default_fit.name])
        data = baecc.merge_multiseries(*params)
        data.index.name = 'datetime'
        if split_date is not None:
            data = split_index(data, date=split_date)
        return data#.sort_index(axis=1) # TODO int col names?

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
        xtra = glob(path.join(self.cache_dir(), '*' + caching.MSGTLD))
        xtra.extend(glob(path.join(self.cache_dir(), '*.h5')))
        super().clear_cache(extra_files=xtra)
        self.reset()
