"""Tools for estimating density and other properties of falling snow"""
import numpy as np
import pandas as pd
import read
from datetime import datetime
from scipy.optimize import minimize
from scipy.special import gamma
from glob import glob
import matplotlib.pyplot as plt
import copy
import locale
import os

locale.setlocale(locale.LC_ALL, 'en_GB.UTF-8')

TAU = 2*np.pi
RHO_W = 1000

def batch_import(dtstr, datadir='../DATA'):
    """Read ASCII data according to a datestring pattern."""
    pipv_files = glob(os.path.join(datadir, 'PIP/a_Velocity_Tables/004%s/*2.dat' % dtstr))
    dsd_files = glob(os.path.join(datadir, 'PIP/a_DSD_Tables/004%s_a_d.dat' % dtstr))
    pluvio200_files = glob(os.path.join(datadir, 'Pluvio200/pluvio200_??_%s*.txt' % dtstr))
    pluvio400_files = glob(os.path.join(datadir, 'Pluvio400/pluvio400_??_%s*.txt' % dtstr))
    pluvio200 = read.Pluvio(pluvio200_files)
    pluvio400 = read.Pluvio(pluvio400_files)
    pipv = read.PipV(pipv_files)
    dsd = read.PipDSD(dsd_files)
    return {'vel': pipv, 'dsd': dsd,
            'pluvio200': pluvio200, 'pluvio400': pluvio400}

def batch_create_hdf(datadir='../DATA', outname='baecc.h5',
                     dtstr='20140[2-3]??'):
    """Read ASCII data and export to hdf."""
    instrdict = batch_import(dtstr, datadir)
    hdf_file = os.path.join(datadir, outname)
    for key in instrdict:
        instrdict[key].to_hdf(filename=hdf_file)

def scatterplot(x, y, c=None, kind='scatter', **kwargs):
    """scatter plot of two Series objects"""
    plotdata = pd.merge(pd.DataFrame(x), pd.DataFrame(y),
                        right_index=True, left_index=True)
    if c is not None:
        kwargs['c'] = c
    return plotdata.plot(kind=kind, x=x.name, y=y.name, **kwargs)

class EventsCollection:
    """Manage multiple events."""
    def __init__(self, csv, dtformat='%d %B %H UTC'):
        """Read event metadata from a csv file."""
        self.dtformat = dtformat
        self.events = pd.read_csv(csv, parse_dates=['start', 'end'],
                                  date_parser=self.parse_datetime)
        self.events.sort(columns=['start', 'end'], inplace=True)
        self.events.start += pd.datetools.timedelta(seconds=1)

    def parse_datetime(self, dtstr):
        date = datetime.strptime(dtstr, self.dtformat)
        date = date.replace(year=2014)
        return date

    def add_data(self, data, autoshift=True, autobias=True):
        """Add data from a Case object."""
        cases = []
        for (i, e) in self.events.iterrows():
            cases.append(data.between_datetime(e.start, e.end,
                                               autoshift=autoshift,
                                               autobias=autobias))
        self.events[data.pluvio.name] = cases

    def autoimport_data(self, datafile=['../DATA/baecc.h5'],
                        autoshift=False, autobias=False, **casekwargs):
        """Import data from a hdf file."""
        timemargin = pd.datetools.timedelta(hours=3)
        dt_start = self.events.iloc[0].start - timemargin
        dt_end = self.events.iloc[-1].end + timemargin
        data = Case.from_hdf(dt_start, dt_end, autoshift=False,
                             filenames=datafile, **casekwargs)
        for d in data:
            self.add_data(d, autoshift=autoshift, autobias=autobias)
        return

class Case(read.PrecipMeasurer, read.Cacher):
    """Calculate snowfall rate from particle size and velocity data."""
    def __init__(self, dsd, pipv, pluvio, varinterval=True, unbias=False,
                 autoshift=False, liquid=False, quess=(0.01, 2.1),
                 bnd=((0, 0.1), (1, 3)), rule='15min', use_cache=True):
        self._use_cache = use_cache
        self.dsd = dsd
        self.pipv = pipv
        self.pluvio = pluvio
        self._varinterval = varinterval
        self.pluvio.varinterval = varinterval
        self.quess = quess
        self.bnd = bnd
        if varinterval:
            self._rule = None
        else:
            self._rule = rule
        self.liquid = liquid
        self._ab = None # alpha, beta
        for instr in [self.dsd, self.pipv, self.pluvio]:
            instr.case = self
        if autoshift:
            self.autoshift()
        if unbias:
            self.noprecip_bias()

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

    @property
    def use_cache(self):
        return self._use_cache

    @use_cache.setter
    def use_cache(self, use_cache):
        self._use_cache = use_cache
        for instr in [self.dsd, self.pipv, self.pluvio]:
            instr.use_cache = use_cache

    @property
    def varinterval(self):
        return self._varinterval

    @varinterval.setter
    def varinterval(self, varinterval):
        self._varinterval = varinterval
        self.pluvio.varinterval = varinterval
        self.reset()

    @property
    def rule(self):
        if self.varinterval and self._rule is None:
            self._rule = self.pluvio.grouper() # TODO: needs to be reset on changes for pluvio data
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
    def from_hdf(cls, dt_start, dt_end, filenames=['../DATA/baecc.h5'],
                 **kwargs):
        """Create Case object from a hdf file."""
        for dt in [dt_start, dt_end]:
            dt = pd.datetools.to_datetime(dt)
        pluvio200 = read.Pluvio(filenames, hdf_table='pluvio200')
        pluvio400 = read.Pluvio(filenames, hdf_table='pluvio400')
        dsd = read.PipDSD(filenames, hdf_table='pip_dsd')
        pipv = read.PipV(filenames, hdf_table='pip_vel')
        for instr in [pluvio200, pluvio400, dsd, pipv]:
            instr.set_span(dt_start, dt_end)
        m200 = cls(dsd, pipv, pluvio200, **kwargs)
        m400 = cls(dsd, pipv, pluvio400, **kwargs)
        return m200, m400

    def between_datetime(self, dt_start, dt_end, inplace=False,
                         autoshift=False, autobias=False):
        """Select data only in chosen time frame."""
        dt_start = pd.datetools.to_datetime(dt_start)
        dt_end = pd.datetools.to_datetime(dt_end)
        if inplace:
            m = self
        else:
            m = copy.deepcopy(self)
        for instr in [m.dsd, m.pipv, m.pluvio]:
            instr.between_datetime(dt_start, dt_end, inplace=True)
            instr.case = m
        m.pluvio.bias = 0
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
        """Calculate precipitation intensity using given or saved parameters."""
        if params is None and not self.liquid:
            params = self.ab
        if self.liquid:
            fits = self.series_nans()
            fits.loc[:] = read.gunn_kinzer
            fits.name = read.gunn_kinzer.name
            self.pipv.fits = pd.DataFrame(fits)
            r = self.sum_over_d(self.r_rho, rho=RHO_W)
        elif simple:
            r = self.sum_over_d(self.r_rho, rho=params[0])
        else:
            r = self.sum_over_d(self.r_ab, alpha=params[0], beta=params[1])
        if self.varinterval:
            return r
        return r.reindex(self.pluvio.amount(rule=self.rule).index).fillna(0)

    def amount(self, **kwargs):
        """Calculate precipitation in mm using given or saved parameters."""
        i = self.intensity(**kwargs)
        if self.varinterval:
            delta = self.pluvio.tdelta()
        else:
            delta = i.index.freq.delta
        return i*(delta/pd.datetools.timedelta(hours=1))

    def sum_over_d(self, func, **kwargs):
        """numerical integration over particle diameter"""
        dD = self.dsd.d_bin
        result = self.series_zeros()
        for d in self.dsd.good_data().columns:
            result = result.add(func(d, **kwargs)*dD, fill_value=0)
        return result

    def r_ab(self, d, alpha, beta):
        """(mm/h)/(m/s)*kg/mg / kg/m**3 * mg/mm**beta * mm**beta * m/s * 1/(mm*m**3)
        """
        return 3.6/RHO_W*alpha*d**beta*self.v(d)*self.n(d)

    def r_rho(self, d, rho):
        """(mm/h)/(m/s)*m**3/mm**3 * kg/m**3 / (kg/m**3) * mm**3 * m/s * 1/(mm*m**3)
        """
        return 3.6e-3*TAU/12*rho/RHO_W*d**3*self.v(d)*self.n(d)

    def v(self, d):
        """velocity wrapper"""
        return self.pipv.v(d, varinterval=self.varinterval, rule=self.rule)

    def n(self, d):
        """N wrapper"""
        return self.dsd.n(d, varinterval=self.varinterval, rule=self.rule)

    def v_fall(self, d, how=np.median):
        """v(D) m/s for every timestep, query is slow"""
        vel = self.pipv.dbin(d, self.dsd.d_bin).vel_v
        if vel.empty:
            return self.series_nans()
        return vel.resample(self.rule, how=how, closed='right', label='right')

    def v_fall_all(self):
        """v(D) in m/s for every timestep and diameter bin"""
        v_d = []
        for d in self.dsd.good_data().columns:
            vel = self.v_fall(d)
            vel.name = d
            v_d.append(vel)
        return pd.concat(v_d, axis=1)

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
        return super().cache_dir(dt_start, dt_end, self.pluvio.name)

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
            idxd = self.dsd.good_data().columns
            dd = pd.Series(idxd)
            dD = self.dsd.d_bin
            d3n = lambda d: d**3*self.n(d)*dD
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
            idxd = self.dsd.good_data().columns
            dd = pd.Series(idxd)
            nd = dd.apply(self.n).T
            nd.columns = idxd
            dmax = nd[nd > 0.0001].T.apply(pd.Series.last_valid_index).fillna(0)
            dmax.name = name
            return dmax
        return self.msger(name, func)

    def n_moment(self, n):
        moment = lambda d: d**n*self.n(d)
        return self.sum_over_d(moment)

    def eta(self):
        return self.n_moment(4)**2/(self.n_moment(6)*self.n_moment(2))

    def mu(self):
        eta = self.eta()
        return ((7-11*eta)-np.sqrt(eta**2+14*eta+1))/(2*(eta-1))

    def lam(self):
        mu = self.mu()
        return np.sqrt(self.n_moment(2)*gamma(mu+5)/(self.n_moment(4)*gamma(mu+3)))

    def n_0(self):
        mu = self.mu()
        return self.n_moment(2)*self.lam()**(mu+3)/gamma(mu+3)

    def d_0_gamma(self):
        return (3.67+self.mu())/self.lam()

    def partcount(self):
        return self.pipv.partcount(rule=self.rule, varinterval=self.varinterval)

    def series_zeros(self):
        """Return series of zeros of the shape of timestep averaged data."""
        return self.pluvio.acc(rule=self.rule)*0

    def series_nans(self):
        """Return series of nans of the shape of timestep averaged data."""
        return self.series_zeros()*np.nan

    def noprecip_bias(self, inplace=True):
        """Wrapper to unbias pluvio using LWC calculated from PIP data."""
        return self.pluvio.noprecip_bias(self.pipv.lwc(), inplace=inplace)

    def pluvargs(self):
        args = {}
        if not self.varinterval:
            args['rule'] = self.rule
        return args

    def cost(self, c, use_accum=True):
        """Cost function for minimization"""
        if use_accum:
            pip_precip = self.acc(params=c)
            cost_method = self.pluvio.acc
        else:
            pip_precip = self.intesity(params=c)
            cost_method = self.pluvio.intensity()
        return abs(pip_precip.add(-1*cost_method(**self.pluvargs())).sum())

    def cost_lsq(self, beta):
        """Single variable cost function using lstsq to find linear coef."""
        alpha = self.alpha_lsq(beta)
        return self.cost([alpha, beta])

    def const_lsq(self, c, simple):
        acc_arr = self.acc(params=c, simple=simple).values
        A = np.vstack([acc_arr, np.ones(len(acc_arr))]).T
        y = self.pluvio.acc(**self.pluvargs()).values
        return np.linalg.lstsq(A, y)[0][0]

    def alpha_lsq(self, beta):
        """Wrapper for const_lsq to calculate alpha"""
        return self.const_lsq(c=[1, beta], simple=False)

    def density_lsq(self):
        """Wrapper for const_lsq to calculate least square particle density"""
        return self.const_lsq(c=[1], simple=True)

    def density(self, pluvio_filter=True, pip_filter=False):
        """Calculates mean density estimate for each timeframe."""
        name = 'density'
        def func():
            rho_r_pip = self.amount(params=[1], simple=True)
            if pluvio_filter: #filter
                rho_r_pip[self.pluvio.intensity() < 0.1] = np.nan
            if pip_filter and self.ab is not None:
                rho_r_pip[self.intensity() < 0.1] = np.nan
            rho = self.pluvio.amount(rule=self.rule)/rho_r_pip
            rho.name = name
            return rho.replace(np.inf, np.nan)
        return self.msger(name, func)

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
        t = self.time_range()
        return (t[0], t[-1])

    def time_range(self):
        """data time ticks on minute interval"""
        return pd.date_range(self.pluvio.acc().index[0],
                             self.pluvio.acc().index[-1], freq='1min')

    def plot(self, axarr=None, kind='line', label_suffix='', pip=True, **kwargs):
        """Plot calculated (PIP) and pluvio intensities."""
        if axarr is None:
            f, axarr = plt.subplots(4, sharex=True, dpi=120)
        if pip:
            self.intensity().plot(label='PIP ' + label_suffix, kind=kind, ax=axarr[0], **kwargs)
        self.pluvio.intensity(rule=self.rule).plot(label=self.pluvio.name + ' ' + label_suffix,
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

    def plot_v_binned(self, ax=None, **kwargs):
        """Plot velocity in diameter bins."""
        if ax is None:
            ax = plt.gca()
        diam = []
        vel = []
        for d in self.dsd.good_data().columns:
            v_new = self.v_fall(d).values
            d_new = [d]*len(v_new)
            vel.extend(v_new)
            diam.extend(d_new)
        ax.plot(diam, vel, 'h', **kwargs)
        return ax

    def plot_v_stuff(self, ax=None, **kwargs):
        """Plot a lot of velocity related stuff in a single figure."""
        if ax is None:
            ax = plt.gca()
        self.plot_v_binned(label='%s bin median' % self.rule, ax=ax, zorder=3,
                           **kwargs)
        self.v_fall_all().mean().plot(label='timestep mean', ax=ax, zorder=4,
                                      **kwargs)
        self.pipv.plot(ax=ax, zorder=2, **kwargs)
        ax.legend(loc='lower right')
        return ax

    def plot_velfitcoefs(self, fig=None, ax=None, rhomin=None, rhomax=None, countmin=1, **kwargs):
        rho = self.density()
        params = self.pipv.fits.polfit.apply(lambda fit: fit.params)
        selection = pd.DataFrame([rho.notnull(), self.partcount() > countmin]).all()
        rho = rho[selection]
        params = params[selection]
        a = params.apply(lambda p: p[0]).values
        b = params.apply(lambda p: p[1]).values
        if fig is None:
            fig = plt.figure(dpi=120)
        if ax is None:
            ax = plt.gca()
        if rhomin is None:
            vmin = rho.min()
        if rhomax is None:
            vmax = rho.max()
        choppa = ax.scatter(a, b, c=rho.values, vmin=rhomin, vmax=rhomax,
                            **kwargs)
        cb = fig.colorbar(choppa, label='bulk density')
        ax.set_xlabel('$a_u$', fontsize=15)
        ax.set_ylabel('$b_u$', fontsize=15)
        return ax

    def plot_d0_bv(self, rhomin=None, rhomax=None, countmin=1, **kwargs):
        rho = self.density()
        params = self.pipv.fits.polfit.apply(lambda fit: fit.params)
        selection = pd.DataFrame([rho.notnull(), self.partcount() > countmin]).all()
        rho = rho[selection]
        params = params[selection]
        dmax = self.d_max()[selection]
        a = params.apply(lambda p: p[0])
        b = params.apply(lambda p: p[1])
        b.name = 'b'
        if rhomin is None:
            vmin = rho.min()
        if rhomax is None:
            vmax = rho.max()
        return scatterplot(x=dmax, y=b, c=rho, **kwargs)

    def xcorr(self, rule='1min', ax=None, **kwargs):
        """Plot cross-correlation between lwc estimate and pluvio intensity.
        Extra arguments are passed to pyplot.xcorr.
        """
        if ax is None:
            ax = plt.gca()
        r = self.pluvio.intensity(rule=rule, unbias=False)
        lwc = self.pipv.lwc(rule).reindex(r.index).fillna(0)
        return ax.xcorr(lwc, r, **kwargs)

    def autoshift(self, rule='1min', inplace=False):
        """Find and correct pluvio time shift using cross correlation."""
        if self.pluvio.shift_periods != 0:
            print('Pluvio already timeshifted, resetting.')
            self.pluvio.shift_reset()
        xc = self.xcorr(rule=rule)
        imaxcorr = xc[1].argmax()
        periods = xc[0][imaxcorr]
        if inplace:
            self.pluvio.shift_periods = periods
            self.pluvio.shift_freq = rule
            print('Pluvio timeshift set to %s*%s.'
                  % (str(self.pluvio.shift_periods), self.pluvio.shift_freq))
        return periods

class Snow2:
    """UNTESTED.
    Calculate snowfall rate using Szyrmer Zawadski's method from Snow Study II.
    """
    def __init__(self):
        return

    @staticmethod
    def best(re, mh=True):
        if mh: # MH05
            cl = np.array([3.8233, -1.5211, 0.30065, -0.06104, 0.13074,
                           -0.073429, 0.016006, -0.0012483])
        else: # KC05
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

