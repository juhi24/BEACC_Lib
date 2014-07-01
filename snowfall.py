"""Tools for estimating density and other properties of falling snow"""
import numpy as np
import pandas as pd
import read
from scipy.optimize import minimize
import matplotlib.pyplot as plt

TAU = 2*np.pi

class Method1:
    """Calculate snowfall rate from particle size and velocity data."""
    def __init__(self, dsd, pipv, pluvio, unbias=False,
                 quess=(0.005, 2.1), bnd=((0, 0.1), (1, 3)), rule='15min'):
        self.dsd = dsd
        self.pipv = pipv
        self.pluvio = pluvio
        self.rho_w = 1000
        self.quess = quess
        self.bnd = bnd
        self.rule = rule
        self.result = None
        self.ab = None
        if unbias:
            self.noprecip_bias()
        
    @classmethod
    def from_hdf(cls, dt_start, dt_end, filenames=['../DATA/baecc.h5'], **kwargs):
        """Create Method1 object from a hdf file."""
        pluvio200 = read.Pluvio(filenames, hdf_table='pluvio200')
        pluvio400 = read.Pluvio(filenames, hdf_table='pluvio400')
        dsd = read.PipDSD(filenames, hdf_table='pip_dsd')
        pipv = read.PipV(filenames, hdf_table='pip_vel')
        for instr in [pluvio200, pluvio400, dsd, pipv]:
            instr.data = instr.data[dt_start:dt_end]
        m200 = cls(dsd, pipv, pluvio200, **kwargs)
        m400 = cls(dsd, pipv, pluvio400, **kwargs)
        return m200, m400
        
    def rainrate(self, consts=None, simple=False):
        """Calculate rainrate using given or saved constants."""
        if self.ab is not None and consts is None:
            consts = self.ab
        dD = self.dsd.d_bin
        R = None
        for D in self.dsd.bin_cen():
            vcond = 'Wad_Dia > %s and Wad_Dia < %s' % (D-0.5*dD, D+0.5*dD)
            vel = self.pipv.data.query(vcond).vel_v
            if vel.empty:
                continue
            # V(D_i) m/s, query is slow
            vel_down = vel.resample(self.rule, how=np.mean, closed='right', label='right')
            # N(D_i) 1/(mm*m**3)
            N = self.dsd.data[str(D)].resample(self.rule, how=self.dsd._sum, closed='right', label='right') 
            if simple:
                addition = 3.6*TAU/12*consts[0]*D**3*vel_down*N*dD
            else:
                addition = 3.6/self.rho_w*consts[0]*D**consts[1]*vel_down*N*dD
                # (mm/h)/(m/s) * kg/m**3 * mg/m**beta * m**beta * m/s * 1/(mm*m**3) * mm == mm/h
            if R is None:
                R = addition
            else:
                R = R.add(addition, fill_value=0)
        return R.fillna(0)
        
    def noprecip_bias(self):
        """Wrapper to unbias pluvio using LWC calculated from PIP data."""
        self.pluvio.noprecip_bias(self.pipv.lwc(), inplace=True)
        
    def cost(self, c):
        """Cost function for minimization"""
        pip_acc = self.rainrate(c).cumsum()
        return abs(pip_acc.add(-1*self.pluvio.acc().dropna()).sum())
        
    def cost_lsq(self, beta):
        """Single variable cost function using lstsq to find linear coef."""
        alpha = self.alpha_lsq(beta)
        return self.cost([alpha, beta])
    
    def const_lsq(self, c, simple):
        acc_arr = self.rainrate(consts=c, simple=simple).cumsum().values
        A = np.vstack([acc_arr, np.ones(len(acc_arr))]).T
        y = self.pluvio.acc(self.rule).values
        return np.linalg.lstsq(A, y)[0][0]
        
    def alpha_lsq(self, beta):
        """Wrapper for const_lsq to calculate alpha"""
        return self.const_lsq(c=[1, beta], simple=False)
        
    def density_lsq(self):
        """Wrapper for const_lsq to calculate mean particle density"""
        return self.const_lsq(c=[1], simple=True)
        
    def density(self):
        """Calculates mean density estimate for each timeframe."""
        rho_r_pip = self.rainrate([1], True)
        rho_r_pip[rho_r_pip < 1000] = np.nan # filter
        return self.pluvio.rainrate(self.rule)/rho_r_pip
        
    def minimize(self, method='SLSQP'):
        """Find constants for calculating particle masses. Save and return results."""
        print('Optimizing constants...')
        self.result = minimize(self.cost, self.quess, method=method, bounds=self.bnd)
        self.ab = self.result.x
        return self.result
        
    def minimize_lsq(self):
        """Find beta by minimization and alpha by linear least square."""
        print('Optimizing constants...')
        self.result = minimize(self.cost_lsq, self.quess[1], method='Nelder-Mead')
        #self.result = minimize(self.cost_lsq, self.quess[1], method='SLSQP', bounds=self.bnd[1])
        print(self.result.message)
        beta = self.result.x[0]
        alpha = self.alpha_lsq(beta)
        self.ab = [alpha, beta]
        return self.result
        
    def time_range(self):
        """data time ticks on minute interval"""
        return pd.date_range(self.pluvio.data.index[0], self.pluvio.data.index[-1], freq='1min')
        
    def plot(self, kind='line', **kwargs):
        """Plot calculated (PIP) and pluvio rainrates."""
        if self.ab is None:
            print('Constants not defined. Will now find them via minimization.')
            self.minimize_lsq()
        f, axarr = plt.subplots(2, sharex=True)
        self.rainrate().plot(label='PIP', kind=kind, ax=axarr[0], **kwargs)
        self.pluvio.rainrate(self.rule).plot(label=self.pluvio.name, kind=kind, ax=axarr[0], **kwargs)
        axarr[0].set_ylabel('mm/%s' % self.rule)
        axarr[0].set_title(r'%s rainrate, $\alpha=%s, \beta=%s$' % (self.rule, self.ab[0], self.ab[1]))
        rho = 1e6*self.density()
        rho.plot(label='mean density', ax=axarr[1])
        for ax in axarr:
            ax.legend(loc='upper right')
            ax.set_xlabel('time')
        axarr[1].set_ylabel(r'$\rho_{part}$')
        plt.show()
    
    def plot_cost(self, resolution=20, ax=None, cmap='binary', **kwargs):
        """The slowest plot you've made"""
        if self.ab is None:
            return
        if ax is None:
            ax = plt.gca()
        alpha0 = self.ab[0]
        alpha = np.linspace(0.4*alpha0, 1.4*alpha0, num=resolution)
        beta = np.linspace(self.bnd[1][0], self.bnd[1][1],num=resolution)
        z = np.zeros((alpha.size, beta.size))
        for i, a in enumerate(alpha):
            for j, b in enumerate(beta):
                z[i][j] = self.cost((a, b))
        ax = plt.gca()
        heat = ax.pcolor(beta, alpha, z, cmap=cmap, **kwargs)
        ax.colorbar()
        ax.xlabel(r'$\beta$')
        ax.ylabel(r'$\alpha$')
        ax.axis('tight')
        ax.title('cost function value')
        return z, heat, ax.plot(self.ab[1], self.ab[0], 'ro')
        
    def plot_cost_lsq(self, resolution, ax=None, *args, **kwargs):
        """Plot cost function value."""
        if ax is None:
            ax = plt.gca()
        beta = np.linspace(self.bnd[1][0],self.bnd[1][1],num=resolution)
        cost = np.array([self.cost_lsq(b) for b in beta])
        ax =  plt.gca()
        ax.xlabel(r'$\beta$')
        ax.ylabel('cost')
        ax.title('cost function value')
        return ax.plot(beta, cost, *args, **kwargs)
        
    def xcorr(self, rule='1min', ax=None, **kwargs):
        """Plot cross-correlation between lwc estimate and pluvio rainrate. 
        Extra arguments are passed to pyplot.xcorr.
        """
        if ax is None:
            ax = plt.gca()
        r = self.pluvio.rainrate(rule)
        lwc = self.pipv.lwc(rule).reindex(self.pluvio.data.index).fillna()
        return plt.xcorr(lwc, r, **kwargs)

class Snow2:
    """UNTESTED. 
    Calculate snowfall rate using Szyrmer Zawadski's method from Snow Study II.
    """
    def __init__():
        return

    @staticmethod
    def best(re, mh=True):
        if mh: # MH05
            cl = np.array([3.8233, -1.5211, 0.30065, -0.06104, 0.13074, -0.073429, 0.016006, -0.0012483])
        else: # KC05
            cl = np.array([3.8816, -1.4579, 0.27749, -0.41521, 0.57683, -0.29220, 0.06467, -0.0053405])    
        logx = 0    
        
        for l, c in enumerate(cl):
            logx += c*np.log(re)**l
        
        return np.exp(logx)
        
    @staticmethod
    def mass(u, ar, D):
        g = 9.81
        fa = 1
        rho_a = 1.275
        nu_a = 1.544e-5
        
        re = u*D/nu_a
        return np.pi*rho_a*nu_a**2/(8*g)*Snow2.best(re)*ar*fa

