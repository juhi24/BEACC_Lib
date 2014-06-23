"""Tools for estimating density and other properties of falling snow"""
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

TAU = 2*np.pi

class Method1:
    """Calculate snowfall rate from particle size and velocity data."""
    def __init__(self, dsd, pipv, pluvio, quess=(0.005, 2.1), bnd=((0, 0.1), (1, 3)), rule='15min'):
        self.dsd = dsd
        self.pipv = pipv
        self.pluvio = pluvio
        self.rho_w = 1000
        self.quess = quess
        self.bnd = bnd
        self.rule = rule
        self.result = None
        self.ab = None
        
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
        print('Optimizing constants...')
        self.result = minimize(self.cost_lsq, self.quess[1], method='Nelder-Mead')
        #self.result = minimize(self.cost_lsq, self.quess[1], method='SLSQP', bounds=self.bnd[1])
        beta = self.result.x[0]
        alpha = self.alpha_lsq(beta)
        self.ab = [alpha, beta]
        return self.result
        
    def plot(self):
        """Plot calculated (PIP) and pluvio rainrates."""
        if self.ab is None:
            print('Constants not defined. Will now find them via minimization.')
            self.minimize_lsq()
        kind = 'line'
        f, axarr = plt.subplots(2, sharex=True)
        self.rainrate().plot(label='PIP', kind=kind,ax=axarr[0])
        self.pluvio.rainrate(self.rule).plot(label=self.pluvio.name, kind=kind, ax=axarr[0])
        axarr[0].set_ylabel('mm/%s' % self.rule)
        axarr[0].set_title(r'%s rainrate, $\alpha=%s, \beta=%s$' % (self.rule, self.ab[0], self.ab[1]))
        rho = 1e6*self.density()
        rho.plot(label='mean density', ax=axarr[1])
        for ax in axarr:
            ax.legend(loc='upper right')
            ax.set_xlabel('time')
        axarr[1].set_ylabel(r'$\rho_{part}$')
    
    def plot_cost(self, resolution=20):
        """The slowest plot you've made"""
        if self.ab is None:
            return
        alpha0 = self.ab[0]
        alpha = np.linspace(0.4*alpha0, 1.4*alpha0, num=resolution)
        beta = np.linspace(self.bnd[1][0], self.bnd[1][1],num=resolution)
        z = np.zeros((alpha.size, beta.size))
        for i, a in enumerate(alpha):
            for j, b in enumerate(beta):
                z[i][j] = self.cost((a, b))
        plt.clf()
        plt.pcolor(beta, alpha, z, cmap='binary')
        plt.colorbar()
        plt.xlabel(r'$\beta$')
        plt.ylabel(r'$\alpha$')
        plt.axis('tight')
        plt.title('cost function value')
        plt.plot(self.ab[1], self.ab[0], 'ro')
        return z
        
    def plot_cost_lsq(self, resolution):
        """Plot cost function value."""
        beta = np.linspace(self.bnd[1][0],self.bnd[1][1],num=resolution)
        cost = np.array([self.cost_lsq(b) for b in beta])
        plt.clf()        
        plt.plot(beta,cost)
        plt.xlabel(r'$\beta$')
        plt.ylabel('cost')
        plt.title('cost function value')

class Snow2:
    """UNTESTED. Calculate snowfall rate using Szyrmer Zawadski's method from Snow Study II."""
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

