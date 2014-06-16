import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

TAU = 2*np.pi

class Method1:
    """Calculate snowfall rate from particle size and velocity data."""
    def __init__(self, dsd, pipv, pluvio, quess=(0.005,2.1), bnd=((0,0.1),(1,3)), rule='30min'):
        self.dsd = dsd
        self.pipv = pipv
        self.pluvio = pluvio
        self.rho_w = 1000
        self.quess = quess
        self.bnd = bnd
        self.rule = rule
        self.result = None
        
    def rainrate(self, consts=None, simple=False):
        """Calculate rainrate using given or saved constants."""
        if self.result is not None and consts is None:
            consts = self.result.x
        dD = self.dsd.d_bin
        R = None
        for D in self.dsd.bin_cen():
            vcond = 'Wad_Dia > %s and Wad_Dia < %s' % (D-0.5*dD, D+0.5*dD)
            vel = self.pipv.data.query(vcond).vel_v.mean() # V(D_i) m/s, query is slow
            N = self.dsd.data[str(D)].resample(self.rule, how=self.dsd._sum) # N(D_i) 1/(mm*m**3)
            if simple:
                addition = TAU/12*consts[0]*D**3*vel*N*dD
            else:
                addition = 3.6/self.rho_w*consts[0]*D**consts[1]*vel*N*dD
                # (mm/h)/(m/s) * mg/m**3 * mg/m**beta * m**beta * m/s * 1/(mm*m**3) * mm == mm/h
            if R is None:
                R = addition
            else:
                R = R.add(addition, fill_value=0)
        return R
        
    def cost(self,c):
        """Cost function for minimization"""
        dsd_acc = self.rainrate(c).cumsum()
        return abs(dsd_acc.add(-1*self.pluvio.acc().dropna()).sum())
        
    def minimize(self):
        """Find constants for calculating particle masses. Save and return results."""
        print('Optimizing constants...')
        self.result = minimize(self.cost, self.quess, method='SLSQP', bounds=self.bnd)
        return self.result
        
    def plot(self):
        """Plot calculated (PIP) and pluvio rainrates."""
        kind = 'line'
        ax = self.rainrate().plot(label='PIP',kind=kind)
        self.pluvio.rainrate(self.rule).plot(label=self.pluvio.name,kind=kind,ax=ax)
        ax.legend()
        ax.set_xlabel('time')
        ax.set_ylabel('mm')
        ax.set_title(r'%s rainrate, $\alpha=%s, \beta=%s$' % (self.rule, self.result.x[0], self.result.x[1]))
    
    def plot_cost(self,resolution=20):
        """The slowest plot you've made"""
        if self.result is None:
            return
        alpha = np.linspace(0,2*self.result.x[0],num=resolution)
        beta = np.linspace(self.bnd[1][0],self.bnd[1][1],num=resolution)
        z = np.zeros((alpha.size,beta.size))
        for i,a in enumerate(alpha):
            for j,b in enumerate(beta):
                z[i][j] = self.cost((a,b))
        plt.clf()
        plt.pcolor(beta,alpha,z,cmap='binary')
        plt.colorbar()
        plt.xlabel(r'$\beta$')
        plt.ylabel(r'$\alpha$')
        plt.axis('tight')
        plt.title('Cost function value')
        plt.plot(self.result.x[1],self.result.x[0],'ro')

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
    def mass(u,ar,D):
        g = 9.81
        fa = 1
        rho_a = 1.275
        nu_a = 1.544e-5
        
        re = u*D/nu_a
        return np.pi*rho_a*nu_a**2/(8*g)*Snow2.best(re)*ar*fa

