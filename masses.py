from glob import glob
import read
import numpy as np
from scipy.optimize import minimize
import copy
import matplotlib.pyplot as plt

class Method1:
    def __init__(self, dsd, pipv, pluvio, quess=(0.005,2.1), bnd=((0,0.1),(1,3)), rule='30min'):
        self.dsd = dsd
        self.pipv = pipv
        self.pluvio = pluvio
        self.rho_w = 1000
        self.quess = quess
        self.bnd = bnd
        self.rule = rule
        self.result = None
        
    def rainrate(self, alphabeta=None):
        if self.result is not None and alphabeta is None:
            alphabeta = self.result.x
        dD = self.dsd.d_bin
        R = None
        for D in self.dsd.bin_cen():
            vcond = 'Wad_Dia > %s and Wad_Dia < %s' % (D-0.5*dD, D+0.5*dD)
            vel = self.pipv.data.query(vcond).vel_v.mean() # V(D_i) m/s
            N = self.dsd.data[str(D)].resample(self.rule, how=self.dsd._sum) # N(D_i) 1/(mm*m**3)
            addition = 3.6/self.rho_w*alphabeta[0]*D**alphabeta[1]*vel*N*dD
            # (mm/h)/(m/s) * kg/m**3 * kg/m**beta * m**beta * m/s * 1/(mm*m**3) * mm == mm/h
            if R is None:
                R = addition
            else:
                R = R.add(addition, fill_value=0)
        return R
        
    def cost(self,c):
        dsd_acc = self.rainrate(c).cumsum()
        return abs(dsd_acc.add(-1*self.pluvio.acc().dropna()).sum())
        
    def minimize(self):
        self.result = minimize(self.cost, self.quess, method='SLSQP', bounds=self.bnd)
        return self.result
        
    def plot(self):
        kind = 'line'
        ax = self.rainrate().plot(label='PIP',kind=kind)
        self.pluvio.rainrate(self.rule).plot(label=self.pluvio.name,kind=kind,ax=ax)
        ax.legend()
        ax.set_xlabel('time')
        ax.set_ylabel('mm')
        ax.set_title(r'%s rainrate, $\alpha=%s, \beta=%s$' % (self.rule, self.result.x[0], self.result.x[1]))

class Snow2:
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
        
pipv_files = glob('/home/jussitii/DATA/PIP/a_Velocity_Tables/00420140212/*.dat')
dsd_files = glob('/home/jussitii/DATA/PIP/a_DSD_Tables/00420140212_a_d.dat')
pluvio200_files = glob('/home/jussitii/DATA/Pluvio200/pluvio200_0?_20140212*')
pluvio400_files = glob('/home/jussitii/DATA/Pluvio400/pluvio400_0?_20140212*')

#dsd_files = glob('/home/jussitii/DATA/PIP/a_DSD_Tables/0042014021[5-6]_a_d.dat')
#pipv_files = glob('/home/jussitii/DATA/PIP/a_Velocity_Tables/0042014021[5-6]/*.dat')
#pluvio400_files = glob('/home/jussitii/DATA/Pluvio400/pluvio400_0?_2014021[5-6]*')
#pluvio200_files = glob('/home/jussitii/DATA/Pluvio200/pluvio200_0?_2014021[5-6]*')

pluvio200 = read.Pluvio(pluvio200_files)
pluvio400 = read.Pluvio(pluvio400_files)
pipv = read.PipV(pipv_files)
dsd = read.PipDSD(dsd_files)

time_start = '4:00'
time_end = '8:00'

dt_start = '20140215T21:00'
dt_end = '20140216T01:00'
pluvio200_lim = copy.deepcopy(pluvio200)
pluvio400_lim = copy.deepcopy(pluvio400)
dsd_lim = dsd
for instr_lim in [pluvio200_lim,pluvio400_lim,dsd_lim]:
    instr_lim.data = instr_lim.data.between_time(time_start,time_end)
    #instr_lim.data = instr_lim.data[dt_start:dt_end]

m_lim200 = Method1(dsd_lim,pipv,pluvio200_lim,rule='30min')
m_lim400 = Method1(dsd_lim,pipv,pluvio400_lim,rule='30min')

#res1 = minimize(m.cost, quess, method='SLSQP', bounds=bnd)
#m.rule = '30min'
#res30 = minimize(m.cost, quess, method='SLSQP', bounds=bnd)
#res30_lim = minimize(m_lim.cost, quess, method='SLSQP', bounds=bnd)
#m.pluvio = pluvio400
#res400 = minimize(m.cost, quess, method='SLSQP', bounds=bnd)
#res400_lim = minimize(m_lim.cost, quess, method='SLSQP', bounds=bnd)