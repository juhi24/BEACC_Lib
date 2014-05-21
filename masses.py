from glob import glob
import read
import numpy as np

V = 0.5
rho_w = 1000

pipv_files = glob('/home/jussitii/DATA/PIP/a_Velocity_Tables/00420140212/*.dat')
dsd_files = glob('/home/jussitii/DATA/PIP/a_DSD_Tables/00420140212_a_d.dat')
pluvio_files = glob('/home/jussitii/DATA/Pluvio200/pluvio200_02_20140212*')

pluvio = read.Pluvio(pluvio_files)
pipv = read.PipV(pipv_files)
dsd = read.PipDSD(dsd_files)

class Simple:
    def __init__():
        return
    
    def mass():
        return
        
    

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
        
