# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 09:22:58 2015

@author: dori
"""

import numpy as np
from scipy.special import cbrt
import matplotlib.pyplot as plt

plt.ioff()
#plt.ion()

def rotateX(theta):
    """ Rotation matrix aroun X axis """
    st = np.sin(theta)
    ct = np.cos(theta)
    return np.matrix(((1.0, 0.0, 0.0),
                      (0.0,  ct, -st),
                      (0.0,  st,  ct)))

def rotateZ(phi):
    """ Rotation matrix aroun Z axis """
    sp = np.sin(phi)
    cp = np.cos(phi)
    return np.matrix((( cp, -sp, 0.0),
                      ( sp,  cp, 0.0),
                      (0.0, 0.0, 1.0)))

# Define rotation of the ellipsoid
gauss = np.random.normal(0.0, np.pi*9.0/180.0, 200)
uni = np.random.uniform(low=0.0, high=2*np.pi, size=300)

# Containers for partial results
phif = np.zeros((len(gauss),len(uni)))
phiw = np.zeros((len(gauss),len(uni)))
phiw_min = np.zeros((len(gauss),len(uni)))
phiec = np.zeros((len(gauss),len(uni)))
ec2vol = np.zeros((len(gauss),len(uni)))

# Define list of horizontal aspect ratios b/a
arh = np.linspace(0.1,1,10)
# Define list of vertical aspect ratios c/a
arv = np.linspace(0.1,1,19)

# containers for definitive results (as a function of axis ratios)
mean_phif = np.zeros((len(arh),len(arv)))
mean_phiw = np.zeros((len(arh),len(arv)))
mean_phiec = np.zeros((len(arh),len(arv)))
mean_ec2vol = np.zeros((len(arh),len(arv)))
std_phif = np.zeros((len(arh),len(arv)))
std_phiw = np.zeros((len(arh),len(arv)))
std_phiec = np.zeros((len(arh),len(arv)))
std_ec2vol = np.zeros((len(arh),len(arv)))

for l in range(len(arh)):
    for m in range(len(arv)):
        if not (arv[m] > arh[l]):
            # Define semi-axis of the ellipsoid
            a = 0.5
            b = a*arh[l]
            c = a*arv[m]
            
            # Find 3d dimensions
            Dmax = 2.0*a
            Dvol = 2.0*cbrt(a*b*c)
            
            Aa = 1.0/(a*a)
            Ab = 1.0/(b*b)
            Ac = 1.0/(c*c)
            A = np.matrix((( Aa, 0.0, 0.0), 
                           (0.0,  Ab, 0.0), 
                           (0.0, 0.0,  Ac)))
            
            for i in range(len(gauss)):
                for j in range(len(uni)): 
                    R = rotateZ(uni[j])*rotateX(gauss[i])
                    Ar = R.T*A*R
                    
                    # Matrix of quadratic terms of projected ellipse on zx plane
                    B = np.matrix(((Ar[0,0]-Ar[0,1]*Ar[0,1]/Ar[1,1], Ar[0,2]-Ar[0,1]*Ar[1,2]/Ar[1,1]),
                                   (Ar[0,2]-Ar[0,1]*Ar[1,2]/Ar[1,1], Ar[2,2]-Ar[1,2]*Ar[1,2]/Ar[1,1])))
                    eigval, eigvec = np.linalg.eigh(B)
                    
                    Dsvif = 2.0/(np.sqrt(min(eigval)))
                    Dsvif_min = 2.0/(np.sqrt(max(eigval)))
                    Dsviec = 2.0/np.sqrt(np.sqrt(eigval[0]*eigval[1]))
                    Dsviw = max(Dsvif*np.cos(np.arccos(eigvec[1,1])),Dsvif*np.sin(np.arccos(eigvec[1,1])))
                    phif[i,j] = Dsvif/Dmax
                    phiw[i,j] = Dsviw/Dmax
                    phiw_min[i,j] = min(Dsvif*np.cos(np.arccos(eigvec[1,1])),Dsvif*np.sin(np.arccos(eigvec[1,1])))
                    phiec[i,j] = Dsviec/Dmax
                    ec2vol[i,j] = Dsviec/Dvol
            mean_phif[l,m] = np.mean(phif)
            mean_phiw[l,m] = np.mean(phiw)
            mean_phiec[l,m] = np.mean(phiec)
            mean_ec2vol[l,m] = np.mean(ec2vol)
            std_phif[l,m] = np.std(phif)
            std_phiw[l,m] = np.std(phiw)
            std_phiec[l,m] = np.std(phiec)
            std_ec2vol[l,m] = np.std(ec2vol)
            print(arv[m],arh[l],np.mean(phif),np.mean(phiw),np.mean(phiw_min),np.mean(phiec))
        else:
            mean_phif[l,m] = np.nan
            mean_phiw[l,m] = np.nan
            mean_phiec[l,m] = np.nan
            mean_ec2vol[l,m] = np.nan
            std_phif[l,m] = np.nan
            std_phiw[l,m] = np.nan
            std_phiec[l,m] = np.nan
            std_ec2vol[l,m] = np.nan
#%%
plf = plt.figure()
for k in range(len(arh)):
    plt.errorbar(arv,mean_phif[k],yerr=std_phif[k],marker='*',label=str(arh[k]))
plt.legend(title='b/a',loc=0)
plt.title('D_{svi,f}')
plt.xlabel('c/a')
plt.ylabel('\phi')
plt.savefig('phif.eps')
plw = plt.figure()
for k in range(len(arh)):
    plt.errorbar(arv,mean_phiw[k],yerr=std_phiw[k],marker='*',label=str(arh[k]))
plt.legend(title='b/a',loc=0)
plt.title('D_{svi,w}')
plt.xlabel('c/a')
plt.ylabel('\phi')
plt.savefig('phiw.eps')
plec = plt.figure()
for k in range(len(arh)):
    plt.errorbar(arv,mean_phiec[k],yerr=std_phiec[k],marker='*',label=str(arh[k]))
plt.legend(title='b/a',loc=0)
plt.title('D_{svi,ec}')
plt.xlabel('c/a')
plt.ylabel('\phi')
plt.savefig('phiec.eps')
plec2vol = plt.figure()
for k in range(len(arh)):
    plt.errorbar(arv,mean_ec2vol[k],yerr=std_ec2vol[k],marker='*',label=str(arh[k]))
plt.legend(title='b/a',loc=0)
plt.title('D_{svi,ec to equal volume}')
plt.xlabel('c/a')
plt.ylabel('\phi')
plt.savefig('ec2vol.eps')