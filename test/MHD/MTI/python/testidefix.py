"""
Created in July 2021

@author: Jean Kempf, Francois Rincon

The test is inspired from the following paper:
Parrish, Ian J., et al. "The effects of anisotropic viscosity on turbulence and heat transport in the intracluster medium." Monthly Notices of the Royal Astronomical Society 422.1 (2012): 704-718.
"""

import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))
import numpy as np
import inifix

n = 1 # fondamental mode

conf = inifix.load("../idefix.ini")
Pr=conf["Setup"]["pr"]
ksi=conf["Setup"]["ksi"]

L = 0.1
H = 3.
wbuoy = 1/np.sqrt(3)
gamma = 5/3

kn = 2*np.pi*n*np.sqrt(2)/L

fid=open('../average.dat',"r")
# read the first line to get data names
varnames=fid.readline().split()
fid.close()
# load the bulk of the file
data=np.loadtxt('../average.dat',skiprows=1)
# store this in our data structure
V={}
i=0
for name in varnames:
    V[name]=data[:,i]
    i=i+1

sigma_vx = 0.5*(np.log(V['kinx'][20]) - np.log(V['kinx'][15]))/(V['t'][20] - V['t'][15])
sigma_vy = 0.5*(np.log(V['kiny'][20]) - np.log(V['kiny'][15]))/(V['t'][20] - V['t'][15])
sigma_by = 0.5*(np.log(V['by_2'][20]) - np.log(V['by_2'][15]))/(V['t'][20] - V['t'][15])
sigma_num = (sigma_vx + sigma_vy + sigma_by)/3.

def dispMTI(wcond, k, Pr):
    cosk_2 = 0.5 # here cosk_2 = sqrt(bhat) dot khat i.e. xhat dot khat
    k_2 = k**2
    ksi = 2.5/k_2/cosk_2*wcond
    kx_2 = k_2*cosk_2
    K = kx_2/k_2

    bk_2 = k_2*cosk_2
    bkhat_2 = cosk_2
    nu = Pr*ksi
    wvisc = 3*nu*bk_2
    v = wvisc*(1 - bkhat_2)
    N_2 = wbuoy**2/gamma*((H-1)*gamma-H)
    return np.array([1, wcond+v, wcond*v+N_2*kx_2/k_2, K*wcond*-wbuoy**2])

def sigmaMTI(wcond, k, Pr):
    roots = np.roots(dispMTI(wcond,k,Pr))
    return np.real(roots[np.real(roots) > 0])[0]

wcd = 0.4*ksi*kn**2/2
sigma_analytic = sigmaMTI(wcd, kn, Pr)

if np.abs(np.log(sigma_analytic/sigma_num)/np.log(sigma_num)) >= 0.1:
    print("Failed")
    sys.exit(1)
else:
    print("SUCCESS")
    sys.exit(0)
