import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.vtk_io import readVTK
import numpy as np
import inifix
from scipy.special import jn

conf = inifix.load("../idefix.ini")
amplitude = conf["Setup"]["amplitude"]
time_step = conf["Output"]["vtk"]

nfile = 10
list_time = [round(i*time_step, 7) for i in range(nfile)]
list_VTK = list()
for k in range(nfile):
    current_VTK = readVTK('../data.00{0:02d}.vtk'.format(k), 'spherical')
    list_VTK.append(current_VTK)

list_VX = list()
list_VY = list()
list_VZ = list()
list_RHO = list()
list_BX = list()
list_BY = list()
list_BZ = list()
for VTK in list_VTK:
    list_VX.append(VTK.data['VX1'])
    list_VY.append(VTK.data['VX2'])
    list_VZ.append(VTK.data['VX3'])
    list_BX.append(VTK.data['BX1'])
    list_BY.append(VTK.data['BX2'])
    list_BZ.append(VTK.data['BX3'])
    list_RHO.append(VTK.data['RHO'])
X = VTK.r
Y = VTK.theta
Z = VTK.phi

mu = conf["Hydro"]["bragViscosity"][-1]

def analytic_sol(r, th, phi, t):
    lambd2 = 2.
    return amplitude*jn(1,r)/r*np.exp(-t*lambd2*mu)


eps = 7e-7
for k, Vx in enumerate(list_VX):
    t = list_time[k]
    current_sol = np.zeros((len(X),len(Y),len(Z)))
    for ix,x in enumerate(X):
        for iy,y in enumerate(Y):
            for iz,z in enumerate(Z):
                current_sol[ix,iy,iz] = analytic_sol(x,y,z,t)
    if np.mean(abs(current_sol - Vx)) > eps:
        print("Failed")
        sys.exit(1)

print("SUCCESS")
sys.exit(0)
