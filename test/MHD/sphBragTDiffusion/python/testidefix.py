import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.vtk_io import readVTK
import numpy as np
import inifix

conf = inifix.load("../idefix.ini")
amplitude = conf["Setup"]["amplitude"]
gamma = conf["Hydro"]["gamma"]
kappa = conf["Hydro"]["bragTDiffusion"][-1]
time_step = conf["Output"]["vtk"]

nfile = 6
list_time = [round(i*time_step, 3) for i in range(nfile)]
list_VTK = list()
for k in range(nfile):
    if k >= 10:
        current_number = str(k)
    else:
        current_number = '0' + str(k)
    current_VTK = readVTK('../data.00' + current_number + '.vtk', geometry='spherical')
    list_VTK.append(current_VTK)

list_VX = list()
list_RHO = list()
list_PRS = list()
list_K = list()
list_T = list()
for VTK in list_VTK:
    list_VX.append(VTK.data['VX1'])
    list_RHO.append(VTK.data['RHO'])
    list_PRS.append(VTK.data['PRS'])
    list_K.append(VTK.data['PRS']/(VTK.data['RHO']**gamma))
    list_T.append(VTK.data['PRS']/VTK.data['RHO'])
VTK = list_VTK[-1]
R = VTK.r
TH = VTK.theta
PHI = VTK.phi
ir = len(R)//5
ith = len(TH)//14
iphi = len(PHI)*5//7

rho0 = 1.
c = 1./(gamma - 1.)
D = kappa/(rho0*c)

def analytic_sol(r, th, phi, t):
    lambd2 = 1./r**2/np.sin(th)**2
    return amplitude*((3/r**2 - 1)*np.sin(r)/r - 3*np.cos(r)/r**2)*(-3*np.cos(th)*np.sin(th))*(np.sin(phi))*np.exp(-t*D*lambd2) + 1

def analytic_rad_sol(r, t):
    th = TH[ith]
    phi = PHI[iphi]
    return analytic_sol(r,th,phi,t)

def analytic_th_sol(th, t):
    r = R[ir]
    phi = PHI[iphi]
    return analytic_sol(r,th,phi,t)


def analytic_phi_sol(phi, t):
    r = R[ir]
    th = TH[ith]
    return analytic_sol(r,th,phi,t)


success = True
eps = 6.4e-8
for it,t in enumerate(list_time):
    analytic_T = np.zeros((R.shape[0], TH.shape[0], PHI.shape[0]))
    for ir,r in enumerate(R):
        for ith,th in enumerate(TH):
            for iphi,phi in enumerate(PHI):
                analytic_T[ir,ith,iphi] = analytic_sol(r,th,phi,t)
    TEMP = list_T[it]
    if np.mean(np.fabs(TEMP - analytic_T)) > eps:
      success = False

if success:
    print("SUCCESS")
    print("Error: {0}".format(np.mean(np.fabs(TEMP - analytic_T))))
    sys.exit(0)
else:
    print("Failed")
    print("Error: {0}".format(np.mean(np.fabs(TEMP - analytic_T))))
    sys.exit(1)
