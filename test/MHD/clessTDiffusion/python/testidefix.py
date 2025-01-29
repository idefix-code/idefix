import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.vtk_io import readVTK
import numpy as np
import inifix

conf = inifix.load("../idefix.ini")
gamma = conf["Hydro"]["gamma"]
kappa = conf["Hydro"]["bragTDiffusion"][-1]

current_VTK = readVTK('../data.0003.vtk', geometry='spherical')

rho = current_VTK.data['RHO'].squeeze()
prs = current_VTK.data['PRS'].squeeze()
r   = current_VTK.r

idx1=np.where(r > 20)[0][0]
idx2=np.where(r > 30)[0][0]

def gamma_prime(gamma, beta):
    return (gamma+beta*(gamma-1))/(1+beta*(gamma-1))

success = True
eps = 3e-4

gamma_eff = np.gradient(prs, r)/np.gradient(rho, r)*rho/prs

Error=np.abs(gamma_eff[idx1:idx2]-gamma_prime(gamma, 1.5)).mean()
if Error > eps:
    success=False

if success:
    print("SUCCESS")
    print("Error: {0}".format(Error))
    sys.exit(0)
else:
    print("Failed")
    print("Error: {0}".format(Error))
    sys.exit(1)
