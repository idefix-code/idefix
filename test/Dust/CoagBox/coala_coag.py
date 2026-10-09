from coala_py.src import *
import numpy as np
import os
import utils
import inifix

import argparse
parser = argparse.ArgumentParser()
parser.add_argument(
    "-ini", required=True, help="Idefix inifile in use for the run", type=str
)
args, unknown = parser.parse_known_args()

input_file = args.ini
c = utils.SetupParams(inifix.load(input_file))

#========== Input parameters required for Coala ==========#
nbins = c.nSpecies
massmin = c.massmin/(c.density*c.length**3)
massmax = c.massmax/(c.density*c.length**3)
kernel = c.kernel
K0 = 1.0
kpol = c.kpol
Q = c.Q
eps = c.eps
coeff_CFL = c.coagCFL
dthydro = c.fixed_dt
ndthydro = int(c.tstop/c.fixed_dt)
#=========================================================#

# init grid
massbins, massgrid = utils.Compute_dustgrid(nbins, massmin, massmax)

# Run coala to solve coagulation equation
if kernel != 3:
    gij_init, gij, time_coag = iterate_coag(kernel,K0,nbins,kpol,dthydro,ndthydro,coeff_CFL,Q,eps,massgrid,massbins)
else:
    _, _, _, dv = utils.Compute_test_velocity(nbins)
    gij_init, gij, time_coag = iterate_coag_kdv(kernel,K0,nbins,kpol,dthydro,ndthydro,coeff_CFL,Q,eps,massgrid,massbins,dv)
    
# save data
path_data = './data'
os.makedirs(path_data, exist_ok=True)

np.savetxt(path_data+"/test-coala-massbins_kernel-%d.txt"%(kernel),massbins)
np.savetxt(path_data+"/test-coala-gij_init_kernel-%d.txt"%(kernel),gij_init)
np.savetxt(path_data+"/test-coala-gij_end_kernel-%d.txt"%(kernel),gij)
print(f"Have generated COALA outputs in {path_data}")





