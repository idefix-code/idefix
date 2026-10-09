import utils
import inifix
import numpy as np 
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.getenv("IDEFIX_DIR"))
from pytools.vtk_io import readVTK
import argparse
parser = argparse.ArgumentParser()
parser.add_argument(
    "-noplot", default=False, help="disable plotting", action="store_true"
)
parser.add_argument(
    "-ini", required=True, help="Idefix inifile in use for the run", type=str
)
args, unknown = parser.parse_known_args()

input_file = args.ini
c = utils.SetupParams(inifix.load(input_file))
data_path = "./data"

match c.kernel:
    case 0:
        error_val = 2.e-4
        kernel_name = "constant"
    case 1:
        error_val = 7.e-4
        kernel_name = "additive"
    case 2:
        error_val = 2.e-4
        kernel_name = "brownian"
    case 3:
        error_val = 6.e-5
        kernel_name = "physical"
    case _:
        print("ERROR: chose a valid kernel (can be 0, 1, 2 or 3)")
 
massbins, massgrid = utils.Compute_dustgrid(c.nSpecies, c.massmin/(c.density*c.length**3), c.massmax/(c.density*c.length**3))
h = massgrid[1:] - massgrid[:-1]

#========== IDEFIX PART ==========#
VTK0 = readVTK(data_path+'/data.0000.vtk')
VTKt = readVTK(data_path+'/data.0001.vtk')

rho_dust_0 = np.empty(c.nSpecies)
rho_dust_t = np.empty(c.nSpecies)
for n in range(c.nSpecies):
    rho_dust_0[n] = np.mean(VTK0.data[f'Dust{n}_RHO'])
    rho_dust_t[n] = np.mean(VTKt.data[f'Dust{n}_RHO'])  
    
#========== COALA PART ===========#
coala_massbins = np.loadtxt(data_path+"/test-coala-massbins_kernel-%d.txt"%(c.kernel))
gij_t0 = np.loadtxt(data_path+"/test-coala-gij_init_kernel-%d.txt"%(c.kernel))
gij_tend = np.loadtxt(data_path+"/test-coala-gij_end_kernel-%d.txt"%(c.kernel))

coala_rho_dust_0 = gij_t0*h
coala_rho_dust_t = gij_tend*h
      
# plot the test
if not args.noplot:
    plt.figure(figsize=(7,5))
    ax = plt.gca()
    #========== IDEFIX PART ==========#
    ax.plot(massbins, rho_dust_0, 'o-', color='skyblue', label=r'$t = t_{0}$')
    ax.plot(massbins, rho_dust_t, 'o-', color='blue', label=r'$t = t_{f}$', alpha=1)

    #========== COALA PART ===========#
    marker_style = dict(marker='o', markersize=8, markerfacecolor='none', linestyle='', markeredgewidth=2)

    plt.plot(coala_massbins, coala_rho_dust_0, markeredgecolor='green', **marker_style, alpha=0.7)
    plt.plot(coala_massbins, coala_rho_dust_t, markeredgecolor='green', label="coala", **marker_style, alpha=1)

    #========== END OF THE PLOT ======#
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('mass  [c.u]', fontsize=16)
    ax.set_ylabel(r'$\langle \rho_{\mathrm{d}} \rangle_{\mathrm{box}}$  [c.u]', fontsize=16)
    ax.minorticks_on()
    ax.tick_params(
        axis='both',
        which='both',
        direction='in',
        top=True,
        right=True,
    )
    ax.set_ylim(1e-15, 1e1)
    ax.grid(True, which='both', alpha=0.2)
    ax.legend(fontsize=12)
    ax.set_title(f'kernel = {c.kernel} ({kernel_name})', fontsize=16)
    plt.tight_layout()
    plt.ioff()
    plt.show()

error = 1/c.nSpecies*np.sum(np.abs(rho_dust_t - coala_rho_dust_t))
print("Error=%e" % error)
if error<error_val:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
