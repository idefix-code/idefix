from pydefix import *
import utils
import numpy as np
import inifix
import matplotlib.pyplot as plt

def initflow(data): 
    input_file = "idefix-constantKernel.ini"  # To modified : we need to collect directly the parameters of the input file via the input class of idefix 
    c = utils.SetupParams(inifix.load(input_file))
    rho0 = c.rho0/c.density
    nSpecies = c.nSpecies
    massmin = c.massmin/(c.density*c.length**3)
    massmax = c.massmax/(c.density*c.length**3)
       
    #CI for the gas
    data.Vc[RHO,:,:,:] = rho0
    data.Vc[VX1,:,:,:] = 0.0
    data.Vc[VX2,:,:,:] = 0.0
    data.Vc[VX3,:,:,:] = 0.0
    
    #CI for the dust
    if nSpecies is not None:
        #Same init than in the Coala test
        _, massgrid = utils.Compute_dustgrid(nSpecies, massmin, massmax)
        dust_mass_fraction_init = utils.Compute_Coala_mass(massgrid, massgrid[0], massgrid[-1], c.eps)
        h = massgrid[1:] - massgrid[:-1]
        dust_mass_fraction = dust_mass_fraction_init*h
        dust_mass_fraction /= np.sum(dust_mass_fraction)
        if c.kernel != 3:
            vx1, vx2, vx3 = np.zeros((nSpecies)), np.zeros((nSpecies)), np.zeros((nSpecies))
        else:
            vx1, vx2, vx3, _ = utils.Compute_test_velocity(nSpecies)   # WARNING: do not change values in utils.Compute_test_velocity() without uploading the COALA solution
        for n in range(nSpecies):
            data.dustVc[n][RHO,:,:,:] = dust_mass_fraction[n]
            data.dustVc[n][VX1,:,:,:] = vx1[n]
            data.dustVc[n][VX2,:,:,:] = vx2[n]
            data.dustVc[n][VX3,:,:,:] = vx3[n]