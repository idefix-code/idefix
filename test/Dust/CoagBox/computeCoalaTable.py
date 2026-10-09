import utils
import numpy as np
import inifix
import os
import sys

import argparse
parser = argparse.ArgumentParser()
parser.add_argument(
    "-ini", required=True, help="Idefix inifile in use for the run", type=str
)
args, unknown = parser.parse_known_args()

input_file = args.ini
c = utils.SetupParams(inifix.load(input_file))

if c.coag == 'True':
    # Collect the input parameters
    nSpecies = c.nSpecies
    kernel = c.kernel
    kpol = c.kpol
    Q = c.Q
    massmin = c.massmin/(c.density*c.length**3)
    massmax = c.massmax/(c.density*c.length**3)  
    print(f'\nHave dust coagulation with:\n  nSpecies = {nSpecies}\n  kernel = {kernel}\n  kpol = {kpol}\n  Q = {Q}')
    print(f'Massgrid and massbins are computed with:\n  massmin = {massmin}, smax = {massmax}\n')
    
    # output directory for coagulation tables
    dir_table = './data'
    os.makedirs(dir_table, exist_ok=True)
    
    # init grid 
    massbins, massgrid = utils.Compute_dustgrid(nSpecies, massmin, massmax)
    
    # compute the coagulation table
    coagTabFlux, coagTabintFlux = utils.Compute_coag_precalc(kernel, 1.0, Q, nSpecies, kpol, massgrid)
    print('Coagulation precalculation done!')
    fname_massgrid = 'test-massgrid-kernel-%d'%(kernel)
    fname_coagTabFlux = 'test-coagTabFlux-kernel-%d'%(kernel)  
    np.save(dir_table+'/'+fname_massgrid+'.npy', massgrid)
    np.save(dir_table+'/'+fname_coagTabFlux+'.npy', coagTabFlux)
    print('-->'+dir_table+'/'+fname_massgrid+'.npy saved')
    print('-->'+dir_table+'/'+fname_coagTabFlux+'.npy saved')
    
    # In our specific choice of units, rho = cs = 1.0 so the size drag
    # gives beta_size = tau_size the stopping time
    # Here we compute the beta_size bins to put in the idefix.ini file 
    # according to the massgrid and a beta_size_max (to give here) which
    # is the maximum stopping time in our distribution
    beta_size_max = c.betaSizeMax/(c.density*c.length)
    beta_size_bins = beta_size_max*(massbins/massbins[-1])**(1.0/3.0)
    print(f'\nParameter size to put in {input_file} (with betaSizeMax = {c.betaSizeMax}):\n {beta_size_bins*(c.density*c.length)}\n')
    
    if coagTabintFlux is not None:
        fname_coagTabIntFlux = 'coagIntTabFlux'
        np.save(dir_table+'/'+fname_coagTabIntFlux+'.npy', coagTabintFlux)
        print('-->'+dir_table+'/'+fname_coagTabIntFlux+'.npy saved')
        print('WARNING : Only kpol = 0 is available to use Coala in Idefix')
else:
    print('No dust coagulation')