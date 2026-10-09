#from coala_py.src import *                 # Only if we need to call COALA, in function Compute_coag_precalc()
#from numba_progress import ProgressBar     # Only if we need to call COALA, in function Compute_coag_precalc()
import numpy as np

class SetupParams:
    """
    Get the parameters from the .ini file 
    If it is not found in the .ini file, parameter is set to None
    """ 
    def __init__(self, conf):
        self.tstop        = self._get(conf, 'TimeIntegrator', 'tstop', cast=float)
        self.fixed_dt     = self._get(conf, 'TimeIntegrator', 'fixed_dt', cast=float)
        self.nstages      = self._get(conf, 'TimeIntegrator', 'nstages', cast=int)
        self.nSpecies     = self._get(conf, 'Dust', 'nSpecies', cast=int)
        self.size         = self._get(conf, 'Dust', 'drag', cast=lambda v: [float(x) for x in v[1:]])
        self.coag         = self._get(conf, 'Coala', 'coag', cast=str)
        self.massmin      = self._get(conf, 'Coala', 'massmin', cast=float)
        self.massmax      = self._get(conf, 'Coala', 'massmax', cast=float)
        self.kernel       = self._get(conf, 'Coala', 'kernel', cast=int)
        self.kpol         = self._get(conf, 'Coala', 'kpol', cast=int)
        self.Q            = self._get(conf, 'Coala', 'Q', cast=int)
        self.eps          = self._get(conf, 'Coala', 'eps', cast=float)
        self.coagCFL      = self._get(conf, 'Coala', 'coagCFL', cast=float)
        self.rho0         = self._get(conf, 'Setup', 'rho0')
        self.dtg_ratio    = self._get(conf, 'Setup', 'dtg_ratio')
        self.betaSizeMax  = self._get(conf, 'Setup', 'betaSizeMax')
        self.length       = self._get(conf, 'Units', 'length')
        self.density      = self._get(conf, 'Units', 'density')
        self.velocity     = self._get(conf, 'Units', 'velocity')
        
    def _get(self, conf, block_name, name, default=None, cast=float):
        if block_name not in conf or name not in conf[block_name]:
            return default
        try:
            return cast(conf[block_name][name])
        except Exception:
            return default
        
def Compute_dustgrid(nbins, mini, maxi):
    '''
    Compute a grid for mass or size in a same way than Coala
    '''
    r = (maxi/mini)**(1/nbins)
    bins = np.zeros(nbins)
    grid = np.zeros(nbins+1)
    grid[0] = mini
    
    for j in range(nbins):
        grid[j+1] = r*grid[j]
        bins[j] = 0.5*(grid[j]+grid[j+1])
        
    return bins, grid

def Compute_coag_precalc(kernel, K0, Q, nbins, kpol, massgrid):
    """
    Precompute coagulation flux tables
    RETURN
    tensor_tabflux_coag : ndarray
    tensor_tabintflux_coag : ndarray or None
    """

    vecnodes, vecweights = np.polynomial.legendre.leggauss(Q)
    mat_coeffs_leg = legendre_coeffs(kpol)

    if kpol == 0:
        tensor_tabflux_coag = np.zeros(
            (nbins, nbins, nbins)
        )
        with ProgressBar(total=nbins) as progress:
            compute_coagtabflux_k0_numba(
                kernel, K0, Q,
                vecnodes, vecweights,
                nbins, massgrid,
                mat_coeffs_leg,
                tensor_tabflux_coag,
                progress
            )
        tensor_tabintflux_coag = None
    else:
        tensor_tabflux_coag = np.zeros(
            (nbins, nbins, nbins, kpol+1, kpol+1)
        )
        with ProgressBar(total=nbins) as progress:
            compute_coagtabflux_numba(
                kernel, K0, Q,
                vecnodes, vecweights,
                nbins, kpol, massgrid,
                mat_coeffs_leg,
                tensor_tabflux_coag,
                progress
            )
        tensor_tabintflux_coag = np.zeros(
            (nbins, kpol+1, nbins, nbins, kpol+1, kpol+1)
        )
        with ProgressBar(total=nbins) as progress:
            compute_coagtabintflux_numba(
                kernel, K0, Q,
                vecnodes, vecweights,
                nbins, kpol, massgrid,
                mat_coeffs_leg,
                tensor_tabintflux_coag,
                progress
            )
    return tensor_tabflux_coag, tensor_tabintflux_coag

def Compute_MRN_mass(grid, q, mini, maxi, eps):
    """
    Initialize the MRN distribution from the bin mini to the bin maxi, in mass
    Normalized to 1
    """
    argmin = np.argmin(np.abs(grid - mini))
    argmax = np.argmin(np.abs(grid - maxi))

    p = (4-q)/3
    norm = 1/(grid[argmax]**p - grid[argmin]**p)
    gij = np.zeros(len(grid)-1)

    gij[argmin:argmax] = (grid[argmin+1:argmax+1]**p - grid[argmin:argmax]**p)/norm
    gij[gij<eps] = eps
    gij /= np.sum(gij)
    return gij

def Compute_mexp_mass(grid, mini, maxi, eps):
    """
    Initialize the distribution m*exp(-m) from the bin mini to the bin max, in mass
    Normalized to 1
    """
    argmin = np.argmin(np.abs(grid - mini))
    argmax = np.argmin(np.abs(grid - maxi))
    
    norm = 1/((1 + grid[argmin]) * np.exp(-grid[argmin]) - (1 + grid[argmax]) * np.exp(-grid[argmax]))
    gij = np.zeros(len(grid) - 1)
    
    gij[argmin:argmax] = ((1 + grid[argmin:argmax]) * np.exp(-grid[argmin:argmax]) - (1 + grid[argmin+1:argmax+1]) * np.exp(-grid[argmin+1:argmax+1])) * norm
    gij[gij < eps] = eps
    gij /= np.sum(gij)
    return gij

def Compute_Coala_mass(grid, mini, maxi, eps): 
    """
    The same initialisation than in Coala (for kpol = 0)
    """
    argmin = np.argmin(np.abs(grid - mini))
    argmax = np.argmin(np.abs(grid - maxi))
    
    gij = np.zeros(len(grid)-1)
    gij[argmin:argmax] = (
        (1 + grid[argmin:argmax]) * np.exp(-grid[argmin:argmax])
        - (1 + grid[argmin+1:argmax+1]) * np.exp(-grid[argmin+1:argmax+1])
        ) / (
        grid[argmin+1:argmax+1] - grid[argmin:argmax]
        )
    gij[gij < eps] = eps
    return gij

def Compute_test_velocity(nbins):
    vmin, vmax = 1e-2, 1e-1
    vx1 = np.linspace(vmin, vmax, nbins)
    vx2 = 1.3*vx1
    vx3 = 1.7*vx1
    dvx1 = vx1[:,None] - vx1[None,:]
    dvx2 = vx2[:,None] - vx2[None,:]
    dvx3 = vx3[:,None] - vx3[None,:]
    dv = np.sqrt(dvx1**2 + dvx2**2 + dvx3**2)
    return vx1, vx2, vx3, dv
