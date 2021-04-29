# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

### Upcoming
- Improved VTK write speeds on non-cartesian geometries
- New `DumpImage` class to load and use the data of dump files without restarting
- New Runge-Kutta-Legendre scheme to speed up parabolic term computation. Compatible with Viscosity, Ambipolar & Ohmic diffusions

## [0.7.1] - 2021-04-26
### Changed
- fixed a bug in the viscous flux computation in spherical coordinates (Thx F. Rincon).

## [0.7.0] - 2021-04-07
### Added
- New orbital advection scheme (aka Fargo) available in 2D and 3D HD/MHD. This implementation fixes the bug found in Pluto 4.4.
- New axis boundary condition allowing one to look at 3D problems up to the spherical axis. This option can also be used for axisymmetric problems or 3D problem with phi ranging on a fraction of 2pi.
- New 2D HLL Riemann solver for MHD evolution based on Londrillo & Del Zanna. Activated using ```#define EMF_AVERAGE  UCT_HLL```
- New python tools to read idefix debug dumps (not enabled by default) in ```pytools/idfx_io.py```
- New option to use a fixed time step during the computation. Use ```fixed_dt``` in your ini file.
- New regression tests on all of the possible EMF reconstruction schemes.
- New optional ```max_runtime``` parameter to stop the code when a given runtime is reached. This can be useful in situations where signaling is not available.

### Changed
- Improved completion logging message, including total memory used by each MPI process in each memory space (requires Kokkos>3.2), and fix performance indication overflow.
- The parameter ```first_dt``` now has a small default value, if not set in the ini file.
- EMF-related methods are now all included in the ElectroMotiveForce class.
- Refactored configure script.
- Fix a bug found in MHD+MPI which led to a slow deviation of adjacent BXs at domain boundaries. The bug fix implies an additional MPI call to synchronise EMFs at domain boundaries.
- Fix a minor bug in Hall-MHD fluxes in the 1D HLL Riemann solver
- Synchronise all of the MPI processes when a SIGUSR signal is received by one process to ensure that they all stop simultaneously.
- Fix a bug in periodic boundary conditions if only one cell is present in said direction.
- Fix an offset bug in the non-ideal EMF computation in 3D.
- Fix a bug in the python VTK Spherical reader leading to unconsistent phi coordinates when the axis is included in the simulation domain.
- Added a changelog

### Removed
- nothing
