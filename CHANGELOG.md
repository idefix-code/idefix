# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Upcoming
### Changed
- enforce positivity of the Limo3 reconstruction scheme for density and pressure (when applicable) by reverting to second order in extreme cases. This makes LimO3 more stable.
- fixed a bug which resulted in a failure at detecting NaNs on some GPU architectures
- Fargo module has been moved to a class that belongs to the datablock, not to hydro (this is for future applications with dust+fargo). Use the [[Fargo]] block in your input file to define properties of the fargo module.
- gravity is now handled in a specific class, so that new gravity modules (e.g. self-gravity) can be handled automatically. Use the [[Gravity]] block in the input file to define properties of the gravity class (this includes user-defined potential and bodyforce).
- Nan detection is now explicit on all MPI processes
- fixed a bug which resulted in the generation of output files at each timestep when the output frequency was reduced at a restart.
- fixed a bug in VTK outputs which produced wrong grids in 1D spherical geometry.
- fixed a bug in VTK and dump outputs with MPI which resulted in garbage at the end of some files when vtk and dmp were overwritten.
- reconstruction is now set by cmake and not in definitions.hpp (ORDER parameter). For backward compatibility, if definitions.hpp sets an ORDER, it supersedes the user choice in cmake.
- the examples in the ``test`` directory that can't be validated by the standard CI test (because of the lack of a quantitative validation test) are now handled separatly in an `Examples` queue: they are only compiled and run for a few cycles, looking for errors.

### Added
- piecewise parabolic reconstruction (PPM)
- coding style guidelines in the documentation
- it is now possible to automatically enable cmake options (e.g. MHD) in each problem directory using set_idefix_property and enable_idefix_property in the problem's CMakeLists.txt.
- the Gravity class now handles automatically central potential wihtout needing to define your own user-defined potential. Use ``potential   central`` in the [[Gravity]] block of your input file.
- new -nolog, -nowrite, -maxcycles and -Werror command line arguments. Check the documentation for their usage.
- new ``Idefix_DEBUG`` cmake option to trigger debugging features (live call stack+array bound checks).

### Removed
- configure.py support and related functions.

## [0.9.1] 2021-10-27
### Changed
- fixed a bug in stretch grid, which led to incorrect grid spacing in s+ grids. This might break the restart of MHD runs from dumps created from previous versions.
- fixed a bug in boundary conditions which could lead to memory corruption when COMPONENTS < DIMENSIONS in MHD

## [0.9.0] 2021-10-19
### Added
- new `-autotune` runtime option, which tests and chooses the best loop unrolling strategy
- 3rd order reconstruction scheme (following Cada \& Torrilhon 2009). Can be enabled in definitions.hpp
- ppm reconstruction in Fargo. Automatically used with ORDER=3, otherwise can be enabled with Idefix_HIGH_ORDER_FARGO in cmake configuration
- new EMF averaging scheme using 2D HLLD Riemann solves
- VTK files now includes TIME, GEOMETRY and PERIODICITY fields which can be read with VTK_io python routines
- new python routines to read Idefix dump files
- new 3D Shearing box boundary conditions
- new command file to stop Idefix when running remotely
- new boundary loop wrappers in `boundaryloop.hpp`

### Changed
- code configuration with `cmake` instead of configure.py
- EMF averaging scheme is now set at run time, not compile time.
- VTK I/O python functions have been refactored: a single VTK_Read routine can now be called for all geometries
- `-restart` is ignored if no restart dump is found
- auto-detect HIP GPU offloading (used for AMD GPUs)

### Removed
- deprecated the `configure.py` script. Cmake should now be used instead (this version is the last one that still supports `configure.py`).

## [0.8.1] - 2021-06-24
### Changed
- Fixed a bug in the hydro HLL solver which used the wrong sound speed for flux computations.

## [0.8.0] - 2021-06-13
### Added
- New `DumpImage` class to load and use the data of dump files without restarting
- New Runge-Kutta-Legendre scheme to speed up parabolic term computation. Compatible with Viscosity, Ambipolar & Ohmic diffusions.

### Changed
- Optimisation: merge ExtrapolatePrimVar and Riemann solves
- Optimisation: Limit array accesses in nonIdeal MHD flux computations
- Optimisation: Improved VTK write speeds on non-cartesian geometries
- rotation now works as it should in polar & spherical coordinates (in this case, it includes both Coriolis and centrifugal acceleration)
- fix a bug in fargo which broke axisymmetric symmetry in some circumstances
- fix a bug when -restart was used without any number (use by default latest restart dump file)
- properly check dependency (uses -M option of the compiler). This generates a series of dependency files (.d) during compilation
- documentation includes on the fly doxygen generated API.

### Removed
- deprecate the `-gpu` option in `configure.py`. The GPU mode is now automatically activated if a GPU architecture is requested.
- add support for persistent options for `configure.py`, specifically for `-arch` and `-cxx`, in a user-defined `idefix.cfg` file.

## [0.7.2] - 2021-05-11
### Changed
- fixed a bug in MakeVsFromAmag which initialised a B from a user-provided vector potential with
non-zero div(B) in spherical coordinates in some circumstances

## [0.7.1] - 2021-04-26
### Changed
 - fixed a bug in the viscous flux computation in spherical coordinates (Thx F. Rincon)

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
