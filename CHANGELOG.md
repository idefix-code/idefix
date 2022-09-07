# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] 2022-09-07
### Changed
- use buffers for mpi axis exchanges to improve performances on GPUs (!195)
- slight optimisation of the cfl estimation for parabolic terms by using the maximum diffusion coefficient instead of the sum of all of the diffusion coefficients (makes a difference when several explicit parabolic terms are used simultaneously) (!176)
- use input::Get<T> and input::GetOrSet<T> instead of the old input:GetInt, input:GetReal... the new functions have a better error handling, and also allows explicit default values. (!208, !179)
- ensure that error messages are sent to std::cerr using a dedicated stream (!179)
- added a parameter check_nan to control the periodicity of Nan checks in the time integration (!191)
- added exception handling for the time integration, which saves a final vtk when a error is detected in the time integration loop (!190)
- fixed a bug in the test shearing box setup which led to memory corruption and incorrect pressure when ISOTHERMAL approximation was disabled (!190)
- make RKL faster when running in 3D with MHD diffusion terms by skipping the evolution of cell-centered fields (!215)
- fixed the many warning messages when compiling on CUDA (!229)
- improved spherical axis regularisation in full 3D (!245)
- improved PPM scheme by using Peterson & Hammet (2013) formulation (!251)

### Added
- single precision version is validated and fully operational. Can be enabled from cmake. (!197)
- isotropic thermal diffusion (anisotropic diffusion in MHD will come later) (!176)
- fixed a bug in 1D+1D in spherical geometry (!175)
- fixed a bug in calcCurrent which led to incorrectly computed currents on non-uniform grids (!206)
- Increase the efficiency of abort checks using MPI_Bcast instead of MPI_Allreduce (!174)
- allow the user to integrate the magnetic vector potential instead of the field to reduce the accumulation of roundoff errors on div(B) (experimental feature, can be enabled at config time) (!177)
- allow for domain decomposition along X3 with axis boundary condition and 2pi azimuthal domains (!182).
- addition of stateContainers to automatically evolve variables in multi-step time integration (!191)
- added a class to use easily lookup tables in idefix_for constructs (class LookupTable) (!213, !198, !178)
- added a shock flattening module (!219)
- added a grid coarsening module to increase the explicit timestep in heterogeneous grids (!246)
- added the linear wave test of Gardiner & Stone (2005) (!251)

## [1.0.0] 2022-01-13
### Changed
- enforce positivity of the Limo3 reconstruction scheme for density and pressure (when applicable) by reverting to second order in extreme cases. This makes LimO3 more stable.
- fixed a bug which resulted in a failure at detecting NaNs on some GPU architectures
- Fargo module has been moved to a class that belongs to the datablock, not to hydro (this is for future applications with dust+fargo). Use the [Fargo] block in your input file to define properties of the fargo module.
- gravity is now handled in a specific class, so that new gravity modules (e.g. self-gravity) can be handled automatically. Use the [Gravity] block in the input file to define properties of the gravity class (this includes user-defined potential and bodyforce).
- Nan detection is now explicit on all MPI processes
- fixed a bug which resulted in the generation of output files at each timestep when the output frequency was reduced at a restart.
- fixed a bug in VTK outputs which produced wrong grids in 1D spherical geometry.
- fixed a bug in VTK and dump outputs with MPI which resulted in garbage at the end of some files when vtk and dmp were overwritten.
- fixed a bug in Datablock initialisation which could lead to memory corruption
- fixed a bug in axis regularisation which could lead to memory corruption.
- reconstruction is now set by cmake and not in definitions.hpp (ORDER parameter). For backward compatibility, if definitions.hpp sets an ORDER, it supersedes the user choice in cmake.
- the examples in the ``test`` directory that can't be validated by the standard CI test (because of the lack of a quantitative validation test) are now handled separatly in an `Examples` queue: they are only compiled and run for a few cycles, looking for errors.
- refactor the MPI class so that it is more general than just for hydro objects.

### Added
- piecewise parabolic reconstruction (PPM)
- coding style guidelines in the documentation
- it is now possible to automatically enable cmake options (e.g. MHD) in each problem directory using set_idefix_property and enable_idefix_property in the problem's CMakeLists.txt.
- the Gravity class now handles automatically central potential wihtout needing to define your own user-defined potential. Use ``potential   central`` in the [[Gravity]] block of your input file.
- new -nolog, -nowrite, -maxcycles and -Werror command line arguments. Check the documentation for their usage.
- new ``Idefix_DEBUG`` cmake option to trigger debugging features (live call stack+array bound checks).
- Fargo now supports domain decomposition in the azimuthal direction.

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
