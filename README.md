[![Documentation Status](https://readthedocs.org/projects/idefix/badge/?version=latest)](https://idefix.readthedocs.io/en/latest/?badge=latest)
[![Idefix CIs](https://github.com/idefix-code/idefix/actions/workflows/idefix-ci.yml/badge.svg?branch=master)](https://github.com/idefix-code/idefix/actions/workflows/idefix-ci.yml)

<!-- toc -->

- [What is Idefix?](#what-is-idefix)
- [Documentation](#documentation)
- [Download:](#download)
- [Installation:](#installation)
- [Compile an example:](#compile-an-example)
- [Running](#running)
  * [serial (gpu/cpu), openMP (cpu)](#serial-gpucpu-openmp-cpu)
  * [With MPI (cpu)](#with-mpi-cpu)
  * [With MPI (gpu)](#with-mpi-gpu)
- [Profiling](#profiling)
- [Debugging](#debugging)
- [Code Validation](#code-validation)
- [Contributing](#contributing)

<!-- tocstop -->

What is Idefix?
---------------
Idefix is a computational fluid dynamics code based on a finite-volume high-order Godunov method, originally designed for astrophysical fluid dynamics applications.  Idefix is designed to be performance-portable, and uses the [Kokkos](https://github.com/kokkos/kokkos) framework to achieve this goal. This means that it can run both on your laptop's cpu and on the largest GPU Exascale clusters. More technically, Idefix can run in serial, use OpenMP and/or MPI (message passing interface) for parallelization, and use GPU acceleration when available (based on Nvidia Cuda, AMD HIP, etc...). All these capabilities are embedded within one single code, so the code relies on relatively abstracted classes and objects available in C++17, which are not necessarily
familiar to astrophysicists. A large effort has been devoted to simplify this level of abstraction so that the code can be modified by researchers and students familiar with C and who are aware of basic object-oriented concepts.


Idefix currently supports the following physics:

* Compressible hydrodynamics in 1D, 2D, 3D
* Compressible magnetohydrodynamics using constrained transport in 1D, 2D, 3D
* Multiple geometry (cartesian, polar, spherical)
* Variable mesh spacing
* Multiple parallelisation strategies (OpenMP, MPI, GPU offloading, etc...)
* Full non-ideal MHD (Ohmic, ambipolar, Hall)
* Viscosity and thermal diffusion
* Super-timestepping for all parabolic terms
* Orbital advection (Fargo-like)
* Self-gravity
* Multi dust species modelled as pressureless fluids
* Multiple planets interraction

Documentation
-------------

A full online documentation is available on [readTheDoc](https://idefix.readthedocs.io/latest/).


Download:
---------

Assuming you want to use https to get idefix (easiest option):

```shell
git clone --recurse-submodules https://github.com/idefix-code/idefix.git idefix
cd idefix
```

This will create and deploy Idefix in the directory `idefix`.

Installation:
-------------

Set the `IDEFIX_DIR` environment variable to the absolute path of the directory

```shell
export IDEFIX_DIR=<idefix main folder>
```

Add this line to `~/.<shell_rc_file>` for a permanent install.


Compile an example:
-------------------
Go to the example directory.
For instance:

```shell
cd test/HD/sod
```

Configure the code launching cmake (version >= 3.16) in the example directory:

```shell
cmake $IDEFIX_DIR
```

Several options can be enabled from the command line (a complete list is available with `cmake $IDEFIX_DIR -LH`). For instance: `-DIdefix_RECONSTRUCTION=Parabolic` (enable PPM reconstruction), `-DIdefix_MPI=ON` (enable mpi), `-DKokkos_ENABLE_OPENMP=ON` (enable openmp parallelisation), etc... For more complex target architectures, it is recommended to use cmake GUI launching `ccmake $IDEFIX_DIR` in place of `cmake` and then switching on the required options. See the [online documentation](https://idefix.readthedocs.io/latest/) for details.


One can then compile the code:

```shell
make -j8
```

Running
-------------------
### serial (gpu/cpu), openMP (cpu)
simple launch the executable

```shell
./idefix
```

### With MPI (cpu)
`-dec` can be used to specify a domain decomposition manually.

It can be omitted for 1D problems, or if `NX`, `NY`, `NZ` and `nproc` are **all** powers of 2.
Otherwise, `-dec` is mandatory. For instance, in 2D, using a 2x2 domain decomposition:

```shell
mpirun -np 4 ./idefix -dec 2 2
```

or in 3D, using a 1x2x4 decomposition:

```shell
mpirun -np 8 ./idefix -dec 1 2 4
```

### With MPI (gpu)
The same rules for cpu domain decomposition applies for gpus. In addition, one should manually specify how many GPU devices one wants to use **per node**. Example, in a run with 2 nodes, 4 gpu per node, one would launch idefix with

```shell
mpirun -np 8 ./idefix -dec 1 2 4 --kokkos-num-devices=4
```

Profiling
-------------------
use the embedded profiling tool by adding "-profile" when calling idefix (no need to recompile)

```shell
./idefix -profile
```

Debugging
-------------------
Add `-DIdefix_DEBUG=ON` when calling cmake, or activate the `Idefix_DEBUG` option in ccmake, and recompile.
Note that this option triggers a lot of outputs and memory access checks which significantly slow down the code.

Code Validation
---------------

Most of tests provided in the `test/` directory can be validated against analytical solution (standard test)
and/or pre-computed solutions (non-regression tests). Note that the validation relies on large reference
files that are stored in the separate `idefix-code/reference` repository that is cloned as a submodule.

Ensure that reference files
were properly downloaded (in the reference/ directory of the root of idefix) before attempting to validate the code.

In order to do a full validation of a particular test
(with all of the possible combination of algorithms), use the script `testme.py`
with the `-all` option, as in e.g.:
```shell
cd $IDEFIX_DIR/test/HD/sod
./testme.py -all
```

Tests require Python 3 along with some third party dependencies to be installed.
To install those deps, run
```shell
pip install -r test/python_requirements.txt
```

Contributing
-------------------
Idefix is developed with the help of the [pre-commit](https://pre-commit.com) framework.
We use [cpplint](https://en.wikipedia.org/wiki/Cpplint) to validate code style, mostly
following the Google standards for C++, and several pre-commit hooks to automatically fix
some coding bad practices.

It is recommended (though not mandatory) to install pre-commit by running the following
script from the top level of the repo
```shell
python3 -m pip install pre-commit
pre-commit install
```

Then, as one checks in their contribution with `git commit`, pre-commit hooks may perform
changes in situ. One then needs to re-add and enter the `git commit` command again for the
commit to be validated.
Note that an important hook that does _not_ perform auto-fixes is `cpplint`, so contributors
need to accomodate for this one by hand.

> Note that if for any reason you do not wish, or are unable to install pre-commit in your
> environment, formatting errors will be caught by our CI after you open a merge-request.
