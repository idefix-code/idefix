The latest documentation is available
[here](https://discs.gricad-pages.univ-grenoble-alpes.fr/idefix/index.html).

<!-- toc -->

- [Download:](#download)
- [Installation:](#installation)
- [Compile an example:](#compile-an-example)
- [Running](#running)
  * [serial (gpu/cpu), openMP (cpu)](#serial-gpucpu-openmp-cpu)
  * [With MPI (cpu)](#with-mpi-cpu)
  * [With MPI (gpu)](#with-mpi-gpu)
- [Profiling](#profiling)
- [Debugging](#debugging)
- [Running tests](#running-tests)

<!-- tocstop -->

Download:
---------

using either ssh or https url as `<address>`
```shell
git clone <address>
cd idefix
git submodule init
git submodule update
```

Installation:
-------------

Set the `IDEFIX_DIR` environment variable to the absolute path of the directory

```shell
export IDEFIX_DIR=<idefix main folder>
```

Add this line to `~/.<shell_rc_file>` for a permanent install.


Compile an example:
-------------------
Go to the example directory:
for exmaple :

```shell
cd test/HD/sod
```

Configure the code:

```shell
python3 $IDEFIX_DIR/configure.py
```

Several options can be enabled (complete list can be accessed with the help `-h` option). For instance: `-mhd` (enable MHD, required in MHD tests), `-mpi` (enable mpi), `-openmp` (enable openmp parallelisation) `-gpu` (use gpu in place of cpu), etc...

One can then compile the code:

```shell
make clean; make -j8
```

Running
-------------------
### serial (gpu/cpu), openMP (cpu)
simple launch the executable

```shell
./idefix
```

### With MPI (cpu)
When `NX`, `NY`, `NZ` and `nproc` are **all** powers of 2, Idefix can guess the best domain
decomposition automatically. In the other cases, the -dec option for idefix is mandatory. Exemple, in 2D, using a 2x2 domain decomposition:

```shell
mpirun -np 4 ./idefix -dec 2 2
```

in 3D, using a 1x2x4 decomposition:

```shell
mpirun -np 8 ./idefix -dec 1 2 4
```

### With MPI (gpu)
The same rules for cpu domain decomposition applies for gpus. In addition, one should manually specify how many GPU devices one wants to use **per node**. Example, in a run with 2 nodes, 4 gpu per node, one would launch idefix with

```shell
mpirun -np 8 ./idefix -dec 1 2 4 --kokkos-ndevices=4
```

Profiling
-------------------
use the kokkos profiler tool, which can be downloaded from kokkos-tools (on github)
Then set the environement variable:

```shell
export KOKKOS_PROFILE_LIBRARY=~/Kokkos/kokkos-tools/src/tools/space-time-stack/kp_space_time_stack.so
````

and then simply run the code (no need to recompile)

Debugging
-------------------
Add `#define DEBUG` in definitions.hpp and recompile. More to come...

Running tests
-------------------
Tests require Python 3 along with some third party dependencies to be installed.
To install those deps, run
```shell
pip install -r test/python_requirements.txt
```

The test suite itself is then run with
```shell
bash test/checks.sh
```
