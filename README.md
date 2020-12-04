Download:
---------
    git clone
    git submodule init
    git submodule update

Compile an example:
-------------------
Go to the example directory:
for exmaple : 
    
`cd Test/HD/sod`

set the `IDEFIX_DIR` environment variable
    
`export IDEFIX_DIR=(idefix main folder)`

Configure the code:

`python3 $IDEFIX_DIR/configure.py`

Several options can be enabled (complete list can be accessed with the help `-h` option). For instance: `-mhd` (enable MHD, required in MHD tests), `-mpi` (enable mpi), `-openmp` (enable openmp parallelisation) `-gpu` (use gpu in place of cpu), etc...

One can then compile the code:

`make clean; make -j8`

Running
-------------------
### serial (gpu/cpu), openMP (cpu)
simple launch the executable

`./idefix`

### With MPI (cpu)
When `NX`, `NY`, `NZ` and `nproc` are **all** powers of 2, Idefix can guess the best domain
decomposition automatically. In the other cases, the -dec option for idefix is mandatory. Exemple, in 2D, using a 2x2 domain decomposition:

`mpirun -np 4 ./idefix -dec 2 2`

in 3D, using a 1x2x4 decomposition:

`mpirun -np 8 ./idefix -dec 1 2 4`

### With MPI (gpu)
The same rules for cpu domain decomposition applies for gpus. In addition, one should manually specify how many GPU devices one wants to use **per node**. Example, in a run with 2 nodes, 4 gpu per node, one would launch idefix with

`mpirun -np 8 ./idefix -dec 1 2 4 --kokkos-ndevices=4`

Profiling
-------------------
use the kokkos profiler tool, which can be downloaded from kokkos-tools (on github)
Then set the environement variable:

`export KOKKOS_PROFILE_LIBRARY=~/Kokkos/kokkos-tools/src/tools/space-time-stack/kp_space_time_stack.so`

and then simply run the code (no need to recompile)

Debugging
-------------------
Add `#define DEBUG` in definitions.hpp and recompile. More to come...

Running tests
-------------------
Tests require Python 3 along with some third party dependencies to be installed.
To install those deps, run
```
pip install -r test/python_requirements.txt
```

The test suite itself is then run with
```
bash test/checks.sh
```
