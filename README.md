Download:
---------
git clone
git submodule init
git submodule update

Compile an example:
-------------------
Go to the example directory:
for exmaple : cd Test/HD/sod 

export IDEFIX_DIR=(idefix main folder)

Configure the code:
python3 $IDEFIX_DIR/configure.py
can add -mpi (enable mpi), -openmp (enable openmp parallelisation) -gpu (use gpu in place of cpu), -mhd (enable mhd), etc...
make clean; make -j8

Running
-------------------
simple launch the executable
./idefix

With mpi:
Launch the executable and provide the domain decomposition. Exemple, in 2D, using a 2x2 domain decomposition
mpirun -np 4 ./idefix -dec 2 2
in 3D, using a 1x2x4 decomposition:
mpirunt -np 8 ./idefix -dec 1 2 4

When NX, NY, NZ and nproc are *all* powers of 2, the -dec option is not mandatory as Idefix can find the best domain
decomposition automatically.

In all the other cases, the -dec option is mandatory

With MPI+GPU
Specify that kokkos should use a different gpu for each MPI process
mpirun -np 8 ./idefix -dec 1 2 4 --kokkos-ndevices=8

Profiling
-------------------
use the kokkos profiler tool, which can be downloaded from kokkos-tools (on github)
set the environement variable:
export KOKKOS_PROFILE_LIBRARY=~/Kokkos/kokkos-tools/src/tools/space-time-stack/kp_space_time_stack.so 
and then simply run the code (no need to recompile)


