The makefile
=====================

Because the code can be configured for many architectures, it relies on a Python configuration script ``$IDEFIX_DIR/configure.py`` to generate the makefile needed. This script accepts
many options to adapt the generated makefile to the architecture on which one wants to run. A complete list of options can be obtained by running ``$IDEFIX_DIR/configure.py -h``. These options are:

``-h, --help``
    Display the help message and exit
``-mhd``
    Enable MHD in the code
``-gpu``
    Enable GPU target architecture (otherwise default target is cpu).
``-arch=xxx``
    Compile for a specific CPU or GPU target. These corresponds to Kokkos target, so user can report to Kokkos documentation to get an up-to-date list of targets. At the time of writing, valid options are
     + Intel CPUs:    KNC, KNL, SNB, HSW, BDW, SKX
     + NVIDIA GPUs :  Kepler, Kepler30, Kepler32, Kepler35, Kepler37, Maxwell, Maxwell50, Maxwell52, Maxwell53, Pascal60, Pascal61, Volta70, Volta72, Turing75, Ampere80
     + ARM CPUS:      ARMv80, ARMv81, ARMv8-ThunderX, ARMv8-TX2
     + IBM:      BGQ, Power7, Power8, Power9
     + AMD-GPUS: Vega900, Vega906
     + AMD-CPUS: AMDAVX, Zen, Zen2
``-openmp``
    Enable OpenMP parallelisation on supported compilers (not available on GPUs for obvious reasons).

``-mpi``
    Enable MPI (message passing interface) when available. Note that this option can be used in conjonction with -gpu to run idefix simultaneously on several GPUs. This feature requires a CUDA-aware installation of MPI, such as OpenMPI.

  .. tip::
    Note that when a source file in the ``makefile`` directory has the same filename as one of the original source file of your *Idefix* distribution, then
    ``make`` will compile the former in place of the original source file. This allows one to easily test a modification of your *Idefix* distribution
    by copying the original file and making your modification in your workdir.
