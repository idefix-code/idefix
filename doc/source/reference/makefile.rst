Code configuration
==================
Configuring with cmake
----------------------

Basic usage
+++++++++++
`Cmake <https://cmake.org>`_ is a tool to control code generation on diverse platforms. It is the default tool used to configure *Idefix*. *Idefix* (and Kokkos)
requires ``cmake`` version >= 3.16. It is also recommended to use the graphical frontend ``ccmake`` to configure *Idefix*, as it allows one to have a rapid
overview of all of the configuration options and switch them according to the target architecture.

To configure *Idefix* with ``Cmake``, simply launch ``cmake $IDEFIX_DIR`` with the desired options **in a problem directory** (that is a directory containing at least ``definitions.hpp`` and ``setup.cpp``).
Alternatively, you can replace ``cmake`` by ``ccmake`` to get a more user-friendly graphical interface).


.. _configurationOptions:

Main configuration options
++++++++++++++++++++++++++

Several options can be enabled from the command line (or are accessible with ``ccmake`` GUI):

``-LH``
    List all of the available configure options and exit

``-D Idefix_MHD=ON``
    Enable MHD in the code (default disabled)

``-D Idefix_ENABLE_MPI=ON``
    Enable MPI parallelisation. Requires an MPI library. When used in conjonction with CUDA (Nvidia GPUs), a CUDA-aware MPI library is required by *Idefix*.

``-D Idefix_DEFS=foo.hpp``
    Specify a particular filename to be used in place of the default problem file ``definitions.hpp``

``-D Idefix_RECONSTRUCTION=x``
    Specify the type of reconstruction scheme (replaces the old "ORDER" parameter in ``definitions.hpp``). Accepted values for ``x`` are:
      + ``Constant``: first order, donor cell reconstruction,
      + ``Linear``: second order, piecewise linear reconstruction (PLM) using Van-leer slope limiter.
      + ``LimO3``: third order, Cada \& Torrilhon 2009
      + ``Parabolic``: fourth order piecewise parabolic reconstruction (PPM, Colella \& Woodward 1984)

.. note::

    The number of ghost cells is automatically adjusted as a function of the order of the reconstruction scheme.
    *Idefix* uses 2 ghost cells when ``ORDER < 4`` and 3 ghost cells when ``ORDER = 4``

``-D Kokkos_ENABLE_OPENMP=ON``
    Enable OpenMP parallelisation on supported compilers. Note that this can be enabled simultaneously with MPI, resulting in a hybrid MPI+OpenMP compilation.

``-D Kokkos_ENABLE_CUDA=ON``
    Enable Nvidia Cuda (for GPU targets). When enabled, ``cmake`` will attempt to auto-detect the target GPU architecture. If this fails, one needs to specify
    the target architecture adding ``-DKokkos_ARCH_{..}=ON`` (see below).

``-D Kokkos_ARCH_{...}=ON``
    Enable architecture-specific optimisation. A complete list can be obtained with the ``-LH`` option. Note that several host and target architecture can be enabled
    simulatenously (e.g for a CPU and a GPU). For instance:
      + Intel CPUs: BDW (Broadwell), HSW (Haswell), KNL (Knights Landing Xeon phi), SKX (Skylake Xeon with AVX512), SNB (Sandy/Ivy bridge), WSM (Westmere) ...
      + NVIDIA GPUs: PASCAL60, PASCAL61, VOLTA70, VOLTA72, AMPERE80, AMPERE86, ...
      + IBM CPUs: POWER7, POWER8, POWER9, BOOL (Blue gene Q)
      + ARM CPUs: ARMV80, ARMV81, ARMV8_THUNDERX, ARMV8_THUNDERX2, A64FX...
      + AMD CPUs: AMDAVX, ZEN, ZEN2, ZEN3...
      + AMD GPUs: VEGA906, VEGA908...



``-D CMAKE_CXX_COMPILER=foo``
    Request a specific ``foo`` compiler for the compilation and link. Alternatively, it is also possible to export the ``CXX`` environement variable with a valid compiler name
    before calling ``cmake``.

.. tip::

    Note that ``cmake`` keeps a cache of the previous configuration performed in a particular problem directory. To reset the configuration and start from scratch,
    delete the file `CMakeCache.txt`.

.. warning::

    *Idefix* ``cmake`` configuration expects the build directory to be a problem directory (that is a directory containing at least ``definitions.hpp`` and ``setup.cpp``).
    Launching ``cmake`` from a problem directory ensures that ``cmake`` will use that directory as its build directory. Note that it is also possible to use the ``-B``
    option to explictely tell ``cmake`` a path to a build=*Idefix* problem directory.


Setup-specific options
++++++++++++++++++++++

Some physical setup might require some ``cmake`` options to be set to specific value (e.g. an MHD setup will surely require MHD to be enable).
To avoid mistakes, it is then recommended to enforce this choice by creating a custom ``CMakeLists.txt`` in your setup directory, and setting
explicitely the options as they are required, using the functions ``set_idefix_property`` (for string properties) and ``enable_idefix_propery``/
``disable_idefix_propery`` (for boolean properties), as in the example below:

.. code-block::
    :caption: CMakeLists.txt

    set_idefix_property(Idefix_RECONSTRUCTION LimO3)
    enable_idefix_property(Idefix_MHD)



.. _customSourceFiles:

Add/replace custom source files
+++++++++++++++++++++++++++++++

It is possible to add custom source files to be compiled and linked against *Idefix*. This can be useful
if your setup requires complex functions and analysis in separate source files. To do so, add a ``CMakeLists.txt`` in your
problem directory, which adds to the ``idefix`` target  *all* the additional source files (i.e cpp *and* hpp headers). For instance,
say you want to add source files for an analysis, your ``CMakeLists.txt`` should look like:

.. code-block::
    :caption: CMakeLists.txt

    add_idefix_source(analysis.cpp)
    add_idefix_source(analysis.hpp)


*Idefix* also allows one to replace a source file in `$IDEFIX_DIR` by your own implementation. This is useful when developping new functionnalities without touching
the main directory of your *Idefix* repository. For instance, say one wants to replace the implementation of viscosity in `$IDEFIX_SRC/src/hydro/viscosity.cpp`,
with a customised `myviscosity.cpp` in the problem directory, one should add a ``CMakeLists.txt`` in the problem directory reading

.. code-block::
    :caption: CMakeLists.txt

    replace_idefix_source(hydro/viscosity.cpp myviscosity.cpp)


Note that the first parameter of ``replace_idefix_source`` is used as a search pattern in `$IDEFIX_DIR`. Hence it is possible to ommit the parent directory
of the file being replaced if there is only one file with that name in the *Idefix* source directory, which is not guaranteed (some classes may implement
methods with the same name). It is therefore recommended to add the parent directory in the first argument of ``replace_idefix_source``.


.. tip::

    Don't forget to delete `CMakeCache.txt` before attempting to reconfigure the code when adding a problem-specific
    ``CmakeLists.txt``.

Using GNU makefile and python configuration script (deprecated)
---------------------------------------------------------------
.. warning::

  Using the ``configure.py`` is deprecated and will be removed in a future version of *Idefix*. In particular, new target architectures
  will *not* be added to ``configure.py``. Use the ``cmake`` procedure instead.


The configure script
++++++++++++++++++++

Because the code can be configured for many architectures, it relies on a Python configuration script ``$IDEFIX_DIR/configure.py`` to generate the makefile needed. This script accepts
many options to adapt the generated makefile to the architecture on which one wants to run. A complete list of options can be obtained by running ``$IDEFIX_DIR/configure.py -h``. These options are:

``-h, --help``
    Display the help message and exit
``-mhd``
    Enable MHD in the code
``-arch=xxx``
    Compile for a specific CPU or GPU target. These corresponds to Kokkos target, so user can report to Kokkos documentation to get an up-to-date list of targets. At the time of writing, valid options are
     + Intel CPUs:    KNC, KNL, SNB, HSW, BDW, SKX
     + NVIDIA GPUs :  Kepler, Kepler30, Kepler32, Kepler35, Kepler37, Maxwell, Maxwell50, Maxwell52, Maxwell53, Pascal60, Pascal61, Volta70, Volta72, Turing75, Ampere80
     + ARM CPUS:      ARMv80, ARMv81, ARMv8-ThunderX, ARMv8-TX2
     + IBM:      BGQ, Power7, Power8, Power9
     + AMD-GPUS: Vega900, Vega906
     + AMD-CPUS: AMDAVX, Zen, Zen2
``-cxx=xxx``
    compile the code with the ``xxx`` C++ compiler. This option is ignored in GPU mode.
``-openmp``
    Enable OpenMP parallelisation on supported compilers (not available on GPUs for obvious reasons).
``-mpi``
    Enable MPI (message passing interface) when available. Note that this option is supported with CPU and GPU architectures as well, though GPUs require a CUDA-aware installation of MPI, such as OpenMPI.
``-defs=filename``
    Specify a particular ``filename`` to be used in place of the default ``definitions.hpp``

  .. tip::
    Note that when a source file in the ``makefile`` directory has the same filename as one of the original source file of your *Idefix* distribution, then
    ``make`` will compile the former in place of the original source file. This allows one to easily test a modification of your *Idefix* distribution
    by copying the original file and making your modification in your workdir.


Persistent configuration options
++++++++++++++++++++++++++++++++

System architecture (``-arch``) and custom compiler (``-cxx``) options can be
saved to a ``idefix.cfg`` file. Such a file can be stored locally, i.e. in the
directory of the physics problem, or globally in ``$HOME/.config/``[#]_ (or
``C:\Users\%USERNAME%\AppData`` on Windows). If both files exist, the global one is
ignored.

Here's an example ``idefix.cfg`` configuration file

.. code-block::

    [compilation]
    GPU = Volta70
    CPU = HSW
    CXX = icx

None of the parameters, or the configuration file itself, are mandatory.
Command line arguments take priority over options stored in ``idefix.cfg``.

.. [#] On POSIX systems, we follow `the XDG specification
<https://specifications.freedesktop.org/basedir-spec/basedir-spec-latest.html>`_,
and use ``$XDG_CONFIG_HOME``. On Windows, we use ``%APPDATA%`` instead.
