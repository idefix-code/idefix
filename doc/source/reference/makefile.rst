Code configuration with Cmake
=============================


Introduction to Cmake
+++++++++++++++++++++
`Cmake <https://cmake.org>`_ is a tool to control code generation on diverse platforms. It is the default tool used to configure *Idefix*. *Idefix* (and Kokkos)
requires ``cmake`` version >= 3.16. It is also recommended to use the graphical frontend ``ccmake`` to configure *Idefix*, as it allows one to have a rapid
overview of all of the configuration options and switch them according to the target architecture.

To configure *Idefix* with ``Cmake``, simply launch ``cmake $IDEFIX_DIR`` with the desired options **in a problem directory** (that is a directory containing at least ``definitions.hpp`` and ``setup.cpp``).
Alternatively, you can replace ``cmake`` by ``ccmake`` to get a more user-friendly graphical interface).

.. warning::

  The old configuration script ``configure.py`` is not supported in this version of *Idefix* and will likely
  break compilation.


.. _configurationOptions:

Main configuration options
++++++++++++++++++++++++++

Several options can be enabled from the command line (or are accessible with ``ccmake`` GUI):

``-LH``
    List all of the available configure options and exit

``-D Idefix_MHD=ON``
    Enable MHD. If your setup is MHD-only, making this option mandatory, it is recommended to use a local ``CMakeLists.txt`` file instead, as described in :ref:`setupSpecificOptions`

``-D Idefix_MPI=ON``
    Enable MPI parallelisation. Requires an MPI library. When used in conjonction with CUDA (Nvidia GPUs), a CUDA-aware MPI library is required by *Idefix*.

``-D Idefix_DEFS=foo.hpp``
    Specify a particular filename to be used in place of the default problem file ``definitions.hpp``

``-D Idefix_DEBUG=ON``
    Enable debug options in *Idefix*. This triggers a lot of outputs, and automatic bound checks of array accesses. As a result, this
    option makes the code very slow.

``-D Idefix_HDF5=ON``
    Enable HDF5 outputs. Requires the HDF5 library on the target system. Required for *Idefix* XDMF outputs.

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


.. _setupExamples:

Configuration examples for selected clusters
++++++++++++++++++++++++++++++++++++++++++++


AdAstra at CINES, AMD Mi250X GPUs
---------------------------------

We recommend the following modules and environement variables on AdAstra:

.. code-block:: bash

    module load PrgEnv-cray-amd
    module load cray-mpich
    module load craype-network-ofi
    module load cce
    module load cpe
    module load rocm/5.2.0
    export LDFLAGS="-L${ROCM_PATH}/lib -lamdhip64 -lstdc++fs"

The last line being there to guarantee the link to the HIP library and the access to specific
C++17 <filesystem> functions.

Finally, *Idefix* can be configured to run on Mi250 by enabling HIP and the desired architecture with the following options to ccmake:
``-DKokkos_ENABLE_HIP=ON -DKokkos_ENABLE_HIP_MULTIPLE_KERNEL_INSTANTIATION=ON -DKokkos_ARCH_VEGA90A=ON``

MPI (multi-GPU) can be enabled by adding ``-DIdefix_MPI=ON`` as usual.

Jean Zay at IDRIS, Nvidia V100 and A100 GPUs
--------------------------------------------

We recommend the following modules and environement variables on Jean Zay:

.. code-block:: bash

    module load cuda/12.1.0
    module load gcc/12.2.0
    module load openmpi/4.1.1-cuda
    module load cmake/3.18.0

*Idefix* can then be configured to run on Nvidia V100 with the following options to ccmake:

``-DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_VOLTA70=ON``

While Ampere A100 GPUs are enabled with
``-DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_AMPERE80=ON``

MPI (multi-GPU) can be enabled by adding ``-DIdefix_MPI=ON`` as usual.

.. _setupSpecificOptions:

Setup-specific options
++++++++++++++++++++++

Some physical setup might require some ``cmake`` options to be set to specific value (e.g. an MHD setup will surely require MHD to be enabled).
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
