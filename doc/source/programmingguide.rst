
.. _programmingGuide:

======================
Programming guide
======================

Because *Idefix* is designed to run on hybrid architecture, there are some subtilities which
are not found in other classical CPU-only codes. These subtilities are described in this section, along
with some basic coding rules which can be handful.

Data types
===========

Because *Idefix* can run on GPUs, and since GPUs experience a significant speedup when working
with single precision arithmetic, a specific ``real`` datatype is used for all floating point
operations in *Idefix*. This is by default aliased to ``double`` in `idefix.hpp`, but it can easily be modified
to ``float`` for specific problems. Note however that this loss of precision might have a strong
impact on the quality of the solution and is therefore not recommended.

Host and device
===============

*Idefix* relies on the Kokkos framework, and therefore assume that system it's running is made
of two sub-systems: a host and a device. The host is traditionnaly the CPU, and is taking care
of inputs and outputs, initialisation and allocation, MPI data exchanges. The device is usually an
accelerator (e.g. a GPU) and is actually performing the computation (or most of it).

Note that while *Idefix* assumes there is a host and a device, the device can be the host (think
of the code running only on your laptop CPU). In this case, Kokkos performs several optimisations,
so that everything effectively runs on the host smoothly.

By construction, the host doesn't have direct access to the device memory and vice-versa. This means
that accessing an array on the device from the host will inevitably lead to a segmentation fault.
This is a very common mistake, so keep this in mind when coding in *Idefix*.

Arrays
======
The fact that most of the computations are performed on the device implies that specific
allocations on the device need to be performed. To simplify the programmer's life, array allocations
should always use one of the ``IdefixArraynD`` where *n* is between 1 and 4. These arrays are
an alias for ``Kokkos::View`` which are fine-tuned for idefix, and are templated by the datatype,
which for most applications should be ``real``.

By definition, an IdefixArray is always allocated on the device, and is therefore not accessible
from the host. To define an array on the host, one use instead ``IdefixHostArraynD``.
Note that ``IdefixHostArraynD`` can be defined as mirrors of an ``IdefixArraynD``, to simplify data
transfer.

It should be noted that these
arrays are not just a contiguous memory zone as one would expect from a C array. Instead, these
arrays are C++ objects, with several properties. In particular, pointer arithmetic, which is a common
(bad) practice on C arrays, is not allowed on ``IdefixArray``. Allocation is performed as an instantiation
of the ``IdefixArraynD`` class. For instance, the following code block

.. code-block:: c++

  int arraySizeX1 = 10;   // 1st dimension of the arrays
  int arraySizeX2 = 10;   // 2nd dimension of the arrays.
  IdefixArray2D<real> myOldrArray = IdefixArray2D<real>("ArrayName", arraySizeX1, arraySizeX2);  // allocation
  IdefixArray2D<real> myNewArray = myOldArray;     // Shallow copy of myOldArray to myNewArray

will allocate myOldArray and perform a shallow copy of ``myOldArray``. No data is copied,
``myNewArray`` is merely a new reference to the same memory block on the device. Similarly,
accessing an element of an ``IdefixArraynD`` should always be with the accessor ``(...)``. In
the example above, one should use for instance ``myNewArray(1,2)`` (note the round brackets).

It is possible to copy data to/from the host/device manually using ``Kokkos::deep_copy`` (see examples
in ``DataBlockHost::syncToDevice()``). However, *Idefix* provides higher level functions
which should be sufficient for most uses through the classes ``DataBlockHost`` and ``GridHost``.

.. tip::
  ``IdefixArray`` internally contains a reference counter, so that when the array is not referenced
  anymore, the memory is automatically freed. Hence there is no equivalent of C ``free`` for
  ``IdefixArray``.

Execution space and loops
=========================
Just like arrays, code can be executed on the host or on the device. Unless otherwise mentionned, code
is by default executed on the host. Since the device is supposed to be performing the computation itself
and since this computation is usually performed using loops, *Idefix* provides z special way to handle
loops which are to be executed on the device, with the function ``idefix_for``.

``idefix_for`` is just a way to write a for loop, with some caveats. Depending on the kind of device
*Idefix* is running, the loop can be unrolled in an arbitrary order, and some iterations might
be executed simultaneously. So care should be taken when writting loops.

A typical loop on three indices looks like

.. code-block:: c++

  // Allocate an Idefix Array
  IdefixArray3D myArray<real> = IdefixArray3D<real>("MyArray", nx1, nx2, nx3);

  idefix_for("LoopName",
             kbeg,kend,
             jbeg,jend,
             ibeg,ieng,
             KOKKOS_LAMBDA (int k, int j, int i) {
                myArray(k,j,i) = 0.0;
              });

This loop will be executed on the device, and will perform a loop on indices i,j,k ranging from
ibeg to iend for i, etc. Note that as expected, we access the data stored in an ``IdefixArray``
inside an ``idefix_for``, i.e. code which is executed on the device.

The string "LoopName" should be descriptive of the loop (i.e. avoid "loop1", "loop2"...).
It is used when profiling or debugging the code and it names the execution kernels.

Note finally that the last argument of ``idefix_for`` relies on the ``KOKKOS_LAMBDA`` construct,
which implies that *Idefix* is actually making a C++ lambda when a loop is called.
While this should be transparent to most users, It should be kept in mind that these lambdas
capture their variables by value [=]. To avoid too much overhead, one should therefore avoid capturing
complex structures. Moreover, a bug in the Nvidia Cuda compiler ``nvcc`` prevents Cuda lambdas
from capturing class members (`see  this post <https://github.com/kokkos/kokkos/issues/695>`_). While
this bug is tightly linked to the C++11 norm and will be addressed in C++17, one should always
make local copies of the class members before using them in loops, to keep compatibility with Cuda
in C++11.

.. warning::
  As stated above, to avoid compatibility issues with nvcc, *always* make local copies (references)
  of the arrays and variables you intend to use before calling ``idefix_loop``. This ensures that
  these variables will be properly captured by device lambdas. It is the most common reason for
  GPU specific segmentation faults.

.. _classes:

Useful classes
==============

.. _inputClass:

The ``Input`` class
-------------------

``Input`` is a class which holds all of the information regarding command line and input file data. It provides accessors such as

.. code-block:: c++

  // Accessor to input parameters
  // the parameters are always: BlockName, EntryName, ParameterNumber (starting from 0)
  std::string GetString(std::string, std::string, int); // Read a string from the input file
  real GetReal(std::string, std::string, int);          // Read a real number from the input file
  int GetInt(std::string, std::string, int);            // Read an integer from the input file
  int CheckEntry(std::string, std::string);             // Check that a block/entry is present in the
                                                        // input file

Note that ``Input`` doesn't really read the input file each time an accessor is called. Internally,
``Input`` reads everything when constructed in a C++ container with all the data coming from the command line and the input file.
Hence there is no read overhead when one calls one of these accessor.

For instance, considering a .ini file::

  [MyBlock]
  myentry   1.0    0.0

It is possible to fetch the entry ``myentry`` using the ``Input`` accessors. Assuming an
instance of ``Input`` is allocated in ``myInput``:

.. code-block:: c++

  real firstParameter = myInput.GetReal("MyBlock","myentry",0)  // firstParameter=1.0
  real secondParameter = myInput.GetReal("MyBlock","myentry",1)  // secondParameter=0.0

If a parameter is not found, *Idefix* will print an error and exit. One can use the ``CheckEntry``
method to check if a parameter is set in the ini file before trying to access it.

.. tip::
  Command line options are also parsed by the ``Input`` class. These options are stored in a
  specific block named ``CommandLine``.

.. _gridClass:

``Grid`` class
------------------

``Grid`` is essentially a datastructure which represents the full computational domain (i.e. without domain decomposition,
if MPI has been enabled). It is useful when one needs to have access to the full grid coordinates for instance. Some of the useful arrays stored
by the grid are:

.. code-block:: c++

  IdefixArray1D<real> x[3];    // geometrical central points
  IdefixArray1D<real> xr[3];   // cell right interface
  IdefixArray1D<real> xl[3];   // cell left interface
  IdefixArray1D<real> dx[3];   // cell width

  real xbeg[3];           // Beginning of grid
  real xend[3];           // End of grid

  int np_tot[3];          // total number of grid points (including ghosts)
  int np_int[3];          // internal number of grid points (excluding ghosts)

.. _datablockClass:

``DataBlock`` class
-----------------------

``DataBlock`` contains all of the data structures that belongs to that particular process (i.e. if MPI is enabled, it contains data
specific to this subprocess, in contrast to ``Grid``). In particular, the DataBlocks have the local grid coordinates, stored
in arrays having the same name as ``Grid``. ``DataBlock`` also contains instances of the physical modules. Currently,
it only contains an instance of the ``Hydro`` class, but future physical modules will follow the same path.

.. _hydroClass:

``Hydro`` class
---------------------
The ``Hydro`` class (and its sub-classes) contains all of the fields and methods specific to (magneto) hydrodynamics. While
interested users may want to read in details the implementation of this class, we provide below a list of the most important
members

.. code-block:: c++

  IdefixArray4D<real> Vc;      // Main cell-centered primitive variables index
  IdefixArray4D<real> Vs;      // Main face-centered varariables
  IdefixArray4D<real> Uc;      // Main cell-centered conservative variables

  // Enroll user-defined gravitational potential
  void EnrollGravPotential(GravPotentialFunc);

  // Enroll user-defined boundary conditions
  void EnrollUserDefBoundary(UserDefBoundaryFunc);

  // Enroll user-defined ohmic, ambipolar and Hall diffusivities
  void EnrollOhmicDiffusivity(DiffusivityFunc);
  void EnrollAmbipolarDiffusivity(DiffusivityFunc);
  void EnrollHallDiffusivity(DiffusivityFunc);

  // Enroll user-defined isothermal sound speed
  void EnrollIsoSoundSpeed(IsoSoundSpeedFunc);


The first two IdefixArrays are the ones storing the primitive variable fields. These arrays
are 4D, the first dimension being the field number. *Idefix* defines aliases for these numbers,
so that one can call ``Vc(VX1,k,j,i)`` in place of ``Vc(1,k,j,i)`` to get the first velocity component.
These aliases are defined in ``idefix.hpp``

Because the code uses contrained transport, the field defined on cell faces is stored in the ``Vs``
array. Just like for ``Vc``, there are aliases, with "s" suffixes defined to simplify the adressing
of the magnetic field components, as ``Vs(BX2s,k,j,i)``.

.. _datablockhostClass:

``DataBlockHost`` class
---------------------------
This class is a *mirror* class, which is designed to be a (partial) copy of the ``DataBlock`` class,
but in which all of the arrays are stored on the *host*. Obviously, ``DataBlockHost`` comes handy
when one has to deal with input/outputs, debugging and initialisation.

The ``DataBlockHost`` should always be constructed with a ``DataBlock`` in argument. This ensures
that the ``DataBlockHost`` knows where its parent ``DataBlock`` is located. When created, a ``DataBlockHost``
fills the following arrays (essentially grid information) with data from its parent ``DataBlock``:

.. code-block:: c++

  IdefixArray1D<real>::HostMirror x[3];   // geometrical central points
  IdefixArray1D<real>::HostMirror xr[3];  // cell right interface
  IdefixArray1D<real>::HostMirror xl[3];  // cell left interface
  IdefixArray1D<real>::HostMirror dx[3];  // cell width

  IdefixArray3D<real>::HostMirror dV;     // cell volume
  IdefixArray3D<real>::HostMirror A[3];   // cell right interface area


Note however that the physics arrays are not automatically synchronized when ``DataBlockHost`` is
created, that is:

.. code-block:: c++

  IdefixArray4D<real>::HostMirror Vc;     // Main cell-centered primitive variables index
  IdefixArray4D<real>::HostMirror Vs;     // Main face-centered primitive variables index
  IdefixArray4D<real>::HostMirror J;      // Current (only when haveCurrent is enabled)
  IdefixArray4D<real>::HostMirror Uc;     // Main cell-centered conservative variables
  IdefixArray3D<real>::HostMirror InvDt;

  IdefixArray3D<real>::HostMirror Ex1;    // x1 electric field
  IdefixArray3D<real>::HostMirror Ex2;    // x2 electric field
  IdefixArray3D<real>::HostMirror Ex3;    // x3 electric field

need to be synchronized *manually*. These IdefixArrays are all defined as ``HostMirror``, implying that they are accessible
from the host only. If modifications are performed on the arrays of the
parent ``DataBlock``, one can call ``DataBlockHost::SyncFromDevice()`` to refresh the host arrays,
and inversely one can call ``DataBlockHost::SyncToDevice()`` to send data from ``DataBlockHost``
to the parent ``DataBlock``.

Finally, ``DataBlockHost`` provides a useful method ``DataBlockHost::MakeVsFromAmag(IdefixHostArray4D<real> &)``
which can be used to initialise the face-centered magnetic field stored in ``DataBlockHost::Vs`` from a user-defined
magnetic potential. See :ref:`setupInitflow`.



Debugging and profiling
=======================