
Problem Setup ``setup.cpp``
===========================
The source code ``setup.cpp`` contains the code specific to the physical setup at hand. It serves several purposes:
  - Get the setup parameters from the problem input file.
  - Init arrays and variables you will need
  - Set the initial conditions
  - Define and enroll setup-specific functions which will be called from the integration loop
  - Define and enroll setup-specific outputs

Most of these initialisations relies on the class ``Setup`` which has to be implemented in your
``setup.cpp``. At this stage, if you are not familliar with the ``idefix_loop`` structures,
``IdefixArray`` data arrays and *Idefix* class structure, it is strongly recommended you have a
look at the :ref:`programmingGuide`.

The ``Setup`` class
--------------------
The ``Setup`` class is declared as follows:

.. code-block:: c++

  class Setup {
    public:
      Setup(Input &, Grid &, DataBlock &, Output &);
      void InitFlow(DataBlock &);
    };

As it can be seen, this class consist of a constructor and one mandatory method: ``InitFlow``, which will handle
the initial condition.

The ``Setup`` constructor
-------------------------
Let us start with the constructor `Setup(Input &, Grid &, DataBlock &, Output &)`.
This constructor is called on the code startup and allows the user to load and set setup-specific parameters. The parameters are four objects
which have already been initialised when ``Setup`` is called: ``Input``, ``Grid``, ``DataBlock`` and ``Output`` (see :ref:`classes`).

A typical constructor first loads the setup parameters calling accessors from the ``Input`` object (see :ref:`inputClass`). Then,
if there are some user-defined functions (for instance a user-defined potential, boundary condition or output),
the constructor also *enrolls* these functions before returning.

.. _functionEnrollment:

Function enrollment
*******************

The enrollment of user functions is required in *Idefix* whenever a parameter is set to "userdef" in
the input file, and for user-defined outputs. This can be seen as a way to link the user
setup to the main code at runtime, and avoid the need to pre-define tens of empty functions. Function enrollment
is achieved by calling one of the ``EnrollXXX`` function of the class associated to it.

For instance, the ``Fluid`` class provide the following list of enrollment functions (declared in hydro.hpp):

.. code-block:: c++

  // Enroll user-defined boundary conditions
  void EnrollUserDefBoundary(UserDefBoundaryFunc);
  void EnrollInternalBoundary(InternalBoundaryFunc);
  void EnrollEmfBoundary(EmfBoundaryFunc);

  // Enroll user-defined gravitational potential
  void EnrollGravPotential(GravPotentialFunc);

  // Enroll user source terms
  void EnrollUserSourceTerm(SrcTermFunc);

  // Enroll user-defined ohmic, ambipolar and Hall diffusivities
  void EnrollOhmicDiffusivity(DiffusivityFunc);
  void EnrollAmbipolarDiffusivity(DiffusivityFunc);
  void EnrollHallDiffusivity(DiffusivityFunc);

  // Enroll user-defined isothermal sound speed
  void EnrollIsoSoundSpeed(IsoSoundSpeedFunc);

When called, these function expects the address of the user-defined function. These user-defined
function should have the following signatures (declared in fluid_defs.hpp and boundary.hpp):

.. code-block:: c++

  using UserDefBoundaryFunc = void (*) (Fluid<Phys> *, int dir, BoundarySide side,
                                      const real t);
  using GravPotentialFunc = void (*) (DataBlock &, const real t, IdefixArray1D<real>&,
                                    IdefixArray1D<real>&, IdefixArray1D<real>&,
                                    IdefixArray3D<real> &);

  using SrcTermFunc = void (*) (Fluid<Phys> *, const real t, const real dt);
  using InternalBoundaryFunc = void (*) (Fluid<Phys>*, const real t);
  using EmfBoundaryFunc = void (*) (DataBlock &, const real t);
  using DiffusivityFunc = void (*) (DataBlock &, const real t, IdefixArray3D<real> &);
  using IsoSoundSpeedFunc = void (*) (DataBlock &, const real t, IdefixArray3D<real> &);


Note that some of these functions involve the template class ``Fluid<Phys>``. The ``Fluid`` class
is indeed capable of handling several types of fluids (described by the template parameter ``Phys``):
MHD, HD, pressureless, etc... Hence, depending on the type of fluid to which the user-defined
function applies, the signature of the function would be different. For instance, a User-defined
boundary condition for a system that would solve for a gas+dust mixture would read

.. code-block:: c++

  void MyBoundary(Fluid<DefaultPhysics> * fluid, int dir, BoundarySide side, const real t) {
  // Here comes the code for the Gas boundary condition
  }

  void MyBoundaryDust(Fluid<DustPhysics> * fluid, int dir, BoundarySide side, const real t) {
  // Here comes the code for the dust boundary condition
  }

Note that *Idefix* defines an alias for the default fluid which is often found in the example provided:

.. code-block:: c++

  using Hydro = Fluid<DefaultPhysics>;


Example
*******

The following example have a user-defined gravitational potential, and defines a ``Setup``
constructor which reads a parameter from the .ini file and enroll the user-defined potential.

.. code-block:: c++

  // a global variable which stores the mass of some object
  real Mass;

  // user-defined potential
  void Potential(DataBlock& data, const real t, IdefixArray1D<real>& x1, IdefixArray1D<real>& x2, IdefixArray1D<real>& x3, IdefixArray3D<real>& phi) {
    idefix_for("Potential",0,data.np_tot[KDIR], 0, data.np_tot[JDIR], 0, data.np_tot[IDIR],
               KOKKOS_LAMBDA (int k, int j, int i) {
                  phi(k,j,i) = -Mass/x1(i);
              });

  }

  // Setup constructor
  Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
    // Read some parameter from the ini file
    Mass = input.Get<real>("Setup","mass",0);

    // Enroll the user-defined potential
    data.gravity->EnrollGravPotential(&Potential);
  }


.. _userdefBoundaries:

User-defined boundaries
-----------------------

If one (or several) boundaries are set to ``userdef`` in the input file, the user needs to
enroll a user-defined boundary function in the ``Setup`` constructor as for the other user-def functions  (see :ref:`functionEnrollment`).
Note that even if several boundaries are ``userdef`` in the input file, only one user-defined function
is required per fluid type. When *Idefix* calls the user defined boundary function, it sets the direction of the boundary (``dir=IDIR``, ``JDIR``,
or ``KDIR``) and the side of the bondary (``side=left`` or ``side=right``). If conveninent, one can use
the ``BoundaryFor`` wrapper functions to automatically loop on the boundary specified by ``dir`` and ``side``.
A typical user-defined boundary condition function looks like this:

.. code-block:: c++

  void UserdefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {
    IdefixArray4D<real> Vc = hydro->Vc;
    IdefixArray4D<real> Vs = hydro->Vs;
    if(dir==IDIR) {
      hydro->boundary->BoundaryFor("UserDefBoundary", dir, side,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vc(RHO,k,j,i) = 1.0;
          Vc(VX1,k,j,i) = 0.0;
          Vc(VX2,k,j,i) = 0.0;
          Vc(VX3,k,j,i) = 0.0;
        });
      // For magnetic field (defined on cell sides), we need specific wrapper functions
      // Note that we don't need to initialise the field component parallel to dir, as it is
      // automatically reconstructed from the solenoidal condition and the tangential components
      hydro->boundary->BoundaryForX2s("UserDefBoundaryBX2s", dir, side,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vs(BX2s,k,j,i) = 0.0;
        });
      hydro->boundary->BoundaryForX3s("UserDefBoundaryBX3s", dir, side,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vs(BX3s,k,j,i) = 0.0;
        });
    }
  }

.. warning::

  Only the tangential field components should be initialised by user-defined boundary conditions.
  *Idefix* automatically reconstruct (and overwrite!) the normal field component from the
  divergence-free condition on B and the user-defined tangential magnetic field components.



.. _setupInitflow:

``Setup::InitFlow`` method
--------------------------

Basics of the Initflow method
*****************************

``InitFlow`` is a method of the ``Setup`` class and is called by *Idefix* after the ``Setup`` constructor.
Its role is to define the initial conditions for the flow, initializing the ``Vc`` (and ``Vs`` in MHD)
arrays of the ``Hydro`` class, for instance. Because this method does not have to be fast, since
it is called only once, it is customary to initialise the flow on the host, and then send it to the
device.

For this, it is useful to first define a mirror ``DataBlockHost`` (see :ref:`datablockhostClass`)
of the ``DataBlock`` given in argument and initialse the flow in ``DataBlockHost`` using a standard
C loop on the host, as in the example below.

.. code-block:: c++

  void Setup::InitFlow(DataBlock &data) {
    // Create a host copy of the DataBlock given in argument
    DataBlockHost dataHost(data);

    // Because we initialise the arrays in DataBlockHost,
    // we can execute the loop on the host
    for(int k = 0; k < dataHost.np_tot[KDIR] ; k++) {
        for(int j = 0; j < dataHost.np_tot[JDIR] ; j++) {
            for(int i = 0; i < dataHost.np_tot[IDIR] ; i++) {
                real x = dataHost.x[IDIR](i);
                real y = dataHost.x[JDIR](j);
                real z = dataHost.x[KDIR](k);

                dataHost.Vc(RHO,k,j,i) = 1.0;
                dataHost.Vc(PRS,k,j,i) = 1.0;
                dataHost.Vc(VX1,k,j,i) = -sin(y);
                dataHost.Vc(VX2,k,j,i) = sin(x)+cos(z);
                dataHost.Vc(VX3,k,j,i) = cos(x);

                dataHost.Vs(BX1s,k,j,i) = -sin(y);
                dataHost.Vs(BX2s,k,j,i) = sin(x);
                dataHost.Vs(BX3s,k,j,i) = 0.0;
            }
        }
    }
    // Do not forget to send our initialisation to the parent dataBlock!
    dataHost.SyncToDevice();
  }

.. warning::
  Do not forget to sync your DataBlockHost to its parent DataBlock using the
  ``DataBlockHost::SyncToDevice()`` method!

Initialising the magnetic field
*******************************

When MHD is used, the face-centered magnetic field stored in ``Vs`` should be initialised with a divergence-free
field *at machine precision*. This might not always be straightforward for some complex field geometry,
so *dataBlockHost* can also be initialised with a vector potential, from which the face-centered field
can be automatically derived using ``DataBlockHost::MakeVsFromAmag`` as in the example below:

.. code-block:: c++

  void Setup::InitFlow(DataBlock &data) {
    // Create a host copy of the DataBlock given in argument
    DataBlockHost dataHost(data);

    // Allocate an array on host to store the vector potential (3 components are expected)
    IdefixHostArray4D<real> A = IdefixHostArray4D<real>("Setup_VectorPotential", 3,
                                                        data.np_tot[KDIR],
                                                        data.np_tot[JDIR],
                                                        data.np_tot[IDIR]);

    for(int k = 0; k < dataHost.np_tot[KDIR] ; k++) {
      for(int j = 0; j < dataHost.np_tot[JDIR] ; j++) {
        for(int i = 0; i < dataHost.np_tot[IDIR] ; i++) {
          real x = dataHost.x[IDIR](i);
          real y = dataHost.x[JDIR](j);
          real z = dataHost.x[KDIR](k);

          // Initialise Vc field (not shown)
          // ...

          // Initialise the 3 components of the vector potential
          A(IDIR,k,j,i) = 0.0;
          A(JDIR,k,j,i) = 0.0;
          A(KDIR,k,j,i) = -y*B0;
        }
      }
    }

    // Compute the face centered Vs from the vector potential
    dataHost.MakeVsFromAmag(A);

    // Do not forget to send our initialisation to the parent dataBlock!
    dataHost.SyncToDevice();
  }

Initialising passive tracers
********************************

Idefix provides the possibility to use passive tracers, that are scalars passively advected by
the associated fluid. Tracers are enabled by setting on non-zero integer to the parameter ``tracer`` in
a ``[Hydro]`` (for hydro tracers) or a ``[Dust]`` block in your input file, so that the tracers
will follow the main fluid or the dust fluids, with the given number of tracers.

Each tracer can be accessed in the dataBlockHost as a particular field in the array Vc as in the example below

.. code-block:: c++

  void Setup::InitFlow(DataBlock &data) {
      // Create a host copy of the DataBlock given in argument
      DataBlockHost dataHost(data);

      for(int k = 0; k < dataHost.np_tot[KDIR] ; k++) {
        for(int j = 0; j < dataHost.np_tot[JDIR] ; j++) {
          for(int i = 0; i < dataHost.np_tot[IDIR] ; i++) {
            real x = dataHost.x[IDIR](i);
            real y = dataHost.x[JDIR](j);
            // First tracer
            dataHost.Vc(TRG,k,j,i) = (x < 0 ? 0 : 1); // TRG for gas tracers
            // second tracer
            dataHost.Vc(TRG+1,k,j,i) = (y < 0 ? 0 : 1); // For gas tracers

.. note::

  Note that when using dust tracers, one should use the field ``TRD`` instead of ``TRG``.

.. _setupInitDump:

Initialising from a restart dump
********************************

In some cases, it can be useful to initialise the flow from a dump taken from a previous
simulation. While one can simply use the ``-restart`` option on the commandline to resume
a simulation (see :ref:`commandLine`), there are some situation when one needs to create
a new initial condition by extrapolating or extanding a restart dump (such as in a resolution
test or a dimension change). In this case, one should use the ``DumpImage`` class which provides
all the tools needed to read a restart dump (see also :ref:`dumpImageClass`).

One typically first construct an instance of ``DumpImage`` in ``Setup::InitFlow``, and then
use this instance to initialise the flow. The procedure is examplified below,
assuming we want to create a dump from ``mydump.dmp``:

.. code-block:: c++

  #include "dumpImage.hpp"

  // Flow initialisation, read directly from the DumpImage
  void Setup::InitFlow(DataBlock &data) {

    // Create a host copy
    DataBlockHost d(data);

    DumpImage image("mydump.dmp", &data);

    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {

          // Note that the restart dump array only contains the full (global) active domain
          // (i.e. it excludes the boundaries, but it is not decomposed accross MPI procs)
          int iglob=i-2*d.beg[IDIR]+d.gbeg[IDIR];
          int jglob=j-2*d.beg[JDIR]+d.gbeg[JDIR];
          int kglob=k-2*d.beg[KDIR]+d.gbeg[KDIR];

          d.Vc(RHO,k,j,i) = image.arrays["Vc-RHO"](kglob,jglob,iglob);
          d.Vc(PRS,k,j,i) = image.arrays["Vc-PRS"](kglob,jglob,iglob);
          d.Vc(VX1,k,j,i) = image.arrays["Vc-VX1"](kglob,jglob,iglob);
  }}}

    // For magnetic variable, we should fill the entire active domain, hence an additional
    // point in the field direction
    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
          for(int i = d.beg[IDIR]; i < d.end[IDIR]+IOFFSET ; i++) {
            int iglob=i-2*d.beg[IDIR]+d.gbeg[IDIR];
            int jglob=j-2*d.beg[JDIR]+d.gbeg[JDIR];
            int kglob=k-2*d.beg[KDIR]+d.gbeg[KDIR];
            d.Vs(BX1s,k,j,i) = image->arrays["Vs-BX1s"](kglob,jglob,iglob);
    }}}

    // And so on for the other components
    // ..


    delete image;   // don't forget to free the memory allocated for dumpImage!

    // Send our datablock to the device
    d.SyncToDevice();
  }


.. note::

  Note that the naming convention in ``DumpImage::arrays`` combines the original array and variable names.
  It is generically written ``XX-YYY`` where ``XX`` is the array name in the ``dataBlock`` (e.g.
  ``Vc`` or ``Vs``) and ``YYY`` is the variable name (e.g. ``VX2`` or ``BX3s``).


User-defined analysis
---------------------

User-defined analysis and outputs can be coded in the ``setup.cpp`` file. Follow the
guidelines in :ref:`output`.

I need a global IdefixArray for my Setup
-----------------------------------------

There are situation where you will need one or several global IdefixArrays that can be accessed from different
functions, e.g. the ``Initflow`` method and the user-defined boundary conditions.

It is important to understand that IdefixArrays (equivalent to ``Kokkos::view`` that are references to memory chunks)
are automatically dealocated when all of the IdefixArrays refereing to that memory chunk have been deleted. This deletion happens either
implicitly (by a closing scope) in which case the objects contained in the scope are all deleted automatically,
or explicitly (through a new/delete pair).

If you define an IdefixArray in the global scope, it is deleted when the program terminates. Hence deallocation should happen then.
Except that, according to `Kokkos documentation <https://kokkos.github.io/kokkos-core-wiki/ProgrammingGuide/Initialization.html#finalization>`_,  we
need to call ``Kokkos::finalize`` before the program terminates and this ``finalize`` should be done once all of
the Kokkos objects have been deleted (including IdefixArray). While *Idefix* makes sure that all of its objects (including the user's ``Setup``) are being deleted before calling ``finalize``,
a simple IdefixArray in the global scope will not be explicitely deleted, and will typically lead to the following error:

.. code-block:: bash

  terminate called after throwing an instance of 'std::runtime_error'
  what():  Kokkos allocation "MyAwesomeArray" is being deallocated after Kokkos::finalize was called

The way to avoid this is to explicitely delete the object when you don't need it anymore. The cleanest way to do this for a setup is to define a "container" class,
containing all of the arrays you will need in the global scope, and just have a global pointer to an instance of this class, that you eventually delete
(and which deletes all of the arrays it contains automatically). More explicitely:

#. start with a declaration of a class container (that we name MyGlobalClass in this example) and a global pointer to a class instance (note that you can put as many arrays as you want in the class)

    .. code-block:: c++

      // Class declaration
      class MyGlobalClass {
      public:
        // Class constructor
        MyGlobalClass(DataBlock &data) {
        //allocate some memory for the array the class contains
          this->array1 = IdefixArray3D<real>("MyAwesomeArray",data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
        }

        // array1, member of the class
        IdefixArray3D<real> array1;
      };

      // A global class instance named "myGlobals"
      MyGlobalClass *myGlobals;

#. initialise your global object in the Setup constructor (this will aumatically allocate the array it contains thanks to the class constructor we have defined):

    .. code-block:: c++

      Setup::Setup(....) {
        ...
        myGlobals = new MyGlobalClass(data);
        ...
      }

#. to avoid the error message above, don't forget to delete the object on exit in the Setup destructor

    .. code-block:: c++

      Setup::~Setup(....) {
        ...
        delete myGlobals;
        ...
      }

#. and finally, use your array when you need it:

    .. code-block:: c++

      MyXXXXFunction(....) {
        // Shallow copy the global array
        IdefixArray3D<real> array = myGlobals->array1;
        // Do stuff
        ....
      }
