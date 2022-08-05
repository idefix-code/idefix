
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

For instance, the ``Hydro`` class provide the following list of enrollment functions (declared in hydro.hpp):

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
function should have the following signatures (declared in hydro_defs.hpp):

.. code-block:: c++

  using UserDefBoundaryFunc = void (*) (DataBlock &, int dir, BoundarySide side,
                                      const real t);
  using GravPotentialFunc = void (*) (DataBlock &, const real t, IdefixArray1D<real>&,
                                    IdefixArray1D<real>&, IdefixArray1D<real>&,
                                    IdefixArray3D<real> &);

  using SrcTermFunc = void (*) (DataBlock &, const real t, const real dt);
  using InternalBoundaryFunc = void (*) (DataBlock &, const real t);
  using EmfBoundaryFunc = void (*) (DataBlock &, const real t);
  using DiffusivityFunc = void (*) (DataBlock &, const real t, IdefixArray3D<real> &);
  using IsoSoundSpeedFunc = void (*) (DataBlock &, const real t, IdefixArray3D<real> &);

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
    data.hydro.EnrollGravPotential(&Potential);
  }


.. _userdefBoundaries:

User-defined boundaries
-----------------------

If one (or several) boundaries are set to ``userdef`` in the input file, the user needs to
enroll a user-defined boundary function in the ``Setup`` constructor as for the other user-def functions  (see :ref:`functionEnrollment`).
Note that even if several boundaries are ``userdef`` in the input file, only one user-defined function
is required. When *Idefix* calls the user defined boundary function, it sets the direction of the boundary (``dir=IDIR``, ``JDIR``,
or ``KDIR``) and the side of the bondary (``side=left`` or ``side=right``). If conveninent, one can use
the ``BoundaryFor`` wrapper functions to automatically loop on the boundary specified by ``dir`` and ``side``.
A typical user-defined boundary condition function looks like this:

.. code-block:: c++

  void UserdefBoundary(DataBlock& data, int dir, BoundarySide side, real t) {
    IdefixArray4D<real> Vc = data.hydro.Vc;
    IdefixArray4D<real> Vs = data.hydro.Vs;
    if(dir==IDIR) {
      data.hydro.boundary.BoundaryFor("UserDefBoundary", dir, side,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vc(RHO,k,j,i) = 1.0;
          Vc(VX1,k,j,i) = 0.0;
          Vc(VX2,k,j,i) = 0.0;
          Vc(VX3,k,j,i) = 0.0;
        });
      // For magnetic field (defined on cell sides), we need specific wrapper functions
      // Note that we don't need to initialise the field component parallel to dir, as it is
      // automatically reconstructed from the solenoidal condition and the tangential components
      data.hydro.boundary.BoundaryForX2s("UserDefBoundaryBX2s", dir, side,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vs(BX2s,k,j,i) = 0.0;
        });
      data.hydro.boundary.BoundaryForX3s("UserDefBoundaryBX3s", dir, side,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vs(BX3s,k,j,i) = 0.0;
        });
    }
  }



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

.. _setupInitDump:

Initialising from a restart dump
********************************

In some cases, it can be useful to initialise the flow from a dump taken from a previous
simulation. While one can simply use the ``-restart`` option on the commandline to resume
a simulation (see :ref:`commandLine`), there are some situation when one needs to create
a new initial condition by extrapolating or extanding a restart dump (such as in a resolution
test or a dimension change). In this case, one should use the ``DumpImage`` class which provides
all the tools needed to read a restart dump (see also :ref:`dumpImageClass`).

One typically first construct an instance of ``DumpImage`` in the ``Setup`` constructor, and then
use this instance to initialise the flow in ``Setup::InitFlow``. The procedure is examplified below,
assuming we want to create a dump from ``mydump.dmp``:

.. code-block:: c++

  DumpImage *image;       // Global pointer to our DumpImage

  // Setup constructor
  Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
    image = new DumpImage("mydump.dmp",output);   // load the dump file and store it in a DumpImage
  }

  // Flow initialisation, read directly from the DumpImage
  void Setup::InitFlow(DataBlock &data) {

    // Create a host copy
    DataBlockHost d(data);

    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {

          // Note that the restart dump array only contains the full (global) active domain
          // (i.e. it excludes the boundaries, but it is not decomposed accross MPI procs)
          int iglob=i-2*d.beg[IDIR]+d.gbeg[IDIR];
          int jglob=j-2*d.beg[JDIR]+d.gbeg[JDIR];
          int kglob=k-2*d.beg[KDIR]+d.gbeg[KDIR];

          d.Vc(RHO,k,j,i) = image->arrays["Vc-RHO"](kglob,jglob,iglob);
          d.Vc(PRS,k,j,i) = image->arrays["Vc-PRS"](kglob,jglob,iglob);
          d.Vc(VX1,k,j,i) = image->arrays["Vc-VX1"](kglob,jglob,iglob);
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
