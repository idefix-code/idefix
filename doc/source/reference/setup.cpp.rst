
Problem Setup ``setup.cpp``
===========================
The source code ``setup.cpp`` contains the code specific to physical setup at hand. It should serve several purposes:
  - Get the setup parameters from the problem input file.
  - Init arrays and variables you would need
  - Set the initial conditions
  - Define setup-specific functions which will be called from the inegration loop
  - Define setup-specific outputs

Most of these initialisations relies on the class ``Setup`` which has to be implemented in your ``setup.cpp``. At this stage it is useful
to have in mind some of the classes defined in *Idefix* to which you can have access to.

Data structures and execution space in *Idefix*
===============================================



Some useful classes of *Idefix*
===================================

*Idefix* being written in c++, it defines a full set of classes can be called from the setup


The ``Setup`` class
-------------------
The ``Setup`` class is declared as follows:

.. code-block:: c++

  class Setup {
    public:
      Setup();
      Setup(Input &, Grid &, DataBlock &);
      void InitFlow(DataBlock &);
      void MakeAnalysis(DataBlock&);
    };

As it can be seen this class consist of constructor and two methods: ``InitFlow`` and ``MakeAnalysis`` which will handle
the initial condition and the setup-specific outputs mentionned above.

Let us start with the constructors. The default constructor ``Setup()`` is not used by *Idefix*. The code only use the constructor `Setup(Input &, Grid &, DataBlock &)`.
This constructor is called on the code startup and allows the user to load and set setup-specific parameters. The parameters are three objects
which have already been initialised when ``Setup`` is called: ``Input``, ``Grid`` and ``DataBlock``.

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

Note that ``Input`` doesn't really read the input file each time the accessor is called. Internally,
Input stores everything when constructed in a c++ container with all the data coming from the command line and the input file.
Hence there is no read overhead when one calls one of these accessor.

.. tip::
  Command line options are also parsed by the ``Input`` class they are stored in the specific block named ``CommandLine``.

The ``Grid`` class
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


The ``DataBlock`` class
-----------------------

``DataBlock`` contains all of the data structures that belongs to that particular process (i.e. if MPI is enabled, it contains data
specific to this subprocess, in contrast to ``Grid``). DataBlock have the same arrays as the ones defined









Function enrollment
-------------------




Writing tips
------------
