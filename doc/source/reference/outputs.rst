.. _output:

Outputs
=======

Output formats
--------------

*Idefix* uses several types of outputs you may want for your setup. By default, *Idefix* allows
for 3 kinds of outputs:

* logs which essentially tells the user what *Idefix* is currently doing. When running in serial, logs are sent to stdout, but when
  MPI is enabled, only the logs of the rank 0 process is sent to stdout, and each process (including rank 0) simultaneously writes a
  log file `idefix.n.log` where *n* is the process MPI rank.
* dump files (.dmp) which are *Idefix* specific binary files containing all of the data at machine precision to restart your run.
  These files are therefore the ones which are read when *Idefix* is restarted.
* VTK files (.vtk) are Visualation Toolkit files, which are easily readable by visualisation softwares such as `Paraview <https://www.paraview.org/>`_
  or `Visit <https://wci.llnl.gov/simulation/computer-codes/visit>`_. A set of python methods is also provided to read vtk file from your
  python scripts in the `pytools` directory.
* user-defined analysis files. These are totally left to the user. They usually consist of ascii tables defined by the user, but they can
  be anything.

The output periodicity and the userdef variables should all be declared in the input file, as described in :ref:`outputSection`.

Defining your own outputs
-------------------------

*Idefix* provides two ways to define your own outputs: analysis, which are used to make your
own output file (e.g. an ascii-tabulated file); and user variables, which are written by *Idefix* output routines.

Both analysis and uservar requires the definition of a user function which needs to be enrolled following the procedure described
in :ref:`functionEnrollment` and using the function signatures declared in `output.hpp`.

We provide below an example of a setup using both analysis outputs and uservar outputs


.. code-block::
  :caption: Input file `idefix.ini`

  [Output]
    analysis  0.01                # A user-defined analysis will be performed every 0.01
    vtk       1.0                 # VTK files will be written every 1.0
    uservar   rhovx rhovy         # Two user variables will be defined and written in vtk files

.. code-block:: c
  :caption: Setup file `setup.cpp`

  // Analyse data to produce an ascii output
  void Analysis(DataBlock & data) {
    // Mirror data on Host
    DataBlockHost d(data);

    // Sync it
    d.SyncFromDevice();

    // Get the field at some specific location
    real by = d.Vc(BX2,0,0,0);

    // Write the data in ascii to our file
    std::ofstream f;
    f.open("timevol.dat",std::ios::app);
    f.precision(10);
    f << std::scientific << data.t << "\t" << by << std::endl;
    f.close();
  }

  // Compute user variables which will be written in vtk files
  void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {
    // Mirror data on Host
    DataBlockHost d(data);

    // Sync it
    d.SyncFromDevice();

    // Make references to the user-defined arrays (variables is a container of IdefixHostArray3D)
    // Note that the labels should match the variable names in the input file
    IdefixHostArray3D<real> rhovx = variables["rhovx"];
    IdefixHostArray3D<real> rhovy = variables["rhovy"];

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
      for(int j = 0; j < d.np_tot[JDIR] ; j++) {
        for(int i = 0; i < d.np_tot[IDIR] ; i++) {
          rhovx(k,j,i) = d.Vc(VX1,k,j,i)*d.Vc(RHO,k,j,i);
          rhovy(k,j,i) = d.Vc(VX2,k,j,i)*d.Vc(RHO,k,j,i);
        }
      }
    }
  }

  // Setup constructor, initialises our setup
  Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
    // Enroll our analysis function
    output.EnrollAnalysis(&Analysis);

    // Enroll our user-defined variables
    output.EnrollUserDefVariables(&ComputeUserVars);
  }

  void Setup::InitFlow(DataBlock &data) {
  // Not shown here
  }
