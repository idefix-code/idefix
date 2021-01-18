
Problem Setup ``setup.cpp``
===========================
The source code ``setup.cpp`` contains the code specific to physical setup at hand. It should serve several purposes:
  - Get the setup parameters from the problem input file.
  - Init arrays and variables you would need
  - Set the initial conditions
  - Define setup-specific functions which will be called from the inegration loop
  - Define setup-specific outputs

Most of these initialisations relies on the class ``Setup`` which has to be implemented in your ``setup.cpp``. At this stage it is useful
to have in mind some of the classes defined in *Idefix* to which you can have access to, which are documented in :ref:`classes`.

The ``Setup`` class
--------------------
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
which have already been initialised when ``Setup`` is called: ``Input``, ``Grid`` and ``DataBlock`` (see :ref:`classes`).








Function enrollment
-------------------




Writing tips
------------
