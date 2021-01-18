
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









Function enrollment
-------------------




Writing tips
------------
