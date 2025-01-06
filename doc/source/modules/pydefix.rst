.. _pydefixModule:

Pydefix module
==============

The Pydefix module allows Idefix to talk directly to a Python interpreter while running. It can be used to create your own initial condition
and/or for custom outputs produced from Python. Pydefix relies on the pybind11 python package

The essence of Pydefix is to allows the user to have a direct access to Idefix data structure from Python without writing/accessing any file. In particular, IdefixArrays are viewed as numpy arrays in Python.
Note however that to keep things simple, Pydefix works on the host memory space only, and hence sync data to/from the GPU (if used) before invoking Python functions. Hence, using Pydefix for outputs induces
a significant loss of performances.


Before you start
----------------
Pybind11 installation
+++++++++++++++++++++

In order to use pydefix, you need to be working in a python>=3.8 environement that includes `pybind11 <https://pybind11.readthedocs.io>`_. Follow the instruction of your package manager to install pybind11>=2.12.

Pydefix usage
-------------
Idefix Configuration
++++++++++++++++++++

In order to use Pydefix, you need to switch on ``Idefix_PYTHON`` in cmake. This will auto-detect Python and check that pybind11 can be used effectively.


Run Idefix with Pydefix
+++++++++++++++++++++++

Pydefix is typically enabled from your input file `idefix.ini` in the block ``[Python]``. The following parameters are available:

+------------------------+-----------------------+-----------------------------------------------------------------------------------------------------------+
|  Entry name            | Parameter type        | Comment                                                                                                   |
+========================+=======================+===========================================================================================================+
| script                 | string                | | (Mandatory) Filename (*without ".py"!*) of the python script that Idefix should use.                    |
|                        |                       | | The script should be in the same location as the Idefix executable file.                                |
+------------------------+-----------------------+-----------------------------------------------------------------------------------------------------------+
| output_function        | string                | | (Optional) Name of the function that will be called for each output event (the function should be       |
|                        |                       | | defined in the  python script above). When ommited, pydefix output functions are disabled.              |
|                        |                       | | The periodicity of the pydefix output routine is set in the block:entry `[Output]:python`               |
+------------------------+-----------------------+-----------------------------------------------------------------------------------------------------------+
| initflow_function      | string                | | (Optional) Name of the python function that will be called to initialize the flow in place of the C++   |
|                        |                       | | function ``Setup::InitFlow``. Revert to ``Setup::Initflow`` when ommited.                               |
+------------------------+-----------------------+-----------------------------------------------------------------------------------------------------------+

Python script
+++++++++++++

When using Pydefix, idefix expects a python script to be specified in the input file (see ``script`` parameter above). To be fully functionnal, you should import the ``pydefix`` module at the beginning
of your python script (you can also import other python module, as any python script).

Your python script should define functions that will be called while Idefix is running:
* The signature of the ``initflow`` function should be ``(data)`` where ``data`` is a python structure matching Idefix's ``DataBlockHost`` class.
* The signature of the ``output`` function should be ``(data,grid,n)`` where ``data`` is a python structure matching Idefix's ``DataBlockHost`` class, ``grid`` is Idefix's ``GridHost`` class, and ``n`` is an integer representing the current number of the output

MPI parallelism
+++++++++++++++
When Idefix runs with MPI parallelism enabled, a python interpreter and script is launched by each MPI process. Each of these script is independent
and have access to its local ``dataBlockHost``. The `pydefix` module however gives access to the local rank ``prank`` and total MPI size ``psize``. In addition,
pydefix provides the function ``GatherIdefixArray`` to gather the data distributed among each process without invoking MPI directly in python. This function
expects a 3D distributed IdefixArray in entry following the signature

.. code-block:: c++

  GatherIdefixArray(IdefixHostArray3D<real> in, // 3D distributed array
                    DataBlockHost dataHost,     // dataBlock structure
                    bool keepBoundaries = true, // Whether we keep the ghost zones in the returned array
                    bool broadcast = true)      // Whether the returned array is available only in proc #0 or in every proc (caution! possibly requires lots of memory)

This function is used as follows:

.. code-block:: python

  import pydefix as pdfx # mandatory
  import numpy as np
  import matplotlib.pyplot as plt

  def output(data,grid,n):
    # Gather pressure field from every process in process #0 (set broadcast to True to distribute the result to all processes)
    prs = pdfx.GatherIdefixArray(data.Vc[pdfx.PRS,:,:,:],data,broadcast=False)

    # Only root process performs this
    if pdfx.prank==0:
      x=grid.x[pdfx.IDIR] # The grid contains the full domain size, while the datablock contains only local information
      y=grid.x[pdfx.JDIR]
      plt.figure()
      plt.pcolormesh(x,y,prs[0,:,:],cmap='plasma')


.. note::
  For more advanced usage, it is also possibly to directly call MPI routines from python using the `Mpi4py <https://pypi.org/project/mpi4py/>`_ module.

Example
+++++++

An example is provided in the directory `test/IO/python`. This example shows how to use Idefix with pure python initial conditions and outputs.
It reproduces the 2D OrszagTang vortex available in MHD/OrszagTang without requiring any single line of C++ code from the user.

The python script `pydefix_example.py` initializes the initial condition of the OT test (``initflow``) and produces a series of PNG files through matplotlib (`output`).

Troubleshooting
---------------

It during configuration stage, you get::

  CMake Error at CMakeLists.txt:122 (find_package):
    By not providing "Findpybind11.cmake" in CMAKE_MODULE_PATH this project has
    asked CMake to find a package configuration file provided by "pybind11",
    but CMake did not find one.

It means that cmake cannot find the location of pybind11 (this typically happens on MacOs). In order to locate pybind11, open a python interpreter, and get pybind11 install dir through:

.. code-block:: python

  import pybind11
  print(pybind11.__file__)

You can then exit the interpreter and set the pybind11_DIR environement variable to the right path:

.. code-block:: bash

  export pybind11_DIR=env/lib/python3.10/site-packages/pybind11

you can then run cmake which should be able to find pybind11, and compile the code.
