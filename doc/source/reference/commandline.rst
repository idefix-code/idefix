Command line & Signal handling
==============================

.. _commandLine:

Command line options
--------------------

Several options can be provided at command line when running the code. These are listed below

+--------------------+-------------------------------------------------------------------------------------------------------------------------+
| Option name        | Comment                                                                                                                 |
+====================+=========================================================================================================================+
| -dec n1 n2 n3      | | Specify the MPI domain decomposition. Idefix will decompose the domain with n1 MPI processes in X1,                   |
|                    | | n2 MPI processes in X2 and n3 processes in X3. Note the number of arguments to -dec should be equal to ``DIMENSIONS``.|
+--------------------+-------------------------------------------------------------------------------------------------------------------------+
| -restart n         | | Restart from the ``n``^th dump file. By default, ``n`` matches the highest value from existing dump files.            |
|                    | | When used, the initial conditions from ``Setup::InitFlow()`` are ignored.                                             |
+--------------------+-------------------------------------------------------------------------------------------------------------------------+
| -i                 |   specify the name of the input file to be used (default ``idefix.ini``)                                                |
+--------------------+-------------------------------------------------------------------------------------------------------------------------+

In addition to these specific *Idefix* options, several Kokkos Options can be used on the command
line when running *Idefix*:

+--------------------------+-------------------------------------------------------------------------------------------------------------------+
| Option name              | Comment                                                                                                           |
+==========================+===================================================================================================================+
| --kokkos-num-devices=x   | | Specify the number of devices (eg CUDA GPU) Kokkos should expect. This option is useful when each MPI           |
|                          | | process should be attached to a different GPU. This option replace --kokkos-ndevices which is now deprecated    |
+--------------------------+-------------------------------------------------------------------------------------------------------------------+


.. _signalHandling:

Signal Handling
---------------

By default, *Idefix* is designed to capture the ``SIGUSR2`` UNIX signal sent by the host operating system when it is running. When such a signal is captured, *Idefix* finishes
its current integration step, dumps a restart file and ends. This signal handling is therefore useful to tell *Idefix* that its allocated time for the current
job is ending, and it is therefore time to stop the integration. Many modern jobs schedulers, such as OAR and SLURM, allow users to send UNIX signals
before killing the jobs. Read the documentation of your job scheduler for more information.

If you are not using a scheduler, it is possible to interupt *Idefix* by directly sending it a ``SIGUSR2`` UNIX
signal with the ``kill`` command. For instance:

.. code-block:: bash

  kill -s SIGUSR2  123456

where 123456 is *Idefix* pid (obtainable with ``ps -a``). When runnning with MPI, the pid used in the ``kill`` command can be any *Idefix* process, as *Idefix*
automatically broadcast the abort signal to all of its running processes.
