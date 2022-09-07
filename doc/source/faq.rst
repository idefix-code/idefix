.. _faq:

=========================
Frequently asked question
=========================

Configuration with Cmake
------------------------

Cmake fails with "CMake wants to use -std=c++1z which is not supported by NVCC" when configuring with Cuda
  This error happens when cmake detects an old version of gcc (<=7). These old version uses
  ``-std=c++1z`` to enable c++17, but this flag is not recognised by Cuda ``nvcc`` compiler. Use
  a more recent version of gcc (>=8).

Cmake fails with "C++17-compliant compiler detected, but unable to compile C++14 or later program.  Verify that Intel:xxxxxxxxx is set up correctly" when using intel compilers
  The error arises with some error on the "std" namespace,. This is because the intel compiler
  tries to use the default GNU C++ Headers, and these are too old. Try to install a more recent
  release of the gcc compilation suite.

How do I use my favourite XXXXX compiler?
  set the ``CXX`` variables to your favourite C++ compiler before calling Cmake (see :ref:`configurationOptions`).
  Do not forget to delete any CmakeCache.txt which could have been produced in your problem directory.

I have a complex setup, and have written some functions in separate .cpp files. How can I add these files to *Idefix* build tree?
  Add a ``CmakeLists.txt`` in your problem directory and use the function `add_idefix_source` (see :ref:`customSourceFiles`).

Compilation
-----------

Is there a way to see explicitely the compilation commands with ``make``?
  Yes, just add ``VERBOSE=1`` after the ``make`` command.

The compilation stops while compiling Kokkos with ``/usr/include/stdlib.h(58): error: expected a ";"``
  This happens on Gricad machines when LIBDL is activated (wrong glibc). Simply disable LIBDL.

Execution
---------

How can I stop the code without loosing the current calculation?
  Two options: the simplest one is to create an empty file named `stop` in the running directory
  (e.g. with ``touch stop``). It is also possible to send the POSIX signal SIGUSR2 to one of the
  *idefix* processes. More details here: :ref:`signalHandling`.

I'm doing performance measures. How do I disable all outputs in *Idefix*?
  Add ``-nowrite`` when you call *Idefix* executable.

Developement
------------

I have a serious bug (e.g. segmentation fault), in my setup, how do I proceed?
  Add ``-DIdefix_DEBUG=ON`` to ``cmake`` and recompile to find out exactly where the code crashes.

I want to test a modification to *Idefix* source code specific to my problem without modifying files in `$IDEFIX_DIR`. How can I proceed?
  Add a ``CmakeLists.txt`` in your problem directory and use the function `replace_idefix_source` (see :ref:`customSourceFiles`).

I want to use a lookup table from a CSV file in my idefix_loop. How could I proceed?
  Use the ``LookupTable`` class (see :ref:`LookupTableClass`)
