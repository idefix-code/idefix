
.. _faq:

=========================
Frequently asked question
=========================

Configuration with Cmake
------------------------

How do I compile and link problem-specific source file with *Idefix*?
  Use a custom ``CMakeLists.txt`` in your problem directory as described in :ref:`customSourceFiles`.

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

Compilation
-----------

Coming soon :-)
