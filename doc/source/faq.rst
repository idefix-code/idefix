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

I want to run on the GPUs of xxx machine, how do I proceed?
  Check the examples in :ref:`setupExamples`

Compilation
-----------

Is there a way to see explicitely the compilation commands with ``make``?
  Yes, just add ``VERBOSE=1`` after the ``make`` command.

The compilation stops while compiling Kokkos with ``/usr/include/stdlib.h(58): error: expected a ";"``
  This happens on Gricad machines when LIBDL is activated (wrong glibc). Simply disable LIBDL.

When using the Intel compiler on a Mac Intel, I get a linking error involving the ``SharedAllocationRecordIvvE18t_tracking_enabledE`` symbol.
  This is a known bug of the Intel Mac compiler with Kokkos. Apparently Intel has decided not to fix it. Check the issue on the `Kokkos git page <https://github.com/kokkos/kokkos/issues/1959>`_.

I get an error during *Idefix* link that says there are undefined symbols to std::filesystem
  This is a known bug/limitation of the stdlibc++ provided with gcc8 that does not include C++17 filesystem extensions.
  While *Idefix* auto-detects gcc8 when it is used as the main compiler, it still misses the cases when another compiler
  (like Clang or Intel) is used with gcc8 as a backend.
  You should try to clear up CMakeCache.txt and explicitely add the required link library when calling cmake as in
  ``LDFLAGS=-lstdc++fs cmake $IDEFIX_DIR ...``

At the end of the compilation phase, during link on MacOS, I get an error ``ld: Assertion failed: (resultIndex < sectData.atoms.size()), function findAtom, ...``.
  This is a known bug of the new linker provided by Apple with Xcode 15. Revert to the old linker:
  ``LDFLAGS=-ld_classic cmake $IDEFIX_DIR ...``

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
  Add ``-DIdefix_DEBUG=ON`` to ``cmake`` and recompile to find out exactly where the code crashes (see :ref:`debugging`).

I want to test a modification to *Idefix* source code specific to my problem without modifying files in `$IDEFIX_DIR`. How can I proceed?
  Add a ``CmakeLists.txt`` in your problem directory and use the function `replace_idefix_source` (see :ref:`customSourceFiles`).

I want to use a lookup table from a CSV file in my idefix_loop. How could I proceed?
  Use the ``LookupTable`` class (see :ref:`LookupTableClass`)
