======================
Performances
======================

We report below the performances obtained on various architectures using Idefix. The reference test
is the 3D MHD Orszag-Tang test problem with 2nd order reconstruction and uct_contact EMFS bundled in
Idefix test suite, disabling passive tracers. The test is computed with a 128\ :sup:`3` resolution per
MPI sub-domain on GPUs or 32\ :sup:`3` per MPI sub-domain on CPUs. All of the performances measures
have been obtained enabling MPI and we reporte here the performance *per GPU*, *per GCD* (on Mi250)
or *per core* (on CPU).

The complete scalability tests are available in Idefix `method paper <https://ui.adsabs.harvard.edu/abs/2023A%26A...677A...9L/abstract>`_.
The performances mentionned below are updated for each major revision of Idefix, so they might slightly differ from the method paper.

.. note::

    You might expect
    slower performances with lower resolution when using GPUs. The overall performances also depends on
    the physical modules activated, the reconstruction scheme, and the efficiency of the parallel network
    on which you are running. The performances reported below are therefore purely indicative. We encourage
    you to use the embedded profiler (see :ref:`commandLine` ) when performances are smaller than expected.


CPU performances
================

+---------------------+--------------------+----------------------------------------------------+
| Cluster name        | Processor          | Performances (in 10\ :sup:`6` cell/s/core)         |
+=====================+====================+====================================================+
| TGCC/Irene Rome     | AMD EPYC Rome      | 0.29                                               |
+---------------------+--------------------+----------------------------------------------------+
| IDRIS/Jean Zay      | Intel Cascade Lake | 0.62                                               |
+---------------------+--------------------+----------------------------------------------------+

GPU performances
================

.. plot::

   import plot_idefix_bench
   plot_idefix_bench.do_plot('Performance on NVidia and AMD GPUs', 'bench.json', ['v100','a100','h100','mi250x'])

.. note::

    The inter-node communication on Jean Zay is not optimal on A100 nodes. A ticket is opened with IDRIS support to fix this issue.
