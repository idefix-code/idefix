======================
Performances
======================

We report below the performances obtained on various architectures using Idefix. The reference test
is the 3D MHD Orszag-Tang test problem with 2nd order reconstruction and uct_contact EMFS bundled in
Idefix test suite, disabling passive tracers. The test is computed with a 128\ :sup:`3` resolution per
MPI sub-domain on GPUs or 32\ :sup:`3` per MPI sub-domain on CPUs. All of the performances measures
have been obtained enabling MPI on *one full node*, but we report here the performance *per GPU*
(i.e. with 2 GCDs on AMD Mi250) or *per core* (on CPU), i.e. dividing the node performance by the number of GPU/core
to simplify the comparison with other clusters.

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

+----------------------+--------------------+----------------------------------------------------+
| Cluster name         | GPU                | Performances (in 10\ :sup:`6` cell/s/GPU)          |
+======================+====================+====================================================+
| IDRIS/Jean Zay       | NVIDIA V100        | 110                                                |
+----------------------+--------------------+----------------------------------------------------+
| IDRIS/Jean Zay       | NVIDIA A100        | 194                                                |
+----------------------+--------------------+----------------------------------------------------+
| CINES/Adastra        | AMD Mi250          | 250                                                |
+----------------------+--------------------+----------------------------------------------------+
