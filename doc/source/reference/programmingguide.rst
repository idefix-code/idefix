
======================
Programming guide
======================

Because *Idefix* is designed to run on hybrid architecture, there are some subtilities which
are not found in other classical CPU-only codes. These subtilities are described in this section, along
with some basic coding rules which can be handful.

Host and target
===============

*Idefix* relies on the Kokkos framework, and therefore assume that system it's running is made
of two sub-systems: a host and a target. The host is traditionnaly the CPU, and is taking care
of inputs and outputs, initialisation and allocation, MPI data exchanges. The target is usually an
accelerator (e.g. a GPU) and is actually performing the computation (or most of it).

Note that while *Idefix* assumes there is a host and a target, the target can be the host (think
of the code running only on your laptop CPU). In this case, Kokkos performs several optimisations,
so that everything effectively runs on the host smoothly.

By construction, the host doesn't have direct access to the target memory and vice-versa. This means
that accessing an array on the target from the host will inevitably lead to a segmentation fault.
This is the most common mistake, so keep this in mind.

Arrays
======



Execution space and loops
=========================
