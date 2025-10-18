Continuous Integration (CI) tests
================================

This document describes the GitHub Actions continuous-integration setup used to run the Idefix
test-suite. The CI is implemented by two workflows checked in .github/workflows:

- .github/workflows/idefix-ci.yml
- .github/workflows/idefix-ci-jobs.yml

Overview
--------

The CI is split in two layers:

- A top-level workflow (.github/workflows/idefix-ci.yml) that:

  - runs a Linter job (pre-commit) on push / PR / manual dispatch,
  - then calls a reusable workflow for different compiler/backends (intel, gcc, cuda)
    providing two inputs: TESTME_OPTIONS and IDEFIX_COMPILER.

- A reusable workflow (.github/workflows/idefix-ci-jobs.yml) that:

  - defines the actual test jobs grouped by physics domain (ShocksHydro, ParabolicHydro,
    ShocksMHD, ParabolicMHD, Fargo, ShearingBox, SelfGravity, Planet, Dust, Braginskii,
    Examples, Utils),
  - runs test scripts on self-hosted runners,
  - expects the repository to be checked out with submodules,
  - invokes the repository-provided CI helper scripts to configure / build / run tests.

Key configuration points
------------------------

- Inputs passed from the top-level workflow:

  - TESTME_OPTIONS (string): flags forwarded to the per-test runner (examples: -cuda, -Werror,
    -intel, -all).
  - IDEFIX_COMPILER (string): which compiler the tests should use (e.g. icc, gcc, nvcc).

- Environment variables set by the reusable workflow:

  - IDEFIX_COMPILER, TESTME_OPTIONS, PYTHONPATH, IDEFIX_DIR

- Linter job:

  - Runs only when repository is the main project (not arbitrary forks).
  - Uses actions/setup-python and runs pre-commit (pre-commit/action@v3 and pre-commit-ci/lite).
  - Prevents regressions in style and common mistakes before running heavy test jobs.

- Test execution:

  - All test jobs call the repository script scripts/ci/run-tests with a test directory
    and the TESTME_OPTIONS flags. Example invocation (from the workflows):
      scripts/ci/run-tests $IDEFIX_DIR/test/HD/sod -all $TESTME_OPTIONS

  - The reusable workflow is written to execute many test directories in separate job steps,
    so each physics group is kept logically separated in CI logs.

Runners and prerequisites
-------------------------

- The heavy numerical tests run on self-hosted runners (see runs-on: self-hosted).
  The CI assumes appropriate hardware and dependencies are available on those runners
  (compilers, MPI, GPUs when CUDA/HIP flags are used, required system libraries).

- The workflows check out the repository and its submodules. Submodules must be available
  on the CI machines.

How tests are driven (testme scripts)
-------------------------------------

Each test directory contains a small Python "testMe" driver that uses the helper Python
class documented in the repository:

- See the test helper documentation: :doc:`idfxTest <testing/idfxTest>`

That helper (idfxTest) is responsible for:

- parsing TESTME_OPTIONS-like flags (precision, MPI, CUDA, reconstruction, vector potential, etc.),
- calling configure / compile / run,
- performing standard python checks and non-regression (RMSE) comparisons against
  reference dumps,
- optionally creating / updating reference dumps (init mode).

Practical examples
------------------

- Example of a CI invocation (triggered by workflows):

  - Top-level workflow calls the reusable jobs workflow for each compiler/back-end, e.g.
    TESTME_OPTIONS="-cuda -Werror" IDEFIX_COMPILER=nvcc

- Running tests locally (developer machine)
  - You can mimic what CI does by calling the repository helper script directly. Example:
    scripts/ci/run-tests /path/to/idefix/test/HD/sod -all -mpi -dec 2 2 -reconstruction 3 -single

Notes for maintainers
---------------------

- The reusable jobs workflow contains a commented concurrency block for optional cancellation
  of in-flight runs — consider enabling it if you want to auto-cancel redundant CI runs.
- Because tests are run on self-hosted runners, ensure the pools have the required compilers,
  MPI stacks and GPU drivers for the requested TESTME_OPTIONS.
- Keep TESTME_OPTIONS in sync with the options understood by the test helper documented in
  :doc:`idfxTest <testing/idfxTest>`.

Relevant files
--------------

- Workflow entry point: .github/workflows/idefix-ci.yml
- Reusable jobs: .github/workflows/idefix-ci-jobs.yml
- Test helper documentation: :doc:`idfxTest <testing/idfxTest>`

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   testing/idfxTest.rst