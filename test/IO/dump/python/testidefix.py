#!/usr/bin/env python3
######################################################################################
# Idefix MHD astrophysical code
#
# Source file test/IO/dump/testidefix.py
#
# Last modified : 11/2026
#
# Copyright(C) by :
# - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
# and other code contributors
#
# Licensed under CeCILL 2.1 License, see COPYING for more information
######################################################################################

"""
Independent check of the dump data type support.

This script reads a dump file produced by test/IO/dump using the pure-Python
reader (pytools/dump_io.py) and verifies that every data type supported by the
dump I/O (int, int64, float, double, bool and complex float/double, both
std::complex and Kokkos::complex) is read back with the expected value.

The expected values are the initial values set in test/IO/dump/setup.cpp (in
the `DumpScalars` structure), which are the values written by the first dump
(before the values are randomised and reloaded for the restart round-trip
checked inside Analysis()).

Usage:
    python3 testidefix.py [dumpfile]

If no dump file is given, dump.0000.dmp in the current directory is used.
"""

import os
import sys

import numpy as np

sys.path.append(os.getenv("IDEFIX_DIR"))

from pytools.dump_io import readDump  # noqa: E402

# Expected scalar values, mirroring the initial values in setup.cpp (DumpScalars).
# key -> (expected value, expected numpy dtype)
EXPECTED = {
    "testInt": (42, np.int32),
    "testInt64": (-1234567890123, np.int64),
    "testFloat": (np.float32(1.25), np.float32),
    "testDouble": (3.14159265358979, np.float64),
    "testBool": (True, np.bool_),
    "testCFloat": (np.complex64(complex(0.5, -3.0)), np.complex64),
    "testCDouble": (np.complex128(complex(1.5, -2.25)), np.complex128),
    "testKCFloat": (np.complex64(complex(-7.0, 8.5)), np.complex64),
    "testKCDouble": (np.complex128(complex(4.0, 5.5)), np.complex128),
}


def check_dump(filename, expected=EXPECTED, rtol=0.0, atol=0.0):
    """Read `filename` and check every field in `expected`.

    Returns the number of failures (0 means success).
    """
    dataset = readDump(filename)
    errors = 0
    for name, (value, dtype) in expected.items():
        if name not in dataset.data:
            print(f"FAIL: field {name!r} not found in {filename}")
            errors += 1
            continue
        field = dataset.data[name]
        if field.dtype != dtype:
            print(f"FAIL: field {name!r} has dtype {field.dtype}, expected {dtype}")
            errors += 1
            continue
        got = field.ravel()[0]
        if isinstance(value, (bool, np.bool_)):
            ok = bool(got) == bool(value)
        else:
            ok = np.isclose(got, value, rtol=rtol, atol=atol)
        if not ok:
            print(f"FAIL: field {name!r} = {got!r}, expected {value!r}")
            errors += 1
        else:
            print(f"OK:   field {name!r} = {got!r}")
    return errors


def validate_dump():
    """Entry point used as a test.py 'callPyFunctionAfter' hook.

    Runs in the test directory; raises if any field is wrong so that the
    test suite reports a failure.
    """
    filename = "../dump.0001.dmp"
    if not os.path.exists(filename):
        raise RuntimeError(f"dump file {filename!r} not found")
    errors = check_dump(filename)
    if errors:
        raise RuntimeError(f"{errors} dump field(s) failed validation in {filename}")
    print(f"All dump data types validated successfully in {filename}")


if __name__ == "__main__":
    validate_dump()
    sys.exit(0)
