# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 22:07:47 2020

@author: lesurg
"""
import os
import re
import struct

import numpy as np

__all__ = ["readDump"]


header_regexp = re.compile(
    r"Idefix (?P<version>[\w\.-]+) Dump Data"
    r"((?P<byteorder>(little|big)) endian)?"
)
NAME_SIZE = 16

DOUBLE_SIZE = 8
FLOAT_SIZE = 4
INT_SIZE = 4
BOOL_SIZE = 1

HEADER_SIZE = 128


class DumpField(object):
    def __init__(self, fh, byteorder="little"):
        # read entry name
        q = fh.read(NAME_SIZE)
        # cut it at the first 0 (cstring format)
        n = q.index(b"\x00")
        self.name = q[:n].decode("utf-8")
        # read datatype, and allocate its size
        self.type = int.from_bytes(fh.read(INT_SIZE), byteorder)
        if self.type == 0:
            mysize = DOUBLE_SIZE
            stringchar = "d"
            dtype = "float64"
        elif self.type == 1:
            mysize = FLOAT_SIZE
            stringchar = "f"
            dtype = "float32"
        elif self.type == 2:
            mysize = INT_SIZE
            stringchar = "i"
            dtype = "int32"
        elif self.type == 3:
            mysize = BOOL_SIZE
            stringchar = "?"
            dtype = bool
        else:
            raise RuntimeError(
                "Found unknown data type %d for field %s" % (self.type, self.name)
            )
        self.ndims = int.from_bytes(fh.read(INT_SIZE), byteorder)
        dims = []
        ntot = 1
        for dim in range(self.ndims):
            dims.append(int.from_bytes(fh.read(INT_SIZE), byteorder))
            ntot = ntot * dims[-1]
        raw = struct.unpack(str(ntot) + stringchar, fh.read(mysize * ntot))
        self.array = np.asarray(raw, dtype=dtype).reshape(dims[::-1]).T


class DumpDataset(object):
    def __init__(self, filename):
        self.filename = os.path.abspath(filename)
        self.metadata = {}
        with open(filename, "rb") as fh:
            self._read_header(fh)
            self._read_fields(fh)

    def _read_header(self, fh):
        q = fh.read(HEADER_SIZE)
        n = q.index(b"\x00")
        self.header = q[:n].decode("utf-8")

        match = re.match(header_regexp, self.header)
        # note that "version" is the only field that is expected to be always present
        # byteorder is set to None if not found in the header
        self.metadata["version"] = match.group("version")
        self.metadata["byteorder"] = match.group("byteorder")

    def _read_field(self, fh):
        if self.metadata["byteorder"] is None:
            # "little" is a safe bet. If anyone ever *needs* to analyze big-endian data produced
            # with old versions of Idefix, then we could offer some flexibility here.
            byteorder = "little"
        else:
            byteorder = self.metadata["byteorder"]

        return DumpField(fh, byteorder)

    def _read_fields(self, fh):
        # read coordinates

        self.x1 = self._read_field(fh).array
        self.x1l = self._read_field(fh).array
        self.x1r = self._read_field(fh).array
        self.x2 = self._read_field(fh).array
        self.x2l = self._read_field(fh).array
        self.x2r = self._read_field(fh).array
        self.x3 = self._read_field(fh).array
        self.x3l = self._read_field(fh).array
        self.x3r = self._read_field(fh).array

        # read remaining fields and store them
        self.data = {}
        while True:
            field = self._read_field(fh)
            if field.name == "eof":
                break
            self.data[field.name] = field.array

    def __repr__(self):
        return "DumpDataset('%s')" % self.filename

# public API
def readDump(filename):
    return DumpDataset(filename)
