import struct

import numpy as np


__all__ = ["readIdfxFile"]
# Read .idfx files which are created
# for debug purposes by DataBlock.DumpToFile(std::string&)

INT_SIZE = 4
NAME_SIZE = 16
DOUBLE_SIZE = 8
FLOAT_SIZE = 4

HEADER_SIZE = 128

# There is one .idfx file per processor
class IdfxFileField(object):

    def __init__(self, fh, byteorder="little"):
        # read entry name
        q = fh.read(NAME_SIZE)
        # cut it at the first 0 (cstring format)
        n = q.index(b"\x00")
        self.name = q[:n].decode("utf-8")
        if self.name == "eof":
            return
        self.ndims = int.from_bytes(fh.read(INT_SIZE), byteorder)
        dims = []
        for dim in range(self.ndims):
            dims.append(int.from_bytes(fh.read(INT_SIZE), byteorder))
        ntot = int(np.prod(dims))
        raw = struct.unpack(str(ntot) + "d", fh.read(DOUBLE_SIZE * ntot))
        self.array = np.asarray(raw).reshape(dims[::-1])

class IdfxFileDataset(object):
    def __init__(self, filename):
        # identical to DumpDataset
        self.filename = filename
        self.metadata = {}
        with open(filename, "rb") as fh:
            self._read_header(fh)
            self._read_fields(fh)

    def _read_header(self, fh):
        # could easily be identical to DumpDataset
        headerSize = 128
        q = fh.read(headerSize)
        n = q.index(b'\x00')
        self.header = q[:n].decode('utf-8')

        self.metadata["byteorder"] = "little"


    def _read_field(self, fh):
        # identical to DumpDataset
        if self.metadata["byteorder"] is None:
            # "little" is a safe bet. If anyone ever *needs* to analyze big-endian data produced
            # with old versions of Idefix, then we could offer some flexibility here.
            byteorder = "little"
        else:
            byteorder = self.metadata["byteorder"]

        return IdfxFileField(fh, byteorder)

    def _read_fields(self, fh):
        self.data = {}
        while True:
            field = self._read_field(fh)
            if field.name == "eof":
                break
            self.data[field.name] = field.array

    def __repr__(self):
        return "IdfxFileDataset('%s')" % self.filename


# public API
def readIdfxFile(filename):
    return IdfxFileDataset(filename)
