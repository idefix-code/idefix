import numpy as np
import struct
# Read .idfx files which are created
# for debug purposes by DataBlock.DumpToFile(std::string&)

# There is one .idfx file per processor

def readIdfxHeader(fileHandler):
    headerSize = 128
    q=fileHandler.read(headerSize)
    n=q.index(b'\x00')
    header=q[:n].decode('utf-8')
    return(header)

def readIdfxEntry(fileHandler):
    class dataStruct:
        pass
    byteorder='little'
    nameSize = 16
    intSize = 4
    dblSize = 8

    data=dataStruct()
    # read entry name
    q=fileHandler.read(nameSize)
    # cut it at the first 0 (cstring format)
    n=q.index(b'\x00')
    data.name=q[:n].decode('utf-8')

    #read dimensions
    data.ndims=int.from_bytes(fileHandler.read(intSize), byteorder)
    dims=[]
    ntot=1
    for dim in range(data.ndims):
        dims.append(int.from_bytes(fileHandler.read(intSize), byteorder))
        ntot=ntot*dims[-1]
    raw=struct.unpack(str(ntot)+"d",fileHandler.read(dblSize*ntot))
    data.array=np.asarray(raw).reshape(dims[::-1])
    return data

def readIdfxFile(filename):
    f=open(filename,"rb")
    class dataStruct:
        pass

    data=dataStruct()
    # read the header
    data.header=readIdfxHeader(f)

    # read remaining fields and store them
    data.V={}
    field=readIdfxEntry(f)
    data.V[field.name]=field.array

    while(field.name!="eof"):

        field=readIdfxEntry(f)
        data.V[field.name]=field.array

    f.close()
    return(data)
