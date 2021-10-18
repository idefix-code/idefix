# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 22:07:47 2020

@author: lesurg
"""
import struct
import numpy as np


class dataStruct:
    pass

def readHeader(fileHandler):
    headerSize = 128
    q=fileHandler.read(headerSize)
    n=q.index(b'\x00')
    header=q[:n].decode('utf-8')
    return(header)

def readEntry(fileHandler):
    intSize = 4
    nameSize = 16
    dblSize = 8

    byteorder='little'
    data=dataStruct()
    # read entry name
    q=fileHandler.read(nameSize)
    # cut it at the first 0 (cstring format)
    n=q.index(b'\x00')
    data.name=q[:n].decode('utf-8')
    #read datatype, and allocate its size
    data.type=int.from_bytes(fileHandler.read(intSize),byteorder)
    if(data.type==0):
        mysize=dblSize
        stringchar="d"
    if(data.type==1):
        mysize=dblSize
        stringchar="f"
    if(data.type==2):
        mysize=intSize
        stringchar="i"
    data.ndims=int.from_bytes(fileHandler.read(intSize), byteorder)
    dims=[]
    ntot=1
    for dim in range(data.ndims):
        dims.append(int.from_bytes(fileHandler.read(intSize), byteorder))
        ntot=ntot*dims[-1]
    raw=struct.unpack(str(ntot)+stringchar,fileHandler.read(mysize*ntot))
    data.array=(np.asarray(raw).reshape(dims[::-1])).T
    return data

def readDump(filename):

    f=open(filename,"rb")
    data=dataStruct()
    # read the header
    data.header=readHeader(f)

    #read coordinates
    field=readEntry(f)
    data.x1=field.array
    field=readEntry(f)
    data.x2=field.array
    field=readEntry(f)
    data.x3=field.array

    # read remaining fields and store them
    data.data={}
    while(field.name!="eof"):

        field=readEntry(f)
        data.data[field.name]=field.array

    f.close()
    return(data)
