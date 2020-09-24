# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 15:23:00 2014

@author: glesur
"""


import numpy as np
class DataStructure:
    pass

# Read a vtk file
def readVTKCart(filename):
    try:
        fid=open(filename,"rb")
    except:
        print("Can't open file")
        return 0

    # define our datastructure
    V=DataStructure()

    # raw data which will be read from the file
    V.data={}

    #print("Hello")
    # datatype we read
    dt=np.dtype(">f")   # Big endian single precision floats

    s=fid.readline()    # VTK DataFile Version x.x
    s=fid.readline()    # Comments

    s=fid.readline()    # BINARY
    s=fid.readline()    # DATASET RECTILINEAR_GRID
    slist=s.split()
    grid_type=str(slist[1],'utf-8')
    if(grid_type != "RECTILINEAR_GRID"):
        print("ERROR: Wrong VTK file type.")
        print("This routine can only open Cartesian or Cylindrical VTK files.")
        fid.close()
        return 0

    s=fid.readline()    # DIMENSIONS NX NY NZ
    slist=s.split()
    #s=fid.readline()    # Extre line feed
    V.nx=int(slist[1])
    V.ny=int(slist[2])
    V.nz=int(slist[3])

    s=fid.readline()    # X_COORDINATES NX float

    x=np.fromfile(fid,dt,V.nx)
    s=fid.readline()    # Extra line feed added by pluto


    s=fid.readline()    # X_COORDINATES NX float

    y=np.fromfile(fid,dt,V.ny)
    s=fid.readline()    # Extra line feed added by pluto

    s=fid.readline()    # X_COORDINATES NX float

    z=np.fromfile(fid,dt,V.nz)
    s=fid.readline()    # Extra line feed added by pluto



    s=fid.readline()    # POINT_DATA NXNYNZ

    slist=s.split()
    point_type=str(slist[0],'utf-8')
    npoints=int(slist[1])
    s=fid.readline()    # EXTRA LINE FEED

    if(point_type == "CELL_DATA"):
        # The file contains face coordinates, so we extrapolate to get the cell center coordinates.
        if V.nx>1:
            V.nx=V.nx-1
            V.x=0.5*(x[1:]+x[:-1])
        else:
            V.x=x
        if V.ny>1:
            V.ny=V.ny-1
            V.y=0.5*(y[1:]+y[:-1])
        else:
            V.y=y
        if V.nz>1:
            V.nz=V.nz-1
            V.z=0.5*(z[1:]+z[:-1])
        else:
            V.z=z
    elif(point_type == "POINT_DATA"):
        V.x=x
        V.y=y
        V.z=z

    if V.nx*V.ny*V.nz != npoints:
        print("ERROR: Grid size incompatible with number of points in the data set")

    while 1:
        s=fid.readline()        # SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
        #print repr(s)
        if len(s)<2:         # leave if end of file
            break
        slist=s.split()
        datatype=str(slist[0],'utf-8')
        varname=str(slist[1],'utf-8')
        if datatype == "SCALARS":
            fid.readline()  # LOOKUP TABLE
            V.data[varname] = np.transpose(np.fromfile(fid,dt,V.nx*V.ny*V.nz).reshape(V.nz,V.ny,V.nx))
        elif datatype == "VECTORS":
            Q=np.fromfile(fid,dt,3*V.nx*V.ny*V.nz)

            V.data[varname+'_X']=np.transpose(Q[::3].reshape(V.nz,V.ny,V.nx))
            V.data[varname+'_Y']=np.transpose(Q[1::3].reshape(V.nz,V.ny,V.nx))
            V.data[varname+'_Z']=np.transpose(Q[2::3].reshape(V.nz,V.ny,V.nx))

        else:
            print("ERROR: Unknown datatype %s" % datatype)
            break;

        fid.readline()  #extra line feed
    fid.close()

    return V

# Read a vtk file
def readVTKPolar(filename):
    try:
        fid=open(filename,"rb")
    except:
        print("Can't open file")
        return 0

    # define our datastructure
    V=DataStructure()

    # raw data which will be read from the file
    V.data={}

    #print("Hello")
    # datatype we read
    dt=np.dtype(">f")   # Big endian single precision floats

    s=fid.readline()    # VTK DataFile Version x.x
    s=fid.readline()    # Comments

    s=fid.readline()    # BINARY
    s=fid.readline()    # DATASET RECTILINEAR_GRID
    print(s)
    slist=s.split()
    grid_type=str(slist[1],'utf-8')
    if(grid_type != "STRUCTURED_GRID"):
        print("ERROR: Wrong VTK file type.")
        print("Current type is: %s"%(grid_type))
        print("This routine can only open Polar VTK files.")
        fid.close()
        return 0

    s=fid.readline()    # DIMENSIONS NX NY NZ
    slist=s.split()
    V.nx=int(slist[1])
    V.ny=int(slist[2])
    V.nz=int(slist[3])

    print("nx=%d, ny=%d, nz=%d"%(V.nx,V.ny,V.nz))

    s=fid.readline()    # POINTS NXNYNZ float
    slist=s.split()
    npoints=int(slist[1])
    points=np.fromfile(fid,dt,3*npoints)
    s=fid.readline()    # EXTRA LINE FEED

    V.points=points

    if V.nx*V.ny*V.nz != npoints:
        print("ERROR: Grid size incompatible with number of points in the data set")
        return 0

    # Reconstruct the polar coordinate system
    x1d=points[::3]
    y1d=points[1::3]
    z1d=points[2::3]

    xcart=np.transpose(x1d.reshape(V.nz,V.ny,V.nx))
    ycart=np.transpose(y1d.reshape(V.nz,V.ny,V.nx))
    zcart=np.transpose(z1d.reshape(V.nz,V.ny,V.nx))

    r=np.sqrt(xcart[:,0,0]**2+ycart[:,0,0]**2)
    theta=np.unwrap(np.arctan2(ycart[0,:,0],xcart[0,:,0]))
    z=zcart[0,0,:]

    s=fid.readline()    # CELL_DATA (NX-1)(NY-1)(NZ-1)
    slist=s.split()
    data_type=str(slist[0],'utf-8')
    if(data_type != "CELL_DATA"):
        print("ERROR: this routine expect CELL DATA as produced by PLUTO.")
        fid.close()
        return 0
    s=fid.readline()    # Line feed
    # Perform averaging on coordinate system to get cell centers
    # The file contains face coordinates, so we extrapolate to get the cell center coordinates.
    if V.nx>1:
        V.nx=V.nx-1
        V.x=0.5*(r[1:]+r[:-1])
    else:
        V.x=r
    if V.ny>1:
        V.ny=V.ny-1
        V.y=(0.5*(theta[1:]+theta[:-1])+np.pi)%(2.0*np.pi)-np.pi
    else:
        V.y=theta
    if V.nz>1:
        V.nz=V.nz-1
        V.z=0.5*(z[1:]+z[:-1])
    else:
        V.z=z


    while 1:
        s=fid.readline()        # SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
        #print repr(s)
        if len(s)<2:         # leave if end of file
            break
        slist=s.split()
        datatype=str(slist[0],'utf-8')
        varname=str(slist[1],'utf-8')
        if datatype == "SCALARS":
            fid.readline()  # LOOKUP TABLE
            V.data[varname] = np.transpose(np.fromfile(fid,dt,V.nx*V.ny*V.nz).reshape(V.nz,V.ny,V.nx))
        elif datatype == "VECTORS":
            Q=np.fromfile(fid,dt,3*V.nx*V.ny*V.nz)

            V.data[varname+'_X']=np.transpose(Q[::3].reshape(V.nz,V.ny,V.nx))
            V.data[varname+'_Y']=np.transpose(Q[1::3].reshape(V.nz,V.ny,V.nx))
            V.data[varname+'_Z']=np.transpose(Q[2::3].reshape(V.nz,V.ny,V.nx))

        else:
            print("ERROR: Unknown datatype %s" % datatype)
            break;

        fid.readline()  #extra line feed
    fid.close()

    return V

# Read a vtk file
def readVTKSpherical(filename):
    try:
        fid=open(filename,"rb")
    except:
        print("Can't open file")
        return 0

    # define our datastructure
    V=DataStructure()

    # raw data which will be read from the file
    V.data={}

    #print("Hello")
    # datatype we read
    dt=np.dtype(">f")   # Big endian single precision floats

    s=fid.readline()    # VTK DataFile Version x.x
    s=fid.readline()    # Comments

    s=fid.readline()    # BINARY
    s=fid.readline()    # DATASET RECTILINEAR_GRID
    slist=s.split()
    grid_type=str(slist[1],'utf-8')
    if(grid_type != "STRUCTURED_GRID"):
        print("ERROR: Wrong VTK file type.")
        print("This routine can only open Spherical VTK files.")
        fid.close()
        return 0

    s=fid.readline()    # DIMENSIONS NX NY NZ
    slist=s.split()
    V.nx=int(slist[1])
    V.ny=int(slist[2])
    V.nz=int(slist[3])

    if(V.nz==1):
        is2d=1
    else:
        is2d=0

    s=fid.readline()    # POINTS NXNYNZ float
    slist=s.split()
    npoints=int(slist[1])
    points=np.fromfile(fid,dt,3*npoints)
    s=fid.readline()    # EXTRA LINE FEED

    V.points=points

    if V.nx*V.ny*V.nz != npoints:
        print("ERROR: Grid size incompatible with number of points in the data set")
        return 0

    # Reconstruct the spherical coordinate system


    x1d=points[::3]
    y1d=points[1::3]
    z1d=points[2::3]

    xcart=np.transpose(x1d.reshape(V.nz,V.ny,V.nx))
    ycart=np.transpose(y1d.reshape(V.nz,V.ny,V.nx))
    zcart=np.transpose(z1d.reshape(V.nz,V.ny,V.nx))

    if(is2d):
        r=np.sqrt(xcart[:,0,0]**2+ycart[:,0,0]**2)
        phi=np.unwrap(np.arctan2(zcart[0,0,:],xcart[0,0,:]))
        theta=np.arccos(ycart[0,:,0]/np.sqrt(xcart[0,:,0]**2+ycart[0,:,0]**2))
    else:
        r=np.sqrt(xcart[:,0,0]**2+ycart[:,0,0]**2+zcart[:,0,0]**2)
        phi=np.unwrap(np.arctan2(ycart[0,0,:],xcart[0,0,:]))
        theta=np.arccos(zcart[0,:,0]/np.sqrt(xcart[0,:,0]**2+ycart[0,:,0]**2+zcart[0,:,0]**2))


    s=fid.readline()    # CELL_DATA (NX-1)(NY-1)(NZ-1)
    slist=s.split()
    data_type=str(slist[0],'utf-8')
    if(data_type != "CELL_DATA"):
        print("ERROR: this routine expect CELL DATA as produced by PLUTO.")
        fid.close()
        return 0
    s=fid.readline()    # Line feed
    # Perform averaging on coordinate system to get cell centers
    # The file contains face coordinates, so we extrapolate to get the cell center coordinates.
    if V.nx>1:
        V.nx=V.nx-1
        V.r=0.5*(r[1:]+r[:-1])
    else:
        V.x=r
    if V.ny>1:
        V.ny=V.ny-1
        V.theta=0.5*(theta[1:]+theta[:-1])
    else:
        V.y=theta
    if V.nz>1:
        V.nz=V.nz-1
        V.phi=(0.5*(phi[1:]+phi[:-1])+np.pi)%(2.0*np.pi)-np.pi
    else:
        V.phi=phi


    while 1:
        s=fid.readline()        # SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
        #print repr(s)
        if len(s)<2:         # leave if end of file
            break
        slist=s.split()
        datatype=str(slist[0],'utf-8')
        varname=str(slist[1],'utf-8')
        if datatype == "SCALARS":
            fid.readline()  # LOOKUP TABLE
            V.data[varname] = np.transpose(np.fromfile(fid,dt,V.nx*V.ny*V.nz).reshape(V.nz,V.ny,V.nx))
        elif datatype == "VECTORS":
            Q=np.fromfile(fid,dt,3*V.nx*V.ny*V.nz)

            V.data[varname+'_X']=np.transpose(Q[::3].reshape(V.nz,V.ny,V.nx))
            V.data[varname+'_Y']=np.transpose(Q[1::3].reshape(V.nz,V.ny,V.nx))
            V.data[varname+'_Z']=np.transpose(Q[2::3].reshape(V.nz,V.ny,V.nx))

        else:
            print("ERROR: Unknown datatype %s" % datatype)
            break;

        fid.readline()  #extra line feed
    fid.close()

    return V
