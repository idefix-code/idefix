"""
Created on Mon Nov  3 15:23:00 2014

@author: glesur
"""


import numpy as np
class DataStructure:
    pass

# Read a VTK file for any geometry
def readVTK(filename, geometry="unknown"):
    try:
        fid=open(filename,"rb")
    except:
        print("Can't open file")
        return 0

    # define our datastructure
    V=DataStructure()

    # raw data which will be read from the file
    V.data={}

    # initialize geometry
    if geometry not in ("unknown", "cartesian", "polar", "spherical"):
      raise ValueError("Received unknown geometry: '{}'.".format(geometry))
    V.geometry = geometry

    # datatype we read
    dt=np.dtype(">f")   # Big endian single precision floats
    dint=np.dtype(">i4")   # Big endian integer

    s=fid.readline()    # VTK DataFile Version x.x
    s=fid.readline()    # Comments

    s=fid.readline()    # BINARY
    s=fid.readline()    # DATASET RECTILINEAR_GRID or STRUCTURED_GRID
    slist=s.split()

    s=fid.readline()    # DIMENSIONS NX NY NZ or FIELD
    slist=s.split()
    entry=str(slist[0].decode('utf-8'))
    if(entry == "FIELD"):
      nfield=int(slist[2])
      for field in range(nfield):
        s=fid.readline()
        slist=s.split()
        entry=str(slist[0].decode('utf-8'))
        if(entry == "TIME"):
          V.t=np.fromfile(fid,dt,1)
        elif(entry == "GEOMETRY"):
          g=np.fromfile(fid,dint,1)
          if g == 0:
            thisgeometry="cartesian"
          elif g == 1:
            thisgeometry="polar"
          elif g == 2:
            thisgeometry="spherical"
          elif g == 3:
            thisgeometry="cylindrical"
          else:
            raise ValueError("Unknown value for GEOMETRY flag ('{}') was found in the VTK file.".format(g))

          if(V.geometry != "unknown"):
            # We already have a proposed geometry, check that what is read from the file matches
            if thisgeometry != V.geometry:
              raise ValueError("geometry argument ('{}') is inconsistent with GEOMETRY flag from the VTK file ('{}')".format(V.geometry, thisgeometry))
          V.geometry=thisgeometry
        elif entry == "PERIODICITY":
          periodicity = np.fromfile(fid,dint,3).astype(bool)
          V.periodicity = tuple(periodicity)
        else:
          raise ValueError("Received unknown field: '{}'.".format(entry))

        s=fid.readline() # extra linefeed

      # finished reading the field entry
      # read next line
      s=fid.readline() #DIMENSIONS...

    if(V.geometry == "unknown"):
      raise RuntimeError(
        "Geometry couldn't be determined from data. "
        "Try to set the geometry keyword argument explicitely."
      )

    slist=s.split() # DIMENSIONS....
    #s=fid.readline()    # Extre line feed
    V.nx=int(slist[1])
    V.ny=int(slist[2])
    V.nz=int(slist[3])

    if V.geometry in ("cartesian", "cylindrical"):
      # CARTESIAN geometry
      # NOTE: cylindrical geometry is meant to be only used in 2D
      #       so the expected coordinates (R, z) are never curvilinear,
      #       which means we can treat them as cartesian
      s=fid.readline()    # X_COORDINATES NX float
      x=np.fromfile(fid,dt,V.nx)

      s=fid.readline()    # Extra line feed added by idefix


      s=fid.readline()    # X_COORDINATES NX float

      y=np.fromfile(fid,dt,V.ny)
      s=fid.readline()    # Extra line feed added by idefix

      s=fid.readline()    # X_COORDINATES NX float

      z=np.fromfile(fid,dt,V.nz)
      s=fid.readline()    # Extra line feed added by idefix



      s=fid.readline()    # POINT_DATA NXNYNZ

      slist=s.split()
      point_type=str(slist[0].decode('utf-8'))
      npoints=int(slist[1])
      s=fid.readline()    # EXTRA LINE FEED

      if(point_type == "CELL_DATA"):
          # The file contains face coordinates, so we extrapolate to get the cell center coordinates.
          if V.nx>1:
              V.nx=V.nx-1
              V.x=0.5*(x[1:]+x[:-1])
              # left side of the cell
              V.xl=x
          else:
              V.x=x
              V.xl=x
          if V.ny>1:
              V.ny=V.ny-1
              V.y=0.5*(y[1:]+y[:-1])
              # left side of the cell
              V.yl=y
          else:
              V.y=y
              V.yl=y
          if V.nz>1:
              V.nz=V.nz-1
              V.z=0.5*(z[1:]+z[:-1])
              V.zl=z
          else:
              V.z=z
              V.zl=z
      elif(point_type == "POINT_DATA"):
          V.x=x
          V.y=y
          V.z=z

      grid_size = V.nx*V.ny*V.nz
      if grid_size != npoints:
          raise RuntimeError("Grid size ({}) is incompatible with number of points in the dataset ({})".format(grid_size, npoints))
    else:
      # POLAR or SPHERICAL coordinates
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

      grid_size = V.nx*V.ny*V.nz
      if grid_size != npoints:
          raise RuntimeError("Grid size ({}) is incompatible with number of points in the dataset ({})".format(grid_size, npoints))

      x1d=points[::3]
      y1d=points[1::3]
      z1d=points[2::3]

      xcart=np.transpose(x1d.reshape(V.nz,V.ny,V.nx))
      ycart=np.transpose(y1d.reshape(V.nz,V.ny,V.nx))
      zcart=np.transpose(z1d.reshape(V.nz,V.ny,V.nx))

      # Reconstruct the polar coordinate system
      if V.geometry == "polar":

        r=np.sqrt(xcart[:,0,0]**2+ycart[:,0,0]**2)
        theta=np.unwrap(np.arctan2(ycart[0,:,0],xcart[0,:,0]))
        z=zcart[0,0,:]

        s=fid.readline()    # CELL_DATA (NX-1)(NY-1)(NZ-1)
        slist=s.split()
        data_type=str(slist[0].decode('utf-8'))
        if(data_type != "CELL_DATA"):
            print("ERROR: this routine expect CELL DATA as produced by idefix.")
            fid.close()
            return 0
        s=fid.readline()    # Line feed
        # Perform averaging on coordinate system to get cell centers
        # The file contains face coordinates, so we extrapolate to get the cell center coordinates.
        if V.nx>1:
            V.nx=V.nx-1
            V.x=0.5*(r[1:]+r[:-1])
            V.xl=r
        else:
            V.x=r
            V.xl=r
        if V.ny>1:
            V.ny=V.ny-1
            V.y=(0.5*(theta[1:]+theta[:-1])+np.pi)%(2.0*np.pi)-np.pi
            V.yl=(theta+np.pi)%(2.0*np.pi)-np.pi
        else:
            V.y=theta
            V.yl=theta
        if V.nz>1:
            V.nz=V.nz-1
            V.z=0.5*(z[1:]+z[:-1])
            V.zl=z
        else:
            V.z=z
            V.zl=z

      # Reconstruct the spherical coordinate system
      if V.geometry == "spherical":
        if(is2d):
          r=np.sqrt(xcart[:,0,0]**2+ycart[:,0,0]**2)
          phi=np.unwrap(np.arctan2(zcart[0,V.ny//2,:],xcart[0,V.ny//2,:]))
          theta=np.arccos(ycart[0,:,0]/np.sqrt(xcart[0,:,0]**2+ycart[0,:,0]**2))
        else:
            r=np.sqrt(xcart[:,0,0]**2+ycart[:,0,0]**2+zcart[:,0,0]**2)
            phi=np.unwrap(np.arctan2(ycart[V.nx//2,V.ny//2,:],xcart[V.nx//2,V.ny//2,:]))
            theta=np.arccos(zcart[0,:,0]/np.sqrt(xcart[0,:,0]**2+ycart[0,:,0]**2+zcart[0,:,0]**2))


        s=fid.readline()    # CELL_DATA (NX-1)(NY-1)(NZ-1)
        slist=s.split()
        data_type=str(slist[0].decode('utf-8'))
        if(data_type != "CELL_DATA"):
            print("ERROR: this routine expect CELL DATA as produced by IDEFIX.")
            fid.close()
            return 0
        s=fid.readline()    # Line feed
        # Perform averaging on coordinate system to get cell centers
        # The file contains face coordinates, so we extrapolate to get the cell center coordinates.
        if V.nx>1:
            V.nx=V.nx-1
            V.r=0.5*(r[1:]+r[:-1])
            V.rl=r
        else:
            V.r=r
            V.rl=r
        if V.ny>1:
            V.ny=V.ny-1
            V.theta=0.5*(theta[1:]+theta[:-1])
            V.thetal=theta
        else:
            V.theta=theta
            V.thetal=theta
        if V.nz>1:
            V.nz=V.nz-1
            V.phi=(0.5*(phi[1:]+phi[:-1]))
            V.phil=phi
        else:
            V.phi=phi
            V.phil=phi

    ## From that point, the coordinates system is known.
    while 1:
        s=fid.readline()        # SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
        #print repr(s)
        if len(s)<2:         # leave if end of file
            break
        slist=s.split()
        datatype=str(slist[0].decode('utf-8'))
        varname=str(slist[1].decode('utf-8'))
        if datatype == "SCALARS":
            fid.readline()  # LOOKUP TABLE
            V.data[varname] = np.transpose(np.fromfile(fid,dt,V.nx*V.ny*V.nz).reshape(V.nz,V.ny,V.nx))
        elif datatype == "VECTORS":
            Q=np.fromfile(fid,dt,3*V.nx*V.ny*V.nz)

            V.data[varname+'_X']=np.transpose(Q[::3].reshape(V.nz,V.ny,V.nx))
            V.data[varname+'_Y']=np.transpose(Q[1::3].reshape(V.nz,V.ny,V.nx))
            V.data[varname+'_Z']=np.transpose(Q[2::3].reshape(V.nz,V.ny,V.nx))

        else:
            raise RuntimeError("Unknown datatype '{}'".format(datatype))

        fid.readline()  #extra line feed
    fid.close()

    return V


# Former geometry-specific readers
def readVTKCart(filename):
  Warning("the use of readVTKCart is discouraged. Use the generic readVTK function with geometry='cartesian'")
  return readVTK(filename, geometry="cartesian")

# Read a vtk file
def readVTKPolar(filename):
  Warning("the use of readVTKPolar is discouraged. Use the generic readVTK function with geometry='polar'")
  return readVTK(filename, geometry="polar")


# Read a vtk file
def readVTKSpherical(filename):
  Warning("the use of readVTKSpherical is discouraged. Use the generic readVTK function with geometry='spherical'")
  return readVTK(filename, geometry="spherical")
