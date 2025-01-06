from pydefix import *
import numpy as np
import matplotlib.pyplot as plt

# The output function
# the only argument is dataBlockHost python object, wrapping a dataBlockHost Idefix object
def output(data,grid,n):
  # Note: if MPI is not enabled, GatherIdefixArray still works (but merely does a local copy)
  pressure = GatherIdefixArray(data.Vc[PRS,:,:,:],data,broadcast=False,keepBoundaries=True)
  # only process #0 performs the output
  if prank==0:
    plt.close()
    plt.figure()
    plt.pcolormesh(grid.x[IDIR],grid.x[JDIR],pressure[0,:,:],label='PRS',vmin=0.02,vmax=0.5,cmap='plasma')
    plt.title("t=%.2f"%data.t)
    plt.colorbar()
    plt.savefig("PRS.%.4d.png"%n)


def initflow(data):
  # Field amplitude
  B0 = 1/np.sqrt(4*np.pi)

  [z,y,x] = np.meshgrid(data.x[KDIR], data.x[JDIR], data.x[IDIR], indexing='ij')

  # Initialize the flow
  data.Vc[RHO,:,:,:] = 25/(36*np.pi)
  data.Vc[PRS,:,:,:] = 5/(12*np.pi)
  data.Vc[VX1,:,:,:] = -np.sin(2*np.pi*y)
  data.Vc[VX2,:,:,:] = np.sin(2*np.pi*x)

  data.Vs[BX1s,:,:-1,:-1] = -B0*np.sin(2*np.pi*y)
  data.Vs[BX2s,:,:-1,:-1] = B0*np.sin(4*np.pi*x)
