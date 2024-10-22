import pydefix as pdfx
import numpy as np
import matplotlib.pyplot as plt

# The output function
# the only argument is dataBlockHost python object, wrapping a dataBlockHost Idefix object
def output(data,n):
  plt.close()
  plt.figure()
  plt.pcolormesh(data.x[pdfx.IDIR],data.x[pdfx.JDIR],data.Vc[pdfx.PRS,0,:,:],label='PRS',vmin=0.02,vmax=0.5,cmap='plasma')
  plt.title("t=%.2f"%data.t)
  plt.colorbar()
  plt.savefig("PRS.%.4d.png"%n)


def initflow(data):
  # Field amplitude
  B0 = 1/np.sqrt(4*np.pi)

  [z,y,x] = np.meshgrid(data.x[pdfx.KDIR], data.x[pdfx.JDIR], data.x[pdfx.IDIR], indexing='ij')

  # Initialize the flow
  data.Vc[pdfx.RHO,:,:,:] = 25/(36*np.pi)
  data.Vc[pdfx.PRS,:,:,:] = 5/(12*np.pi)
  data.Vc[pdfx.VX1,:,:,:] = -np.sin(2*np.pi*y)
  data.Vc[pdfx.VX2,:,:,:] = np.sin(2*np.pi*x)

  data.Vs[pdfx.BX1s,:,:-1,:-1] = -B0*np.sin(2*np.pi*y)
  data.Vs[pdfx.BX2s,:,:-1,:-1] = B0*np.sin(4*np.pi*x)
