import pydefix as pdfx
import numpy as np
import matplotlib.pyplot as plt

def output(data):

  plt.figure()
  plt.plot(data.x[pdfx.IDIR],data.Vc[pdfx.VX1,0,0,:],label='VX1')
  plt.plot(data.x[pdfx.IDIR],data.Vc[pdfx.RHO,0,0,:],label='RHO')
  plt.legend()
  plt.show()
  #data.Vc[0,0,0,10] = 2.0

def initflow(data):

  # Initialize the flow
  data.Vc[pdfx.RHO,0,0,:] = 1.0
  data.Vc[pdfx.PRS,0,0,:] = 2.0
  data.Vc[pdfx.VX1,0,0,:] = np.sin(2.0*np.pi*data.x[0][:])
