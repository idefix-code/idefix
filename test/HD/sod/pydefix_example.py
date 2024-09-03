
def output(data):
  print("data=")
  print(data.Vc[0,0,0,:])
  data.Vc[0,0,0,10] = 2.0

def initflow(data):
  # Initialize the flow
  print(data.Vc.shape)
  data.Vc[0,0,0,:] = 1+0*data.x[0][:]
  data.Vc[0,0,0,1] = 1.0
  data.Vc[1,0,0,:] = data.x[0][:]
  print(data.Vc[0,0,0,:])
  print(data.Vc[1,0,0,:])
