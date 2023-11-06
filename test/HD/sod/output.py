
def myfunc(n, data):
  print("Python n=%d"%n)
  print(dir(data))
  print("Dimensions="+str(data.Vc.shape))
  print("x1=")
  print(data.x[0][:])
  print("data=")
  print(data.Vc[0,0,0,:])
