#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 15:33:28 2020

@author: lesurg
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 17:09:32 2020

@author: lesurg
"""
import os
import sys
TESTDIR_PATH = os.path.join(os.getenv("IDEFIX_DIR"), "test")
sys.path.append(TESTDIR_PATH)
from idefix_testing.framework import readVTKCart
import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=False,
                    help="disable plotting",
                    action="store_true")


args=parser.parse_args()

rep='../'
nend=1000

# values for the field strength to compute theoretical eigenfrequency values
lH=1.0
B=1.0
Va=1.0
k=2.0*np.pi


n=0
dt=0.01

V=readVTKCart(rep+'/data.'+'%0*d'%(4,n)+'.vtk')
nx=V.x.size*V.y.size*V.z.size
nx=V.x.size

Bx=np.zeros((nend,nx))
By=np.zeros((nend,nx))
Bz=np.zeros((nend,nx))

Jx=np.zeros((nend,nx))
Jy=np.zeros((nend,nx))
Jz=np.zeros((nend,nx))

Vx=np.zeros((nend,nx))
Vy=np.zeros((nend,nx))
Vz=np.zeros((nend,nx))

t=dt*np.arange(0,nend)
for n in range(nend):
    V=readVTKCart(rep+'/data.'+'%0*d'%(4,n)+'.vtk')


    Bx[n,:]=V.data['BX1'][:,0,0]
    By[n,:]=V.data['BX2'][:,0,0]
    Bz[n,:]=V.data['BX3'][:,0,0]





# Theoretical speedes
f_w=Va/(2*np.pi)*k*(np.sqrt(1+(k*lH/2)**2)+k*lH/2)
f_i=Va/(2*np.pi)*k*(np.sqrt(1+(k*lH/2)**2)-k*lH/2)
f_A=k*Va/(2*np.pi)



# Compute temporal spectrum and get frequency
by=By[:,0]
sp=np.abs(np.fft.rfft(by))
f=np.arange(sp.size)/(t[-1]-t[0])

r=10**np.arange(np.floor(np.log10(np.amin(sp))),np.ceil(np.log10(np.amax(sp))),0.1)

imax=argrelextrema(sp,np.greater,order=10)[0][0]
f_wnum=f[imax]
error=np.fabs(f_w-f_wnum)/f_w
print("Theoretical whistler frequency=%g, numerical=%g, error=%g"%(f_w,f_wnum,error))


if(not args.noplot):

    plt.close('all')
    plt.figure()
    plt.contourf(t,V.x,By.T,64)
    plt.title('By')
    plt.colorbar()

    plt.figure()
    plt.contourf(t,V.x,Bz.T,64)
    plt.title('Bz')
    plt.colorbar()

    plt.figure()
    plt.loglog(f,sp,label='signal')
    plt.loglog(f_w+0*r,r,'--',label='whistler')
    plt.loglog(f_i+0*r,r,'-.',label='ion cyclotron')
    plt.loglog(f_A+0*r,r,':',label='Alfven')
    plt.xlabel('frequency')
    plt.ylabel('By amplitude')
    plt.legend()

    plt.ioff()
    plt.show()



if(error<0.05):
    print("SUCCESS")
    sys.exit(0)
else:
    print("Failed")
    sys.exit(1)
