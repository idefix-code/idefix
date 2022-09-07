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
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=False,
                    help="disable plotting",
                    action="store_true")


args, unknown=parser.parse_known_args()


# values for the field strength to compute theoretical eigenfrequency values
lH=1.0
B=1.0
Va=1.0
k=2.0*np.pi


# load the dat file produced by the setup
raw=np.loadtxt('../timevol.dat',skiprows=1)
t=raw[:,0]
by=raw[:,1]


# Theoretical speedes
f_w=Va/(2*np.pi)*k*(np.sqrt(1+(k*lH/2)**2)+k*lH/2)
f_i=Va/(2*np.pi)*k*(np.sqrt(1+(k*lH/2)**2)-k*lH/2)
f_A=k*Va/(2*np.pi)



# Compute temporal spectrum and get frequency
#by=By[:,0]
sp=np.abs(np.fft.rfft(by))
f=np.arange(sp.size)/(t[-1]-t[0])

r=10**np.arange(np.floor(np.log10(np.amin(sp))),np.ceil(np.log10(np.amax(sp))),0.1)




if(not args.noplot):

    plt.close('all')


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


imax=argrelextrema(sp,np.greater,order=2)[0][0]
f_wnum=f[imax]
error=np.fabs(f_w-f_wnum)/f_w
print("Theoretical whistler frequency=%g, numerical=%g, error=%g"%(f_w,f_wnum,error))

if(error<0.06):
    print("SUCCESS")
    sys.exit(0)
else:
    print("Failed")
    sys.exit(1)
