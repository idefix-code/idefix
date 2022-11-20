import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-noplot",
                    default=False,
                    help="disable plotting",
                    action="store_true")


args, unknown=parser.parse_known_args()


# values for the field strength to compute theoretical decay rate
eta=0.05
k=2.0*np.pi
B0=1.0

# decay rate of the energy, assuming k^2eta^2 < 2B0^2
gamma=-k**2*eta/2

raw=np.loadtxt('../timevol.dat',skiprows=1)
t=raw[:,0]
e=raw[:,1]

eth=e[0]*np.exp(2*gamma*t)
error=np.mean(e/eth-1.0)

if(not args.noplot):

    plt.close('all')
    plt.figure()
    plt.semilogy(t,e)
    plt.semilogy(t,eth)

    plt.ioff()
    plt.show()

if(error<0.03):
    print("SUCCESS")
    sys.exit(0)
else:
    print("Failed")
    sys.exit(1)
