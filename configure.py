import argparse
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument("-mhd", 
                    default=False,
                    help="enable MHD",
                    action="store_true")
parser.add_argument("-gpu",
                    default=False,
                    help="Enable KOKKOS+CUDA",
                    action="store_true")

parser.add_argument("-cxx",
                    help="Override default compiler")

cpuarch = ["HSW"
          ,"BDW"
          ,"SKX"
          ,"EPYC"]

gpuarch = ["Kepler30"
          ,"Maxwell50"
          ,"Pascal60"
          ,"Pascal61"
          ,"Volta70"
          ,"Volta72"
          ,"Turing75"]

parser.add_argument("-arch",
                    nargs='+',
                    choices = cpuarch + gpuarch,
                     help="Target Kokkos architecture")

parser.add_argument("-openmp",
                    help="enable OpenMP parallelism",
                    action="store_true")

parser.add_argument("-mpi",
                    help="enable MPI parallelism",
                    action="store_true")

args=parser.parse_args()

idefixDir = os.getenv("IDEFIX_DIR")
if idefixDir is None:
    print("Environment variable IDEFIX_DIR is not set")
    print("Please export IDEFIX_DIR=/path/to/idefix")
    exit()

iMakefile = idefixDir+"/Makefile.in"
oMakefile = "Makefile"

# List of objects
makefileOptions = {}
makefileOptions['extraIncludeDir']=""
makefileOptions['extraVpath']=""
makefileOptions['extraHeaders']=""
makefileOptions['extraObj']=""
makefileOptions['extraLine']=""
makefileOptions['cxxflags']=""
makefileOptions['ldflags']=""


# extract cpu & gpu architectures from args.arch
cpu=""
gpu=""

if args.arch is not None:
    for archItem in args.arch:
        if archItem in cpuarch:
            cpu = archItem
        if archItem in gpuarch:
            gpu = archItem

if cpu == "":
    cpu="BDW"

if gpu == "":
    gpu="Pascal60"

if args.gpu:
    makefileOptions['cxx'] = '${KOKKOS_PATH}/bin/nvcc_wrapper'
    makefileOptions['extraLine'] += '\nKOKKOS_CUDA_OPTIONS = "enable_lambda"'
    makefileOptions['kokkosDevices'] = '"Cuda"'
    makefileOptions['kokkosArch'] = cpu+","+gpu
    makefileOptions['cxxflags'] = "-O3 "

    # This assumes openmpi. TODO: do a more general routine for all compilers
    if(args.mpi):
        stream=os.popen('mpicxx --showme:compile')
        makefileOptions['cxxflags'] += stream.read().strip()
        stream=os.popen('mpicxx --showme:link')
        makefileOptions['ldflags'] += stream.read().strip()
else:
    if(args.mpi):
        makefileOptions['cxx'] = "mpicxx"
    else:
        makefileOptions['cxx'] = "g++"
    makefileOptions['kokkosArch'] = cpu
    makefileOptions['cxxflags'] = "-O3"
    if args.openmp:
         makefileOptions['kokkosDevices'] = '"OpenMP"'
    else:
        makefileOptions['kokkosDevices'] = '"Serial"'

if(args.mpi):
    makefileOptions['cxxflags'] += " -DWITH_MPI"

if args.mhd:
    makefileOptions['extraIncludeDir'] += " -I$(SRC)/hydro/MHDsolvers"
    makefileOptions['extraVpath'] += "$(SRC)/hydro/MHDsolvers"
    makefileOptions['extraHeaders'] += " solversMHD.hpp"
    makefileOptions['extraObj'] += " hlldMHD.o hllMHD.o roeMHD.o tvdlfMHD.o"
    makefileOptions['cxxflags'] += " -DMHD=YES"
else:
    makefileOptions['extraIncludeDir'] += " -I$(SRC)/hydro/HDsolvers"
    makefileOptions['extraVpath'] += "$(SRC)/hydro/HDsolvers"
    makefileOptions['extraHeaders'] += " solversHD.hpp"
    makefileOptions['extraObj'] += " hllcHD.o hllHD.o roeHD.o tvdlfHD.o"
    makefileOptions['cxxflags'] += " -DMHD=NO"

if args.cxx:
    makefileOptions['cxx'] = args.cxx



# Makefile substitution
with open(iMakefile, 'r') as file:
    makefile = file.read()

for key, val in makefileOptions.items():
    makefile = re.sub(r'@{0}@'.format(key), val, makefile)

with open(oMakefile,'w') as file:
    file.write(makefile)



