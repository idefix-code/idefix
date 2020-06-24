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

parser.add_argument("-arch", 
                    choices = ["HSW"
                             ,"BDW"
                             ,"SKX"
                             ,"Kepler30"
                             ,"Maxwell50"
                             ,"Pascal60"
                             ,"Pascal61"
                             ,"Volta70"
                             ,"Volta72"],
                     help="Target Kokkos architecture")

parser.add_argument("-openmp",
                    help="enable OpenMP parallelism",
                    action="store_true")

parser.add_argument("-mpi",
                    help="enable MPI parallelism",
                    action="store_true")

args=parser.parse_args()

idefixDir = os.getenv("IDEFIX_DIR")
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

if args.gpu:
    if(args.mpi):
        raise NotImplementedError('MPI+Cuda compiler is not implemented yet')
    makefileOptions['cxx'] = '${KOKKOS_PATH}/bin/nvcc_wrapper'
    makefileOptions['extraLine'] += '\nKOKKOS_CUDA_OPTIONS = "enable_lambda"'
    makefileOptions['kokkosDevices'] = '"Cuda"'
    makefileOptions['kokkosArch'] = '"Pascal60"'
    makefileOptions['cxxflags'] = "-O3"
else:
    if(args.mpi):
        makefileOptions['cxx'] = "mpicxx"
    else:
        makefileOptions['cxx'] = "g++"
    makefileOptions['kokkosArch'] = '"BDW"'
    makefileOptions['cxxflags'] += " -O3"
    if args.openmp:
         makefileOptions['kokkosDevices'] = '"OpenMP"'
    else:
        makefileOptions['kokkosDevices'] = '"Serial"'

if(args.mpi):
    makefileOptions['cxxflags'] += " -D WITH_MPI"

if args.mhd:
    makefileOptions['extraIncludeDir'] += " -I$(SRC)/hydro/MHDsolvers"
    makefileOptions['extraVpath'] += "$(SRC)/hydro/MHDsolvers"
    makefileOptions['extraHeaders'] += " solversMHD.hpp"
    makefileOptions['extraObj'] += " hlldMHD.o hllMHD.o roeMHD.o tvdlfMHD.o"
    makefileOptions['cxxflags'] += " -D MHD=YES"
else:
    makefileOptions['extraIncludeDir'] += " -I$(SRC)/hydro/HDsolvers"
    makefileOptions['extraVpath'] += "$(SRC)/hydro/HDsolvers"
    makefileOptions['extraHeaders'] += " solversHD.hpp"
    makefileOptions['extraObj'] += " hllcHD.o hllHD.o roeHD.o tvdlfHD.o"
    makefileOptions['cxxflags'] += " -D MHD=NO"

if args.cxx:
    makefileOptions['cxx'] = args.cxx



# Makefile substitution
with open(iMakefile, 'r') as file:
    makefile = file.read()

for key, val in makefileOptions.items():
    makefile = re.sub(r'@{0}@'.format(key), val, makefile)

with open(oMakefile,'w') as file:
    file.write(makefile)



