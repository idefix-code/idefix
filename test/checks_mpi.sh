#!/bin/bash

rep_2D_mpi_list="HD/MachReflection HD/ViscousFlowPastCylinder HD/FargoPlanet MHD/OrszagTang"
rep_3D_mpi_list="MHD/AxisFluxTube MHD/LinearWaveTest MHD/FargoMHDSpherical MHD/OrszagTang3D"
rep_3D_noX3_list="MHD/AmbipolarCshock3D SelfGravity/RandomSphere"

# refer to the parent dir of this file, wherever this is called from
# a python equivalent is e.g.
#
# import pathlib
# TEST_DIR = pathlib.Path(__file__).parent
TEST_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
TMP_DIR="$(mktemp -d)"

function resolve_path {
    # resolve relative paths
    # work around the fact that `realpath` is not bundled with every UNIX distro
    echo "`cd "$1";pwd`"
}

function finish ()
{
  echo $1
  echo "Cleaning directory $TMP_DIR"
  cd $TEST_DIR
  rm -rf $TMP_DIR
  exit 1
}

target_dir=$(resolve_path $TEST_DIR/..)

if [ -z ${var+IDEFIX_DIR} ] & [ -d "$IDEFIX_DIR" ] ; then
    global_dir=$(resolve_path $IDEFIX_DIR)
    if [ $target_dir != $global_dir ] ; then
        echo \
        "Warning: IDEFIX_DIR is set globally to $global_dir" \
        ", but is redefined to $target_dir"
    fi
fi

export IDEFIX_DIR=$target_dir
echo $IDEFIX_DIR

set -e
options=$@

# 2D MPI tests
for rep in $rep_2D_mpi_list; do
    cp -R $TEST_DIR/$rep/* $TMP_DIR
    cd $TMP_DIR
    echo "***********************************************"
    echo "Configuring  $rep"
    echo "Using $TMP_DIR/$rep as working directory"
    echo "***********************************************"
    cmake $IDEFIX_DIR -DIdefix_MPI=ON  $options || finish "!!!! MPI $rep failed during configuration"
    echo "***********************************************"
    echo "Making  $rep"
    echo "***********************************************"
    make -j 8 || finish "!!!! MPI $rep failed during compilation"

    ini_files=$(ls *.ini)
    for ini in $ini_files; do
        echo "***********************************************"
        echo "Running  $rep with $ini"
        echo "***********************************************"
        mpirun -np 4 ./idefix -i $ini -dec 2 2 -nolog || finish "!!!! MPI $rep failed during runtime with $ini"

        cd python
        echo "***********************************************"
        echo "Testing  $rep with $ini"
        echo "***********************************************"
        python3 testidefix.py -noplot -i ../$ini || finish "!!!! MPI $rep failed during validation"
        cd ..
    done
   rm -rf *.vtk *.dbl *.dmp *.ini python CMakeLists.txt
done

# MHD tests
for rep in $rep_3D_noX3_list; do
    cp -R $TEST_DIR/$rep/* $TMP_DIR
    cd $TMP_DIR
    echo "***********************************************"
    echo "Configuring  $rep"
    echo "Using $TMP_DIR/$rep as working directory"
    echo "***********************************************"
    cmake $IDEFIX_DIR -DIdefix_MPI=ON $options || finish "!!!! MPI $rep failed during configuration"
    echo "***********************************************"
    echo "Making  $rep"
    echo "***********************************************"
    make -j 8 || finish "!!!! MPI $rep failed during compilation"

    ini_files=$(ls *.ini)
    for ini in $ini_files; do
        echo "***********************************************"
        echo "Running  $rep with $ini"
        echo "***********************************************"
        mpirun -np 4 ./idefix -i $ini -dec 2 2 1 -nolog || finish "!!!! MPI $rep failed during runtime with $ini"

        cd python
        echo "***********************************************"
        echo "Testing  $rep with $ini"
        echo "***********************************************"
        python3 testidefix.py -noplot || finish "!!!! MPI $rep failed during validation"
        cd ..
    done
   rm -rf *.vtk *.dbl *.dmp *.ini python CMakeLists.txt
    cd $TEST_DIR
done

for rep in $rep_3D_mpi_list; do
    cp -R $TEST_DIR/$rep/* $TMP_DIR
    cd $TMP_DIR
    echo "***********************************************"
    echo "Configuring  $rep"
    echo "Using $TMP_DIR/$rep as working directory"
    echo "***********************************************"
    cmake $IDEFIX_DIR -DIdefix_MPI=ON $options || finish "!!!! MPI $rep failed during configuration"
    echo "***********************************************"
    echo "Making  $rep"
    echo "***********************************************"
    make -j 8 || finish "!!!! MPI $rep failed during compilation"

    ini_files=$(ls *.ini)
    for ini in $ini_files; do
        echo "***********************************************"
        echo "Running  $rep with $ini"
        echo "***********************************************"
        mpirun -np 8 ./idefix -i $ini -dec 2 2 2 -nolog || finish "!!!! MPI $rep failed during runtime with $ini"

        cd python
        echo "***********************************************"
        echo "Testing  $rep with $ini"
        echo "***********************************************"
        python3 testidefix.py -noplot -i ../$ini  || finish "!!!! MPI $rep failed during validation"
        cd ..
    done
   rm -rf *.vtk *.dbl *.dmp *.ini python CMakeLists.txt
    cd $TEST_DIR
done

# Blast wave test
rep="HD/SedovBlastWave"

## Cartesian blast
cp -R $TEST_DIR/$rep/* $TMP_DIR
cd $TMP_DIR
echo "***********************************************"
echo "Configuring  $rep"
echo "Using $TMP_DIR/$rep as working directory"
echo "***********************************************"
cmake $IDEFIX_DIR -DIdefix_MPI=ON -DIdefix_MHD=OFF $options || finish "!!!! MPI $rep failed during configuration"
echo "***********************************************"
echo "Making  $rep"
echo "***********************************************"
make -j 8 || finish "!!!! MPI $rep failed during compilation"

ini="idefix.ini"
echo "***********************************************"
echo "Running  $rep with $ini"
echo "***********************************************"
mpirun -np 8 ./idefix -i $ini -dec 2 2 2 -nolog || finish "!!!! MPI $rep failed during runtime with $ini"

cd python
echo "***********************************************"
echo "Testing  $rep with $ini"
echo "***********************************************"
python3 testidefix.py -noplot -i ../$ini  || finish "!!!! MPI $rep failed during validation"
cd ..

## Spherical blast
echo "***********************************************"
echo "Configuring  $rep in spherical geometry"
echo "Using $TMP_DIR/$rep as working directory"
echo "***********************************************"

cmake $IDEFIX_DIR -DIdefix_MPI=ON -DIdefix_MHD=OFF -DIdefix_DEFS=definitions-spherical.hpp $options || finish "!!!! MPI $rep failed during configuration"
echo "***********************************************"
echo "Making  $rep"
echo "***********************************************"
make -j 8 || finish "!!!! MPI $rep failed during compilation"

ini="idefix-spherical.ini"
echo "***********************************************"
echo "Running  $rep with $ini"
echo "***********************************************"
mpirun -np 8 ./idefix -i $ini -dec 2 2 2 -nolog || finish "!!!! MPI $rep failed during runtime with $ini"

cd python
echo "***********************************************"
echo "Testing  $rep with $ini"
echo "***********************************************"
python3 testidefix.py -noplot -i ../$ini  || finish "!!!! MPI $rep failed during validation"

## done
echo "Test was successfull"
cd $TEST_DIR
rm -rf $TMP_DIR
