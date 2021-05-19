#!/bin/bash

rep_HD_2D_mpi_list="MachReflection ViscousFlowPastCylinder"
rep_MHD_2D_mpi_list="OrszagTang"
rep_MHD_3D_mpi_list="OrszagTang3D"

# refer to the parent dir of this file, wherever this is called from
# a python equivalent is e.g.
#
# import pathlib
# TEST_DIR = pathlib.Path(__file__).parent
TEST_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


function resolve_path {
    # resolve relative paths
    # work around the fact that `realpath` is not bundled with every UNIX distro
    echo "`cd "$1";pwd`"
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

# HD MPI tests
for rep in $rep_HD_2D_mpi_list; do
    cd $TEST_DIR/HD/$rep
    echo "***********************************************"
    echo "Configuring  $rep"
    echo "***********************************************"
    python3 $IDEFIX_DIR/configure.py -mpi $options
    echo "***********************************************"
    echo "Making  $rep"
    echo "***********************************************"
    make clean; make -j 4

    ini_files=$(ls *.ini)
    for ini in $ini_files; do
        echo "***********************************************"
        echo "Running  $rep with $ini"
        echo "***********************************************"
        mpirun -np 4 ./idefix -i $ini -dec 2 2

        cd python
        echo "***********************************************"
        echo "Testing  $rep with $ini"
        echo "***********************************************"
        python3 testidefix.py -noplot
        cd ..
    done
    make clean
    rm -f *.vtk *.dbl
done

# MHD tests
for rep in $rep_MHD_2D_mpi_list; do
    cd $TEST_DIR/MHD/$rep
    echo "***********************************************"
    echo "Configuring  $rep"
    echo "***********************************************"
    python3 $IDEFIX_DIR/configure.py -mhd -mpi $options
    echo "***********************************************"
    echo "Making  $rep"
    echo "***********************************************"
    make clean; make -j 4

    ini_files=$(ls *.ini)
    for ini in $ini_files; do
        echo "***********************************************"
        echo "Running  $rep with $ini"
        echo "***********************************************"
        mpirun -np 4 ./idefix -i $ini -dec 2 2

        cd python
        echo "***********************************************"
        echo "Testing  $rep with $ini"
        echo "***********************************************"
        python3 testidefix.py -noplot
        cd ..
    done

    cd $TEST_DIR
done

# MHD tests
for rep in $rep_MHD_3D_mpi_list; do
    cd $TEST_DIR/MHD/$rep
    echo "***********************************************"
    echo "Configuring  $rep"
    echo "***********************************************"
    python3 $IDEFIX_DIR/configure.py -mhd -mpi $options
    echo "***********************************************"
    echo "Making  $rep"
    echo "***********************************************"
    make clean; make -j 4

    ini_files=$(ls *.ini)
    for ini in $ini_files; do
        echo "***********************************************"
        echo "Running  $rep with $ini"
        echo "***********************************************"
        mpirun -np 8 ./idefix -i $ini -dec 2 2 2

        cd python
        echo "***********************************************"
        echo "Testing  $rep with $ini"
        echo "***********************************************"
        python3 testidefix.py -noplot
        cd ..
    done

    cd $TEST_DIR
done
