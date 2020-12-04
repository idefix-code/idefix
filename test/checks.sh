#!/bin/sh   

rep_HD_list="sod-iso sod MachReflection"
rep_MHD_list="sod-iso sod AmbipolarCshock HallWhistler OrszagTang"

# refer to the parent dir of this file, wherever this is called from
# a python equivalent is e.g.
#
# import pathlib
# TEST_DIR = pathlib.Path(__file__).parent
TEST_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# todo: add a warning here if the env var is already defined _and_ different
export IDEFIX_DIR=$TEST_DIR/..
echo $IDEFIX_DIR

set -e
options=""

# HD tests
for rep in $rep_HD_list; do
    cd $TEST_DIR/HD/$rep
    echo "***********************************************"
    echo "Configuring  $rep"
    echo "***********************************************"
    python3 $IDEFIX_DIR/configure.py $options
    echo "***********************************************"
    echo "Making  $rep"
    echo "***********************************************"
    make clean; make -j 4

    ini_files=$(ls *.ini)
    for ini in $ini_files; do
        echo "***********************************************"
        echo "Running  $rep with $ini"
        echo "***********************************************"
        ./idefix -i $ini

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
for rep in $rep_MHD_list; do
    cd $TEST_DIR/MHD/$rep
    echo "***********************************************"
    echo "Configuring  $rep"
    echo "***********************************************"
    python3 $IDEFIX_DIR/configure.py -mhd $options
    echo "***********************************************"
    echo "Making  $rep"
    echo "***********************************************"
    make clean; make -j 4

    ini_files=$(ls *.ini)
    for ini in $ini_files; do
        echo "***********************************************"
        echo "Running  $rep with $ini"
        echo "***********************************************"
        ./idefix -i $ini

        cd python
        echo "***********************************************"
        echo "Testing  $rep with $ini"
        echo "***********************************************"
        python3 testidefix.py -noplot
        cd ..
    done
    
    cd $TEST_DIR
done



