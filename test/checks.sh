#!/bin/sh   

rep_HD_list="sod-iso sod"
rep_MHD_list="sod-iso sod AmbipolarCshock HallWhistler OrszagTang"

cwd=$(pwd)
export IDEFIX_DIR=$cwd/..
echo $IDEFIX_DIR

set -e

# HD tests
for rep in $rep_HD_list; do
    cd HD/$rep
    echo "***********************************************"
    echo "Configuring  $rep"
    echo "***********************************************"
    python3 $IDEFIX_DIR/configure.py
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

    cd $cwd
done

# MHD tests
for rep in $rep_MHD_list; do
    cd MHD/$rep
    echo "***********************************************"
    echo "Configuring  $rep"
    echo "***********************************************"
    python3 $IDEFIX_DIR/configure.py -mhd
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
    
    cd $cwd
done



