#!/bin/sh   

rep_list="sod-iso sod"

cwd=$(pwd)
export IDEFIX_DIR=$cwd/../..
echo $IDEFIX_DIR

set -e
for rep in $rep_list; do
    cd $rep
    echo "***********************************************"
    echo "Configuring  $rep"
    echo "***********************************************"
    python $IDEFIX_DIR/configure.py
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
        python testidefix.py -noplot
        cd ..
    done
    
    cd $cwd
done


