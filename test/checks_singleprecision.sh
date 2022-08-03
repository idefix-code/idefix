#!/bin/bash

rep_list="HD/sod-iso HD/sod MHD/sod MHD/sod-iso MHD/OrszagTang"
rep_MPI_list="MHD/OrszagTang3D"

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

# High order tests
for rep in $rep_list; do
    TMP_DIR="$(mktemp -d)"
    cp -R $TEST_DIR/$rep/* $TMP_DIR
    cd $TMP_DIR
    echo "***********************************************"
    echo "Configuring  $rep"
    echo "Using $TMP_DIR as working directory"
    echo "***********************************************"
    rm -f CMakeCache.txt

    cmake $IDEFIX_DIR $options -DIdefix_PRECISION=Single || { echo "!!!!$rep in single precision failed during configuration"; exit 1; }
    echo "***********************************************"
    echo "Making  $rep in single precision"
    echo "***********************************************"
    make clean; make -j 10 || { echo "!!!! $rep in single precision failed during compilation with"; exit 1; }

    ini_files=$(ls *.ini)
    for ini in $ini_files; do
        echo "***********************************************"
        echo "Running  $rep in single precision and $ini"
        echo "***********************************************"
        ./idefix -i $ini -nolog || { echo "!!!! $rep in single precision failed running with $ini"; exit 1; }

        cd python
        echo "***********************************************"
        echo "Testing  $rep in single precision and $ini"
        echo "***********************************************"
        python3 testidefix.py -noplot || { echo "!!!! $rep in single precision failed validation with $ini"; exit 1; }
        cd ..
    done
    make clean
    rm -f *.vtk *.dbl *.dmp

    echo "***********************************************"
    echo "Cleaning  $rep in $TMP_DIR"
    echo "***********************************************"
    rm -rf $TMP_DIR
done

#do it with MPI (only the default .ini files though)
# High order tests
for rep in $rep_MPI_list; do
    TMP_DIR="$(mktemp -d)"
    cp -R $TEST_DIR/$rep/* $TMP_DIR
    cd $TMP_DIR
    echo "***********************************************"
    echo "Configuring  $rep"
    echo "Using $TMP_DIR as working directory"
    echo "***********************************************"
    rm -f CMakeCache.txt

    cmake $IDEFIX_DIR $options -DIdefix_MPI=ON -DIdefix_PRECISION=Single || { echo "!!!!$rep in single precision failed during configuration"; exit 1; }
    echo "***********************************************"
    echo "Making  $rep in single precision"
    echo "***********************************************"
    make clean; make -j 10 || { echo "!!!! $rep in single precision and MPI failed during compilation with"; exit 1; }

    ini_files=$(ls *.ini)
    for ini in $ini_files; do
      echo "***********************************************"
      echo "Running  $rep in single precision, $ini and MPI"
      echo "***********************************************"
      mpirun -np 4 ./idefix -i $ini -nolog || { echo "!!!! $rep in single precision and MPI failed running "; exit 1; }

      cd python
      echo "***********************************************"
      echo "Testing  $rep in single precision and MPI"
      echo "***********************************************"
      python3 testidefix.py -noplot || { echo "!!!! $rep in single precision and MPI failed validation"; exit 1; }
      cd ..
    done
    make clean
    rm -f *.vtk *.dbl *.dmp

    echo "***********************************************"
    echo "Cleaning  $rep in $TMP_DIR"
    echo "***********************************************"
    rm -rf $TMP_DIR
done
