#!/bin/bash

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

# Validate LookupTable
rep="utils/lookupTable"

TMP_DIR="$(mktemp -d)"
cp -R $TEST_DIR/$rep/* $TMP_DIR
cd $TMP_DIR
echo "***********************************************"
echo "Making numpy test files"
echo "***********************************************"
python3 makeNpy.py
echo "***********************************************"
echo "Configuring  $rep"
echo "Using $TMP_DIR/$rep as working directory"
echo "***********************************************"
rm -f CMakeCache.txt
cmake $IDEFIX_DIR -DIdefix_MPI=ON  $options
echo "***********************************************"
echo "Making  $rep"
echo "***********************************************"
make clean; make -j 10
mpirun -np 4 ./idefix
cd $TEST_DIR
rm -rf $TMP_DIR


# Validate DumpImage
rep="utils/dumpImage"

TMP_DIR="$(mktemp -d)"
cp -R $TEST_DIR/$rep/* $TMP_DIR
cd $TMP_DIR
echo "***********************************************"
echo "Configuring  $rep"
echo "Using $TMP_DIR/$rep as working directory"
echo "***********************************************"
rm -f CMakeCache.txt
cmake $IDEFIX_DIR  $options
echo "***********************************************"
echo "Making  $rep"
echo "***********************************************"
make clean; make -j 10
./idefix
cd $TEST_DIR
rm -rf $TMP_DIR
