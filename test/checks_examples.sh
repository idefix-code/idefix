#!/bin/bash

# check that the provided examples (which do not have a reference case) compile and run for a few loops

rep_example_list="HD/KHI HD/RWI-cavity HD/VSI MHD/RotorCartesian MHD/RotorPolar MHD/advectionOBlique MHD/AmbipolarShearingBox MHD/AmbipolarWind MHD/FieldLoop MHD/HallDisk MHD/disk MHD/diskSpherical Dust/StreamingInstability Dust/FargoPlanet"

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


TMP_DIR="$(mktemp -d)"
function finish ()
{
  echo $1
  echo "Cleaning directory $TMP_DIR"
  cd $TEST_DIR
  rm -rf $TMP_DIR
  exit 1
}

for rep in $rep_example_list; do

    cp -R $TEST_DIR/$rep/* $TMP_DIR
    cd $TMP_DIR
    echo "***********************************************"
    echo "Configuring  $rep"
    echo "Using $TMP_DIR as working directory"
    echo "***********************************************"

    cmake $IDEFIX_DIR $options || finish "!!!! Example $rep failed during configuration"
    echo "***********************************************"
    echo "Making  $rep"
    echo "***********************************************"
    make -j 8 || finish  "!!!! Example $rep failed during compilation"


    echo "***********************************************"
    echo "Running  $rep"
    echo "***********************************************"
    ./idefix -maxcycles 10 -nowrite -Werror || finish "!!!! Example $rep failed running"

    echo "***********************************************"
    echo "Cleaning  $rep in $TMP_DIR"
    echo "***********************************************"
   rm -rf *.vtk *.dbl *.dmp *.ini python CMakeLists.txt CMakeCache.txt

done

echo "Test was successfull"
cd $TEST_DIR
rm -rf $TMP_DIR
