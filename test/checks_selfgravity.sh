#!/bin/bash

rep_SG_list="JeansInstability UniformCollapse RandomSphere RandomSphereCartesian"

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

# HD tests
for rep in $rep_SG_list; do

    cp -R $TEST_DIR/SelfGravity/$rep/* $TMP_DIR
    cd $TMP_DIR
    echo "***********************************************"
    echo "Configuring  $rep"
    echo "Using $TMP_DIR as working directory"
    echo "***********************************************"

    def_files=$(ls definitions*.hpp)
    for def in $def_files; do

        cmake $IDEFIX_DIR $options -DIdefix_DEFS=$def|| finish "!!!! Self Gravity $rep failed during configuration"
        echo "***********************************************"
        echo "Making  $rep with $def"
        echo "***********************************************"
        make -j 8 || finish "!!!! Self Gravity $rep failed during compilation with $def"

        ini_files=$(ls *.ini)
        for ini in $ini_files; do
            echo "***********************************************"
            echo "Running  $rep with $ini"
            echo "***********************************************"
            ./idefix -i $ini -nolog || finish "!!!! Self Gravity $rep failed running with $def and $ini"

            cd python
            echo "***********************************************"
            echo "Testing  $rep with $ini and $def"
            echo "***********************************************"
            python3 testidefix.py -noplot || finish "!!!! Self Gravity $rep failed validation with $def and $ini"
            cd ..
        done
        rm -f *.vtk *.dbl
    done
    echo "***********************************************"
    echo "Cleaning  $rep in $TMP_DIR"
    echo "***********************************************"
   rm -rf *.vtk *.dbl *.dmp *.ini python CMakeLists.txt
done
echo "Test was successfull"
cd $TEST_DIR
rm -rf $TMP_DIR
