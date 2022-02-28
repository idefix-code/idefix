#!/bin/bash
rep_MHD_list="AmbipolarCshock3D FargoMHDSpherical ShearingBox OrszagTang OrszagTang3D"

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


# MHD+Vector potential tests
for rep in $rep_MHD_list; do
    cp -R $TEST_DIR/MHD/$rep $TMP_DIR
    cd $TMP_DIR/$rep
    echo "***********************************************"
    echo "Configuring  $rep"
    echo "Using $TMP_DIR/$rep as working directory"
    echo "***********************************************"
    rm -f CMakeCache.txt
    def_files=$(ls definitions*.hpp)
    for def in $def_files; do
        cmake $IDEFIX_DIR $options -DIdefix_EVOLVE_VECTOR_POTENTIAL=ON -DIdefix_DEFS=$def|| { echo "!!!! MHD $rep failed during configuration"; exit 1; }
        echo "***********************************************"
        echo "Making  $rep with $def"
        echo "***********************************************"
        make clean; make -j 10 || { echo "!!!! MHD $rep failed during compilation with $def"; exit 1; }

        ini_files=$(ls *.ini)
        for ini in $ini_files; do
            echo "***********************************************"
            echo "Running  $rep with $ini"
            echo "***********************************************"
            ./idefix -i $ini -nolog || { echo "!!!! MHD $rep failed running with $def and $ini"; exit 1; }

            cd python
            echo "***********************************************"
            echo "Testing  $rep with $ini and $def"
            echo "***********************************************"
            python3 testidefix.py -noplot || { echo "!!!! MHD $rep failed validation with $def and $ini"; exit 1; }
            cd ..
        done
    done
    cd $TEST_DIR
done

echo "Cleaning temporary directory $TMP_DIR"
rm -rf $TMP_DIR
