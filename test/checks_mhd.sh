#!/bin/bash
rep_MHD_list="sod-iso sod AxisFluxTube AmbipolarCshock HallWhistler ResistiveAlfvenWave FargoMHDSpherical ShearingBox OrszagTang OrszagTang3D"

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


# MHD tests
for rep in $rep_MHD_list; do
    cd $TEST_DIR/MHD/$rep
    echo "***********************************************"
    echo "Configuring  $rep"
    echo "***********************************************"
    rm -f CMakeCache.txt
    def_files=$(ls definitions*.hpp)
    for def in $def_files; do
        cmake $IDEFIX_DIR $options -DIdefix_DEFS=$def -DKokkos_ENABLE_DEBUG_BOUNDS_CHECK=ON|| { echo "!!!! MHD $rep failed during configuration"; exit 1; }
        echo "***********************************************"
        echo "Making  $rep with $def"
        echo "***********************************************"
        make clean; make -j 4 || { echo "!!!! MHD $rep failed during compilation with $def"; exit 1; }

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

# Test restart functions with OT3D which have generated a dump during the first pass
rep=OrszagTang3D
cd $TEST_DIR/MHD/$rep
# remove generated vtk from previous run
mv data.0001.vtk data.0001.old.vtk
echo "***********************************************"
echo "Running  $rep with restart dump # 1"
echo "***********************************************"
./idefix -restart 1 || { echo "!!!! MHD $rep failed running restart dump validation"; exit 1; }
cd python
echo "***********************************************"
echo "Testing  $rep with restart dump # 1"
echo "***********************************************"
python3 testidefix.py -noplot || { echo "!!!! MHD $rep failed checking restart dump validation"; exit 1; }
cd ..
echo "***********************************************"
echo "Checking bitwise compatibility of output from restarts"
echo "***********************************************"
diff data.0001.vtk data.0001.old.vtk || { echo "!!!! MHD $rep failed: restart dumps do not produce exactly the same results"; exit 1; }
echo "Success"
