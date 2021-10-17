#!/usr/bin/env sh
# make two tarballs for public releases:
# - a full tarball including sources and tests.
# - a "lite" one (barebones source files, no tests)
# Project management orientied files are excluded in all cases, as well as
# non essential kokkos files with significant footprints.


set -euxo pipefail

VERSION=$(git describe --tags --always)
BLACK_LIST=".git|.pre-commit-config.yaml|CPPLINT.cfg|publish.sh|src/kokkos/core|src/kokkos/example/"

function make_tarball () {
    local black_list=${1:?Must provide a black_list}
    local output=${2:?Must provide an output}.tar

    # ensure kokkos version is clean
    git submodule update --init

    if [[ $(git diff --name-only) ]]; then
        echo "local copy is dirty. Commit or stash changes before publishing."
        exit 1
    fi


    git ls-files | egrep -v $black_list | xargs tar -cvf $output

    echo "#define GITVERSION $VERSION" > src/gitversion.hpp
    tar -rvf $output src/gitversion.hpp

    gzip $output
}

make_tarball "($BLACK_LIST)" /tmp/idefix_$VERSION
make_tarball "($BLACK_LIST|test)" /tmp/idefix_"$VERSION"_lite
