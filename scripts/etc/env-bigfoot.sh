#!/usr/bin/env bash
IDEFIX_CMAKE_OPTIONS+=( -DKokkos_ENABLE_HIP=ON -DKokkos_ARCH_VEGA90A=ON -DCMAKE_CXX_COMPILER=hipcc -DIdefix_MPI=ON )

function in_env() {
    local -a cmd=( "$@" )
    NIX_BUILD_SHELL=$(which bash) nix-shell ~/idefix.drv --run "$(declare -p cmd); "'"${cmd[@]}"'
}
