#!/usr/bin/env bash
source /applis/site/nix.sh

# IDEFIX_CMAKE_OPTIONS+=( -DKokkos_ENABLE_HIP=ON -DKokkos_ARCH_VEGA90A=ON -DCMAKE_CXX_COMPILER=hipcc )
mkdir -p "$HOME/.nix-shell"
ln -fs "$SHELL" "$HOME/.nix-shell/sh"

function in_env() {
    local -a cmd=( "$@" )
    NIX_BUILD_SHELL="$HOME/.nix-shell/sh" nix-shell ~/idefix.drv --run "$(declare -p cmd); "'"${cmd[@]}"'
}
