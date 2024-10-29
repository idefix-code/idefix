#!/usr/bin/env bash
function setup_env() {
    set +ue
    source /applis/site/nix.sh >/dev/null 2>&1
    set -ue
    
    mkdir -p "$HOME/.nix-shell"
    ln -fs "$(which bash)" "$HOME/.nix-shell/bash"
}

function in_env() {
    local -a cmd=( "$@" )
    NIX_BUILD_SHELL="$HOME/.nix-shell/bash" nix-shell ~/idefix.drv --run "$(declare -p cmd); "'"${cmd[@]}"'
}
