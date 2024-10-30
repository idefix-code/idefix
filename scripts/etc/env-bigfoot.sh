#!/usr/bin/env bash
function setup_env() {
    set +ue
    source /applis/site/nix.sh >/dev/null 2>&1
    set -ue
    
    mkdir -p "$HOME/.nix-shell"
    ln -fs "$(which bash)" "$HOME/.nix-shell/bash"
}

function set_gpu_options() {
    local model="$1"; shift
    case "$model" in
        V100)
            IDEFIX_CMAKE_OPTIONS+=( -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_VOLTA70=ON );;
        A100)
            IDEFIX_CMAKE_OPTIONS+=( -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_AMPERE80=ON );;
        H100)
            IDEFIX_CMAKE_OPTIONS+=( -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_HOPPER90=ON );;
        '') ;;
        *)
            printf "Error: unknown gpu architecture '%s'\n" "$model"
            return 1
    esac
}

function in_env() {
    local -a cmd=( "$@" )
    local drvfile="$IDEFIX_DIR/scripts/etc/env-bigfoot.drv"
    if [ ! -e "$drvfile" ]; then
        printf "Cacheing an Idefix shell derivation in %s\n" "$drvfile"
        nix-instantiate --add-root "$drvfile" "$IDEFIX_DIR/scripts/etc/env-bigfoot.nix"
    fi
    NIX_BUILD_SHELL="$HOME/.nix-shell/bash" nix-shell "$drvfile" --run "$(declare -p cmd); "'"${cmd[@]}"'
}
