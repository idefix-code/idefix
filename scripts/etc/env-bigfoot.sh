#!/usr/bin/env bash
function setup_env() {
    set +ue
    source /applis/site/nix.sh >/dev/null 2>&1
    set -ue

    local cache=~/.nix-cache/nix-shell-idefix
    if [ ! -e "$cache" ]; then
        mkdir -p ~/.nix-cache
	nix-instantiate --add-root "$cache.drv" --expr 'with import ./scripts/etc/env-bigfoot-nixpkgs.nix; bashInteractive'
	nix-store --realise --add-root "$cache" "$cache.drv" 
    fi
}

function set_gpu_options() {
    local model="$1"; shift
    case "$model" in
        V100)
            IDEFIX_CMAKE_OPTIONS+=( -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_VOLTA70=ON "-DCUDAToolkit_INCLUDE_DIR=$(in_env_raw 'echo $IDEFIX_CUDA_INCLUDE')" );;
        A100)
            IDEFIX_CMAKE_OPTIONS+=( -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_AMPERE80=ON "-DCUDAToolkit_INCLUDE_DIR=$(in_env_raw 'echo $IDEFIX_CUDA_INCLUDE')" );;
        H100)
            IDEFIX_CMAKE_OPTIONS+=( -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_HOPPER90=ON "-DCUDAToolkit_INCLUDE_DIR=$(in_env_raw 'echo $IDEFIX_CUDA_INCLUDE')" );;
        '') ;;
        *)
            printf "Error: unknown gpu architecture '%s'\n" "$model"
            return 1
    esac
}

function in_env_raw() {
    local -a cmd="$1"
    local drvfile="$IDEFIX_DIR/scripts/etc/env-bigfoot.drv"
    if [ ! -e "$drvfile" ]; then
        printf "Cacheing an Idefix shell derivation in %s\n" "$drvfile" >&2
        nix-instantiate --add-root "$drvfile" "$IDEFIX_DIR/scripts/etc/env-bigfoot.nix" >/dev/null
    fi
    printf "Running command: %s\n" "$cmd" >&2
    NIX_BUILD_SHELL="$HOME/.nix-cache/nix-shell/bin/bash" nix-shell "$drvfile" --run "$cmd"
}
function in_env() {
    local -a cmd=( "$@" )
    in_env_raw "$(declare -p cmd); "'"${cmd[@]}"'
}
