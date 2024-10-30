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
    local drvfile="$IDEFIX_DIR/scripts/etc/env-bigfoot.drv"
    if [ ! -e "$drvfile" ]; then
        printf "Cacheing an Idefix shell derivation in %s\n" "$drvfile"
        nix-instantiate --add-root "$drvfile" "$IDEFIX_DIR/scripts/etc/env-bigfoot.nix"
    fi
    NIX_BUILD_SHELL="$HOME/.nix-shell/bash" nix-shell "$drvfile" --run "$(declare -p cmd); "'"${cmd[@]}"'
}
