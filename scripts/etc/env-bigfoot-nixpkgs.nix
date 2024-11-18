let # doit correspondre a la glibc du systeme (2.31 a compter du Fri Nov 15 11:27:37 CET 2024)
    nixpkgs_bigfoot = import (builtins.fetchGit {
      name = "nixpkgs_glibc_2_31";
      url = "https://github.com/NixOS/nixpkgs/";
      ref = "refs/heads/nixpkgs-unstable";
      rev = "3913f6a514fa3eb29e34af744cc97d0b0f93c35c";
    }) {};

    glibc_bigfoot = nixpkgs_bigfoot.glibc;
    gcc_bigfoot = nixpkgs_bigfoot.gcc-unwrapped;

    getCustomGccStdenv = customGcc: customGlibc: pkgs:
      with pkgs; let
        compilerWrapped = wrapCCWith {
          cc = customGcc;
          bintools = wrapBintoolsWith {
            bintools = binutils-unwrapped;
            libc = customGlibc;
          };
        };
      in
        overrideCC stdenv compilerWrapped;
        
    withBigfootStdenv = _final: prev: {
      config.replaceStdenv = {pkgs}: getCustomGccStdenv gcc_bigfoot glibc_bigfoot pkgs;
      config.allowAliases = true;
      config.rocmSupport = true;
      config.cudaSupport = true;
      config.packageOverrides = pkgs: {
        nur = import (builtins.fetchTarball "https://github.com/nix-community/NUR/archive/master.tar.gz") {
          inherit pkgs;
        };
      };
    };
in import <nixpkgs> {
  overlays = [ withBigfootStdenv ];
};
