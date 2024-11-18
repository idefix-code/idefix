let nixpkgs = import ./env-bigfoot-nigpkgs.nix;
    
    cudatoolkit = nixpkgs.cudaPackages_12.cudatoolkit;

    libInputs = with nixpkgs; [
      stdenv.cc.cc
      (nur.repos.gricad.openmpi4.override {
        cudaSupport = true;
        inherit cudatoolkit;
      })
      nur.repos.gricad.ucx
      rocmPackages.rocm-smi
      rocmPackages.clr
      rocmPackages.rocthrust
      rocmPackages.rocprim

      zlib
      cudatoolkit
    ];
in nixpkgs.mkShell {
  buildInputs = with nixpkgs; [
    cmake
    pkg-config
  ];
  shellHook = ''
    export LD_LIBRARY_PATH="${nixpkgs.lib.makeLibraryPath libInputs}:/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH"
  '';
  NIX_SHELL_PROMPT_TAG = "idefix";
  IDEFIX_CUDA_INCLUDE = "${nixpkgs.lib.getDev cudatoolkit}/include";
}
