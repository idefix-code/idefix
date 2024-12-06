let nixpkgs = import ./env-bigfoot-nixpkgs.nix;
    
    cudatoolkit = nixpkgs.cudaPackages_12.cudatoolkit;

    libInputs = with nixpkgs; [
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
    rocmPackages.hipcc
  ];
  NIX_SHELL_PROMPT_TAG = "idefix";
  IDEFIX_CUDA_INCLUDE = "${nixpkgs.lib.getDev cudatoolkit}/include";
}
