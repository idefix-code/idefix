let nixpkgs = import <nixpkgs> {};
    cudart = nixpkgs.lib.getDev nixpkgs.cudaPackages_12.cuda_cudart;

    inputs = with nixpkgs; [
      nur.repos.gricad.openmpi4
      nur.repos.gricad.ucx
      rocmPackages.rocm-smi
      rocmPackages.clr
      rocmPackages.rocthrust
      rocmPackages.rocprim
      
      cudaPackages_12.cuda_nvcc
      cudart

      cmake
      pkg-config
    ];
in nixpkgs.mkShell {
  buildInputs = inputs;
  LD_LIBRARY_PATH = nixpkgs.lib.makeLibraryPath inputs;
  NIX_SHELL_PROMPT_TAG = "idefix";
  IDEFIX_CUDA_INCLUDE = "${cudart}/include";
}
