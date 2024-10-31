let nixpkgs = import <nixpkgs> {};

    inputs = with nixpkgs; [
      nur.repos.gricad.openmpi4
      nur.repos.gricad.ucx
      rocmPackages.rocm-smi
      rocmPackages.clr
      rocmPackages.rocthrust
      rocmPackages.rocprim
      
      cudaPackages_12.cudatoolkit

      cmake
      pkg-config
    ];
in nixpkgs.mkShell {
  buildInputs = inputs;
  LD_LIBRARY_PATH = nixpkgs.lib.makeLibraryPath inputs;
  NIX_SHELL_PROMPT_TAG = "idefix";
  IDEFIX_CUDA_INCLUDE = "${cudaPackages_12.cudatoolkit}/include";
}
