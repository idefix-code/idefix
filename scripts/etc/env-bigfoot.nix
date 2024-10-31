let nixpkgs = import <nixpkgs> { };
    cudatoolkit = nixpkgs.cudaPackages_12.cudatoolkit;
    
    inputs = with nixpkgs; [
      (nur.repos.gricad.openmpi4.override {
        cudaSupport = true;
        inherit cudatoolkit;
      })
      nur.repos.gricad.ucx
      rocmPackages.rocm-smi
      rocmPackages.clr
      rocmPackages.rocthrust
      rocmPackages.rocprim
      
      cudatoolkit

      cmake
      pkg-config
    ];
in nixpkgs.mkShell {
  buildInputs = inputs;
  LD_LIBRARY_PATH = nixpkgs.lib.makeLibraryPath inputs;
  NIX_SHELL_PROMPT_TAG = "idefix";
  IDEFIX_CUDA_INCLUDE = "${nixpkgs.lib.getDev cudatoolkit}/include";
}
