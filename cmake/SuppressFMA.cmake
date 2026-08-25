#[[============================================================================
  SuppressFMA.cmake

  Provides an option and a function to optionally disable fused
  multiply-add (FMA) code generation / contraction for a Kokkos-based CXX
  target.

  Kokkos wraps the real device compiler behind nvcc_wrapper (CUDA) or
  hipcc (HIP), which makes CMAKE_CXX_COMPILER_ID report the *underlying
  host* compiler (e.g. "GNU" or "Clang") instead of "NVIDIA" or "Clang
  as HIP". Backend detection therefore relies on the Kokkos_ENABLE_*
  variables exported by KokkosConfig.cmake / set by the Kokkos build,
  and only falls back to CMAKE_CXX_COMPILER_ID for the plain host
  compilers (no CUDA/HIP backend active):

    - Kokkos_ENABLE_HIP  ON -> AMD HIP  (hipcc, clang-based)  : -ffp-contract=off
    - Kokkos_ENABLE_CUDA ON -> NVIDIA nvcc (via nvcc_wrapper)  : --fmad=false
    - otherwise, CMAKE_CXX_COMPILER_ID selects among:
        GNU        : gcc/g++
        Intel      : classic icc/icpc
        IntelLLVM  : Intel oneAPI icx/icpx
        Clang      : LLVM clang++
        AppleClang : Xcode clang++
        NVHPC      : NVIDIA HPC SDK (nvc++), e.g. for OpenMPTarget/OpenACC

  Usage (after find_package(Kokkos) so Kokkos_ENABLE_* are defined):
    include(SuppressFMA.cmake)
    add_library(mylib source.cpp)
    target_link_libraries(mylib PUBLIC Kokkos::kokkos)
    target_suppress_fma(mylib)

============================================================================]]

include_guard(GLOBAL)

# Determine the CXX FMA-suppression flags for the active Kokkos backend /
# CXX compiler. Returns the list of flags (possibly empty) via out_var.
function(_fma_suppression_flags_cxx out_var)
    set(flags "")
    set(id "${CMAKE_CXX_COMPILER_ID}")

    # Kokkos backend takes priority: nvcc_wrapper/hipcc hide the real
    # device compiler from CMAKE_CXX_COMPILER_ID.
    if(Kokkos_ENABLE_HIP)
        set(flags "-ffp-contract=off")

    elseif(Kokkos_ENABLE_CUDA)
        set(flags "--fmad=false")

    elseif(id STREQUAL "GNU")
        set(flags "-ffp-contract=off" "-mno-fma")

    elseif(id MATCHES "^(Clang|AppleClang)$")
        set(flags "-ffp-contract=off")

    elseif(id STREQUAL "Intel")
        # Intel classic compiler
        set(flags "-fp-model=precise" "-no-fma")

    elseif(id STREQUAL "IntelLLVM")
        # Intel oneAPI compiler (clang-based)
        set(flags "-ffp-contract=off" "-fp-model=strict")

    elseif(id STREQUAL "NVHPC")
        # NVIDIA HPC SDK (formerly PGI), e.g. OpenMPTarget/OpenACC backend
        set(flags "-Mnofma") # untested
    endif()

    set(${out_var} "${flags}" PARENT_SCOPE)
endfunction()

# target_suppress_fma(<target>)
#
# Applies compiler-specific FMA-suppression flags to <target>'s CXX
# sources, but only when the SUPPRESS_FMA option is ON. Safe to call
# unconditionally.
function(target_suppress_fma target)

    if(NOT TARGET ${target})
        message(FATAL_ERROR "target_suppress_fma: '${target}' is not a target")
    endif()

    _fma_suppression_flags_cxx(cxx_flags)

    if(cxx_flags)
        foreach(flag IN LISTS cxx_flags)
            target_compile_options(${target} PRIVATE
                $<$<COMPILE_LANGUAGE:CXX>:${flag}>
            )
        endforeach()
    else()
        message(VERBOSE
            "target_suppress_fma: no FMA-suppression flag known for "
            "CXX compiler '${CMAKE_CXX_COMPILER_ID}' "
            "(Kokkos_ENABLE_CUDA=${Kokkos_ENABLE_CUDA}, "
            "Kokkos_ENABLE_HIP=${Kokkos_ENABLE_HIP}) (target ${target})")
    endif()
endfunction()
