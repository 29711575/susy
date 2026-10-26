#!/usr/bin/env bash
# CMake-free fallback build that mirrors the GNU branch of CMakeLists.txt.
# Use this only when cmake / the Intel classic toolchain is unavailable.
# Requires: g++ (C++17), gfortran (for the PFPACK Fortran sources), MKL
# (mkl_intel_lp64 + mkl_sequential + mkl_core), Eigen and Boost headers.
#
#   EIGEN3_INCLUDE_DIR=... BOOST_INCLUDE_DIR=... MKLROOT=... \
#     bash scripts/build_nodeps.sh [build_dir] [target ...]
#
# With no target arguments every kitaev_ssh program and validation binary in
# this source tree is built. The legacy reproduction/ targets are available
# through CMake and are not selected by this SSH-focused fallback script.
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MKL_DYNAMIC=FALSE
root=$(cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"$root/build"}; shift || true
: "${EIGEN3_INCLUDE_DIR:?Set EIGEN3_INCLUDE_DIR to the directory containing Eigen/}"
: "${BOOST_INCLUDE_DIR:?Set BOOST_INCLUDE_DIR to the directory containing boost/}"
: "${MKLROOT:?Source oneAPI setvars.sh or set MKLROOT}"
CXX=${CXX:-g++}; CC=${CC:-gcc}; FC=${FC:-gfortran}
mkdir -p "$build/pfapack"; build=$(cd "$build" && pwd)

# 1) PFPACK static libraries (serial build, one process at a time).
fortran_objects=()
for file in "$root"/inc/pfapack/fortran/*.f; do
  object="$build/pfapack/$(basename "${file%.f}").o"
  if [ ! -s "$object" ] || [ "$file" -nt "$object" ]; then "$FC" -O3 -c "$file" -o "$object"; fi
  fortran_objects+=("$object")
done
(cd "$build/pfapack" && ar rcs libpfapack.a "${fortran_objects[@]}")
c_objects=()
for file in "$root"/inc/pfapack/c_interface/*.c; do
  object="$build/pfapack/c_$(basename "${file%.c}").o"
  if [ ! -s "$object" ] || [ "$file" -nt "$object" ]; then "$CC" -O3 -c "$file" -o "$object"; fi
  c_objects+=("$object")
done
(cd "$build/pfapack" && ar rcs libcpfapack.a "${c_objects[@]}")

# 2) GNU/MKL flags, kept identical in spirit to the CMakeLists GNU branch.
flags=(-O3 -DNDEBUG -std=c++17 -static-libstdc++ -static-libgcc -DPFQMC_SCALE_SAFE_UDT
       -no-pie -fno-fast-math -ffp-contract=off
       -I"$EIGEN3_INCLUDE_DIR" -I"$BOOST_INCLUDE_DIR" -I"$root/inc" -I"$MKLROOT/include")
libs=(-L"$MKLROOT/lib" -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
      -lpthread -lm -ldl -lgfortran)

# 3) program -> source mapping (only files present in this tree).
declare -A PROGS=(
  [kitaev_ssh_center_green_validate]=kitaev_ssh/center_green_validate.cpp
  [kitaev_ssh_guided_validate]=kitaev_ssh/guided_proposal_validate.cpp
  [kitaev_ssh_pfqmc_free_phonon_validate]=kitaev_ssh/pfqmc_free_phonon_validate.cpp
  [kitaev_ssh_interacting_fock_validate]=kitaev_ssh/interacting_fock_validate.cpp
  [kitaev_ssh_total_energy_validate]=kitaev_ssh/total_energy_validate.cpp
  [kitaev_ssh_pfqmc_fock_validate]=kitaev_ssh/pfqmc_fock_validate.cpp
  [kitaev_ssh_pfqmc_cache_validate]=kitaev_ssh/pfqmc_cache_validate.cpp
  [kitaev_ssh_vbs_field_validate]=kitaev_ssh/vbs_field_validate.cpp
  [kitaev_ssh_green_precision_benchmark]=kitaev_ssh/green_precision_benchmark.cpp
  [kitaev_ssh_green_refinement_validate]=kitaev_ssh/green_refinement_validate.cpp
  [kitaev_ssh_vbs_susceptibility_validate]=kitaev_ssh/vbs_susceptibility_validate.cpp
  [kitaev_ssh_vbs_amplitude_validate]=kitaev_ssh/vbs_amplitude_validate.cpp
  [kitaev_ssh_real_magnitude_validate]=kitaev_ssh/real_magnitude_validate.cpp
  [kitaev_ssh_physical_trial_validate]=kitaev_ssh/physical_matter_trial_validate.cpp
  [kitaev_ssh_ported_validate]=kitaev_ssh/ported_spectral_validate.cpp
  [kitaev_ssh_free]=kitaev_ssh/free_projector.cpp
  [kitaev_ssh_phonon_example]=kitaev_ssh/phonon_example.cpp
  [kitaev_ssh_mc]=kitaev_ssh/coupled_mc.cpp
  [kitaev_ssh_profile]=kitaev_ssh/profile_chain.cpp
  [kitaev_ssh_production_profile]=kitaev_ssh/production_profile.cpp
  [kitaev_ssh_weight_benchmark]=kitaev_ssh/weight_benchmark.cpp
  [kitaev_ssh_overlap_svd_validate]=kitaev_ssh/overlap_svd_validate.cpp
  [kitaev_ssh_bipartite_validate]=kitaev_ssh/bipartite_spectrum_validate.cpp
  [kitaev_ssh_fermion_validate]=kitaev_ssh/fermion_validate.cpp
  [kitaev_ssh_weight_validate]=kitaev_ssh/weight_optimization_validate.cpp
  [kitaev_ssh_validate]=kitaev_ssh/validate.cpp
  [kitaev_ssh_phonon_validate]=kitaev_ssh/phonon_validate.cpp
  [kitaev_ssh_numerical_failure_validate]=kitaev_ssh/numerical_failure_validate.cpp
  [kitaev_ssh_deferred_phase_validate]=kitaev_ssh/deferred_phase_validate.cpp
  [kitaev_ssh_gamma_validate]=kitaev_ssh/gamma_validate.cpp
  [kitaev_ssh_q_mode_validate]=kitaev_ssh/q_mode_validate.cpp
  [kitaev_ssh_mc_validate]=kitaev_ssh/mc_validate.cpp
  [kitaev_ssh_shifted_validate]=kitaev_ssh/shifted_trial_validate.cpp
  [kitaev_ssh_vbs_validate]=kitaev_ssh/vbs_validate.cpp
  [kitaev_ssh_vbs_q_validate]=kitaev_ssh/vbs_q_validate.cpp
  [kitaev_ssh_sheared_pbc_validate]=kitaev_ssh/sheared_pbc_validate.cpp
)
requested=("$@")
[ ${#requested[@]} -eq 0 ] && requested=("${!PROGS[@]}")
for name in "${requested[@]}"; do
  src="$root/${PROGS[$name]:-}"
  [ -n "$src" ] || { echo "unknown target $name" >&2; exit 2; }
  [ -f "$src" ] || { echo "skip $name (missing $src)" >&2; continue; }
  echo "building $name"
  "$CXX" "${flags[@]}" "$src" "$root/src/skewMatUtils.cpp" \
    "$build/pfapack/libcpfapack.a" "$build/pfapack/libpfapack.a" "${libs[@]}" -o "$build/$name"
done
"$CXX" --version > "$build/compiler.txt"
echo "Built in $build"
