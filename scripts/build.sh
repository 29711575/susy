#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MKL_DYNAMIC=FALSE
root=$(cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"$root/build"}
mkdir -p "$build"
build=$(cd "$build" && pwd)
: "${EIGEN3_INCLUDE_DIR:?Set EIGEN3_INCLUDE_DIR to the directory containing Eigen/}"
: "${BOOST_INCLUDE_DIR:?Set BOOST_INCLUDE_DIR to the directory containing boost/}"
CXX=${CXX:-g++}; CC=${CC:-icx}; FC=${FC:-ifx}
mkdir -p "$build/pfapack"
cd "$build/pfapack"
# Serial compilation: at most one compiler process at a time.
fortran_objects=()
for file in "$root"/inc/pfapack/fortran/*.f; do
 object="$(basename "${file%.f}").o"
 "$FC" -O3 -c "$file" -o "$object"
 fortran_objects+=("$object")
done
rm -f libpfapack.a libcpfapack.a
ar rcs libpfapack.a "${fortran_objects[@]}"
objects=()
for file in "$root"/inc/pfapack/c_interface/*.c; do
 object="c_$(basename "${file%.c}").o"
 "$CC" -O3 -c "$file" -o "$object"
 objects+=("$object")
done
ar rcs libcpfapack.a "${objects[@]}"
libs=()
case "$(basename "$CXX")" in
 g++*) flags=(-no-pie -fno-fast-math -ffp-contract=off -I"${MKLROOT:?Source oneAPI setvars.sh}/include"); libs=(-L"$MKLROOT/lib" -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl "$("$FC" -print-file-name=libirc.a)");;
 icpc) flags=(-mkl=sequential -fp-model=precise -fimf-arch-consistency=true -no-fma);;
 *) flags=(-Wno-enum-constexpr-conversion -qmkl=sequential -fp-model=precise -ffp-contract=off);;
esac
for name in driven_qmc density_sign_regression cdw_initial_test hs1_review_test; do
 target=$name
 if [[ $name == driven_qmc ]]; then target=pure_state_driven_qmc; fi
 "$CXX" -O3 -DNDEBUG -std=c++17 -DPFQMC_SCALE_SAFE_UDT "${flags[@]}" \
  -I"$EIGEN3_INCLUDE_DIR" -I"$BOOST_INCLUDE_DIR" -I"$root/inc" \
  "$root/reproduction/pure_state_driven_stage1/$name.cpp" "$root/src/skewMatUtils.cpp" \
  "$build/pfapack/libcpfapack.a" "$build/pfapack/libpfapack.a" "${libs[@]}" -o "$build/$target"
done
"$CXX" --version > "$build/compiler.txt"
echo "Built in $build"
