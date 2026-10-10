cmake_minimum_required(VERSION 3.18)

# Configure-only regression for the MKL discovery cache. The fake files prove
# path selection and cache invalidation; they are deliberately never linked.
get_filename_component(_repo "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef _suffix)
set(_work "/tmp/pfqmc-mklroot-regression-${_suffix}")
file(MAKE_DIRECTORY "${_work}")

set(_deps "${_work}/deps")
file(MAKE_DIRECTORY "${_deps}/Eigen" "${_deps}/boost/serialization")
file(WRITE "${_deps}/Eigen/Core" "// fake Eigen include for configure-only testing\n")
file(WRITE "${_deps}/boost/serialization/nvp.hpp" "// fake Boost include for configure-only testing\n")

function(make_fake_root root layout)
  file(MAKE_DIRECTORY "${root}/include")
  if(NOT "${ARGN}" STREQUAL "missing-header")
    file(WRITE "${root}/include/mkl_types.h" "/* fake MKL header */\n")
  endif()
  if("${layout}" STREQUAL "intel64")
    set(_libdir "${root}/lib/intel64")
  else()
    set(_libdir "${root}/lib")
  endif()
  file(MAKE_DIRECTORY "${_libdir}")
  foreach(_component intel_lp64 sequential core)
    if(NOT "${ARGN}" STREQUAL "missing-${_component}")
      file(WRITE "${_libdir}/libmkl_${_component}.so" "fake ${_component}\n")
    endif()
  endforeach()
endfunction()

set(_a "${_work}/mkl-A")
set(_b "${_work}/mkl-B")
set(_missing_header "${_work}/mkl-missing-header")
set(_missing_intel "${_work}/mkl-missing-intel")
set(_missing_sequential "${_work}/mkl-missing-sequential")
set(_missing_core "${_work}/mkl-missing-core")
make_fake_root("${_a}" direct)
make_fake_root("${_b}" intel64 intel64-marker)
make_fake_root("${_missing_header}" direct missing-header)
make_fake_root("${_missing_intel}" direct missing-intel_lp64)
make_fake_root("${_missing_sequential}" direct missing-sequential)
make_fake_root("${_missing_core}" direct missing-core)
set(_a_link "${_work}/mkl-A-symlink")
file(CREATE_LINK "${_a}" "${_a_link}" SYMBOLIC)

function(configure_ok build root explicit env_root)
  set(_args
    "${CMAKE_COMMAND}" -S "${_repo}" -B "${build}"
    "-DEIGEN3_INCLUDE_DIR=${_deps}"
    "-DBOOST_INCLUDE_DIR=${_deps}"
    "-DBUILD_TESTING=OFF")
  if(explicit)
    list(APPEND _args "-DMKLROOT=${root}")
  endif()
  if("${env_root}" STREQUAL "")
    execute_process(COMMAND ${_args} RESULT_VARIABLE _result
      OUTPUT_VARIABLE _stdout ERROR_VARIABLE _stderr)
  else()
    execute_process(COMMAND "${CMAKE_COMMAND}" -E env "MKLROOT=${env_root}" ${_args}
      RESULT_VARIABLE _result OUTPUT_VARIABLE _stdout ERROR_VARIABLE _stderr)
  endif()
  if(NOT _result EQUAL 0)
    message(FATAL_ERROR "configure failed for ${root}:\n${_stdout}\n${_stderr}")
  endif()
endfunction()

function(configure_fail build root)
  set(_args
    "${CMAKE_COMMAND}" -S "${_repo}" -B "${build}"
    "-DEIGEN3_INCLUDE_DIR=${_deps}"
    "-DBOOST_INCLUDE_DIR=${_deps}"
    "-DBUILD_TESTING=OFF"
    "-DMKLROOT=${root}")
  execute_process(COMMAND ${_args} RESULT_VARIABLE _result
    OUTPUT_VARIABLE _stdout ERROR_VARIABLE _stderr)
  if(_result EQUAL 0)
    message(FATAL_ERROR "invalid MKL root unexpectedly configured: ${root}")
  endif()
  string(CONCAT _diagnostic "${_stdout}\n${_stderr}")
  if(NOT _diagnostic MATCHES "MKL")
    message(FATAL_ERROR "MKL failure did not identify the requested root: ${root}\n${_diagnostic}")
  endif()
endfunction()

function(assert_cache_root build expected)
  get_filename_component(_expected_real "${expected}" REALPATH)
  file(READ "${build}/CMakeCache.txt" _cache)
  string(FIND "${_cache}" "MKLROOT:PATH=${_expected_real}" _root_pos)
  if(_root_pos EQUAL -1)
    message(FATAL_ERROR "cache root mismatch in ${build}; expected ${_expected_real}")
  endif()
  foreach(_component INTEL_LP64 SEQUENTIAL CORE)
    string(FIND "${_cache}" "MKL_${_component}_LIBRARY:FILEPATH=${_expected_real}/" _lib_pos)
    if(_lib_pos EQUAL -1)
      message(FATAL_ERROR "cache library ${_component} did not resolve below ${_expected_real}")
    endif()
  endforeach()
  string(FIND "${_cache}" "MKL_INCLUDE_DIR:PATH=${_expected_real}/include" _include_pos)
  if(_include_pos EQUAL -1)
    message(FATAL_ERROR "cache include path did not resolve below ${_expected_real}")
  endif()
endfunction()

set(_main "${_work}/build-main")
configure_ok("${_main}" "${_a}" TRUE "")
assert_cache_root("${_main}" "${_a}")
configure_ok("${_main}" "${_a}" TRUE "")
assert_cache_root("${_main}" "${_a}")

# Same build directory, both directions. This is the stale-library regression.
configure_ok("${_main}" "${_b}" TRUE "")
assert_cache_root("${_main}" "${_b}")
configure_ok("${_main}" "${_a}" TRUE "")
assert_cache_root("${_main}" "${_a}")

set(_env_build "${_work}/build-environment")
configure_ok("${_env_build}" "${_b}" FALSE "${_b}")
assert_cache_root("${_env_build}" "${_b}")
set(_explicit_build "${_work}/build-explicit-priority")
configure_ok("${_explicit_build}" "${_b}" TRUE "${_a}")
assert_cache_root("${_explicit_build}" "${_b}")

set(_link_build "${_work}/build-symlink")
configure_ok("${_link_build}" "${_a_link}" TRUE "")
assert_cache_root("${_link_build}" "${_a}")

configure_fail("${_work}/build-missing-header" "${_missing_header}")
configure_fail("${_work}/build-missing-intel" "${_missing_intel}")
configure_fail("${_work}/build-missing-sequential" "${_missing_sequential}")
configure_fail("${_work}/build-missing-core" "${_missing_core}")
configure_fail("${_work}/build-nonexistent" "${_work}/mkl-does-not-exist")

message(STATUS "MKLROOT_CACHE_REGRESSION_OK")
file(REMOVE_RECURSE "${_work}")
