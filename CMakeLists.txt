cmake_minimum_required(VERSION 3.18)
project(picsp3d VERSION 3.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(PICSP_NATIVE "Optimise for the build machine (-march=native)" ON)
option(PICSP_TESTS "Build tests" ON)
option(PICSP_WERROR "Treat compiler warnings as errors (CI)" OFF)
option(PICSP_SANITIZE "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(PICSP_FAST_MATH "Enable -ffast-math style optimisations in particle kernels" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(MPI_DETERMINE_LIBRARY_VERSION ON)   # fills MPI_CXX_LIBRARY_VERSION_STRING, checked below
find_package(MPI REQUIRED COMPONENTS CXX)
find_package(OpenMP REQUIRED COMPONENTS CXX)

# Identify the MPI library and the launcher separately: with conda or modules on PATH they can come
# from different implementations, and an Open MPI binary started by MPICH's mpiexec (or vice versa)
# runs every rank as rank 0 of 1.
function(picsp_mpi_impl text out)
  if(text MATCHES "Open MPI|OpenRTE")
    set(${out} openmpi PARENT_SCOPE)
  elseif(text MATCHES "MPICH|HYDRA|Intel")
    set(${out} mpich PARENT_SCOPE)
  else()
    set(${out} unknown PARENT_SCOPE)
  endif()
endfunction()
picsp_mpi_impl("${MPI_CXX_LIBRARY_VERSION_STRING}" PICSP_MPILIB_IMPL)
execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version OUTPUT_VARIABLE _mpiexec_version
                ERROR_VARIABLE _mpiexec_version RESULT_VARIABLE _mpiexec_result TIMEOUT 10)
picsp_mpi_impl("${_mpiexec_version}" PICSP_MPIEXEC_IMPL)
message(STATUS "MPI library: ${PICSP_MPILIB_IMPL} (${MPI_CXX_COMPILER}), launcher: ${PICSP_MPIEXEC_IMPL} (${MPIEXEC_EXECUTABLE})")
if(NOT PICSP_MPILIB_IMPL STREQUAL "unknown" AND NOT PICSP_MPIEXEC_IMPL STREQUAL "unknown"
   AND NOT PICSP_MPILIB_IMPL STREQUAL PICSP_MPIEXEC_IMPL)
  message(WARNING "The MPI library (${PICSP_MPILIB_IMPL}) and mpiexec (${PICSP_MPIEXEC_IMPL}) come from "
                  "different MPI implementations. Point both at the same installation, e.g. "
                  "-DMPI_CXX_COMPILER=/usr/bin/mpicxx -DMPIEXEC_EXECUTABLE=/usr/bin/mpirun")
endif()

file(GLOB_RECURSE PICSP_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/*.cpp)
list(REMOVE_ITEM PICSP_SOURCES ${CMAKE_SOURCE_DIR}/src/main.cpp)

add_library(picsp_core STATIC ${PICSP_SOURCES})
target_include_directories(picsp_core PUBLIC ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/extern)
target_link_libraries(picsp_core PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX)
target_compile_options(picsp_core PUBLIC
  $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wno-unused-parameter -Wno-sign-compare>
  $<$<CONFIG:Release>:-O3>)
if(PICSP_NATIVE)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-march=native HAS_MARCH_NATIVE)
  if(HAS_MARCH_NATIVE)
    target_compile_options(picsp_core PUBLIC -march=native)
  endif()
endif()
if(PICSP_FAST_MATH)
  target_compile_options(picsp_core PUBLIC -fno-math-errno -fno-trapping-math)
endif()
if(PICSP_WERROR)
  target_compile_options(picsp_core PUBLIC -Werror)
endif()
if(PICSP_SANITIZE)
  target_compile_options(picsp_core PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(picsp_core PUBLIC -fsanitize=address,undefined)
endif()
target_compile_definitions(picsp_core PUBLIC PICSP_VERSION="${PROJECT_VERSION}")

add_executable(picsp3d ${CMAKE_SOURCE_DIR}/src/main.cpp)
target_link_libraries(picsp3d PRIVATE picsp_core)

if(PICSP_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

install(TARGETS picsp3d RUNTIME DESTINATION bin)
