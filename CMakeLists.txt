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

find_package(MPI REQUIRED COMPONENTS CXX)
find_package(OpenMP REQUIRED COMPONENTS CXX)

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
