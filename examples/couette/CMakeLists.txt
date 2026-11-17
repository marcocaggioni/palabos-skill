# Couette example — Palabos case.
# Build (serial, no MPI):
#   mkdir -p build && cd build && cmake .. -DPALABOS_ROOT=$HOME/palabos && make -j4
# Run: ./couette --gamma 10 --out runs/gamma_10
# Follows the standalone pattern of Palabos' own showcase examples
# (examples/showCases/cavity2d/CMakeLists.txt): the Palabos library is
# built from its sources (no installed Palabos package exists).
cmake_minimum_required(VERSION 3.10)
project(couette)
enable_language(CXX)

if(NOT DEFINED PALABOS_ROOT)
  set(PALABOS_ROOT "$ENV{PALABOS_ROOT}" CACHE PATH "Palabos root directory")
endif()
if(NOT PALABOS_ROOT)
  message(FATAL_ERROR "Set -DPALABOS_ROOT=<path-to-palabos> or PALABOS_ROOT env var")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(CMAKE_CXX_STANDARD 20)
if(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
  set(CMAKE_CXX_FLAGS "-std=c++20 -Wall -Wextra -Wnon-virtual-dtor")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
endif()

option(ENABLE_MPI "Enable MPI" OFF)
if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  include_directories(${MPI_CXX_INCLUDE_PATH})
  add_definitions(-DPLB_MPI_PARALLEL)
endif()

add_definitions(-DPLB_USE_POSIX)

include_directories("${PALABOS_ROOT}/src")
include_directories("${PALABOS_ROOT}/externalLibraries")
include_directories("${PALABOS_ROOT}/externalLibraries/Eigen3")

file(GLOB_RECURSE PALABOS_SRC "${PALABOS_ROOT}/src/*.cpp")
file(GLOB_RECURSE EXT_SRC "${PALABOS_ROOT}/externalLibraries/tinyxml/*.cpp")
add_library(palabos STATIC ${PALABOS_SRC} ${EXT_SRC})

# Binary lands in build/ (couette_sweep.py expects build/couette).
add_executable(couette couette.cpp)
target_link_libraries(couette palabos)
if(ENABLE_MPI)
  target_link_libraries(couette ${MPI_CXX_LIBRARIES})
endif()
find_package(TBB REQUIRED)
target_link_libraries(couette TBB::tbb)
