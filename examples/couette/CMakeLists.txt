# Couette example — Palabos case.
# Build: mkdir build && cd build && cmake .. -DPALABOS_ROOT=<path> && make
cmake_minimum_required(VERSION 3.10)
project(couette)
set(CMAKE_CXX_STANDARD 14)
if(NOT DEFINED PALABOS_ROOT)
  set(PALABOS_ROOT "$ENV{PALABOS_ROOT}" CACHE PATH "Palabos root directory")
endif()
if(NOT PALABOS_ROOT)
  message(FATAL_ERROR "Set -DPALABOS_ROOT=<path-to-palabos> or PALABOS_ROOT env var")
endif()
include(${PALABOS_ROOT}/cmake/Palabos.cmake)
add_executable(couette couette.cpp)
target_link_libraries(couette palabos)
