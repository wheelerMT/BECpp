cmake_minimum_required(VERSION 3.21)
project(BECpp)

set(CMAKE_CXX_STANDARD 20)

set(SOURCES src/main.cpp src/grid.cpp src/wavefunction.cpp)
set(INCLUDES include/constants.h include/grid.h include/wavefunction.h)

add_executable(BECpp ${SOURCES} ${INCLUDES})

target_include_directories(BECpp PRIVATE include/)
