cmake_minimum_required(VERSION 3.21)
project(BECpp)

set(CMAKE_CXX_STANDARD 20)

set(SOURCES src/main.cpp src/grid.cpp src/wavefunction.cpp)
set(INCLUDES include/constants.h include/grid.h include/wavefunction.h)
set(FFTW_PATH C:/dev/cpp/libs/fftwcmake)

add_executable(${PROJECT_NAME} ${SOURCES} ${INCLUDES})

target_include_directories(${PROJECT_NAME} PRIVATE include/ ${FFTW_PATH})
target_link_directories(${PROJECT_NAME} PRIVATE ${FFTW_PATH})
target_link_libraries(${PROJECT_NAME} libfftw3-3)
