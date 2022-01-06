cmake_minimum_required(VERSION 3.21)
project(BECpp)

set(CMAKE_CXX_STANDARD 20)

set(SOURCES src/main.cpp src/grid.cpp src/wavefunction.cpp src/data.cpp)
set(INCLUDES include/constants.h include/grid.h include/wavefunction.h include/data.h)

add_executable(${PROJECT_NAME} ${SOURCES} ${INCLUDES})

target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include/
        ${PROJECT_SOURCE_DIR}/lib/fftw ${PROJECT_SOURCE_DIR}/lib/HighFive/include
        ${PROJECT_SOURCE_DIR}/lib/HDF5/1.13.0/include)
target_link_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/lib/fftw
        ${PROJECT_SOURCE_DIR}/lib/HDF5/1.13.0/lib)
target_link_libraries(${PROJECT_NAME} libfftw3-3 hdf5)
