cmake_minimum_required(VERSION 3.21..3.26)
project(BECpp)

set(CMAKE_CXX_STANDARD 20)

set(SOURCES src/grid.cpp src/wavefunction.cpp src/data.cpp src/evolution.cpp)
set(INCLUDES include/constants.h include/grid.h include/wavefunction.h
  include/data.h include/evolution.h include/BECpp.h)

find_package(OpenMP REQUIRED)
find_package(HDF5 REQUIRED)

set(USE_BOOST OFF CACHE BOOL "Enable Boost Support")
option(HIGHFIVE_EXAMPLES "Compile examples" OFF)
option(HIGHFIVE_BUILD_DOCS "Enable documentation building" OFF)
add_subdirectory(lib/HighFive)

configure_file(cmake/downloadFindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download)
if (result)
    message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
else ()
    message("CMake step for findFFTW completed (${result}).")
endif ()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download)
if (result)
    message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
endif ()

set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")
find_package(FFTW REQUIRED)

add_library(${PROJECT_NAME} STATIC ${SOURCES} ${INCLUDES})

target_include_directories(${PROJECT_NAME} PUBLIC
        ${PROJECT_SOURCE_DIR}/include/
        FFTW_INCLUDE_DIRS)

target_link_libraries(${PROJECT_NAME}
        OpenMP::OpenMP_CXX
        hdf5::hdf5
        HighFive
        FFTW::Double)

add_subdirectory(examples)

enable_testing()
add_subdirectory(test)
