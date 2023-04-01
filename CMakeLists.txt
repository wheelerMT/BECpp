cmake_minimum_required(VERSION 3.21)
project(BECpp)

set(CMAKE_CXX_STANDARD 20)

set(SOURCES src/grid.cpp)
set(INCLUDES include/constants.h include/grid.h)

find_package(OpenMP REQUIRED)

add_library(${PROJECT_NAME} STATIC ${SOURCES} ${INCLUDES})

target_include_directories(${PROJECT_NAME} PRIVATE
        ${PROJECT_SOURCE_DIR}/include/)
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)

enable_testing()
add_subdirectory(test)
