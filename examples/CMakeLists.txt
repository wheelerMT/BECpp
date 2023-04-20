cmake_minimum_required(VERSION 3.21..3.26)

project(2d_example)

set(CMAKE_CXX_STANDARD 20)

add_executable(${PROJECT_NAME} 2d_ground_state.cpp)

target_link_libraries(${PROJECT_NAME} BECpp)

