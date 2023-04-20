cmake_minimum_required(VERSION 3.21)
set(CMAKE_CXX_STANDARD 20)

include(FetchContent)
FETCHCONTENT_DECLARE(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.13.0.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set(SOURCE_FILES test_grid.cpp test_wavefunction.cpp test_data.cpp)

add_executable(tests
        ${SOURCE_FILES}
)
target_link_libraries(tests
        BECpp
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(tests)
