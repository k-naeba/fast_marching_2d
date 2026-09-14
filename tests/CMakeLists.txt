include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.14.0
)
FetchContent_MakeAvailable(googletest)

add_executable(test_fast_marching_2d test_fast_marching_2d.cpp)
target_link_libraries(test_fast_marching_2d
    PRIVATE
        fast_marching_2d
        GTest::gtest_main
)
add_test(NAME test_fast_marching_2d COMMAND test_fast_marching_2d)
