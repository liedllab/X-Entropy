cmake_minimum_required(VERSION 3.10)
project(entropy_tests)

set(CMAKE_CXX_STANDARD 17)


include_directories(src)

add_subdirectory(src)

ENABLE_TESTING()
add_subdirectory(tests_cpp)
