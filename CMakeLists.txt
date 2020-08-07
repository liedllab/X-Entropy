cmake_minimum_required(VERSION 3.10)
project(entropy_tests)

set(CMAKE_CXX_STANDARD 17)

include(ExternalProject)

set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR})

ExternalProject_Add(googletest
	GIT_REPOSITORY https://github.com/google/googletest.git
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
)

include_directories(src ${EXTERNAL_INSTALL_LOCATION}/include)
link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)

add_subdirectory(src)

ENABLE_TESTING()
add_subdirectory(tests_cpp)
