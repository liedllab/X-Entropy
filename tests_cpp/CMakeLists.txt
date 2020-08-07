include(GoogleTest)

add_executable(entropy_tests Gce_Tester.cpp  Integrator_Tester.cpp  pch.h)
add_dependencies(entropy_tests googletest)

target_link_libraries(entropy_tests PRIVATE gtest_main gtest pthread fftw3 gomp entropy_tests_lib)

gtest_discover_tests(
	entropy_tests
)
