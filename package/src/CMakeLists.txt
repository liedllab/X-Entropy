set(CMAKE_REQUIRED_LIBRARIES "fftw3" "omp")

add_library(entropy_tests_lib Exceptions.h  Integrators.cpp  Integrators.h  kde.cpp  kde.h)

target_link_libraries(entropy_tests_lib PRIVATE fftw3 gomp)
