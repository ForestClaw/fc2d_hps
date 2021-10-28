cmake_minimum_required(VERSION 3.13...3.20)
project(fc2d_hps)

set(CMAKE_CXX_STANDARD 11)

# -- unit tests
include(FetchContent)
FetchContent_Declare(
	googletest
	URL https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip
)
FetchContent_MakeAvailable(googletest)
# configure_file(CMakeLists.txt.in ${CMAKE_BINARY_DIR}/src/solvers/fc2d_hps/googletest-download/CMakeLists.txt)
# execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" . WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/src/solvers/fc2d_hps/googletest-download )
# execute_process(COMMAND ${CMAKE_COMMAND} --build . WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/src/solvers/fc2d_hps/googletest-download )
# add_subdirectory(${CMAKE_BINARY_DIR}/src/solvers/fc2d_hps/googletest-src ${CMAKE_BINARY_DIR}/src/solvers/fc2d_hps/googletest-build)
add_subdirectory(tests)

# -- external libraries
add_subdirectory(external)

# -- source files
add_subdirectory(src)