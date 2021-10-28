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
add_subdirectory(tests)

# -- external libraries
add_subdirectory(external)

# -- source files
add_subdirectory(src)