cmake_minimum_required(VERSION 3.13...3.24)
project(HPS LANGUAGES C CXX Fortran)
set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

# --=== External Packages ===--
add_subdirectory(external)

# --=== ForestClaw ===--
# Check if this being included as a subdirectory
# list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/ForestClawExternalProject)

# Get installed ForestClaw package

# Set to True if you want to download ForestClaw and build from scratch
set(external_build False)
if(NOT external_build)
    # 'clawpack' here refers to the "-Dclawpack=on" cmake flag
    # This uses -DFORESTCLAW_ROOT to find ForestClaw installation
    find_package(FORESTCLAW REQUIRED CONFIG COMPONENTS clawpack mpi)
else()
    # Fall back to external project build
    include(ForestClawExternalProject)
    ForestClawExternalProject(TAG develop COMPONENTS clawpack mpi)
endif()

# --=== Source ===-- 
add_subdirectory(src)

# --=== Testing ===--
if (BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

# --=== Examples and Applications --===
add_subdirectory(examples)



