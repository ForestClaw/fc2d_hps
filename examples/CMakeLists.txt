# --------------------------------------------
# List of examples
#
# --------------------------------------------

cmake_minimum_required(VERSION 3.13...3.20)

project(HPSApps LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 17)

# -- check if this being included as a subdirectory
if(NOT TARGET FORESTCLAW::FORESTCLAW)

    set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ON)
    find_package(MPI)

    # --- get installed ForestClaw package
    find_package(FORESTCLAW REQUIRED)

endif()

# -- set environment string used for regression tests
set(FCLAW_TEST_ENVIRONMENT "FCLAW_APPLICATIONS_BUILD_DIR=${PROJECT_BINARY_DIR};FCLAW_APPLICATIONS_SRC_DIR=${PROJECT_SOURCE_DIR}")
if(TARGET MPI::MPI_C)
    list(APPEND FCLAW_TEST_ENVIRONMENT "FCLAW_MPIRUN=${MPIEXEC_EXECUTABLE}")
    list(APPEND FCLAW_TEST_ENVIRONMENT "FCLAW_MPI_TEST_FLAGS=${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}")
endif()



# ----------------------------------
# HPS library and examples
# ----------------------------------
if(TARGET FORESTCLAW::HPS)

    # Hello world program
    add_subdirectory(hello)

    # Simple
    add_subdirectory(simple)

    # Poisson's equation (timing and error testing)
    add_subdirectory(poisson)

endif()

