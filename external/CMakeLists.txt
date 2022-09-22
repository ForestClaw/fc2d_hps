# PETSc
# Sets the following variables:
#   PETSC_INCLUDE   : Directory with PETSc header files
#   PETSC_LIB       : Built PETSc library
#   PETSC_FOUND     : Flag for if PETSc was found
if(PETSC_PATH)

    if(EXISTS ${PETSC_PATH}/include/petsc.h)
        set(PETSC_INCLUDE ${PETSC_PATH}/include CACHE INTERNAL "PETSc include directory")
    else()
        message(FATAL_ERROR "Could not find petsc.h")
    endif()

    if(EXISTS ${PETSC_PATH}/lib/libpetsc.dylib)
        set(PETSC_LIB ${PETSC_PATH}/lib CACHE INTERNAL "PETSc library directory")
        # add_library(PETSC_LIB SHARED IMPORTED GLOBAL)
        # set_target_properties(PETSC_LIB PROPERTIES IMPORTED_LOCATION ${PETSC_PATH}/lib/libpetsc.dylib)
    elseif(EXISTS ${PETSC_PATH}/lib/libpetsc.a)
        set(PETSC_LIB ${PETSC_PATH}/lib CACHE INTERNAL "PETSc library directory")
        # add_library(PETSC_LIB SHARED IMPORTED GLOBAL)
        # set_target_properties(PETSC_LIB PROPERTIES IMPORTED_LOCATION ${PETSC_PATH}/lib/libpetsc.a)
    elseif(EXISTS ${PETSC_PATH}/lib/libpetsc.so)
        set(PETSC_LIB ${PETSC_PATH}/lib CACHE INTERNAL "PETSc library directory")
        # add_library(PETSC_LIB SHARED IMPORTED GLOBAL)
        # set_target_properties(PETSC_LIB PROPERTIES IMPORTED_LOCATION ${PETSC_PATH}/lib/libpetsc.so)
    else()
        message(FATAL_ERROR "Could not find libpetsc")
    endif()

    set(PETSC_FOUND TRUE CACHE INTERNAL "Flag for if PETSc was found")

    message(STATUS "Found external PETSc:")
    message(STATUS "  PETSC_INCLUDE = " ${PETSC_INCLUDE})
    message(STATUS "  PETSC_LIB = " ${PETSC_LIB})

else()

    message(FATAL_ERROR "NEED TO BUILD PETSC")

endif()

# SLEPc
# Sets the following variables:
#   SLEPC_INCLUDE   : Directory with SLEPc header files
#   SLEPC_LIB       : Built SLEPc library
#   SLEPC_FOUND     : Flag for if SLEPc was found
# if(SLEPC_PATH)
#     if(EXISTS ${SLEPC_PATH}/include/slepc.h)
#         set(SLEPC_INCLUDE ${SLEPC_PATH}/include CACHE INTERNAL "SLEPc include directory")
#     else()
#         message(FATAL_ERROR "Could not find slepc.h")
#     endif()

#     if(EXISTS ${SLEPC_PATH}/lib/libslepc.dylib)
#         set(SLEPC_LIB ${SLEPC_PATH}/lib CACHE INTERNAL "SLEPc library directory")
#         # add_library(SLEPC_LIB SHARED IMPORTED)
#         # set_target_properties(SLEPC_LIB PROPERTIES IMPORTED_LOCATION ${SLEPC_PATH}/libslepc.dylib)
#     elseif(EXISTS ${SLEPC_PATH}/lib/libslepc.a)
#         set(SLEPC_LIB ${SLEPC_PATH}/lib CACHE INTERNAL "SLEPc library directory")
#         # add_library(SLEPC_LIB SHARED IMPORTED)
#         # set_target_properties(SLEPC_LIB PROPERTIES IMPORTED_LOCATION ${SLEPC_PATH}/libslepc.a)
#     elseif(EXISTS ${SLEPC_PATH}/lib/libslepc.so)
#         set(SLEPC_LIB ${SLEPC_PATH}/lib CACHE INTERNAL "SLEPc library directory")
#         # add_library(SLEPC_LIB SHARED IMPORTED)
#         # set_target_properties(SLEPC_LIB PROPERTIES IMPORTED_LOCATION ${SLEPC_PATH}/libslepc.so)
#     else()
#         message(FATAL_ERROR "Could not find libslepc")
#     endif()

#     set(SLEPC_FOUND TRUE CACHE INTERNAL "Flag for if SLEPc was found")

#     message(STATUS "Found external SLEPc:")
#     message(STATUS "  SLEPC_INCLUDE = " ${SLEPC_INCLUDE})
#     message(STATUS "  SLEPC_LIB = " ${SLEPC_LIB})
# else()
#     message(FATAL_ERROR "NEED TO BUILD SLEPC")
# endif()

# FISHPACK90
if (FISHPACK_PATH)



else()

    message(STATUS "Building FISHPACK90...")
    add_subdirectory(fishpack90)

endif()

# Matrix Market
if (MMIO_PATH)

else()
    
    message(STATUS "Building MatrixMarket...")
    add_subdirectory(mmio)

endif()

add_subdirectory(matplotlibcpp)

FetchContent_GetProperties(matplotlibcpp)

if(${matplotlibcpp_POPULATED} AND ${PETSC_FOUND})
    
    set(EXTERNAL_PKGS_INCLUDE
        ${PYTHON_INCLUDE_DIR}
        ${NUMPY_INCLUDE_DIR}
        ${matplotlibcpp_SOURCE_DIR}
        ${PETSC_INCLUDE}
        CACHE INTERNAL "HydroForest external pacakges' include directories"
    )

    set(EXTERNAL_PKGS_LIB
        ${PYTHON_LIBRARY_DIR}
        ${PETSC_LIB}
        CACHE INTERNAL "HydroForest external packages' library directories"
    )

else()

    message(FATAL_ERROR "External packages not ready...")

endif()