# PETSc
# Sets the following variables:
#   PETSC_INCLUDE   : Directory with PETSc header files
#   PETSC_LIB       : Built PETSc library
#   PETSC_FOUND     : Flag for if PETSc was found
if(PETSC_PATH)

    if(EXISTS ${PETSC_PATH}/include/petsc.h)
        set(PETSC_INCLUDE ${PETSC_PATH}/include PARENT_SCOPE)
    else()
        message(FATAL_ERROR "Could not find petsc.h")
    endif()

    if(EXISTS ${PETSC_PATH}/lib/libpetsc.dylib)
        add_library(PETSC_LIB SHARED IMPORTED)
        set_target_properties(PETSC_LIB PROPERTIES IMPORTED_LOCATION ${PETSC_PATH}/libpetsc.dylib)
    elseif(EXISTS ${PETSC_PATH}/lib/libpetsc.a)
        add_library(PETSC_LIB SHARED IMPORTED)
        set_target_properties(PETSC_LIB PROPERTIES IMPORTED_LOCATION ${PETSC_PATH}/libpetsc.a)
    elseif(EXISTS ${PETSC_PATH}/lib/libpetsc.so)
        add_library(PETSC_LIB SHARED IMPORTED)
        set_target_properties(PETSC_LIB PROPERTIES IMPORTED_LOCATION ${PETSC_PATH}/libpetsc.so)
    else()
        message(FATAL_ERROR "Could not find libpetsc")
    endif()

    set(PETSC_FOUND TRUE PARENT_SCOPE)

else()

    message(FATAL_ERROR "NEED TO BUILD PETSC")

endif()

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
