# Check if EllipticForest path is set
if(ELLIPTICFOREST_ROOT)
    message(STATUS ${ELLIPTICFOREST_ROOT})

    # Check if EllipticForest libraries exist
    if(EXISTS ${ELLIPTICFOREST_ROOT}/include/EllipticForest.hpp)
        set(ELLIPTICFOREST_INCLUDE ${ELLIPTICFOREST_ROOT}/include CACHE INTERNAL "EllipticForest include directory")
    else()
        message(FATAL_ERROR "Could not find EllipticForest.hpp")
    endif()

    if(EXISTS ${ELLIPTICFOREST_ROOT}/lib/libEllipticForest.a)
        set(ELLIPTICFOREST_LIB ${ELLIPTICFOREST_ROOT}/lib CACHE INTERNAL "EllipticForest library directory")
    else()
        message(FATAL_ERROR "Could not find libEllipticForest.a")
    endif()

    set(ELLIPTICFOREST_FOUND TRUE CACHE INTERNAL "Flag for if EllipticForest was found")

    add_library(fc2d_hps)
    set_target_properties(fc2d_hps PROPERTIES EXPORT_NAME FC2D_HPS)
    target_include_directories(fc2d_hps PUBLIC
        ${ELLIPTICFOREST_INCLUDE}
        ${PETSC_ROOT}/include
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
    target_link_directories(fc2d_hps PUBLIC ${ELLIPTICFOREST_LIB} ${PETSC_ROOT}/lib)
    target_link_libraries(fc2d_hps PUBLIC forestclaw clawpatch EllipticForest fishpack90 petsc)
    target_sources(fc2d_hps PRIVATE fc2d_hps_interface.cpp)

    install(FILES fc2d_hps_interface.hpp DESTINATION include)
    install(TARGETS fc2d_hps
        EXPORT ${PROJECT_NAME}-targets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

    add_library(FORESTCLAW::FC2D_HPS INTERFACE IMPORTED GLOBAL)
    target_link_libraries(FORESTCLAW::FC2D_HPS INTERFACE fc2d_hps)

else()

    message(FATAL_ERROR "HPS solver is turned on, but no `ELLIPTICFOREST_ROOT` is provided")

endif()