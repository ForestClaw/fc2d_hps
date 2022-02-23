# -- fortran library
add_library(hps_f OBJECT
  fortran_source/fc2d_hps_tag4refinement.f90
  fortran_source/fc2d_hps_tag4coarsening.f90
  fortran_source/fc2d_hps_apply_bc_default.f90  
  fortran_source/fc2d_hps_eval_bc_default.f90
)

target_include_directories(hps_f INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)

# -- add library
add_library(hps)
set_target_properties(hps PROPERTIES EXPORT_NAME HPS)

set(HPS_sources
  HPS/fc2d_hps.cpp
)
set(Methods_sources
  Methods/fc2d_hps_setup.cpp
  Methods/fc2d_hps_merge.cpp
  Methods/fc2d_hps_build.cpp
  Methods/fc2d_hps_upwards.cpp
  Methods/fc2d_hps_split.cpp
  Methods/fc2d_hps_solve.cpp
  Methods/fc2d_hps_copy.cpp
)
set(Structures_sources
  Structures/fc2d_hps_patch.cpp
  Structures/fc2d_hps_patchgrid.cpp
  Structures/fc2d_hps_patchsolver.cpp
  Structures/fc2d_hps_poissonproblem.cpp
)
set(Util_sources
  Util/fc2d_hps_diagnostics.c
  Util/fc2d_hps_log.cpp
  Util/fc2d_hps_options.c
  Util/fc2d_hps_output_ascii.c
  Util/fc2d_hps_physical_bc.c
)

target_sources(hps PRIVATE 
  ${HPS_sources}
  ${Methods_sources}
  ${Structures_sources}
  ${Util_sources}
  $<TARGET_OBJECTS:hps_f>
)

target_include_directories(hps
  PUBLIC
  ${PROJECT_SOURCE_DIR}/include
)

target_include_directories(hps
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

target_include_directories(hps PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../external/mmio)
# target_link_directories(hps PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../external/mmio)

target_link_libraries(hps PUBLIC clawpatch blas lapack mmio fishpack90)

# -- install

# install(FILES
#   fc2d_hps.h       
#   fc2d_hps_fort.h  
#   fc2d_hps_options.h  
#   fc2d_hps_physical_bc.h  
#   fc2d_hps_solve.h
#   fc2d_hps_vector.hpp
#   fc2d_hps_matrix.hpp
#   DESTINATION include
# )

# install(TARGETS hps
#   EXPORT ${PROJECT_NAME}Targets
#   ARCHIVE DESTINATION lib
#   LIBRARY DESTINATION lib)


# -- imported target, for use from FetchContent
add_library(FORESTCLAW::HPS INTERFACE IMPORTED GLOBAL)
target_link_libraries(FORESTCLAW::HPS INTERFACE hps)
