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

target_sources(hps PRIVATE 
  fc2d_hps.cpp
  fc2d_hps_diagnostics.c  
  fc2d_hps_merge.cpp  
  fc2d_hps_options.c  
  fc2d_hps_output_ascii.c  
  fc2d_hps_patch.cpp
  fc2d_hps_patchgrid.cpp
  fc2d_hps_patchsolver.cpp
  fc2d_hps_physical_bc.c  
  fc2d_hps_poissonproblem.cpp
  fc2d_hps_solve.cpp
  fc2d_hps_split.cpp
  $<TARGET_OBJECTS:hps_f>
)

target_link_libraries(hps PUBLIC clawpatch blas lapack fishpack90)

target_include_directories(hps
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

# -- test bench executable (delete later)
add_executable(
  test_bench test_bench.cpp
)
target_link_libraries(
  test_bench hps
)

# -- install

install(FILES
  fc2d_hps.h       
  fc2d_hps_fort.h  
  fc2d_hps_options.h  
  fc2d_hps_physical_bc.h  
  fc2d_hps_solve.h
  fc2d_hps_vector.hpp
  fc2d_hps_matrix.hpp
  DESTINATION include
)

install(TARGETS hps
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)


# -- imported target, for use from FetchContent
add_library(FORESTCLAW::HPS INTERFACE IMPORTED GLOBAL)
target_link_libraries(FORESTCLAW::HPS INTERFACE hps)