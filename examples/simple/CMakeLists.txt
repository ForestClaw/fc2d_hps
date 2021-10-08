add_library(laplace_f OBJECT
  fclaw2d_clawpatch_user_exceeds_th.f90
  hsmooth.f90 
  fortran/laplace_setprob.f90 
  fortran/laplace_fort_rhs.f90 
  fortran/laplace_compute_error.f 
  fortran/laplace_qexact.f90 
  fortran/laplace_apply_bc.f90 
  fortran/laplace_eval_bc.f90 
  fortran/laplace_output.f90
)

add_executable(laplace 
  laplace.cpp 
  laplace_user.cpp 
  laplace_options.c 
  laplace_diagnostics.c 
  $<TARGET_OBJECTS:laplace_f>
)

target_include_directories(laplace PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(laplace PRIVATE 
  FORESTCLAW::HPS
)