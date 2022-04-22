set(include_dir ${PROJECT_SOURCE_DIR}/include)

add_library(simple_f OBJECT
  fclaw2d_clawpatch_user_exceeds_th.f90
  hsmooth.f90 
  fortran/simple_setprob.f90 
  fortran/simple_fort_rhs.f90 
  fortran/simple_compute_error.f 
  fortran/simple_qexact.f90 
  fortran/simple_apply_bc.f90 
  fortran/simple_eval_bc.f90 
  fortran/simple_output.f90
)

add_executable(simple 
  main.cpp 
  simple_user.cpp 
  simple_options.c 
  simple_diagnostics.c 
  $<TARGET_OBJECTS:simple_f>
)

target_link_libraries(simple PRIVATE
  mmio
)


target_include_directories(simple PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${include_dir}
)

target_link_libraries(simple PRIVATE 
  FORESTCLAW::FORESTCLAW
)

configure_file(fclaw_options.ini fclaw_options.ini)
configure_file(scripts/run_timing.py run_timing.py)
configure_file(scripts/hps_simple_analysis.py hps_simple_analysis.py)