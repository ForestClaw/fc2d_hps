# -- unit tests
configure_file(CMakeLists.txt.in ${CMAKE_BINARY_DIR}/src/solvers/fc2d_hps/googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" . WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/src/solvers/fc2d_hps/googletest-download )
execute_process(COMMAND ${CMAKE_COMMAND} --build . WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/src/solvers/fc2d_hps/googletest-download )
add_subdirectory(${CMAKE_BINARY_DIR}/src/solvers/fc2d_hps/googletest-src ${CMAKE_BINARY_DIR}/src/solvers/fc2d_hps/googletest-build)
add_subdirectory(tests)

# -- external libraries
add_subdirectory(external)

# -- source files
add_subdirectory(src)