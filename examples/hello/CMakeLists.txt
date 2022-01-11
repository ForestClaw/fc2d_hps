add_executable(hello 
  main.cpp
)

target_include_directories(hello PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(hello PRIVATE 
  FORESTCLAW::HPS
)

configure_file(fclaw_options.ini fclaw_options.ini)