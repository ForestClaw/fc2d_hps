add_executable(svd
    main.cpp
)

target_include_directories(svd PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(svd PUBLIC
    hps
)