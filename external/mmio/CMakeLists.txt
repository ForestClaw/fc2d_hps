add_library(mmio STATIC mmio.c)
target_compile_options(mmio PUBLIC "-w")
install(TARGETS mmio EXPORT ${PROJECT_NAME}Targets)

target_include_directories(mmio PUBLIC 
     ${CMAKE_CURRENT_BINARY_DIR}
)
