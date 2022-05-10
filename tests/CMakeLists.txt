include(CTest)
include(FetchContent)
FetchContent_Declare(googletest URL https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

file(GLOB tests "./*.cpp")
list(REMOVE_ITEM tests "./main.cpp")

foreach(file ${tests})
	set(name)
	get_filename_component(name ${file} NAME_WE)
	add_executable("${name}_tests" ${file} "./main.cpp")
	target_link_libraries("${name}_tests" hps gtest_main fishpack90 lapack mmio)
	gtest_discover_tests("${name}_tests")
endforeach()

