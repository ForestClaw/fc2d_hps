enable_testing()
find_package(GTest REQUIRED)
include(GoogleTest)

file(
	GLOB tests "./*.cpp"
)
list(REMOVE_ITEM tests "./main.cpp")

foreach(file ${tests})
	set(name)
	get_filename_component(name ${file} NAME_WE)
	add_executable(
		"${name}_tests"
		${file}
		"./main.cpp"
	)
	target_link_libraries(
		"${name}_tests"
		FORESTCLAW::HPS
		gtest_main
	)
	gtest_add_tests(
		TARGET "${name}_tests"
	)
endforeach()

