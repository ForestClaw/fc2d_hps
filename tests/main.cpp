#include "gtest/gtest.h"
#include <string>

int main(int argc, char* argv[]) {
	::testing::InitGoogleTest(&argc, argv);
	testing::internal::CaptureStdout();
	std::string output = testing::internal::GetCapturedStdout();
	std::cout << "Running tests..." << std::endl;
	std::cout << output << std::endl;
	return RUN_ALL_TESTS();
}