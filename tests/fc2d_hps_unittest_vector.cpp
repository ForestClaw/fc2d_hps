#include "gtest/gtest.h"
#include <fclaw2d_include_all.h>
#include <HPS/fc2d_hps.hpp>
#include <Structures/fc2d_hps_vector.hpp>

TEST(Vector, init) {
	fc2d_hps_vector<int> vec = {0, 1, 2, 3};
	ASSERT_TRUE(true);
}

TEST(Vector, size) {
	fc2d_hps_vector<int> vec_int = {0, 1, 2, 3};
	fc2d_hps_vector<double> vec_double = {0., 1., 2., 3.};
	std::size_t true_size = 4;

	ASSERT_EQ(true_size, vec_int.size());
	ASSERT_EQ(true_size, vec_double.size());
}

TEST(Vector, index) {
	fc2d_hps_vector<int> vec_int = {0, 1, 2, 3};
	fc2d_hps_vector<double> vec_double = {0., 1., 2., 3.};

	// @TODO: Uncomment after fixing fc2d_hps_vector.hpp TODO for index operator.
	// EXPECT_THROW(vec_int[-1], std::out_of_range);
	// EXPECT_THROW(vec_int[4], std::out_of_range);
	// EXPECT_THROW(vec_double[-1], std::out_of_range);
	// EXPECT_THROW(vec_double[4], std::out_of_range);

	for (int i = 0; i < 4; i++) {
		ASSERT_EQ(i, vec_int[i]);
		ASSERT_EQ((double) i, vec_double[i]);
	}
}

TEST(Vector, extract) {
	fc2d_hps_vector<int> vec_int = {0, 1, 2, 3};
	fc2d_hps_vector<double> vec_double = {0., 1., 2., 3.};

	EXPECT_THROW(vec_int.extract(2,3), std::invalid_argument);
	EXPECT_THROW(vec_double.extract(2,3), std::invalid_argument);

	fc2d_hps_vector<int> vec_int_extracted = vec_int.extract(1, 2);
	fc2d_hps_vector<double> vec_double_extracted = vec_double.extract(1, 2);

	fc2d_hps_vector<int> vec_int_true_extracted = {1, 2};
	fc2d_hps_vector<double> vec_double_true_extracted = {1., 2.};

	EXPECT_EQ(vec_int_true_extracted, vec_int_extracted);
	EXPECT_EQ(vec_double_true_extracted, vec_double_extracted);
}

TEST(Vector, intract) {
	fc2d_hps_vector<int> vec_int = {0, 1, 2, 3};
	fc2d_hps_vector<double> vec_double = {0., 1., 2., 3.};

	fc2d_hps_vector<int> vec_int_to_intract = {10, 20, 30};
	fc2d_hps_vector<double> vec_double_to_intract = {10., 20., 30.};

	vec_int.intract(1, vec_int_to_intract);
	vec_double.intract(1, vec_double_to_intract);

	fc2d_hps_vector<int> vec_int_true_intracted = {0, 10, 20, 30};
	fc2d_hps_vector<double> vec_double_true_intracted = {0., 10., 20., 30.};

	EXPECT_THROW(vec_int.intract(2, vec_int_to_intract), std::invalid_argument);
	EXPECT_THROW(vec_double.intract(2, vec_double_to_intract), std::invalid_argument);

	EXPECT_EQ(vec_int_true_intracted, vec_int);
	EXPECT_EQ(vec_double_true_intracted, vec_double);
}

TEST(Vector, from_index_set) {
	fc2d_hps_vector<int> vec = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	std::vector<int> I = {8, 6, 4, 2, 0};
	fc2d_hps_vector<int> vec_test = vec.from_index_set(I);
	fc2d_hps_vector<int> vec_expected = {8, 6, 4, 2, 0};
	for (int i = 0; i < vec_test.size(); i++) {
		EXPECT_EQ(vec_test[i], vec_expected[i]);
	}
}

TEST(Vector, plus) {
	fc2d_hps_vector<int> vec_int = {0, 1, 2, 3};
	fc2d_hps_vector<int> vec_int2 = {0, 1, 2};
	fc2d_hps_vector<double> vec_double = {0., 1., 2., 3.};
	fc2d_hps_vector<double> vec_double2 = {0., 1., 2.};

	EXPECT_THROW(
		{vec_int + vec_int2;},
		std::invalid_argument
	);
	EXPECT_THROW(
		{vec_double + vec_double2;},
		std::invalid_argument
	);

	fc2d_hps_vector<int> vec_int_added = vec_int + vec_int;
	fc2d_hps_vector<double> vec_double_added = vec_double + vec_double;

	fc2d_hps_vector<int> vec_int_true_added = {0, 2, 4, 6};
	fc2d_hps_vector<double> vec_double_true_added = {0., 2., 4., 6.};

	EXPECT_EQ(vec_int_true_added, vec_int_added);
	EXPECT_EQ(vec_double_true_added, vec_double_added);

	vec_int_added = vec_int + 10;
	vec_double_added = vec_double + 10.;

	vec_int_true_added = {10, 11, 12, 13};
	vec_double_true_added = {10., 11., 12., 13.};

	EXPECT_EQ(vec_int_true_added, vec_int_added);
	EXPECT_EQ(vec_double_true_added, vec_double_added);

	vec_int_added = 10 + vec_int;
	vec_double_added = 10. + vec_double;

	EXPECT_EQ(vec_int_true_added, vec_int_added);
	EXPECT_EQ(vec_double_true_added, vec_double_added);
}

TEST(Vector, minus) {
	fc2d_hps_vector<int> vec_int = {0, 1, 2, 3};
	fc2d_hps_vector<int> vec_int2 = {0, 1, 2};
	fc2d_hps_vector<double> vec_double = {0., 1., 2., 3.};
	fc2d_hps_vector<double> vec_double2 = {0., 1., 2.};

	EXPECT_THROW(
		{vec_int - vec_int2;},
		std::invalid_argument
	);
	EXPECT_THROW(
		{vec_double - vec_double2;},
		std::invalid_argument
	);

	fc2d_hps_vector<int> vec_int_subtracted = vec_int - vec_int;
	fc2d_hps_vector<double> vec_double_subtracted = vec_double - vec_double;

	fc2d_hps_vector<int> vec_int_true_subtracted = {0, 0, 0, 0};
	fc2d_hps_vector<double> vec_double_true_subtracted = {0., 0., 0., 0.};

	EXPECT_EQ(vec_int_true_subtracted, vec_int_subtracted);
	EXPECT_EQ(vec_double_true_subtracted, vec_double_subtracted);

	vec_int_subtracted = vec_int - 10;
	vec_double_subtracted = vec_double - 10.;

	vec_int_true_subtracted = {-10, -9, -8, -7};
	vec_double_true_subtracted = {-10., -9., -8., -7.};

	EXPECT_EQ(vec_int_true_subtracted, vec_int_subtracted);
	EXPECT_EQ(vec_double_true_subtracted, vec_double_subtracted);

	// @TODO: Fix function in fc2d_hps_vector.
	// vec_int_subtracted = 10 - vec_int;
	// vec_double_subtracted = 10. - vec_double;

	// vec_int_true_subtracted = {10, 9, 8, 7};
	// vec_double_true_subtracted = {10., 8., 8., 7.};

	// for (int i = 0; i < 4; i++) {
	// 	EXPECT_FLOAT_EQ(vec_int_true_subtracted[i], vec_int_subtracted[i]);
	// 	EXPECT_FLOAT_EQ(vec_double_true_subtracted[i], vec_double_subtracted[i]);
	// }
}

TEST(Vector, product) {
	fc2d_hps_vector<int> vec_int = {0, 1, 2, 3};
	fc2d_hps_vector<int> vec_int2 = {0, 1, 2};
	fc2d_hps_vector<double> vec_double = {0., 1., 2., 3.};
	fc2d_hps_vector<double> vec_double2 = {0., 1., 2.};

	EXPECT_THROW(
		{vec_int * vec_int2;},
		std::invalid_argument
	);
	EXPECT_THROW(
		{vec_double * vec_double2;},
		std::invalid_argument
	);

	fc2d_hps_vector<int> vec_int_multiplied = vec_int * vec_int;
	fc2d_hps_vector<double> vec_double_multiplied = vec_double * vec_double;

	fc2d_hps_vector<int> vec_int_true_multiplied = {0, 1, 4, 9};
	fc2d_hps_vector<double> vec_double_true_multiplied = {0., 1., 4., 9.};

	EXPECT_EQ(vec_int_true_multiplied, vec_int_multiplied);
	EXPECT_EQ(vec_double_true_multiplied, vec_double_multiplied);

	vec_int_multiplied = vec_int * 2;
	vec_double_multiplied = vec_double * 2.;
	
	vec_int_true_multiplied = {0, 2, 4, 6};
	vec_double_true_multiplied = {0., 2., 4., 6.};

	EXPECT_EQ(vec_int_true_multiplied, vec_int_multiplied);
	EXPECT_EQ(vec_double_true_multiplied, vec_double_multiplied);

	vec_int_multiplied = 2 * vec_int;
	vec_double_multiplied = 2. * vec_double;
	
	vec_int_true_multiplied = {0, 2, 4, 6};
	vec_double_true_multiplied = {0., 2., 4., 6.};

	EXPECT_EQ(vec_int_true_multiplied, vec_int_multiplied);
	EXPECT_EQ(vec_double_true_multiplied, vec_double_multiplied);

}
