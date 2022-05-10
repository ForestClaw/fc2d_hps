#include "gtest/gtest.h"
#include <fclaw2d_include_all.h>
#include <HPS/fc2d_hps.hpp>
#include <Structures/fc2d_hps_vector.hpp>
#include <Structures/fc2d_hps_matrix.hpp>

TEST(Matrix, init) {
	int n_rows = 3;
	int n_cols = 4;
	std::vector<int> data = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12
	};

	fc2d_hps_matrix<int> mat1();
	
	fc2d_hps_matrix<int> mat2(n_rows, n_cols, 1);

	EXPECT_EQ(mat2.size(), 12);
	EXPECT_EQ(mat2.rows, n_rows);
	EXPECT_EQ(mat2.cols, n_cols);

	fc2d_hps_matrix<int> mat3(n_rows, n_cols, data);

	EXPECT_EQ(mat3.size(), 12);
	EXPECT_EQ(mat3.rows, n_rows);
	EXPECT_EQ(mat3.cols, n_cols);
	for (int i = 0; i < n_rows; i++) {
		for (int j = 0; j < n_cols; j++) {
			EXPECT_EQ(data[j + i*n_cols], mat3[j + i*n_cols]);
		}
	}

	fc2d_hps_vector<int> vec = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12
	};
	fc2d_hps_matrix<int> mat4(n_rows, n_cols, vec);

	EXPECT_EQ(mat4.size(), 12);
	EXPECT_EQ(mat4.rows, n_rows);
	EXPECT_EQ(mat4.cols, n_cols);
	for (int i = 0; i < n_rows; i++) {
		for (int j = 0; j < n_cols; j++) {
			EXPECT_EQ(data[j + i*n_cols], mat4[j + i*n_cols]);
		}
	}
}

TEST(Matrix, index) {
	int n_rows = 3;
	int n_cols = 4;
	std::vector<int> data = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12
	};

	fc2d_hps_matrix<int> mat(n_rows, n_cols, data);

	for (int i = 0; i < n_rows; i++) {
		for (int j = 0; j < n_cols; j++) {
			EXPECT_EQ(data[j + i*n_cols], mat(i, j));
		}
	}
	for (int n = 0; n < n_rows*n_cols; n++) {
		EXPECT_EQ(data[n], mat[n]);
	}
}

TEST(Matrix, transpose) {
	int n_rows = 3;
	int n_cols = 4;
	std::vector<int> data = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12
	};
	fc2d_hps_matrix<int> mat(n_rows, n_cols, data);

	std::vector<int> data_true = {
		1, 5, 9,
		2, 6, 10,
		3, 7, 11,
		4, 8, 12
	};
	fc2d_hps_matrix<int> mat_true(n_cols, n_rows, data_true);

	fc2d_hps_matrix<int> mat_transposed = mat.transpose();

	for (int i = 0; i < n_cols; i++) {
		for (int j = 0; j < n_rows; j++) {
			EXPECT_EQ(mat_true(i,j), mat_transposed(i,j));
		}
	}
}

TEST(Matrix, from_index_set) {
	int n_rows = 4;
	int n_cols = 4;
	std::vector<int> data = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12,
		13, 14, 15, 16
	};
	fc2d_hps_matrix<int> mat(n_rows, n_cols, data);

	std::vector<int> I = {1, 3};
	std::vector<int> J = {2, 3};

	fc2d_hps_matrix<int> extracted = mat.from_index_set(I, J);
	std::vector<int> data_true = {
		7, 8,
		15, 16
	};
	fc2d_hps_matrix<int> expected(I.size(), J.size(), data_true);
	
	for (int i = 0; i < extracted.rows; i++) {
		for (int j = 0; j < extracted.cols; j++) {
			EXPECT_EQ(extracted(i,j), expected(i,j));
		}
	}

	std::vector<int> I2 = {1, 2, 3, 4, 5};
	std::vector<int> J2 = {1, 2, 3, 4, 5};
	// EXPECT_THROW(mat.from_index_set(I, J2), std::invalid_argument);

	std::vector<int> I3 = {-1};
	std::vector<int> J3 = {7};
	// EXPECT_THROW(mat.from_index_set(I3, J), std::invalid_argument);
	// EXPECT_THROW(mat.from_index_set(I, J3), std::invalid_argument);

}

TEST(Matrix, permute) {
	int n_rows = 4;
	int n_cols = 4;
	std::vector<int> data = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12,
		13, 14, 15, 16
	};
	fc2d_hps_matrix<int> mat(n_rows, n_cols, data);

	std::vector<int> pi_I1 = {0, 1, 2, 3,};
	std::vector<int> pi_J1 = {1, 0, 3, 2};

	fc2d_hps_matrix<int> permuted1 = mat.from_index_set(pi_I1, pi_J1);
	std::vector<int> data_true1 = {
		2, 1, 4, 3,
		6, 5, 8, 7,
		10, 9, 12, 11,
		14, 13, 16, 15
	};
	fc2d_hps_matrix<int> expected1(pi_I1.size(), pi_J1.size(), data_true1);

	for (int i = 0; i < permuted1.rows; i++) {
		for (int j = 0; j < permuted1.cols; j++) {
			EXPECT_EQ(permuted1(i,j), expected1(i,j));
			printf("%3i", permuted1(i,j));
		}
		printf("\n");
	}

	std::vector<int> pi_I2 = {1, 0, 3, 2};
	std::vector<int> pi_J2 = {0, 1, 2, 3};

	fc2d_hps_matrix<int> permuted2 = mat.from_index_set(pi_I2, pi_J2);
	std::vector<int> data_true2 = {
		5, 6, 7, 8,
		1, 2, 3, 4,
		13, 14, 15, 16,
		9, 10, 11, 12
	};
	fc2d_hps_matrix<int> expected2(pi_I2.size(), pi_J2.size(), data_true2);

	for (int i = 0; i < permuted2.rows; i++) {
		for (int j = 0; j < permuted2.cols; j++) {
			EXPECT_EQ(permuted2(i,j), expected2(i,j));
			printf("%3i", permuted2(i,j));
		}
		printf("\n");
	}

	std::vector<int> pi_I3 = {0, 2, 1, 3};
	std::vector<int> pi_J3 = {0, 2, 1, 3};

	fc2d_hps_matrix<int> permuted3 = mat.from_index_set(pi_I3, pi_J3);
	std::vector<int> data_true3 = {
		1, 3, 2, 4,
		9, 11, 10, 12,
		5, 7, 6, 8,
		13, 15, 14, 16
	};
	fc2d_hps_matrix<int> expected3(pi_I3.size(), pi_J3.size(), data_true3);

	for (int i = 0; i < permuted3.rows; i++) {
		for (int j = 0; j < permuted3.cols; j++) {
			EXPECT_EQ(permuted3(i,j), expected3(i,j));
			printf("%3i", permuted3(i,j));
		}
		printf("\n");
	}

}

TEST(Matrix, permute_block) {
	int n_rows1 = 3;
	int n_cols1 = 4;
	std::vector<int> data1 = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12
	};
	fc2d_hps_matrix<int> mat1(n_rows1, n_cols1, data1);

	std::vector<int> pi_I1 = {1, 0};
	std::vector<int> pi_J1 = {1, 0};
	std::vector<int> R1 = {2, 1};
	std::vector<int> C1 = {3, 1};

	fc2d_hps_matrix<int> permuted1 = mat1.block_permute(pi_I1, pi_J1, R1, C1);
	std::vector<int> data_true1 = {
		12, 9, 10, 11,
		4, 1, 2, 3,
		8, 5, 6, 7
	};
	fc2d_hps_matrix<int> expected1(n_rows1, n_cols1, data_true1);

	for (int i = 0; i < expected1.rows; i++) {
		for (int j = 0; j < expected1.cols; j++) {
			EXPECT_EQ(permuted1(i,j), expected1(i,j));
		}
	}

	int n_rows2 = 6;
	int n_cols2 = 8;
	std::vector<int> data2 = {
		1, 2, 1, 1, 2, 3, 1, 2,
		3, 4, 2, 4, 5, 6, 3, 4,
		1, 2, 1, 1, 2, 3, 1, 2,
		3, 4, 2, 4, 5, 6, 3, 4,
		5, 6, 3, 7, 8, 9, 5, 6,
		1, 2, 1, 1, 2, 3, 1, 2
	};
	fc2d_hps_matrix<int> mat2(n_rows2, n_cols2, data2);

	std::vector<int> pi_I2 = {2, 1, 0};
	std::vector<int> pi_J2 = {2, 3, 0, 1};
	std::vector<int> R2 = {2, 3, 1};
	std::vector<int> C2 = {2, 1, 3, 2};

	fc2d_hps_matrix<int> permuted2 = mat2.block_permute(pi_I2, pi_J2, R2, C2);
	std::vector<int> data_true2 = {
		1, 2, 3, 1, 2, 1, 2, 1,
		1, 2, 3, 1, 2, 1, 2, 1,
		4, 5, 6, 3, 4, 3, 4, 2,
		7, 8, 9, 5, 6, 5, 6, 3,
		1, 2, 3, 1, 2, 1, 2, 1,
		4, 5, 6, 3, 4, 3, 4, 2
	};
	fc2d_hps_matrix<int> expected2(n_rows2, n_cols2, data_true2);

	for (int i = 0; i < expected2.rows; i++) {
		for (int j = 0; j < expected2.cols; j++) {
			EXPECT_EQ(permuted2(i,j), expected2(i,j));
		}
	}

}

TEST(Matrix, extract) {
	int n_rows = 3;
	int n_cols = 4;
	std::vector<int> data = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12
	};
	fc2d_hps_matrix<int> mat(n_rows, n_cols, data);

	fc2d_hps_matrix<int> extracted = mat.extract(0, 1, 3, 2);

	std::vector<int> data_true = {
		2, 3,
		6, 7,
		10, 11
	};
	fc2d_hps_matrix<int> extracted_true(3, 2, data_true);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {
			EXPECT_EQ(extracted_true(i,j), extracted(i,j));
		}
	}
}

TEST(Matrix, extract_row) {
	int n_rows = 3;
	int n_cols = 4;
	std::vector<int> data = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12
	};
	fc2d_hps_matrix<int> mat(n_rows, n_cols, data);

	fc2d_hps_vector<int> row = mat.extract_row(1);
	fc2d_hps_vector<int> row_true = {5, 6, 7, 8};

	EXPECT_EQ(row_true, row);
}

TEST(Matrix, extract_col) {
	int n_rows = 3;
	int n_cols = 4;
	std::vector<int> data = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12
	};
	fc2d_hps_matrix<int> mat(n_rows, n_cols, data);

	fc2d_hps_vector<int> col = mat.extract_col(1);
	fc2d_hps_vector<int> col_true = {2, 6, 10};

	EXPECT_EQ(col_true, col);
}

TEST(Matrix, intract) {
	int n_rows = 3;
	int n_cols = 4;
	std::vector<int> data = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12
	};
	fc2d_hps_matrix<int> mat(n_rows, n_cols, data);

	std::vector<int> data_true = {
		1, 2, 3, 4,
		5, 60, 70, 8,
		9, 100, 110, 12
	};
	fc2d_hps_matrix<int> mat_intracted_true(n_rows, n_cols, data_true);

	std::vector<int> data_to_intract = {
		60, 70,
		100, 110
	};
	fc2d_hps_matrix<int> to_intract(2, 2, data_to_intract);

	mat.intract(1, 1, to_intract);
	
	for (int i = 0; i < n_rows; i++) {
		for (int j = 0; j < n_cols; j++) {
			EXPECT_EQ(mat_intracted_true(i,j), mat(i,j));
		}
	}
}

TEST(Matrix, intract_row) {
	int n_rows = 3;
	int n_cols = 4;
	std::vector<int> data = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12
	};
	fc2d_hps_matrix<int> mat(n_rows, n_cols, data);

	std::vector<int> data_true = {
		1, 2, 3, 4,
		50, 60, 70, 80,
		9, 10, 11, 12
	};
	fc2d_hps_matrix<int> mat_true(n_rows, n_cols, data_true);

	fc2d_hps_vector<int> data_to_intract = {50, 60, 70, 80};

	mat.intract_row(1, data_to_intract);

	for (int i = 0; i < n_rows; i++) {
		for (int j = 0; j < n_cols; j++) {
			EXPECT_EQ(mat_true(i,j), mat(i,j));
		}
	}
}

TEST(Matrix, intract_column) {
	int n_rows = 3;
	int n_cols = 4;
	std::vector<int> data = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12
	};
	fc2d_hps_matrix<int> mat(n_rows, n_cols, data);

	std::vector<int> data_true = {
		1, 20, 3, 4,
		5, 60, 7, 8,
		9, 100, 11, 12
	};
	fc2d_hps_matrix<int> mat_true(n_rows, n_cols, data_true);

	fc2d_hps_vector<int> data_to_intract = {20, 60, 100};

	mat.intract_column(1, data_to_intract);

	for (int i = 0; i < n_rows; i++) {
		for (int j = 0; j < n_cols; j++) {
			EXPECT_EQ(mat_true(i,j), mat(i,j));
		}
	}
}

TEST(Matrix, plus) {
	int n_rows = 3;
	int n_cols = 4;
	std::vector<int> data = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12
	};
	fc2d_hps_matrix<int> mat(n_rows, n_cols, data);

	fc2d_hps_matrix<int> res = mat + mat;

	std::vector<int> data_true = {
		2, 4, 6, 8,
		10, 12, 14, 16,
		18, 20, 22, 24
	};
	fc2d_hps_matrix<int> mat_true(n_rows, n_cols, data_true);

	for (int i = 0; i < n_rows; i++) {
		for (int j = 0; j < n_cols; j++) {
			EXPECT_EQ(mat_true(i,j), res(i,j));
		}
	}
}

TEST(Matrix, minus) {
	int n_rows = 3;
	int n_cols = 4;
	std::vector<int> data = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12
	};
	fc2d_hps_matrix<int> mat(n_rows, n_cols, data);

	fc2d_hps_matrix<int> res = mat - mat;

	fc2d_hps_matrix<int> mat_true(n_rows, n_cols, 0);

	for (int i = 0; i < n_rows; i++) {
		for (int j = 0; j < n_cols; j++) {
			EXPECT_EQ(mat_true(i,j), res(i,j));
		}
	}
}

TEST(Matrix, mat_vec) {
	int n_rows = 3;
	int n_cols = 4;
	std::vector<double> data = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12
	};
	fc2d_hps_matrix<double> mat(n_rows, n_cols, data);
	fc2d_hps_vector<double> vec = {1, 1, 1, 1};

	fc2d_hps_vector<double> vec_true = {10, 26, 42};

	fc2d_hps_vector<double> res = mat * vec;

	EXPECT_EQ(vec_true, res);
}

TEST(Matrix, mat_mat) {
	int n_rows = 3;
	int n_cols = 4;
	std::vector<double> data = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12
	};
	fc2d_hps_matrix<double> mat(n_rows, n_cols, data);
	fc2d_hps_matrix<double> mat_transpose = mat.transpose();

	std::vector<double> data_true = {
		30, 70, 110,
		70, 174, 278,
		110, 278, 446
	};
	fc2d_hps_matrix<double> mat_true(n_rows, n_rows, data_true);

	fc2d_hps_matrix<double> res = mat * mat_transpose;

	for (int i = 0; i < n_rows; i++) {
		for (int j = 0; j < n_rows; j++) {
			EXPECT_FLOAT_EQ(mat_true(i,j), res(i,j));
		}
	}
}

TEST(Matrix, vec_solve) {
	int n_vars = 3;
	std::vector<double> data = {
		1, 2, 3,
		4, 5, 6,
		7, 8, 10
	};
	fc2d_hps_matrix<double> A(n_vars, n_vars, data);

	fc2d_hps_vector<double> b = {1, 1, 1};

	fc2d_hps_vector<double> x = solve(A, b);

	fc2d_hps_vector<double> x_true = {-1, 1, 0};

	for (int i = 0; i < n_vars; i++) {
		EXPECT_NEAR(x_true[i], x[i], 1e-15);
	}
}

TEST(Matrix, mat_solve) {
	int n_vars = 3;
	std::vector<double> data = {
		1, 2, 3,
		4, 5, 6,
		7, 8, 10
	};
	fc2d_hps_matrix<double> A(n_vars, n_vars, data);

	int n_rhss = 4;
	std::vector<double> rhs_data = {
		1, 2, 3, 4,
		1, 2, 3, 4,
		1, 2, 3, 4,
	};
	fc2d_hps_matrix<double> B(n_vars, n_rhss, rhs_data);

	fc2d_hps_matrix<double> X = solve(A, B);

	std::vector<double> data_true = {
		-1, -2, -3, -4,
		1, 2, 3, 4,
		0, 0, 0, 0
	};
	fc2d_hps_matrix<double> X_true(n_vars, n_rhss, data_true);

	for (int i = 0; i < n_vars; i++) {
		for (int j = 0; j < n_rhss; j++) {
			EXPECT_NEAR(X_true(i,j), X(i,j), 1e-14);
		}
	}
}

TEST(Matrix, mmio) {
	int n_vars = 3;
	std::vector<double> data = {
		1, 2, 3,
		4, 5, 6,
		7, 8, 10
	};
	fc2d_hps_matrix<double> A(n_vars, n_vars, data);

	std::string filename = "test_mmio.mmio";
    std::FILE* file = fopen(filename.c_str(), "w");
	A.write_to_mmio(file, DataType::real, "%16.8e");
    fclose(file);

}