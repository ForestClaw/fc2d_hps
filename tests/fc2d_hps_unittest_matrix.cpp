#include "gtest/gtest.h"
#include <fclaw2d_include_all.h>
#include <fc2d_hps.h>
#include <fc2d_hps_vector.hpp>
#include <fc2d_hps_matrix.hpp>

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