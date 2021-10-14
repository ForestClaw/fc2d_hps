/*
  Copyright (c) 2019-2021 Carsten Burstedde, Donna Calhoun, Damyn Chipman
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef FC2D_HPS_MATRIX_HPP
#define FC2D_HPS_MATRIX_HPP

#include <vector>
#include <string>
#include <iostream>
#include "fc2d_hps_vector.hpp"

extern "C" {
    void dgemv_(char* TRANS, int* M, int* N, double* ALPHA, double* A, int* LDA, double* X, int* INCX, double* BETA, double* Y, int* INCY);
    void dgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* ALPHA, double* A, int* LDA, double* B, int* LDB, double* BETA, double* C, int* LDC);
    void dgesv_(int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
}

template<class T>
class fc2d_hps_matrix : public std::vector<T> {

public:

	// Inherit and use std::vector class
	using std::vector<T>::vector;
	std::size_t rows{};
	std::size_t cols{};

	fc2d_hps_matrix() :
		rows(0), cols(0)
			{}

	fc2d_hps_matrix(std::size_t n_rows, std::size_t n_cols) :
		std::vector<T>(n_rows * n_cols), rows(n_rows), cols(n_cols)
			{}

	fc2d_hps_matrix(std::size_t n_rows, std::size_t n_cols, T value) :
		std::vector<T>(n_rows * n_cols, value), rows(n_rows), cols(n_cols)
			{}

	fc2d_hps_matrix(std::size_t n_rows, std::size_t n_cols, std::vector<T> data) :
		std::vector<T>(data), rows(n_rows), cols(n_cols)
			{}

	// TODO: Write copy and move constructors
	// fc2d_hps_matrix(const fc2d_hps_matrix<T>& to_copy) :
	// 	std::vector<T>(to_copy.rows * to_copy.cols), rows(to_copy.rows), cols(to_copy.cols)
	// 	{
	// 		for (int i = 0; i < rows; i++) {
	// 			for (int j = 0; j < cols; j++) {
	// 				this->operator()(i, j) = to_copy.operator()(i, j);
	// 			}
	// 		}
	// 	}

	T& operator()(std::size_t i, std::size_t j) {
		return this->operator[](flatten_index(i,j));
	}

	fc2d_hps_matrix<T> transpose() {
		fc2d_hps_vector<T> data;
		data.reserve(cols*rows);

		for (int j = 0; j < this->cols; j++) {
			for (int i = 0; i < this->rows; i++) {
				data.emplace_back(this->operator()(i,j));
			}
		}

		return fc2d_hps_matrix<T>{this->cols, this->rows, std::move(data)};
	}

	fc2d_hps_matrix<T> extract(std::size_t row_index, std::size_t col_index, std::size_t row_length, std::size_t col_length) {
		if (row_index + row_length > rows) {
			throw std::invalid_argument("[fc2d_hps_matrix<T>::extract] Row size exceeds matrix size");
		}
		if (col_index + col_length > cols) {
			throw std::invalid_argument("[fc2d_hps_matrix<T>::extract] Column size exceeds matrix size");
		}

		std::vector<T> out_data;
		out_data.reserve(row_length * col_length);
		for (int i = row_index; i < row_index + row_length; i++) {
			for (int j = col_index; j < col_index + col_length; j++) {
				out_data.emplace_back(this->operator()(i,j));
			}
		}

		return {row_length, col_length, std::move(out_data)};
	}

	fc2d_hps_vector<T> extract_row(std::size_t row_index) {
		if (row_index > rows) {
			throw std::invalid_argument("[fc2d_hps_matrix<T>::extract_row] Row size exceeds matrix size");
		}

		fc2d_hps_vector<T> out_data(cols);
		for (int j = 0; j < cols; j++) {
			out_data[j] = this->operator()(row_index, j);
		}
		return out_data;
	}

	fc2d_hps_vector<T> extract_col(std::size_t col_index) {
		if (col_index > cols) {
			throw std::invalid_argument("[fc2d_hps_matrix<T>::extract_col] Column size exceeds matrix size");
		}

		fc2d_hps_vector<T> out_data(rows);
		for (int i = 0; i < rows; i++) {
			out_data[i] = this->operator()(i, col_index);
		}
		return out_data;
	}

	void intract(std::size_t row_index, std::size_t col_index, fc2d_hps_matrix<T>& mat) {
		if (row_index + mat.rows > this->rows) {
			throw std::invalid_argument("[fc2d_hps_matrix<T>::intract] Row size exceeds matrix size");
		}
		if (col_index + mat.cols > this->cols) {
			throw std::invalid_argument("[fc2d_hps_matrix<T>::intract] Column size exceeds matrix size");
		}

		for (int i = row_index; i < row_index + mat.rows; i++) {
			for (int j = col_index; j < col_index + mat.cols; j++) {
				this->operator()(i,j) = mat(i - row_index, j - col_index);
			}
		}
	}

	void intract_row(std::size_t row_index, fc2d_hps_vector<T>& vec) {
		if (row_index > this->rows) {
			throw std::invalid_argument("[fc2d_hps_matrix<T>::intract_row] Row index exceeds matrix size");
		}
		if (vec.size() != this->cols) {
			throw std::invalid_argument("[fc2d_hps_matrix<T>::intract_row] Size of vector to intract does not match matrix size");
		}

		for (int j = 0; j < this->cols; j++) {
			this->operator()(row_index, j) = vec[j];
		}
	}

	void intract_column(std::size_t col_index, fc2d_hps_vector<T>& vec) {
		if (col_index > this->cols) {
			throw std::invalid_argument("[fc2d_hps_matrix<T>::intract_column] Column index exceeds matrix size");
		}
		if (vec.size() != this->rows) {
			throw std::invalid_argument("[fc2d_hps_matrix<T>::intract_column] Size of vector to intract does not match matrix size");
		}

		for (int i = 0; i < this->rows; i++) {
			this->operator()(i, col_index) = vec[i];
		}
	}

private:

	std::size_t flatten_index(std::size_t i, std::size_t j) {
		assert(i < rows);
		assert(j < cols);
		std::size_t index = i*cols + j;
		assert(index < rows*cols);
		return index;
	}

};

template<class T>
fc2d_hps_matrix<T> operator+(fc2d_hps_matrix<T>& A, fc2d_hps_matrix<T>& B) {
	if (A.rows != B.rows && A.cols != B.cols) {
		throw std::invalid_argument("[fc2d_hps_matrix<T> operator+] Invalid matrix sizes; `LHS` and `RHS` are not the same size");
	}

	fc2d_hps_matrix<T> R(A.rows, A.cols, 0);
	for (int i = 0; i < R.rows; i++) {
		for (int j = 0; j < R.cols; j++) {
			R(i,j) = A(i,j) + B(i,j);
		}
	}

	return R;
}

template<class T>
fc2d_hps_matrix<T> operator-(fc2d_hps_matrix<T>& A, fc2d_hps_matrix<T>& B) {
	if (A.rows != B.rows && A.cols != B.cols) {
		throw std::invalid_argument("[fc2d_hps_matrix<T> operator+] Invalid matrix sizes; `LHS` and `RHS` are not the same size");
	}

	fc2d_hps_matrix<T> R(A.rows, A.cols, 0);
	for (int i = 0; i < R.rows; i++) {
		for (int j = 0; j < R.cols; j++) {
			R(i,j) = A(i,j) - B(i,j);
		}
	}

	return R;
}

template<class T>
fc2d_hps_vector<T> operator*(fc2d_hps_matrix<T>& A, fc2d_hps_vector<T>& x) {
	if (A.cols != x.size()) {
		throw std::invalid_argument("[fc2d_hps_matrix<T> operator*] Invalid matrix and vector dimensions. Matrix `A` must have same number of columns as entries in `x`");
	}

	fc2d_hps_vector<T> b(A.rows);

	// Setup BLAS call
    // void dgemv_(int* TRANS, int* M, int* N, double* ALPHA, double* A, int* LDA, double* X, int* INCX, double* BETA, double* Y, int* INCY);
    char TRANS_ = 'C';
    int M_ = A.cols;
    int N_ = A.rows;
    double ALPHA_ = 1.0;
    double* A_ = A.data();
    int LDA_ = A.cols;
    double* X_ = x.data();
    int INCX_ = 1;
    double BETA_ = 0.0;
    double* Y_ = b.data();
    int INCY_ = 1;
    dgemv_(&TRANS_, &M_, &N_, &ALPHA_, A_, &LDA_, X_, &INCX_, &BETA_, Y_, &INCY_);

	return b;
}

template<class T>
fc2d_hps_matrix<T> operator*(fc2d_hps_matrix<T>& A, fc2d_hps_matrix<T>& B) {
	if (A.cols != B.rows) {
		throw std::invalid_argument("[fc2d_hps_matrix<T> operator*] Invalid matrix dimensions. Matrix `A` must have same number of columns as rows of `B`");
	}

	fc2d_hps_matrix<T> CT(B.cols, A.rows, 0);

    // Setup call
    // void dgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* ALPHA, double* A, int* LDA, double* B, int* LDB, double* BETA, dobule* C, int* LDC);
    char TRANSA_ = 'C';
    char TRANSB_ = 'C';
    int M_ = A.rows;
    int N_ = B.cols;
    int K_ = A.cols; // B.cols();
    double ALPHA_ = 1.0;
    double* A_ = A.data();
    int LDA_ = K_;
    double* B_ = B.data();
    int LDB_ = N_;
    double BETA_ = 0.0;
    double* C_ = CT.data();
    int LDC_ = M_;
    dgemm_(&TRANSA_, &TRANSB_, &M_, &N_, &K_, &ALPHA_, A_, &LDA_, B_, &LDB_, &BETA_, C_, &LDC_);

	return CT.transpose();
}

template<class T>
fc2d_hps_vector<T> solve(fc2d_hps_matrix<T>& A, fc2d_hps_vector<T>& b) {
    if (A.rows != A.cols) {
        throw std::invalid_argument("[fc2d_hps_matrix<T> solve] Matrix must be square (rows != cols).");
    }

    // Setup output vector and utility variables
    fc2d_hps_vector<T> x(b);
    fc2d_hps_matrix<T> AT = A.transpose();
    fc2d_hps_vector<int> p(b.size());

    // Setup call
    int N_ = AT.rows;
    int NRHS_ = 1;
    double* A_ = AT.data();
    int LDA_ = AT.rows;
    int* IPIV_ = p.data();
    double* B_ = x.data();
    int LDB_ = b.size();
    int INFO_;
    dgesv_(&N_, &NRHS_, A_, &LDA_, IPIV_, B_, &LDB_, &INFO_);

    // Check output
    if (INFO_ != 0) {
        std::cerr << "[fc2d_hps_matrix<T> solve] Fortran call to `dgesv_` returned non-zero flag of: " << INFO_ << std::endl;
    }

    return x;
}

template<class T>
fc2d_hps_matrix<T> solve(fc2d_hps_matrix<T>& A, fc2d_hps_matrix<T>& B) {

    // Check inputs
    if (A.rows != A.cols) {
        throw std::invalid_argument("[fc2d_hps_matrix<T> solve] Matrix must be square (rows != cols).");
    }

    // Setup output matrix and utility variables
    fc2d_hps_matrix<T> X = B.transpose();
    fc2d_hps_matrix<T> AT = A.transpose();
    fc2d_hps_vector<int> p(A.rows);

    // Setup call
    int N_ = AT.rows;
    int NRHS_ = B.cols;
    double* A_ = AT.data();
    int LDA_ = AT.rows;
    int* IPIV_ = p.data();
    double* B_ = X.data();
    int LDB_ = B.rows;
    int INFO_;
    dgesv_(&N_, &NRHS_, A_, &LDA_, IPIV_, B_, &LDB_, &INFO_);

    // Check output
    if (INFO_ != 0) {
        std::cerr << "[fc2d_hps_matrix<T> solve] Fortran call to `dgesv_` returned non-zero flag of: " << INFO_ << std::endl;
    }

    return X.transpose();
}

#endif // FC2D_HPS_MATRIX_HPP