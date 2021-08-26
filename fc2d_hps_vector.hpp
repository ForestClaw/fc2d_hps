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

#ifndef FC2D_HPS_VECTOR_HPP
#define FC2D_HPS_VECTOR_HPP

#include <vector>

template<class T>
class fc2d_hps_vector : public std::vector<T> {

public:

	// Inherit and use std::vector class
	using std::vector<T>::vector;

	/**
	 * @brief Returns size of vector.
	 * 
	 * @return std::size_t Size of vector
	 */
	std::size_t size() const {
		return this->std::vector<T>::size();
	}

	/**
	 * @brief Index operator; returns object at location `index`.
	 * 
	 * @param index Index of object
	 * @return T& Object at index
	 */
	T& operator[](std::size_t index) {
		return this->std::vector<T>::operator[](index);
	}

	/**
	 * @brief Extracts and returns a subset of the vector.
	 * 
	 * @param start_index Starting index of subset to extract
	 * @param length Length of subset to extract
	 * @return fc2d_hps_vector<T> Extracted vector as new fc2d_hps_vector instance
	 */
	fc2d_hps_vector<T> extract(int start_index, int length) const {

		if (start_index + length > this->size()) {
			throw std::invalid_argument("[fc2d_hps_vector<T>::extract] Index mismatch. `start_index` + `length` is greater than size of vector.");
		}
			
		fc2d_hps_vector<T> extracted_vector;
		extracted_vector.reserve(length);
		for (int i = start_index; i < start_index + length; i++) {
			extracted_vector.emplace_back(this->operator[](i));
		}
		return extracted_vector;

	}

	/**
	 * @brief Intracts (i.e., inserts, opposite of extract) `vec` into this vector starting at `start_index`.
	 * 
	 * @param start_index Starting index to insert `vec` into this vector
	 * @param vec Vector to intract into this vector
	 * @return None
	 */
	void intract(int start_index, const fc2d_hps_vector<T>& vec) {

		if (start_index + vec.size() > this->size()) {
			throw std::invalid_argument("[fc2d_hps_vector<T>::intract] Index mismatch. `start_index` + `vec.size()` is greater than size of host vector");
		}

		for (int i = start_index; i < start_index + vec.size(); i++) {
			this->operator[](i) = vec[i - start_index];
		}

	}

private:



};

#endif // FC2D_HPS_VECTOR_HPP