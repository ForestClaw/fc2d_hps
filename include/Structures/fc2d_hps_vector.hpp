#ifndef FC2D_HPS_VECTOR_HPP
#define FC2D_HPS_VECTOR_HPP

#include <vector>
#include <string>
#include <stdexcept>


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

	fc2d_hps_vector<T> from_index_set(std::vector<int> I) {
		if (I.size() > this->size()) {
			throw std::invalid_argument("[fc2d_hps_vector<T>::from_index_set] Size of index `I` is greater than size of vector.");
		}

		fc2d_hps_vector<T> output(I.size());
		for (std::size_t i = 0; i < I.size(); i++) {
			if (I[i] > this->size() || I[i] < 0) {
				throw std::invalid_argument("[fc2d_hps_vector<T>::from_index_set] Index in `I` is out of range.");
			}
			output[i] = this->operator[](I[i]);
		}
		return output;
	}

	fc2d_hps_vector<T> block_permute(std::vector<int> I, std::vector<int> S) {
		/**
		 * std::vector<int> I : Index set of permuted row indices for each block
		 * std::vector<int> S : Vector containing size of each block
		 */

		// Error checks
		std::size_t size_check = 0;
		for (auto& s : S) size_check += s;
		if (size_check != this->size()) {
			throw std::invalid_argument("[fc2d_hps_vector<T>::block_permute] Sizes in `S` do not add up to size of `this`.");
		}

		std::vector<int> S_global(this->size());

		// Build global index set
		std::size_t I_counter = 0;
		for (auto& i : I) {
			// Get starting index for i-th block
			std::size_t s = 0;
			for (std::size_t ii = 0; ii < i; ii++) s += S[ii];

			// Create increasing index set from starting index to ending index of i-th block and put into S_global
			for (std::size_t iii = s; iii < (s + S[i]); iii++) {
				S_global[I_counter++] = iii;
			}
		}

		return this->from_index_set(S_global);
	}

};

template<class T>
fc2d_hps_vector<T> operator+(fc2d_hps_vector<T> & lhs, fc2d_hps_vector<T> & rhs) {

	if (lhs.size() != rhs.size()) {
		throw std::invalid_argument("[fc2d_hps_vector<T>::operator+] Size mismatch. Size of `lhs` and `rhs` are not the same");
	}

	fc2d_hps_vector<T> res = lhs;
	for (int i = 0; i < res.size(); i++) {
		res[i] += rhs[i];
	}
	return res;

}

template<class T>
fc2d_hps_vector<T> operator+(fc2d_hps_vector<T> & lhs, T rhs) {

	fc2d_hps_vector<T> res = lhs;
	for (int i = 0; i < res.size(); i++) {
		res[i] += rhs;
	}
	return res;

}

template<class T>
fc2d_hps_vector<T> operator+(T lhs, fc2d_hps_vector<T> & rhs) {

	fc2d_hps_vector<T> res = rhs;
	for (int i = 0; i < res.size(); i++) {
		res[i] += lhs;
	}
	return res;

}

template<class T>
fc2d_hps_vector<T> operator-(fc2d_hps_vector<T> & lhs, fc2d_hps_vector<T> & rhs) {

	if (lhs.size() != rhs.size()) {
		throw std::invalid_argument("[fc2d_hps_vector<T>::operator+] Size mismatch. Size of `lhs` and `rhs` are not the same");
	}

	fc2d_hps_vector<T> res = lhs;
	for (int i = 0; i < res.size(); i++) {
		res[i] -= rhs[i];
	}
	return res;

}

template<class T>
fc2d_hps_vector<T> operator-(fc2d_hps_vector<T> & lhs, T rhs) {

	fc2d_hps_vector<T> res = lhs;
	for (int i = 0; i < res.size(); i++) {
		res[i] -= rhs;
	}
	return res;

}

template<class T>
fc2d_hps_vector<T> operator*(fc2d_hps_vector<T> & lhs, fc2d_hps_vector<T> & rhs) {

	if (lhs.size() != rhs.size()) {
		throw std::invalid_argument("[fc2d_hps_vector<T>::operator+] Size mismatch. Size of `lhs` and `rhs` are not the same");
	}

	fc2d_hps_vector<T> res = lhs;
	for (int i = 0; i < res.size(); i++) {
		res[i] *= rhs[i];
	}
	return res;

}

template<class T>
fc2d_hps_vector<T> operator*(fc2d_hps_vector<T> & lhs, T rhs) {

	fc2d_hps_vector<T> res = lhs;
	for (int i = 0; i < res.size(); i++) {
		res[i] *= rhs;
	}
	return res;

}

template<class T>
fc2d_hps_vector<T> operator*(T lhs, fc2d_hps_vector<T> & rhs) {

	fc2d_hps_vector<T> res = rhs;
	for (int i = 0; i < res.size(); i++) {
		res[i] *= lhs;
	}
	return res;

}

#endif // FC2D_HPS_VECTOR_HPP