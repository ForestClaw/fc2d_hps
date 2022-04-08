#ifndef FC2D_HPS_QUADTREE_HPP
#define FC2D_HPS_QUADTREE_HPP
#pragma once

#include <vector>
#include <string>
#include <forestclaw2d.h>
#include <p4est.h>
#include <p4est_wrap.h>
#include <p4est_iterate.h>
#include <p4est_search.h>
#include "fc2d_hps_patch.hpp"

#define FC2D_HPS_NUMBER_CHILDREN 		4
#define FC2D_HPS_QUADTREE_LOWER_LEFT 	0
#define FC2D_HPS_QUADTREE_LOWER_RIGHT 	1
#define FC2D_HPS_QUADTREE_UPPER_LEFT 	2
#define FC2D_HPS_QUADTREE_UPPER_RIGHT 	3
#define FC2D_HPS_QUADTREE_MAX_HEIGHT    40 // Arbitrary right now, change if needed

typedef std::vector<std::vector<int>> level_array;

/**
 * DOCS
 * 
 * PATTERN: Singleton
 * 
 * NOTE: How to access nodes in the tree:
 * ```
 *	for (int l = 0; l < global_indices.size(); l++) {
 *		for (int i = 0; i < global_indices[l].size(); i++) {
 *			int mID = global_indices[l][i];
 *			int pID, cID;
 *			int p_index = parent_indices[l][i];
 *			int c_index = child_indices[l][i];
 *			if (p_index != -1) {
 *				pID = global_indices[l-1][p_index];
 *			}
 *			else {
 *				pID = -1;
 *			}
 *			if (c_index != -1) {
 *				cID = global_indices[l+1][c_index];
 *			}
 *			else {
 *				cID = -1;
 *			}
 *			data = this->data[mID];
 *		}
 *	}
 * ```
 */
template<class T>
class fc2d_hps_quadtree {

public:

	p4est_t* p4est;																	// Access to the underlying `p4est_t`
	level_array global_indices;														// Index set of global IDs, by level
	level_array parent_indices;														// Index set of parent IDs of that node in the global index set, by level
	level_array child_indices;														// Index set of child IDs of that node in the global index set, by level
	std::vector<T> data;															// Data container to be accessed via iterators, store by the global index
	std::function<T(fc2d_hps_quadtree<T>*, int, int, void*)> init_fn;				// Data initialization function
	void* user;																		// User data pointer

	typedef struct quadtree_data_wrapper {

		level_array* g;
		level_array* p;
		level_array* c;
		std::vector<T>* d;
		int global_counter;

	} quadtree_data_wrapper_t;

	typedef struct quadtree_traversal_wrapper {

		std::function<void(T&)> visit;
		T& data;

	} quadtree_traversal_wrapper_t;

	fc2d_hps_quadtree(fc2d_hps_quadtree& other) = delete;
	void operator=(const fc2d_hps_quadtree&) = delete;
	~fc2d_hps_quadtree() {
		delete instance;
	}

	static fc2d_hps_quadtree* get_instance(p4est_t* p4est, std::function<T(fc2d_hps_quadtree<T>*, int, int, void*)> init_fn, void* user) {

		if (instance == nullptr) {
			instance = new fc2d_hps_quadtree(p4est, init_fn, user);
		}
		return instance;

	}

	static fc2d_hps_quadtree* get_instance() {

		if (instance == nullptr) {
			throw std::invalid_argument("[fc2d_hps_quadtree<T>::get_instance] Instance of quadtree is not set; call with arguments first.");
		}
		return instance;

	}

	void traverse_preorder(std::function<void(T&)> visit) {

		for (auto& d : this->data) visit(d);

	}

	void traverse_postorder(std::function<void(T&)> visit) {

		this->_traverse_postorder(visit, 0, 0);

	}

	void merge(std::function<void(T&, T&, T&, T&, T&)> visit) {

		for (int l = this->global_indices.size()-1; l > 0; l--) {
			for (int i = 0; i < this->global_indices[l].size(); i++) {
				// printf("l = %i,  i = %i\n", l, i);
				if (i % 4 == 0) {
					int c0_idx = i;
					int c1_idx = i+1;
					int c2_idx = i+2;
					int c3_idx = i+3;
					// int p_idx = this->parent_indices[l][i];
					// printf("gID = %i, c0_idx = %i, c1_idx = %i, c2_idx = %i, c3_idx = %i\n", this->global_indices[l][i], c0_idx, c1_idx, c2_idx, c3_idx);

					int c0ID = this->global_indices[l][c0_idx];
					int c1ID = this->global_indices[l][c1_idx];
					int c2ID = this->global_indices[l][c2_idx];
					int c3ID = this->global_indices[l][c3_idx];
					int pID = this->parent_indices[l][i];

					// printf("merging: pID = %i, c0ID = %i, c1ID = %i, c2ID = %i, c3ID = %i\n", pID, c0ID, c1ID, c2ID, c3ID);
					visit(
						this->data[pID],
						this->data[c0ID],
						this->data[c1ID],
						this->data[c2ID],
						this->data[c3ID]
					);
				}
			}
		}

	}

	void split(std::function<void(T&, T&, T&, T&, T&)> visit) {

		for (int l = 0; l < this->global_indices.size(); l++) {
			for (int i = 0; i < this->global_indices[l].size(); i++) {

				int pID = this->global_indices[l][i];
				int c0ID = this->child_indices[l][i];

				if (c0ID != -1) {
					auto it = std::find(this->global_indices[l+1].begin(), this->global_indices[l+1].end(), c0ID);
					int c_idx = it - this->global_indices[l+1].begin();
					int c1ID = this->global_indices[l+1][c_idx+1];
					int c2ID = this->global_indices[l+1][c_idx+2];
					int c3ID = this->global_indices[l+1][c_idx+3];

					// printf("splitting: pID = %i, c0ID = %i, c1ID = %i, c2ID = %i, c3ID = %i\n", pID, c0ID, c1ID, c2ID, c3ID);
					visit(
						this->data[pID],
						this->data[c0ID],
						this->data[c1ID],
						this->data[c2ID],
						this->data[c3ID]
					);
				}

				// int p_idx = i;
				// int c0_idx = this->child_indices[l][i];
				// int c1_idx = c0_idx + 1;
				// int c2_idx = c0_idx + 2;
				// int c3_idx = c0_idx + 3;

				// if (c0_idx != -1) {
				// 	int pID = this->global_indices[l][p_idx];
				// 	int c0ID = this->global_indices[l+1][c0_idx];
				// 	int c1ID = this->global_indices[l+1][c1_idx];
				// 	int c2ID = this->global_indices[l+1][c2_idx];
				// 	int c3ID = this->global_indices[l+1][c3_idx];

				// 	visit(
				// 		this->data[pID],
				// 		this->data[c0ID],
				// 		this->data[c1ID],
				// 		this->data[c2ID],
				// 		this->data[c3ID]
				// 	);
				// }
			}
		}

	}

	friend std::ostream& operator<<(std::ostream& os, const fc2d_hps_quadtree& quadtree) {

		os << "Global Indices:" << std::endl;
		for (auto& g : quadtree.global_indices) {
			for (auto& i : g) {
				os << i << ",  ";
			}
			os << std::endl;
		}
		os << std::endl;

		os << "Parent Indices:" << std::endl;
		for (auto& p : quadtree.parent_indices) {
			for (auto& i : p) {
				os << i << ",  ";
			}
			os << std::endl;
		}
		os << std::endl;

		os << "Child Indices:" << std::endl;
		for (auto& c : quadtree.child_indices) {
			for (auto& i : c) {
				os << i << ",  ";
			}
			os << std::endl;
		}
		os << std::endl;

		return os;

	}

private:
	
	static fc2d_hps_quadtree<T>* instance; // Single instance of quadtree pointer

	fc2d_hps_quadtree(); // Delete default constructor

	// Constructor from p4est object; populates level arrays
	fc2d_hps_quadtree(p4est_t* p4est, std::function<T(fc2d_hps_quadtree<T>*, int, int, void*)> init_fn, void* user) : p4est(p4est), init_fn(init_fn), user(user) { 

		_build_level_arrays();
		_build_data();

	}

	void _build_level_arrays() {

		// Reserve memory for level arrays
		p4est_tree_t* p4est_tree = p4est_tree_array_index(p4est->trees, 0); // TODO: What if this is part of a forest of trees?
		int n_levels = p4est_tree->maxlevel + 1;
		global_indices.reserve(n_levels);
		parent_indices.reserve(n_levels);
		child_indices.reserve(n_levels);

		for (int l = 0; l < n_levels; l++) {
			global_indices.push_back(std::vector<int>{});
			parent_indices.push_back(std::vector<int>{});
			child_indices.push_back(std::vector<int>{});

			global_indices[l].reserve((std::size_t) pow(2, 2*l));
			parent_indices[l].reserve((std::size_t) pow(2, 2*l));
			child_indices[l].reserve((std::size_t) pow(2, 2*l));
		}

		// Set up user pointer
		void* save_user_pointer = p4est->user_pointer;

		// Create wrapper and store in p4est user pointer
		quadtree_data_wrapper_t wrapper;
		wrapper.g = &global_indices;
		wrapper.p = &parent_indices;
		wrapper.c = &child_indices;
		wrapper.d = &data;
		wrapper.global_counter = 0;
		p4est->user_pointer = &wrapper;

		// Call `p4est_search_reorder` to traverse tree and populate level arrays
		int skip_levels = 0;
		p4est_search_reorder(p4est, skip_levels, NULL, this->p4est_visit_pre, this->p4est_visit_post, NULL, NULL);

		// Restore original p4est user pointer
		p4est->user_pointer = save_user_pointer;

	}

	void _build_data() {

		// Count total number of nodes
		int node_counter = 0;
		for (int l = 0; l < this->global_indices.size(); l++) {
			for (int i = 0; i < this->global_indices[l].size(); i++) {
				node_counter++;
			}
		}
		this->data.resize(node_counter); // Assumes a default constructor is built for T

		// Iterate through tree via levels
		for (int l = 0; l < this->global_indices.size(); l++) {
			for (int i = 0; i < this->global_indices[l].size(); i++) {
				// Call init function and place in data array
				this->data[this->global_indices[l][i]] = init_fn(this, l, i, user);
			}
		}

	}

	void _traverse_postorder(std::function<void(T&)> visit, int level, int idx) {
		int gID = this->global_indices[level][idx];
		int cID = this->child_indices[level][idx];
		if (cID == -1) {
			visit(this->data[gID]);
			return;
		}
		else {
			auto it = std::find(this->global_indices[level+1].begin(), this->global_indices[level+1].end(), cID);
			int c_idx = it - this->global_indices[level+1].begin();
			for (int i = 0; i < 4; i++) {
				_traverse_postorder(visit, level+1, c_idx+i);
			}
			visit(this->data[gID]);
		}
	}

	void _merge(std::function<void(T&, T&, T&, T&, T&)> visit, int level, int idx) {
		int gID = this->global_indices[level][idx];
		int cID = this->child_indices[level][idx];
		if (cID == -1) {

		}
	}

	static int p4est_visit_pre(p4est_t* p4est, p4est_topidx_t which_tree, p4est_quadrant_t* quadrant, p4est_locidx_t local_num, void* point) {

		// Get access to level arrays
		quadtree_data_wrapper_t* wrapper = (quadtree_data_wrapper_t*) p4est->user_pointer;
		level_array& global_indices = *(wrapper->g);
		level_array& parent_indices = *(wrapper->p);
		level_array& child_indices = *(wrapper->c);
		std::vector<T> data = *(wrapper->d);

		// Populate global index array
		global_indices[quadrant->level].push_back(wrapper->global_counter++);

		// Populate parent index array
		int pID;
		if (quadrant->level == 0) {
			pID = -1;
		}
		else {
			pID = global_indices[quadrant->level-1][global_indices[quadrant->level-1].size() - 1];
		}
		parent_indices[quadrant->level].push_back(pID);

		return 1;
	}

	static int p4est_visit_post(p4est_t* p4est, p4est_topidx_t which_tree, p4est_quadrant_t* quadrant, p4est_locidx_t local_num, void* point) {
		
		// Get access to level arrays
		quadtree_data_wrapper_t* wrapper = (quadtree_data_wrapper_t*) p4est->user_pointer;
		level_array& global_indices = *(wrapper->g);
		level_array& parent_indices = *(wrapper->p);
		level_array& child_indices = *(wrapper->c);
		std::vector<T> data = *(wrapper->d);

		// Populate child array
		int cID;
		if (local_num < 0) {
			// Current patch is not a leaf
			cID = global_indices[quadrant->level+1][global_indices[quadrant->level+1].size()-4];
		}
		else {
			// Patch is a leaf
			cID = -1;
		}
		child_indices[quadrant->level].push_back(cID);

		return 1;
	}

};

// Set initial instance to nullptr
template<class T>
fc2d_hps_quadtree<T>* fc2d_hps_quadtree<T>::instance = nullptr;

#endif // FC2D_HPS_QUADTREE_HPP
