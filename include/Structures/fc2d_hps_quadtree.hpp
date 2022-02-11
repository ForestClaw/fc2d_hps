#ifndef FC2D_HPS_QUADTREE_HPP
#define FC2D_HPS_QUADTREE_HPP

#include <vector>
#include <forestclaw2d.h>
#include <p4est.h>
#include <p4est_wrap.h>
#include <p4est_iterate.h>
#include "fc2d_hps_patch.hpp"

#define FC2D_HPS_NUMBER_CHILDREN 		4
#define FC2D_HPS_QUADTREE_LOWER_LEFT 	0
#define FC2D_HPS_QUADTREE_LOWER_RIGHT 	1
#define FC2D_HPS_QUADTREE_UPPER_LEFT 	2
#define FC2D_HPS_QUADTREE_UPPER_RIGHT 	3
#define FC2D_HPS_QUADTREE_MAX_HEIGHT    40 // Arbitrary right now, change if needed

template<class T>
class fc2d_hps_quadtree {

public:

	struct fc2d_hps_quadnode {
		T data;
		std::size_t level;
		fc2d_hps_quadnode* children[FC2D_HPS_NUMBER_CHILDREN];
		fc2d_hps_quadnode* parent;
		fc2d_hps_quadnode() : level(0) {}
		fc2d_hps_quadnode(T& node_data) : data(node_data), level(0), parent(nullptr) {
			for (std::size_t i = 0; i < FC2D_HPS_NUMBER_CHILDREN; i++) {
				children[i] = nullptr;
			}
		}
	};

	fc2d_hps_quadnode* root;	// Root of tree
	std::size_t height;			// Height of tree, i.e., number of levels

	fc2d_hps_quadtree() :
		root(nullptr), height(0)
			{}

	fc2d_hps_quadtree(T& root_data) :
		root(new fc2d_hps_quadnode(root_data)), height(1)
			{}

	fc2d_hps_quadtree(const fc2d_hps_quadtree& to_copy) {
		this->root = to_copy.root;
		this->height = to_copy.height;
	}

	fc2d_hps_quadtree& operator=(const fc2d_hps_quadtree& rhs) {
		return *this;
	}

	~fc2d_hps_quadtree() {
		// printf("QUADTREE DESTRUCTOR CALLED\n");
		this->remove();
		this->root = nullptr;
	}

	void build(std::function<bool(T&)> bigger, std::function<std::vector<T>(T&)> init) {
		// Build a quadtree by calling `bigger` which returns true or false to determine if the tree gets bigger.
		// Any new leaves are built by calling `init`

		if (this->root == nullptr) {
			throw std::invalid_argument("[fc2d_hps_quadtree::build] Tree's `root` is null. Build tree from filled in root.");
		}

		std::size_t current_height = 0;
		build_(this->root, bigger, init, current_height);
	}

	void grow(fc2d_hps_quadnode* node, T children_data[FC2D_HPS_NUMBER_CHILDREN]) {
		// Grows a leaf node with supplied data
		// Check if node is a leaf
		if (node->children[0] != nullptr) {
			throw std::invalid_argument("[fc2d_hps_quadtree::grow] `node` is not a leaf");
		}
		// Get node's level
		// std::size_t node_level = node->level;
		
		// Create four children nodes of contiguous data
		fc2d_hps_quadnode* children_nodes[FC2D_HPS_NUMBER_CHILDREN];
		for (std::size_t i = 0; i < FC2D_HPS_NUMBER_CHILDREN; i++) {
			// Create new quadnode with data
			fc2d_hps_quadnode* new_quadnode = new fc2d_hps_quadnode(children_data[i]);
			children_nodes[i] = new_quadnode;
			
			// Point children nodes to parent
			children_nodes[i]->parent = node;

			// Point parent to children
			node->children[i] = children_nodes[i];

			// Set children height
			node->children[i]->level = node->level + 1;

			// Set tree's height
			if (node->children[i]->level > this->height - 1) {
				this->height++;
			}
		}
	}

	// TODO: Write grow function that takes a function that takes the parent node and returns 4 children

	void trim(fc2d_hps_quadnode* siblings[FC2D_HPS_NUMBER_CHILDREN]) {
		// Trims a siblings of leaves

		// Check if all siblings are leaves
		if (
			is_leaf_(siblings[FC2D_HPS_QUADTREE_LOWER_LEFT]) == false ||
			is_leaf_(siblings[FC2D_HPS_QUADTREE_LOWER_RIGHT]) == false ||
			is_leaf_(siblings[FC2D_HPS_QUADTREE_UPPER_LEFT]) == false ||
			is_leaf_(siblings[FC2D_HPS_QUADTREE_UPPER_RIGHT]) == false
		) {
			throw std::invalid_argument("[fc2d_hps_quadtree::trim] `siblings` are not all leaves");
		}

		// Get node's level
		// std::size_t node_level = siblings[0]->level;

		// Update tree height if trimming a level
		// if ((node_level + 1)

		for (std::size_t i = 0; i < FC2D_HPS_NUMBER_CHILDREN; i++) {
			// Point siblings' parent's children to null
			siblings[i]->parent->children[i] = nullptr;

			// Point siblings' parent to null
			// siblings[i]->parent = nullptr;

			// Update tree hieght
			

			// Delete siblings
			delete siblings[i];
		}
	}

	void remove() {
		// Recursively trim the tree until removed
		remove_(this->root);
	}

	void traverse_inorder(std::function<void(T&)> visit) {
		// Traverse the tree by visiting all children prior to visiting node (Inorder traversal)
		traverse_inorder_(this->root, visit);
		visit(this->root->data);
	}

	void traverse_preorder(std::function<void(T&)> visit) {
		// Traverse the tree by visitiing all nodes with parents visited first (Preorder traversal)
		// visit(this->root->data);
		traverse_preorder_(this->root, visit);
	}

	void merge(std::function<void(T&, T&, T&, T&, T&)> visit) {
		// Traverse the tree and merge all children up the tree
		// `visit` function should !replace! the data in the parent with a merged version of the children data
		merge_(this->root, visit);
	}

	void split(std::function<void(T&, T&, T&, T&, T&)> visit) {
		// Traverse the tree in preorder and split all parents into children
		// `visit` function should !replace! the data in the children with a split version of the parent data
		split_(this->root, visit);
	}

	// TODO: Create delete function

private:

	void build_(fc2d_hps_quadnode* node, std::function<bool(T&)> bigger, std::function<std::vector<T>(T&)> init, std::size_t current_height) {
		if (current_height >= FC2D_HPS_QUADTREE_MAX_HEIGHT) {
			std::cerr << "[fc2d_hps_quadtree::build] WARNING: Tree exceeds `FC2D_HPS_QUADTREE_MAX_HEIGHT`. Change max hieght or adjust `bigger` function." << std::endl;
			return;
		}
		if (bigger(node->data)) {
			std::vector<T> children = init(node->data);
			grow(node, children.data());
			for (std::size_t i = 0; i < FC2D_HPS_NUMBER_CHILDREN; i++) {
				build_(node->children[i], bigger, init, current_height++);
			}
			return;
		}
	}

	bool is_leaf_(fc2d_hps_quadnode* node) {
		for (std::size_t i = 0; i < FC2D_HPS_NUMBER_CHILDREN; i++) {
			if (node->children[i] != nullptr) {
				return false;
			}
		}
		return true;
	}

	void traverse_inorder_(fc2d_hps_quadnode* node, std::function<void(T&)> visit) {
		for (std::size_t i = 0; i < FC2D_HPS_NUMBER_CHILDREN; i++) {
			if (node->children[i] != nullptr) {
				traverse_inorder_(node->children[i], visit);
				visit(node->children[i]->data);
			}
		}
	}

	void traverse_preorder_(fc2d_hps_quadnode* node, std::function<void(T&)> visit) {
		for (std::size_t i = 0; i < FC2D_HPS_NUMBER_CHILDREN; i++) {
			if (node->children[i] != nullptr) {
				visit(node->children[i]->data);
				traverse_preorder_(node->children[i], visit);
			}
		}
	}

	void merge_(fc2d_hps_quadnode* node, std::function<void(T&, T&, T&, T&, T&)> visit) {
		if (is_leaf_(node) == true) {
			return;
		}
		else {
			for (std::size_t i = 0; i < FC2D_HPS_NUMBER_CHILDREN; i++) {
				if (is_leaf_(node->children[i]) == false) {
					merge_(node->children[i], visit);
				}
			}
			visit(
				node->data,
				node->children[FC2D_HPS_QUADTREE_LOWER_LEFT]->data,
				node->children[FC2D_HPS_QUADTREE_LOWER_RIGHT]->data,
				node->children[FC2D_HPS_QUADTREE_UPPER_LEFT]->data,
				node->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]->data
			);
		}
	}

	void split_(fc2d_hps_quadnode* node, std::function<void(T&, T&, T&, T&, T&)> visit) {
		if (is_leaf_(node) == true) {
			return;
		}
		else {
			visit(
				node->data,
				node->children[FC2D_HPS_QUADTREE_LOWER_LEFT]->data,
				node->children[FC2D_HPS_QUADTREE_LOWER_RIGHT]->data,
				node->children[FC2D_HPS_QUADTREE_UPPER_LEFT]->data,
				node->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]->data
			);
			for (std::size_t i = 0; i < FC2D_HPS_NUMBER_CHILDREN; i++) {
				split_(node->children[i], visit);
			}
		}
	}

	void remove_(fc2d_hps_quadnode* node) {
		if (is_leaf_(node) == true) {
			return;
		}
		else if (
			is_leaf_(node->children[FC2D_HPS_QUADTREE_LOWER_LEFT]) == true &&
			is_leaf_(node->children[FC2D_HPS_QUADTREE_LOWER_RIGHT]) == true &&
			is_leaf_(node->children[FC2D_HPS_QUADTREE_UPPER_LEFT]) == true &&
			is_leaf_(node->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]) == true
		) {
			trim(node->children);
			return;
		}
		else {
			for (std::size_t i = 0; i < FC2D_HPS_NUMBER_CHILDREN; i++) {
				remove_(node->children[i]);
			}
		}
	}

};

#endif // FC2D_HPS_QUADTREE_HPP
