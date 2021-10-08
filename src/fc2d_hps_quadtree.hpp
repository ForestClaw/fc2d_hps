#ifndef FC2D_HPS_QUADTREE_HPP
#define FC2D_HPS_QUADTREE_HPP

#include <forestclaw2d.h>
#include <p4est.h>
#include <p4est_wrap.h>
#include <p4est_iterate.h>
#include <fc2d_hps_patch.hpp>

#define FC2D_HPS_NUMBER_CHILDREN 		4
#define FC2D_HPS_QUADTREE_LOWER_LEFT 	0
#define FC2D_HPS_QUADTREE_LOWER_RIGHT 	1
#define FC2D_HPS_QUADTREE_UPPER_LEFT 	2
#define FC2D_HPS_QUADTREE_UPPER_RIGHT 	3

// Forward declaration
void p4est_iterate_fn(p4est_iter_volume_info_t* info, void* user_data);

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
	// std::size_t height;			// Height of tree, i.e., number of levels

	fc2d_hps_quadtree() :
		root(nullptr)
			{}

	fc2d_hps_quadtree(T& root_data) :
		root(new fc2d_hps_quadnode(root_data))
			{}
	
	// fc2d_hps_quadtree(fclaw2d_domain_t* domain) :
	// 	root(nullptr) {
	// 	// Build a quadtree from a forestclaw domain object

	// 	// Variables
	// 	p4est_wrap_t* wrap = (p4est_wrap_t*) domain->pp;
	// 	p4est_t* p4est = wrap->p4est;

	// 	// Begin growth algorithm
	// 	fc2d_hps_quadnode* temp_node = root;


	// }

	~fc2d_hps_quadtree() {}

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
		}

		// Update tree height if adding another level
		// if ((node_level + 1) >= this->height) {
		// 	this->height = node_level + 2;
		// }
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
			throw std::invalid_argument("[fc2d_hps_quadtree::trim] `siblings` are note all leaves");
		}

		// Get node's level
		// std::size_t node_level = siblings[0]->level;

		// Update tree height if trimming a level
		// if ((node_level + 1)

		for (std::size_t i = 0; i < FC2D_HPS_NUMBER_CHILDREN; i++) {
			// Point siblings' parent's children to null
			siblings[i]->parent->children[i] = nullptr;

			// Delete siblings
			delete siblings[i];
		}
		// Delete or point array of siblings...?


	}

	void traverse(std::function<void(T&)> visit) {
		// Traverse the tree by visiting all children prior to visiting node (Inorder traversal)
		traverse_(this->root, visit);
		visit(this->root->data);
	}

	void merge(std::function<void(T&, T&, T&, T&, T&)> visit) {
		// Traverse the tree and merge all children up the tree
		// `visit` function should !replace! the data in the parent with a merged version of the children data
		merge_(this->root, visit);
	}

private:

	bool is_leaf_(fc2d_hps_quadnode* node) {
		for (std::size_t i = 0; i < FC2D_HPS_NUMBER_CHILDREN; i++) {
			if (node->children[i] != nullptr) {
				return false;
			}
		}
		return true;
	}

	void traverse_(fc2d_hps_quadnode* node, std::function<void(T&)> visit) {
		for (std::size_t i = 0; i < FC2D_HPS_NUMBER_CHILDREN; i++) {
			if (node->children[i] != nullptr) {
				traverse_(node->children[i], visit);
				visit(node->children[i]->data);
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


};

fc2d_hps_quadtree<fc2d_hps_patch> fc2d_hps_create_quadtree_from_domain(fclaw2d_domain_t* domain) {
	// WORKING HERE!
	// My current approach doesn't look like it will work... I believe I need to come up with a recursive approach.


	// Create and initialize variables
	fc2d_hps_patch root_patch;
	fc2d_hps_quadtree<fc2d_hps_patch> tree(root_patch);
	fc2d_hps_quadtree<fc2d_hps_patch>::fc2d_hps_quadnode** temp_node = &(tree.root);
	p4est_wrap_t* wrap = (p4est_wrap_t*) domain->pp;
	p4est_t* p4est = wrap->p4est;
	p4est_tree_t* p4est_tree = p4est_tree_array_index(p4est->trees, 0);

	// Iterate through quadrants in tree in Morton order
	for (std::size_t ID = 0; ID < p4est->local_num_quadrants; ID++) {

		// Get access to quadrant
		p4est_quadrant_t* p4est_quad = p4est_quadrant_array_index(&(p4est_tree->quadrants), ID);

		// // Get levels of node and quad
		// int node_level = temp_node->level;
		// int quad_level = p4est_quad->level;

		while ((*temp_node)->level != p4est_quad->level) {

			// Create four new nodes of patches
			fc2d_hps_patch child_patches[FC2D_HPS_NUMBER_CHILDREN];
			child_patches[0] = *(new fc2d_hps_patch());
			child_patches[1] = *(new fc2d_hps_patch());
			child_patches[2] = *(new fc2d_hps_patch());
			child_patches[3] = *(new fc2d_hps_patch());
			
			// Grow tree
			tree.grow(*temp_node, child_patches);

			// Move temp_node to point to first child of new family
			temp_node = &((*temp_node)->children[0]);

		}

		// temp_node should point to a leaf corresponding to the correct level p4est_quadrant
		// Store RHS data on leaf patches
		// Access RHS data from forestclaw domain

	}


	return tree;

}

#endif // FC2D_HPS_QUADTREE_HPP