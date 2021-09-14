#ifndef FC2D_HPS_QUADTREE_HPP
#define FC2D_HPS_QUADTREE_HPP

#define FC2D_HPS_NUMBER_CHILDREN 		4
#define FC2D_HPS_QUADTREE_LOWER_LEFT 	0
#define FC2D_HPS_QUADTREE_LOWER_RIGHT 	1
#define FC2D_HPS_QUADTREE_UPPER_LEFT 	2
#define FC2D_HPS_QUADTREE_UPPER_RIGHT 	3

template<class T>
class fc2d_hps_quadtree {

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

public:

	fc2d_hps_quadnode* root;	// Root of tree
	// std::size_t height;			// Height of tree, i.e., number of levels

	fc2d_hps_quadtree() :
		root(nullptr)
			{}

	fc2d_hps_quadtree(T& root_data) :
		root(new fc2d_hps_quadnode(root_data))
			{}

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

#endif // FC2D_HPS_QUADTREE_HPP