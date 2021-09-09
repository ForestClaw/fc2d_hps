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
		fc2d_hps_quadnode* children[FC2D_HPS_NUMBER_CHILDREN];
		fc2d_hps_quadnode* parent;
		fc2d_hps_quadnode() {}
		fc2d_hps_quadnode(T& node_data) : data(node_data), parent(nullptr) {
			for (std::size_t i = 0; i < FC2D_HPS_NUMBER_CHILDREN; i++) {
				children[i] = nullptr;
			}
		}
	};

public:

	fc2d_hps_quadnode* root;

	fc2d_hps_quadtree() :
		root(nullptr)
			{}

	fc2d_hps_quadtree(T& root_data) :
		root(new fc2d_hps_quadnode(root_data))
			{}

	void grow(fc2d_hps_quadnode* node, T children_data[FC2D_HPS_NUMBER_CHILDREN]) {
		// Grows a leaf node with supplied data

		// Check if node is a leaf
		if (node->children[0] != nullptr) {
			throw std::invalid_argument("[fc2d_hps_quadtree::grow] `node` is not a leaf");
		}

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
		}
	}

	void trim(fc2d_hps_quadnode* siblings[FC2D_HPS_NUMBER_CHILDREN]) {
		// Trims a siblings of leaves

		// Check if all siblings are leaves
		if (
			siblings[0]->children[0] != nullptr ||
			siblings[1]->children[0] != nullptr ||
			siblings[2]->children[0] != nullptr ||
			siblings[3]->children[0] != nullptr
		) {
			throw std::invalid_argument("[fc2d_hps_quadtree::trim] `siblings` are note all leaves");
		}

		for (std::size_t i = 0; i < FC2D_HPS_NUMBER_CHILDREN; i++) {
			// Point siblings' parent's children to null
			siblings[i]->parent->children[i] = nullptr;

			// Delete siblings
			delete siblings[i];
		}
		// Delete or point array of siblings...?
	}

};

#endif // FC2D_HPS_QUADTREE_HPP