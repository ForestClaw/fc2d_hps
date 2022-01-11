#include "gtest/gtest.h"
#include <fclaw2d_include_all.h>
#include <fc2d_hps.h>
#include <fc2d_hps_quadtree.hpp>

TEST(QuadTree, init) {
	int root = 4;
	fc2d_hps_quadtree<int> tree(root);
	EXPECT_EQ(tree.root->data, root);
}

TEST(QuadTree, grow_root) {
	int root = 4;
	fc2d_hps_quadtree<int> tree(root);
	int children[4] = {3, 2, 1, 0};
	tree.grow(tree.root, children);

	EXPECT_EQ(tree.root->children[0]->data, children[0]);
	EXPECT_EQ(tree.root->children[1]->data, children[1]);
	EXPECT_EQ(tree.root->children[2]->data, children[2]);
	EXPECT_EQ(tree.root->children[3]->data, children[3]);

	EXPECT_EQ(tree.height, 2);
}

TEST(QuadTree, grow_tree) {
	/**
	 * Let's make a tree that looks like this:
	 * 					   0
	 * 				/	 |	 |	 \	 
	 * 				1  	 2 	 3 	  4 
	 * 			 / | | \        / | | \
	 *           5 6 7 8       9 10 11 12 
	 * and check to make sure it got built correctly.
	 */

	int root = 0;
	fc2d_hps_quadtree<int> tree(root);

	int children_data1[FC2D_HPS_NUMBER_CHILDREN] = {1, 2, 3, 4};
	int children_data2[FC2D_HPS_NUMBER_CHILDREN] = {5, 6, 7, 8};
	int children_data3[FC2D_HPS_NUMBER_CHILDREN] = {9, 10, 11, 12};
	tree.grow(tree.root, children_data1);
	
	EXPECT_EQ(tree.height, 2);

	tree.grow(tree.root->children[FC2D_HPS_QUADTREE_LOWER_LEFT], children_data2);

	EXPECT_EQ(tree.height, 3);

	tree.grow(tree.root->children[FC2D_HPS_QUADTREE_UPPER_RIGHT], children_data3);

	EXPECT_EQ(tree.height, 3);

	EXPECT_EQ(tree.root->data, 0);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_LOWER_LEFT]->data, 1);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_LOWER_RIGHT]->data, 2);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_UPPER_LEFT]->data, 3);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]->data, 4);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_LOWER_LEFT]->children[FC2D_HPS_QUADTREE_LOWER_LEFT]->data, 5);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_LOWER_LEFT]->children[FC2D_HPS_QUADTREE_LOWER_RIGHT]->data, 6);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_LOWER_LEFT]->children[FC2D_HPS_QUADTREE_UPPER_LEFT]->data, 7);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_LOWER_LEFT]->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]->data, 8);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]->children[FC2D_HPS_QUADTREE_LOWER_LEFT]->data, 9);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]->children[FC2D_HPS_QUADTREE_LOWER_RIGHT]->data, 10);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]->children[FC2D_HPS_QUADTREE_UPPER_LEFT]->data, 11);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]->data, 12);
}

TEST(QuadTree, grow_trim_tree) {
	/**
	 * Let's make a tree that looks like this:
	 * 					   0
	 * 				/	 |	 |	 \	 
	 * 				1  	 2 	 3 	  4 
	 * 			 / | | \        / | | \
	 *           5 6 7 8       9 10 11 12 
	 * and then trim {9, 10, 11, 12}.
	 */

	int root = 0;
	fc2d_hps_quadtree<int> tree(root);

	int children_data1[FC2D_HPS_NUMBER_CHILDREN] = {1, 2, 3, 4};
	int children_data2[FC2D_HPS_NUMBER_CHILDREN] = {5, 6, 7, 8};
	int children_data3[FC2D_HPS_NUMBER_CHILDREN] = {9, 10, 11, 12};
	tree.grow(tree.root, children_data1);

	EXPECT_EQ(tree.height, 2);

	tree.grow(tree.root->children[FC2D_HPS_QUADTREE_LOWER_LEFT], children_data2);

	EXPECT_EQ(tree.height, 3);

	tree.grow(tree.root->children[FC2D_HPS_QUADTREE_UPPER_RIGHT], children_data3);

	EXPECT_EQ(tree.height, 3);

	tree.trim(tree.root->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]->children);

	EXPECT_EQ(tree.height, 3);
	
	EXPECT_EQ(tree.root->data, 0);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_LOWER_LEFT]->data, 1);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_LOWER_RIGHT]->data, 2);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_UPPER_LEFT]->data, 3);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]->data, 4);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_LOWER_LEFT]->children[FC2D_HPS_QUADTREE_LOWER_LEFT]->data, 5);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_LOWER_LEFT]->children[FC2D_HPS_QUADTREE_LOWER_RIGHT]->data, 6);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_LOWER_LEFT]->children[FC2D_HPS_QUADTREE_UPPER_LEFT]->data, 7);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_LOWER_LEFT]->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]->data, 8);

	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]->children[FC2D_HPS_QUADTREE_LOWER_LEFT], nullptr);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]->children[FC2D_HPS_QUADTREE_LOWER_RIGHT], nullptr);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]->children[FC2D_HPS_QUADTREE_UPPER_LEFT], nullptr);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]->children[FC2D_HPS_QUADTREE_UPPER_RIGHT], nullptr);
}

void visit1(int& data) {
	std::cout << data << " ";
	return;
}

TEST(QuadTree, traverse) {
	int root = 0;
	fc2d_hps_quadtree<int> tree(root);

	int children_data1[FC2D_HPS_NUMBER_CHILDREN] = {1, 2, 3, 4};
	int children_data2[FC2D_HPS_NUMBER_CHILDREN] = {5, 6, 7, 8};
	int children_data3[FC2D_HPS_NUMBER_CHILDREN] = {9, 10, 11, 12};
	tree.grow(tree.root, children_data1);
	tree.grow(tree.root->children[FC2D_HPS_QUADTREE_LOWER_LEFT], children_data2);
	tree.grow(tree.root->children[FC2D_HPS_QUADTREE_UPPER_RIGHT], children_data3);

	tree.traverse(visit1);
}

void visit2(int& data, int& data0, int& data1, int& data2, int& data3) {
	std::cout
		<< "data  = " << data  << ", "
		<< "data0 = " << data0 << ", "
		<< "data1 = " << data1 << ", "
		<< "data2 = " << data2 << ", "
		<< "data3 = " << data3 << ", "
	<< std:: endl;
	
	data = data + data0 + data1 + data2 + data3;
	// return data0 + data1 + data2 + data3;
}

TEST(QuadTree, merge) {
	/**
	 * Let's make a tree that looks like this:
	 * 					   0
	 * 				/	 |	 |	 \	 
	 * 				1  	 2 	 3 	  4 
	 * 			 / | | \        / | | \
	 *           5 6 7 8       9 10 11 12
	 *                      /  |  |  \
	 *                     13 14 15 16 
	 * and then merge by adding the leaf values.
	 * 
	 * The result at the root should be (5 + 6 + 7 + 8) + 1 + 2 + 3 + (9 + 10 + 11 + 12) + 4 = 78.
	 * The result at the (1) node should be (5 + 6 + 7 + 8) + 1 = 27.
	 * The result at the (4) node should be (9 + 10 + 11 + 12) + 4 = 46.
	 */
	int root = 0;
	fc2d_hps_quadtree<int> tree(root);

	int children_data1[FC2D_HPS_NUMBER_CHILDREN] = {1, 2, 3, 4};
	int children_data2[FC2D_HPS_NUMBER_CHILDREN] = {5, 6, 7, 8};
	int children_data3[FC2D_HPS_NUMBER_CHILDREN] = {9, 10, 11, 12};
	tree.grow(tree.root, children_data1);
	tree.grow(tree.root->children[FC2D_HPS_QUADTREE_LOWER_LEFT], children_data2);
	tree.grow(tree.root->children[FC2D_HPS_QUADTREE_UPPER_RIGHT], children_data3);

	tree.merge(visit2);

	EXPECT_EQ(tree.root->data, 78);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_LOWER_LEFT]->data, 27);
	EXPECT_EQ(tree.root->children[FC2D_HPS_QUADTREE_UPPER_RIGHT]->data, 46);
}

bool bigger1(int& data) {
	printf("in bigger1: data = %i\n", data); 
	if (data < 2) return true;
	else return false;
}

std::vector<int> init1(int& data) {
	printf("in init1: data = %i\n", data);
	std::vector<int> children(FC2D_HPS_NUMBER_CHILDREN, data+1);
	for (int i = 0; i < FC2D_HPS_NUMBER_CHILDREN; i++) {
		printf("in init1: children[%i] = %i\n", i, children[i]);
	}
	return children;
}

TEST(QuadTree, build) {
	/**
	 * Testing the build function. We'll build a tree that looks like this:
	 * 					   0
	 * 				/	 |	 |	 \	 
	 * 				1  	 1 	 1 	  1 
	 * 			 / | | \ |   |  / | | \
	 *           2 2 2 2  ...  2  2 2 2
	 *
	 */

	int root = 0;
	fc2d_hps_quadtree<int> tree(root);
	tree.build(bigger1, init1);
	tree.traverse(visit1);

	EXPECT_EQ(tree.root->data, 0);
	for (int i = 0; i < FC2D_HPS_NUMBER_CHILDREN; i++) {
		EXPECT_EQ(tree.root->children[i]->data, 1);
		for (int j = 0; j < FC2D_HPS_NUMBER_CHILDREN; j++) {
			EXPECT_EQ(tree.root->children[i]->children[j]->data, 2);
		}
	}

}