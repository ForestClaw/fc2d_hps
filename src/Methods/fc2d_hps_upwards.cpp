#include <Methods/fc2d_hps_upwards.hpp>

void visit_upwards(fc2d_hps_patch& patch) {
    patch.print_info();
}

void fc2d_hps_upwards(fclaw2d_global_t* glob) {
    fclaw_global_essentialf("Begin HPS upwards pass\n");

    // Traverse inorder to set non-homogeneous RHS and BC data into solution vector
    // quadtree.traverse_inorder(visit_upwards);
    fclaw_global_essentialf("End HPS upwards pass\n");
}