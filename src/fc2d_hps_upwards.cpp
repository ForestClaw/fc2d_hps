#include "fc2d_hps_upwards.hpp"

void upwards_assign(fc2d_hps_patch& patch) {

    if (patch.f.size() == 0) {
        throw std::invalid_argument("[fc2d_hps_upwards::upwards_assign] Patch's RHS data `f` is not set.");
    }

    if (patch.is_leaf == true) {
        // Compute derivative of local particular solution
        fc2d_hps_vector<double> g_zeros(2*patch.grid.Nx + 2*patch.grid.Ny, 0);
        fc2d_hps_FISHPACK_solver solver;
        patch.h = solver.dtn(patch.grid, g_zeros, patch.f);
    }
    else {
        
    }

}