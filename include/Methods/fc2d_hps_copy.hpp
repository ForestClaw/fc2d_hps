#ifndef FC2D_HPS_COPY_HPP_
#define FC2D_HPS_COPY_HPP_

#include <iostream>

#include "fclaw2d_global.h"
#include "fclaw2d_patch.h"
#include "fclaw2d_physical_bc.h"
#include "fclaw2d_clawpatch.h"

#include <HPS/fc2d_hps.hpp>
#include <Util/fc2d_hps_physical_bc.h>
#include <Util/fc2d_hps_options.h>
#include <Structures/fc2d_hps_patchgrid.hpp>
#include <Structures/fc2d_hps_patch.hpp>
#include <Structures/fc2d_hps_quadtree.hpp>
#include <Structures/fc2d_hps_patchsolver.hpp>

#include "fc2d_hps_methods.hpp"

// Globals
extern patch_tree quadtree; // use static
extern int current_ID; // use static
extern std::vector<fc2d_hps_matrix<double>> T_cache;

void visit_copy_data(fc2d_hps_patch& patch);
void fc2d_hps_clawpatch_data_move(fclaw2d_global* glob);

#endif // FC2D_HPS_COPY_HPP_