#ifndef FC2D_HPS_SETUP_HPP_
#define FC2D_HPS_SETUP_HPP_

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

// Globals
typedef fc2d_hps_quadtree<fc2d_hps_patch> patch_tree;
extern patch_tree quadtree; // use static
extern int current_ID; // use static

bool build_from_p4est_callback_bigger(fc2d_hps_patch& patch);
std::vector<fc2d_hps_patch> build_from_p4est_callback_init(fc2d_hps_patch& parent);
// fc2d_hps_quadtree<fc2d_hps_patch> fc2d_hps_create_quadtree_from_domain(fclaw2d_global_t* glob);
void fc2d_hps_create_quadtree_from_domain(fclaw2d_global_t* glob);
void fc2d_hps_setup(struct fclaw2d_global* glob);

#endif // FC2D_HPS_SETUP_HPP_