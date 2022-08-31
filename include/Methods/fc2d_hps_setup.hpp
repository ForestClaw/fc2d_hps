#ifndef FC2D_HPS_SETUP_HPP_
#define FC2D_HPS_SETUP_HPP_

#include <iostream>

#include <p4est_wrap.h>
#include <p4est_iterate.h>

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
#include <Util/GenericSingleton.hpp>
#include <Structures/fc2d_hps_datacache.hpp>

// Globals
typedef fc2d_hps_quadtree<fc2d_hps_patch> patch_tree;

fc2d_hps_patch init_fn(fc2d_hps_quadtree<fc2d_hps_patch>* quadtree, int level, int idx, void* user);
void visit_set_ID(fc2d_hps_patch& patch);
void fc2d_hps_setup(struct fclaw2d_global* glob);

#endif // FC2D_HPS_SETUP_HPP_