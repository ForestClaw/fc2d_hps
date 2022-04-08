/*
  Copyright (c) 2019-2021 Carsten Burstedde, Donna Calhoun, Damyn Chipman
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef FC2D_HPS_SOLVE_HPP_
#define FC2D_HPS_SOLVE_HPP_

// #include <iostream>

// #include "fclaw2d_global.h"
// #include "fclaw2d_patch.h"
// #include "fclaw2d_physical_bc.h"
// #include "fclaw2d_clawpatch.h"

// #include <HPS/fc2d_hps.hpp>
// #include <Util/fc2d_hps_physical_bc.h>
// #include <Util/fc2d_hps_options.h>
// #include <Structures/fc2d_hps_patchgrid.hpp>
// #include <Structures/fc2d_hps_patch.hpp>
// #include <Structures/fc2d_hps_quadtree.hpp>
// #include <Structures/fc2d_hps_patchsolver.hpp>
// #include <Methods/fc2d_hps_merge.hpp>
// #include <Methods/fc2d_hps_split.hpp>
#include "fc2d_hps_methods.hpp"

extern patch_tree quadtree; // use static
extern int current_ID; // use static

// HPS Setup Routines
// bool build_from_p4est_callback_bigger(fc2d_hps_patch& patch);
// std::vector<fc2d_hps_patch> build_from_p4est_callback_init(fc2d_hps_patch& parent);
// // fc2d_hps_quadtree<fc2d_hps_patch> fc2d_hps_create_quadtree_from_domain(fclaw2d_global_t* glob);
// void fc2d_hps_create_quadtree_from_domain(fclaw2d_global_t* glob);
// void fc2d_hps_setup(struct fclaw2d_global* glob);

// HPS Build Routines
// static void cb_merge(fclaw2d_global_t *glob, fclaw2d_patch_t *fine_patches, int blockno, int fine0_patchno, void *user);
// void visit_leaves(fc2d_hps_patch& patch);
// void visit_merge(fc2d_hps_patch& tau, fc2d_hps_patch& alpha, fc2d_hps_patch& beta, fc2d_hps_patch& gamma, fc2d_hps_patch& omega);
// void fc2d_hps_build(struct fclaw2d_global* glob);

// HPS Upwards Pass Routines
// void visit_upwards(fc2d_hps_patch& patch);
// void fc2d_hps_upwards(struct fclaw2d_global* glob);

// HPS Solve Routine
void visit_split(fc2d_hps_patch& tau, fc2d_hps_patch& alpha, fc2d_hps_patch& beta, fc2d_hps_patch& gamma, fc2d_hps_patch& omega);
void visit_patchsolver(fc2d_hps_patch& patch);
void set_root_boundary_data(fc2d_hps_patch& root_patch);
void fc2d_hps_solve(struct fclaw2d_global* glob);

// Patch Solver Routine
void fc2d_hps_fishpack_solve(fclaw2d_global_t* glob);

// ForestClaw - HPS Interface Routines
// void visit_copy_data(fc2d_hps_patch& patch);
// void fc2d_hps_clawpatch_data_move(fclaw2d_global* glob);

#endif /* !FC2D_HPS_SOLVE_HPP_ */

