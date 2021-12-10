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

#include <iostream>
#include "fclaw2d_global.h"
#include "fc2d_hps.h"
#include "fc2d_hps_options.h"
#include "fc2d_hps_patchgrid.hpp"
#include "fc2d_hps_patch.hpp"
#include "fc2d_hps_patchsolver.hpp"

// struct fclaw2d_global; 

static void cb_merge(fclaw2d_global_t *glob, fclaw2d_patch_t *fine_patches, int blockno, int fine0_patchno, void *user);
void fc2d_hps_build(struct fclaw2d_global* glob);
void fc2d_hps_solve(struct fclaw2d_global* glob);

#endif /* !FC2D_HPS_SOLVE_HPP_ */

