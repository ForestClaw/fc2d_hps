/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun, Damyn Chipman
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

#ifndef FC2D_HPS_PHYSICAL_BC_H
#define FC2D_HPS_PHYSICAL_BC_H

#include <fclaw2d_elliptic_solver.h>
#include <fclaw2d_global.h>
#include <fclaw2d_domain.h>
#include <fclaw2d_patch.h>
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_physical_bc.h>

#include <HPS/fc2d_hps.hpp>
#include <Util/fc2d_hps_options.h>

#ifdef __cplusplus
extern "C"
{
#endif

struct fclaw2d_global;
struct fclaw2d_domain;
struct fclaw2d_patch;


typedef struct fc2d_hps_time_info
{
    double t;
} fc2d_hps_time_info_t;


/* Inhomogeneous boundary conditions */
void fc2d_hps_physical_bc(struct fclaw2d_global *glob);


#ifdef __cplusplus
}
#endif

#endif