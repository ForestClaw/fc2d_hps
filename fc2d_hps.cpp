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

#include "fc2d_hps.h"
#include "fc2d_hps_options.h"
#include "fc2d_hps_physical_bc.h"
#include "fc2d_hps_fort.h"

#include <fclaw2d_elliptic_solver.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_output_ascii.h> 
#include <fclaw2d_clawpatch_output_vtk.h>

#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_output.h>

#include <fclaw2d_domain.h>

#include "fc2d_hps_solve.h"



static fc2d_hps_vtable_t s_hps_vt;

/* --------------------- Hps solver (required) ------------------------- */

static
void hps_setup_solver(fclaw2d_global_t *glob)
{
	//fc2d_hps_vtable_t*  hps_vt = fc2d_hps_vt();
}


static
void hps_rhs(fclaw2d_global_t *glob,
             fclaw2d_patch_t *patch,
             int blockno,
             int patchno)
{
	fc2d_hps_vtable_t*  hps_vt = fc2d_hps_vt();

	FCLAW_ASSERT(hps_vt->fort_rhs != NULL); /* Must be initialized */

    int mx,my,mbc;
    double dx,dy,xlower,ylower;
	fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

    int mfields;
    double *rhs;
	fclaw2d_clawpatch_rhs_data(glob,patch,&rhs,&mfields);
	FCLAW_ASSERT(mfields == 1);

	/* Compute right hand side */
	hps_vt->fort_rhs(&blockno,&mbc,&mx,&my,&mfields,
                    &xlower,&ylower,&dx,&dy,rhs);
}

static
void hps_solve(fclaw2d_global_t* glob)
{
    // Apply non-homogeneous boundary conditions 
    fc2d_hps_physical_bc(glob);

    fc2d_hps_vtable_t  *hps_vt  = fc2d_hps_vt();  
    fc2d_hps_options_t *hps_opt = fc2d_hps_get_options(glob);

    /* Should the solvers be part of the hps library? Yes, for now, at least */
    switch (hps_opt->patch_solver)
    {
        case FISHPACK:
            // hps_vt->patch_solver = fc2d_hps_fishpack_solve;
            break;
#if 0
        case DST:
            hps_vt->patch_solver = fc2d_hps_dst_solve;
            break;
        case VARPOISSON:
            hps_vt->patch_solver = fc2d_hps_varpoisson_solve;
            break;
        case HEAT:
            hps_vt->patch_solver = fc2d_hps_heat_solve;
            break;
        case USER_solver:
            if (hps_vt->patch_solver == NULL)
            {
                fclaw_global_essentialf("hps_solve : User specified solver not set\n");
                exit(0);
            }
#endif            
        default:
            break;
            /* user has specified something, hopefully */
    }
    
    FCLAW_ASSERT(hps_vt->patch_solver != NULL);

    hps_vt->patch_solver(glob);
}

/* ---------------------------------- Output functions -------------------------------- */

static
void hps_output(fclaw2d_global_t *glob, int iframe)
{
	const fc2d_hps_options_t* hps_opt;
	hps_opt = fc2d_hps_get_options(glob);

	if (hps_opt->ascii_out != 0)
		fclaw2d_clawpatch_output_ascii(glob,iframe);

	if (hps_opt->vtk_out != 0)
		fclaw2d_clawpatch_output_vtk(glob,iframe);
}



/* ------------------------------ Virtual functions  ---------------------------------- */

static
fc2d_hps_vtable_t* hps_vt_init()
{
	FCLAW_ASSERT(s_hps_vt.is_set == 0);
	return &s_hps_vt;
}

void fc2d_hps_solver_initialize()
{
	int claw_version = 4; /* solution data is organized as (i,j,m) */
	fclaw2d_clawpatch_vtable_initialize(claw_version);


    //fclaw2d_clawpatch_vtable_t*      clawpatch_vt = fclaw2d_clawpatch_vt();

	/* ForestClaw vtable items */
	fclaw2d_vtable_t*   fclaw_vt = fclaw2d_vt();
	fclaw_vt->output_frame      = hps_output;

	/* These could be over-written by user specific settings */
	fclaw2d_patch_vtable_t*   patch_vt = fclaw2d_patch_vt();  
	patch_vt->rhs            = hps_rhs;  /* Calls FORTRAN routine */
	patch_vt->setup          = NULL;
    
    fclaw2d_elliptic_vtable_t *elliptic_vt = fclaw2d_elliptic_vt();
    elliptic_vt->setup = hps_setup_solver;
    elliptic_vt->solve = hps_solve;    
    elliptic_vt->apply_bc = fc2d_hps_physical_bc;

	fc2d_hps_vtable_t*  hps_vt = hps_vt_init();	
    hps_vt->fort_apply_bc = &HPS_FORT_APPLY_BC_DEFAULT;
    hps_vt->fort_eval_bc  = &HPS_FORT_EVAL_BC_DEFAULT;

#if 0
    /* solver is specified in solve routine, above*/
    hps_vt->patch_solver = NULL;
#endif    

	hps_vt->is_set = 1;
}


/* ----------------------------- User access to solver functions --------------------------- */

fc2d_hps_vtable_t* fc2d_hps_vt()
{
	FCLAW_ASSERT(s_hps_vt.is_set != 0);
	return &s_hps_vt;
}





