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
#include "fc2d_hps_solve.hpp"
#include "fc2d_hps_diagnostics.h"
#include "fc2d_hps_vector.hpp"

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
    int mx,my,mbc;
    double dx,dy,xlower,ylower;
	fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

    int mfields;
    double *rhs;
	fclaw2d_clawpatch_rhs_data(glob,patch,&rhs,&mfields);
	FCLAW_ASSERT(mfields == 1);

	/* Compute right hand side */
    fc2d_hps_vtable_t*  hps_vt = fc2d_hps_vt();
    FCLAW_ASSERT(hps_vt->fort_rhs != NULL); /* Must be initialized */

	hps_vt->fort_rhs(&blockno,&mbc,&mx,&my,&mfields,
                    &xlower,&ylower,&dx,&dy,rhs);
}

void visit_print(fc2d_hps_patch& patch) {
    patch.print_info();
}

// static fc2d_hps_quadtree<fc2d_hps_patch> quadtree;
// patch_tree quadtree;

static
void hps_solve(fclaw2d_global_t *glob)
{
    fclaw_global_essentialf("----- Begin HPS solver -----\n");

    // HPS setup stage
    fc2d_hps_physical_bc(glob); // These could be done independently of each other
    fc2d_hps_setup(glob);       // These could be done independently of each other

    // HPS build stage
    fc2d_hps_build(glob);

    // Upwards pass for non-homogeneous RHS and BCs
    // fc2d_hps_upwards(glob);

    // HPS solve stage
    fc2d_hps_solve(glob);

    // Copy data from leaf patches into clawpatches
    fc2d_hps_clawpatch_data_move(glob);

    fclaw_global_essentialf("----- End of HPS solver -----\n");
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


/* ---------------------------------- Tagging functions ------------------------------- */

int hps_tag4refinement(fclaw2d_global_t *glob,
                       fclaw2d_patch_t *this_patch,
                       int blockno, int patchno,
                       int initflag)
{
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

    int tag_patch;
    double refine_threshold;

    refine_threshold = fclaw_opt->refine_threshold;

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *rhs;
    int mfields;
    fclaw2d_clawpatch_rhs_data(glob,this_patch,&rhs,&mfields);

    fclaw2d_clawpatch_vtable_t* clawpatch_vt = fclaw2d_clawpatch_vt();
    FCLAW_ASSERT(clawpatch_vt->fort_tag4refinement != NULL);


    /* Use default fortran tagging routines.  Choose refinement based on criteria
       set in configuration files (clawpatch:refinement-criteria) */
    tag_patch = 0;
    clawpatch_vt->fort_tag4refinement(&mx,&my,&mbc,&mfields,&xlower,&ylower,&dx,&dy,
                                      &blockno, rhs, &refine_threshold,
                                      &initflag, &tag_patch);
    return tag_patch;
}


static
int hps_tag4coarsening(fclaw2d_global_t *glob,
                       fclaw2d_patch_t *fine_patches,
                       int blockno,
                       int patchno,
                       int initflag)
{
    fclaw2d_patch_t *patch0 = &fine_patches[0];
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch0,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *rhs[4];
    int mfields;
    for (int igrid = 0; igrid < 4; igrid++)
    {
        fclaw2d_clawpatch_rhs_data(glob,&fine_patches[igrid],&rhs[igrid],&mfields);
    }

    fclaw2d_clawpatch_vtable_t* clawpatch_vt = fclaw2d_clawpatch_vt();
    FCLAW_ASSERT(clawpatch_vt->fort_tag4coarsening != NULL);

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    double coarsen_threshold = fclaw_opt->coarsen_threshold;

    int tag_patch = 0;
    clawpatch_vt->fort_tag4coarsening(&mx,&my,&mbc,&mfields,&xlower,&ylower,&dx,&dy,
                                      &blockno, rhs[0],rhs[1],rhs[2],rhs[3],
                                      &coarsen_threshold,&initflag,&tag_patch);
    return tag_patch == 1;
}

/* -------------------------------- Diagnostic functions ------------------------------ */
static
void hps_compute_error(fclaw2d_global_t *glob,
                          fclaw2d_patch_t *patch,
                          int blockno,
                          int patchno,
                          void *user)
{
    fc2d_hps_error_info_t* error_data = (fc2d_hps_error_info_t*) user;

    fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    if (fclaw_opt->compute_error)
    {
        fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
        FCLAW_ASSERT(clawpatch_vt->fort_compute_patch_error != NULL);

        int mx, my, mbc;
        double xlower,ylower,dx,dy;
        fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,&xlower,
                                    &ylower,&dx,&dy);

        double *area = fclaw2d_clawpatch_get_area(glob,patch);  /* Might be null */

        /* Computing =olution is stored in the RHS; true solution is stored in soln */
        int mfields;

        /* Computed solution */
        double *rhs;
        fclaw2d_clawpatch_rhs_data(glob,patch,&rhs,&mfields);

        double *err;
        fclaw2d_clawpatch_elliptic_error_data(glob,patch,&err,&mfields);

        /* True solution */
        double *soln;
        fclaw2d_clawpatch_elliptic_soln_data(glob,patch,&soln,&mfields);
        double t = glob->curr_time;
        clawpatch_vt->fort_compute_patch_error(&blockno, &mx,&my,&mbc,
                                               &mfields,&dx,&dy,
                                               &xlower,&ylower, &t, rhs, err, soln);
        /* Accumulate sums and maximums needed to compute error norms */

        FCLAW_ASSERT(clawpatch_vt->fort_compute_error_norm != NULL);
        clawpatch_vt->fort_compute_error_norm(&blockno, &mx, &my, &mbc, &mfields, 
                                              &dx,&dy, area, err,
                                              error_data->local_error);

    }
}


static
void hps_conservation_check(fclaw2d_global_t *glob,
                            fclaw2d_patch_t *patch,
                            int blockno,
                            int patchno,
                            void *user)
{
    fc2d_hps_error_info_t* error_data = (fc2d_hps_error_info_t*) user;
    int mx, my, mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int mfields;
    double *rhs;  /* Solution is stored in the right hand side */ 
    fclaw2d_clawpatch_rhs_data(glob,patch,&rhs,&mfields);

    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    FCLAW_ASSERT(clawpatch_vt->fort_conservation_check != NULL);


    /* Need a better way to determine which diagnostic to do */
    double* area = fclaw2d_clawpatch_get_area(glob,patch);  
    clawpatch_vt->fort_conservation_check(&mx, &my, &mbc, &mfields, &dx,&dy,
                                          area, rhs, error_data->rhs,
                                          error_data->c_kahan);
    fc2d_hps_options_t *hps_opt = fc2d_hps_get_options(glob);

    int intersects_bc[4];
    fclaw2d_physical_get_bc(glob,blockno,patchno,intersects_bc);

    double t = glob->curr_time;
    int cons_check = 1;

    fc2d_hps_vtable_t*  hps_vt = fc2d_hps_vt();
    FCLAW_ASSERT(hps_vt->fort_apply_bc != NULL);

    /* Sum up the normal derivative around the boundary */
    hps_vt->fort_apply_bc(&blockno, &mx, &my, &mbc, &mfields, 
                         &xlower, &ylower, &dx,&dy,&t, intersects_bc,
                         hps_opt->boundary_conditions,rhs, hps_vt->fort_eval_bc,
                         &cons_check, error_data->boundary);
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

	/* Patch : These could be over-written by user specific settings */
	fclaw2d_patch_vtable_t*   patch_vt = fclaw2d_patch_vt();  
	patch_vt->rhs            = hps_rhs;   /* Calls FORTRAN routine */
    patch_vt->initialize     = hps_rhs;   /* Get an initial refinement */
	patch_vt->setup          = NULL;

    /* Tagging functions : Base refinement on the right hand side */
    patch_vt->tag4refinement = hps_tag4refinement;
    patch_vt->tag4coarsening = hps_tag4coarsening;

    /* Clawpatch and ForestClaw : Output functions */
    fclaw2d_vtable_t*   fclaw_vt = fclaw2d_vt();
    fclaw_vt->output_frame = hps_output;

    /* Elliptic specific functions */
    fclaw2d_elliptic_vtable_t *elliptic_vt = fclaw2d_elliptic_vt();
    elliptic_vt->setup = hps_setup_solver;

    /* Solver doesn't do anything so far */
    elliptic_vt->solve = hps_solve;    
    elliptic_vt->apply_bc = fc2d_hps_physical_bc;

    /* BCs : Homogeneous BCs by default */
	fc2d_hps_vtable_t*  hps_vt = hps_vt_init();	
    hps_vt->fort_apply_bc = &FC2D_HPS_FORT_APPLY_BC_DEFAULT;
    hps_vt->fort_eval_bc  = &FC2D_HPS_FORT_EVAL_BC_DEFAULT;

    /* Diagnostics : Error, conservation */
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    clawpatch_vt->compute_error = hps_compute_error;  /* calls user-defined fortran routine */

    /* Conservation check : Compares sum(rhs) with sum of normal fluxes around the boundary
       of the solution.   (uses divergence theorem) */
    clawpatch_vt->conservation_check = hps_conservation_check;        

    /* These are specialized for the elliptic problem */
    fclaw2d_diagnostics_vtable_t *diag_vt = fclaw2d_diagnostics_vt();
    diag_vt->patch_init_diagnostics     = fc2d_hps_diagnostics_initialize;
    diag_vt->patch_reset_diagnostics    = fc2d_hps_diagnostics_reset;
    diag_vt->patch_compute_diagnostics  = fc2d_hps_diagnostics_compute;
    diag_vt->patch_gather_diagnostics   = fc2d_hps_diagnostics_gather;
    diag_vt->patch_finalize_diagnostics = fc2d_hps_diagnostics_finalize;

	hps_vt->is_set = 1;
}


/* ----------------------------- User access to solver functions --------------------------- */

fc2d_hps_vtable_t* fc2d_hps_vt()
{
	FCLAW_ASSERT(s_hps_vt.is_set != 0);
	return &s_hps_vt;
}





