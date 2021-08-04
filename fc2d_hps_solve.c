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

#include "fc2d_hps_solve.h"

#include "fc2d_hps.h"

struct hps_tree_T_storage_t;

typedef struct hps_tree_T_storage
{
    double *Tmat;  /* T Matrix for parent  */
    int is_leaf;
    struct hps_tree_T_storage *children;  /* Pointer to first of four quadrants */
} hpt_tree_T_storage_t;


void cb_init_storage(fclaw2d_domain_t *domain,
                     fclaw2d_patch_t *fine_patches,
                     int blockno, int fine0_patchno,
                     void *user)

{
    fclaw2d_global_iterate_t* g = (fclaw2d_global_iterate_t*) user;

    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(g->glob);
    hps_tree_T_storage_t *hps_store_parent = ((hps_tree_T_storage_t*) g->user);

    fclaw2d_patch_t *patch0 = &fine_patches[0];
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch0,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    hps_store_parent->Tmat = FCLAW_ALLOC(double,2*(mx + my));
    hps_store_parent->is_leaf = 0;  /* Parent is not a leaf */

    for (int igrid = 0; igrid < 4; igrid++)
    {
        hps_tree_T_storage_t *child = hps_store_parent->children[igrid];
        child->Tmat = FCLAW_ALLOC(double,2*(mx + my));
        child->is_leaf = 1;  /* Assume children are leaves for now */
    }
}

static
void cb_merge(fclaw2d_domain_t *domain,
              fclaw2d_patch_t *fine_patches,
              int blockno,
              int fine0_patchno,
              void *user)
{
    fclaw2d_global_iterate_t* g = (fclaw2d_global_iterate_t*) user;

    hps_tree_T_storage_t *hps_store_parent = ((hps_tree_T_storage_t*) g->user);

    for(int igrid = 0; igrid < 4; igrid++)
    {
        int patchno = fine0_patchno+igrid;
        hps_tree_T_storage_t *child = hps_store_parent->children[igrid];
        if (child->is_leaf) 
            fc2d_hps_construct_T_leaf(glob,&fine_patch[igrid],blockno,patchno,child->Tmat);
    }

    /* Merge to four children into single parent T */
    fc2d_hps_merge_children(glob,*fine_patch[0],blockno,fine0_patch,hps_store_parent);
}
    

void fc2d_hps_solve(fclaw2d_global_t* glob)
{
    /* Apply non-homogeneous boundary conditions to any patches on the boundary */
    fc2d_hps_physical_bc(glob);

    /* Check that we only have one level */
    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    FCLAW_ASSERT(fclaw_opt->maxlevel <= 1);   

#if 0
    /* Initialize storage for T matrices */
    hps_tree_T_storage_t Tstore;   /* Pointer to  node zero */
    fclaw2d_global_iterate_families(glob, cb_init_storage, &Tstore);
#endif

    /* Assign patch solver */
    fc2d_hps_vtable_t  *hps_vt  = fc2d_hps_vt();  
    fc2d_hps_options_t *hps_opt = fc2d_hps_get_options(glob);
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

    /* Create T matrices for each leaf and merge to create T for parent leaves */
    fclaw2d_global_iterate_families(glob, cb_merge, &Tstore);


#if 0
    /* We will do something like this ... */
    double *Tmat;   /* Pointer to 
    fclaw2d_global_iterate_families(glob, hps_vt->patch_solver,Tmat);
                                    
#endif                                    

    // hps_vt->patch_solver(glob);
}

#if 0
void fc2d_hps_solve(fclaw2d_global_t *glob)
{
    fclaw2d_clawpatch_options_t *clawpatch_opt =
                               fclaw2d_clawpatch_get_options(glob);

    fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

    fc2d_hps_options_t *hps_opt = fc2d_hps_get_options(glob);

    /* Solver goes here ... */

    /* Factor : Use tree iterator */

    /* Solver : ... */
}
#endif


