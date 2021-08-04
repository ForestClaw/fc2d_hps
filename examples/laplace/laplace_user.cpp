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

#include "laplace_user.h"

static
void laplace_problem_setup(fclaw2d_global_t *glob)
{
    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        const laplace_options_t* user = laplace_get_options(glob);

        fprintf(f,  "%-24d   %s",  user->example,    "% example\n");
        fprintf(f,  "%-24.6f   %s",user->alpha,      "% alpha\n");
        fprintf(f,  "%-24.6f   %s",user->x0,         "% x0\n");
        fprintf(f,  "%-24.6f   %s",user->y0,         "% y0\n");
        fprintf(f,  "%-24.6f   %s",user->a,          "% a\n");
        fprintf(f,  "%-24.6f   %s",user->b,          "% b\n");
        fprintf(f,  "%-24.6f   %s",user->eps_disk,   "% eps_disk\n");
        fprintf(f,  "%-24d   %s",user->m_polar,    "% m_polar\n");
        for(int i = 0; i < user->m_polar; i++)
            fprintf(f,"%-24.6f   %% x0[%d]\n",user->x0_polar[i],i); 

        for(int i = 0; i < user->m_polar; i++)
            fprintf(f,"%-24.6f   %% y0[%d]\n",user->y0_polar[i],i);            

        for(int i = 0; i < user->m_polar; i++)
            fprintf(f,"%-24.6f   %% r0[%d]\n",user->r0_polar[i],i);            

        for(int i = 0; i < user->m_polar; i++)
            fprintf(f,"%-24.6f   %% r1[%d]\n",user->r1_polar[i],i);            

        for(int i = 0; i < user->m_polar; i++)
            fprintf(f,"%-24d   %% n[%d]\n",user->n_polar[i],i);            

        fc2d_hps_options_t*  hps_opt = fc2d_hps_get_options(glob);    
        fprintf(f,  "%-24d   %s",hps_opt->boundary_conditions[0],  "% bc[0]\n");
        fprintf(f,  "%-24d   %s",hps_opt->boundary_conditions[1],  "% bc[1]\n");
        fprintf(f,  "%-24d   %s",hps_opt->boundary_conditions[2],  "% bc[2]\n");
        fprintf(f,  "%-24d   %s",hps_opt->boundary_conditions[3],  "% bc[3]\n");


        fclose(f);
    }
    fclaw2d_domain_barrier (glob->domain);

    SETPROB(); /* This file reads the file just created above */
}


void laplace_link_solvers(fclaw2d_global_t *glob)
{
    /* ForestClaw vtable */
    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt();
    fclaw_vt->problem_setup = &laplace_problem_setup;  

    /* HPS virtual table : Initialize RHS */
    fc2d_hps_vtable_t*  hps_vt = fc2d_hps_vt();
    hps_vt->fort_rhs = &LAPLACE_FORT_RHS;

    /* Clawpatch : Compute the error */
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();

    /* Specialized for this example */
    clawpatch_vt->fort_compute_patch_error = &LAPLACE_COMPUTE_ERROR;

    /* BCs : Include inhomogeneous boundary conditions on the right hand side */
    hps_vt->fort_apply_bc = &LAPLACE_FORT_APPLY_BC;
    hps_vt->fort_eval_bc  = &LAPLACE_FORT_EVAL_BC;

    // Output routines
    clawpatch_vt->time_header_ascii = fc2d_hps_time_header_ascii;
    clawpatch_vt->cb_output_ascii = cb_hps_output_ascii;     

    /* This one will be used only if we are also computing the error.  In this case
       the error and exact solution will be output, along with computed solution */          
    hps_vt->fort_output = &LAPLACE_FORT_OUTPUT_ASCII; 
}

