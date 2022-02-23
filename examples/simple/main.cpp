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

#include "simple_user.h"
#include <Structures/fc2d_hps_quadtree.hpp>
#include <iostream>

static
fclaw2d_domain_t* create_domain(sc_MPI_Comm mpicomm, fclaw_options_t* fclaw_opt)
{
    /* Mapped, multi-block domain */
    p4est_connectivity_t     *conn = NULL;
    fclaw2d_domain_t         *domain;
    fclaw2d_map_context_t    *cont = NULL, *brick = NULL;
 
    int mi = fclaw_opt->mi;
    int mj = fclaw_opt->mj;

    int a = fclaw_opt->periodic_x;
    int b = fclaw_opt->periodic_y;

    /* Map unit square to disk using mapc2m_disk.f */
    conn = p4est_connectivity_new_brick(mi,mj,a,b);
    brick = fclaw2d_map_new_brick(conn,mi,mj);
    cont = fclaw2d_map_new_nomap_brick(brick);

    domain = fclaw2d_domain_new_conn_map (mpicomm, fclaw_opt->minlevel, conn, cont);
    fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);  
    return domain;
}

double init_fn(p4est_t* p4est, void* user) {
    return 1.0;
}

void merge_fn(double& p, double& c0, double& c1, double& c2, double& c3) {
    p = c0 + c1 + c2 + c3;
}

void split_fn(double& p, double& c0, double& c1, double& c2, double& c3) {
    c0 = p / 4;
    c1 = p / 4;
    c2 = p / 4;
    c3 = p / 4;
}

void visit_print(double& d) {
    printf("node data: d = %f\n", d);
}

static
void run_program(fclaw2d_global_t* glob)
{
    std::cout << "[simple.cpp::run_program]  Running program" << std::endl;
    /* ---------------------------------------------------------------
       Set domain data.
       --------------------------------------------------------------- */
    fclaw2d_domain_data_new(glob->domain);


    /* Initialize virtual table for ForestClaw */
    fclaw2d_vtables_initialize(glob);

    /* Test hps solver */
    fc2d_hps_solver_initialize();

    /* set up elliptic solver to use the hps solver */
    simple_link_solvers(glob);

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */

    /* Set up grid and RHS */
    fclaw2d_initialize(glob);

    /* Compute sum of RHS; reset error accumulators */
    int init_flag = 1;  
    // fclaw2d_diagnostics_gather(glob,init_flag);
    init_flag = 0;

    /* Output rhs */
    int Frame = 0;
    // fclaw2d_output_frame(glob,Frame);
 
    /* Solve the elliptic problem */
    fclaw2d_elliptic_solve(glob);

    /* Compute error, compute conservation */
    fclaw2d_diagnostics_gather(glob, init_flag);

    /* Reset to time FISHPACK */
    // init_flag = 1;
    // fclaw2d_diagnostics_gather(glob, init_flag);

    /* Run FISHPACK */


    /* Output solution */
    Frame = 0;
    fclaw2d_output_frame(glob,Frame);

    /* ---------------------------------------------------------------
       Finalize
       --------------------------------------------------------------- */
    fclaw2d_finalize(glob);
}

int
main (int argc, char **argv)
{
    fclaw_app_t *app;
    int first_arg;
    fclaw_exit_type_t vexit;

    /* Options */
    sc_options_t                *options;
    fclaw_options_t             *fclaw_opt;

    fclaw2d_clawpatch_options_t *clawpatch_opt;
    fc2d_hps_options_t *hps_opt;
    simple_options_t *user_opt;

    fclaw2d_global_t *glob;
    fclaw2d_domain_t *domain;
    sc_MPI_Comm mpicomm;

    int retval;

    /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);

    /* Create new options packages */
    fclaw_opt = fclaw_options_register(app,"fclaw_options.ini");
    clawpatch_opt = fclaw2d_clawpatch_options_register(app,"fclaw_options.ini");
    hps_opt = fc2d_hps_options_register(app,"fclaw_options.ini");
    user_opt = simple_options_register(app,"fclaw_options.ini");  

    /* Read configuration file(s) and command line, and process options */
    options = fclaw_app_get_options (app);
    retval = fclaw_options_read_from_file(options);
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    /* Run the program */
    if (!retval & !vexit)
    {
        /* Options have been checked and are valid */

        mpicomm = fclaw_app_get_mpi_size_rank(app, NULL, NULL);
        domain = create_domain(mpicomm, fclaw_opt);
    
        /* Create global structure which stores the domain, timers, etc */
        glob = fclaw2d_global_new();
        fclaw2d_global_store_domain(glob, domain);

        /* Store option packages in glob */
        fclaw2d_options_store (glob, fclaw_opt);
        fclaw2d_clawpatch_options_store (glob, clawpatch_opt);
        fc2d_hps_options_store (glob, hps_opt);
        simple_options_store (glob, user_opt);
        // fclaw_global_essentialf("Running program\n");
        run_program(glob);
        fclaw_global_essentialf("Finished!\n");
        fclaw2d_global_destroy(glob);        
    }
    fclaw_global_essentialf("Destroying app\n");
    fclaw_app_destroy (app);
    printf("Returning...\n");

    return 0;
}
