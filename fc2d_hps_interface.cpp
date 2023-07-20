#include "fc2d_hps_interface.hpp"

using HPSAlgorithm = EllipticForest::HPSAlgorithm<EllipticForest::FISHPACK::FISHPACKFVGrid, EllipticForest::FISHPACK::FISHPACKFVSolver, EllipticForest::FISHPACK::FISHPACKPatch, double>;

/*************************************************/
// fc2d_hps
/*************************************************/

static fc2d_hps_vtable_t s_hps_vt;

fc2d_hps_vtable_t* fc2d_hps_vt() {
    FCLAW_ASSERT(s_hps_vt.is_set != 0);
	return &s_hps_vt;
}

void fc2d_hps_solver_initialize(fclaw2d_global_t* glob) {
    int claw_version = 4; /* solution data is organized as (i,j,m) */
	fclaw2d_clawpatch_vtable_initialize(glob, claw_version);

	/* Patch : These could be over-written by user specific settings */
	fclaw2d_patch_vtable_t*   patch_vt = fclaw2d_patch_vt(glob);  
	patch_vt->rhs            = hps_rhs;   /* Calls FORTRAN routine */
    patch_vt->initialize     = hps_rhs;   /* Get an initial refinement */
	patch_vt->setup          = NULL;

    /* Tagging functions : Base refinement on the right hand side */
    patch_vt->tag4refinement = hps_tag4refinement;
    patch_vt->tag4coarsening = hps_tag4coarsening;

    /* Clawpatch and ForestClaw : Output functions */
    fclaw2d_vtable_t*   fclaw_vt = fclaw2d_vt(glob);
    fclaw_vt->output_frame = hps_output;
    fclaw_vt->hook_regrid = hps_regrid_hook;

    /* Elliptic specific functions */
    fclaw2d_elliptic_vtable_t *elliptic_vt = fclaw2d_elliptic_vt(glob);
    elliptic_vt->setup = hps_setup_solver;
    elliptic_vt->factor = hps_factor_solver;

    /* Solver doesn't do anything so far */
    elliptic_vt->solve = hps_solve;    
    elliptic_vt->apply_bc = fc2d_hps_physical_bc;

    /* BCs : Homogeneous BCs by default */
	fc2d_hps_vtable_t*  hps_vt = hps_vt_init();	
    // hps_vt->fort_apply_bc = &FC2D_HPS_FORT_APPLY_BC_DEFAULT;
    // hps_vt->fort_eval_bc  = &FC2D_HPS_FORT_EVAL_BC_DEFAULT;

    /* Diagnostics : Error, conservation */
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt(glob);
    clawpatch_vt->compute_error = hps_compute_error;  /* calls user-defined fortran routine */

    /* Conservation check : Compares sum(rhs) with sum of normal fluxes around the boundary
       of the solution.   (uses divergence theorem) */
    clawpatch_vt->conservation_check = hps_conservation_check;        

    /* These are specialized for the elliptic problem */
    fclaw2d_diagnostics_vtable_t *diag_vt = fclaw2d_diagnostics_vt(glob);
    diag_vt->patch_init_diagnostics     = fc2d_hps_diagnostics_initialize;
    diag_vt->patch_reset_diagnostics    = fc2d_hps_diagnostics_reset;
    diag_vt->patch_compute_diagnostics  = fc2d_hps_diagnostics_compute;
    diag_vt->patch_gather_diagnostics   = fc2d_hps_diagnostics_gather;
    diag_vt->patch_finalize_diagnostics = fc2d_hps_diagnostics_finalize;

	hps_vt->is_set = 1;
}

void fc2d_hps_setprob(fclaw2d_global_t* glob) {
    return;
}

void fc2d_hps_rhs(fclaw2d_global_t* glob, fclaw2d_patch_t *patch, int blockno, int patchno) {
    return;
}

static double lambdaGlobal = 0.0;

void fc2d_hps_heat_set_lambda(fclaw2d_global_t* glob, double lambda) {
    // EllipticForest::EllipticForestApp& app = EllipticForest::EllipticForestApp::getInstance();
    // app.options.setOption("lambda", lambda);
    lambdaGlobal = lambda;
    return;
}

double fc2d_hps_heat_get_lambda() {
    return lambdaGlobal;
}

void hps_setup_solver(fclaw2d_global_t *glob) {
    // Get EllipticForest app
    EllipticForest::EllipticForestApp& app = EllipticForest::EllipticForestApp::getInstance();
    app.log("Setting up EllipticForest...");
    app.options["homogeneous-rhs"] = false;

    // Get glob from user and get options
	fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
	fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);

    // Create root patch
    int nx = clawpatch_opt->mx;
    int ny = clawpatch_opt->my;
    double x_lower = fclaw_opt->ax;
    double x_upper = fclaw_opt->bx;
    double y_lower = fclaw_opt->ay;
    double y_upper = fclaw_opt->by;
    EllipticForest::FISHPACK::FISHPACKFVGrid root_grid(nx, ny, x_lower, x_upper, y_lower, y_upper);
    EllipticForest::FISHPACK::FISHPACKPatch root_patch(root_grid);
    root_patch.level = 0;
    root_patch.isLeaf = true;

    // Create patch solver
    double lambda = fc2d_hps_heat_get_lambda();
    EllipticForest::FISHPACK::FISHPACKFVSolver solver(lambda);

    // Create new HPS algorithm
    // TODO: delete HPS in clean up function (where...?)
    HPSAlgorithm* HPS = new HPSAlgorithm(root_patch, solver);

    // Save HPS into ForestClaw glob
    // TODO: Should I put this somewhere else?
    glob->user = (HPSAlgorithm*) HPS;

    // Call setup stage
    fclaw2d_domain_t* domain = glob->domain;
    p4est_wrap_t* p4est_wrap = (p4est_wrap_t*) domain->pp;
    p4est_t* p4est = p4est_wrap->p4est;
    HPS->setupStage(p4est);
}

void hps_factor_solver(fclaw2d_global_t* glob) {
    HPSAlgorithm* HPS = (HPSAlgorithm*) glob->user;
    // HPS->patchSolver = EllipticForest::FISHPACK::FISHPACKFVSolver(fc2d_hps_heat_get_lambda());
    HPS->buildStage();
}

void hps_rhs(fclaw2d_global_t *glob, fclaw2d_patch_t *patch, int blockno, int patchno) {
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

void hps_solve(fclaw2d_global_t *glob) {
    // Get EllipticForest app
    EllipticForest::EllipticForestApp& app = EllipticForest::EllipticForestApp::getInstance();
    app.log("Beginning HPS solve...");

    // Get fc2d_hps virtual table
    fc2d_hps_vtable_t* hps_vt = fc2d_hps_vt();
    fc2d_hps_options_t hps_opt = fc2d_hps_get_options(glob);

    // Get HPS algorithm from glob
    // TODO: Should I get this from somewhere else?
    HPSAlgorithm* HPS = (HPSAlgorithm*) glob->user;

    //
    // HPS->quadtree.merge([&](EllipticForest::FISHPACK::FISHPACKPatch& tau, EllipticForest::FISHPACK::FISHPACKPatch& alpha, EllipticForest::FISHPACK::FISHPACKPatch& beta, EllipticForest::FISHPACK::FISHPACKPatch& gamma, EllipticForest::FISHPACK::FISHPACKPatch& omega){

    //     HPS->coarsen_(tau, alpha, beta, gamma, omega);
    //     // HPS->createIndexSets_(tau, alpha, beta, gamma, omega);
    //     // HPS->createMatrixBlocks_(tau, alpha, beta, gamma, omega);

    //     // EllipticForest::Matrix<double>& T_alpha = alpha.matrixT();
    //     // EllipticForest::Matrix<double>& T_beta = beta.matrixT();
    //     // EllipticForest::Matrix<double>& T_gamma = gamma.matrixT();
    //     // EllipticForest::Matrix<double>& T_omega = omega.matrixT();

    //     // int n_alpha = T_alpha.nRows();
    //     // int n_beta = T_beta.nRows();
    //     // int n_gamma = T_gamma.nRows();
    //     // int n_omega = T_omega.nRows();

    //     // app.log("Pre-solve merge:");
    //     // app.log("\tn_alpha = %i, n_beta = %i, n_gamma = %i, n_omega = %i", n_alpha, n_beta, n_gamma, n_omega);



    // });

    // Update solver with lambda
    // TODO: There's gotta be a better way to do this...
    // app.log("Setting lambda to: %e", fc2d_hps_heat_get_lambda());
    HPS->patchSolver = EllipticForest::FISHPACK::FISHPACKFVSolver(fc2d_hps_heat_get_lambda());

    // Call build stage
    HPS->buildStage();

    // Call upwards stage
    HPS->upwardsStage([&](EllipticForest::FISHPACK::FISHPACKPatch& leafPatch){
        EllipticForest::FISHPACK::FISHPACKFVGrid& grid = leafPatch.grid();
        fclaw2d_patch_t* fc_patch = &(glob->domain->blocks->patches[leafPatch.leafID]);

        int mx, my;
        int mbc;
        double x_lower, y_lower, dx, dy;
        int mfields;
        double* rhs;
        fclaw2d_clawpatch_grid_data(glob, fc_patch, &mx, &my, &mbc, &x_lower, &y_lower, &dx, &dy);
        fclaw2d_clawpatch_rhs_data(glob, fc_patch, &rhs, &mfields);
        int nx = mx + 2*mbc;
        int ny = mx + 2*mbc;

        leafPatch.vectorF() = EllipticForest::Vector<double>(grid.nPointsX() * grid.nPointsY());
        for (auto i = 0; i < nx; i++) {
            for (auto j = 0; j < ny; j++) {
                if (i > mbc-1 && i < nx-mbc && j > mbc-1 && j < ny-mbc) {
                    int ii = i - mbc;
                    int jj = j - mbc;
                    int idx_ef = jj + ii*mx;
                    int idx_ef_T = ii + jj*my;
                    int idx_fc = j + i*nx;
                    int idx_fc_T = i + j*ny;

                    // q[idx_fc] = patch.vectorU()[idx_ef];
                    // rhs[idx_fc] = patch.vectorU()[idx_ef];
                    leafPatch.vectorF()[idx_ef] = rhs[idx_fc];
                }
                // int idx = j + i*grid.nPointsY();
                // int idx_T = i + j*grid.nPointsX();
                // leafPatch.vectorF()[idx] = rhs[idx];
                // printf("idx = %i, idx_T = %i, rhs = %f\n", idx, idx_T, rhs[idx_T]);
            }
        }
        return;
    });

    // Call solve stage; provide Dirichlet data via function
    std::vector<int> boundaryTypes = {
        std::get<int>(app.options["boundary-type-west"]),
        std::get<int>(app.options["boundary-type-east"]),
        std::get<int>(app.options["boundary-type-south"]),
        std::get<int>(app.options["boundary-type-north"])
    };
    if (!hps_opt->use_ext_bc) {
        HPS->solveStage([&](int side, double x, double y, double* a, double* b){
            double time = glob->curr_time;
            if (boundaryTypes[side] == 0) {
                *a = 1.0;
                *b = 0.0;
            }
            else if (boundaryTypes[side] == 1) {
                *a = 0.0;
                *b = 1.0;
            }
            else {
                std::cerr << "INVALID BC TYPE" << std::endl;
            }
            return hps_vt->fort_eval_bc(&side, &time, &x, &y);
        });
    }
    else {
        // Can I create an interpolation function by iterating through forestclaw patches, then call that interpolation function for the actual solveStage call?
        // Use Numerical Recipies algorithms!
        // NOTE: Need to determine order of interpolation
        HPS->solveStage([&](EllipticForest::FISHPACK::FISHPACKPatch& rootPatch){
            // Iterate through patches; get boundary intersections; call extended boundary function
            fclaw2d_global_iterate_patches(
                fclaw2d_global_t* glob,
                [](fclaw2d_domain_t* domain, fclaw2d_patch_t* patch, int blockno, int patchno, void* user){
                    // 


                },
                NULL);
            
        });
    }

    // HPS->solveStage([&](EllipticForest::FISHPACK::FISHPACKPatch& rootPatch){
    //     fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
    //     fc2d_hps_vtable_t* hps_vt = fc2d_hps_vt();

    //     EllipticForest::FISHPACK::FISHPACKFVGrid& grid = rootPatch.grid();
    //     rootPatch.vectorG() = EllipticForest::Vector<double>(2*grid.nPointsX() + 2*grid.nPointsY());

    //     // std::vector<int> boundaryTypes = {
    //     //     std::get<int>(app.options["boundary-type-west"]),
    //     //     std::get<int>(app.options["boundary-type-east"]),
    //     //     std::get<int>(app.options["boundary-type-south"]),
    //     //     std::get<int>(app.options["boundary-type-north"])
    //     // };

    //     // for (auto s = 0; s < 4; s++) {
    //     //     double x, y;
    //     //     int nPointsSide;

    //     //     if (s == 0 || s == 1)   nPointsSide = grid.nPointsY();
    //     //     else                    nPointsSide = grid.nPointsX();

    //     //     switch (s) {
    //     //         case 0:
    //     //             x = grid.xLower();
    //     //             break;
    //     //         case 1:
    //     //             x = grid.xUpper();
    //     //             break;
    //     //         case 2:
    //     //             y = grid.yLower();
    //     //             break;
    //     //         case 3:
    //     //             y = grid.yUpper();
    //     //             break;
    //     //     }

    //     //     EllipticForest::Vector<double> boundaryDataSide(nPointsSide);
    //     //     for (auto i = 0; i < nPointsSide; i++) {
    //     //         if (s == 0 || s == 1)   y = grid(1, i);
    //     //         else                    x = grid(0, i);
    //     //         boundaryDataSide[i] = hps_vt->fort_eval_bc(&boundaryTypes[s], &glob->curr_time, &x, &y);
    //     //     }

    //     //     if (boundaryTypes[s] == 0) {
    //     //         // Dirichlet BC
    //     //         rootPatch.vectorG().setSegment(s*nPointsSide, boundaryDataSide);
    //     //     }
    //     //     else if (boundaryTypes[s] == 1) {
    //     //         // Neumann BC
    //     //         EllipticForest::Matrix<double> T = rootPatch.matrixT()(s*nPointsSide, (s+1)*nPointsSide-1, s*nPointsSide, (s+1)*nPointsSide-1);
    //     //         EllipticForest::Vector<double> g = EllipticForest::solve(T, boundaryDataSide);
    //     //         rootPatch.vectorG().setSegment(s*nPointsSide, g);
    //     //     }

    //     // }

    //     EllipticForest::Vector<double> gWest(grid.nPointsY());
    //     EllipticForest::Vector<double> gEast(grid.nPointsY());
    //     EllipticForest::Vector<double> gSouth(grid.nPointsX());
    //     EllipticForest::Vector<double> gNorth(grid.nPointsX());

    //     int dirichletBC = 1;

    //     for (auto j = 0; j < grid.nPointsY(); j++) {
    //         double y = grid(1, j);
    //         double x_lower = grid.xLower();
    //         double x_upper = grid.xUpper();
    //         gWest[j] = hps_vt->fort_eval_bc(&dirichletBC, &glob->curr_time, &x_lower, &y);
    //         gEast[j] = hps_vt->fort_eval_bc(&dirichletBC, &glob->curr_time, &x_upper, &y);
    //     }
    //     for (auto i = 0; i < grid.nPointsX(); i++) {
    //         double x = grid(0, i);
    //         double y_lower = grid.yLower();
    //         double y_upper = grid.yUpper();
    //         gSouth[i] = hps_vt->fort_eval_bc(&dirichletBC, &glob->curr_time, &x, &y_lower);
    //         gNorth[i] = hps_vt->fort_eval_bc(&dirichletBC, &glob->curr_time, &x, &y_upper);
    //     }

    //     rootPatch.vectorG().setSegment(0*grid.nPointsX(), gWest);
    //     rootPatch.vectorG().setSegment(1*grid.nPointsX(), gEast);
    //     rootPatch.vectorG().setSegment(2*grid.nPointsX(), gSouth);
    //     rootPatch.vectorG().setSegment(3*grid.nPointsX(), gNorth);
        
    //     // std::cout << rootPatch.str() << std::endl;
    //     // std::cout << "g = " << rootPatch.vectorG() << std::endl;

    //     return;
    // });

    // Copy data to ForestClaw patch
    HPS->quadtree.traversePreOrder([&](EllipticForest::FISHPACK::FISHPACKPatch& patch){
        if (patch.isLeaf) {
            fclaw2d_patch_t* fc_patch = &(glob->domain->blocks->patches[patch.leafID]);

            int mbc;
            int mx, my;
            double x_lower, y_lower, dx, dy;
            double* q;
            int meqn, mfields;
            double* rhs;
            fclaw2d_clawpatch_grid_data(glob, fc_patch, &mx, &my, &mbc, &x_lower, &y_lower, &dx, &dy);
            fclaw2d_clawpatch_soln_data(glob, fc_patch, &q, &meqn);
            fclaw2d_clawpatch_rhs_data(glob, fc_patch, &rhs, &mfields);

            EllipticForest::FISHPACK::FISHPACKFVGrid& grid = patch.grid();
            int nx = mx + 2*mbc;
            int ny = mx + 2*mbc;
            for (auto i = 0; i < nx; i++) {
                for (auto j = 0; j < ny; j++) {
                    if (i > mbc-1 && i < nx-mbc && j > mbc-1 && j < ny-mbc) {
                        int ii = i - mbc;
                        int jj = j - mbc;
                        int idx_ef = jj + ii*mx;
                        int idx_ef_T = ii + jj*my;
                        int idx_fc = j + i*nx;
                        int idx_fc_T = i + j*ny;

                        // q[idx_fc] = patch.vectorU()[idx_ef];
                        rhs[idx_fc] = patch.vectorU()[idx_ef];
                    }
                }
            }
            // std::cout << patch.str() << std::endl;
            // std::cout << "g = " << patch.vectorG() << std::endl;
            // std::cout << "u = " << patch.vectorU() << std::endl;
            // std::cout << "f = " << patch.vectorF() << std::endl;
        }
        return;
    });

    return;
}

void hps_output(fclaw2d_global_t *glob, int iframe) {
    const fc2d_hps_options_t* hps_opt;
	hps_opt = fc2d_hps_get_options(glob);

	// if (hps_opt->ascii_out != 0)
	// 	fclaw2d_clawpatch_output_ascii(glob,iframe);

	// if (hps_opt->vtk_out != 0)
	// 	fclaw2d_clawpatch_output_vtk(glob,iframe);

    fclaw2d_clawpatch_output_vtk(glob, iframe);
}

int hps_tag4refinement(fclaw2d_global_t *glob, fclaw2d_patch_t *this_patch, int blockno, int patchno, int initflag) {
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

    fclaw2d_clawpatch_vtable_t* clawpatch_vt = fclaw2d_clawpatch_vt(glob);
    FCLAW_ASSERT(clawpatch_vt->fort_tag4refinement != NULL);


    /* Use default fortran tagging routines.  Choose refinement based on criteria
       set in configuration files (clawpatch:refinement-criteria) */
    tag_patch = 0;
    clawpatch_vt->fort_tag4refinement(&mx,&my,&mbc,&mfields,&xlower,&ylower,&dx,&dy,
                                      &blockno, rhs, &refine_threshold,
                                      &initflag, &tag_patch);
    return tag_patch;
}

int hps_tag4coarsening(fclaw2d_global_t *glob, fclaw2d_patch_t *fine_patches, int blockno, int patchno, int initflag) {
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

    fclaw2d_clawpatch_vtable_t* clawpatch_vt = fclaw2d_clawpatch_vt(glob);
    FCLAW_ASSERT(clawpatch_vt->fort_tag4coarsening != NULL);

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    double coarsen_threshold = fclaw_opt->coarsen_threshold;

    int tag_patch = 0;
    clawpatch_vt->fort_tag4coarsening(&mx,&my,&mbc,&mfields,&xlower,&ylower,&dx,&dy,
                                      &blockno, rhs[0],rhs[1],rhs[2],rhs[3],
                                      &coarsen_threshold,&initflag,&tag_patch);
    return tag_patch == 1;
}

void hps_compute_error(fclaw2d_global_t *glob, fclaw2d_patch_t *patch, int blockno, int patchno, void *user) {
    fc2d_hps_error_info_t* error_data = (fc2d_hps_error_info_t*) user;

    fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    if (fclaw_opt->compute_error)
    {
        fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt(glob);
        FCLAW_ASSERT(clawpatch_vt->fort_compute_patch_error != NULL);

        int mx, my, mbc;
        double xlower,ylower,dx,dy;
        fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,&xlower,
                                    &ylower,&dx,&dy);

        double *area = fclaw2d_clawpatch_get_area(glob,patch);  /* Might be null */

        /* Computing solution is stored in the RHS; true solution is stored in soln */
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

void hps_conservation_check(fclaw2d_global_t *glob, fclaw2d_patch_t *patch, int blockno, int patchno, void *user) {
    fc2d_hps_error_info_t* error_data = (fc2d_hps_error_info_t*) user;
    int mx, my, mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int mfields;
    double *rhs;  /* Solution is stored in the right hand side */ 
    fclaw2d_clawpatch_rhs_data(glob,patch,&rhs,&mfields);

    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt(glob);
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

fc2d_hps_vtable_t* hps_vt_init() {
    FCLAW_ASSERT(s_hps_vt.is_set == 0);
	return &s_hps_vt;
}

void hps_regrid_hook(fclaw2d_domain_t * old_domain,
                     fclaw2d_patch_t * old_patch,
                     fclaw2d_domain_t * new_domain,
                     fclaw2d_patch_t * new_patch,
                     fclaw2d_patch_relation_t newsize,
                     int blockno,
                     int old_patchno,
                     int new_patchno,
                     void *user) {

    // Get the EF app, HPS, and quadtree
    EllipticForest::EllipticForestApp& app = EllipticForest::EllipticForestApp::getInstance();
    // app.log("REGRID HOOK");
    // app.log("  OLD_PATCH # = %i", old_patchno);
    // app.log("  NEW_PATCH # = %i", new_patchno);
    fclaw2d_global_t* glob = fclaw2d_global_get_global();
    HPSAlgorithm* HPS = (HPSAlgorithm*) glob->user;
    auto& quadtree = HPS->quadtree;
    auto& solver = HPS->patchSolver;

    // Get current quadtree node
    int target_globalID = 0;
    quadtree.traversePreOrder([&](EllipticForest::Quadtree<EllipticForest::FISHPACK::FISHPACKPatch>::QuadtreeNode node){
        auto& patch = node.data;
        if (newsize == FCLAW2D_PATCH_HALFSIZE) {
            if (node.leafID == new_patchno) {
                // app.log("\tFOUND NODE");
                target_globalID = node.globalID;
                return false;
            }
        }
        else if (newsize == FCLAW2D_PATCH_DOUBLESIZE) {
            if (node.leafID == new_patchno) {
                // app.log("\tFOUND NODE");
                target_globalID = node.parentID;
                return false;
            }
        }
        return true;
    });
    auto target_node = quadtree.getNode(target_globalID);

    if (newsize == FCLAW2D_PATCH_SAMESIZE) {
        // Old patch has same refinement as new patch; copy data
        // app.log("\tNO REGRID");
    }
    else if (newsize == FCLAW2D_PATCH_HALFSIZE) {
        // Old patch is more coarse than new patch; interpolate 1 patch to 4
        // app.log("\tINTERPOLATE 1 TO 4");
        quadtree.refineNode(target_node.globalID, [&](EllipticForest::FISHPACK::FISHPACKPatch& parentPatch){
            std::vector<EllipticForest::FISHPACK::FISHPACKPatch> childrenPatches(4);

            // Build children grids
            auto& parentGrid = parentPatch.grid();
            int nx = parentGrid.nPointsX();
            int ny = parentGrid.nPointsY();
            double xLower = parentGrid.xLower();
            double xUpper = parentGrid.xUpper();
            double xMid = (xLower + xUpper) / 2.0;
            double yLower = parentGrid.yLower();
            double yUpper = parentGrid.yUpper();
            double yMid = (yLower + yUpper) / 2.0;
            EllipticForest::FISHPACK::FISHPACKFVGrid child0Grid(nx, ny, xLower, xMid, yLower, yMid);
            EllipticForest::FISHPACK::FISHPACKFVGrid child1Grid(nx, ny, xMid, xUpper, yLower, yMid);
            EllipticForest::FISHPACK::FISHPACKFVGrid child2Grid(nx, ny, xLower, xMid, yMid, yUpper);
            EllipticForest::FISHPACK::FISHPACKFVGrid child3Grid(nx, ny, xMid, xUpper, yMid, yUpper);

            // Build child patches
            childrenPatches[0] = EllipticForest::FISHPACK::FISHPACKPatch(child0Grid);
            childrenPatches[1] = EllipticForest::FISHPACK::FISHPACKPatch(child1Grid);
            childrenPatches[2] = EllipticForest::FISHPACK::FISHPACKPatch(child2Grid);
            childrenPatches[3] = EllipticForest::FISHPACK::FISHPACKPatch(child3Grid);
            int globalIDCounter = 1;
            int leafIDCounter = 0;
            for (auto& patch : childrenPatches) {
                patch.globalID = parentPatch.globalID + globalIDCounter++; // ??? Maybe...
                patch.leafID = parentPatch.leafID + leafIDCounter++; // ??? Maybe...
                patch.level = parentPatch.level + 1;
                patch.isLeaf = true;
                patch.nCoarsens = 0;
            }
            parentPatch.leafID = -1;
            // parentPatch.nCoarsens = 0;

            // Compute child data
            for (auto& patch : childrenPatches) {
                // Compute T
                patch.matrixT() = solver.buildD2N(patch.grid());

                // Compute f
                auto& grid = patch.grid();
                fclaw2d_patch_t* fc_patch = &(new_domain->blocks->patches[patch.leafID]);

                int mx,my,mbc;
                double dx,dy,xlower,ylower;
                fclaw2d_clawpatch_grid_data(glob,fc_patch,&mx,&my,&mbc,
                                            &xlower,&ylower,&dx,&dy);

                int mfields;
                double *rhs;
                fclaw2d_clawpatch_rhs_data(glob,fc_patch,&rhs,&mfields);

                int meqn;
                double *q;
                fclaw2d_clawpatch_soln_data(glob,fc_patch,&q,&meqn);

                double lambda = fc2d_hps_heat_get_lambda();

                

                // int mfields;
                // double* rhs;
                // global_heat_rhs(glob, fc_patch, blockno, new_patchno);
                // NOTE: Need to call heat_rhs (heat_user.cpp) in order to init this data
                // fclaw2d_clawpatch_rhs_data(glob, fc_patch, &rhs, &mfields);
                // NOTE: Check if rhs is a null pointer
                int nx = mx + 2*mbc;
                int ny = mx + 2*mbc;
                patch.vectorF() = EllipticForest::Vector<double>(grid.nPointsX() * grid.nPointsY());
                for (auto i = 0; i < nx; i++) {
                    for (auto j = 0; j < ny; j++) {
                        if (i > mbc-1 && i < nx-mbc && j > mbc-1 && j < ny-mbc) {
                            int ii = i - mbc;
                            int jj = j - mbc;
                            int idx_ef = jj + ii*mx;
                            int idx_ef_T = ii + jj*my;
                            int idx_fc = j + i*nx;
                            int idx_fc_T = i + j*ny;

                            // rhs[idx_fc] = lambda*q[idx_fc];
                            // printf("rhs[%i] = %10.4f\n", idx_fc, rhs[idx_fc]);

                            // patch.vectorF()[idx_ef] = rhs[idx_fc];
                            patch.vectorF()[idx_ef] = lambda*q[idx_fc];

                            // q[idx_fc] = patch.vectorU()[idx_ef];
                            // rhs[idx_fc] = patch.vectorU()[idx_ef];
                        }
                        // int idx = j + i*grid.nPointsY();
                        // int idx_T = i + j*grid.nPointsX();
                        // patch.vectorF()[idx] = rhs[idx];
                    }
                }

                // Compute h
                patch.vectorH() = solver.particularNeumannData(patch.grid(), patch.vectorF());
            }

            // Merge child data to get new parent data; updates parent node with appropiate data
            HPS->merge4to1(parentPatch, childrenPatches[0], childrenPatches[1], childrenPatches[2], childrenPatches[3]);
            HPS->upwards4to1(parentPatch, childrenPatches[0], childrenPatches[1], childrenPatches[2], childrenPatches[3]);
            // HPS->split1to4(parentPatch, childrenPatches[0], childrenPatches[1], childrenPatches[2], childrenPatches[3]);

            return childrenPatches;
        });
    }
    else if (newsize == FCLAW2D_PATCH_DOUBLESIZE) {
        // Old patch is more fine than new patch; average 4 patches to 1
        // app.log("\tAVERAGE 4 TO 1");
        quadtree.coarsenNode(target_node.globalID, [&](EllipticForest::FISHPACK::FISHPACKPatch& child0Patch, EllipticForest::FISHPACK::FISHPACKPatch& child1Patch, EllipticForest::FISHPACK::FISHPACKPatch& child2Patch, EllipticForest::FISHPACK::FISHPACKPatch& child3Patch){
            EllipticForest::FISHPACK::FISHPACKPatch parentPatch;

            // Build parent grid
            int nx = child0Patch.grid().nPointsX();
            int ny = child0Patch.grid().nPointsY();
            double xLower = child0Patch.grid().xLower();
            double xUpper = child1Patch.grid().xUpper();
            double yLower = child0Patch.grid().yLower();
            double yUpper = child2Patch.grid().yUpper();
            EllipticForest::FISHPACK::FISHPACKFVGrid parentGrid(nx, ny, xLower, xUpper, yLower, yUpper);
            parentPatch.grid() = parentGrid;
            parentPatch.globalID = target_node.globalID;
            parentPatch.leafID = child0Patch.leafID;
            parentPatch.level = target_node.level;
            parentPatch.isLeaf = true;

            // Compute parent data
            //    Compute T
            parentPatch.matrixT() = solver.buildD2N(parentPatch.grid());

            //    Compute f
            fclaw2d_patch_t* fc_patch = &(new_domain->blocks->patches[new_patchno]);
            int mfields;
            double* rhs;
            fclaw2d_clawpatch_rhs_data(glob, fc_patch, &rhs, &mfields);
            parentPatch.vectorF() = EllipticForest::Vector<double>(parentPatch.grid().nPointsX() * parentPatch.grid().nPointsY());
            for (auto i = 0; i < parentPatch.grid().nPointsX(); i++) {
                for (auto j = 0; j < parentPatch.grid().nPointsY(); j++) {
                    int idx = j + i*parentPatch.grid().nPointsY();
                    int idx_T = i + j*parentPatch.grid().nPointsX();
                    parentPatch.vectorF()[idx] = rhs[idx];
                }
            }

            //    Compute h
            parentPatch.vectorH() = solver.particularNeumannData(parentPatch.grid(), parentPatch.vectorF());

            return parentPatch;
        });
    }

    // For debugging: Check if quadtree is valid
    // if (!quadtree.isValid()) {
    //     std::cerr << "[EllipticForest] Quadtree is invalid." << std::endl;
    //     exit;
    // }

    return;

}

/*************************************************/
// fc2d_hps_physical_bc
/*************************************************/

void fc2d_hps_physical_bc(struct fclaw2d_global *glob) {
    fc2d_hps_time_info_t tinfo;
    tinfo.t = glob->curr_time;
    fclaw2d_global_iterate_patches(glob,
                                   cb_fc2d_hps_physical_bc,
                                   (void *) &tinfo);
}

void cb_fc2d_hps_physical_bc(fclaw2d_domain_t *domain, fclaw2d_patch_t *patch, int blockno, int patchno, void *user) {
    fclaw2d_global_iterate_t* s = (fclaw2d_global_iterate_t*) user;
    fc2d_hps_time_info_t *tinfo = (fc2d_hps_time_info_t*) s->user;

    double t = tinfo->t;


    /* Determine which faces are at the physical boundary */
    int intersects_bc[4];
    fclaw2d_physical_get_bc(s->glob,blockno,patchno,intersects_bc);

    int mx, my, mbc;
    double xlower, ylower, dx, dy;
    fclaw2d_clawpatch_grid_data(s->glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    const fc2d_hps_options_t* hps_opt = fc2d_hps_get_options(s->glob);
    int mfields;
    double *rhs;
    fclaw2d_clawpatch_rhs_data(s->glob,patch,&rhs,&mfields);

    fc2d_hps_vtable_t*  hps_vt = fc2d_hps_vt();

    int cons_check = 0;
    double flux_sum[4];

    hps_vt->fort_apply_bc(&blockno, &mx, &my, &mbc, &mfields, 
                         &xlower, &ylower, &dx,&dy,&t, intersects_bc,
                         hps_opt->boundary_conditions,rhs, hps_vt->fort_eval_bc,
                         &cons_check, flux_sum);
}

void fc2d_hps_physical_get_bc(fclaw2d_global_t *glob, int blockno, int patchno, int *intersects_bdry) {
    // const int numfaces = get_faces_per_patch(domain);
    int bdry[4];
    fclaw2d_patch_boundary_type(glob->domain,blockno,patchno,bdry);
    int i;
    for(i = 0; i < 4; i++)
    {
        // Physical boundary conditions
        intersects_bdry[i] = bdry[i] == 1;
    }
}

/*************************************************/
// fc2d_hps_output_ascii
/*************************************************/

void fc2d_hps_time_header_ascii(fclaw2d_global_t* glob, int iframe) {
    const fclaw2d_clawpatch_options_t *clawpatch_opt = 
    fclaw2d_clawpatch_get_options(glob);
    char matname1[20];
    sprintf(matname1,"fort.q%04d",iframe);

    FILE *f1 = fopen(matname1,"w");
    fclose(f1);

    char matname2[20];
    sprintf(matname2,"fort.t%04d",iframe);

    double time = glob->curr_time;

    int ngrids = glob->domain->global_num_patches;

    int mfields = clawpatch_opt->rhs_fields;  
    int maux = clawpatch_opt->maux;

    int mf;
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    if (fclaw_opt->compute_error)
        mf = mfields + 2;  /* Print out error and true solution */
    else
        mf = mfields;  /* Only print out computed solution */


    FILE *f2 = fopen(matname2,"w");
    fprintf(f2,"%12.6f %23s\n%5d %30s\n%5d %30s\n%5d %30s\n%5d %30s\n",time,"time",
            mf,"mfields",ngrids,"ngrids",maux,"num_aux",2,"num_dim");
    fclose(f2);
}

void cb_hps_output_ascii(fclaw2d_domain_t * domain, fclaw2d_patch_t * patch, int blockno, int patchno, void *user) {
    fclaw2d_global_iterate_t* s = (fclaw2d_global_iterate_t*) user;
    fclaw2d_global_t *glob = (fclaw2d_global_t*) s->glob;
    int iframe = *((int *) s->user);

    /* Get info not readily available to user */
    int global_num, local_num;
    int level;
    fclaw2d_patch_get_info(glob->domain,patch,
                           blockno,patchno,
                           &global_num,&local_num, &level);
    
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *rhs;
    int mfields;
    fclaw2d_clawpatch_rhs_data(glob,patch,&rhs,&mfields);

    double *err;
    fclaw2d_clawpatch_elliptic_error_data(glob,patch,&err,&mfields);

    double *soln;
    fclaw2d_clawpatch_elliptic_soln_data(glob,patch,&soln,&mfields);

    char fname[BUFSIZ];
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    snprintf (fname, BUFSIZ, "%s.q%04d", fclaw_opt->prefix, iframe);

    fc2d_hps_vtable_t*  hps_vt = fc2d_hps_vt(); 
    FCLAW_ASSERT(hps_vt->fort_output != NULL);

    if (fclaw_opt->compute_error)
    {        
        hps_vt->fort_output(fname,&mx,&my,&mfields,&mbc,
                            &xlower,&ylower,&dx,&dy,rhs,
                            soln, err, &global_num, &level,&blockno,
                            &glob->mpirank);
    }
    else
    {
        fclaw2d_clawpatch_vtable_t* clawpatch_vt = fclaw2d_clawpatch_vt(glob);
        clawpatch_vt->fort_output_ascii(fname,&mx,&my,&mfields,&mbc,
                                        &xlower,&ylower,&dx,&dy,rhs,
                                        &global_num,&level,&blockno,
                                        &glob->mpirank);        
    }
}

/*************************************************/
// fc2d_hps_options
/*************************************************/

static int s_hps_options_package_id = -1;
static const fclaw_app_options_vtable_t hps_options_vtable = {
    options_register,
    options_postprocess,
    options_check,
    options_destroy,
};

fc2d_hps_options_t*  fc2d_hps_options_register (fclaw_app_t * app, const char *configfile) {
    fc2d_hps_options_t *hps_opt;

    FCLAW_ASSERT (app != NULL);

    hps_opt = FCLAW_ALLOC (fc2d_hps_options_t, 1);
    fclaw_app_options_register (app, "hps", configfile,
                                &hps_options_vtable, hps_opt);
    
    fclaw_app_set_attribute(app,"hps",hps_opt);
    return hps_opt;
}

void fc2d_hps_package_register(fclaw_app_t* app, fc2d_hps_options_t* hps_opt) {
    return;
}

fc2d_hps_options_t* fc2d_hps_get_options(struct fclaw2d_global *glob) {
    int id = s_hps_options_package_id;
    return (fc2d_hps_options_t*) fclaw_package_get_options(glob,id);
}

void fc2d_hps_options_store (struct fclaw2d_global* glob, fc2d_hps_options_t* hps_opt) {
    int id = fclaw_package_container_add_pkg(glob,hps_opt);
    s_hps_options_package_id = id;
}

void* hps_register (fc2d_hps_options_t* hps_opt, sc_options_t * opt) {

    sc_options_add_bool(opt, 0, "use-ext-bc", &hps_opt->use_ext_bc, 0, "Option to use extended boundary condition FORTRAN function");

    /* Array of NumFaces=4 values */
    fclaw_options_add_int_array (opt, 0, "boundary_conditions", 
                                 &hps_opt->bc_cond_string, "1 1 1 1",
                                 &hps_opt->boundary_conditions, 4,
                                 "[hps] Physical boundary condition type [1 1 1 1]");

    sc_options_add_bool (opt, 0, "use-ext-bc", &hps_opt->use_ext_bc, 0, "Option to use extended version of boundary condition function (user supplied)");

    sc_options_add_bool (opt, 0, "ascii-out", &hps_opt->ascii_out, 0,
                           "Output ASCII formatted data [F]");

    sc_options_add_bool (opt, 0, "vtk-out", &hps_opt->vtk_out, 0,
                           "Output VTK formatted data [F]");

    // sc_options_add_bool(opt, 0, "mmio_out", &hps_opt->mmio_out, 0,
    //                         "Output matrices in Matrix Market format [F]");

    sc_options_add_bool(opt, 0, "cache_T", &hps_opt->cache_T, 1,
                            "Cache DtN matrix T by level [T]");

    sc_options_add_bool(opt, 0, "time_setup", &hps_opt->time_setup, 0,
                            "Time setup stage of HPS method [F]");

    sc_options_add_bool(opt, 0, "time_build", &hps_opt->time_build, 0,
                            "Time build stage of HPS method [F]");

    sc_options_add_bool(opt, 0, "time_upwards", &hps_opt->time_upwards, 0,
                            "Time upwards stage of HPS method [F]");

    sc_options_add_bool(opt, 0, "time_solve", &hps_opt->time_solve, 0,
                            "Time solve stage of HPS method [F]");

    // sc_options_add_bool(opt, 0, "time_copy", &hps_opt->time_copy, 0,
    //                         "Time copy stage of HPS method [F]");

    sc_options_add_bool(opt, 0, "nonhomogeneous_rhs", &hps_opt->nonhomogeneous_rhs, 0,
                            "Flag for non-homogeneous RHS; determines if upwards pass is done [F]");

    sc_options_add_bool(opt, 0, "only_patch_solver", &hps_opt->only_patch_solver, 0,
                            "Option to bypass HPS method and just use the patch solver [F]");


    /* Set operator type (laplace, varpoisson, heat, ...) */
    // sc_keyvalue_t *kv_op = hps_opt->kv_operator_type = sc_keyvalue_new ();
    // sc_keyvalue_set_int (kv_op, "laplace",  LAPLACE);     /* Uses FFT or BICG */
    // sc_keyvalue_set_int (kv_op, "varpoisson", VARPOISSON);   /* Uses BICG */
    // sc_keyvalue_set_int (kv_op, "heat",       HEAT);   /* Uses BICG */
    // sc_keyvalue_set_int (kv_op, "user_solver",  USER_SOLVER);   /* Uses BICG */
    // sc_options_add_keyvalue (opt, 0, "operator-type", &hps_opt->operator_type,
    //                          "laplace", kv_op, "Set operator type [laplace]");

    /* Set solver type (FFT, DST, BICG, ...) */
    // sc_keyvalue_t *kv_s = hps_opt->kv_patch_solver = sc_keyvalue_new ();
    // sc_keyvalue_set_int (kv_s, "bicg", BICG);
    // sc_keyvalue_set_int (kv_s, "fishpack", FISHPACK);
    // sc_keyvalue_set_int (kv_s, "dst", DST);
    // sc_keyvalue_set_int (kv_s, "fft",  FFT);     
    // sc_keyvalue_set_int (kv_s, "user_solver",  USER_SOLVER);     
    // sc_options_add_keyvalue (opt, 0, "patch-solver", &hps_opt->patch_solver,
    //                          "fishpack", kv_s, "Set patch solver [fishpack]");

    hps_opt->is_registered = 1;
    return NULL;
}

fclaw_exit_type_t hps_postprocess (fc2d_hps_options_t * hps_opt) {
    fclaw_options_convert_int_array (hps_opt->bc_cond_string, 
                                     &hps_opt->boundary_conditions,4);
    
    return FCLAW_NOEXIT;
}

fclaw_exit_type_t hps_check(fc2d_hps_options_t *hps_opt, fclaw2d_clawpatch_options_t *clawpatch_opt) {
    return FCLAW_NOEXIT;
}

void hps_destroy (fc2d_hps_options_t * hps_opt) {
    fclaw_options_destroy_array (hps_opt->boundary_conditions);

    FCLAW_ASSERT (hps_opt->kv_patch_solver != NULL);
    sc_keyvalue_destroy (hps_opt->kv_patch_solver);

    FCLAW_ASSERT (hps_opt->kv_operator_type != NULL);
    sc_keyvalue_destroy (hps_opt->kv_operator_type);
}

void* options_register (fclaw_app_t * app, void *package, sc_options_t * opt) {
    fc2d_hps_options_t *hps_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);

    hps_opt = (fc2d_hps_options_t*) package;

    return hps_register(hps_opt,opt);
}

fclaw_exit_type_t options_postprocess (fclaw_app_t * app, void *package, void *registered) {
    fc2d_hps_options_t *hps_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    hps_opt = (fc2d_hps_options_t*) package;
    FCLAW_ASSERT (hps_opt->is_registered);

    return hps_postprocess (hps_opt);
}

fclaw_exit_type_t options_check (fclaw_app_t * app, void *package, void *registered) {
    fc2d_hps_options_t *hps_opt;
    fclaw2d_clawpatch_options_t *clawpatch_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    hps_opt = (fc2d_hps_options_t*) package;
    FCLAW_ASSERT (hps_opt->is_registered);

    clawpatch_opt = (fclaw2d_clawpatch_options_t *)
        fclaw_app_get_attribute(app,"clawpatch",NULL);
    FCLAW_ASSERT(clawpatch_opt->is_registered);

    return hps_check(hps_opt,clawpatch_opt);  
}

void options_destroy (fclaw_app_t * app, void *package, void *registered) {
    fc2d_hps_options_t *hps_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    hps_opt = (fc2d_hps_options_t*) package;
    FCLAW_ASSERT (hps_opt->is_registered);

    hps_destroy (hps_opt);

    FCLAW_FREE (hps_opt);
}

/*************************************************/
// fc2d_hps_diagnostics
/*************************************************/

void fc2d_hps_diagnostics_initialize(fclaw2d_global_t *glob, void **acc_patch) {
    const fclaw2d_clawpatch_options_t *clawpatch_opt = 
              fclaw2d_clawpatch_get_options(glob);

    fc2d_hps_error_info_t *error_data;

    int mfields = clawpatch_opt->rhs_fields;

    error_data = FCLAW_ALLOC(fc2d_hps_error_info_t,1);

    /* Allocate memory for 1-norm, 2-norm, and inf-norm errors */
    error_data->local_error  = FCLAW_ALLOC_ZERO(double,3*mfields);
    error_data->global_error  = FCLAW_ALLOC_ZERO(double,3*mfields);
    error_data->mass   = FCLAW_ALLOC_ZERO(double,mfields);
    error_data->mass0  = FCLAW_ALLOC_ZERO(double,mfields);
    error_data->rhs   = FCLAW_ALLOC_ZERO(double,mfields);
    error_data->boundary   = FCLAW_ALLOC_ZERO(double,mfields);
    error_data->area = 0;
    error_data->c_kahan = FCLAW_ALLOC_ZERO(double,mfields);   

    *acc_patch = error_data;
}

void fc2d_hps_diagnostics_reset(fclaw2d_global_t *glob, void* patch_acc) {
    fc2d_hps_error_info_t *error_data = (fc2d_hps_error_info_t*) patch_acc;
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    int mfields = clawpatch_opt->rhs_fields;

    for(int m = 0; m < mfields; m++)
    {
        int i1 = m;            /* 1-norm */
        int i2 = mfields + m;     /* 2-norm */
        int iinf = 2*mfields + m; /* inf-norm */
        error_data->local_error[i1] = 0;
        error_data->local_error[i2] = 0;
        error_data->local_error[iinf] = 0;
        error_data->mass[m]= 0;        

        error_data->rhs[m] = 0;
        error_data->boundary[m] = 0;
        error_data->c_kahan[m] = 0;
    }
    error_data->area = 0;
}

void fc2d_hps_compute(fclaw2d_domain_t *domain, fclaw2d_patch_t *patch, int blockno, int patchno, void* user) {
    fclaw2d_global_iterate_t *s = (fclaw2d_global_iterate_t *) user;
    fc2d_hps_error_info_t *error_data = (fc2d_hps_error_info_t*) s->user; 

    /* Accumulate area for final computation of error */
    int mx, my, mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(s->glob,patch,&mx,&my,&mbc,&xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt(s->glob);
    double *area = fclaw2d_clawpatch_get_area(s->glob,patch);  
    FCLAW_ASSERT(clawpatch_vt->fort_compute_patch_area != NULL);
    error_data->area += clawpatch_vt->fort_compute_patch_area(&mx,&my,&mbc,&dx,&dy,area);

    /* Compute error */
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(s->glob);
    if (fclaw_opt->compute_error)
    {
        clawpatch_vt->compute_error(s->glob, patch, blockno, patchno, error_data);
    }

    if (fclaw_opt->conservation_check)
    {
        clawpatch_vt->conservation_check(s->glob, patch, blockno, patchno, error_data);
    }
}

void fc2d_hps_diagnostics_compute(fclaw2d_global_t* glob, void* patch_acc) {
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    int check = fclaw_opt->compute_error || fclaw_opt->conservation_check;
    if (!check) return;

    fclaw2d_global_iterate_patches(glob, fc2d_hps_compute, patch_acc);
}

void fc2d_hps_diagnostics_gather(fclaw2d_global_t *glob, void* patch_acc, int init_flag) {
    fclaw2d_domain_t *domain = glob->domain;
    
    fc2d_hps_error_info_t *error_data = (fc2d_hps_error_info_t*) patch_acc;
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    
    int mfields = clawpatch_opt->rhs_fields;  /* clawpatch->meqn */

    if (fclaw_opt->compute_error != 0)
    {
        double total_area = fclaw2d_domain_global_sum(domain, error_data->area);
        FCLAW_ASSERT(total_area != 0);

        double *error_norm = FCLAW_ALLOC_ZERO(double,3*mfields);
        for (int m = 0; m < mfields; m++)
        {
            int i1 = m;            /* 1-norm */
            int i2 = mfields + m;     /* 2-norm */
            int i3 = 2*mfields + m; /* inf-norm */

            error_norm[i1]  = fclaw2d_domain_global_sum(domain, error_data->local_error[i1]);
            error_norm[i1] /= total_area;

            error_norm[i2]  = fclaw2d_domain_global_sum(domain, error_data->local_error[i2]);
            error_norm[i2] /= total_area;
            error_norm[i2] = sqrt(error_norm[i2]);

            error_norm[i3] = fclaw2d_domain_global_maximum(domain, 
                                                           error_data->local_error[i3]);

            error_data->global_error[i1] = error_norm[i1];
            error_data->global_error[i2] = error_norm[i2];
            error_data->global_error[i3] = error_norm[i3];

            fclaw_global_essentialf("error[%d] = %16.6e %16.6e %16.6e\n",m,
                                    error_norm[i1], error_norm[i2],error_norm[i3]);
        }
        FCLAW_FREE(error_norm);
    }


    if (fclaw_opt->conservation_check != 0)
    {
        double *total_mass = FCLAW_ALLOC_ZERO(double,mfields);
        for(int m = 0; m < mfields; m++)
        {
            /* Store mass for future checks */
            if (init_flag)
            {
                total_mass[m] = fclaw2d_domain_global_sum(domain, error_data->rhs[m]);
                error_data->mass0[m] = total_mass[m];                
            }
            else
            {
                total_mass[m] = fclaw2d_domain_global_sum(domain, error_data->boundary[m]);
            }
            fclaw_global_essentialf("sum[%d] =  %24.16e  %24.16e\n",m,total_mass[m],
                                    fabs(total_mass[m]-error_data->mass0[m]));
        }
        FCLAW_FREE(total_mass);        
    }
}

void fc2d_hps_diagnostics_finalize(fclaw2d_global_t *glob, void** patch_acc) {
    fc2d_hps_error_info_t *error_data = *((fc2d_hps_error_info_t**) patch_acc);
    FCLAW_FREE(error_data->mass);
    FCLAW_FREE(error_data->mass0);
    FCLAW_FREE(error_data->rhs);
    FCLAW_FREE(error_data->boundary);
    FCLAW_FREE(error_data->c_kahan);
    FCLAW_FREE(error_data->local_error);
    FCLAW_FREE(error_data->global_error);
    FCLAW_FREE(error_data);
    *patch_acc = NULL;
}

void fc2d_hps_compute_diagnostics(fclaw2d_domain_t *domain, fclaw2d_patch_t *patch, int blockno, int patchno, void* user) {
    return;
}