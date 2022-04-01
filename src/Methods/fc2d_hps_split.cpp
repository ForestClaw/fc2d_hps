#include <Methods/fc2d_hps_split.hpp>

void split_vertical(fc2d_hps_patch& tau, fc2d_hps_patch& alpha, fc2d_hps_patch& beta) {

    // Get options
    fclaw2d_global_t* glob = (fclaw2d_global_t*) tau.user;
    fc2d_hps_options_t* hps_opt = fc2d_hps_get_options(glob);

    // Create child grids
    // std::cout << "[split_vertical]  creating child grids" << std::endl;
    fc2d_hps_patchgrid grid_alpha(tau.grid.Nx, tau.grid.Ny/2, tau.grid.x_lower, tau.grid.x_upper, tau.grid.y_lower, (tau.grid.y_lower + tau.grid.y_upper)/2);
    fc2d_hps_patchgrid grid_beta(tau.grid.Nx, tau.grid.Ny/2, tau.grid.x_lower, tau.grid.x_upper, (tau.grid.y_lower + tau.grid.y_upper)/2, tau.grid.y_upper);
    alpha.grid = grid_alpha;
    beta.grid = grid_beta;

    // Apply solution operator to get interior data
    // std::cout << "[split_vertical]  applying solution operator" << std::endl;
    // printf("HERE1\n");
    fc2d_hps_vector<double> u_tau = tau.S * tau.g;
    if (hps_opt->nonhomogeneous_rhs) {
        // printf("HERE1a\n");
        u_tau = u_tau + tau.w;
    }

    // Set child patch data
    //    Extract WESN components of tau.g
    // std::cout << "[split_vertical]  extracting WESN components of tau" << std::endl;
    int N_tau = tau.grid.Nx;
    fc2d_hps_vector<double> g_W_tau = tau.g.extract(0*N_tau, N_tau);
    fc2d_hps_vector<double> g_E_tau = tau.g.extract(1*N_tau, N_tau);
    fc2d_hps_vector<double> g_S_tau = tau.g.extract(2*N_tau, N_tau);
    fc2d_hps_vector<double> g_N_tau = tau.g.extract(3*N_tau, N_tau);

    //    Create WESN components of alpha.g
    // std::cout << "[split_vertical]  creating WESN components of alpha" << std::endl;
    fc2d_hps_vector<double> g_W_alpha = g_W_tau.extract(0, tau.grid.Ny / 2);
    fc2d_hps_vector<double> g_E_alpha = g_E_tau.extract(0, tau.grid.Ny / 2);
    // g_S_alpha = g_S_tau
    // g_N_alpha = u_tau

    //    Create WESN components of beta.g
    // std::cout << "[split_vertical]  creating WESN components of beta" << std::endl;
    fc2d_hps_vector<double> g_W_beta = g_W_tau.extract(tau.grid.Ny / 2, tau.grid.Ny / 2);
    fc2d_hps_vector<double> g_E_beta = g_E_tau.extract(tau.grid.Ny / 2, tau.grid.Ny / 2);
    // g_S_beta = u_tau
    // g_N_beta = g_N_tau

    //    Intract into alpha.g
    // std::cout << "[split_vertical]  creating g_alpha" << std::endl;
    fc2d_hps_vector<double> g_alpha(2*alpha.grid.Nx + 2*alpha.grid.Ny);
    g_alpha.intract(0, g_W_alpha);
    g_alpha.intract(alpha.grid.Ny, g_E_alpha);
    g_alpha.intract(2*alpha.grid.Ny, g_S_tau);
    g_alpha.intract(2*alpha.grid.Ny + alpha.grid.Nx, u_tau);

    //    Intract into beta.g
    // std::cout << "[split_vertical]  creating g_beta" << std::endl;
    fc2d_hps_vector<double> g_beta(2*beta.grid.Nx + 2*beta.grid.Ny);
    g_beta.intract(0, g_W_beta);
    g_beta.intract(beta.grid.Ny, g_E_beta);
    g_beta.intract(2*beta.grid.Ny, u_tau);
    g_beta.intract(2*beta.grid.Ny + beta.grid.Nx, g_N_tau);

    // Set child patch data
    // std::cout << "[split_vertical]  assigning child patch data" << std::endl;
    alpha.g = g_alpha;
    beta.g = g_beta;

    // std::cout << "[split_vertical]  returning" << std::endl;
    return;
}

void split_horizontal(fc2d_hps_patch& tau, fc2d_hps_patch& alpha, fc2d_hps_patch& beta) {

    // Get options
    fclaw2d_global_t* glob = (fclaw2d_global_t*) tau.user;
    fc2d_hps_options_t* hps_opt = fc2d_hps_get_options(glob);
    
    // Create child grids
    // std::cout << "[split_horizontal]  creating child grids" << std::endl;
    fc2d_hps_patchgrid grid_alpha(tau.grid.Nx/2, tau.grid.Ny, tau.grid.x_lower, (tau.grid.x_lower + tau.grid.x_upper)/2, tau.grid.y_lower, tau.grid.y_upper);
    fc2d_hps_patchgrid grid_beta(tau.grid.Nx/2, tau.grid.Ny, (tau.grid.x_lower + tau.grid.x_upper)/2, tau.grid.x_upper, tau.grid.y_lower, tau.grid.y_upper);
    alpha.grid = grid_alpha;
    beta.grid = grid_beta;

    // Apply solution operator to get interior data
    // std::cout << "[split_horizontal]  applying solution operator" << std::endl;
    // printf("HERE2\n");
    fc2d_hps_vector<double> u_tau = tau.S * tau.g;
    if (hps_opt->nonhomogeneous_rhs) {
        // printf("HERE2a\n");
        u_tau = u_tau + tau.w;
    }

    // Set child patch data
    //    Extract WESN components of tau.g
    // std::cout << "[split_horizontal]  extracting WESN components of tau" << std::endl;
    fc2d_hps_vector<double> g_W_tau = tau.g.extract(0, tau.grid.Ny);
    fc2d_hps_vector<double> g_E_tau = tau.g.extract(tau.grid.Ny, tau.grid.Ny);
    fc2d_hps_vector<double> g_S_tau = tau.g.extract(2*tau.grid.Ny, tau.grid.Nx);
    fc2d_hps_vector<double> g_N_tau = tau.g.extract(2*tau.grid.Ny + tau.grid.Nx, tau.grid.Nx);

    //    Create WESN components of alpha.g
    // std::cout << "[split_horizontal]  creating WESN components of alpha" << std::endl;
    // g_W_alpha = g_W_tau
    // g_E_alpha = u_tau
    fc2d_hps_vector<double> g_S_alpha = g_S_tau.extract(0, tau.grid.Nx / 2);
    fc2d_hps_vector<double> g_N_alpha = g_N_tau.extract(0, tau.grid.Nx / 2);

    //    Create WESN components of beta.g
    // std::cout << "[split_horizontal]  creating WESN components of beta" << std::endl;
    // g_W_beta = u_tau
    // g_E_beta = g_E_tau
    fc2d_hps_vector<double> g_S_beta = g_S_tau.extract(tau.grid.Nx / 2, tau.grid.Nx / 2);
    fc2d_hps_vector<double> g_N_beta = g_N_tau.extract(tau.grid.Nx / 2, tau.grid.Nx / 2);

    //    Intract into alpha.g
    // std::cout << "[split_horizontal]  creating g_alpha" << std::endl;
    fc2d_hps_vector<double> g_alpha(2*alpha.grid.Nx + 2*alpha.grid.Ny);
    g_alpha.intract(0, g_W_tau);
    g_alpha.intract(alpha.grid.Ny, u_tau);
    g_alpha.intract(2*alpha.grid.Ny, g_S_alpha);
    g_alpha.intract(2*alpha.grid.Ny + alpha.grid.Nx, g_N_alpha);

    //    Intract into beta.g
    // std::cout << "[split_horizontal]  creating g_beta" << std::endl;
    fc2d_hps_vector<double> g_beta(2*beta.grid.Nx + 2*beta.grid.Ny);
    g_beta.intract(0, u_tau);
    g_beta.intract(beta.grid.Ny, g_E_tau);
    g_beta.intract(2*beta.grid.Ny, g_S_beta);
    g_beta.intract(2*beta.grid.Ny + beta.grid.Nx, g_N_beta);

    // Set child patch data
    // std::cout << "[split_horizontal]  assigning child patch data" << std::endl;
    alpha.g = g_alpha;
    beta.g = g_beta;
    
    // std::cout << "[split_horizontal]  returning" << std::endl;
    return;

}

void uncoarsen_patch(fc2d_hps_patch& patch) {

    // Get options
    fclaw2d_global_t* glob = (fclaw2d_global_t*) patch.user;
    fc2d_hps_options_t* hps_opt = fc2d_hps_get_options(glob);

    // Build L21
	int N = patch.N_cells_leaf;
	// fc2d_hps_matrix<double> L21_west = build_L21(N*patch.coarsened->N_patch_side[WEST], N*patch.N_patch_side[WEST]);
	// fc2d_hps_matrix<double> L21_east = build_L21(N*patch.coarsened->N_patch_side[EAST], N*patch.N_patch_side[EAST]);
	// fc2d_hps_matrix<double> L21_south = build_L21(N*patch.coarsened->N_patch_side[SOUTH], N*patch.N_patch_side[SOUTH]);
	// fc2d_hps_matrix<double> L21_north = build_L21(N*patch.coarsened->N_patch_side[NORTH], N*patch.N_patch_side[NORTH]);
	// std::vector<fc2d_hps_matrix<double>> L21_diagonals = {L21_west, L21_east, L21_south, L21_north};
	// fc2d_hps_matrix<double> L21_patch = block_diag(L21_diagonals);

	// Build L12
	fc2d_hps_matrix<double> L12_west = build_L12(N*patch.N_patch_side[WEST], N*patch.coarsened->N_patch_side[WEST]);
	fc2d_hps_matrix<double> L12_east = build_L12(N*patch.N_patch_side[EAST], N*patch.coarsened->N_patch_side[EAST]);
	fc2d_hps_matrix<double> L12_south = build_L12(N*patch.N_patch_side[SOUTH], N*patch.coarsened->N_patch_side[SOUTH]);
	fc2d_hps_matrix<double> L12_north = build_L12(N*patch.N_patch_side[NORTH], N*patch.coarsened->N_patch_side[NORTH]);
	std::vector<fc2d_hps_matrix<double>> L12_diagonals = {L12_west, L12_east, L12_south, L12_north};
	fc2d_hps_matrix<double> L12_patch = block_diag(L12_diagonals);

    // Interpolate data down to original patch
    // printf("uncoarsening patch\n");
    printf("HERE1\n");
    printf("L12 = [%i, %i], g = [%i]\n", L12_patch.rows, L12_patch.cols, patch.coarsened->g.size());
    patch.g = L12_patch * patch.coarsened->g;
    if (hps_opt->nonhomogeneous_rhs) {
        // printf("w:\n");
        printf("HERE2\n");
        patch.w = L12_south * patch.coarsened->w;
    }

    // patch.print_info();
    // patch.coarsened->print_info();

    return;

}

void split_1to4(fc2d_hps_patch& parent, fc2d_hps_patch& child0, fc2d_hps_patch& child1, fc2d_hps_patch& child2, fc2d_hps_patch& child3) {

    // Assumptions on entry
    // if (parent.S.rows == 0 || parent.S.cols == 0) {
    //     throw std::invalid_argument("[fc2d_hps_split::split_1to4] `parent` patch solution matrix `S` does not have data.");
    // }
    // if (parent.g.size() == 0) {
    //     throw std::invalid_argument("[fc2d_hps_split::split_1to4] `parent` patch does not have Dirichelt data set.");
    // }

    // parent.print_info();
    // child0.print_info();
    // child1.print_info();
    // child2.print_info();
    // child3.print_info();

    // Check for coarsened versions, uncoarsen if exists
    if (parent.has_coarsened) {
        if (parent.coarsened->has_coarsened) {
            printf("HERE A\n");
            parent.print_info();
            parent.coarsened->print_info();
            parent.coarsened->coarsened->print_info();
            parent.coarsened->coarsened->g = parent.g;
            uncoarsen_patch(*parent.coarsened);
        }
        uncoarsen_patch(parent);
    }

    // Vertical split
    // Create patches to merge (either original or coarsened)
    fc2d_hps_patch* alpha;
    fc2d_hps_patch* beta;
    fc2d_hps_patch* gamma;
    fc2d_hps_patch* omega;

    // Set based on if patch has coarsened data
    child0.has_coarsened ? alpha = child0.coarsened : alpha = &child0;
    child1.has_coarsened ? beta = child1.coarsened : beta = &child1;
    child2.has_coarsened ? gamma = child2.coarsened : gamma = &child2;
    child3.has_coarsened ? omega = child3.coarsened : omega = &child3;

    //    Create patches for rectangular pieces
    fc2d_hps_patch alpha_prime;
    fc2d_hps_patch beta_prime;

    // alpha_prime.print_info();

    //    Perform split
    // std::cout << "[split_1to4]  begin vertical split" << std::endl;
    split_vertical(parent, alpha_prime, beta_prime);
    // NOTE: g and grid for alpha_prime and beta_prime are set inside split_vertical


    //    Assign alpha_prime and beta_prime the solution matrices
    alpha_prime.S = alpha->S_prime;
    alpha_prime.w = alpha->w_prime;
    alpha_prime.user = alpha->user;

    beta_prime.S = gamma->S_prime;
    beta_prime.w = gamma->w_prime;
    beta_prime.user = gamma->user;

    // alpha_prime.print_info();
    // beta_prime.print_info();
    // alpha_prime.print_info();


    // Horizontal split
    //    Bottom split
    // std::cout << "[split_1to4]  begin horizontal split 1" << std::endl;
    split_horizontal(alpha_prime, *alpha, *beta);
    
    //    Top split
    // std::cout << "[split_1to4]  begin horizontal split 2" << std::endl;
    split_horizontal(beta_prime, *gamma, *omega);
    return;

}