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
    fc2d_hps_vector<double> u_tau = tau.S * tau.g;
    if (hps_opt->nonhomogeneous_rhs) {
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
    fc2d_hps_vector<double> u_tau = tau.S * tau.g;
    if (hps_opt->nonhomogeneous_rhs) {
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

void split_1to4(fc2d_hps_patch& parent, fc2d_hps_patch& child0, fc2d_hps_patch& child1, fc2d_hps_patch& child2, fc2d_hps_patch& child3) {

    // Assumptions on entry
    if (parent.S.rows == 0 || parent.S.cols == 0) {
        throw std::invalid_argument("[fc2d_hps_split::split_1to4] `parent` patch solution matrix `S` does not have data.");
    }
    if (parent.g.size() == 0) {
        throw std::invalid_argument("[fc2d_hps_split::split_1to4] `parent` patch does not have Dirichelt data set.");
    }

    // Vertical split
    //    Create patches for rectangular pieces
    fc2d_hps_patch alpha_prime;
    fc2d_hps_patch beta_prime;

    // alpha_prime.print_info();

    //    Perform split
    // std::cout << "[split_1to4]  begin vertical split" << std::endl;
    split_vertical(parent, alpha_prime, beta_prime);
    // NOTE: g and grid for alpha_prime and beta_prime are set inside split_vertical

    // alpha_prime.print_info();

    //    Assign alpha_prime and beta_prime the solution matrices
    alpha_prime.S = child0.S_prime;
    alpha_prime.w = child0.w_prime;
    alpha_prime.user = child0.user;

    beta_prime.S = child2.S_prime;
    beta_prime.w = child2.w_prime;
    beta_prime.user = child2.user;

    // alpha_prime.print_info();


    // Horizontal split
    //    Bottom split
    // std::cout << "[split_1to4]  begin horizontal split 1" << std::endl;
    split_horizontal(alpha_prime, child0, child1);
    
    //    Top split
    // std::cout << "[split_1to4]  begin horizontal split 2" << std::endl;
    split_horizontal(beta_prime, child2, child3);
    return;

}