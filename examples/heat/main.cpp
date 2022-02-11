#include <iostream>

#include <fclaw2d_include_all.h>
#include <fc2d_hps.h>
#include <fc2d_hps_options.h>

// Forward declarations
static fclaw2d_domain_t* create_domain(sc_MPI_Comm mpi_comm, fclaw_options_t* fclaw_opt);

int main(int argc, char* argv[]) {

    // Create app
    fclaw_app_t* fclaw_app = fclaw_app_new(&argc, &argv, NULL);

    // Read in options and setup
    int first_arg;
    fclaw_options_t* fclaw_opt = fclaw_options_register(fclaw_app, argv[1]);
    fclaw2d_clawpatch_options_t* clawpatch_opt = fclaw2d_clawpatch_options_register(fclaw_app, argv[1]);
    fc2d_hps_options_t* hps_opt = fc2d_hps_options_register(fclaw_app, argv[1]);
    sc_options_t* sc_opt = fclaw_app_get_options(fclaw_app);
    int return_value = fclaw_options_read_from_file(sc_opt);
    fclaw_exit_type_t vexit = fclaw_app_options_parse(fclaw_app, &first_arg, "fclaw_options_cache.ini");

    // Run program if options are good
    if (!return_value && !vexit) {
        
        // Get MPI information
        int mpi_rank;
        int mpi_size;
        sc_MPI_Comm mpi_comm = fclaw_app_get_mpi_size_rank(fclaw_app, &mpi_size, &mpi_rank);

        // Create ForestClaw domain and glob
        // fclaw2d_domain_t* fc_domain = create_domain(mpi_comm, fclaw_opt);
        fclaw2d_domain_t* fc_domain;
        fclaw2d_domain_data_new(fc_domain);
        fclaw2d_global_t* fc_glob = fclaw2d_global_new();
        // fclaw2d_global_store_domain(fc_glob, fc_domain);

        
        // Destroy glob
        fclaw2d_global_destroy(fc_glob);

    }

    fclaw_app_destroy(fclaw_app);
    return 0;
}

static fclaw2d_domain_t* create_domain(sc_MPI_Comm mpi_comm, fclaw_options_t* fclaw_opt) {

    // Mapped, multi-block domain
    p4est_connectivity_t* conn = NULL;
    fclaw2d_domain_t* domain;
    fclaw2d_map_context_t* cont = NULL;
    fclaw2d_map_context_t* brick = NULL;

    // Problem domain
    int mi = fclaw_opt->mi;
    int mj = fclaw_opt->mj;
    int a = fclaw_opt->periodic_x;
    int b = fclaw_opt->periodic_y;

    // Map unit square to disk
    // conn = p4est_connectivity_new_brick(mi,mj,a,b);
    // brick = fclaw2d_map_new_brick(conn,mi,mj);
    // cont = fclaw2d_map_new_nomap_brick(brick);

    // domain = fclaw2d_domain_new_conn_map (mpi_comm, fclaw_opt->minlevel, conn, cont);
    // fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    // fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);  
    return domain;

}