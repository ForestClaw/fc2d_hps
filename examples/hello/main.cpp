#include <iostream>

#include <fclaw2d_include_all.h>
#include <HPS/fc2d_hps.hpp>
#include <Util/fc2d_hps_options.h>

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
        // Global stuff
        int mpi_comm;
        int mpi_rank;
        fclaw_app_get_mpi_size_rank(fclaw_app, &mpi_comm, &mpi_rank);

        // Say hello!
        if (mpi_rank == 0) {
            std::cout << "[hello::main.cpp::main]  Hello from ForestClaw::HPS! This is the head_rank." << std::endl;
        }
        else {
            std::cout << "[hello::main.cpp::main]  Hello from ForestClaw::HPS! This is rank " << mpi_rank << std::endl;
        }

    }

    fclaw_app_destroy(fclaw_app);
    return 0;
}