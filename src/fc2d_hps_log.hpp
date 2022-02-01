#ifndef FC2D_HPS_LOG_HPP_
#define FC2D_HPS_LOG_HPP_
#pragma once

#include <iostream>
#include <fstream>
#include <string>

enum fc2d_hps_verbosity_options {
    NO_OUTPUT = 0,
    ONLY_LOG = 1,
    ONLY_CONSOLE = 2,
    
};

class fc2d_hps_log {

public:

    std::size_t verbosity;
    std::ofstream log_file;

    fc2d_hps_log(std::size_t verbosity, std::string log_filename);
    ~fc2d_hps_log();

    

private:


};


#endif // FC2D_HPS_LOG_HPP_