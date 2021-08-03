/*
Copyright (c) 2019 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#ifndef FC2D_HPS_OPTIONS_H
#define FC2D_HPS_OPTIONS_H

#include <fclaw_base.h>

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif

struct fclaw2d_global;

typedef struct fc2d_hps_options fc2d_hps_options_t;


typedef enum {
    LAPLACE = 0,      /* Laplacian (no beta) */    
    VARPOISSON,         /* Variable Poisson operator */
    HEAT,               /* Variable Poisson operator */
    USER_OPERATOR
} fc2d_hps_operator_types;

typedef enum {
    FFT = 0,    /* Must use fivepoint */
    FISHPACK,   /* Must use fivepoint */
    DST,        /* Must use fivepoint */
    BICG,       /* Can be used with variable coefficient */
    USER_SOLVER
} fc2d_hps_solver_types;


struct fc2d_hps_options
{
    /* Boundary conditions */
    int *boundary_conditions;
    const char *bc_cond_string;

    /* Output */
    int ascii_out;
    int vtk_out;

    int operator_type;  /* laplace, varpoisson, etc ... */
    sc_keyvalue_t *kv_operator_type;

    int patch_solver;                /* FFT, DST, FISHPACK, ... */
    sc_keyvalue_t *kv_patch_solver;

    int is_registered;
};

#if 0
fclaw_exit_type_t fc2d_hps_postprocess (fc2d_hps_options_t *
                                               mg_opt);

fclaw_exit_type_t fc2d_hps_check (fc2d_hps_options_t * mg_opt);

void fc2d_hps_reset (fc2d_hps_options_t * mg_opt);
#endif

fc2d_hps_options_t*  fc2d_hps_options_register (fclaw_app_t * app,
                                                const char *configfile);

void fc2d_hps_package_register(fclaw_app_t* app,
                               fc2d_hps_options_t* mg_opt);

fc2d_hps_options_t* fc2d_hps_get_options(struct fclaw2d_global *glob);

void fc2d_hps_options_store (struct fclaw2d_global* glob, 
                             fc2d_hps_options_t* mg_opt);

#ifdef __cplusplus
}
#endif

#endif
