/*
Copyright (c) 2019 Carsten Burstedde, Donna Calhoun, Damyn Chipman
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

#ifndef SIMPLE_USER_H
#define SIMPLE_USER_H

/* ForestClaw headers */
#include <fclaw2d_include_all.h>

#include <fclaw2d_output.h>
#include <fclaw2d_diagnostics.h>

#include <fclaw2d_elliptic_solver.h>

#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_clawpatch.h>

#include <HPS/fc2d_hps.hpp>
#include <Util/fc2d_hps_options.h>
#include <Util/fc2d_hps_output_ascii.h>
// #include <Structures/fc2d_hps_patchsolver.hpp>

/* HPS headers */
// #include <fc2d_hps_patch.hpp>
// #include <fc2d_hps_quadtree.hpp>
// #include <fc2d_hps_interface.hpp>

/* Application headers */
#include "simple_options.h"


#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif


/* --------------------------- Problem dependent functions -----------------------------*/

void simple_link_solvers(fclaw2d_global_t *glob);


/* --------------------------- Fortran functions ---------------------------------------*/

#define SETPROB FCLAW_F77_FUNC(setprob,SETPROB)

void SETPROB();

#define SIMPLE_FORT_QEXACT_COMPLETE FCLAW_F77_FUNC(simple_fort_qexact_rhs,SIMPLE_FORT_QEXACT_COMPLETE)

void SIMPLE_FORT_QEXACT_COMPLETE(const int* example,
                                 const double* x, const double* y,
                                 const double* q, const double* qlap, const double* grad,
                                 const int* flag);

#define SIMPLE_FORT_RHS FCLAW_F77_FUNC(simple_fort_rhs,SIMPLE_FORT_RHS)

void SIMPLE_FORT_RHS(const int* blockno, const int* mbc, const int* mx, 
                     const int* my, const int* mfields, 
                     const double *xlower, const double *ylower,
                     const double* dx, const double* dy, double rhs[]);


#define SIMPLE_COMPUTE_ERROR FCLAW_F77_FUNC(simple_compute_error,SIMPLE_COMPUTE_ERROR)

void SIMPLE_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* mfields,
                        double *dx, double *dy, double *xlower,
                        double *ylower, double *t, double q[],
                           double error[], double soln[]);


#define SIMPLE_FORT_APPLY_BC FCLAW_F77_FUNC(simple_fort_apply_bc, \
                                            SIMPLE_FORT_APPLY_BC)

void SIMPLE_FORT_APPLY_BC(const int* blockno, const  int* mx, const  int* my, 
                          const  int* mbc, const  int* mfields, 
                          const double* xlower, const double* ylower,
                          const double* dx, const double* dy, const double* t,
                          int intersects_bc[], int mthbc[], 
                          double rhs[], fc2d_hps_fort_eval_bc_t g_bc, 
                          int* cons_check, double flux_sum[]);


#define SIMPLE_FORT_EVAL_BC FCLAW_F77_FUNC(simple_fort_eval_bc, SIMPLE_FORT_EVAL_BC)

double SIMPLE_FORT_EVAL_BC(const int* iface, const double* t,
                            const double* x, const double* y);


/* ----------------------------- Fortran - output functions --------------------------- */

#define  SIMPLE_FORT_OUTPUT_ASCII \
           FCLAW_F77_FUNC(simple_fort_output_ascii, \
                          SIMPLE_FORT_OUTPUT_ASCII)
void SIMPLE_FORT_OUTPUT_ASCII(const char* matname1,
                              int* mx,        int* my,
                              int* meqn,      int* mbc,
                              double* xlower, double* ylower,
                              double* dx,     double* dy,
                              double q[],double soln[], double error[],
                              int* patch_num, int* level,
                              int* blockno,   int* mpirank);

#if 0
#define SIMPLE_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(simple_fort_header_ascii, \
                        SIMPLE_FORT_HEADER_ASCII)
void SIMPLE_FORT_HEADER_ASCII(char* matname1, char* matname2,
                              double* time, int* meqn, int* maux, 
                              int* ngrids);
#endif




#ifdef __cplusplus
}
#endif

#endif
