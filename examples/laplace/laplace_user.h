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

#ifndef LAPLACE_USER_H
#define LAPLACE_USER_H

/* ForestClaw headers */
#include <fclaw2d_include_all.h>

#include <fclaw2d_output.h>
#include <fclaw2d_diagnostics.h>

#include <fclaw2d_elliptic_solver.h>

#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_clawpatch.h>

#include <fc2d_hps.h>
#include <fc2d_hps_options.h>
#include <fc2d_hps_output_ascii.h>

/* Application headers */
#include "laplace_options.h"


#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif


/* --------------------------- Problem dependent functions -----------------------------*/

void laplace_link_solvers(fclaw2d_global_t *glob);


/* --------------------------- Fortran functions ---------------------------------------*/

#define SETPROB FCLAW_F77_FUNC(setprob,SETPROB)

void SETPROB();



#define LAPLACE_FORT_RHS FCLAW_F77_FUNC(laplace_fort_rhs,LAPLACE_FORT_RHS)

void LAPLACE_FORT_RHS(const int* blockno, const int* mbc, const int* mx, 
                     const int* my, const int* mfields, 
                     const double *xlower, const double *ylower,
                     const double* dx, const double* dy, double rhs[]);


#define LAPLACE_COMPUTE_ERROR FCLAW_F77_FUNC(laplace_compute_error,LAPLACE_COMPUTE_ERROR)

void LAPLACE_COMPUTE_ERROR(const int* blockno, 
                           const int *mx, 
                           const int *my, 
                           const int* mbc, 
                           const int* mfields,
                           const double *dx, 
                           const double *dy, 
                           const double *xlower,
                           const double *ylower, 
                           const double *t, double q[],
                           double error[], double soln[]);


#define LAPLACE_FORT_APPLY_BC FCLAW_F77_FUNC(laplace_fort_apply_bc, \
                                            LAPLACE_FORT_APPLY_BC)

void LAPLACE_FORT_APPLY_BC(const int* blockno, const  int* mx, const  int* my, 
                          const  int* mbc, const  int* mfields, 
                          const double* xlower, const double* ylower,
                          const double* dx, const double* dy, const double* t,
                          int intersects_bc[], int mthbc[], 
                          double rhs[], fc2d_hps_fort_eval_bc_t g_bc, 
                          int* cons_check, double flux_sum[]);


#define LAPLACE_FORT_EVAL_BC FCLAW_F77_FUNC(laplace_fort_eval_bc, LAPLACE_FORT_EVAL_BC)

double LAPLACE_FORT_EVAL_BC(const int* iface, const double* t,
                            const double* x, const double* y);


/* ----------------------------- Fortran - output functions --------------------------- */

#define  LAPLACE_FORT_OUTPUT_ASCII \
           FCLAW_F77_FUNC(laplace_fort_output_ascii, \
                          LAPLACE_FORT_OUTPUT_ASCII)
void LAPLACE_FORT_OUTPUT_ASCII(const char* matname1,
                              int* mx,        int* my,
                              int* meqn,      int* mbc,
                              double* xlower, double* ylower,
                              double* dx,     double* dy,
                              double q[],double soln[], double error[],
                              int* patch_num, int* level,
                              int* blockno,   int* mpirank);

#if 0
#define LAPLACE_FORT_HEADER_ASCII \
         FCLAW_F77_FUNC(laplace_fort_header_ascii, \
                        LAPLACE_FORT_HEADER_ASCII)
void LAPLACE_FORT_HEADER_ASCII(char* matname1, char* matname2,
                              double* time, int* meqn, int* maux, 
                              int* ngrids);
#endif




#ifdef __cplusplus
}
#endif

#endif
