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

#ifndef FC2D_HPS_H
#define FC2D_HPS_H

/* Include headers here that are needed for other HPS routines */
#include <fclaw2d_include_all.h>

#include <fclaw2d_elliptic_solver.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_clawpatch_output_ascii.h> 
#include <fclaw2d_clawpatch_output_vtk.h>

#ifdef __cplusplus
extern "C"
{
#endif

typedef  struct fc2d_hps_vtable  fc2d_hps_vtable_t;


/* --------------------------- Fortran defs solver functions -------------------------- */

typedef void (*fc2d_hps_fort_qexact_t)(const int* example,
                                       const double* x,
                                       const double* y,
                                       const double* q,
                                       const double* qlap,
                                       const double* grad,
                                       const int* flag);

typedef  void (*fc2d_hps_fort_rhs_t)(const int* blockno, 
                                     const int* mbc,
                                     const int* mx, const int* my,
                                     const int* mfields,
                                     const double* xlower, const double* ylower,
                                     const double* dx, const double* dy,
                                     double rhs[]);

typedef  void (*fc2d_hps_fort_beta_t)(const double* x,
                                      const double* y,
                                      const double* beta,
                                      double grad[]);

typedef double (*fc2d_hps_fort_eval_bc_t)(const int *iface, const double *t, 
                                          const double *x, const double *y);

typedef void (*fc2d_hps_fort_apply_bc_t)(const int* blockno, const  int* mx, const  int* my, 
                                         const  int* mbc, const  int* meqn, 
                                         const double* xlower, const double* ylower,
                                         const double* dx, const double* dy, 
                                         const double *t, 
                                         int intersects_bc[], int mthbc[], 
                                         double rhs[], fc2d_hps_fort_eval_bc_t g_bc, 
                                         int* cons_check, double flux_sum[]);

typedef void (*fc2d_hps_patch_solver_t)(struct fclaw2d_global *glob);


typedef void (*fc2d_hps_fort_output_t)(const char* matname1,
                                       int* mx,        int* my,
                                       int* meqn,      int* mbc,
                                       double* xlower, double* ylower,
                                       double* dx,     double* dy,
                                       double q[],double soln[], double error[],
                                       int* patch_num, int* level,
                                       int* blockno,   int* mpirank);


/* -------------------------- Solver and utilities ------------------------------------ */

//void fc2d_hps_solve(struct fclaw2d_global *glob);

fc2d_hps_vtable_t* fc2d_hps_vt();


/* --------------------------------- Virtual table ------------------------------------ */

struct fc2d_hps_vtable
{
    /* Solver that defines patch solver and calls Hps solver */
    fc2d_hps_patch_solver_t   patch_solver;  /* 'solver' is a keyword */

    /* Fortran routines */
    fc2d_hps_fort_rhs_t        fort_rhs;	
    fc2d_hps_fort_beta_t       fort_beta;	
    fc2d_hps_fort_apply_bc_t   fort_apply_bc;
    fc2d_hps_fort_eval_bc_t    fort_eval_bc;
    fc2d_hps_fort_qexact_t     fort_qexact;

    /* Allows us to output error and exact solution, along with computed solution */
    fc2d_hps_fort_output_t     fort_output;   

	  int is_set;
};

void fc2d_hps_solver_initialize(void);

fc2d_hps_vtable_t* fc2d_hps_vt(void);


/* ----------------------------- User access to solver functions ---------------------- */

void fc2d_hps_setprob(fclaw2d_global_t* glob);


void fc2d_hps_rhs(fclaw2d_global_t* glob,
                  fclaw2d_patch_t *patch,
                  int blockno,
                  int patchno);


/* -------------------------------- solver utilities -------------------------------- */

/* Put this here so that user does not have to include fc2d_hps_heat.h */
void fc2d_hps_heat_set_lambda(double lambda);

double fc2d_hps_heat_get_lambda();



#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif /* !FC2D_HPS_H */
