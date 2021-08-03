## Makefile.am in src/solvers/fc2d_hps/examples/laplace

bin_PROGRAMS += src/solvers/fc2d_hps/examples/laplace/laplace

src_solvers_fc2d_hps_examples_laplace_laplace_SOURCES = \
	src/solvers/fc2d_hps/examples/laplace/laplace.cpp \
	src/solvers/fc2d_hps/examples/laplace/laplace_user.cpp \
	src/solvers/fc2d_hps/examples/laplace/laplace_user.h \
	src/solvers/fc2d_hps/examples/laplace/laplace_options.c \
	src/solvers/fc2d_hps/examples/laplace/laplace_diagnostics.c \
	src/solvers/fc2d_hps/examples/laplace/hsmooth.f90 \
	src/solvers/fc2d_hps/examples/laplace/fortran/laplace_setprob.f90 \
	src/solvers/fc2d_hps/examples/laplace/fortran/laplace_fort_rhs.f90 \
	src/solvers/fc2d_hps/examples/laplace/fortran/laplace_compute_error.f \
	src/solvers/fc2d_hps/examples/laplace/fortran/laplace_qexact.f90 \
	src/solvers/fc2d_hps/examples/laplace/fortran/laplace_apply_bc.f90 \
	src/solvers/fc2d_hps/examples/laplace/fortran/laplace_eval_bc.f90 \
	src/solvers/fc2d_hps/examples/laplace/fortran/laplace_output.f90


EXTRA_DIST += \
	src/solvers/fc2d_hps/examples/laplacefclaw_options.ini \
	src/solvers/fc2d_hps/examples/laplace/regression.ini \
	src/solvers/fc2d_hps/examples/laplace/afterframe.m \
	src/solvers/fc2d_hps/examples/laplace/setplot2.m

## Include headers and libraries needed to build this application
##  -- Only really need AM_CPPFLAGS and LDADD (since these include all the other 
##     variables).  But others are included here for completeness

src_solvers_fc2d_hps_examples_laplace_laplace_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_HPS_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)


src_solvers_fc2d_hps_examples_laplace_laplace_LDADD = \
        $(LDADD) \
        $(FCLAW_HPS_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
