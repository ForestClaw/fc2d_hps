# -*- Makefile -*-
libhps_generated_headers =

libhps_installed_headers = \
	src/solvers/fc2d_hps/fc2d_hps.h \
	src/solvers/fc2d_hps/fc2d_hps_options.h \
	src/solvers/fc2d_hps/fc2d_hps_physical_bc.h

libhps_compiled_sources = \
	src/solvers/fc2d_hps/fc2d_hps.cpp \
	src/solvers/fc2d_hps/fc2d_hps_options.c \
	src/solvers/fc2d_hps/fc2d_hps_physical_bc.c \
	src/solvers/fc2d_hps/fc2d_hps_output_ascii.c \
	src/solvers/fc2d_hps/fc2d_hps_diagnostics.c \
 	src/solvers/fc2d_hps/fortran_source/fc2d_hps_apply_bc.f90 \
	src/solvers/fc2d_hps/fortran_source/fc2d_hps_eval_bc.f90 \
	src/solvers/fc2d_hps/fortran_source/fc2d_hps_tag4refinement.f90 \
	src/solvers/fc2d_hps/fortran_source/fc2d_hps_tag4coarsening.f90

## Name of library to build
lib_LTLIBRARIES += src/solvers/fc2d_hps/libfc2d_hps.la

## Named variables that can be referenced from other libraries/apps
FCLAW_HPS_CPPFLAGS = -I@top_srcdir@/src/solvers/fc2d_hps \
	                       -I@top_srcdir@/src/solvers/fc2d_hps/patch_solvers \
						   $(HPS_CFLAGS)

FCLAW_HPS_LDADD    = @top_builddir@/src/solvers/fc2d_hps/libfc2d_hps.la \
						   $(HPS_LIBS)

## Sources needed to build this library
src_solvers_fc2d_hps_libfc2d_hps_la_SOURCES = \
          $(libhps_internal_headers) \
          $(libhps_compiled_sources)

## Headers needed to build this library.  Only AM_CPPFLAGS required; others are 
## included in AM_CPPFLAGS, but are included here for clarity (although makes 
## for lengthy compile lines)
src_solvers_fc2d_hps_libfc2d_hps_la_CPPFLAGS = \
          $(AM_CPPFLAGS) \
          $(FCLAW_CPPFLAGS) \
          $(FCLAW_HPS_CPPFLAGS) \
          $(FCLAW_CLAWPATCH_CPPFLAGS)

## Add headers from this library to general header pool, to be used by other apps
AM_CPPFLAGS  += $(FCLAW_HPS_CPPFLAGS)

LDADD        += $(FCLAW_HPS_LDADD)

include_HEADERS += $(libhps_installed_headers)
