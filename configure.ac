AC_INIT([medys], [0.1], [meitnerium109@gmail.com])
m4_define([_AC_F95_FC], [ifort])
AC_CONFIG_SRCDIR([src/main/main.f90])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE
#AC_SEARCH_LIBS([mkl_lapack])
#AC_SEARCH_LIBS([pthread])
#AC_CHECK_LIBS([pthread])

AC_CHECK_LIB([pthread], [main],[PTHREAD=-lpthread],echo "Error:Required library pthread not found. Install it and try again" && exit 1)
#AC_CHECK_LIB([mkl_intel_lp64], [main],[FFTW_LIBS="-Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -WI, --end-group -liomp5 -lpthread"],echo "Error:Required library mkl not found. Install it and try again" && exit 1)
#AC_CHECK_LIB(mkl_intel_lp64, fftw_plan_dft_1d, [FFTW_LIBS="-Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -WI, --end-group -liomp5 -lpthread"; with_fftw3=yes],AC_MSG_RESULT("MKL library not found.  Octave will use the FFTW3 instead."),[-Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread])   
AC_PROG_FC 
AC_PROG_CC 
#AM_FCFLAGS([-lmkl_intel_lp64])
#AM_FFLAGS = ${PTHREAD}
AC_CONFIG_FILES([Makefile])
AC_OUTPUT


