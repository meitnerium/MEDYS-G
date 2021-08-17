AUTOMAKE_OPTIONS = foreign subdir-objects
bin_PROGRAMS = medys 
medys_SOURCES =   src/all_v001/basics.f90 src/all_v001/variables.f90 src/all_v001/read_input.f90 src/all_v001/math.f90 src/all_v001/champ1.f90 src/all_v001/grille.f90 src/all_v001/NewSubroutines.f90 src/all_v001/Tests.f90 src/all_v001/string.f90 src/io_interface/io_module.f90 src/all_v001/precision.f90 src/core/error.f90 src/all_v001/matrix_modifications.f90 src/all_v001/matrix_temp_lowlvl.f90 src/all_v001/matrix_temp_highlvl.f90 src/all_v001/read_columbus.f90 src/all_v001/ortho.f90 src/all_v001/qp_integrals.f90 src/all_v001/observable.f90 src/all_v001/nouveaudrt.f90 src/all_v001/read_gamess.f90 src/all_v001/dynamique_new.f90 src/main/main.f90


# For debugging
medys_FCFLAGS = -I${MKLROOT}/include/intel64/lp64/ -check bounds -nogen-interface -fpe0 -g -traceback -debug extended -heap-arrays -fp-stack-check -debug all -qopenmp 
medys_LDADD = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -liomp5 -lmkl_blas95_lp64 -lmkl_lapack95_lp64

#For production
#medys_FCFLAGS = -I${MKLROOT}/include/intel64/lp64/ -nogen-interface -fpe0  -heap-arrays -fp-stack-check -qopenmp  -O2
#medys_LDADD = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -liomp5 -lmkl_blas95_lp64 -lmkl_lapack95_lp64
