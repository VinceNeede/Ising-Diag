ifndef FILENAME
    FILENAME=""
endif

PATHINCLUDE="/opt/intel/oneapi/mkl/latest/include"

PATHLIB="/opt/intel/oneapi/mkl/latest/lib/intel64"
FLAGS = -O3 -static #-fast -msse4.2 -axAVX,CORE-AVX2

STATICLIB= /usr/local/lib/libarpackILP64.a ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lmkl_intel_ilp64 -liomp5 -lpthread -lm -ldl
INCLUDEOP=-I${MKLROOT}/include/intel64/ilp64 -i8  -I"${MKLROOT}/include" -I/usr/local/lib

compile: compile_spblas compile_lib
	@echo "Compiling Ising"
	@ifx $(FLAGS) -i8 -warn all -module _results -fpp  -I$(PATHINCLUDE) -w  Ising.f90 -c -o _results/Ising.o
	@ifx $(FLAGS) -i8 -warn all -module _results -fpp -w _results/Ising_Parity.o _results/diag.o _results/mkl_spblas.o -o _results/Ising_Parity ${INCLUDEOP} ${STATICLIB}
	@echo "Ising compiled"
	@echo ""

	@echo "Compiling Ising_Parity"
	@ifx $(FLAGS) -i8 -warn all -module _results -fpp  -I$(PATHINCLUDE) -w  Ising_Parity.f90 -c -o _results/Ising_Parity.o
	@ifx $(FLAGS) -i8 -warn all -module _results -fpp -w _results/Ising.o _results/diag.o _results/mkl_spblas.o -o _results/Ising ${INCLUDEOP} ${STATICLIB}
	@echo "Ising_Parity compiled"
	
	@echo ""
	@echo "Compiled Succesfully"
compile_lib:
	@echo "Compiling Library diagonalization"
	@ifx $(FLAGS) -i8 -warn all -module _results -fpp  -I$(PATHINCLUDE) -w  diag.f90 -c -o _results/diag.o
	@echo "Library diagonalization compiled"
	@echo ""

compile_spblas:
	@echo "compiling MKL_SPBLAS"
	@ifx $(FLAGS) -i8 -warn all -module _results -fpp  -I$(PATHINCLUDE) -w  $(PATHINCLUDE)/mkl_spblas.f90 -c -o _results/mkl_spblas.o 
	@echo "MKL_SPBLAS compiled"
	@echo ""