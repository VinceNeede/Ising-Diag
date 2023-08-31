ifndef FILENAME
    FILENAME=""
endif

ifndef OUTDIR
    OUTDIR="_results"
endif

PATHINCLUDE="/opt/intel/oneapi/mkl/latest/include"

PATHLIB="/opt/intel/oneapi/mkl/latest/lib/intel64"
FLAGS = -O3 -static #-fast -msse4.2 -axAVX,CORE-AVX2

STATICLIB= /usr/local/lib/libarpackILP64.a ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lmkl_intel_ilp64 -liomp5 -lpthread -lm -ldl
INCLUDEOP=-I${MKLROOT}/include/intel64/ilp64 -i8  -I"${MKLROOT}/include" -I/usr/local/lib

compile: check_dir compile_spblas compile_lib
	@echo "Compiling Ising"
	@ifx $(FLAGS) -i8 -warn all -module ${OUTDIR} -fpp  -I$(PATHINCLUDE) -w  Ising.f90 -c -o ${OUTDIR}/Ising.o
	@ifx $(FLAGS) -i8 -warn all -module ${OUTDIR} -fpp -w ${OUTDIR}/Ising.o ${OUTDIR}/diag.o ${OUTDIR}/mkl_spblas.o -o ${OUTDIR}/Ising ${INCLUDEOP} ${STATICLIB}
	@echo "Ising compiled"
	@echo ""

	@echo "Compiling Ising_Parity"
	@ifx $(FLAGS) -i8 -warn all -module ${OUTDIR} -fpp  -I$(PATHINCLUDE) -w  Ising_Parity.f90 -c -o ${OUTDIR}/Ising_Parity.o
	@ifx $(FLAGS) -i8 -warn all -module ${OUTDIR} -fpp -w ${OUTDIR}/Ising_Parity.o ${OUTDIR}/diag.o ${OUTDIR}/mkl_spblas.o -o ${OUTDIR}/Ising_Parity ${INCLUDEOP} ${STATICLIB}
	@echo "Ising_Parity compiled"
	
	@echo ""
	@echo "Compiled Succesfully"
compile_lib:
	@echo "Compiling Library diagonalization"
	@ifx $(FLAGS) -i8 -warn all -module ${OUTDIR} -fpp  -I$(PATHINCLUDE) -w  diag.f90 -c -o ${OUTDIR}/diag.o
	@echo "Library diagonalization compiled"
	@echo ""

compile_spblas:
	@echo "compiling MKL_SPBLAS"
	@ifx $(FLAGS) -i8 -warn all -module ${OUTDIR} -fpp  -I$(PATHINCLUDE) -w  $(PATHINCLUDE)/mkl_spblas.f90 -c -o ${OUTDIR}/mkl_spblas.o 
	@echo "MKL_SPBLAS compiled"
	@echo ""
check_dir:
	@[ -d "${OUTDIR}" ] || mkdir ${OUTDIR}