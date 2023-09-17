ifndef FILENAME
    FILENAME=""
endif

FLAGS = -i8 -r8 -O3 -static#-fast -msse4.2 -axAVX,CORE-AVX2

ifeq ($(OS),Windows_NT)
STATICLIB= mkl_intel_ilp64.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib /link /LIBPATH:"C:/Program Files (x86)/arpack/lib/" arpackILP64.lib /NODEFAULTLIB:MSVCRT
INCLUDEOP= /4I8  -I"%MKLROOT%\include"
OPTIONS=-warn:all -module:${OUTDIR}
VOIDLINE=@echo.
PATHINCLUDE="${MKLROOT}\include"
ifndef OUTDIR
    OUTDIR=win_results/
endif
else
STATICLIB= /usr/local/lib/libarpackILP64.a ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lmkl_intel_ilp64 -liomp5 -lpthread -lm -ldl
INCLUDEOP=-I${MKLROOT}/include/intel64/ilp64 -i8  -I"${MKLROOT}/include" -I/usr/local/lib
OPTIONS=-warn all -module ${OUTDIR}
VOIDLINE=@echo ""
PATHINCLUDE="${MKLROOT}/include"
ifndef OUTDIR
    OUTDIR=Linux_results/
endif
endif


compile: check_dir compile_spblas compile_lib
	@echo "Compiling Quenching"
	@ifx $(FLAGS) ${OPTIONS} -fpp  -I$(PATHINCLUDE) -w  SRC/quenching.f90 -c -o ${OUTDIR}/quenching.o
	@ifx $(FLAGS) ${OPTIONS} -fpp -w ${OUTDIR}/quenching.o ${OUTDIR}/diag.o ${OUTDIR}/mkl_spblas.o -o ${OUTDIR}/quenching ${INCLUDEOP} ${STATICLIB}
	@echo "Quenching compiled"
	$(VOIDLINE)
	@echo "Compiled Succesfully"
compile_lib: check_dir compile_spblas
	@echo "Compiling Library diagonalization"
	@ifx $(FLAGS) ${OPTIONS} -fpp  -I$(PATHINCLUDE) -w  SRC/diag.f90 -c -o ${OUTDIR}/diag.o
	@echo "Library diagonalization compiled"
	$(VOIDLINE)

compile_spblas: 
	@echo "compiling MKL_SPBLAS"
	@ifx $(FLAGS) ${OPTIONS} -fpp  -I$(PATHINCLUDE) -w  $(PATHINCLUDE)/mkl_spblas.f90 -c -o ${OUTDIR}/mkl_spblas.o 
	@echo "MKL_SPBLAS compiled"
	$(VOIDLINE)

check_dir:
	@echo "Checking for directory existence..."
ifeq ($(OS),Windows_NT)
	@if not exist "$(OUTDIR)" ( \
		mkdir $(OUTDIR); \
	) 
else
	@if [ ! -d "$(OUTDIR)" ]; then \
		mkdir $(OUTDIR); \
	fi
endif
#@[ -d "${OUTDIR}" ] || mkdir ${OUTDIR}