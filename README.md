# Ising #
This library diagonalize the 1-D Ising Hamiltonian:
$$\hat H = -J \sum_j \hat X_j \hat X_{j+1} - g \sum_j \hat Z_j +h \sum_j \hat X_j$$
Using the `arpack-ng` library, build with `MKL`, and storing the matrix in a sparse format using `MKL_SPBLAS` library.
# Prerequisites #
These are the libraries needed by the program.
## For Ubuntu ##
### Intel Compiler and MKL ###
  The MKL Library contains:

* Blas:
    * Level 1: Vector and Vector-Vector operations
    * Level 2: Matrix-Vector operations
    * Level 3: Matrix-Matrix operations
* Lapack: Provides routines for solving systems of simultaneous linear equations, least-squares solutions of linear systems of equations, eigenvalue problems, and singular value problems.
* MKL_SPBLAS: An implementation of Blas Library for Sparse Matrices.

Actually Blas and Lapack exists outside MKL, but Intel implementation is more optimized and parallelized.

MKL is part of the Intel oneAPI Base Toolkit, it can be downloaded from the relative [site](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html).

For Ubuntu users the online installer can be downloaded by writing the following lines in the Command Line:

```bash
wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/992857b9-624c-45de-9701-f6445d845359/l_BaseKit_p_2023.2.0.49397.sh

sudo sh ./l_BaseKit_p_2023.2.0.49397.sh
```
and follow the installer instructions.

The Intel compilers (ifort for classic Fortran, ifx for modern Fortran) are part of the Intel oneAPI HPC Toolkit, which can be downloaded from the relative [site](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html).

For Ubuntu users the online installer can be downloaded by writing the following lines in the Command Line:

```bash
wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/0722521a-34b5-4c41-af3f-d5d14e88248d/l_HPCKit_p_2023.2.0.49440.sh

sudo sh ./l_HPCKit_p_2023.2.0.49440.sh
```
and follow the installer instructions.

### ARPACK-NG ###
Arpack contains implementations of Lanczos and Arnoldi algorithms for extreme eigenvalue problems. Such library is no longer mantained and supported, for this reason there exists an open collaboratio called ARPACK-NG, which is still not the
best modern implementation of the algorithms (see the great number of go-to as an example), but the best we have so far. 

Before installing ARPACK-NG make sure you have installed cmake, this can be veriefied by writing in the command shell:
```bash
which cmake
```
if there is no output, then cmake is not installed, if the output looks like `/usr/bin` it is correctly installed. 

You can install cmake from pip:

```bash
python3 -m pip install cmake
which cmake && cmake --version
```

Now we have to set up the Intel environment, this is done by writing
```bash
source /opt/intel/oneapi/setvars.sh
```

We are now going to download the repository of ARPACK-NG from git and move inside the folder

```bash
git clone https://github.com/opencollab/arpack-ng.git
cd ./arpack-ng
```

Inside the arpack-ng folder, create a dir called `build` and move inside that folder:
```bash
mkdir build
cd build
```

We are now going to set the installation:

```bash
cmake -D INTERFACE64=ON -D ITF64SUFFIX="ILP64" -D BUILD_SHARED_LIBS=ON -D BLA_VENDOR=Intel10_64ilp ..
```

Please make sure the Configuration Summary (the end of the output of the last command) looks like this:
```bash
-- Configuration summary for arpack-ng-3.9.0:
   -- prefix: /usr/local
   -- MPI: OFF (ICB provided )
   -- ICB: OFF
   -- INTERFACE64: 1
   -- FC:      /opt/intel/oneapi/compiler/2023.2.0/linux/bin/intel64/ifort
   -- FCFLAGS: -O3  -i8
   -- BLAS:
      -- link:    /opt/intel/oneapi/mkl/2023.2.0/lib/intel64/libmkl_intel_ilp64.so
      -- link:    /opt/intel/oneapi/mkl/2023.2.0/lib/intel64/libmkl_intel_thread.so
      -- link:    /opt/intel/oneapi/mkl/2023.2.0/lib/intel64/libmkl_core.so
      -- link:    /opt/intel/oneapi/compiler/2023.2.0/linux/compiler/lib/intel64_lin/libiomp5.so
      -- link:    -lm
      -- link:    -ldl
   -- LAPACK:
      -- link:    /opt/intel/oneapi/mkl/2023.2.0/lib/intel64/libmkl_intel_ilp64.so
      -- link:    /opt/intel/oneapi/mkl/2023.2.0/lib/intel64/libmkl_intel_thread.so
      -- link:    /opt/intel/oneapi/mkl/2023.2.0/lib/intel64/libmkl_core.so
      -- link:    /opt/intel/oneapi/compiler/2023.2.0/linux/compiler/lib/intel64_lin/libiomp5.so
      -- link:    -lm
      -- link:    -ldl
      -- link:    -lm
      -- link:    -ldl
-- Configuring done
-- Generating done

```
In particular:
  * The compiler is ifort, this can be checked by the flag --FC
  * the BLAS and LAPACK library linked is `libmkl_intel_ilp64.so` and not `libmkl_intel_lp64.so`, this can be seen from the first row after --LAPACK and --BLAS. The difference between the two is the integer representation (8 bytes for ilp64 and
    4 bytes for lp64, we are going to need 8 bytes)

We can build the libraries by running just: 
```bash
make
```
and then install it by running:
```bash
sudo make install
```

We have isntalled the shared libraries, we are going to need the static library also, so execute the following commands:
```bash
cd ../
rm -r -f build
mkdir build
cd build
cmake -D INTERFACE64=ON -D ITF64SUFFIX="ILP64" -D BUILD_SHARED_LIBS=OFF -D BLA_VENDOR=Intel10_64ilp ..
make
sudo make install
```
✨ Congratulations, you have installed `arpack` lib.

# HOW TO USE ISING #
<!-- ## Clone the repository ##
First you have to clone the repository, since this is a private repository it is not easy at all, at first you are going to need your all PAT (Personal Access Token) that you can create
following [this tutorial](https://nira.com/how-to-clone-a-private-repository-in-github/).

For saving your pat on your local pc we are going to create a local variable. First open in writing mode the file `.bashrc`
```bash
nano /home/user/.bashrc
```
where user is your account on linux. At the end of the file, write
```bash
export PAT="<pat>"
```
where <pat> is the PAT you created (for pasting `ctrl+shift+v`). Save the file by pressing `ctrl+O` and then `Enter`, and exit by pressing `ctrl+x`.

Load the variable with the following command:
```bash
source /home/user/.bashrc
```
Finally you can clone the repository:
```bash
git clone https://${PAT}@github.com/VinceNeede/Ising.git
```
Here `${PAT}` outputs the variable we saved earlier.

## Compile the file ## 
-->
First clone the repository:
```bash
git clone https://github.com/VinceNeede/Ising-Diag.git
cd Ising-Diag/
```
If not already done, set the intel environment by running 
```bash
source /opt/intel/oneapi/setvars.sh
```
It is now possible to compile the two files `Ising.f90` and `Ising_Parity.f90`, together with the library `diag.f90`, by just running
```bash
make
```
By default the output will be saved in a directory called `_results/` in the parent directory (if it does not exists it will be created). If you want to specify a different directory you can set the `OUTDIR` variable, for example if we want to save the output files in a directory called `FOO` we would write
```bash
make OUTDIR=FOO
```
If the directory is not already present, it will be created by the makefile.

You can now use the two files by setting the parameters - number of spins, transversal and longitudinal fields (the last one is not considered by Ising_Parity), OBC or PBC - in the `chain.in` file, and then run Ising with the command `./_results/Ising` or Ising_Parity with the command `._results/Ising_Parity`.
> NOTE: you don't have to compile the file each time you modify the file `chain.in`.
