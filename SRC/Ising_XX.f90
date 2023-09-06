MODULE SystemParameters
    !!! This Library contains the variables to be imported 
    LOGICAL  :: PBC, ctrlDav
    DOUBLE PRECISION :: Lambda, gfield
END MODULE SystemParameters

program dsaupd_mkl
    ! Import the library where the function dlanc is defined
    use diagonalization
    ! Import the Sparse BLAS Library
    use mkl_spblas
    ! Import the variables
    use SystemParameters

    implicit none
    integer, allocatable :: iSpin(:,:), spin_z(:)
    integer:: ell, N, nev, ncv, nnz
    real(rk) :: d(3)

    type(SPARSE_MATRIX_T) :: A
    TYPE(MATRIX_DESCR) :: descrA
    integer,allocatable :: row(:), col(:)
    real(8), allocatable :: val(:)

    !   %--------------------------------------%
    !   | Import the parameteres from the      |
    !   | file chain.in                        |
    !   %--------------------------------------%

    write(*,*) 'WARNING: this version is obsolate'


    OPEN (Unit=1,file='chain.in',status='old')
    READ (1,*) ell       ! Chain length
    READ (1,*) gfield    ! transverse magnetic field
    READ (1,*) Lambda    ! longitudinal field
    READ (1,*) PBC       ! Type of boundary condirions (.true. -> PBC,   .false. -> OBC)
    READ (1,*) ctrlDav   ! Type of diagonalization     (.true. -> Davidson,   .false. -> Lapack full diag)
    CLOSE (unit=1)
    N=2**ell
    nnz = (ell-1)*2**(ell-1)+N+ell*2**(ell-1)           ! Interaction XX, gfield over Z, lambdafield over X
    if (PBC) nnz=nnz+2**(ell-1)
    
    allocate(row(nnz), col(nnz), val(nnz))
    
    !call system_clock(count=t0)
    allocate(iSpin(N,ell), spin_z(N))

    !   %---------------------------------------%
    !   | Call the set_variables subroutine     |
    !   | which construct the matrix ispin, used|
    !   | for decimal to binary representation, |
    !   | and spin_z which contains the total   |
    !   | spin along the z axis for a given     |
    !   | basis vector.                         |
    !   %---------------------------------------%

    call set_variables(ispin, spin_z, ell)

    !   %---------------------------------------%
    !   | Call the Build_Ham subroutine which   |
    !   | construct a sparse rappresentation of |
    !   | the Ising Hamiltonian in coo format.  |
    !   %---------------------------------------%

    call Build_Ham(A, descrA, row, col, val, ell, N, nnz)

    !   %---------------------------------------%
    !   | Call the dlanc subroutine which       |
    !   | computes the eigenvalues, and         |
    !   | eigenvectors (if requested), of the   |
    !   | sparse matrix. NEV is the number of   |
    !   | eigenvalues to search for. NCV is the |
    !   | size of the Krylov sub-space which    |
    !   | must obey to the condition:           |
    !   |       NEV + 1 <= NCV <= MAXNCV        |
    !   | Where MAXNCV is set in the            |
    !   | diagonalization module.               |        
    !   | 'I' specifies that the problem to     |
    !   | solve is a standard Eigenvalues       |
    !   | problem, i.e.                         |
    !   |           A*x=lambda*B*x              |
    !   | where x is a vector, lambda the       |
    !   | eigenvalue to search for, and B=I.    |
    !   | 'SA' specifies that we are searching  |
    !   | for the Smallest Algebric eigenvalues.|
    !   | Other options can be select, as 'LM'  |
    !   | which stands for Larges Magnitude.    |
    !   %---------------------------------------%

    nev = 3
    ncv = min(n,20)
    call dlanc(A, descrA, n, nev, ncv, 'I', 'SA', d)

    !   %---------------------------------------%
    !   | TODO: not necessarly all the          |
    !   | eigenvalues searched are going to     |
    !   | converge, add a Flag and organize     |
    !   | the array such that d(1) is the       |
    !   | smallest.                             |
    !   %---------------------------------------%

    print*,d

    contains

    subroutine set_variables(ispin, spin_z, ell)
        implicit none
        integer, intent(in) :: ell
        integer, intent(out) :: ispin(:,:), spin_z(:)

        !   %---------------------------------------%
        !   | Construct the matrix iSpin(N,L) and   |
        !   | the array spin_z.                     |
        !   | The i-th row of the matrix iSpin is   |
        !   | the i-th vector of the computational  |
        !   | basis.                                |
        !   |           |00>  ---> 1                |
        !   |           |10>  ---> 2                |
        !   |           |01>  ---> 3                |
        !   |           |11>  ---> 4                |
        !   | The j-th row of iSpin is the value of |
        !   | the j-th spin. (0 for up, 1 for down).|
        !   | The i-th elent of spin_z is the total |
        !   | spin along the z axis of the i-th     |
        !   | vector of the computational basis.    |
        !   | That is done by couting the number of |
        !   | down spins c_down, so that the total  | 
        !   | spin is L-2*c_down, since the number  |
        !   | of up spins is L-c_down.              |
        !   %---------------------------------------%

        integer :: ii, jj, c_down, NumTot, itemp
        NumTot=2**ell

        iSpin = 0
        do ii = 1,NumTot
            c_down = 0
            itemp = ii-1
            do jj = 1,ell
                iSpin(ii,jj) = Mod(itemp,2)
                c_down = c_down + Mod(itemp,2)
                itemp = itemp/2
            enddo
            spin_z (ii)= ell-2*c_down
        enddo
    end subroutine

    subroutine Build_Ham(A, descrA, row, col, val, ell, N, nnz)

    !   %---------------------------------------%
    !   | Build the Sparse Matrix in COO Format.|
    !   | The valeus of the non zero elements   |
    !   | (in number nnz) are stored in an array|
    !   | val, the i-th element of that array   |
    !   | occupies the row(i), col(i) position  |
    !   | of the matrix. mkl_spblas allows other|
    !   | formats. Furthermore it can be        |
    !   | specified if the matrix is symmetric  |
    !   | (or hermitian) in that way only the   |
    !   | upper (or lower) triangular part can  |
    !   | be saved.                             |
    !   %---------------------------------------%

        implicit none
        type(SPARSE_MATRIX_T), intent(out) :: A
        TYPE(MATRIX_DESCR), intent(out) :: descrA
        integer, intent(out) :: row(:), col(:)
        real(8), intent(out) :: val(:)
        integer, intent(in) :: ell, N, nnz

        integer :: ii, jj, kk, iarr, stat, indx, start_arr


        !   %---------------------------------------%
        !   |                                       |
        !   |      -J\sum_i X_i X_{i+1} Coupling    |
        !   |                                       |
        !   %---------------------------------------%


        !   %---------------------------------------%
        !   | Each element X_ii X_{ii+1} of the sum |
        !   | is carried out separately, since each |
        !   | terms acts different, i.e. if         |
        !   |    <ss|X_ii X_{ii+1}|kk> .neq. 0      |
        !   | then all the other terms will be zero.|
        !   | The only non zero term for a give ii  |
        !   | is the term:                          |
        !   |       <indx|X_ii X_{ii+1}|kk>         |
        !   | furthermore since we can fill only    |
        !   | the upper triangular of the matrix    |
        !   | by imposing:                          |
        !   |               kk .geq. indx           |
        !   | which is equivalent to impose that    |
        !   | the spin value of the ii+1 spin must  |
        !   | always be 1. Such thing is implemented|
        !   | by the further loop over jj.          |
        !   | The iarr variable is an index that    |
        !   | goes from 1 to nnz, that could be     |
        !   | implemented by setting it to zero and |
        !   | adding one (iarr=iarr+1) at each      |
        !   | assigment. Since such implementation  |
        !   | can create some problems in parallel  |
        !   | programming, by preferring            |
        !   | optimization over readability I       |
        !   | decided to adopt another way.         |
        !   %---------------------------------------%

        start_arr=0
        do ii = 1, ell-1
            do jj = 1, 2**(ell-ii-1)
                do kk = (jj-1)*2**(ii+1)+2**ii+1, jj*2**(ii+1)
                    indx = kk+2**(ii)*(1-2*iSpin(kk,ii+1)) + 2**(ii-1)*(1-2*iSpin(kk,ii))
                    iarr=kk-2**(ii)-(jj-1)*2**(ii)+2**(ell-1)*(ii-1)
                    col(iarr)=kk              !!! kk>indx     
                    row(iarr)=indx
                    val(iarr)=-1.d0
                enddo
            enddo
        enddo
        start_arr=(ell-1)*2**(ell-1)
        if (PBC) then
            do kk = 2**(ell-1)+1, 2**ell
                iarr=kk-2**(ell-1)
                indx=kk + (1-2*ispin(kk,1)) + 2**(ell-1)*(1-2*ispin(kk,ell))
                col(start_arr+iarr)=kk
                row(start_arr+iarr)=indx
                val(start_arr+iarr)=-1.d0
            enddo
            start_arr=start_arr+2**(ell-1)
        endif

        !   %---------------------------------%
        !   |                                 |
        !   |  -g\sum_j Z_j Transverse Field  |
        !   |                                 |
        !   %---------------------------------%

        do ii = 1, N
            row(start_arr+ii) = ii
            col(start_arr+ii) = ii
            val(start_arr+ii) = -gfield*spin_z(ii)
        enddo
        start_arr=start_arr+N
                
        !   %-----------------------------------%
        !   |                                   |
        !   |  +h\sum_j X_j Longitudinal Field  |
        !   |                                   |
        !   %-----------------------------------%

        do ii = 1, ell
            do jj = 1, 2**(ell-ii)
                do kk = (jj-1)*2**ii + 2**(ii-1)+1, jj*2**ii
                    indx=kk+2**(ii-1)*(1-2*ispin(kk,ii))
                    iarr=kk-2**(ii-1)-(jj-1)*2**(ii-1)+2**(ell-1)*(ii-1)
                    col(start_arr+iarr) = kk
                    row(start_arr+iarr) = indx
                    val(start_arr+iarr) = Lambda
                enddo
            enddo
        enddo


        stat=mkl_sparse_d_create_coo(A,SPARSE_INDEX_BASE_ONE,N, N, nnz, row,col, val)
        if (stat .ne. 0) print*, 'error in sparse matrix creation: ', stat

        descrA % TYPE = SPARSE_MATRIX_TYPE_SYMMETRIC            !!! col .geq. row
        descrA % Mode = SPARSE_FILL_MODE_UPPER
        descrA % diag = SPARSE_DIAG_NON_UNIT
    end subroutine
end program dsaupd_mkl

