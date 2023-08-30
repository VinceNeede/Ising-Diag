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
    integer, allocatable :: iSpin(:,:), spin_z(:), order(:)
    integer:: ell, N, nev, ncv, nnz
    real(rk) :: d(3),x(2),y(2)
    
    type(SPARSE_MATRIX_T) :: A_p, A_m
    TYPE(MATRIX_DESCR) :: descrA
    integer, allocatable :: row_p(:), col_p(:), row_m(:), col_m(:)
    real(8), allocatable :: val_p(:), val_m(:)
    
    !   %--------------------------------------%
    !   | Import the parameteres from the      |
    !   | file chain.in                        |
    !   %--------------------------------------%
    
    OPEN (Unit=1,file='chain.in',status='old')
    READ (1,*) ell       ! Chain length
    READ (1,*) gfield    ! transverse magnetic field
    READ (1,*) Lambda    ! longitudinal field
    READ (1,*) PBC       ! Type of boundary condirions (.true. -> PBC,   .false. -> OBC)
    READ (1,*) ctrlDav   ! Type of diagonalization     (.true. -> Davidson,   .false. -> Lapack full diag)
    CLOSE (unit=1)
    N=2**ell
    
    nnz = (ell-1)*2**(ell-2)+2**(ell-1)           ! Interaction XX, gfield over Z
    if (PBC) nnz=nnz+2**(ell-2)
    
    allocate(iSpin(N,ell), spin_z(N), order(N))
    allocate(row_p(nnz), row_m(nnz), col_p(nnz), col_m(nnz)&
    , val_p(nnz), val_m(nnz))
    
    call set_variables(iSpin, order, spin_z, ell)    
    
    !   %---------------------------------------%
    !   | Call the Build_Ham subroutine which   |
    !   | construct a sparse rappresentation of |
    !   | the Ising Hamiltonian in coo format.  |
    !   %---------------------------------------%
    
    call Build_Ham(A_p, A_m, descrA, row_p, col_p, val_p, row_m, col_m, val_m, ell, N, nnz)
    
    ! x=[1.d0,0.d0]
    ! print*, 'fine till here'
    ! nev=mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.d0, A_m, descrA, x,0.d0, y)
    ! if (nev .ne. 0) print*, 'error mkl_sparse_d_mv: ', nev
    
    ! print*, y
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
    call dlanc(A_p, descrA, n, nev, ncv, 'I', 'SA', d)
    
    !   %---------------------------------------%
    !   | TODO: not necessarly all the          |
    !   | eigenvalues searched are going to     |
    !   | converge, add a Flag and organize     |
    !   | the array such that d(1) is the       |
    !   | smallest.                             |
    !   %---------------------------------------%
    
    print*,d
    call dlanc(A_m, descrA, n, nev, ncv, 'I', 'SA', d)
    print*, d
    contains
    
    subroutine set_variables(ispin, order, spin_z, ell)
        implicit none
        integer, intent(in) :: ell
        integer, intent(out) :: ispin(:,:), order(:), spin_z(:)
        
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
        
        integer :: kk, tt, ii, jj, c_down, N, itemp
        
        N=2**ell
        
        iSpin = 0
        kk = 1
        tt = 2**(ell-1)+1
        do ii = 1,N
            c_down = 0
            itemp = ii-1
            do jj = 1,ell
                iSpin(ii,jj) = Mod(itemp,2)
                c_down = c_down + Mod(itemp,2)
                itemp = itemp/2
            enddo
            if(Mod(c_down,2)==0) then
                order(ii) = kk
                kk=kk+1
            else
                order(ii) = tt
                tt=tt+1
            end if
            spin_z (ii)= ell-2*c_down
        enddo
    end subroutine
    
    subroutine Build_Ham(A_p, A_m, descrA, row_p, col_p, val_p, row_m, col_m, val_m, ell, N, nnz)
        
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
        type(SPARSE_MATRIX_T), intent(out) :: A_p, A_m
        TYPE(MATRIX_DESCR), intent(out) :: descrA
        integer, intent(in) :: ell, N, nnz
        
        integer :: row_p(:), col_p(:), row_m(:), col_m(:)
        real(8):: val_p(:), val_m(:)
        
        integer :: ii, jj, kk, iarr_p, iarr_m, stat, indx, start_arr, exc
        
        
        
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
        
        iarr_p = 0
        iarr_m = 0
        do ii = 1, ell-1
            do jj = 1, 2**(ell-ii-1)
                do kk = (jj-1)*2**(ii+1)+2**ii+1, jj*2**(ii+1)
                    indx = kk+2**(ii)*(1-2*iSpin(kk,ii+1)) + 2**(ii-1)*(1-2*iSpin(kk,ii))
                    !iarr=kk-2**(ii)-(jj-1)*2**(ii)+2**(ell-1)*(ii-1)
                    !print*, indx, kk
                    
                    exc=order(kk)
                    indx=order(indx)
                    !print*, indx, exc
                    
                    if (indx<=2**(ell-1) .and. exc<=2**(ell-1)) then
                        iarr_p=iarr_p+1
                        row_p(iarr_p) = minval([indx, exc])
                        col_p(iarr_p) = maxval([indx, exc])
                        val_p(iarr_p)=-1.d0
                    else
                        iarr_m = iarr_m+1
                        row_m(iarr_m) = minval([indx, exc])-2**(ell-1)
                        col_m(iarr_m) = maxval([indx, exc])-2**(ell-1)
                        val_m(iarr_m)=-1.d0
                    endif
                enddo
            enddo
        enddo
        
        if (PBC) then
            do kk = 2**(ell-1)+1, 2**ell
                indx=kk + (1-2*ispin(kk,1)) + 2**(ell-1)*(1-2*ispin(kk,ell))
                indx = Order(indx)
                Exc = Order(kk)
                if (indx<=N/2 .and. Exc<=N/2) then
                    iarr_p=iarr_p+1
                    row_p(iarr_p) = minval([indx, Exc])
                    col_p(iarr_p) = maxval([indx, Exc])
                    val_p(iarr_p) = -1.d0
                else
                    iarr_m = iarr_m+1
                    row_m(iarr_m) = minval([indx, Exc])-N/2
                    col_m(iarr_m) = maxval([indx, Exc])-N/2
                    val_m(iarr_m) = -1.d0
                endif
            enddo
        endif
        ! start_arr=(ell-1)*2**(ell-1)
        ! if (PBC) then
        !     do kk = 2**(ell-1)+1, 2**ell
        !         iarr=kk-2**(ell-1)
        !         indx=kk + (1-2*ispin(kk,1)) + 2**(ell-1)*(1-2*ispin(kk,ell))
        !         col(start_arr+iarr)=kk
        !         row(start_arr+iarr)=indx
        !         val(start_arr+iarr)=-1.d0
        !     enddo
        !     start_arr=start_arr+2**(ell-1)
        ! endif
        
        ! !   %---------------------------------%
        ! !   |                                 |
        ! !   |  -g\sum_j Z_j Transverse Field  |
        ! !   |                                 |
        ! !   %---------------------------------%
        
        do ii = 1, N
            kk = order(ii)
            if (kk<=N/2) then
                iarr_p = iarr_p + 1
                row_p(iarr_p) = kk
                col_p(iarr_p) = kk
                val_p(iarr_p) = - gfield*spin_z(ii)
            else
                iarr_m = iarr_m + 1
                row_m(iarr_m) = kk-N/2
                col_m(iarr_m) = kk-N/2
                val_m(iarr_m) = -gfield*spin_z(ii)
            endif
        end do
        ! do ii = 1, N
        !     exc=order(ii)
        !     row(start_arr+ii) = ii
        !     col(start_arr+ii) = ii
        !     val(start_arr+ii) = -gfield*spin_z(ii)
        ! enddo
        start_arr=start_arr+N
        
        ! !   %-----------------------------------%
        ! !   |                                   |
        ! !   |  +h\sum_j X_j Longitudinal Field  |
        ! !   |                                   |
        ! !   %-----------------------------------%
        
        ! do ii = 1, ell
        !     do jj = 1, 2**(ell-ii)
        !         do kk = (jj-1)*2**ii + 2**(ii-1)+1, jj*2**ii
        !             indx=kk+2**(ii-1)*(1-2*ispin(kk,ii))
        !             iarr=kk-2**(ii-1)-(jj-1)*2**(ii-1)+2**(ell-1)*(ii-1)
        !             col(start_arr+iarr) = kk
        !             row(start_arr+iarr) = indx
        !             val(start_arr+iarr) = Lambda
        !         enddo
        !     enddo
        ! enddo
        ! print*, col_p
        ! print*, row_p
        ! print*, val_p
        ! print*, col_m
        ! print*, row_m
        ! print*, val_m
        
        stat=mkl_sparse_d_create_coo(A_p,SPARSE_INDEX_BASE_ONE,N/2,N/2, nnz, row_p,col_p, val_p)
        if (stat .ne. 0) print*, 'error in sparse matrix creation: ', stat
        
        stat=mkl_sparse_d_create_coo(A_m,SPARSE_INDEX_BASE_ONE,N/2, N/2, nnz, row_m,col_m, val_m)
        if (stat .ne. 0) print*, 'error in sparse matrix creation: ', stat
        
        descrA % TYPE = SPARSE_MATRIX_TYPE_SYMMETRIC            !!! col .geq. row
        descrA % Mode = SPARSE_FILL_MODE_UPPER
        descrA % diag = SPARSE_DIAG_NON_UNIT
        
    end subroutine
end program dsaupd_mkl

