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
    integer:: ell, N, nev, ncv, nnz, ii
    real(8) :: d_m(3), d_p(3)
    real(8), allocatable :: evec_m(:,:), evec_p(:,:)
    
    type(SPARSE_MATRIX_T) :: A_p, A_m
    TYPE(MATRIX_DESCR) :: descrA
    integer, allocatable :: row_p(:), col_p(:), row_m(:), col_m(:)
    real(8), allocatable :: val_p(:), val_m(:), mag_x(:), mag_Z0(:), mag_z1(:)
    
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
    !   | If evec is given as an argument to the|
    !   | routine, on the output it will store  |
    !   | the eigenvectors computed. It has to  |
    !   | a NxNEV matrx. The eigenvectors are   |
    !   | stored as columns, and the order is   |
    !   | the same of the eigenvalues.          |
    !   %---------------------------------------%
    
    nev = 3
    ncv = min(n/2,20)
    
    allocate(evec_p(N,nev), evec_m(N,nev), mag_x(ell), mag_z0(ell), mag_z1(ell))
    
    call dlanc(A_p, descrA, n/2, nev, ncv, 'I', 'SA', d_p, evec=evec_p)
    
    
    ! do ii=1, N
    !     if (order(ii) .gt. N/2) then
    !         print '(3(F9.2,x))', 0.d0,0.d0,0.d0
    !     else
    !         print'(3(F9.2,x))', evec_p(order(ii),:)
    !     endif
    ! enddo
    ! print*,''
    call reorder(d_p, evec_p, N, nev)
    !print*, d_p
    ! do ii=1, N
    !     print'(3(F9.2,x))', evec_p(order(ii),:)
    ! enddo
    ! print*,''
    
    call dlanc(A_m, descrA, n/2, nev, ncv, 'I', 'SA', d_m, evec=evec_m)
    ! do ii=1, N
    !     if (order(ii) .le. N/2) then 
    !         print'(3(F9.2,x))', evec_m(order(ii)+N/2,:)
    !     else
    !         print'(3(F9.2,x))', evec_m(order(ii)-N/2,:)
    !     endif
    ! enddo
    ! print*,''
    call reorder(d_m, evec_m, N, nev)
    !print*, d_m
    ! do ii=1, N
    !     if (order(ii) .le. N/2) then 
    !         print'(3(F9.2,x))', evec_m(order(ii)+N/2,:)
    !     else
    !         print'(3(F9.2,x))', evec_m(order(ii)-N/2,:)
    !     endif
    ! enddo
    ! print*,''
    
    call Magnetization(evec_p(:,1),evec_m(:,1), mag_x, mag_Z0, mag_z1)
    !print*, mag_x
    
    call out()
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
        !   | In order to seperate the two parity   |
        !   | sector we have to consider a permuted |
        !   | basis, take L=2 for example:          |
        !   |               1       1               |
        !   |               2       4               |
        !   |               3       2               |
        !   |               4       3               |
        !   | For this reason we introduce a new    |
        !   | array order of size N, which idea is: |
        !   | Calling A the matrix in the           |
        !   | computational basis, and B in the     |
        !   | permuted one, we would have:          |
        !   |                                       |
        !   |  (A)_{i,j} = (B)_{order(i), order(j)} |
        !   |                                       |
        !   | for order(i)<=N/2 we would have a     |
        !   | symmetric parity sector, for          |
        !   | order(i)>N/2 the antisymmetric one.   |
        !   | The array order can be build by takin |
        !   | the array of the permuted basis and   |
        !   | exchange the values with the array    |
        !   | index, i.e.:                          |
        !   |               order(1)=1              |  
        !   |               order(4)=2              |
        !   |               order(2)=3              |
        !   |               order(3)=4              |
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
        !   | The two parity sectors are implemented|
        !   | simultaneously, for the minus parity  |
        !   | sector, since the indexes hav to go   |
        !   | from 1 to N/2, we have to subtract N/2|
        !   %---------------------------------------%
        
        iarr_p = 0
        iarr_m = 0
        do ii = 1, ell-1
            do jj = 1, 2**(ell-ii-1)
                do kk = (jj-1)*2**(ii+1)+2**ii+1, jj*2**(ii+1)
                    indx = kk+2**(ii)*(1-2*iSpin(kk,ii+1)) + 2**(ii-1)*(1-2*iSpin(kk,ii))
                    exc=order(kk)
                    indx=order(indx)
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
        
        stat=mkl_sparse_d_create_coo(A_p,SPARSE_INDEX_BASE_ONE,N/2,N/2, nnz, row_p,col_p, val_p)
        if (stat .ne. 0) print*, 'error in sparse matrix creation: ', stat
        
        stat=mkl_sparse_d_create_coo(A_m,SPARSE_INDEX_BASE_ONE,N/2, N/2, nnz, row_m,col_m, val_m)
        if (stat .ne. 0) print*, 'error in sparse matrix creation: ', stat
        
        descrA % TYPE = SPARSE_MATRIX_TYPE_SYMMETRIC            !!! col .geq. row
        descrA % Mode = SPARSE_FILL_MODE_UPPER
        descrA % diag = SPARSE_DIAG_NON_UNIT
        
    end subroutine
    
    subroutine Magnetization(psi_0, psi_1, mag_x, mag_z0, mag_z1)
        implicit none
        real(8), intent(in) :: psi_0(:)
        real(8), intent(in) :: psi_1(:)
        real(8), intent(out) :: mag_x(:), mag_z0(:), mag_z1(:)
        
        !   %---------------------------------------%
        !   | Here we compute the magnetization of  |
        !   | ground state by computing:            |
        !   |           <psi_0|Mx|psi_1>            |
        !   | where |psi_0> is the ground state with|
        !   | parity +, and psi_1 the ground state  |
        !   | with parity -, since Mx anticommutes  |
        !   | withe the parity operator, his mean   |
        !   | value on psi_0 and psi_1 alone is     |
        !   | zero.                                 |
        !   | In doing the computation we have to   |
        !   | remind that the states are not written|
        !   | in the computational basis but in the |
        !   | permuted one, so we have to permute it|
        !   | back.                                 |
        !   %---------------------------------------%
        
        integer :: ii, jj, jord, exc
        
        do ii =1, ell
            mag_x(ii)=0.d0
            mag_z0(ii)=0.d0
            mag_z1(ii)=0.d0
            do jj =1, N
                jord=order(jj)
                exc=order(jj+2**(ii-1)*(1-2*ispin(jj,ii)))
                if (jord .le. N/2) then
                    mag_Z0(ii)=mag_Z0(ii) + (1-2*iSpin(jj,ii))*abs(psi_0(jord))**2
                    if (exc .gt. N/2) then
                        mag_x(ii)=mag_x(ii)+psi_0(jord)*psi_1(exc-N/2)
                    endif
                else
                    mag_Z1(ii)=mag_Z1(ii) + (1-2*iSpin(jj,ii))*abs(psi_1(jord-N/2))**2
                endif
            enddo
        enddo
    end subroutine
    
    subroutine out()
        implicit none
        ! logical :: exist
        
        ! inquire(file="chain.out", exist=exist)
        ! if (exist) then
        !   open(12, file="chain.out", status="old", position="append", action="write")
        ! else
        !   open(12, file="chain.out", status="new", action="write")
        !   write(12,'(A44)') "L,g,E_p0,E_p1,E_p2,E_m0,E_m1,E_m2,Mx"
        ! end if
        if (PBC) then
            write(*,'(I4,10(",",ES21.14))') ell, gfield, d_p, d_m, mag_x(1), mag_z0(1),mag_z1(1)
        else
            write(*,100) ell, gfield, d_p, d_m, mag_x, mag_z0, mag_z1
100   Format(I4,<7+3*ell>(",",ES21.14))
        endif
        
    end subroutine
    
    subroutine reorder(d, evec, N, nev)
        implicit none
        integer, intent(in) :: N, nev
        real(8), intent(inout) :: d(:)
        real(8), intent(inout) :: evec(:,:)
        
        real(8) :: rtemp, atemp(N)
        logical :: exc
        integer :: i
        
        exc=.true.
        do while (exc)      ! exc signals if along the array there was at least one swap
            exc=.false.
            do i = 1, nev-1
                if (d(i)>d(i+1)) then
                    call dswap(d(i), d(i+1))
                    call aswap(evec(:,i), evec(:,i+1),N)
                    exc=.true.
                endif
            enddo
        enddo
        
    end subroutine
    
    subroutine dswap (x,y)
        implicit none
        real(8),intent(inout) :: x,y
        
        real(8) :: c
        c=x
        x=y
        y=c
    end subroutine
    subroutine aswap (x,y,n)
        implicit none
        real(8),intent(inout) :: x(:),y(:)
        integer, intent(in) :: n
        
        real(8) :: c(n)
        c=x
        x=y
        y=c
    end subroutine
end program dsaupd_mkl

