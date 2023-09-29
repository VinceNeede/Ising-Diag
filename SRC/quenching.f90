module ODE
    implicit none
    
    abstract interface
    Subroutine fun_vC_vCR(z, y, x, N)
        integer, intent(in) :: N
        Real, Intent(in) :: x
        complex, intent(in) :: y(N)
        complex, intent(out) :: z(N)
    end subroutine  fun_vC_vCR
    end interface
    
    contains
    Subroutine Runge_Kutta(fun, y0, x0, dx, y, N)
        ! Compute y at point x0+dx starting from y0
        ! at point x0 given the current fun by Runge_Kutta method
        ! stat = 0 success
        ! stat = 1 y0 and y are not congruent
        ! Adolfo Avella
        ! ver:1.0 Salerno, 2013.05.10
        
        Implicit None
        procedure(fun_vC_vCR) :: fun
        integer, intent(in) :: N
        complex, intent(in) :: y0(N)
        real,intent(in) :: x0, dx
        complex, intent(out) :: y(N)
        
        
        complex :: c1(N)
        
        call fun(c1, y0, x0, N)
        y=y0+dx*c1/(6.d0,0.)
        call fun(c1, y0+dx*(0.5,0.)*c1, x0+0.5*dx, N)
        y=y+dx/(3.d0,0.)*c1
        call fun(c1, y0+dx*(0.5,0.)*c1, x0+0.5*dx, N)
        y=y+dx/(3.d0,0.)*c1
        call fun(c1, y0+dx*c1, x0+dx, N)
        
        y=y+c1*dx/(6.d0,0.)
        !y=y0+c1*cmplx(dx)
    End Subroutine Runge_Kutta
end module ODE

MODULE SystemParameters
    !!! This Library contains the variables to be imported 
    LOGICAL  :: PBC, ctrlDav, verbose
    DOUBLE PRECISION :: Lambda, Lambda0, gfield0, gfield1, theta
    integer :: steps
END MODULE SystemParameters

program dsaupd_mkl
    ! Import the library where the function dlanc is defined
    use diagonalization
    ! Import the Sparse BLAS Library
    use mkl_spblas
    ! Import the variables
    use SystemParameters

    use ODE
    
    implicit none
    integer, allocatable :: iSpin(:,:), spin_z(:)
    integer:: ell, N, nev, ncv, nnz, ii
    real :: broken_mag
    real, allocatable :: evec(:,:), E0(:)
    
    type(SPARSE_MATRIX_T) :: A, A_cmplx
    TYPE(MATRIX_DESCR) :: descrA,descrA_cmplx
    integer,allocatable :: row(:), col(:)
    real, allocatable :: val(:), magX(:), magZ(:)
    complex, allocatable :: val_cmplx(:)
    character(len=15) :: in_c

    real :: x0, dx
    complex, allocatable :: y0(:), y(:)
    !   %--------------------------------------%
    !   | Import the parameteres from the      |
    !   | file chain.in                        |
    !   %--------------------------------------%
    call get_command_argument(1,in_c)
    read(in_c,*) ell
    call get_command_argument(2,in_c)
    read(in_c,*) gfield0
    call get_command_argument(3,in_c)
    read(in_c,*) gfield1
    call get_command_argument(4,in_c)
    read(in_c,*) Lambda0
    call get_command_argument(5,in_c)
    read(in_c,*) Lambda
    call get_command_argument(6,in_c)
    read(in_c,*) theta
    call get_command_argument(7,in_c)
    read(in_c,*) PBC
    call get_command_argument(8,in_c)
    read(in_c,*) steps
    call get_command_argument(9,in_c)
    read(in_c,*) verbose




    
    ! print*, ell
    ! OPEN (Unit=1,file='chain.in',status='old')
    ! READ (1,*) ell       ! Chain length
    ! READ (1,*) gfield0    ! transverse magnetic field
    ! READ (1,*) gfield1    ! transverse magnetic field
    ! READ (1,*) Lambda0    ! starting longitudinal field
    ! READ (1,*) Lambda    ! longitudinal field
    ! READ (1,*) theta    ! transverse magnetic field
    ! READ (1,*) PBC       ! Type of boundary condirions (.true. -> PBC,   .false. -> OBC)
    ! READ (1,*) ctrlDav   ! Type of diagonalization     (.true. -> Davidson,   .false. -> Lapack full diag)
    ! read (1,*) steps    ! Number of steps for the runge-kutta
    ! read (1,*) verbose  ! .true. prints each step, .false. prints only the ending step
    ! CLOSE (unit=1)
    N=2**ell
    nnz = N+ell*2**(ell-1)           ! Interaction XX, gfield over Z, lambdafield over X
    
    allocate(row(nnz), col(nnz), val(nnz), val_cmplx(nnz))
    
    !call system_clock(count=t0)
    allocate(iSpin(N,ell), spin_z(N), magX(ell), magz(ell))
    call set_variables(ispin, spin_z, ell)
    !   %---------------------------------------%
    !   | Call the Build_Ham subroutine which   |
    !   | construct a sparse rappresentation of |
    !   | the Ising Hamiltonian in coo format.  |
    !   %---------------------------------------%
    
    call Build_Ham(A, A_cmplx, descrA,descrA_cmplx, row, col, val, val_cmplx, ell, N, nnz)
    
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
    
    nev = 1
    ncv = min(n,10)
    
    allocate(evec(N,nev), E0(nev))
    call dlanc(A, descrA, n, nev, ncv, 'I', 'SA', E0, evec=evec)
    
    allocate(y0(N), y(N))
    do ii = 1, N
        y0(ii)=cmplx(evec(ii,1))
    enddo
    x0=0.
    
    ! call fun(y,y0,0.,N)
    dx=real(ell)*theta/real(steps)
    ! call Runge_Kutta(fun, y0, x0, dx, y, N)
    ! print*,dot_product(y0,y), exp(cmplx(0.,-E0(1)*dx))
    
    if (verbose) then
        do ii =1, steps
            call Runge_Kutta(fun, y0, x0, dx, y, N)
            y=y/nrm2(y)
            call Magnetization(y, broken_mag, magX, magZ)
            y0=y
            x0=x0+dx
            call out(x0/ell)
        end do
    else
        do ii=1,steps
            call Runge_Kutta(fun, y0, x0, dx, y, N)
            y=y/nrm2(y)
            y0=y
            x0=x0+dx
        end do
        call Magnetization(y, broken_mag, magX, magZ)
        call out(x0/ell)
    endif
    
    contains
    
    ! subroutine fun2(z, y, x, alpha, N)
    !     integer, intent(in) :: N
    !     Real, Intent(in) :: x, alpha
    !     complex, intent(in) :: y(N)
    !     complex, intent(out) :: z(N)
    !     z=-alpha*y
    
    ! end subroutine fun2
    
    
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
    subroutine Build_Ham(A, A_cmplx, descrA,descrA_cmplx, row, col, val, val_cmplx, ell, N, nnz)
        
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
        type(SPARSE_MATRIX_T), intent(out) :: A, A_cmplx
        TYPE(MATRIX_DESCR), intent(out) :: descrA, descrA_cmplx
        integer, intent(out) :: row(:), col(:)
        real, intent(out) :: val(:)
        complex, intent(out) :: val_cmplx(:)
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
        
        do kk=1,N
            col(kk)=kk
            row(kk)=kk
            val(kk)=Lambda0*spin_z(kk)
            val_cmplx(kk)=cmplx(Lambda*spin_z(kk))
            do ii=1,ell-1
                val(kk)=val(kk)+2.d0*real(mod(ispin(kk,ii)+iSpin(kk,ii+1),2),8)-1.d0
                val_cmplx(kk)=val_cmplx(kk)+cmplx(2.d0*real(mod(ispin(kk,ii)+iSpin(kk,ii+1),2),8)-1.d0)
            enddo
        enddo
        if (PBC) then
            do kk = 1,N
                val(kk)=val(kk)+2.d0*real(mod(ispin(kk,ell)+iSpin(kk,1),2),8)-1.d0
                val_cmplx(kk)=val_cmplx(kk)+cmplx(2.d0*real(mod(ispin(kk,ell)+iSpin(kk,1),2),8)-1.d0)
            enddo
        endif
        start_arr=N
        
        !   %---------------------------------%
        !   |                        >        |
        !   |  -g\sum_j Z_j Transverse Field  |
        !   |                                 |
        !   %---------------------------------%
        
        do ii = 1, ell
            do jj = 1, 2**(ell-ii)
                do kk = (jj-1)*2**ii + 2**(ii-1)+1, jj*2**ii
                    indx=kk+2**(ii-1)*(1-2*ispin(kk,ii))
                    iarr=kk-2**(ii-1)-(jj-1)*2**(ii-1)+2**(ell-1)*(ii-1)
                    col(start_arr+iarr) = kk
                    row(start_arr+iarr) = indx
                    val(start_arr+iarr) = -gfield0
                    val_cmplx(start_arr+iarr) = cmplx(-gfield1)
                enddo
            enddo
        enddo
        start_arr=start_arr+ell*2**(ell-1)
        !   %-----------------------------------%
        !   |                                   |
        !   |  +h\sum_j X_j Longitudinal Field  |
        !   |                                   |
        !   %-----------------------------------%
        
        stat=mkl_sparse_d_create_coo(A,SPARSE_INDEX_BASE_ONE,N, N, nnz, row,col, val)
        if (stat .ne. 0) print*, 'error in sparse matrix creation: ', stat
        stat=mkl_sparse_z_create_coo(A_cmplx,SPARSE_INDEX_BASE_ONE,N, N, nnz, row,col, val_cmplx)
        if (stat .ne. 0) print*, 'error in sparse matrix creation Complex: ', stat
        
        descrA_cmplx % TYPE = SPARSE_MATRIX_TYPE_HERMITIAN           !!! col .geq. row
        descrA_cmplx % Mode = SPARSE_FILL_MODE_UPPER
        descrA_cmplx % diag = SPARSE_DIAG_NON_UNIT
        
        descrA % TYPE = SPARSE_MATRIX_TYPE_SYMMETRIC            !!! col .geq. row
        descrA % Mode = SPARSE_FILL_MODE_UPPER
        descrA % diag = SPARSE_DIAG_NON_UNIT
    end subroutine
    
    subroutine fun(z, y, x, N)
        integer, intent(in) :: N
        Real, Intent(in) :: x
        complex, intent(in) :: y(N)
        complex, intent(out) :: z(N)
        
        integer :: stat
        stat=mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, cmplx(0.,-1.), A_cmplx, descrA_cmplx, y,(0.,0.), z)
        if (stat .ne. 0) print*, 'error mkl_sparse_d_mv: ', stat
        
    end subroutine fun
    
    subroutine Magnetization(psi, broken_mag,magx, magz)
        implicit none
        complex, intent(in) :: psi(:)
        real(8), intent(out) ::  broken_mag, magx(:), magz(:)
        
        integer :: ii, jj
        
        broken_mag=0.d0
        do ii = 1,N
            broken_mag=broken_mag+abs(spin_z(ii))*real(conjg(psi(ii))*psi(ii))
        enddo
        do ii=1,ell
            magz(ii)=0.d0
            magx(ii)=0.d0
            do jj=1,N
                magz(ii)=magz(ii)-(1-2*iSpin(jj,ii))*real(conjg(psi(jj))*psi(jj))
                magx(ii)=magx(ii)-real(psi(jj)*conjg(psi(jj+2**(ii-1)*(1-2*ispin(jj,ii)))))
            enddo
        enddo
        broken_mag=broken_mag/ell
    end subroutine Magnetization
    
    subroutine out(t)
        implicit none
        real, intent(in) :: t
        
        if (PBC) then
            write(*,'(I4,",",I8,5(",",ES22.14)",",L2,4(",",ES22.14))') ell,steps, gfield0, gfield1, Lambda0, Lambda, theta, PBC, t, broken_mag, magz(1), magX(1)
        else
            write(*,100) ell,steps, gfield0, gfield1, Lambda0, Lambda, theta, PBC, t, broken_mag, magz, magX
            100   Format((I4,",",I8,5(",",ES22.14)",",L2,<2+2*ell>(",",ES22.14)))
        endif
        
    end subroutine out    
    
    function nrm2(psi)
        complex, intent(in) :: psi(:)
        real :: nrm2

        integer :: i

        nrm2=0.d0
        do i=1,size(psi)
            nrm2=nrm2+conjg(psi(i))*psi(i)
        enddo
        nrm2=sqrt(nrm2)
    end function nrm2
end program dsaupd_mkl

