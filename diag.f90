!! ifx -warn all diag_ex_arpack.f90 -lblas -llapack -larpack

module diagonalization
    !   USE http://rsusu1.rnd.runnet.ru/libraries/ARPACK/node136.html#SECTION001210000000000000000 AS DOCUMENTATION
    !
    !     This example program is intended to illustrate the
    !     simplest case of using ARPACK in considerable detail.
    !     This code may be used to understand basic usage of ARPACK
    !     and as a template for creating an interface to ARPACK.
    !
    !     This code shows how to use ARPACK to find a few eigenvalues
    !     (lambda) and corresponding eigenvectors (x) for the standard
    !     eigenvalue problem:
    !
    !                        A*x = lambda*x
    !
    !     where A is an n by n real symmetric matrix.
    !
    !     The main points illustrated here are
    !
    !        1) How to declare sufficient memory to find NEV
    !           eigenvalues of largest magnitude.  Other options
    !           are available.
    !
    !        2) Illustration of the reverse communication interface
    !           needed to utilize the top level ARPACK routine DSAUPD
    !           that computes the quantities needed to construct
    !           the desired eigenvalues and eigenvectors(if requested).
    !
    !        3) How to extract the desired eigenvalues and eigenvectors
    !           using the ARPACK routine DSEUPD.
    !
    !     The only thing that must be supplied in order to use this
    !     routine on your problem is to change the array dimensions
    !     appropriately, to specify WHICH eigenvalues you want to compute
    !     and to supply a matrix-vector product
    !
    !                         w <-  Av
    !
    !     in place of the call to AV( ) below.
    !
    !     Once usage of this routine is understood, you may wish to explore
    !     the other available options to improve convergence, to solve generalized
    !     problems, etc.  Look at the file ex-sym.doc in DOCUMENTS directory.
    !     This codes implements
    !
    !\Example-1
    !     ... Suppose we want to solve A*x = lambda*x in regular mode.
    !     ... OP = A  and  B = I.
    !     ... Assume "call av (n,x,y)" computes y = A*x
    !     ... Use mode 1 of DSAUPD.
    !
    !\BeginLib
    !
    !\Routines called:
    !     dsaupd  ARPACK reverse communication interface routine.
    !     dseupd  ARPACK routine that returns Ritz values and (optionally)
    !             Ritz vectors.
    !     dnrm2   Level 1 BLAS that computes the norm of a vector.
    !     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
    !
    !\Author
    !     Richard Lehoucq
    !     Danny Sorensen
    !     Chao Yang
    !     Dept. of Computational &
    !     Applied Mathematics
    !     Rice University
    !     Houston, Texas
    !
    !
    !\Remarks
    !     1. None
    !
    !\EndLib
    !
    !-----------------------------------------------------------------------
    !
    !     %------------------------------------------------------%
    !     | Storage Declarations:                                |
    !     |                                                      |
    !     | The maximum dimensions for all arrays are            |
    !     | set here to accommodate a problem size of            |
    !     | N .le. MAXN                                          |
    !     |                                                      |
    !     | NEV is the number of eigenvalues requested.          |
    !     |     See specifications for ARPACK usage below.       |
    !     |                                                      |
    !     | NCV is the largest number of basis vectors that will |
    !     |     be used in the Implicitly Restarted Arnoldi      |
    !     |     Process.  Work per major iteration is            |
    !     |     proportional to N*NCV*NCV.                       |
    !     |                                                      |
    !     | You must set:                                        |
    !     |                                                      |
    !     | MAXN:   Maximum dimension of the A allowed.          |
    !     | MAXNEV: Maximum NEV allowed.                         |
    !     | MAXNCV: Maximum NCV allowed.                         |
    !     %------------------------------------------------------%
    !
    implicit none
    integer, parameter :: max_ell=20, maxn=2**max_ell, maxnev=10, maxncv=25, ldv=maxn
    integer, parameter :: rk=8
    !
    !     %--------------%
    !     | Local Arrays |
    !     %--------------%
    !
    Real(rk) :: v(ldv,maxncv),workl(maxncv*(maxncv+8)),workd(3*maxn), resid(maxn),ax(maxn)
    logical :: select(maxncv)
    integer :: iparam(11), ipntr(11)
    !
    !     %---------------%
    !     | Local Scalars |
    !     %---------------%
    !
    integer :: ido, lworkl, info, ierr, j, ishfts, maxitr, mode1, nconv, itry
    Real(rk) :: tol, sigma, rnorm
    !
    !     %------------%
    !     | Parameters |
    !     %------------%
    !
    Real(rk), parameter :: zero=0.0d+0        !! Set real kind
    !
    !     %-----------------------------%
    !     | BLAS & LAPACK routines used |
    !     %-----------------------------%
    !
    Real(rk)&
    &      dnrm2
    external dnrm2, daxpy, dseupd, dsaupd, dcopy, dgetv0
    !
    !
    !     %-------------------------------------------------%
    !     | The following sets dimensions for this problem. |
    !     %-------------------------------------------------%
    !
    
    abstract interface
        subroutine mvmul(n,x,y)
            import :: rk
            integer,intent(in) :: n
            real(rk),intent(in) :: x(n)
            real(rk),intent(out) :: y(n)

        end subroutine
    end interface
    
    contains
    
    subroutine dlanc(A, descrA,n,nev,ncv,bmat,which,d,evec,v0)
        use mkl_spblas
        implicit none
        type(SPARSE_MATRIX_T),intent(in),target :: A
        TYPE(MATRIX_DESCR) :: descrA
        integer,intent(in) :: n,nev,ncv
        character,intent(in) :: bmat(1), which(2)
        real(8),intent(out) :: d(nev)
        real(8),intent(out),optional :: evec(n,nev)
        real(8),intent(in),optional :: v0(n)
        logical :: rvec = .false.
        integer :: stat=0
        ido = 0
        lworkl = ncv*(ncv+8)
        tol = zero            !! if zero, machine precision is used
        if (present(v0)) then
            info=1
            resid(1:n)=v0
        else
            info=0
        end if

        if (present(evec)) rvec=.true.
        !!! Parameters to set
        !
        !     %-----------------------------------------------%
        !     |                                               |
        !     | Specifications for ARPACK usage are set       |
        !     | below:                                        |
        !     |                                               |
        !     |    1) NEV = 4  asks for 4 eigenvalues to be   |
        !     |       computed.                               |
        !     |                                               |
        !     |    2) NCV = 20 sets the length of the Arnoldi |
        !     |       factorization                           |
        !     |                                               |
        !     |    3) This is a standard problem              |
        !     |         (indicated by bmat  = 'I')            |
        !     |                                               |
        !     |    4) Ask for the NEV eigenvalues of          |
        !     |       largest magnitude                       |
        !     |         (indicated by which = 'LM')           |
        !     |       See documentation in DSAUPD for the     |
        !     |       other options SM, LA, SA, LI, SI.       |
        !     |                                               |
        !     | Note: NEV and NCV must satisfy the following  |
        !     | conditions:                                   |
        !     |              NEV <= MAXNEV                    |
        !     |          NEV + 1 <= NCV <= MAXNCV             |
        !     %-----------------------------------------------%
        !
        
        !
        if ( n .gt. maxn ) then
            stop ' ERROR with _DSIMP: N is greater than MAXN '
        else if ( nev .gt. maxnev ) then
            stop ' ERROR with _DSIMP: NEV is greater than MAXNEV '
        else if ( ncv .gt. maxncv ) then
            stop ' ERROR with _DSIMP: NCV is greater than MAXNCV '
        end if
        !
        !     %-----------------------------------------------------%
        !     |                                                     |
        !     | Specification of stopping rules and initial         |
        !     | conditions before calling SSAUPD                    |
        !     |                                                     |
        !     | TOL  determines the stopping criterion.             |
        !     |                                                     |
        !     |      Expect                                         |
        !     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
        !     |               computed   true                       |
        !     |                                                     |
        !     |      If TOL .le. 0,  then TOL <- macheps            |
        !     |           (machine precision) is used.              |
        !     |                                                     |
        !     | IDO  is the REVERSE COMMUNICATION parameter         |
        !     |      used to specify actions to be taken on return  |
        !     |      from SSAUPD. (See usage below.)                |
        !     |                                                     |
        !     |      It MUST initially be set to 0 before the first |
        !     |      call to SSAUPD.                                |
        !     |                                                     |
        !     | INFO on entry specifies starting vector information |
        !     |      and on return indicates error codes            |
        !     |                                                     |
        !     |      Initially, setting INFO=0 indicates that a     |
        !     |      random starting vector is requested to         |
        !     |      start the ARNOLDI iteration.  Setting INFO to  |
        !     |      a nonzero value on the initial call is used    |
        !     |      if you want to specify your own starting       |
        !     |      vector (This vector must be placed in RESID.)  |
        !     |                                                     |
        !     | The work array WORKL is used in SSAUPD as           |
        !     | workspace.  Its dimension LWORKL is set as          |
        !     | illustrated below.                                  |
        !     |                                                     |
        !     %-----------------------------------------------------%
        !
        
        !
        !     %---------------------------------------------------%
        !     | Specification of Algorithm Mode:                  |
        !     |                                                   |
        !     | This program uses the exact shift strategy        |
        !     | (indicated by setting PARAM(1) = 1).              |
        !     | IPARAM(3) specifies the maximum number of Arnoldi |
        !     | iterations allowed.  Mode 1 of SSAUPD is used     |
        !     | (IPARAM(7) = 1). All these options can be changed |
        !     | by the user. For details see the documentation in |
        !     | SSAUPD.                                           |
        !     %---------------------------------------------------%
        !
        ishfts = 1
        maxitr = 1000
        mode1 = 1
        !
        iparam(1) = ishfts
        !
        iparam(3) = maxitr
        !
        iparam(7) = mode1
        !
        !     %------------------------------------------------%
        !     | M A I N   L O O P (Reverse communication loop) |
        !     %------------------------------------------------%
        !
        !     %------------------------------------------------%
        !     | Repeadtly call dvget0 in order to obtain the   |
        !     | starting vector. THIIS IS ACTUALLY DONE BY     |
        !     |   DSAUPD VIA DSAUPD2                           |
        !     %------------------------------------------------%
        !
        
        ! itry=0
        ! do while (ido/=99)
        !     itry=itry+1
        !     call dgetv0(ido, 'I', itry, .false., n,1, v, ldv, resid, rnorm, ipntr, workd, info)
        !     if (info<0) stop "Couldn't create starting vector"
        !     if (ido==-1) call av(n,workd(ipntr(1)),workd(ipntr(2)))
        
        ! end do
        ! itry=0
        ido=0
        !        %---------------------------------------------%
        !        | Repeatedly call the routine DSAUPD and take |
        !        | actions indicated by parameter IDO until    |
        !        | either convergence is indicated or maxitr   |
        !        | has been exceeded.                          |
        !        %---------------------------------------------%
        !
        do while (ido/= 99)
            !print*, info
            call dsaupd ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl,&
            lworkl, info )
            !
            !print*, ido
            if (ido .eq. -1 .or. ido .eq. 1) then
                !
                !           %--------------------------------------%
                !           | Perform matrix vector multiplication |
                !           |              y <--- OP*x             |
                !           | The user should supply his/her own   |
                !           | matrix vector multiplication routine |
                !           | here that takes workd(ipntr(1)) as   |
                !           | the input, and return the result to  |
                !           | workd(ipntr(2)).                     |
                !           %--------------------------------------%
                !
                stat=mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.d0, A, descrA, workd(ipntr(1)),0.d0, workd(ipntr(2)))
                if (stat .ne. 0) print*, 'error mkl_sparse_d_mv: ', stat
                !
                !           %-----------------------------------------%
                !           | L O O P   B A C K to call DSAUPD again. |
                !           %-----------------------------------------%
                !
            end if
        end do
        !
        !     %----------------------------------------%
        !     | Either we have convergence or there is |
        !     | an error.                              |
        !     %----------------------------------------%
        !
        if ( info .lt. 0 ) then
            print*, 'error in diagonalization dsaupd', info
            stop
        endif
        !
        !        %-------------------------------------------%
        !        | No fatal errors occurred.                 |
        !        | Post-Process using DSEUPD.                |
        !        |                                           |
        !        | Computed eigenvalues may be extracted.    |
        !        |                                           |
        !        | Eigenvectors may be also computed now if  |
        !        | desired.  (indicated by rvec = .true.)    |
        !        |                                           |
        !        | The routine DSEUPD now called to do this  |
        !        | post processing (Other modes may require  |
        !        | more complicated post processing than     |
        !        | mode1.)                                   |
        !        |                                           |
        !        %-------------------------------------------%
        !
        !                                       !! ALL computes NEV eigenvector, "S" computes only the ones specified by the array select
        !! The firs v contains the vector B-Orthonormal in case of Generalized eigenvalues
        call dseupd ( rvec, 'All', select, d, v, ldv, sigma,&
        &         bmat, n, which, nev, tol, resid, ncv, v, ldv,&
        &         iparam, ipntr, workd, workl, lworkl, ierr )
        !
        !         %----------------------------------------------%
        !         | Eigenvalues are returned in the first column |
        !         | of the two dimensional array D and the       |
        !         | corresponding eigenvectors are returned in   |
        !         | the first NCONV (=IPARAM(5)) columns of the  |
        !         | two dimensional array V if requested.        |
        !         | Otherwise, an orthogonal basis for the       |
        !         | invariant subspace corresponding to the      |
        !         | eigenvalues in D is returned in V.           |
        !         %----------------------------------------------%
        !
        if ( ierr .ne. 0) then
            print*, 'error in diagonalization dseupd', ierr
            stop
        endif
        !
        nconv =  iparam(5)
        do j = 1, (nconv+1)/2
            call swap(d(j), d(nconv-(j-1)))
        enddo
        if (rvec) then
            do j=1,nconv
                call dcopy(n,v(:,nconv-(j-1)),1,evec(:,j),1)
            end do
        endif
    end subroutine

    subroutine swap (x,y)
        implicit none
        real(8),intent(inout) :: x,y

        real(8) :: c
        c=x
        x=y
        y=c
    end subroutine

end module

! program prova
!     use diagonalization
!     implicit none
!     real(8) :: d(3)


!     call dlanc(av,maxn,3,10,'I','SA',d)

!     print*,d

    

!     contains

!         ! ------------------------------------------------------------------
!     !     matrix vector subroutine
!     subroutine av (n, v, w)
!         integer, intent(in) :: n
!         real(8), intent(in) :: v(n)
!         real(8), intent(out) :: w(n)
        
!         w(1)=v(2)
!         w(2)=v(1)
!         w(3)=0.5d0*v(3)
!         return
        
!     end subroutine
    
! end program prova