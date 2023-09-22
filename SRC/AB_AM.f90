program ode_adams
    implicit none
    integer (4), parameter :: N=1
    real(8),parameter :: b0=55.d0/24.d0, b1=-59.d0/24.d0
    real(8),parameter :: b2=37.d0/24.d0, b3=-9.d0/24.d0
    real(8),parameter :: b_1=9.d0/24.d0, b0_=19.d0/24.d0
    real(8),parameter :: b1_=-5.d0/24.d0, b2_=1.d0/24.d0
    real(8), pointer, dimension(:) :: d,y0,y,c0,c1,c2,c3
    real(8) :: x0, xf, dx
    integer :: steps, i, j
    
    allocate(y0(N),y(N),c0(N),c1(N),c2(N),c3(N))
    x0=1.d0
    xf=2.d0
    steps=10000
    
    y0=[exp(1.d0)]
    y=y0
    print '(3(ES21.14,x))',x0,abs((y-sol([1.d0,1.d0],x0))/sol([1.d0,1.d0],x0))
    
    
    dx=1.e-10!(xf-x0)/real(steps)
    ! Euler Step
    c0=func(y0,x0,N)
    y=y0+dx*c0
    x0=x0+dx
    print '(3(ES21.14,x))',x0,abs((y-sol([1.d0,1.d0],x0))/sol([1.d0,1.d0],x0))
    
    ! AB(p=1) + AM(p=0)
    d=>y0
    y0=>y
    y=>d
    c1=func(y0,x0,N)
    y=y0+dx*(3.d0*c1-c0)*0.5d0
    c2=func(y,x0+dx,N)
    y=y0+dx*(c1+c2)*0.5d0
    x0=x0+dx
    print '(3(ES21.14,x))',x0,abs((y-sol([1.d0,1.d0],x0))/sol([1.d0,1.d0],x0))
    
    ! AB(p=2) + AM(p=1)
    d=>y0
    y0=>y
    y=>d
    c2=func(y0,x0,N)
    y=y0+dx*(23.d0/12.d0*c2-16.d0/12.d0*c1+5.d0/12.d0*c0)
    c3=func(y,x0+dx,N)
    y=y0+dx*(5.d0/12.d0*c3+8.d0/12.d0*c2-1.d0/12.d0*c1)
    x0=x0+dx
    print '(3(ES21.14,x))',x0,abs((y-sol([1.d0,1.d0],x0))/sol([1.d0,1.d0],x0))
    
    dx=(xf-x0)/real(steps)
    
    do i=1,steps-3
        d=>y0
        y0=>y
        y=>d
        
        d=>c0
        c0=>c1
        c1=>c2
        c2=>c3
        c3=>d
        
        c2=func(y0,x0,N)
        y=y0+dx*(b0*c2+b1*c1+b2*c0+b3*c3)
        do j=1,2
            c3=func(y,x0+dx,N)
            y=y0+dx*(b_1*c3+b0_*c2+b1_*c1+b2_*c0)    
        enddo
        x0=x0+dx
        print '(3(ES21.14,x))',x0,abs((y-sol([1.d0,1.d0],x0))/sol([1.d0,1.d0],x0))
    enddo

    print*,'\n'
    print*, "Runge Kutta"
    x0=1.d0
    xf=2.d0
    steps=5
    dx=(xf-x0)/real(steps)
    
    y0=[exp(1.d0)]
    y=y0
    do i=1,steps
        y0=y
        c0=dx*func(y0,x0,N)
        c1=dx*func(y0+0.5d0*c0,x0+0.5d0*dx,N)
        c2=dx*func(y0+0.5d0*c1,x0+0.5d0*dx,N)
        c3=dx*func(y0+c2,x0+dx,N)

        y=y0+(c0+2.d0*(c1+c2)+c3)/6.d0
        x0=x0+dx
        print '(3(ES21.14,x))',x0,abs((y-sol([1.d0,1.d0],x0))/sol([1.d0,1.d0],x0))

    enddo
    contains
    
    function func(y,x,N) 
        implicit none
        real,dimension(:) :: y
        real :: x
        integer :: N
        
        real :: func(N)
        func=y*(1.d0+0.5d0/x)
    end function func
    
    function sol(y0,x) 
        real :: y0(:)
        real :: x
        real :: sol(size(y0))
        sol=sqrt(x)*exp(x)
        
    end function sol
end program ode_adams