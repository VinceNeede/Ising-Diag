program PC
    implicit none
    complex,parameter :: b0=cmplx(55.d0/24.d0), b1=cmplx(-59.d0/24.d0)
    complex,parameter :: b2=cmplx(37.d0/24.d0), b3=cmplx(-9.d0/24.d0)
    complex,parameter :: b_1=cmplx(9.d0/24.d0), b0_=cmplx(19.d0/24.d0)
    complex,parameter :: b1_=cmplx(-5.d0/24.d0), b2_=cmplx(1.d0/24.d0)
    complex,pointer,dimension(:) :: c0,c1,c2,c3,d
    complex,pointer,dimension(:) :: y, y0
    real :: x0,xf,dx
    integer :: steps,i,N
    
    N=1
    allocate(c0(N),c1(N),c2(N),c3(N))
    allocate(y(N),y0(N))
    steps=10
    x0=0.d0
    xf=1.d0
    dx=(xf-x0)/real(steps)
    y0=cmplx(0.d0)
    y=y0
    
    d=>c0
    call runge_kutta()
    print '(5(ES21.14,x))',x0,y0,sol(x0),abs((y0-sol(x0))/sol(x0))
    
    d=> c1
    call runge_kutta()
    print '(5(ES21.14,x))',x0,y0,sol(x0),abs((y0-sol(x0))/sol(x0))
    
    d=>c2
    call runge_kutta()
    print '(5(ES21.14,x))',x0,y0,sol(x0),abs((y0-sol(x0))/sol(x0))
    
    
    do i=1,steps-3
        c3=func(y0,x0)
        y=y0+dx*(b0*c3+b1*c2+b2*c1+b3*c0)
        c0=func(y,x0+dx)
        y=y0+dx*(b_1*c0+b0_*c3+b1_*c2+b2_*c1)    
        x0=x0+dx

        d=>y0
        y0=>y
        y=>d
        
        d=>c0
        c0=>c1
        c1=>c2
        c2=>c3
        c3=>d
        
        print '(5(ES21.14,x))',x0,y0,sol(x0),abs((y0-sol(x0))/sol(x0))
    enddo
    contains
    subroutine runge_kutta()
        d=func(y0,x0)
        c3=func(y0+0.5d0*dx*d,x0+0.5*dx)
        y = y0 + dx/6.d0*(d+2.d0*c3)
        c3=func(y0+0.5d0*dx*c3,x0+0.5*dx)
        y=y+dx/3.d0*c3
        c3=func(y0+dx*c3,x0+dx)
        y = y + dx/6.d0*c3
        d=>y0
        y0=>y
        y=>d
        x0=x0+dx
        
    end subroutine runge_kutta
    function func(y,x)
        complex, intent(in) :: y(:)
        real,intent(in) :: x
        real :: func(size(y))
        
        integer :: i
        do i =1,size(y)
            func(i)=cmplx(x)-y(i)
        enddo

    end function func
    
    function sol(x)
        real, intent(in) :: x
        real :: sol
        
        sol=exp(-x)+x-1.d0
        
    end function sol
end program PC