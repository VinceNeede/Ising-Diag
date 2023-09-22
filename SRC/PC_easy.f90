program PC
    implicit none
    integer,parameter :: steps=10
    real :: y(steps),c0,c1,c2,c3,y0,x0,dx,xf
    real :: b0,b1,b2,b3
    integer :: i,j

    x0=0.d0
    xf=1.d0
    dx=(xf-x0)/real(steps)
    y0=0.d0
    y(1)=y0

    b0=55.d0/24.d0
    b1=-59.d0/24.d0
    b2=37.d0/24.d0
    b3=- 9.d0/24.d0
    do i=1,steps
        if (i.le.3) then
            c0=dx*func(y(i),x0)
            c1=dx*func(y(i)+0.5d0*c0,x0+0.5d0*dx)
            c2=dx*func(y(i)+0.5d0*c1,x0+0.5d0*dx)
            c3=dx*func(y(i)+c2,x0+dx)
            y(i+1)=y(i)+(c0+2.d0*(c1+c2)+c3)/6.d0
        else
            y(i+1)=y(i)+dx*(b0*func(y(i),x0)&
               +b1*func(y(i-1),x0-dx)+b2*func(y(i-2),x0-2.d0*dx) +b3*func(y(i-3),x0-3.d0*dx))
            ! do j=1,2
            y(i+1)=y(i)+dx/24.d0*(9.d0 * func(y(i+1),x0 + dx) + 19.d0 * func(y(i),x0) &
                - 5.d0 * func(y(i - 1),x0 - dx) + func(y(i - 2),x0 - 2.d0 * dx))
            ! enddo
        endif
        x0=x0+dx
        print '(4(ES21.14,x))',x0,y(i+1),sol(x0),abs((y(i+1)-sol(x0))/sol(x0))
    enddo

    contains
    function func(y,x)
        real, intent(in) :: y,x
        real :: func
        

            func=x-y
    end function func

    function sol(x)
        real, intent(in) :: x
        real :: sol
    
        sol=exp(-x)+x-1.d0
        
    end function sol
end program PC