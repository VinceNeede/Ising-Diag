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
        
        
        complex :: c1(N), c2(N), c3(N), c4(N)
        
        call fun(c1, y0, x0, N)
        call fun(c2, y0+dx*(0.5,0.)*c1, x0+0.5*dx, N)
        call fun(c3, y0+dx*(0.5,0.)*c2, x0+0.5*dx, N)
        call fun(c4, y0+dx*c3, x0+dx, N)
        
        !y=y0+(c1+(2.d0,0.)*(c2+c3)+c4)*dx/(6.d0,0.)
        y=y0+c1*cmplx(dx)
    End Subroutine Runge_Kutta
end module ODE
program main
    use ode
    implicit none
    integer,parameter :: N=2
    real :: x0,dx,x_fin
    integer :: steps, ii
    complex ::  y0(N), y(N),y0_0(N)


    x0=1.d0
    x_fin=2.
    steps=10000
    dx=(x_fin-x0)/real(steps)

    y0_0=[1.,1.]
    y0=y0_0
    do ii =1,steps
        call Runge_Kutta(fun,y0,x0,dx,y,N)
        x0=x0+dx
        print'(13(ES21.14))', x0,abs(y-x0**2)/(x0**2)
        y0=y
    enddo
    contains
    subroutine fun(z,y,x,N)
        implicit none
        integer, intent(in) :: N
        Real, Intent(in) :: x
        complex, intent(in) :: y(N)
        complex, intent(out) :: z(N)

        z=2.*x
    end subroutine
end program main
