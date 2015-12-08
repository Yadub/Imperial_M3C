!Project part 3
!module contains parameters for advection equation and main program tests Euler integration of adv. eqn.
!You may need to compile this code with gfortran -c advmodule.f90 before using f2py
module advmodule
    implicit none
    real(kind=8) :: S_adv,c_adv
    real(kind=8) :: c1_adv,c2_adv !for 2d advection eqn.

end module advmodule
!--------------------

program test_adv
!test program for solving advection eqn with Euler time marching (not required for assignment)
	use ode
    use fdmodule
    use advmodule
	implicit none
	integer :: i1,t1,nt
	real(kind=8) :: dt,pi
	real(kind=8),allocatable, dimension(:) :: x,f0,f

    pi = acos(-1.d0)

!Read input
    open(unit=10,file='data.in')
    read(10,*) nt
    read(10,*) n
    close(10)


!set numerical parameters
    allocate(x(n),f0(n),f(n))
    dx = 1.d0/(n)
    dt = 1.d0/nt

!generate grid and initial condition
    do i1=1,n
        x(i1) = dx*(i1-1)
    end do
    f0 = sin((2*pi)*x)

!set eqn. parameters
    c_adv = 1.d0
    S_adv = 0.d0

!euler integration
    call euler(0.d0,f0,dt,nt,f)

!add code to compare f to exact solution


end program test_adv
