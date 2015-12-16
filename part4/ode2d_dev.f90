!Yadu Bhageria
!CID: 00733164
!Project part 4
!module to use RK4 time marching to solve an initial value problem
!Solves: dy/dt = RHS(y,t) where y is a 2d matrix
module ode2d
    use omp_lib
    integer :: numthreads,par
    save
contains
!-----------------
subroutine rk4(t0,y0,dt,nt,y,time)
    !4th order RK method
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: y0
    real(kind=8), intent(in) :: t0,dt
	integer, intent (in) :: nt
    real(kind=8), dimension(size(y0,1),size(y0,2)), intent(out) :: y
    real(kind=8), dimension(size(y0,1),size(y0,2)) :: f1, f2, f3, f4
    real(kind=8) :: t,halfdt,fac
    real(kind=8), intent(out) :: time
	integer:: k
    integer(kind=8) :: t1,t2,clock_rate
        print *, 'par=',par
        if (par==1) call omp_set_num_threads(numthreads)
        halfdt = 0.5d0*dt
        fac = 1.d0/6.d0

        y = y0
        t = t0
        call system_clock(t1)
        do k = 1, nt

           f1 = dt*RHS(t, y)

           f2 = dt*RHS(t + halfdt, y + 0.5d0*f1)

           f3 = dt*RHS(t + halfdt, y + 0.5d0*f2)

           f4 = dt*RHS(t + dt, y + f3)

           y = y + (f1 + 2*f2  + 2*f3 + f4)*fac

           t = t + dt*dble(k)

        end do
        call system_clock(t2,clock_rate)
        time = dble(t2-t1)/dble(clock_rate)
end subroutine rk4
!-----------------
subroutine euler(t0,y0,dt,nt,y)
    !explicit Euler method
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: y0
    real(kind=8), intent(in) :: t0,dt
	integer, intent (in) :: nt
    real(kind=8), dimension(size(y0,1),size(y0,2)), intent(out) :: y
    real(kind=8) :: t,halfdt,fac
	integer:: k

    y = y0
    t = t0
    do k = 1,nt

        y = y + dt*RHS(t,y)
        t = t + dt

    end do


end subroutine euler
!-----------------
function RHS(t,f)
    use fdmodule2d
    use fdmodule
    use advmodule
    implicit none
    integer :: i1
    real(kind=8), intent(in) :: t
    real(kind=8), dimension(:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1),size(f,2)) :: RHS,df1,df2

    if (par==1) then
        !compute df/dx1
        n = n1
        dx = dx1
        !$OMP parallel do
        do i1=1,n2
            call fd2(f(i1,:),df1(i1,:))
        end do
        !OMP end parallel do

        n = n2
        dx = dx2
        !compute df/dx2, df_amp,df_max
        !$OMP parallel do
        do i1=1,n1
            call fd2(f(:,i1),df2(:,i1))
        end do
        !$OMP end parallel do
    else
        call grad(f,df1,df2)
    end if
    
    RHS = S_adv - c1_adv*df1 - c2_adv*df2

end function RHS


end module ode2d
!-----------------



