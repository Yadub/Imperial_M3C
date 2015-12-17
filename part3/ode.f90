!Yadu Bhageria
!CID: 00733164
!
!Project part 3
!module to use RK4 or Euler time marching to solve an initial value problem
!Solves: dy/dt = RHS(y,t)
module ode

contains
subroutine rk4(t0,y0,dt,nt,y)
    !4th order RK method
    implicit none
    real(kind=8), dimension(:), intent(in) :: y0
    real(kind=8), intent(in) :: t0,dt
	integer, intent (in) :: nt
    real(kind=8), dimension(size(y0)), intent(out) :: y
    real(kind=8), dimension(size(y0)) :: f1, f2, f3, f4
    real(kind=8) :: t,halfdt,fac
	integer:: k

        halfdt = 0.5d0*dt
        fac = 1.d0/6.d0

        y = y0
        t = t0

        do k = 1, nt

           f1 = dt*RHS(t, y)

           f2 = dt*RHS(t + halfdt, y + 0.5d0*f1)

           f3 = dt*RHS(t + halfdt, y + 0.5d0*f2)

           f4 = dt*RHS(t + dt, y + f3)

           y = y + (f1 + 2*f2  + 2*f3 + f4)*fac

           t = t + dt*dble(k)

        end do
end subroutine rk4
!------------------
subroutine euler(t0,y0,dt,nt,y)
    !explicit Euler method
    implicit none
    real(kind=8), dimension(:), intent(in) :: y0
    real(kind=8), intent(in) :: t0,dt
	integer, intent (in) :: nt
    real(kind=8), dimension(size(y0)), intent(out) :: y
    real(kind=8) :: t,halfdt,fac
	integer:: k

    y = y0
    t = t0
    do k = 1,nt

        y = y + dt*RHS(t,y)
        t = t + dt

    end do


end subroutine euler
!--------------------
subroutine euler_omp(t0,y0,dt,nt,numthreads,y)
    !explicit Euler method, parallelized with OpenMP
    use omp_lib
    implicit none
    real(kind=8), dimension(:), intent(in) :: y0
    real(kind=8), intent(in) :: t0,dt
	integer, intent (in) :: nt,numthreads
    real(kind=8), dimension(size(y0)), intent(out) :: y
    real(kind=8) :: t,upbound,lowbound
	integer:: k,threadID,istart,iend,ntot

!$ call omp_set_num_threads(numthreads)

    y = y0
    t = t0
    ntot = size(y)

!$OMP parallel private(istart,iend,threadID,upbound,lowbound)
    threadID = omp_get_thread_num()
    call mpe_decomp1d(ntot,numthreads,threadID,istart,iend) !construct domain decomposition
    print *, 'istart,iend,threadID=',istart,iend,threadID

    do k = 1,nt
        !$OMP barrier
        upbound = y(modulo(iend+1,ntot))
        lowbound = y(modulo(istart-2,ntot)+1)
        !$OMP barrier
        y(istart:iend) = y(istart:iend) + dt*RHS_omp(t,(/ lowbound, y(istart:iend), upbound /))
        !$OMP master
        t = t + dt
        !$OMP end master

    end do
    print *, 'finished loop:',threadID,maxval(abs(y(istart:iend))) !this will slow the code but is included for assessment.
!$OMP end parallel

!The code has been parallelized by having the threads split the work on the inputted array, y0, using mpe_decomp1d
!For each time step, each thread first stores the upper and lower boundary values needed to compute the derivative on private variables upbound and lowbound before they all then update their respective sections of the array, y. This loops over nt times

end subroutine euler_omp

!-----------------------
function RHS_omp(t,f)
    !called by euler_omp
    !RHS_omp = df/dt
    use fdmodule
    use advmodule
    implicit none
    real(kind=8), intent(in) :: t
    real(kind=8), dimension(:), intent(in) :: f
    !add variable declaration for RHS_omp
    real(kind=8), dimension(size(f)-2) :: RHS_omp

    RHS_omp = S_adv-c_adv*0.5d0*(f(3:size(f))-f(1:size(f)-2))/dx

 end function RHS_omp
!--------------------------------------

function RHS(t,f)
    !called by euler and rk4
    !RHS = df/dt
    use fdmodule
    use advmodule
    implicit none
    real(kind=8), intent(in) :: t
    real(kind=8), dimension(:), intent(in) :: f
    real(kind=8), dimension(size(f)) :: RHS

    call fd2(f,RHS)

    RHS = S_adv-c_adv*RHS

end function RHS
!-----------------
end module ode


!--------------------------------------------------------------------
!  (C) 2001 by Argonne National Laboratory.
!      See COPYRIGHT in online MPE documentation.
!  This file contains a routine for producing a decomposition of a 1-d array
!  when given a number of processors.  It may be used in "direct" product
!  decomposition.  The values returned assume a "global" domain in [1:n]
!
subroutine MPE_DECOMP1D( n, numprocs, myid, s, e )
    implicit none
    integer :: n, numprocs, myid, s, e
    integer :: nlocal
    integer :: deficit

    nlocal  = n / numprocs
    s       = myid * nlocal + 1
    deficit = mod(n,numprocs)
    s       = s + min(myid,deficit)
    if (myid .lt. deficit) then
        nlocal = nlocal + 1
    endif
    e = s + nlocal - 1
    if (e .gt. n .or. myid .eq. numprocs-1) e = n

end subroutine MPE_DECOMP1D

