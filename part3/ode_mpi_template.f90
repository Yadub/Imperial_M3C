!mpif90 -O3 -o adv_mpi.exe advection_mpi.f90
!to run mpiexec -n 2 adv_mpi.exe (using 2 processes in this example, could be any other number)
module advmodule
    implicit none
    real(kind=8) :: S_adv,c_adv
end module advmodule
!-------------------------------
program advection_mpi
    use mpi
    use advmodule
    implicit none
    integer :: i1,j1
    integer :: nt,n !number of time steps, number of spatial points
    real(kind=8) :: dt,dx,pi !time step, grid size, advection eqn parameters
    integer :: myid, numprocs, ierr
    real(kind=8), allocatable, dimension(:) :: x,f0,f !grid, initial condition, solution

 ! Initialize MPI
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

!gather input
    open(unit=10,file='data_mpi.in')
        read(10,*) n
        read(10,*) nt
        read(10,*) dt
        read(10,*) c_adv
        read(10,*) S_adv
    close(10)

    allocate(x(n),f0(n),f(n))

!make grid from 0 to 1-dx, dx = 1/n
    do i1=1,n
        x(i1) = dble(i1-1)/dble(n)
    end do
    dx = 1.d0/dble(n)

!generate initial condition
    pi = acos(-1.d0)
    f0 = sin(2.d0*pi*x)


!compute solution
    call euler_mpi(MPI_COMM_WORLD,numprocs,n,0.d0,f0,dt,nt,dx,f)


!output solution
        call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
       if (myid==0) then
        print *, 'max(f)=',maxval(abs(f))
        open(unit=11,file='fmpi.dat')
        do i1=1,n
            write(11,*) f(i1)
        end do
        close(11)
    end if
    !can be loaded in python with: f=np.loadtxt('fmpi.dat')

    call MPI_FINALIZE(ierr)
end program advection_mpi
!-------------------------------


subroutine euler_mpi(comm,numprocs,n,t0,y0,dt,nt,dx,y)
    !explicit Euler method
    use mpi
    use advmodule
    implicit none
    integer, intent (in) :: n,nt
    real(kind=8), dimension(n), intent(in) :: y0
    real(kind=8), intent(in) :: t0,dt,dx
    real(kind=8), dimension(n), intent(out) :: y
    real(kind=8) :: t
	integer :: i1,k,istart,iend
    integer :: comm,myid,ierr,numprocs


    call MPI_COMM_RANK(comm, myid, ierr)
    print *, 'start euler_mpi, myid=',myid

    !set initial conditions
    y = y0
    t = t0

    !generate decomposition and allocate sub-domain variables
    call mpe_decomp1d(size(y),numprocs,myid,istart,iend)

    print *, 'istart,iend,threadID,npart=',istart,iend,myid


    !time marching
    do k = 1,nt


        call RHS_mpi(ADD CODE HERE)
        ylocal= ylocal + dt*Rpart !ylocal must be declared and defined, Rpart must be declared, and 
                                  !should be returned by RHS_mpi
        t = t + dt
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    end do

    print *, 'before collection',myid, maxval(abs(ylocal))
    !collect ylocal from each processor onto myid=0

    if (myid==0) print *, 'finished',maxval(abs(y))



end subroutine euler_mpi
!-------------------------
subroutine RHS_mpi(nn,dx,t,f,rhs)
    use advmodule
    implicit none
    integer, intent(in) :: nn
    real(kind=8), intent(in) :: dx,t
    real(kind=8), dimension(nn), intent(in) :: f
    !output variable rhs must be declared

end subroutine RHS_mpi

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


