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
    real(kind=8) :: t,upbound,lowbound
	integer :: i1,k,istart,iend
    integer :: comm,myid,ierr,numprocs,N_local,sender,receiver
    real(kind=8), dimension(:), allocatable :: ylocal
    real(kind=8), dimension(:), allocatable :: Rpart
    integer, dimension(MPI_STATUS_SIZE) :: status
    logical, parameter :: display = .false.

    !for gather:
    integer, allocatable, dimension(:) :: Nper_proc, disps

    call MPI_COMM_RANK(comm, myid, ierr)
    print *, 'start euler_mpi, myid=',myid

    !set initial conditions
    y = y0
    t = t0

    !generate decomposition and allocate sub-domain variables
    call mpe_decomp1d(size(y),numprocs,myid,istart,iend)

    print *, 'istart,iend,threadID,npart=',istart,iend,myid

    N_local = iend - istart + 1 !set size of yloca for each proc
    allocate (ylocal(N_local),Rpart(N_local)) !allocate ylocal and Rpart based on N_local
    ylocal = y(istart:iend) !set initial condition values for ylocal

    !time marching
    do k = 1,nt
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !-----------------------------------------------------------
        !Send data at top boundary up to next processor
        !-----------------------------------------------------------
        if (myid<numprocs-1) then
            receiver = myid+1
        else
            receiver = 0
        end if

        if (myid>0) then
            sender = myid-1
        else
            sender = numprocs-1
        end if

        call MPI_SEND(ylocal(N_local),1,MPI_DOUBLE_PRECISION,receiver,0,comm,ierr)
        call MPI_RECV(lowbound,1,MPI_DOUBLE_PRECISION,sender,MPI_ANY_TAG,comm,status,ierr)

        call MPI_BARRIER(comm,ierr)

        !-----------------------------------------------------------
        !Send data at bottom boundary down to previous processor
        !-----------------------------------------------------------
        if (myid>0) then
            receiver = myid-1
        else
            receiver = numprocs-1
        end if

        if (myid<numprocs-1) then
            sender = myid+1
        else
            sender = 0
        end if
        call MPI_SEND(ylocal(1),1,MPI_DOUBLE_PRECISION,receiver,0,comm,ierr)
        call MPI_RECV(upbound,1,MPI_DOUBLE_PRECISION,sender,MPI_ANY_TAG,comm,status,ierr)
        !----finished sending/receiving

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)


        call RHS_mpi(N_local+2,dx,t,(/lowbound,ylocal,upbound/),Rpart)
        ylocal= ylocal + dt*Rpart !ylocal must be declared and defined, Rpart must be declared, and
                                  !should be returned by RHS_mpi
        t = t + dt
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    end do

    print *, 'before collection',myid, maxval(abs(ylocal))
    !collect ylocal from each processor onto myid=0

    allocate(Nper_proc(numprocs),disps(numprocs))

    !gather N_local from each proc to array Nper_proc on myid=0
    call MPI_GATHER(N_local,1,MPI_INT,Nper_proc,1,MPI_INT,0,comm,ierr)
    if (myid==0) then
        if (display) print *,'Nper_proc=',Nper_proc
        disps(1)=0
        do i1=2,numprocs
            disps(i1) = disps(i1-1)+Nper_proc(i1-1) !needed for gatherv below
        end do
        if (display) print *, 'disps=', disps
    end if

    call MPI_GATHERV(ylocal,N_local,MPI_DOUBLE_PRECISION,y,Nper_proc,disps,MPI_DOUBLE_PRECISION,0,comm,ierr)

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
    real(kind=8), dimension(size(f)-2), intent(out) :: rhs

    rhs = S_adv-c_adv*0.5d0*(f(3:nn)-f(1:nn-2))/dx

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


