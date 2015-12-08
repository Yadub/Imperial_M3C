!Project part 2
!Module containing routines for differentiating an array of size N with
!2nd order and 4th order-compact finite differences
!Test routine test_fd applies these methods to a Gaussian function and
!returns the error while test_fd_time also returns timing information.

module fdmodule
    implicit none
    integer :: n
    real(kind=8) :: dx
    save

contains
!-------------------
subroutine fd2(f,df)
    !2nd order centered finite difference scheme with
    !periodic boundary conditions
    !Assumes N, dx have been set in calling program
    implicit none
    real(kind=8) :: invtwodx
    real(kind=8), dimension(:), intent(in) :: f
    real(kind=8), dimension(size(f)), intent(out) :: df

    invtwodx = 0.5d0/dx
    df(2:N-1) = invtwodx*(f(3:N)-f(1:N-2))

    !b.c.'s
    df(1) = invtwodx*(f(2) - f(N))
    df(N) = invtwodx*(f(1) - f(N-1))

end subroutine fd2
!-----------------

subroutine cfd4(f,df)
    !4th order centered finite difference scheme with
    !one-sided 3rd order boundary conditions
    !Assumes N, dx have been set in calling program
    implicit none
    real(kind=8) :: dxfac,c1,c2,c3,c4
    real(kind=8), dimension(:), intent(in) :: f
    real(kind=8), dimension(size(f)), intent(out) :: df
!variables for DGTSV
    real(kind=8), dimension(size(f)) :: D
    real(kind=8), dimension(size(f)-1) :: DL,DU
    integer :: INFO

    !needed constants
    dxfac = 1.d0/dx
    c1 = 3.d0
    c2 = 2.5d0
    c3 = 2.d0
    c4 = 0.5d0

    !RHS
    df(2:N-1) = c1*(f(3:N)-f(1:N-2))
    df(1) = (-c2*f(1) + c3*f(2) + c4*f(3))
    df(N) = (c2*f(N) - c3*f(N-1) - c4*f(N-2))

    df = dxfac*df


    !Diagonals for DGTSV
    D =  4.d0
    DL = 1.d0
    DU = 1.d0

    !b.c.'s
    D(1) = 1.d0
    DU(1) = 2.d0
    D(N) = 1.d0
    DL(N-1) = 2.d0

    !solution is returned in df
    call DGTSV(N,1,DL,D,DU,df,N,INFO)


    if (info .ne. 0) then
        print *, 'error in compactfd4, DGTSV gives INFO=',INFO
        STOP
    end if

end subroutine cfd4
!------------------

subroutine test_fd(alpha,error)
    !test finite difference schemes with Gaussian function
    !return L1 error for each scheme
    !Assumes N, dx have been set in calling program
    implicit none
    integer :: i1
    real(kind=8), intent(in) :: alpha
    real(kind=8), intent(out) :: error(2)
    real(kind=8), dimension(N) :: x,fgauss,dfgauss,dfgauss_exact
    real(kind=8) :: x0

    !generate grid
    do i1=0,N-1
        x(i1+1) = dble(i1)*dx
    end do

    !test function and exact solution
    x0 = 0.5d0*x(N)
    fgauss = exp(-alpha*(x-x0)**2)
    dfgauss_exact = (-2.d0*alpha)*(x-x0)*fgauss

    !tests
    call fd2(fgauss,dfgauss)
    error(1) = sum(abs(dfgauss-dfgauss_exact))/dble(N)

    call cfd4(fgauss,dfgauss)
    error(2) = sum(abs(dfgauss-dfgauss_exact))/dble(N)


end subroutine test_fd
!---------------------

subroutine test_fd_time(alpha,error,time)
    !test finite difference schemes with Gaussian function
    !return L1 error for each scheme and wall-clock time
    !for 1000 calls to fd2 and cfd4
    !Assumes N, dx have been set in calling program
    !Assumes N, dx have been set in calling program
    implicit none
    integer :: i1
    real(kind=8), intent(in) :: alpha
    real(kind=8), intent(out) :: error(2),time(2)
    real(kind=8), dimension(N) :: x,fgauss,dfgauss,dfgauss_exact
    real(kind=8) :: x0
    integer(kind=8) :: t1,t2,clock_rate

    !generate grid
    do i1=0,N-1
        x(i1+1) = dble(i1)*dx
    end do

    !test function and exact solution
    x0 = 0.5d0*x(N)
    fgauss = exp(-alpha*(x-x0)**2)
    dfgauss_exact = (-2.d0*alpha)*(x-x0)*fgauss


!Test fd2
    call system_clock(t1)
    do i1=1,1000
        call fd2(fgauss,dfgauss)
    end do
    call system_clock(t2,clock_rate)

    time(1) = dble(t2-t1)/dble(clock_rate)
    error(1) = sum(abs(dfgauss-dfgauss_exact))/dble(N)


!Test cfd4
    call system_clock(t1)

    do i1=1,1000
        call cfd4(fgauss,dfgauss)
    end do

    call system_clock(t2,clock_rate)
    time(2) = dble(t2-t1)/dble(clock_rate)
    error(2) = sum(abs(dfgauss-dfgauss_exact))/dble(N)


end subroutine test_fd_time
!---------------------------
end module fdmodule















