!Yadu Bhageria
!CID:00733164
!Copied from solution for Homework 4, part 1
!module containing routines for differentiating n2 x n1 matrix with
!2nd order and 4th order compact finite differences
!Test routines apply these methods to a simple test function and
!return the error and time.

module fdmodule2d
    !$    use omp_lib
    use fdmodule
    implicit none
    integer :: n1,n2,numthreads
    real(kind=8) :: dx1,dx2
    save
contains
!-------------------
subroutine grad(f,df1,df2)
    !compute df1=df/dx1, df2=df/dx2, df_amp = |df1| + |df2|, df_max = max(df_amp)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1),size(f,2)), intent(out) :: df1,df2
    integer :: i1

    !compute df/dx1
    n = n1
    dx = dx1
    do i1=1,n2
        call fd2(f(i1,:),df1(i1,:))
    end do


    !compute df/dx2
    n = n2
    dx = dx2
    do i1=1,n1
        call fd2(f(:,i1),df2(:,i1))
    end do

end subroutine grad
!---------------------------
subroutine test_grad(error,time)
    !tests accuracy and speed of grad, assumes n1,n2,dx1,dx2 have been set in calling program
    implicit none
    integer :: i1,j1
    real(kind=8), allocatable, dimension(:) :: x1,x2 !coordinates stored in arrays
    real(kind=8), allocatable, dimension(:,:) :: x11,x22 !coordinates stored in matrices
    real(kind=8), allocatable, dimension(:,:) :: ftest,df1,df2,df_amp_exact !test function and results frorom grad
    integer(kind=8) :: t1,t2,clock_rate
    real(kind=8), intent(out) :: time,error(2) !error is two-element array

    allocate(x1(n1),x2(n2),x11(n2,n1),x22(n2,n1),ftest(n2,n1),df1(n2,n1),df2(n2,n1),df_amp_exact(n2,n1))


    !generate mesh
    dx1 = 1.d0/dble(n1-1)
    dx2 = 1.d0/dble(n2-1)

    do i1=1,n1
        x1(i1) = dble(i1-1)*dx1
    end do

    do i1=1,n2
        x2(i1) = dble(i1-1)*dx2
    end do

    do i1=1,n2
        x11(i1,:) = x1
    end do

    do i1=1,n1
        x22(:,i1) = x2
    end do

    ftest = sin(x11)*cos(x22) !test function

    call system_clock(t1)
    call grad(ftest,df1,df2)
    call system_clock(t2,clock_rate)

    time = dble(t2-t1)/dble(clock_rate)

    error(1) = sum(abs(df1-cos(x11)*cos(x22)))/(n1*n2)
    error(2) = sum(abs(df2+sin(x11)*sin(x22)))/(n1*n2)
end subroutine test_grad
!-------------------
subroutine grad_omp(f,df1,df2,df_amp,df_max)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1),size(f,2)), intent(out) :: df1,df2,df_amp
    real(kind=8), intent(out) :: df_max
    integer :: i1

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
    !$OMP parallel do reduction(max:df_max)
    do i1=1,n1
    call fd2(f(:,i1),df2(:,i1))
    df_amp(:,i1) = abs(df1(:,i1))+ abs(df2(:,i1))
    df_max = max(df_max,maxval(df_amp(:,i1)))
    end do
    !OMP end parallel do


end subroutine grad_omp
!---------------------------
subroutine test_grad_omp(error,time)
    !tests accuracy and speed of grad_omp, assumes n1,n2,dx1,dx2, and numthreads have been set in calling program
    implicit none
    integer :: i1,j1
    real(kind=8), allocatable, dimension(:) :: x1,x2 !coordinates stored in arrays
    real(kind=8), allocatable, dimension(:,:) :: x11,x22 !coordinates stored in matrices
    real(kind=8), allocatable, dimension(:,:) :: ftest,df1,df2,df_amp,df_amp_exact !test function and results frorom grad
    real(kind=8) :: df_max
    integer(kind=8) :: t1,t2,clock_rate
    real(kind=8), intent(out) :: time,error(2) !error is two-element array

    allocate(x1(n1),x2(n2),x11(n2,n1),x22(n2,n1),ftest(n2,n1),df1(n2,n1),df2(n2,n1), &
    df_amp(n2,n1),df_amp_exact(n2,n1))

    !$ call omp_set_num_threads(numthreads)

    !generate mesh
    dx1 = 1.d0/dble(n1-1)
    dx2 = 1.d0/dble(n2-1)

    do i1=1,n1
    x1(i1) = dble(i1-1)*dx1
    end do

    do i1=1,n2
    x2(i1) = dble(i1-1)*dx2
    end do

    do i1=1,n2
    x11(i1,:) = x1
    end do

    do i1=1,n1
    x22(:,i1) = x2
    end do

    ftest = sin(x11)*cos(x22) !test function

    call system_clock(t1)
    call grad_omp(ftest,df1,df2,df_amp,df_max)
    call system_clock(t2,clock_rate)

    time = dble(t2-t1)/dble(clock_rate)

    error(1) = sum(abs(df1-cos(x11)*cos(x22)))/(n1*n2)
    error(2) = sum(abs(df2+sin(x11)*sin(x22)))/(n1*n2)

end subroutine test_grad_omp
!---------------------------
end module fdmodule2d


    !test program, not required for assignment
    program test
    !$    use omp_lib
    use fdmodule2d
    implicit none
    integer :: i1,j1
    real(kind=8), allocatable, dimension(:) :: x1,x2 !coordinates stored in arrays
    real(kind=8), allocatable, dimension(:,:) :: x11,x22 !coordinates stored in matrices
    real(kind=8), allocatable, dimension(:,:) :: ftest,df1,df2,df_amp_exact !test function and results from grad
    real(kind=8) :: time,error(2)
    integer(kind=8) :: t1,t2,clock_rate

    open(unit=10,file='data.in')
    read(10,*) n1
    read(10,*) n2
    read(10,*) numThreads
    close(10)

    !$ call omp_set_num_threads(numthreads)


    dx1 = 1.d0/dble(n1-1)
    dx2 = 1.d0/dble(n2-1)

    call test_grad(error,time)

    print *, 'test, n1,n2,numThreads=',n1,n2,numThreads
    print *, 'error=',error
    print *, 'time=',time



end program test