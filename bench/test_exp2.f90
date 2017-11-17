program test_exp2
use iso_c_binding, only: c_double
use iso_fortran_env, only: qp => real128
implicit none
integer, parameter :: dp = kind(0.d0)

interface
    real(c_double) function c_exp2(x) bind(c, name="exp2")
    ! The range of x is (-1075, 1024).
    import :: c_double
    real(c_double), value, intent(in) :: x
    end function
end interface

real(dp), allocatable :: x(:), r(:), r_exact(:), err_abs(:)
integer, allocatable :: err_ulp(:)
real(dp) :: t1, t2, x_min, x_max
integer :: i, j, n, n2
integer, allocatable :: seed(:)

call random_seed(size=n)
allocate(seed(n))
seed = 1
call random_seed(put=seed)

! Used for testing:
!n = 100000
! Used for benchmarking:
n = 1000
n2 = 100000
allocate(x(n), r(n), r_exact(n), err_abs(n), err_ulp(n))

call random_number(x)
x_min = -1075
x_max = 1024
x = x*(x_max - x_min) + x_min
print *, x(:10)
print *, "Running c_exp2()"
call cpu_time(t1)
do i = 1, n
    do j = 1, n2
        r(i) = c_exp2(x(i))
    end do
end do
call cpu_time(t2)
print *, "Time:", t2-t1
print *, "Computing Errors"
do i = 1, n
    r_exact(i) = real(2**real(x(i), qp), dp)
    err_abs(i) = abs(r(i)-r_exact(i))
    err_ulp(i) = ulp(r(i), r_exact(i))
end do
print *, "Done"
print *, "Max abs error:", maxval(err_abs)
print *, "Max ulp error:", maxval(err_ulp), "ulp"
call histogram(err_ulp)

contains

    integer recursive function ulp(a, b) result(r)
    real(dp), intent(in) :: a, b
    r = int(abs(a - b) / spacing(abs(b)))
    !if (r /= ulp2(a, b)) error stop
    end function

    integer recursive function ulp2(a, b) result(r)
    real(dp), intent(in) :: a, b
    real(dp) :: x
    if (a > b) then
        r = ulp2(b, a)
        return
    end if
    x = nearest(a, 1.0_dp)
    r = 0
    if (x > b) return
    r = 1
    do while (x < b)
        x = nearest(x, 1.0_dp)
        r = r + 1
    end do
    end function

    subroutine histogram(a)
    integer, intent(in) :: a(:)
    integer :: i
    print *, "Histogram"
    do i = 0, 4
        print *, i, count(a == i), 100._dp * count(a == i) / size(a), "%"
    end do
    end subroutine

end program
