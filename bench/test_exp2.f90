program test_exp2
use iso_c_binding, only: c_double
implicit none
integer, parameter :: dp = kind(0.d0)

interface
    real(c_double) function c_exp2(x) bind(c, name="exp2")
    import :: c_double
    real(c_double), value, intent(in) :: x
    end function
end interface

real(dp) :: y

y = c_exp2(0.7_dp)
print *, y

end program
