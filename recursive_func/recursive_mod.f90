module recursive_mod
  implicit none
  
  contains

  recursive function factorial(n) result(res)
    integer, intent(in) :: n
    integer :: res
    if (n == 0) then
        res = 1
    elseif(n > 0) then
        res = n * factorial(n - 1)
    else
        stop 'Error: Factorial not defined for negative integers'
    end if
  end function factorial

  recursive function legendre(n,x) result(y)
    integer, intent(in) :: n
    real, intent(in) :: x
    real :: y

    if ( n == 0 ) then
      y = 1.0
    else if ( n == 1 ) then
      y = x
    else
      y = ( (2*n-1) * x * legendre(n-1,x) - &
            (n-1) * legendre(n-2,x) ) / n
    endif
  end function legendre

end module recursive_mod
