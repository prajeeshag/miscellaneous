module recursive_mod
	implicit none
	
	contains

	recursive function factorial(n) result(res)
	  integer :: res, n
	  if (n == 0) then
	      res = 1
	  elseif(n > 0) then
	      res = n * factorial(n - 1)
		else
				stop 'Error: Factorial not defined for negative integers'
	  end if
	end

end module recursive_mod
