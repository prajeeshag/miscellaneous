program main
	use recursive_mod, only: factorial
	implicit none
	
	integer :: n=0

  print *, 'enter a positive integer :'
	read(*,*) n

	write(*,*) factorial(n)
end program main
