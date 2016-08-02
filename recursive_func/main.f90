program main
  use recursive_mod, only: factorial, legendre
  implicit none
  
  integer :: n=0, i
  real, dimension(-10:10) :: x, y

  forall(i=-10:10) x(i)=real(i)/10.0

  print *, x

  print *, 'enter a positive integer :'
  read(*,*) n

  do i=-10,10
    y(i) = legendre(n,x(i))
  enddo

  print *, y
end program main
