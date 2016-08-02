program main
  use recursive_mod, only: factorial, legendre
  implicit none
  
  integer :: n=0, i
  real, dimension(-100:100) :: x, y
	character (len=32) :: carg

  forall(i=-100:100) x(i)=real(i)/100.0

	call getarg(1,carg)

  read(carg,*) n

  open(10,file='output.dat')
  do i=-100,100
    write(10,*) x(i), legendre(n,x(i))
  enddo
	close(10)
end program main
