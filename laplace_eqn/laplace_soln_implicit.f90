program form_matrix
  use netcdf_write, only: write_nc
  
  implicit none
  integer :: ni=10, nj=10
  real :: delx=1.0, dely=1.0
  real :: b, a
  real, allocatable :: amat(:,:)
  real, allocatable :: bndx(:,:), bndy(:,:), l_AP(:), l_B(:,:)
  integer, allocatable :: ipiv(:)
  integer, allocatable :: xaxis(:), yaxis(:)
  integer ::  ij, ninj, i, j, dim_ap, istrt, nele, ldb, nrhs=1, info
  character :: uplo='L'
  character (len=512) :: fmt, cnj
  logical :: bndopt(4)=(/.true.,.true.,.true.,.true./), exist

  namelist/laplace_nml/bndopt, ni, nj

  inquire(file='input.nml',exist=exist)

  if(exist) then
     open(10, file='input.nml')
     read(10, nml=laplace_nml)
     close(10)
  endif

  write(*,nml=laplace_nml)
  
  ninj=ni*nj
  dim_ap = (ninj*(ninj+1)/2)
  
  b=(delx/dely)**2
  a = -2.0*(1+b)

  allocate(bndx(2,nj),bndy(2,ni))
  allocate(l_AP(dim_ap))
  allocate(l_B(ninj,nrhs))
  allocate(amat(ni,nj))
  allocate(yaxis(nj))
  allocate(xaxis(ni))

  forall(i=1:ni) xaxis(i)=i
  forall(i=1:nj) yaxis(i)=i
  
  bndx(:,:) = 0.0; bndy(:,:) = 0.0
  
  if(bndopt(1)) forall(i=1:nj) bndx(1,i) = sin(3.1414 * (i-1)/(nj-1))
  if(bndopt(3)) forall(i=1:nj) bndx(2,i) = sin(3.1414 * (i-1)/(nj-1))
  if(bndopt(2)) forall(i=1:ni) bndy(1,i) = sin(3.1414 * (i-1)/(ni-1))
  if(bndopt(4)) forall(i=1:ni) bndy(2,i) = sin(3.1414 * (i-1)/(ni-1))
  
  l_AP(:) = 0.0
  l_B(:,:) = 0.0
  
  nele = ninj
  istrt = 1
  
  do i = 1, ni
     do j = 1, nj
        ij = (i-1)*nj+j
        l_AP(istrt) = a
        if (j<nj) l_AP(istrt+1) = b
        if (istrt+nj<=dim_ap) l_AP(istrt+nj) = 1.
        istrt=istrt+nele
        nele=nele-1
        if (i==1) l_B(ij,1) = l_B(ij,1) - bndx(1,j)
        if (i==ni) l_B(ij,1) = l_B(ij,1) - bndx(2,j)
        if (j==1) l_B(ij,1) = l_B(ij,1) - bndy(1,i) * b
        if (j==nj) l_B(ij,1) = l_B(ij,1) - bndy(2,i) * b
     end do
  end do
  
!  call dppsv(uplo, ninj, nrhs, l_AP, l_B, ninj, info)
  
  allocate(ipiv(ninj))
  
  call dspsv(uplo, ninj, nrhs, l_AP, ipiv, l_B, ninj, info)

  if (info/=0) then
     print *, 'error while running dppsv, error code:', info
     stop
  endif

  amat(:,:) = reshape(l_B(:,1),(/ni,nj/))

  call write_nc(real(xaxis), real(yaxis), amat, 'temp')

end program form_matrix
