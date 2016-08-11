program form_matrix
  use netcdf_write, only: write_nc
  
  implicit none
  integer :: ni=10, nj=10, ntstep=10000, outstep
  real :: delx=1.0, dely=1.0, delt=1.0
  real :: b, a, p, q, r
  real, allocatable :: amat(:,:)
  real, allocatable :: bndx(:,:), bndy(:,:), l_AP(:), l_B(:,:), tmp_AP(:)
  integer, allocatable :: ipiv(:)
  integer, allocatable :: xaxis(:), yaxis(:)
  integer ::  ij, ninj, i, j, dim_ap, istrt, nele, ldb, nrhs=1, info, n
  character :: uplo='L'
  character (len=512) :: fmt, cnj
  logical :: bndopt(4)=(/.true.,.true.,.true.,.true./), exist

  namelist/heat_eqn_nml/bndopt, ni, nj, delt, ntstep, outstep

  outstep = ntstep

  inquire(file='input.nml',exist=exist)

  if(exist) then
     open(10, file='input.nml')
     read(10, nml=heat_eqn_nml)
     close(10)
  endif

  write(*,nml=heat_eqn_nml)
  
  ninj=ni*nj

  a = delt/(delx**2)
  b = delt/(dely**2)

  p = -1.*a; q = -1.*b; r = 1+2*a+2*b 
  
  dim_ap = (ninj*(ninj+1)/2)
  
  allocate(bndx(2,nj),bndy(2,ni))
  
  allocate(l_AP(dim_ap))
  allocate(l_B(ninj,nrhs))
  allocate(tmp_AP(dim_ap))
  
  allocate(amat(ni,nj))
  allocate(yaxis(nj))
  allocate(xaxis(ni))
  allocate(ipiv(ninj))
  
  forall(i=1:ni) xaxis(i)=i
  forall(i=1:nj) yaxis(i)=i
  
  bndx(:,:) = 0.0; bndy(:,:) = 0.0
  
  if(bndopt(1)) forall(i=1:nj) bndx(1,i) = sin(3.1414 * (i-1)/(nj-1))
  if(bndopt(3)) forall(i=1:nj) bndx(2,i) = sin(3.1414 * (i-1)/(nj-1))
  if(bndopt(2)) forall(i=1:ni) bndy(1,i) = sin(3.1414 * (i-1)/(ni-1))
  if(bndopt(4)) forall(i=1:ni) bndy(2,i) = sin(3.1414 * (i-1)/(ni-1))
  
  l_AP(:) = 0.0
  
  l_B(:,:) = 0.0 ! put initial condition here

  nele = ninj
  istrt = 1
  
  do i = 1, ni
     do j = 1, nj
        l_AP(istrt) = r
        if (j<nj) l_AP(istrt+1) = q
        if (i<ni) l_AP(istrt+nj) = p
        istrt=istrt+nele
        nele=nele-1
     enddo
  enddo

  do n = 1, ntstep
     tmp_AP = l_AP
     do i = 1, ni
        do j =1, nj
           ij = (i-1)*nj+j
           if (i==1) l_B(ij,1) = l_B(ij,1) - bndx(1,j) * p
           if (i==ni) l_B(ij,1) = l_B(ij,1) - bndx(2,j) * p
           if (j==1) l_B(ij,1) = l_B(ij,1) - bndy(1,i) * q
           if (j==nj) l_B(ij,1) = l_B(ij,1) - bndy(2,i) * q
        end do
     end do
    
     call dspsv(uplo, ninj, nrhs, tmp_AP, ipiv, l_B, ninj, info)
     ! send_data(l_B)
     
     if (info/=0) then
        print *, 'error while running dppsv, error code:', info
        stop "ERROR"
     endif
     
     amat(:,:) = reshape(l_B(:,1),(/ni,nj/))     
     if (mod(n,outstep)==0) call write_nc(real(xaxis), real(yaxis), amat, 'temp', n/outstep)
  enddo !time loop ends





end program form_matrix
