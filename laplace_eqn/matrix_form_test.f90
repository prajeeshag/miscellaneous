program form_matrix
  implicit none
  integer :: ni=3, nj=3
  real :: delx=1.0, dely=1.0
  real :: b, a
  real, allocatable :: amat(:,:), bmat(:)
  real, allocatable :: bndi(:,:), bndj(:,:), l_AP(:)
  character (len=1), allocatable :: amatc(:,:), cl_ap(:)
  integer :: im1j, ijp1, ijm1, ij, ip1j, ninj, i, j, dim_ap, istrt, nele
  character (len=512) :: cnjni, fmt
  
  ninj=ni*nj
  dim_ap = (ninj*(ninj+1)/2)
  
  b=(delx/dely)**2
  a = 2*(1-b)

  allocate(amatc(ni*nj,ni*nj))
  allocate(amat(ni*nj,ni*nj),bmat(ni*nj))
  allocate(bndi(2,nj),bndj(2,ni))
  allocate(l_AP(dim_ap))
  allocate(cl_ap(dim_ap))
  
  forall(i=1:ninj,j=1:ninj) amatc(i,j)='0'
  forall(i=1:dim_ap) cl_ap(i)='0'
  forall(i=1:ninj,j=1:ninj) amat(i,j)=0.0

  nele = ninj
  istrt = 1
  
  do i = 1, ni
     do j = 1, nj
        im1j=(i-2)*nj+j
        ijm1=(i-1)*nj+j-1
        ij=(i-1)*nj+j
        ijp1=(i-1)*nj+j+1
        ip1j=i*nj+j
!        amatc(ij,im1j) = '0'
!        amatc(ij,ijm1) = '0'
        amatc(ij,ijp1) = '0'
        amatc(ij,ip1j) = '0'
!        if ((i-1)>0) amatc(ij,im1j)='1'
!       if ((j-1)>0) amatc(ij,ijm1)='b'
        amatc(ij,ij)='a'
        if ((j+1)<=nj) amatc(ij,ijp1)='b'
        if ((i+1)<=ni) amatc(ij,ip1j)='1'
        cl_ap(istrt) = 'a'
        if (j<nj) cl_ap(istrt+1) = 'b'
        if (istrt+nj<=dim_ap) cl_ap(istrt+nj) = '1'
        istrt=istrt+nele
        nele=nele-1
     end do
  end do
  
  write(cnjni,*)ninj
  write(fmt,*)'(i2.2,1x,i2.2,'//trim(cnjni)//'(1x,a1))'
  write(*,*) trim(fmt)
  do i = 1, ni
     do j = 1, nj
        ij = (i-1)*nj+j
        write (*,trim(fmt)) i,j,amatc(ij,:)
     enddo
  enddo
  write(*,*)
  write(*,*)cl_ap

end program form_matrix