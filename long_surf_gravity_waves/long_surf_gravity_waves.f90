program main
  use netcdf_write, only : write_nc
  
  implicit none
  real, parameter :: GRAV = 9.8
  real :: delt = 0.1, delx = 10.0, dely = 10.0
  integer :: ntstep=2000
  integer :: ni = 300, nj = 100, nip1, njp1
  real, allocatable :: uvel(:,:,:), vvel(:,:,:), eta(:,:,:), h(:,:), h0(:,:), topo(:,:)
  real, allocatable :: xaxis(:), yaxis(:), tmp(:,:)
  logical, allocatable :: wet(:,:)
  integer :: i, j, t, taup1, tau, ki=81
  real :: rdelx, rdely, cfl, temp_b1, temp_b2, temp_a1, temp_a2, epsl=0.0
  real :: hmin = 0.05, tmpvel
  logical :: exist
  
  namelist/main_nml/ delt, delx, dely, ni, nj, ntstep, epsl

  inquire(file='input.nml',exist=exist)
  
  if (exist) then
     open(10, file='input.nml')
     read(10, nml=main_nml)
     close(10)
  endif
  
  write(*, nml=main_nml)
  
  nip1 = ni + 1; njp1 = nj + 1
  rdelx = 1./delx ; rdely = 1./dely
  temp_a1 = delt * GRAV * rdelx
  temp_a2 = delt * GRAV * rdely 
  temp_b1 = delt * rdelx * 0.5
  temp_b2 = delt * rdely * 0.5
  
  
  allocate(uvel(0:1,1:nip1,1:njp1))
  allocate(vvel(0:1,1:nip1,1:njp1))
  allocate(eta(0:1,0:nip1,0:nip1))
  allocate(h(0:nip1,0:njp1))
  allocate(h0(1:ni,1:nj))
  allocate(tmp(ni,nj))
  allocate(topo(1:ni,1:nj))
  allocate(wet(0:nip1,0:njp1))
  allocate(xaxis(nip1))
  allocate(yaxis(njp1))

  forall(i=1:nip1) xaxis(i) = i
  forall(i=1:njp1) yaxis(i) = i

  topo(:,:) = 100.0
!  topo(1,:) = -2.0
  topo(ni,:) = -2.0
  do i = 50, 250
     topo(i,:) = 50.0 - ((real(i)-50.0)/200.0) * 50.3
  enddo
  
  topo(251:300,:) = -0.3
  topo(:,1) = -2.0
  topo(:,nj) = -2.0

  uvel(:,:,:) = 0.0; vvel = 0.0

  h0(:,:) = topo(:,:)

  eta(:,:,:) = 0.0
  
  where (h0<0.0) eta(0,1:ni,1:nj) = -1.0*h0(:,:)

  do i = 1, ki
     eta(:,i,:) = 2.0*sin(3.1414*(i-1)/ki)
  enddo
  
  cfl = delt * sqrt(GRAV * maxval(h0)) / delx

  print *, 'cfl=',cfl
  if ( cfl >= 1.0 ) then
     print *, 'CFL Violation'
     stop 'ERROR'
  endif
  
  h(0,:) = hmin
  h(nip1,:) = hmin
  h(:,0) = hmin
  h(:,njp1) = hmin

  eta(:,0,:) = 0.0
  eta(:,nip1,:) = 0.0
  eta(:,:,0) = 0.0
  eta(:,:,njp1) = 0.0
     
  h(1:ni,1:nj) = h0(:,:) + eta(0,1:ni,1:nj)
  
  wet(:,:) = .true.
  
  where(h(1:ni,1:nj) < hmin) wet(1:ni,1:nj)=.false.
  
  do t = 1, ntstep
     taup1 = mod(t,2)
     tau = mod(t+1,2)

     do j = 1, njp1
        do i = 1, nip1
           tmpvel = uvel(tau,i,j) - temp_a1 * (eta(tau,i,j)-eta(tau,i-1,j))
           if (.not. wet(i-1,j) .and. tmpvel > 0.0) tmpvel = 0.0
           if (.not. wet(i,j) .and. tmpvel < 0.0) tmpvel = 0.0
           uvel(taup1,i,j) = tmpvel
           
           tmpvel = vvel(tau,i,j) - temp_a2 * (eta(tau,i,j)-eta(tau,i,j-1))
           if (.not. wet(i,j-1) .and. tmpvel > 0.0) tmpvel = 0.0
           if (.not. wet(i,j) .and. tmpvel < 0.0) tmpvel = 0.0
           vvel(taup1,i,j) = tmpvel
        enddo
     enddo
     
     do j = 1, nj
        do i = 1, ni
           eta(taup1,i,j) = eta(tau,i,j) - temp_b1 * ( (uvel(taup1,i+1,j) + abs(uvel(taup1,i+1,j))) * h(i,j) &
                                                     + (uvel(taup1,i+1,j) - abs(uvel(taup1,i+1,j))) * h(i+1,j) &
                                                     - (uvel(taup1,i,j) + abs(uvel(taup1,i,j))) * h(i-1,j) &
                                                     - (uvel(taup1,i,j) - abs(uvel(taup1,i,j))) * h(i,j) ) &
                                                     
                                         - temp_b2 * ( (vvel(taup1,i,j+1) + abs(vvel(taup1,i,j+1))) * h(i,j) &
                                                     + (vvel(taup1,i,j+1) - abs(vvel(taup1,i,j+1))) * h(i,j+1) &
                                                     - (vvel(taup1,i,j) + abs(vvel(taup1,i,j))) * h(i,j-1) &
                                                     - (vvel(taup1,i,j) - abs(vvel(taup1,i,j))) * h(i,j) )
        enddo
     enddo
     
     h(1:ni,1:nj) = h0(:,:) + eta(taup1,1:ni,1:nj)
      
     wet(:,:) = .true.
      
     where(h(1:ni,1:nj) < hmin) wet(1:ni,1:nj)=.false.

     call shapiro_filter_2d(eta(taup1,1:ni,1:nj),wet(1:ni,1:nj))

     h(1:ni,1:nj) = h0(:,:) + eta(taup1,1:ni,1:nj)
     
     tmp = -999.0
     where(wet(1:ni,1:nj)) tmp(:,:) = eta(taup1, 1:ni,1:nj)
     
     call write_nc(xaxis(1:ni), yaxis(1:nj), tmp, 'eta', t, missing_value=-999.0)
     call write_nc(xaxis(1:ni), yaxis(1:nj), h0, 'h0', t)
!     call write_nc(xaxis(1:ni), yaxis(1:nj), uvel(taup1,1:ni,1:nj), 'uvel', t, missing_value=-999.0)
!     call write_nc(xaxis(1:ni), yaxis(1:nj), vvel(taup1,1:ni,1:nj), 'vvel', t, missing_value=-999.0)

  enddo

  deallocate(uvel)
  deallocate(vvel)
  deallocate(eta)
  deallocate(h)
  deallocate(h0)
  deallocate(topo)
  deallocate(wet)

contains

  subroutine shapiro_filter_2d (field, mask)
    real, intent(inout) :: field(:,:)
    logical, intent(in) :: mask(:,:)
    integer :: is, ie, js, je, i, j
    
    is = 2; ie = size(field,1) - 1
    js = 2; je = size(field,2) - 1

    do i = is, ie
       do j = js, je, 1
          if (.not.mask(i-1,j)) cycle
          if (.not.mask(i+1,j)) cycle
          if (.not.mask(i,j-1)) cycle
          if (.not.mask(i,j+1)) cycle
          field(i,j) = (1-epsl)*field(i,j) + 0.25 * epsl * &
               (field(i-1,j)+field(i+1,j)+field(i,j-1)+field(i,j+1))
       end do
    end do
  end subroutine shapiro_filter_2d

end program main
