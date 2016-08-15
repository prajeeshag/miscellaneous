program main
  use netcdf_write, only : write_nc
  
  implicit none
  real, parameter :: GRAV = 9.8
  real :: delt = 1.0, delx = 10.0, dely = 10.0
  integer :: ntstep=2000
  integer :: ni = 100, nj = 100, nip1, njp1
  real, allocatable :: uvel(:,:,:), vvel(:,:,:), eta(:,:,:), h(:,:), h0(:,:), topo(:,:)
  real, allocatable :: xaxis(:), yaxis(:), tmp(:,:)
  logical, allocatable :: maskt(:,:)
  integer :: i, j, t, taup1, tau 
  real :: rdelx, rdely, cfl, temp_b1, temp_b2, temp_a1, temp_a2
  logical :: exist
  
  namelist/main_nml/ delt, delx, dely, ni, nj

  inquire(file='input.nml',exist=exist)
  
  if (exist) then
     open(10, file='input.nml')
     read(10, nml=main_nml)
     close(10)
  endif
  
  write(*, nml=main_nml)
  
  nip1 = ni + 1; njp1 = nj + 1
  rdelx = 1./delx ; rdely = 1./dely
  
  
  allocate(uvel(0:1,1:nip1,1:njp1))
  allocate(vvel(0:1,1:nip1,1:njp1))
  allocate(eta(0:1,1:ni,1:nj))
  allocate(h(0:nip1,0:njp1))
  allocate(h0(1:ni,1:nj))
  allocate(tmp(ni,nj))
  allocate(topo(1:ni,1:nj))
  allocate(maskt(0:nip1,0:njp1))
  allocate(xaxis(nip1))
  allocate(yaxis(njp1))

  forall(i=1:nip1) xaxis(i) = i
  forall(i=1:njp1) yaxis(i) = i

  topo(:,:) = -10.0
  topo(1,:) = 5.0
  topo(ni,:) = 5.0
  topo(:,1) = 5.0
  topo(:,nj) = 5.0

  maskt(:,:) = .false.

  h0(:,:) = 0.0; eta(:,:,:) = 0.0; uvel(:,:,:) = 0.0; h(:,:) = 0.0; vvel = 0.0
  
  where(topo(:,:) < 0.0) maskt(1:ni,1:nj)=.true.
  
  eta(:,int(ni/2)-5:int(ni/2)+5,int(nj/2)-5:int(nj/2)+5) = 0.5
  
  where (maskt(1:ni,1:nj)) h0(:,:) = -1.0*topo(:,:)
  
  temp_a1 = delt * GRAV * rdelx
  temp_a2 = delt * GRAV * rdely 
  temp_b1 = delt * rdelx * 0.5
  temp_b2 = delt * rdely * 0.5
  
  h(1:ni,1:nj) = h0(:,:) + eta(0,:,:)
  
  cfl = delt * sqrt(GRAV * maxval(h)) / delx

  print *, 'cfl=',cfl
  if ( cfl >= 1.0 ) then
     print *, 'CFL Violation'
     stop 'ERROR'
  endif
  
  do t = 1, ntstep
     taup1 = mod(t,2)
     tau = mod(t+1,2)

     do j = 1, njp1
        do i = 1, nip1
           if (maskt(i,j).and.maskt(i-1,j)) then
              uvel(taup1,i,j) = uvel(tau,i,j) - temp_a1 * (eta(tau,i,j)-eta(tau,i-1,j))
           endif
           if (maskt(i,j).and.maskt(i,j-1)) then
              vvel(taup1,i,j) = vvel(tau,i,j) - temp_a2 * (eta(tau,i,j)-eta(tau,i,j-1))
           endif
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
     
     h(1:ni,1:nj) = h0(:,:) + eta(taup1,:,:)
     tmp = -999.0
     where(h(1:ni,1:nj)>0.0) tmp(:,:) = h(1:ni,1:nj)
     call write_nc(xaxis(1:ni), yaxis(1:nj), tmp, 'h', t, missing_value=-999.0)
     call write_nc(xaxis(1:ni), yaxis(1:nj), uvel(taup1,1:ni,1:nj), 'uvel', t, missing_value=-999.0)
     call write_nc(xaxis(1:ni), yaxis(1:nj), vvel(taup1,1:ni,1:nj), 'vvel', t, missing_value=-999.0)

  enddo

  deallocate(uvel)
  deallocate(vvel)
  deallocate(eta)
  deallocate(h)
  deallocate(h0)
  deallocate(topo)
  deallocate(maskt)

end program main
