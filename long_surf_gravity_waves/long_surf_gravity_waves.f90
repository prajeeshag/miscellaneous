program shallow_water
  
  use mpp_mod, only: mpp_npes, mpp_pe, mpp_error, stdout, FATAL, WARNING, NOTE, mpp_init, mpp_exit, mpp_max
  use mpp_io_init, only: mpp_io_init, mpp_open, mpp_close, MPP_RDONLY, MPP_ASCII, MPP_MULTI
  use mpp_domains_mod, only: domain2d, domain1d, mpp_define_layout, mpp_define_domains, &
       mpp_get_compute_domain, mpp_get_domain_components, mpp_update_domains
  use diag_manager_mod, only: diag_manager_init, register_diag_field, register_static_field, &
       diag_axis_init, send_data, diag_manager_end
  use time_manager_mod, only: set_calendar_type, NO_CALENDAR, time_type, set_time, operator(+), assignment(=), &
       print_time
  
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
  real :: ibuf
  
  namelist/shallow_water_nml/ delt, ntstep, epsl

  call init_model()
  
  cfl = delt * sqrt(GRAV * maxval(h0)) / delx

  print *, 'cfl=',cfl
  if ( cfl >= 1.0 ) then
     print *, 'CFL Violation'
     stop 'ERROR'
  endif
  
  
  temp_a1 = delt * GRAV * rdelx
  temp_a2 = delt * GRAV * rdely
  
  temp_b1 = delt * rdelx * 0.5
  temp_b2 = delt * rdely * 0.5
  
  do t = 1, ntstep
     taup1 = mod(t,2)
     tau = mod(t+1,2)

     do i = isc, iec
        do j = jsc, jec
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
     
     do i = isc, iec 
        do j = jsc, jec
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

  subroutine init_model()
    integer :: siz(2)

    call mpp_init()
    call mpp_io_init()
    call set_calendar_type(NO_CALENDAR)
    call diag_manager_init()
     
    call mpp_open(ibuf, 'input.nml', action=MPP_RDONLY,form=MPP_ASCII)
    rewind(ibuf)
    read(ibuf, nml=shallow_water_nml)
    call mpp_close(ibuf)
     
    write(stdout(), nml=shallow_water_nml)
     
    Time = set_time(seconds=0)
  
    if (.not.field_exist(grid_file, "topo")) call mpp_error(fatal, 'topo does not exist in '//trim(grid_file))
    call field_size(grid_file, "topo", siz)
    ni = siz(1)
    nj = siz(2)
    
    ! Make some checks
    if ( ni <= 0 .OR. nj <= 0) then
       call mpp_error(FATAL,'=>Error ni or nj is <= 0')
    endif

    allocate(xaxis(ni), yaxis(nj))

    delx = xaxis(2)-xaxis(1)
    dely = yaxis(2)-yaxis(1)
    rdelx = 1./delx ; rdely = 1./dely
    
    call read_data(grid_file, 'x', xaxis, no_domain=.true.)
    call read_data(grid_file, 'y', yaxis, no_domain=.true.)
    
    ! Defining the domain decomposition layout. 
    call mpp_define_layout ((/1,ni,1,nj/),mpp_npes(),domain_layout)
     
    write(msg,*) 'domain layout',domain_layout
    call mpp_error(NOTE,trim(msg))
     
    ! Defining the local domain, with halo (for this example halo is set to 1).
    call mpp_define_domains( (/1,ni,1,nj/), domain_layout, local_domain, &
         xhalo=halo, yhalo=halo )

    call mpp_get_compute_domain(local_domain, isc, iec, jsc, jec)
    isd = isc - halo; ied = iec + halo
    jsc = jsc - halo; jed = jec + halo

    allocate(topo(isc:iec, jsc:jec))
    
    call read_data(grid_file, 'topo', topo, local_domain)

    allocate( uvel(0:1,isd:ied,jsd:jed), vvel(0:1,isd:ied,jsd:jed) )
    allocate( eta(0:1,isd:ied,jsd:jed),  h(isd:ied,isd:ied) )
    allocate( h0(isc:iec,jsc:jec),       tmp(isc:iec,jsc:jec) )
    allocate( wet(isc:iec,jsc:jec) )
 
    uvel(:,:,:) = 0.0 ; vvel(:,:,:) = 0.0 ; eta(:,:,:) = 0.0
    h(:,:) = 0.0 ; h0(:,:) = topo(:,:) ; tmp(:,:) = 0.0
    wet(:,:) = .false.
    
    where(h0>hmin) wet = .true.

    where (h0(:,:) < 0.0) eta(0,isc:iec,jsc:jec) = -1.0*h0(:,:)
  
    call get_initial_conditions()

    h(1:ni,1:nj) = h0(:,:) + eta(0,1:ni,1:nj)
  
    wet(:,:) = .true.
  
    where(h(1:ni,1:nj) < hmin) wet(1:ni,1:nj)=.false.
    
  end subroutine init_model
    
end program shallow_water
