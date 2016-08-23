module shallow_water_mod
  
  use mpp_mod, only: mpp_npes, mpp_pe, mpp_error, stdout, FATAL, WARNING, NOTE, mpp_init, mpp_exit, mpp_max
  use mpp_io_mod, only: mpp_io_init, mpp_open, mpp_close, MPP_RDONLY, MPP_ASCII, MPP_MULTI
  use fms_mod,                    only: field_exist, field_size, read_data
  use mpp_domains_mod, only: domain2d, domain1d, mpp_define_layout, mpp_define_domains, &
       mpp_get_compute_domain, mpp_get_domain_components, mpp_update_domains
  use diag_manager_mod, only: diag_manager_init, register_diag_field, register_static_field, &
       diag_axis_init, send_data, diag_manager_end
  use time_manager_mod, only: set_calendar_type, NO_CALENDAR, time_type, set_time, operator(+), assignment(=), &
       print_time, set_ticks_per_second
  use constants_mod, only: GRAV, PI
  
  implicit none

  private 

  real :: delt, delx = 0.0, dely = 0.0
  
  integer :: ni = 0, nj = 0, halo = 1
  real :: dt_model=1.0 !seconds
  integer :: i, j, taup1, tau
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: tstep = 0
  integer :: id_x=0, id_y=0, id_uvel=0, id_vvel=0, id_eta=0, id_thick=0
  
  real :: rdelx, rdely,  b1,  b2,  a1,  a2, epsl_shapiro =0.0
  real :: hmin = 0.1
  type(time_type) :: dt_Time

  integer :: domain_layout(2)=(/1,1/)

  type(domain2D) :: domain2 
  
  real, allocatable :: uvel(:,:,:), vvel(:,:,:), eta(:,:,:), &
                       h(:,:), h0(:,:), topo(:,:)
  real, allocatable :: xaxis(:), yaxis(:), tmp(:,:)
  logical, allocatable :: wet(:,:)
  logical :: do_shapiro_filter = .false.

  character (len=1024) :: msg

  character(len=256) , parameter :: module_name='shallow_water_mod'
  
  namelist/shallow_water_nml/ dt_model, epsl_shapiro, do_shapiro_filter

  public :: init_model, update_model, end_model

contains

  !shapiro_filter_2d
  subroutine shapiro_filter_2d (field, mask)
    real, intent(inout) :: field(:,:)
    logical, intent(in) :: mask(:,:)
    integer :: is, ie, js, je, i, j

    if (.not.do_shapiro_filter) return
    
    is = 2; ie = size(field,1) - 1
    js = 2; je = size(field,2) - 1

    do i = is, ie
       do j = js, je, 1
          if (.not.mask(i-1,j)) cycle
          if (.not.mask(i+1,j)) cycle
          if (.not.mask(i,j-1)) cycle
          if (.not.mask(i,j+1)) cycle
          field(i,j) = (1-epsl_shapiro)*field(i,j) + 0.25 * epsl_shapiro * &
               (field(i-1,j)+field(i+1,j)+field(i,j-1)+field(i,j+1))
       end do
    end do
  end subroutine shapiro_filter_2d
  !----------------------------------------------------------------------------

  !Init_model
  subroutine init_model(Time)
    type(time_type), intent(in) :: Time
    integer :: siz(4), ibuf
    real :: cfl
    character (len=256) :: grid_file='INPUT/grid_spec.nc'
    integer :: id_topo = 0
    logical :: used
    integer :: seconds, ticks
    
    call mpp_init()
    call mpp_io_init()
    call diag_manager_init()
     
    call mpp_open(ibuf, 'input.nml', action=MPP_RDONLY,form=MPP_ASCII)
    rewind(ibuf)
    read(ibuf, nml=shallow_water_nml)
    call mpp_close(ibuf)
     
    write(stdout(), nml=shallow_water_nml)

    seconds = 0; ticks = 1
    
    if(int(dt_model)>0) then
       seconds = int(dt_model)
       call set_ticks_per_second(ticks)
       dt_Time = set_time(seconds=seconds)
       delt = real(int(dt_model))
    else
       ticks = int(1.0/dt_model)
       call set_ticks_per_second(ticks)
       dt_Time = set_time(seconds=seconds, ticks=1)
       delt = 1.0/ticks
    endif

    write(msg,*) delt 
    call mpp_error(NOTE, '=> Time step = '//trim(msg)//' seconds')
     
        
    if (.not.field_exist(grid_file, "topo")) call mpp_error(fatal, 'topo does not exist in '//trim(grid_file))
    call field_size(grid_file, "topo", siz)
    ni = siz(1)
    nj = siz(2)
    
    ! Make some checks
    if ( ni <= 0 .OR. nj <= 0) then
       call mpp_error(FATAL,'=>Error ni or nj is <= 0')
    endif

    allocate(xaxis(ni), yaxis(nj))
    
    call read_data(grid_file, 'x', xaxis, no_domain=.true.)
    call read_data(grid_file, 'y', yaxis, no_domain=.true.)
    call mpp_error(NOTE, '=>Read axis from grid_spec.')
    
    ! Inititializing some coeficients
    delx = xaxis(2)-xaxis(1)
    dely = yaxis(2)-yaxis(1)
    rdelx = 1./delx ; rdely = 1./dely
    
    a1 = delt * GRAV * rdelx
    a2 = delt * GRAV * rdely
    b1 = delt * rdelx * 0.5
    b2 = delt * rdely * 0.5

    
    ! Defining the domain decomposition layout. 
    call mpp_define_layout ((/1,ni,1,nj/),mpp_npes(),domain_layout)
     
    write(msg,*) 'domain layout',domain_layout
    call mpp_error(NOTE,trim(msg))
     
    ! Defining the local domain, with halo (for this example halo is set to 1).
    call mpp_define_domains( (/1,ni,1,nj/), domain_layout, domain2, &
         xhalo=halo, yhalo=halo )

    call mpp_get_compute_domain(domain2, isc, iec, jsc, jec)
    isd = isc - halo; ied = iec + halo
    jsd = jsc - halo; jed = jec + halo

    allocate(topo(isc:iec, jsc:jec))
    
    call read_data(grid_file, 'topo', topo, domain2)
    call mpp_error(NOTE, '=>Read topo from grid_spec')

    allocate( uvel(0:1,isd:ied,jsd:jed), vvel(0:1,isd:ied,jsd:jed) )
    allocate( eta(0:1,isd:ied,jsd:jed),  h(isd:ied,isd:ied) )
    allocate( h0(isc:iec,jsc:jec),       tmp(isc:iec,jsc:jec) )
    allocate( wet(isc:iec,jsc:jec) )
 
    uvel(:,:,:) = 0.0 ; vvel(:,:,:) = 0.0 ; eta(:,:,:) = 0.0
    h(:,:) = 0.0 ; h0(:,:) = topo(:,:) ; tmp(:,:) = 0.0
    wet(:,:) = .false.
    
    where (h0>hmin) wet = .true.

    where (h0(:,:) < 0.0) eta(0,isc:iec,jsc:jec) = -1.0*h0(:,:)
  
    call get_initial_conditions()

    h(1:ni,1:nj) = h0(:,:) + eta(0,1:ni,1:nj)
  
    wet(:,:) = .true.
  
    where(h(1:ni,1:nj) < hmin) wet(1:ni,1:nj)=.false.

    cfl = delt * sqrt(GRAV * maxval(h0)) / delx

    if ( cfl >= 1.0 ) then
       write(msg, *)'delx = ', delx, 'dely = ',dely, 'h0_max=', maxval(h0), 'delt=',delt, 'cfl=', cfl
       call mpp_error(WARNING, '=>'//trim(msg))
       call mpp_error(fatal, '=> CFL Violation')
    endif
    

    !diagnostics initializations
    
    id_x = diag_axis_init(name='X', data=xaxis, units='m' , cart_name='X', &
         long_name='X-axis', domain2=domain2)
    
    id_y = diag_axis_init(name='Y', data=yaxis, units='m' , cart_name='Y', &
         long_name='Y-axis', domain2=domain2)
    
    id_uvel = register_diag_field(module_name, 'uvel', (/id_x, id_y/), Time, &
         'U-velocity', 'm/s', interp_method="conserve_order1")
    
    id_vvel = register_diag_field(module_name, 'vvel', (/id_x, id_y/), Time, &
         'V-velocity', 'm/s', interp_method="conserve_order1")
    
    id_eta = register_diag_field(module_name, 'eta', (/id_x, id_y/), Time, &
         'Sea surface elevation', 'm', interp_method="conserve_order1", &
         mask_variant=.true., missing_value=1.e20)
    
    id_thick = register_diag_field(module_name, 'thick', (/id_x, id_y/), Time, &
         'Column Thickness', 'm', interp_method="conserve_order1")
    
    id_topo = register_static_field(module_name, 'topo', (/id_x, id_y/), &
         'Topography', 'm', interp_method="conserve_order1")

    used = send_data(id_topo, topo(isc:iec,jsc:jec), Time)    
  
  end subroutine init_model
  !----------------------------------------------------------------------------

  ! update_model
  subroutine update_model(Time)
    type(time_type), intent(inout) :: Time
    real :: tmpvel
    logical :: used

    tstep = tstep + 1
    Time = Time + dt_Time
    
    taup1 = mod(tstep,2)
    tau = mod(tstep+1,2)

    do i = isc, iec
       do j = jsc, jec
          tmpvel = uvel(tau,i,j) -  a1 * (eta(tau,i,j)-eta(tau,i-1,j))
          if (.not. wet(i-1,j) .and. tmpvel > 0.0) tmpvel = 0.0
          if (.not. wet(i,j) .and. tmpvel < 0.0) tmpvel = 0.0
          uvel(taup1,i,j) = tmpvel
          
          tmpvel = vvel(tau,i,j) -  a2 * (eta(tau,i,j)-eta(tau,i,j-1))
          if (.not. wet(i,j-1) .and. tmpvel > 0.0) tmpvel = 0.0
          if (.not. wet(i,j) .and. tmpvel < 0.0) tmpvel = 0.0
          vvel(taup1,i,j) = tmpvel
       enddo
    enddo
    
    do i = isc, iec 
       do j = jsc, jec
          eta(taup1,i,j) = eta(tau,i,j) -  b1 * ( (uvel(taup1,i+1,j) + abs(uvel(taup1,i+1,j))) * h(i,j) &
                                                    + (uvel(taup1,i+1,j) - abs(uvel(taup1,i+1,j))) * h(i+1,j) &
                                                    - (uvel(taup1,i,j) + abs(uvel(taup1,i,j))) * h(i-1,j) &
                                                    - (uvel(taup1,i,j) - abs(uvel(taup1,i,j))) * h(i,j) ) &
                                                    
                                        -  b2 * ( (vvel(taup1,i,j+1) + abs(vvel(taup1,i,j+1))) * h(i,j) &
                                                    + (vvel(taup1,i,j+1) - abs(vvel(taup1,i,j+1))) * h(i,j+1) &
                                                    - (vvel(taup1,i,j) + abs(vvel(taup1,i,j))) * h(i,j-1) &
                                                    - (vvel(taup1,i,j) - abs(vvel(taup1,i,j))) * h(i,j) )
       enddo
    enddo
    
    h(isc:iec,jsc:jec) = h0(isc:iec,jsc:jec) + eta(taup1,isc:iec,jsc:jec)
     
    wet(:,:) = .true.
     
    where(h(isc:iec,jsc:jec) < hmin) wet(isc:iec,jsc:jec)=.false.

    call shapiro_filter_2d (eta(taup1,isc:iec,jsc:jec),wet(isc:iec,jsc:jec))

    h(isc:iec,jsc:jec) = h0(isc:iec,jsc:jec) + eta(taup1,isc:iec,jsc:jec)

    !Diagnostics
        
    used = send_data(id_eta, eta(taup1,isc:iec, jsc:jec), Time, &
         mask=wet(isc:iec,jsc:jec))
    
    used = send_data(id_uvel, uvel(taup1, isc:iec, jsc:iec), Time)
    
    used = send_data(id_vvel, vvel(taup1, isc:iec, jsc:iec), Time)

    used = send_data(id_thick, h(isc:iec,jsc:jec), Time) 
  
  end subroutine update_model
  !----------------------------------------------------------------------------------------------------------

  
  ! get_initial_conditions 
  subroutine get_initial_conditions()
    integer :: i, j
    do i = 1, 21
       do j = 1, 21
          eta(:,i+20,j+20) = sin(PI*(i-1)/21)*sin(PI*(j-1)/21)
       enddo   
    enddo
    
  end subroutine get_initial_conditions
  !---------------------------------------------------------------------------------------------------------

  !end_model
  subroutine end_model()
    
    deallocate(uvel)
    deallocate(vvel)
    deallocate(eta)
    deallocate(h)
    deallocate(h0)
    deallocate(topo)
    deallocate(wet)

  end subroutine end_model
  !---------------------------------------------------------------------------------------------------------
    
end module shallow_water_mod

#ifdef test_shallow_water_mod
program main
  
  use mpp_mod, only: mpp_init, mpp_exit 
  use diag_manager_mod, only: diag_manager_init, diag_manager_end
  use shallow_water_mod, only : init_model, update_model, end_model
  use time_manager_mod, only: set_calendar_type, NO_CALENDAR, time_type, set_time, operator(+), assignment(=), &
       print_time

  implicit none

  integer :: ntstep = 1000, t
  type(time_type) :: Time

  call mpp_init()

  call diag_manager_init()

  call set_calendar_type(NO_CALENDAR)
  
  Time = set_time(seconds=0)
  
  call init_model(Time)

  do t = 1, ntstep
     call update_model(Time)
  end do

  call end_model()

  call diag_manager_end(Time)

  call mpp_exit()
  
end program main
#endif 
