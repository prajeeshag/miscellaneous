! Swathi..Test program to solve u<sub>t = delsq(u)
!
! with zero bc's and sin(pi*x)*sin(pi*y) initial condition

PROGRAM  heat_eqn
  use mpp_mod, only: mpp_npes, mpp_pe, mpp_error, stdout, FATAL, WARNING, NOTE, mpp_init, mpp_exit
  use mpp_mod, only: mpp_max

  use mpp_io_mod, only: mpp_io_init, mpp_open, mpp_close, MPP_RDONLY, MPP_ASCII, MPP_MULTI
                     
  use mpp_domains_mod, only: domain2d, domain1d, mpp_define_layout, mpp_define_domains, &
                             mpp_get_compute_domain, mpp_get_domain_components, mpp_update_domains

  use diag_manager_mod, only: diag_manager_init, register_diag_field, register_static_field, &
                              diag_axis_init, send_data, diag_manager_end

  use time_manager_mod, only: set_calendar_type, NO_CALENDAR, time_type, set_time, operator(+), assignment(=), &
                              print_time

  implicit none

  integer :: ni=100, nj=100, halo = 1, n, unit, ntsteps=10000
  real :: epsl = 1e-5, depsl
  real, dimension(:),allocatable :: xt0,yt0
  real, dimension(:),allocatable :: xt1,yt1
  real, allocatable,dimension(:,:,:) :: u
  integer :: dtts = 1
  real :: delx=2.0, dely=2.0, cfl
  integer :: ioun,npes,pe, isc,iec,jsc,jec,isd,ied,jsd,jed
  integer, dimension(2) :: domain_layout=(/1,1/)
        integer :: tau=0,taup1=1

  type(domain2D) :: local_domain   ! domain with halos on this processor

  type(time_type) :: Time, timestep

  integer :: i,j,io_status,io_int, id_temp, id_y, id_x

  logical :: used

  integer :: bnd_option = 1, init_option = 2 

  namelist /heat_eqn_nml/ ni, nj, dtts, ntsteps, delx, dely, bnd_option, epsl

  depsl = huge(depsl)
  call mpp_init()
  call mpp_io_init()
  call set_calendar_type(NO_CALENDAR)
  call diag_manager_init()

  Time = set_time(seconds=0)
 
  npes = mpp_npes()
  pe = mpp_pe()


! Read some namelist inputs
  call mpp_open(ioun,'input.nml',action=MPP_RDONLY,form=MPP_ASCII)
  read(ioun,heat_eqn_nml,iostat=io_status)
  write(stdout(),heat_eqn_nml)
  call mpp_close(ioun)

! Make some checks
  if ( nj <= 0 .OR. nj <= 0 .OR. dtts <= 0 .OR. ntsteps <= 0) then
          write(stdout(),*)'nj,nj,dtts,ntsteps',nj,nj,dtts,ntsteps
          call mpp_error(FATAL,'=>Error one of the above is <= 0')
  endif

! Do a cfl check
  cfl = 0.5*delx**2*dely**2/(delx**2+dely**2)
  write(stdout(),*)'cfl=',cfl
  if ( dtts > cfl)  call mpp_error(FATAL,'CFL violation')

! timestep is an derived data type defined by time_manager_mod of fms, various
! mathematical operation such as + - > < = is also defined for this time data type.
  timestep = set_time(seconds=dtts)

  allocate (xt0(ni))
  allocate (yt0(nj))

  allocate (xt1(ni))
  allocate (yt1(nj))

  do i = 1,ni
    xt0(i) = i*delx
    xt1(i) = real(i)/(ni+1)
  enddo

  do j = 1,nj
    yt0(j) = j*dely
    yt1(j) = real(j)/(nj+1)
  enddo

! Defining the domain decomposition layout.  mpp_define_layout will give the decomposition layout into a most optimal configuaration. 
  call mpp_define_layout ((/1,ni,1,nj/),npes,domain_layout)

  if (pe ==0) then
    write(stdout(),*)'domain layout',domain_layout
    write(stdout(),*)'ni,nj',ni,nj
  endif

! Defining the local domain, with halo (for this example halo is set to 1).
  call mpp_define_domains( (/1,ni,1,nj/), domain_layout, local_domain, &
    xhalo=halo, yhalo=halo )

! isc = start index of local compute domain in i direction. (compute domain : local domain excluding the halo boundaries)
! iec = end index of local compute domain in i direction.
! jsc = start index of local compute domain in j direction.
! jec = end index of local compute domain in j direction.

! isd = start index of local full domain in i direction. (full domain : local domain including the halo boundaries)
! ied = end index of local full domain in i direction.
! jsd = start index of local full domain in j direction.
! jed = end index of local full domain in j direction.
  call mpp_get_compute_domain (local_domain, isc, iec, jsc, jec)
        isd = isc - halo; ied = iec + halo
        jsd = jsc - halo; jed = jec + halo

  write(stdout(),*)'pe,isc,iec,jsc,jec',pe,isc,iec,jsc,jec

  allocate(u(isd:ied,jsd:jed,0:1)) ! allocated as local domain size including halos, 
                                   ! but the computation happens only for the compute domain.

  u(:,:,:) = 0.0

! Seting Up Output using diag_manager
  ! register the netcdf axis X and Y
  id_x = diag_axis_init(name='X', data=xt0(isc:iec), units='m' , cart_name='X', long_name='X-axis', domain2=local_domain)
  id_y = diag_axis_init(name='Y', data=yt0(jsc:jec), units='m' , cart_name='Y', long_name='Y-axis', domain2=local_domain)

  ! register a variable (here name is temp)
  id_temp = register_diag_field ( module_name='heat_eqn', field_name='temp', axes=(/id_x,id_y/), init_time=Time, &
              long_name='Temperature', units='Degree Celsius', interp_method = "conserve_order1" )

! Initial conditions
  call get_initial_conditions()

! halo update using mpp_update_domains
  call mpp_update_domains(u,local_domain)

! Prescribing boundary conditions at the boundary of full domain.

  call generate_boundary_conditions()

! Time Step loop
  do n = 1,ntsteps 
    tau = mod(n-1,2)
    taup1 = mod(n,2)
    ! Time is incremented with the timestep
    Time = Time + timestep

    do j = jsc,jec
      do i = isc,iec
        u(i,j,taup1) = u(i,j,tau)*(1-2*dtts*(delx**2 + dely**2)/(delx**2*dely**2)) + & 
             dtts*((u(i-1,j,tau) + u(i+1,j,tau))/delx**2 + (u(i,j-1,tau) + u(i,j+1,tau))/dely**2)
        depsl = max(u(i,j,taup1)-u(i,j,tau),depsl)
      enddo
   enddo
   call mpp_max(depsl)
   
    ! Halo update after each time step computation.
    call mpp_update_domains(u,local_domain)
    
    ! data for the netcdf variable temp is being send at every time-step
    used = send_data(id_temp, u(isc:iec,jsc:jec,taup1), Time)
    call print_time(Time)
    if(depsl<=epsl) exit
  enddo

  ! It is very important to end the diag_manager, other wise the data at the
  ! last few times won't get flushed on to the disk.
  call diag_manager_end(Time)

  !
  call mpp_exit()
  contains
    subroutine generate_boundary_conditions ()
      select case (bnd_option)
         case(1)
            if (jsc == 1) u(:,jsc-1,0) = 0.0
            if (jec == nj) u(:,jec+1,0) = 0.0
            if (isc == 1) u(isc-1,:,0) = 0.0
            if (iec == ni) u(iec+1,:,0) = 0.0
         case(2)
            if (jsc == 1) u(:,jsc-1,0) = 1.0
            if (jec == nj) u(:,jec+1,0) = 0.0
            if (isc == 1) u(isc-1,:,0) = 0.0
            if (iec == ni) u(iec+1,:,0) = 0.0
         case(3)
            if (jsc == 1) u(:,jsc-1,0) = 1.0
            if (jec == nj) u(:,jec+1,0) = 0.0
            if (isc == 1) u(isc-1,:,0) = 1.0
            if (iec == ni) u(iec+1,:,0) = 0.0
         case(4)
            if (jsc == 1) u(:,jsc-1,0) = 1.0
            if (jec == nj) u(:,jec+1,0) = 1.0
            if (isc == 1) u(isc-1,:,0) = 1.0
            if (iec == ni) u(iec+1,:,0) = 0.0
         case(5)
            if (jsc == 1) u(:,jsc-1,0) = 1.0
            if (jec == nj) u(:,jec+1,0) =1.0
            if (isc == 1) u(isc-1,:,0) = 0.0
            if (iec == ni) u(iec+1,:,0) = 0.0
         case default
            call mpp_error(FATAL,'Invalid option for boundary value.')
       end select

    end subroutine generate_boundary_conditions
    
    subroutine get_initial_conditions()
      select case (init_option) 
         case (1)
            u(isc:iec,jsc:jec,0) = 0.0
         case default
            do j = jsc,jec
               do i = isc, iec
                  u(i,j,0) = sin(3.1414*xt1(i))*sin(3.1414*yt1(j))
               enddo
            enddo
      end select
    end subroutine get_initial_conditions

END PROGRAM heat_eqn
