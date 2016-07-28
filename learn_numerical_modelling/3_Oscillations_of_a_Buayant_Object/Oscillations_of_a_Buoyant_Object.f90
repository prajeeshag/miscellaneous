program main
	implicit none
	
	real :: depth=-100.0, rho_surf=1025.0, rho_bot=1026.0, rho_obj=1025.5
	real :: dt=0.1, w=0.0, z=-80.0, fric_coef=0.001, t=0, g=9.8, rho_amb
	real :: tot_t=600.0, alfa=1.0, w1
	logical :: file_opened=.false., two_layer=.false.

	namelist/input_nml/depth, rho_surf, rho_bot, rho_obj, dt, z, fric_coef, tot_t, alfa, two_layer

	open(10,file='input.nml')
		read(10,nml=input_nml)
		write(*,nml=input_nml)
	close(10)

	do while (t <= tot_t)
		t = t + dt

		call update_w()

		call update_z()

		call print_out()
		
	end do

	close(10)

	call system('./plot.py')

	contains
	
	subroutine update_w()
		if (two_layer) then
			rho_amb = rho_surf
			if ((z/depth) > 0.5) rho_amb = rho_bot
		else
			rho_amb = rho_surf + ((rho_bot - rho_surf) / depth) * z
		end if
		if (quadratic_fric) then
			w = (w - dt * (g * ( 1.0 - rho_amb/rho_obj ) - fric_coef * (1-alfa) * w )) / (1 + fric_coef * dt * alfa)
		else
			w = (w - dt * (g * ( 1.0 - rho_amb/rho_obj ) - fric_coef * (1-alfa) * w)) / (1 + fric_coef * dt * alfa)
		end if
	end subroutine

	subroutine update_z()
	  z = z + dt * w
		z = min(z,0.0)
		z = max(z,depth)
	end subroutine

	subroutine print_out()
		
		if (.not.file_opened) then 
			open(10,file='output.dat')
			file_opened = .true.
		endif
		
		write(10,*) t, z, w
	
	end subroutine 
		
end program main
