!======================================!
!  Program bubble_test                 !
!                                      !
!  Computes effective buoyancy of 3D   !
!  rho field and writes to file        !
!                                      !
!  Input: .nc filepath,                !
!  Output:                             !
!  rho recalculated and appended to    !
! original netcdf file.                !
!                                      !
!  Copyright (c) Nadir Jeevanjee       !
!  http://www.romps.org                !
!======================================!

program bubble_test

   use buoyancy_mod
   use netcdf
   use netcdf_mod
   use, intrinsic :: iso_c_binding

   implicit none

   include '/home1/02291/jeevanje/domain/fftw-3.3.3/install/include/fftw3.f03'

   ! In
   character(200) :: filepath

   ! Privates
   integer :: nx, ny, nz, r
   integer :: xdimid, ydimid, zdimid, tdimid, ncid, betaid, status, deltarhoid, deltarho_ftid, deltarho_altid
   real(c_double), dimension(:,:,:), allocatable :: rho, beta, deltarho, deltarho_alt, deltarhohat, deltarho_ft
   real(c_double), dimension(:), allocatable :: x,y,z
   real(c_double) :: dx, dy
   logical :: file_exist
   ! For FFT
   real(c_double), dimension(:,:), allocatable :: in
   complex(c_double_complex), dimension(:,:), allocatable :: out 
   type(c_ptr) :: planf, planb

   ! Read in filepath and t index
   call getarg(1,filepath)

   ! Check that file exists
   inquire(file=trim(filepath),exist=file_exist)
         if (.not.file_exist) then
            print *, 'Error in bubble_test: Can not find '//trim(filepath)
          stop
         end if
         if (file_exist) then
            print *, 'file '//trim(filepath)//' found.'
         end if

  ! Open file 
   call handle_err(nf90_open(path = filepath, mode = nf90_write, ncid = ncid))

   ! Retrieve dimensions
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'x', dimid = xdimid))
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'y', dimid = ydimid))
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'z', dimid = zdimid))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = xdimid, len = nx))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = ydimid, len = ny))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = zdimid, len = nz))
   write(*,*) 'nx,ny,nz =', nx, ny, nz

   allocate(rho(nx,ny,nz),&
        beta(nx,ny,nz), deltarho(nx,ny,nz),&
        deltarhohat(nx/2+1,ny,nz), & 
        deltarho_ft(nx,ny,nz), & 
        x(nx), y(ny), z(nz), &
        in(nx,ny), out(nx/2 +1,ny) )

   print *, 'Arrays allocated'

   ! Get input data
   x = get_netCDF1(trim(filepath),'x')
   y = get_netCDF1(trim(filepath),'y')
   z = get_netCDF1(trim(filepath),'z')
   dx = x(2)-x(1)
   dy = y(2)-y(1)
   rho = get_netCDF3(trim(filepath),'rho') 

   ! compute delta rho
   call laplacian2d(rho,dx,dy,nz,deltarho)

   print *, "Laplacian of rho computed"

   ! Plans for FFTW
   planf = fftw_plan_dft_r2c_2d(ny,nx,in,out,fftw_measure) 
   planb = fftw_plan_dft_c2r_2d(ny,nx,out,in,fftw_measure) 

   ! Compute FT of deltarho for check
   do r = 1,nz
      in = deltarho(:,:,r)
      call fftw_execute_dft_r2c(planf,in,out)
      deltarhohat(:,:,r) = out
   end do

   ! Backwards FT
   do r = 1,nz
      out = deltarhohat(:,:,r)
      call fftw_execute_dft_c2r(planb,out,in)
      deltarho_ft(:,:,r) = in/(nx*ny)
   end do

   ! Solve poisson eqn with g*deltarho as source
   call solve_poisson(-gr*deltarho,x,y,z,beta)

   print *, "Poisson equation solved"

   ! Compute 3D Laplacian of beta as check
   call laplacian3d(beta,dx,dy,z,deltarho_alt)
   deltarho_alt = -1/gr*deltarho_alt
   
   !=================
   ! Write output
   !=================

   print *, "Calculation complete, writing output to netCDF file"

   ! Write variables to nc file
   call write_netCDF3(ncid,(/ xdimid, ydimid, zdimid/),'beta','kg/(m^2 s^2)',&
        'Effective Buoyancy', beta)
   call write_netCDF3(ncid,(/ xdimid, ydimid, zdimid/),'deltarho','kg/m^5',&
        'Horizontal laplacian of rho', deltarho)
   call write_netCDF3(ncid,(/ xdimid, ydimid, zdimid/),'deltarho_ft','kg/m^5',&
        'Doubly transformed horizontal laplacian of rho', deltarho_ft)
   call write_netCDF3(ncid,(/ xdimid, ydimid, zdimid/),'deltarho_alt','kg/m^5',&
        'Horizontal laplacian of rho,recomputed from beta', deltarho_alt)

   ! Close nc file
   call handle_err(nf90_close(ncid = ncid))      

end program bubble_test
