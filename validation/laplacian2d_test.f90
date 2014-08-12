!======================================!
!  Program fftw_test                   !
!                                      !
! Code to calculate FFT of 2-D arrays  !
! and check proper operation of FFTW   !
! code                                 !
!                                      !
!  Input: .nc filepath,                !
!     z and t indices for 2-D slab     !
!  Output:                             !
!  rho recalculated and appended to    !
! original netcdf file.                !
!                                      !
!  Copyright (c) Nadir Jeevanjee       !
!  http://www.romps.org                !
!======================================!

program fftw_test

   use buoyancy_mod
   use netcdf
   use netcdf_mod
   use, intrinsic :: iso_c_binding

   implicit none

   include '/home1/02291/jeevanje/domain/fftw-3.3.3/install/include/fftw3.f03'

   ! In
   character(200) :: filepath, zarg, targ 

   ! Privates
   integer :: nx, ny, nz, nt, zind, tind
   integer :: xdimid, ydimid, zdimid, tdimid, ncid, rhoid, status, deltarhoid
   real(c_double), dimension(:,:), allocatable :: rho2d, rho2d_alt, deltarho
   real(c_double), dimension(:,:,:,:), allocatable :: rho4d
   real(c_double), dimension(:), allocatable :: x,y
   real(c_double) :: dx, dy
   complex(c_double_complex), dimension(:,:), allocatable :: out 
   type(c_ptr) :: planf, planb
   logical :: file_exist


   ! Read in filepath and z, t indices
   call getarg(1,filepath)
   call getarg(2,zarg)
   call getarg(3,targ)
   read(zarg,*), zind
   read(targ,*), tind


   ! Check that file exists
   inquire(file=trim(filepath),exist=file_exist)
         if (.not.file_exist) then
            print *, 'Error in fftw_test: Can not find '//trim(filepath)
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
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'time', dimid = tdimid))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = xdimid, len = nx))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = ydimid, len = ny))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = zdimid, len = nz))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = tdimid, len = nt))
   write(*,*) 'nx,ny,nz,nt =', nx, ny, nz, nt


   allocate(rho4d(nx,ny,nz,nt),rho2d(nx,ny), out(ny,nx/2+1),rho2d_alt(nx,ny), & 
            deltarho(nx,ny), x(nx), y(ny))
   !planf = fftw_plan_dft_r2c_2d(ny,nx,rho2d,out,fftw_measure)
   !planb = fftw_plan_dft_c2r_2d(ny,nx,out,rho2d_alt,fftw_measure)
   print *, 'Arrays allocated'

   ! Get input data
   x = get_netCDF1(trim(filepath),'x')
   y = get_netCDF1(trim(filepath),'y')
   dx = x(2)-x(1)
   dy = y(2)-y(1)
   rho4d = get_netCDF4(trim(filepath),'rho') 
   rho2d = rho4d(1:nx,1:ny,zind,tind)

   ! compute delta rho
   call laplacian2d(rho2d,dx,dy,deltarho)

   !call fftw_execute_dft_r2c(planf,rho2d,out)
   !call fftw_destroy_plan(planf)

   !write (*,*) 'FFT =', out 

   !call fftw_execute_dft_c2r(planb,out,rho2d_alt)
   !rho2d_alt = rho2d_alt/(nx*ny) ! Renormalize
   !call fftw_destroy_plan(planb)

   ! Write deltarho to nc file
   call handle_err(nf90_redef(ncid = ncid)) ! Put file in define mode
   call handle_err(nf90_def_var(ncid = ncid, name = 'deltarho', &
        xtype = nf90_float, dimids = (/ xdimid, ydimid/), &
        varid = deltarhoid ) ) ! Get varid

   ! Add attributes
   call handle_err(nf90_put_att(ncid, deltarhoid, 'units', 'kg/m^5'))
   call handle_err(nf90_put_att(ncid, deltarhoid, 'long_name', 'Horizontal Laplacian of density at z=475 m'))
   call handle_err(nf90_enddef(ncid = ncid)) ! End define mode
 
   ! Write data
   call handle_err(nf90_put_var(ncid,deltarhoid,deltarho))
   print *, "deltarho written"

   call handle_err(nf90_close(ncid = ncid))      

end program fftw_test
