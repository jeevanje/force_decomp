!======================================!
!  Program poisson_solve_test          !
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

program poisson_solve_test

   use buoyancy_mod
   use calculus_mod
   use netcdf
   use netcdf_mod
   use, intrinsic :: iso_c_binding

   implicit none

   include '/home1/02291/jeevanje/domain/fftw-3.3.3/install/include/fftw3.f03'

   ! In
   character(200) :: filepath, zarg, targ 

   ! Privates
   integer :: nx, ny, nz, nt, zind, tind, r
   integer :: xdimid, ydimid, zdimid, tdimid, ncid, betaid, status, deltarhoid, deltarho_ftid, deltarho_altid
   real(c_double), dimension(:,:,:), allocatable :: rho, beta, deltarho, deltarho_alt, deltarho_ft
   complex(c_double_complex), dimension(:,:,:), allocatable :: deltarhohat
   real(c_double), dimension(:,:,:,:), allocatable :: rho4d
   real(c_double), dimension(:), allocatable :: x,y,z
   real(c_double) :: dx, dy
   logical :: file_exist

   ! For FFT
   real(c_double), dimension(:,:), allocatable :: in
   complex(c_double_complex), dimension(:,:), allocatable :: out 
   type(c_ptr) :: planf, planb

   ! Read in filepath and t index
   call getarg(1,filepath)
   call getarg(2,zarg)
   read(zarg,*), zind
   call getarg(3,targ)
   read(targ,*), tind

   ! Check that file exists
   inquire(file=trim(filepath),exist=file_exist)
         if (.not.file_exist) then
            print *, 'Error in poisson_solve_test: Can not find '//trim(filepath)
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

   allocate(rho4d(nx,ny,nz,nt),rho(nx,ny,nz),&
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
   rho4d = get_netCDF4(trim(filepath),'rho') 
   rho = rho4d(:,:,:,tind)

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

   call handle_err(nf90_redef(ncid = ncid)) ! Put file in define mode

   ! Write beta to nc file
   call handle_err(nf90_def_var(ncid = ncid, name = 'beta', &
        xtype = nf90_float, dimids = (/ xdimid, ydimid, zdimid/), &
        varid = betaid ) ) ! Get varid

   ! Add attributes
   call handle_err(nf90_put_att(ncid, betaid, 'units', 'kg/(m^2 s^2)'))
   call handle_err(nf90_put_att(ncid, betaid, 'long_name', 'Effective Buoyancy'))
   call handle_err(nf90_enddef(ncid = ncid)) ! End define mode
 
   ! Write beta
   call handle_err(nf90_put_var(ncid,betaid,beta))
   print *, "beta written"

!!!!!!!!!!!

   ! Write deltarho to nc file
   call handle_err(nf90_redef(ncid = ncid)) ! Put file in define mode

   call handle_err(nf90_def_var(ncid = ncid, name = 'deltarho', &
        xtype = nf90_float, dimids = (/ xdimid, ydimid, zdimid/), &
        varid = deltarhoid ) ) ! Get varid

   ! Add attributes
   call handle_err(nf90_put_att(ncid, deltarhoid, 'units', 'kg/m^5'))
   call handle_err(nf90_put_att(ncid, deltarhoid, 'long_name', 'Horizontal Laplacian of density at z=475 m'))
   call handle_err(nf90_enddef(ncid = ncid)) ! End define mode
 
   ! Write deltarho
   call handle_err(nf90_put_var(ncid,deltarhoid,deltarho))
   print *, "deltarho written"

!!!!!!!!!!!

   ! Write deltarho_ft to nc file
   call handle_err(nf90_redef(ncid = ncid)) ! Put file in define mode

   call handle_err(nf90_def_var(ncid = ncid, name = 'deltarho_ft', &
        xtype = nf90_float, dimids = (/ xdimid, ydimid, zdimid/), &
        varid = deltarho_ftid ) ) ! Get varid

   ! Add attributes
   call handle_err(nf90_put_att(ncid, deltarho_ftid, 'units', 'kg/m^5'))
   call handle_err(nf90_put_att(ncid, deltarho_ftid, 'long_name', 'Double transformed horizontal Laplacian of density at z=475 m'))
   call handle_err(nf90_enddef(ncid = ncid)) ! End define mode
 
   ! Write deltarho_ft
   call handle_err(nf90_put_var(ncid,deltarho_ftid,deltarho_ft))
   print *, "deltarho_ft written"

!!!!!!!!!!!!

   ! Write deltarho_alt to nc file
   call handle_err(nf90_redef(ncid = ncid)) ! Put file in define mode

   call handle_err(nf90_def_var(ncid = ncid, name = 'deltarho_alt', &
        xtype = nf90_float, dimids = (/ xdimid, ydimid, zdimid/), &
        varid = deltarho_altid ) ) ! Get varid

   ! Add attributes
   call handle_err(nf90_put_att(ncid, deltarho_altid, 'units', 'kg/m^5'))
   call handle_err(nf90_put_att(ncid, deltarho_altid, 'long_name', 'Recomputed from beta Horizontal Laplacian of density at z=475 m'))
   call handle_err(nf90_enddef(ncid = ncid)) ! End define mode
 
   ! Write deltarho_alt
   call handle_err(nf90_put_var(ncid,deltarho_altid,deltarho_alt))
   print *, "deltarho_alt written"

   call handle_err(nf90_close(ncid = ncid))      

end program poisson_solve_test
