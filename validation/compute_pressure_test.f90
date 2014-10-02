!======================================!
!  Program compute_pressure_test       !
!                                      !
! Code to calculate beta and pdyn      !
! and sources via                       !
! subroutine compute_pressure, for     !
! single 3D slice from netcdf file.    !
!                                      !
!  Input: .nc filepath,                !
!     t index for 3-D slab             ! 
!  Output:                             !
!  beta and pdyn calculated and        !
!  appended to original netcdf file.   !
!                                      !
!  Copyright (c) Nadir Jeevanjee       !
!  http://www.romps.org                !
!======================================!

program compute_pressure_test

   use calculus_mod 
   use buoyancy_mod
   use netcdf
   use netcdf_mod
   use, intrinsic :: iso_c_binding

   implicit none

   ! In
   character(200) :: filepath, targ 

   ! Privates
   integer :: nx, ny, nz, nt, tind, ndim
   integer :: xdimid, ydimid, zdimid, tdimid,  ncid, status
   real(c_double), dimension(:,:,:), allocatable :: u,v,w, rho, tmp 
   real(c_double), dimension(:,:,:), allocatable :: s_beta, &
        s_divh, s_eh, s_omega3, s_dh3, s_rhobar, s_dyn, beta, pdyn, Fdyn, & 
       s_dyn_alt
   real(c_double), dimension(:,:,:,:), allocatable :: u4d, v4d, w4d, rho4d
   real(c_double), dimension(:), allocatable :: x,y,z
   real(c_double) :: dx, dy
   logical :: file_exist
   character :: C_pos = 'sss'

   ! Read in filepath and t index
   call getarg(1,filepath)
   call getarg(2,targ,status=status)
   if (status == -1) then
      print *, 'Error: No targ given. Stop.'
      stop
   end if
   read(targ,*), tind 

   ! Check that file exists
   inquire(file=trim(filepath),exist=file_exist)
         if (.not.file_exist) then
            print *, 'Error in compute_pressure_test: Can not find '//trim(filepath)
          stop
         end if
         if (file_exist) then
            print *, 'file '//trim(filepath)//' found.'
         end if

  ! Open file 
   call handle_err(nf90_open(path = filepath, mode = nf90_write, ncid = ncid))

   ! Retrieve dimensions, necessary to allocate arrays
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'x', dimid = xdimid))
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'y', dimid = ydimid))
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'z', dimid = zdimid))
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'time', dimid = tdimid))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = xdimid, len = nx))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = ydimid, len = ny))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = zdimid, len = nz))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = tdimid, len = nt))
   write(*,*) 'nx,ny,nz,nt =', nx, ny, nz, nt

   allocate(rho4d(nx,ny,nz,nt), u4d(nx,ny,nz,nt), &
        v4d(nx,ny,nz,nt), w4d(nx,ny,nz,nt), &
        rho(nx,ny,nz), u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), &
        x(nx), y(ny), z(nz), tmp(nx,ny,nz), s_dyn_alt(nx,ny,nz), &
        Fdyn(nx,ny,nz) )

   print *, 'Arrays allocated'

   ! Get input data
   x = get_netCDF1(trim(filepath),'x')
   y = get_netCDF1(trim(filepath),'y')
   z = get_netCDF1(trim(filepath),'z')
   dx = x(2)-x(1)
   dy = y(2)-y(1)
   rho4d = get_netCDF4(trim(filepath),'rho') 
   u4d   = get_netCDF4(trim(filepath),'u')  
   v4d   = get_netCDF4(trim(filepath),'v')  
   w4d   = get_netCDF4(trim(filepath),'w')  
   
   rho = rho4d(:,:,:,tind)
   u   = u4d(:,:,:,tind)
   v   = v4d(:,:,:,tind)
   w   = w4d(:,:,:,tind)

   ! compute beta, pdyn
  call compute_pressure(x,y,z,rho,u,v,w,s_beta, & 
        s_divh, s_eh, s_omega3, s_dh3, s_rhobar, s_dyn, beta, pdyn )

  ! Compute Fdyn
  Fdyn = -partialder_s2i(3,z,pdyn,'ns')

  ! Check pdyn calculation
  !call laplacian3d(-pdyn,dx,dy,z,'s','d',s_dyn_alt)

   !=================!
   ! Write output    !
   !=================!

   print *, "Calculation complete, writing output to netCDF file"

   ! Write beta_source to nc file
   call write_netCDF3(ncid, (/ xdimid, ydimid, zdimid/), 's_beta', &
        'kg/(m^4 s^2)', 'Source of Effective Buoyancy (sss)',s_beta)

   ! Write s_dyn to nc file
   call write_netCDF3(ncid, (/ xdimid, ydimid, zdimid/), 's_dyn', &
        'kg/(m^3 s^2)', 'Source for Dynamic Pressure ('//C_pos//')',s_dyn)

   ! Write Fdyn to nc file
   call write_netCDF3(ncid, (/ xdimid, ydimid, zdimid/), 'Fdyn', &
       'kg/(m^2 s^2)', 'Dynamic Pressure Force (ssi)', Fdyn)

   ! Write pdyn_source fields to nc file
   call write_netCDF3(ncid, (/ xdimid, ydimid, zdimid/), 's_divh', &
        'kg/(m^3 s^2)', 'x-y convergence component of s_pdyn', s_divh )

   call write_netCDF3(ncid, (/ xdimid, ydimid, zdimid/), 's_eh', &
        'kg/(m^3 s^2)', 'x-y strain component of s_pdyn', s_eh )

   call write_netCDF3(ncid, (/ xdimid, ydimid, zdimid/), 's_omega3', &
        'kg/(m^3 s^2)', 'vertical vorticity component of s_dyn', s_omega3 )

   call write_netCDF3(ncid, (/ xdimid, ydimid, zdimid/), 's_dh3', &
        'kg/(m^3 s^2)', '2*D_13*D_31+2*D_23*D_32 component of s_dyn', s_dh3 )

   call write_netCDF3(ncid, (/ xdimid, ydimid, zdimid/), 's_rhobar', &
        'kg/(m^3 s^2)', 'rhobar component of s_pdyn', s_rhobar )

   ! Write beta to nc file
   call write_netCDF3(ncid, (/ xdimid, ydimid, zdimid/), 'beta', &
        'kg/(m^2 s^2)', 'Effective Buoyancy (sss)',beta)

   ! Write pdyn to nc file
   tmp=pdyn(:,:,1:nz)   ! ??
   call write_netCDF3(ncid, (/ xdimid, ydimid, zdimid/), 'pdyn', &
        'kg/(m s^2)', 'Dynamic Pressure ('//C_pos//')',tmp)

   call handle_err(nf90_close(ncid = ncid))      

end program compute_pressure_test
