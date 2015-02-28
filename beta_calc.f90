!======================================!
!  Program beta_calc                   !
!                                      !
!  Computes effective buoyancy of 3D   !
!  or 4D rho field and writes to file  !
!                                      !
!  Input: .nc filepath,                !
!  Output: beta and appended to        !
!          original netcdf file.       !
!                                      !
!  Copyright (c) Nadir Jeevanjee 2015  !
!  http://www.romps.org                !
!======================================!

program beta_calc

   use buoyancy_mod
   use netcdf
   use netcdf_mod
   use, intrinsic :: iso_c_binding

   implicit none

   include '/home1/02291/jeevanje/domain/fftw-3.3.3/install/include/fftw3.f03'

   ! In
   character(200) :: filepath

   ! Privates
   integer :: nx, ny, nz, nt, r, l, file_dim
   integer :: xdimid, ydimid, zdimid, tdimid, ncid, betaid, status
   integer, dimension(:), allocatable :: dimids

   real(c_double) :: dx, dy
   real(c_double), dimension(:), allocatable :: x,y,z
   real(c_double), dimension(:,:,:), allocatable :: rho, beta, deltarho
   real(c_double), dimension(:,:,:,:), allocatable :: tmp

   logical :: file_exist

   ! For FFT
   type(c_ptr) :: planf, planb
   real(c_double), dimension(:,:), allocatable :: in
   complex(c_double_complex), dimension(:,:), allocatable :: out 

   ! Read in filepath 
   call getarg(1,filepath)

   ! Check that file exists
   inquire(file=trim(filepath),exist=file_exist)
         if (.not.file_exist) then
            print *, 'Error in beta_calc: Can not find '//trim(filepath)
          stop
         end if
         if (file_exist) then
            print *, 'file '//trim(filepath)//' found.'
         end if

  ! Open file 
   call handle_err(nf90_open(path = filepath, mode = nf90_write, ncid = ncid))

   ! Retrieve spatial dimensions
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'x', dimid = xdimid))
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'y', dimid = ydimid))
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'z', dimid = zdimid))
   call handle_err(nf90_inquire_dimension(ncid = ncid, &
                   dimid = xdimid, len = nx))
   call handle_err(nf90_inquire_dimension(ncid = ncid, &
                   dimid = ydimid, len = ny))
   call handle_err(nf90_inquire_dimension(ncid = ncid, &
                   dimid = zdimid, len = nz))

   ! Handle time, determine file_dim
   status = nf90_inq_dimid(ncid = ncid, name = 'time', dimid = tdimid)
   if (status == nf90_noerr) then
      file_dim = 4
      allocate( dimids(file_dim), tmp(nx,ny,nz,1) )
      call handle_err(nf90_inquire_dimension(ncid = ncid, &
                       dimid = tdimid, len = nt)) 
      dimids = (/ xdimid, ydimid, zdimid, tdimid /)
   else
      file_dim = 3
      allocate( dimids(file_dim) )
      dimids = (/ xdimid, ydimid, zdimid /)
      nt = 0
   end if

   write(*,*) 'nx,ny,nz,nt =', nx, ny, nz, nt

   allocate(rho(nx,ny,nz), beta(nx,ny,nz), deltarho(nx,ny,nz), &
            x(nx), y(ny), z(nz), in(nx,ny), out(nx/2 +1,ny) )

   print *, 'Arrays allocated'

   ! Get dimension data
   x = get_netCDF1(trim(filepath),'x')
   y = get_netCDF1(trim(filepath),'y')
   z = get_netCDF1(trim(filepath),'z')
   dx = x(2)-x(1)
   dy = y(2)-y(1)

   ! Check to make sure beta isn't already written. Exit if not.
   if ( nf90_inq_varid(ncid = ncid, name = 'rhobeta', varid = betaid) &
        .eq. nf90_noerr) then 
      write(*,*) 'rhobeta already written. Exiting.'
      stop
   end if

   ! Preliminaries for writing output
   ! Put open file in define mode
   call handle_err(nf90_redef(ncid = ncid)) 

   ! Get variable ids
!  call handle_err(nf90_def_var(ncid = ncid, name = 's_rhobeta', &
!         xtype = nf90_float, dimids = dimids, varid = s_betaid ) ) 
   call handle_err(nf90_def_var(ncid = ncid, name = 'rhobeta', &
         xtype = nf90_float, dimids = dimids, varid = betaid ) ) 

   ! Add attributes
!  call handle_err(nf90_put_att(ncid, s_betaid, 'units', 'kg/(m^4 s^2)') ) 
!  call handle_err(nf90_put_att(ncid, s_betaid, &
!        'long_name', 'Source of Effective Buoyancy ('//C_pos//')'))
   call handle_err(nf90_put_att(ncid, betaid, 'units', 'kg/(m^2 s^2)') )
   call handle_err(nf90_put_att(ncid, betaid, &
        'long_name',  'Effective Buoyancy (sss)' ) )

   call handle_err(nf90_enddef(ncid = ncid)) ! End define mode

   ! Begin time loop for calculating an writing variables
   if (file_dim == 3) then
      nt = 1    ! For loop
   end if 

   do l=1,nt
      write(*,*) 'Computing rhobeta for slice', l
      if (file_dim == 4 ) then
           tmp = get_netCDF4(trim(filepath),'rho', &
                 start=(/ 1,1,1,l /), count=(/ nx,ny,nz,1 /) ) 
           rho = tmp(:,:,:,1)
      else
           rho = get_netCDF3(trim(filepath),'rho') 
      end if

      ! Compute delta rho
      deltarho = laplacian2d(rho,dx,dy)
      print *, "Laplacian of rho computed"

      ! Plans for FFTW
      planf = fftw_plan_dft_r2c_2d(ny,nx,in,out,fftw_measure) 
      planb = fftw_plan_dft_c2r_2d(ny,nx,out,in,fftw_measure) 

      ! Solve poisson eqn with gr*deltarho as source
      call solve_poisson(gr*deltarho,x,y,z,beta,'s','d')
      print *, "Poisson equation solved"
   
      ! Write output
      if (file_dim == 4 ) then
         call handle_err(nf90_put_var(ncid,betaid,beta, &
                 start=(/ 1,1,1,l /), count=(/ nx,ny,nz,1 /) ))
      else if  (file_dim == 3) then
         call handle_err(nf90_put_var(ncid,betaid,beta, &
                 start=(/ 1,1,1 /), count=(/ nx,ny,nz/) ))
      end if
   end do

   write(*,*) 'Calculation complete. Closing ncdf file.'

   ! Close nc file
   call handle_err(nf90_close(ncid = ncid))      

end program beta_calc
