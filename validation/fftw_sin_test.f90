!======================================!
!  Program fftw_sin_test                   !
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
   character(200) :: filepath

   ! Privates
   integer :: nx, ny
   integer :: xdimid, ydimid, ncid, varid, status, field_altid
   real(c_double), dimension(:,:), allocatable :: field, field_alt, fieldhat 
   real(c_double), dimension(:), allocatable :: x,y
   real(c_double), dimension(:,:), allocatable :: in
   real(c_double) :: dx, dy
   complex(c_double_complex), dimension(:,:), allocatable :: out 
   type(c_ptr) :: planf, planb
   logical :: file_exist


   ! Read in filepath and z, t indices
   call getarg(1,filepath)

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
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = xdimid, len = nx))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = ydimid, len = ny))
   write(*,*) 'nx,ny =', nx, ny


   allocate(field(nx,ny), in(nx,ny), &
         out(nx/2+1,ny),  fieldhat(nx/2+1,ny), & 
         field_alt(nx,ny), x(nx), y(ny) )
   planf = fftw_plan_dft_r2c_2d(ny,nx,in,out,fftw_measure)
   planb = fftw_plan_dft_c2r_2d(ny,nx,out,in,fftw_measure)
   print *, 'Arrays allocated'

   ! Get input data
   x = get_netCDF1(trim(filepath),'x')
   y = get_netCDF1(trim(filepath),'y')
   dx = x(2)-x(1)
   dy = y(2)-y(1)
   field = get_netCDF2(trim(filepath),'field') 

      ! Compute FT of deltarho for check
   in = field(:,:)
   call fftw_execute_dft_r2c(planf,in,out)
   fieldhat(:,:) = out
      
   call fftw_destroy_plan(planf)

   ! Backwards FT
   out = fieldhat(:,:)
   call fftw_execute_dft_c2r(planb,out,in)
   field_alt(:,:) = in/(nx*ny)
      
   call fftw_destroy_plan(planb)


   ! Write field_alt to nc file
   call handle_err(nf90_redef(ncid = ncid)) ! Put file in define mode
   call handle_err(nf90_def_var(ncid = ncid, name = 'field_alt', &
        xtype = nf90_float, dimids = (/ xdimid, ydimid /), &
        varid = field_altid ) ) ! Get varid

   ! Add attributes
   call handle_err(nf90_put_att(ncid, field_altid, 'units', 'None'))
   call handle_err(nf90_put_att(ncid, field_altid, 'long_name', 'Sinusoidal field'))
   call handle_err(nf90_enddef(ncid = ncid)) ! End define mode
 
   ! Write data
   call handle_err(nf90_put_var(ncid,field_altid,field_alt))
   print *, "field_alt written"

   call handle_err(nf90_close(ncid = ncid))      

end program fftw_test
