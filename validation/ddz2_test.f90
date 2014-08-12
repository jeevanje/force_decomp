program ddz2_test

   use calculus_mod
   use buoyancy_mod, only : pi
   use netcdf
   use netcdf_mod
   use, intrinsic :: iso_c_binding

   implicit none

   include '/home1/02291/jeevanje/domain/fftw-3.3.3/install/include/fftw3.f03'

   ! In
   character(200) :: filepath

   ! Privates
   integer :: nz,  zind, k
   integer :: zdimid, tdimid, ncid,  status, ddz2sinzid, sinzid
   real, dimension(:,:), allocatable ::  A
   real, dimension(:), allocatable :: x, y, z, zhalo, z2, sinz, ddz2sinz
   real :: dx, dy
   complex(c_double_complex), dimension(:,:), allocatable :: out 
   logical :: file_exist


   ! Read in filepath 
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
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'z', dimid = zdimid))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = zdimid, len = nz))
   write(*,*) 'nz = ',nz

   allocate(zhalo(0:(nz+1)),z(1:nz), &
        z2(0:(nz+1)), sinz(0:(nz+1)), ddz2sinz(nz) )
   print *, 'Arrays allocated'

   ! Get input data
   z = get_netCDF1(trim(filepath),'z')
   zhalo(1:nz) = z
   zhalo(0) = -z(1)
   zhalo(nz+1) = 2*z(nz)-z(nz-1)
   z2(1:nz) = zint(z)
   z2(nz+1) = 2*z(nz)-z2(nz)
   z2(0) = -z2(2)
   sinz = sin(2*pi*z2/(2.0e3)) ! 2 km wavelength

   ! Compute A
   call ddz2matrix(z,A,'i')
   ddz2sinz = matmul(A,sinz)

   write(*,*) 'ddz2sinz =', ddz2sinz

   ! Write sinz, ddz2sinz to nc file
   call write_netCDF1(ncid, (/ zdimid/), 'ddz2sinz','dimensionless', &
        '(ddz)^2 of sin(2 pi z/2e3)', ddz2sinz)
   call write_netCDF1(ncid, (/ zdimid/), 'sinz','dimensionless', &
        'sin(2 pi z/2e3)', sinz(1:nz))

   call handle_err(nf90_close(ncid = ncid))      

end program ddz2_test
