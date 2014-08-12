program partialder_test

   use calculus_mod
   use netcdf
   use netcdf_mod
   use buoyancy_mod, only: pi
   use, intrinsic :: iso_c_binding

   implicit none

   ! In
   character(200) :: filepath

   ! Privates
   integer :: nz, ny, nx, i, j, k
   integer :: xdimid, ydimid, zdimid, tdimid, ncid, status
   integer :: ppxsinxid, ppysinyid, ppzsinzid, sinxid, sinyid, sinzid
   real, dimension(:,:,:), allocatable :: sinfield, ppxsinfield, ppysinfield, ppzsinfield,nabla_l2d_sinfield,nabla_pp_sinfield, tmp
   real, dimension(:), allocatable :: x, y, z, z2
   real :: dx, dy, Lx, Ly, Lz
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
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'x', dimid = xdimid))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = xdimid, len = nx))
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'y', dimid = ydimid))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = ydimid, len = ny))
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'z', dimid = zdimid))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = zdimid, len = nz))
   write(*,*) 'nx = ', nx, 'ny= ', ny, 'nz = ', nz

   ! Allocate arrays
   allocate(x(nx),y(ny),z(nz), z2(nz), sinfield(nx,ny,nz), tmp(nx,ny,nz), &
        ppxsinfield(nx,ny,nz), ppysinfield(nx,ny,nz), ppzsinfield(nx,ny,nz), &
        nabla_l2d_sinfield(nx,ny,nz), nabla_pp_sinfield(nx,ny,nz) )
   print *, 'Arrays allocated'

   ! Get coordinates and other parameters
   x = get_netCDF1(trim(filepath),'x')
   y = get_netCDF1(trim(filepath),'y')
   z = get_netCDF1(trim(filepath),'z')
   dx = x(2) - x(1)
   dy = y(2) - y(1)
   Lx = maxval(x)+dx
   Ly = maxval(y)+dy
   Lz = 2.0e3 ! Wavelength in z (not total length in z!)

   z2(1) = 0
   do k=2,nz
      z2(k) = 1./2.*(z(k)+z(k-1))
   end do

   ! Define sinfield w/ horiz wavelength = L/4
   do i=1,nx 
      do j=1,ny
         do k=1,nz
            sinfield(i,j,k) = sin(2.*pi*x(i)*4./Lx)*sin(2.*pi*y(j)*4./Ly)*sin(2.*pi*z(k)/(2.0e3))
         end do
      end do
   end do

   ! Compute partial derivatives
   ppxsinfield = partialder_s2i(1,x,sinfield)
   ppysinfield = partialder_s2i(2,y,sinfield)
   ppzsinfield = partialder_s2i(3,z,sinfield)

   ! Compute 2D Laplacian in two ways
   ! Via laplacian2d
   nabla_l2d_sinfield = laplacian2d(sinfield,dx,dy)
   ! Via partialder2
   nabla_pp_sinfield = partialder2(1,x,sinfield) + partialder2(2,y,sinfield)

   ! Write field, partial derivatives to nc file
   call write_netCDF3(ncid, (/ xdimid, ydimid, zdimid/), 'sinfield','dimensionless', 'Product of sin functions', sinfield)
   call write_netCDF3(ncid, (/ xdimid, ydimid, zdimid/), 'ppxsinfield','dimensionless', 'ppx of sinfield', ppxsinfield)
   call write_netCDF3(ncid, (/ xdimid, ydimid, zdimid/), 'ppysinfield','dimensionless', 'ppy of sinfield', ppysinfield)
   call write_netCDF3(ncid, (/ xdimid, ydimid, zdimid/), 'ppzsinfield','dimensionless', 'ppz of sinfield', ppzsinfield)

   ! Write laplacians to nc file
   call write_netCDF3(ncid, (/ xdimid, ydimid, zdimid/), 'nabla_l2d_sinfield','dimensionless', '2D Laplacian of sinfield', nabla_l2d_sinfield)
   call write_netCDF3(ncid, (/ xdimid, ydimid, zdimid/), 'nabla_pp_sinfield','dimensionless', 'ddx^2+ddy^2 of sinfield', nabla_pp_sinfield)
   call handle_err(nf90_close(ncid = ncid))      

end program partialder_test
