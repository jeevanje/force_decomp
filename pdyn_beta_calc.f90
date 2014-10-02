!======================================!
!  Program pdyn_beta_calc              !
!                                      !
! Code to calculate beta and pdyn via  !
! subroutine compute_pressure, and     !
! append to  netcdf file.              !
!                                      !
!  Input: .nc filepath,                !
!  Output:                             !
!  beta and pdyn calculated and        !
!  appended to original netcdf file.   !
!                                      !
!  Copyright (c) Nadir Jeevanjee       !
!  http://www.romps.org                !
!======================================!

program pdyn_beta_calc

   use buoyancy_mod
   use calculus_mod
   use netcdf
   use netcdf_mod
   use, intrinsic :: iso_c_binding

   implicit none

   ! In
   character(200) :: filepath

   ! Privates
   integer :: nx, ny, nz, nt, l
   integer ::  betaid, pdynid, s_betaid,s_dynid, Fdynid
   integer, dimension(4) :: dimids
   integer :: xdimid, ydimid, zdimid, tdimid,  ncid, status
   real(c_double), dimension(:,:,:), allocatable :: u3d,v3d,w3d, rho3d, &
        s_beta3d , s_divh, s_eh, s_omega3, s_dh3, s_rhobar, s_dyn3d, beta3d, pdyn3d, &
        Fdyn3d
   real, dimension(:,:,:,:), allocatable :: tmp
   real(c_double), dimension(:), allocatable :: x,y,z, time
   real(c_double) :: dx, dy
   logical :: file_exist
   character(3) :: C_pos='sss'   ! ssi or sss

   ! Read in filepath and t index
   call getarg(1,filepath)

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

   ! Retrieve dimensions, necessary to allocate arrays
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'x', dimid = xdimid))
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'y', dimid = ydimid))
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'z', dimid = zdimid))
   call handle_err(nf90_inq_dimid(ncid = ncid, name = 'time', dimid = tdimid))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = xdimid, len = nx))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = ydimid, len = ny))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = zdimid, len = nz))
   call handle_err(nf90_inquire_dimension(ncid = ncid, dimid = tdimid, len = nt))
   dimids = (/ xdimid, ydimid, zdimid, tdimid /)
   write(*,*) 'nx,ny,nz,nt =', nx, ny, nz, nt

   ! Allocate. Vars not allocated here are allocated in subroutines
   allocate( rho3d(nx,ny,nz), u3d(nx,ny,nz), v3d(nx,ny,nz), w3d(nx,ny,nz), &
        tmp(nx,ny,nz,1), Fdyn3d(nx,ny,nz), &
        x(nx), y(ny), z(nz), time(nt) )

   print *, 'Arrays allocated'

   ! Get input data
   x = get_netCDF1(trim(filepath),'x')
   y = get_netCDF1(trim(filepath),'y')
   z = get_netCDF1(trim(filepath),'z')
   time = get_netCDF1(trim(filepath),'time')
   dx = x(2)-x(1)
   dy = y(2)-y(1)

   ! Preliminaries for writing output
   ! Check to make sure beta isn't already written. Exit if not.
   if ( nf90_inq_varid(ncid = ncid, name = 'beta', varid = betaid) &
        .eq. nf90_noerr) then 
      write(*,*) 'Beta already written. Exiting.'
      stop
   end if

   ! Put open file in define mode
   call handle_err(nf90_redef(ncid = ncid)) 

   ! Get variable ids
   call handle_err(nf90_def_var(ncid = ncid, name = 's_beta', &
            xtype = nf90_float, dimids = dimids, varid = s_betaid ) ) 
   call handle_err(nf90_def_var(ncid = ncid, name = 's_dyn', &
            xtype = nf90_float, dimids = dimids, varid = s_dynid ) ) 
   call handle_err(nf90_def_var(ncid = ncid, name = 'beta', &
            xtype = nf90_float, dimids = dimids, varid = betaid ) ) 
   call handle_err(nf90_def_var(ncid = ncid, name = 'pdyn', &
            xtype = nf90_float, dimids = dimids, varid = pdynid ) ) 
   call handle_err(nf90_def_var(ncid = ncid, name = 'Fdyn', &
            xtype = nf90_float, dimids = dimids, varid = Fdynid ) ) 

   ! Add attributes
   call handle_err(nf90_put_att(ncid, s_betaid, 'units', 'kg/(m^4 s^2)') ) 
   call handle_err(nf90_put_att(ncid, s_betaid, &
        'long_name', 'Source of Effective Buoyancy ('//C_pos//')'))

   call handle_err(nf90_put_att(ncid, s_dynid, 'units', 'kg/(m^3 s^2)'))
   call handle_err(nf90_put_att(ncid, s_dynid, &
        'long_name',  'Source for Dynamic Pressure ('//C_pos//')' ) )

   call handle_err(nf90_put_att(ncid, betaid, 'units', 'kg/(m^2 s^2)') )
   call handle_err(nf90_put_att(ncid, betaid, &
        'long_name',  'Effective Buoyancy ('//C_pos//')') )

   call handle_err(nf90_put_att(ncid, pdynid, 'units','kg/(m s^2)') )
   call handle_err(nf90_put_att(ncid, pdynid, & 
        'long_name',  'Dynamic Pressure ('//C_pos//')'))

   call handle_err(nf90_put_att(ncid, Fdynid, 'units', 'kg/(m^2 s^2)') )
   call handle_err(nf90_put_att(ncid, Fdynid, &
        'long_name',  'Dynamic Pressure Force (ssi)' ))

   call handle_err(nf90_enddef(ncid = ncid)) ! End define mode

   ! time loop for calculating and writing variables   
   do l=1,nt
      write(*,*) 'Computing beta, pdyn, Fdyn for snapshot', l
      tmp = get_netCDF4(trim(filepath),'rho', &
                 start=(/ 1,1,1,l /), count=(/ nx,ny,nz,1 /) ) 
      rho3d = tmp(:,:,:,1)

      tmp   = get_netCDF4(trim(filepath),'u', &
                 start=(/ 1,1,1,l /), count=(/ nx,ny,nz,1 /) ) 
      u3d   = tmp(:,:,:,1)      

      tmp   = get_netCDF4(trim(filepath),'v', &
                 start=(/ 1,1,1,l /), count=(/ nx,ny,nz,1 /) ) 
      v3d   = tmp(:,:,:,1)

      tmp   = get_netCDF4(trim(filepath),'w', &
                 start=(/ 1,1,1,l /), count=(/ nx,ny,nz,1 /) ) 
      w3d   = tmp(:,:,:,1)

      ! Compute beta, pdyn, Fdyn 
      call compute_pressure(x,y,z,rho3d,u3d,v3d,w3d,s_beta3d, &
            s_divh, s_eh, s_omega3, s_dh3, s_rhobar, s_dyn3d, beta3d, pdyn3d )
      Fdyn3d  = - partialder_s2i(3,z,pdyn3d,'ns')

      ! Write output
      call handle_err(nf90_put_var(ncid,s_betaid,s_beta3d, &
                 start=(/ 1,1,1,l /), count=(/ nx,ny,nz,1 /) ))
      call handle_err(nf90_put_var(ncid,s_dynid,s_dyn3d, &
                 start=(/ 1,1,1,l /), count=(/ nx,ny,nz,1 /) ))
      call handle_err(nf90_put_var(ncid,betaid,beta3d, &
                 start=(/ 1,1,1,l /), count=(/ nx,ny,nz,1 /) ))
      call handle_err(nf90_put_var(ncid,pdynid,pdyn3d, &
                 start=(/ 1,1,1,l /), count=(/ nx,ny,nz,1 /) ))
      call handle_err(nf90_put_var(ncid,Fdynid,Fdyn3d, &
                 start=(/ 1,1,1,l /), count=(/ nx,ny,nz,1 /) ))
   end do

   print *, "beta, pdyn, Fdyn computed"
   
   ! Close netcdf file
   call handle_err(nf90_close(ncid = ncid))      

end program pdyn_beta_calc
