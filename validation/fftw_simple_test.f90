!======================================!
!  Program fftw_simple_test                   !
!                                      !
! Code to calculate FFT of 2-D arrays  !
! and check proper operation of FFTW   !
! code                                 !
!                                      !
!  Copyright (c) Nadir Jeevanjee       !
!  http://www.romps.org                !
!======================================!

program fftw_simple_test

   use, intrinsic :: iso_c_binding

   implicit none

   include '/home1/02291/jeevanje/domain/fftw-3.3.3/install/include/fftw3.f03'

   ! In
   character(200) :: filepath

   ! Privates
   integer :: n, nx, ny, i,j 
   real(c_double), dimension(:,:), allocatable :: field, field_alt
   real(c_double), dimension(:), allocatable :: field1d, field_alt1d
   real(c_double), dimension(:,:), allocatable :: in
   real(c_double) :: dx, dy
   complex(c_double_complex), dimension(:,:), allocatable :: fieldhat,out
   complex(c_double_complex), dimension(:), allocatable :: fieldhat1d
   type(c_ptr) :: planf, planb, planf1d, planb1d
   logical :: file_exist

   n = 16
   nx = 6
   ny = 6


   allocate(field(nx,ny), in(nx,ny), out(nx/2+1,ny), &
        fieldhat(nx/2+1,ny), field_alt(nx,ny), &
        field1d(n), fieldhat1d(n/2+1), field_alt1d(n) )

   planf1d = fftw_plan_dft_r2c_1d(n,field1d,fieldhat1d,fftw_measure)
   planb1d = fftw_plan_dft_c2r_1d(n,fieldhat1d,field_alt1d,fftw_measure)
   planf = fftw_plan_dft_r2c_2d(ny,nx,in,out,fftw_measure)
   planb = fftw_plan_dft_c2r_2d(ny,nx,out,in,fftw_measure)
   print *, 'Arrays allocated'

   write(*,*) 'n =', n
   call random_seed
   call random_number(field1d)
   write(*,*) 'field1d = ', field1d

   call fftw_execute_dft_r2c(planf1d,field1d,fieldhat1d)
   call fftw_execute_dft_c2r(planb1d,fieldhat1d,field_alt1d)
   field_alt1d = field_alt1d/n
   write(*,*) 'Doubly-transformed field1d = ', field_alt1d

   !stop

   write(*,*) 'nx,ny =', nx, ny
   call random_number(field)
   do i = 1,nx
      do j = 1,ny
         write(*,*) 'field',i,j,' = ', field(i,j)
      end do
   end do

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

   do i = 1,nx
      do j = 1,ny
         write(*,*) 'Doubly transformed field',i,j,' = ', field_alt(i,j)
      end do
   end do

end program fftw_simple_test
