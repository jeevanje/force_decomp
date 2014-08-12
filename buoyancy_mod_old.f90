! Module buoyancy_mod
! Tools for calculating effective buoyancy and dynamic pressure

module buoyancy_mod

! Use statements
  use calculus_mod

  implicit none

  real, parameter :: pi = 3.14159265359
  real, parameter :: gr = 9.81  ! Consistent with dam

  contains

      include 'invert_tridiagonal.f90'

  !==============================!
  ! Subroutine compute_pressure  !
  !                              !
  ! Computes effective buoyancy, !
  ! dynamic pressure, and their  !
  ! sources from 3D              !
  ! density and velocity fields  !
  !==============================!

  subroutine compute_pressure(x,y,z,rho,u,v,w,beta_source,pdyn_source_tensor,pdyn_source,beta,pdyn)

    use, intrinsic :: iso_c_binding

    ! In
    real, dimension(:), allocatable, intent(in) :: x, y, z
    real, dimension(:,:,:), allocatable, intent(in) :: rho, u, v, w 


    ! Private
    integer :: nx, ny, nz, i, j, k, a, b
    real :: dx, dy
    real, dimension(:), allocatable :: rhobar, rhobar_i
    real, dimension(:,:,:), allocatable :: deltarho, tmp, tmpu, tmpv,tmpw

    ! Out
    real, dimension(:,:,:), allocatable, intent(out) :: beta, pdyn, beta_source, pdyn_source
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: pdyn_source_tensor

    ! Get dimension size and grid spacing
    dx = x(2) - x(1)
    dy = y(2) - y(1)
    nx = size(x,dim=1)
    ny = size(y,dim=1)
    nz = size(z,dim=1)

    ! Allocate arrays
    allocate( rhobar(nz), rhobar_i(nz), deltarho(nx,ny,nz),  &
         beta(nx,ny,nz), pdyn(nx,ny,nz), &
         beta_source(nx,ny,nz), pdyn_source(nx,ny,nz), &
         pdyn_source_tensor(3,3,nx,ny,nz), &
         tmp(nx,ny,nz), tmpu(nx,ny,nz), tmpv(nx,ny,nz), tmpw(nx,ny,nz) )

    ! Compute beta, beta_source
    tmp = s2i(3,z,rho)  ! Interpolate rho to ssi
    beta_source = -gr*laplacian2d(tmp,dx,dy)  ! ssi
    call solve_poisson(beta_source,x,y,z,beta,'i')  ! ssi

    ! Compute rhobar, rhobar_i
    tmp = s2i(3,z,rho) ! move rho to ssi
    do k=1,nz
       rhobar(k) = 1/(real(nx)*real(ny))*sum(rho(:,:,k))
       rhobar_i(k) = 1/(real(nx)*real(ny))*sum(tmp(:,:,k))
    end do
    
    ! Compute pdyn_source_tensor on ssi positions
    pdyn_source_tensor = 0.
    pdyn_source = 0.

    ! (a,b) = (1,1)
    tmp = i2s(1,x,u) ! Move u to sss
    do k=1,nz
       tmp(:,:,k) = rhobar(k)*tmp(:,:,k)**2  
    end do
    tmp = partialder2(1,x,tmp,'s') ! ppx ppx rhobar*u^2 on sss
    pdyn_source_tensor(1,1,:,:,:) = partialder_s2i(3,z,tmp) 

    ! (a,b) = (2,2)
    tmp = i2s(2,y,v) ! Move v to sss
    do k=1,nz
       tmp(:,:,k) = rhobar(k)*tmp(:,:,k)**2  
    end do
    tmp = partialder2(2,y,tmp,'s') ! ppy ppy rhobar*v^2 on sss
    pdyn_source_tensor(2,2,:,:,:) = partialder_s2i(3,z,tmp) 

    ! (a,b) = (3,3)
    do k=1,nz
       tmp(:,:,k) = rhobar_i(k)*w(:,:,k)**2  ! ssi
    end do
    tmp = partialder2(3,z,tmp,'i') ! ppz^2 rhobar*w^2 on ssi
    tmp = partialder_i2s(3,z,tmp) ! ppz^3 rhobar*w^2 on sss 
    pdyn_source_tensor(3,3,:,:,:) = s2i(3,z,tmp) ! move to ssi

    ! (a,b) = (1,2)=(2,1)
    tmpu = s2i(2,y,u) ! Move u to iis
    tmpv = s2i(1,x,v) ! Move v to iis 
    do k=1,nz
       tmp(:,:,k) = rhobar(k)*tmpu(:,:,k)*tmpv(:,:,k)  ! rhobar*u*v, iis
    end do
    tmp = partialder_i2s(1,x,tmp) ! ppx rhobar*u*v, on sis
    tmp = partialder_i2s(2,y,tmp) ! ppy ppx rhobar*u*v, sss
    pdyn_source_tensor(1,2,:,:,:) = partialder_s2i(3,z,tmp) ! ppz ppx ppx rhobar*u*v, ssi
    pdyn_source_tensor(2,1,:,:,:) = pdyn_source_tensor(1,2,:,:,:)

    ! (a,b) = (1,3)=(3,1)
    tmpw = s2i(1,x,w) ! move w to isi
    tmpu = s2i(3,z,u) ! move u to isi
    do k=1,nz
       tmp(:,:,k) = rhobar_i(k)*tmpu(:,:,k)*tmpw(:,:,k)  
    end do
    tmp = partialder_i2s(1,x,tmp) ! ppx rhobar*u*w, on ssi
    pdyn_source_tensor(1,3,:,:,:) = partialder2(3,z,tmp,'i') ! ppz ppz ppx rhobar*u*w, ssi
    pdyn_source_tensor(3,1,:,:,:) = pdyn_source_tensor(1,3,:,:,:)

    ! (a,b) = (2,3)=(3,2)
    tmpw = s2i(2,y,w) ! move w to sii
    tmpv = s2i(3,z,v) ! move v to sii
    do k=1,nz
       tmp(:,:,k) = rhobar_i(k)*tmpv(:,:,k)*tmpw(:,:,k)  
    end do
    tmp = partialder_i2s(2,y,tmp) ! ppy rhobar*v*w, on ssi
    pdyn_source_tensor(2,3,:,:,:) = partialder2(3,z,tmp,'i') ! ppz ppz ppy rhobar*v*w, ssi
    pdyn_source_tensor(3,2,:,:,:) = pdyn_source_tensor(2,3,:,:,:)
   
    ! Sum on a,b to get pdyn
    do a=1,3
       do b=1,3
          pdyn_source = pdyn_source + pdyn_source_tensor(a,b,:,:,:) 
       end do
    end do
    
    ! Compute pdyn
    call solve_poisson(pdyn_source,x,y,z, pdyn,'i')
          
  end subroutine compute_pressure


  !===========================!
  ! Subroutine solve_poisson  !
  !                           !
  ! Solves poisson equation   !
  ! for source field g        !
  !===========================!

  subroutine solve_poisson(g,x,y,z,beta,s_or_i)

      use, intrinsic :: iso_c_binding
      include '/home1/02291/jeevanje/domain/fftw-3.3.3/install/include/fftw3.f03'
      ! In
      real, dimension(:,:,:), intent(in) :: g ! source field 
      real, dimension(:), allocatable, intent(in) :: x,y,z 
      character(1) :: s_or_i   ! scalar or interface position of z-coord

      ! Private
      integer :: nx, ny, nz, i,j,k,l,r
      real :: dx, dy
      real, dimension(:,:), allocatable :: A
      complex(c_double_complex), dimension(:,:,:), allocatable :: ghat, betahat ! 2D FT of g
      complex(c_double_complex), dimension(:,:,:), allocatable :: tmpa, tmpb, tmpc, tmprhs ! For invert_tridiagonal

      ! For FFT
      real(c_double), dimension(:,:), allocatable :: in
      complex(c_double_complex), dimension(:,:), allocatable :: out 
      type(c_ptr) :: planf, planb

      ! Out
      real, dimension(:,:,:), intent(out) :: beta ! Solution Poisson eqn

      ! Get array size
      nx = size(x,dim=1)
      ny = size(y,dim=1)
      nz = size(z,dim=1)
      dx = x(2)- x(1)
      dy = y(2) - y(1)

      allocate( in(nx,ny), out(nx/2+1,ny), &
           ghat(nx/2 +1,ny,nz), &
           tmpa(nx/2+1,ny,nz), tmpb(nx/2+1,ny,nz), & 
           tmpc(nx/2+1,ny,nz), tmprhs(nx/2+1,ny,nz), & 
           betahat(nx/2+1,ny,nz) )

      ! Plans for FFTW
      planf = fftw_plan_dft_r2c_2d(ny,nx,in,out,fftw_measure)
      planb = fftw_plan_dft_c2r_2d(ny,nx,out,in,fftw_measure)

      ! Compute FT of source field g
      do r = 1,nz
         in = g(:,:,r)
         call fftw_execute_dft_r2c(planf,in,out)
         ghat(:,:,r) = out
      end do

      ! Get matrix for (d/dz)^2
      call ddz2matrix(z,A,s_or_i)
            
      ! Invert tridiagonal
      do k= 1, nx/2 +1
         do l = 1,ny
            do r = 1, nz
               tmpb(k,l,r) = (2./(dx)**2*(cos(2.*pi*(k-1)/nx)-1) + 2./(dy)**2*(cos(2.*pi*(l-1)/ny)-1)) + A(r,r)
               tmpa(k,l,r) = A(r,r-1)
               tmpc(k,l,r) = A(r,r+1)
             end do
         end do
      end do
      tmprhs = ghat
      call invert_tridiagonal(nx/2+1,ny,nz,tmpa,tmpb,tmpc,tmprhs)
      betahat = tmprhs

      ! Backwards FT
      do r = 1,nz
         out = betahat(:,:,r)
         call fftw_execute_dft_c2r(planb,out,in)
         beta(:,:,r) = in/(nx*ny)
      end do

   end subroutine solve_poisson


end module buoyancy_mod
