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
    beta_source = gr*laplacian2d(rho,dx,dy)  ! sss
    call solve_poisson(beta_source,x,y,z,beta,'s','d')  ! sss, Dirichlet

    ! Compute rhobar, rhobar_i
    tmp = s2i(3,z,rho) ! move rho to ssi. OK to have rho(z=0)=0 by dp/dz=-rho*g

    do k=1,nz
       rhobar(k) = 1/(real(nx)*real(ny))*sum(rho(:,:,k))
       rhobar_i(k) = 1/(real(nx)*real(ny))*sum(tmp(:,:,k))
    end do
    
    ! Compute pdyn_source_tensor on ssi positions
    pdyn_source_tensor = 0.
    pdyn_source = 0.

    ! (a,b) = (1,1)
    tmpu = s2i(3,z,u) ! u on isi
    tmp = partialder_i2s(1,x,tmpu) ! ppx u on ssi
    do k=1,nz
       tmp(:,:,k) = rhobar_i(k)*tmp(:,:,k)**2 ! rhobar*(ppx u)^2 on ssi  
    end do
    pdyn_source_tensor(1,1,:,:,:) = tmp

    ! (a,b) = (2,2)
    tmpv = s2i(3,z,v) ! v on sii
    tmp = partialder_i2s(2,y,tmpv) ! ppv v on ssi
    do k=1,nz
       tmp(:,:,k) = rhobar_i(k)*tmp(:,:,k)**2 ! rhobar*(ppy v)^2 on ssi  
    end do
    pdyn_source_tensor(2,2,:,:,:) = tmp

    ! (a,b) = (3,3)
    tmpw = i2s(3,z,w) ! w on sss
    tmp = partialder_s2i(3,z,tmpw) ! ppz w on ssi
    do k=1,nz
       tmp(:,:,k) = rhobar_i(k)*tmp(:,:,k)**2 ! rhobar*(ppz w)^2 on ssi  
    end do
    pdyn_source_tensor(3,3,:,:,:) = tmp

    ! (a,b) = (1,2)=(2,1)    
    tmpu = s2i(2,y,i2s(1,x,s2i(3,z,u))) ! Move u to sii
    tmpv = s2i(1,x,i2s(2,y,s2i(3,z,v))) ! Move v to isi 
    tmpu = partialder_i2s(2,y,tmpu) ! ppy u on ssi
    tmpv = partialder_i2s(1,x,tmpv) ! ppx v on ssi
    do k=1,nz
       tmp(:,:,k) = rhobar_i(k)*tmpu(:,:,k)*tmpv(:,:,k) ! rhobar*ppy u*ppx v, ssi
    end do
    pdyn_source_tensor(1,2,:,:,:) = tmp 
    pdyn_source_tensor(2,1,:,:,:) = pdyn_source_tensor(1,2,:,:,:)

    ! (a,b) = (1,3)=(3,1)
    tmpu = i2s(1,x,u) ! move u to sss
    tmpw = s2i(1,x,w) ! move w to isi
    tmpu = partialder_s2i(3,z,u) ! ppz u on ssi
    tmpw = partialder_i2s(1,x,w) ! ppx w on ssi
    do k=1,nz
       tmp(:,:,k) = rhobar_i(k)*tmpu(:,:,k)*tmpw(:,:,k) ! rhobar*ppz u *ppx w 
    end do
    pdyn_source_tensor(1,3,:,:,:) = tmp
    pdyn_source_tensor(3,1,:,:,:) = pdyn_source_tensor(1,3,:,:,:)

    ! (a,b) = (2,3=(3,2)
    tmpv = i2s(2,y,v) ! move v to sss
    tmpw = s2i(2,y,w) ! move w to sii
    tmpv = partialder_s2i(3,z,v) ! ppz v on ssi
    tmpw = partialder_i2s(2,y,w) ! ppy w on ssi
    do k=1,nz
       tmp(:,:,k) = rhobar_i(k)*tmpv(:,:,k)*tmpw(:,:,k) ! rhobar*ppz v *ppy w 
    end do
    pdyn_source_tensor(2,3,:,:,:) = tmp
    pdyn_source_tensor(3,2,:,:,:) = pdyn_source_tensor(2,3,:,:,:)
   
    ! Sum on a,b to get pdyn
    do a=1,3
       do b=1,3
          pdyn_source = pdyn_source + pdyn_source_tensor(a,b,:,:,:) 
       end do
    end do
    
    ! Compute pdyn
    call solve_poisson(pdyn_source,x,y,z, pdyn,'i','n') ! ssi, Neumann
          
  end subroutine compute_pressure


  !===========================!
  ! Subroutine solve_poisson  !
  !                           !
  ! Solves poisson equation   !
  ! -nabla^2 beta = g         !
  ! for source field g        !
  !===========================!

  subroutine solve_poisson(g,x,y,z,beta,s_or_i,d_or_n)

      use, intrinsic :: iso_c_binding
      include '/home1/02291/jeevanje/domain/fftw-3.3.3/install/include/fftw3.f03'
      ! In
      real, dimension(:,:,:), intent(in) :: g ! source field 
      real, dimension(:), allocatable, intent(in) :: x,y,z 
      character(1) :: s_or_i   ! scalar or interface position of z-coord
      character(1) :: d_or_n   ! Dirichlet or Neumann BCs

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
      call ddz2matrix(z,A,s_or_i,d_or_n)  
            
      ! Invert tridiagonal
      do k= 1, nx/2 +1
         do l = 1,ny
            do r = 1, nz
               tmpb(k,l,r) = (2./(dx)**2*(cos(2.*pi*(k-1)/nx)-1) + 2./(dy)**2*(cos(2.*pi*(l-1)/ny)-1)) + A(r,r)
               tmpa(k,l,r) = A(r,r-1)
               tmpc(k,l,r) = A(r,r+1)
               ! These are matrix elements for +nabla^2
             end do
         end do
      end do
      tmprhs = ghat
      ! Use - of tmp matrix elements since inverting -nabla^2
      call invert_tridiagonal(nx/2+1,ny,nz,-tmpa,-tmpb,-tmpc,tmprhs)
      betahat = tmprhs

      ! Backwards FT
      do r = 1,nz
         out = betahat(:,:,r)
         call fftw_execute_dft_c2r(planb,out,in)
         beta(:,:,r) = in/(nx*ny)
      end do

   end subroutine solve_poisson


end module buoyancy_mod
