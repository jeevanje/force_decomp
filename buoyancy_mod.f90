! Module buoyancy_mod
! Routines for calculating effective buoyancy and dynamic pressure

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

  subroutine compute_pressure(x,y,z,rho,u,v,w,s_beta, & 
        s_eh, s_omega3, s_dh3, s_rhobar, s_dyn, beta, pdyn )

    use, intrinsic :: iso_c_binding

    ! In
    real, dimension(:), allocatable, intent(in) :: x, y, z
    real, dimension(:,:,:), allocatable, intent(in) :: rho, u, v, w 


    ! Private
    integer :: nx, ny, nz, i, j, k
    real :: dx, dy
    real, dimension(:), allocatable :: rhobar, tmp1d
    real, dimension(:,:), allocatable :: A
    real, dimension(:,:,:), allocatable :: deltarho, tmp, tmpu, tmpv,tmpw
    real, dimension(:,:,:,:,:), allocatable :: D ! deformation tensor

    ! Out
    real, dimension(:,:,:), allocatable, intent(out) :: s_beta, &
        s_eh, s_omega3, s_dh3, s_rhobar, s_dyn, beta, pdyn 

    ! Get dimension size and grid spacing
    dx = x(2) - x(1)
    dy = y(2) - y(1)
    nx = size(x,dim=1)
    ny = size(y,dim=1)
    nz = size(z,dim=1)

    ! Allocate arrays
    allocate( rhobar(nz),  deltarho(nx,ny,nz),  &
         beta(nx,ny,nz), pdyn(nx,ny,nz), &
         s_beta(nx,ny,nz), s_dyn(nx,ny,nz), &
         s_eh(nx,ny,nz), s_omega3(nx,ny,nz), &
         s_dh3(nx,ny,nz), s_rhobar(nx,ny,nz), &
         D(3,3,nx,ny,nz), &
         tmpu(nx,ny,nz), tmpv(nx,ny,nz), tmpw(nx,ny,nz), tmp(nx,ny,nz), &
         tmp1d(nz) )

    ! Compute beta, s_beta
    s_beta = gr*laplacian2d(rho,dx,dy) ! sss
    call solve_poisson(s_beta,x,y,z,beta,'s','d')  ! sss, Dirichlet
    
    write(*,*) 'beta, s_beta computed'

    ! Compute rhobar
    do k=1,nz
       rhobar(k) = 1/(real(nx)*real(ny))*sum(rho(:,:,k))
    end do
    
    ! Compute deformation tensor D on sss positions
    D(1,1,:,:,:) =  partialder_i2s(1,x,u)                   ! ppx u on sss
    tmpv         = s2i(1,x,i2s(2,y,v))                      ! Move v to iss 
    D(1,2,:,:,:) = partialder_i2s(1,x,tmpv)                 ! ppx v on sss
    tmpw         = s2i(1,x,i2s(3,z,w,'ns'))                 ! Move w to iss
    D(1,3,:,:,:) = partialder_i2s(1,x,tmpw)                 ! ppx w on sss
    tmpu         = s2i(2,y,i2s(1,x,u))                      ! Move u to sis
    D(2,1,:,:,:) = partialder_i2s(2,y,tmpu)                 ! ppy u on sss
    D(2,2,:,:,:) = partialder_i2s(2,y,v)                    ! ppv v on sss
    tmpw         = s2i(2,y,i2s(3,z,w,'ns'))                 ! Move w to sis
    D(2,3,:,:,:) = partialder_i2s(2,y,tmpw)                 ! ppy w on sss
    tmpu         = i2s(1,x,s2i(3,z,u,'fs'))                 ! Move u to ssi
    D(3,1,:,:,:) = partialder_s2i(3,z,tmpu,'fs')            ! ppz u on sss
    tmpv         = i2s(2,y,s2i(3,z,v,'fs'))                 ! Move v to ssi
    D(3,2,:,:,:) = partialder_s2i(3,z,tmpv,'fs')            ! ppz v on sss
    D(3,3,:,:,:) = partialder_i2s(3,z,w,'ns')               ! ppz w on sss

    ! Compute components of sdyn (ssi)
    call ddz2matrix(z,A,'s','fs')  ! for s_pdyn_rho calc
    tmp1d=matmul(A,log(rhobar))
    tmp1d(1) = 2*tmp1d(2)-tmp1d(3) ! linearly interpolate to bottom level
    do k=1,nz
       s_eh(:,:,k)     = rhobar(k)*( D(1,1,:,:,k)**2 &
                             + 1./2.*( D(1,2,:,:,k)+D(2,1,:,:,k) )**2 &
                             + D(2,2,:,:,k)**2 + D(3,3,:,:,k)**2 )
       s_omega3(:,:,k) = -1/2.*rhobar(k)*( D(1,2,:,:,k)-D(2,1,:,:,k) )**2
       s_dh3(:,:,k)    = 2.*rhobar(k)*( D(1,3,:,:,k)*D(3,1,:,:,k) &
                             + D(2,3,:,:,k)*D(3,2,:,:,k) )
       s_rhobar(:,:,k) = -rhobar(k)*tmp1d(k)*w(:,:,k)**2
    end do
   
    ! Sum to get pdyn
    s_dyn = s_eh + s_omega3 + s_dh3 +s_rhobar
    
    ! Compute pdyn
    call solve_poisson(s_dyn,x,y,z,pdyn,'s','d') ! sss, Dirichlet

    write(*,*) 'pdyn, pdyn_source computed'
          
  end subroutine compute_pressure


  !===========================!
  ! Subroutine solve_poisson  !
  !                           !
  ! Solves poisson equation   !
  ! -nabla^2 beta = g         !
  ! for source field g.       !
  ! Requires specification of !
  ! BCs and scalar or         !
  ! interface levels.A that   !
  !===========================!

  subroutine solve_poisson(g,x,y,z,beta,s_or_i,d_or_n)

      use, intrinsic :: iso_c_binding
      include '/home1/02291/jeevanje/domain/fftw-3.3.3/install/include/fftw3.f03'
      ! In
      real, dimension(:,:,:), intent(in) :: g ! source field 
      real, dimension(:), allocatable, intent(in) :: x,y,z 
      character*(*), intent(in) :: s_or_i, d_or_n

      ! Private
      integer :: nx, ny, nz, i,j,k,l,r
      real :: dx, dy
      real, dimension(:,:), allocatable :: A ! ddz2 matrix

      complex(c_double_complex), dimension(:,:,:), allocatable :: ghat, betahat ! 2D FT of g
      complex(c_double_complex), dimension(:,:,:), allocatable :: tmpa, tmpb, tmpc, tmprhs ! For invert_tridiagonal

      ! For FFT
      real(c_double), dimension(:,:), allocatable :: in
      complex(c_double_complex), dimension(:,:), allocatable :: out 
      type(c_ptr) :: planf, planb

      ! Out
      real, dimension(:,:,:), intent(out) :: beta ! Solution Poisson eqn

      ! Get array size
      nx = size(g,dim=1)
      ny = size(g,dim=2)
      nz = size(g,dim=3) 
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
            
      ! Construct and invert tridiagonal matrix for +nabla^2
      call ddz2matrix(z,A,s_or_i,d_or_n)
      do k= 1, nx/2 +1
         do l = 1,ny
            do r = 1, nz
               tmpb(k,l,r) = (2./(dx)**2*(cos(2.*pi*(k-1)/nx)-1) + 2./(dy)**2*(cos(2.*pi*(l-1)/ny)-1)) + A(r,r)
               if (r .ne. 1) then
                  tmpa(k,l,r) = A(r,r-1)
               end if
               if (r .ne. nz) then
                  tmpc(k,l,r) = A(r,r+1)
               end if
             end do
         end do
      end do
      ! Set beta(kx=0,ky=0)=0 at top for N,i calculation to remove degeneracy
      if ( (d_or_n .eq. 'n').and.(s_or_i .eq. 'i') ) then
         tmpb(1,1,nz) = tmpb(1,1,nz)+2.*A(nz,nz) 
         tmpa(1,1,nz) = tmpa(1,1,nz)+1./2.*A(nz,nz-1)
      end if
      ! Alternate to above: Impose D BC on beta(kx=0,ky=0) for N, i calculation
      !if ( (d_or_n .eq. 'n').and.(s_or_i .eq. 'i') ) then
      !   call ddz2matrix(z,A,'i','d')
      !   do r = 1, nz
      !      tmpb(1,1,r) = A(r,r)
      !      if (r .ne. 1) then
      !         tmpa(1,1,r) = A(r,r-1)
      !      end if
      !      if (r .ne. nz) then
      !         tmpc(1,1,r) = A(r,r+1)
      !      end if
      !   end do
      !end if
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
