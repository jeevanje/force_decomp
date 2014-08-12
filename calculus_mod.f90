! Module calculus_mod
! Tools for calculating finite differences

module calculus_mod

! Use statements

  implicit none

 
  contains

  !================================!
  ! function zint                  !
  !                                !
  ! Computes interface values of z !
  ! from scalar values.            !
  !================================!

  function zint(z)

    implicit none
    
    ! In
    real, dimension(:), allocatable, intent(in) :: z  ! z vector   

    ! Private
    integer :: nz, k  

    ! Out
    real, dimension(:), allocatable :: zint ! interface z vector

    ! Get dimension sizes
    nz = size(z,dim=1)

    allocate( zint(nz) )

    zint(1) = 0
    do k=2,nz
      zint(k) = 1./2.*(z(k)+z(k-1))
    end do

  end function zint

  !================================!
  ! function i2s                   !
  !                                !
  ! Moves pth coordinate position  !
  ! of 3d field from interface (i) !
  ! to scalar (s).                 !
  ! Requires optional argument 'bc'!
  ! if moving along z-axis         !
  !================================!

  function i2s(p,coord,f,bc)

    implicit none
    
    ! In
    integer :: p                                          ! Coordinate rank
    real, dimension(:), allocatable, intent(in) :: coord  ! Coordinate vector   
    real, dimension(:,:,:), allocatable, intent(in) :: f  ! 3d  field
    character,  optional :: bc ! ns (no-slip) or fs (free-slip)

    ! Private
    real, dimension(:,:,:), allocatable :: fhalo ! halo field
    real, dimension(:), allocatable :: zhalo
    integer :: n, nx, ny, nz, i, j, k  

    ! Out
    real, dimension(:,:,:), allocatable :: i2s ! translated field

    ! Get dimension sizes
    nx = size(f,dim=1)
    ny = size(f,dim=2)
    nz = size(f,dim=3)

    allocate( fhalo(1:(nx+1),1:(ny+1),1:(nz+1)), & 
         i2s(nx,ny,nz), zhalo(0:(nz+1)) )

    if (p .eq.1) then   ! convert f from ixx to sxx 

       ! Construct f-halo on x-direction
       fhalo(1:nx,1:ny,1:nz) = f
       fhalo(nx+1,1:ny,1:nz) = f(1,:,:)

       ! Compute i2s
       do i=1,nx
          i2s(i,:,:) = 1./2.*( fhalo(i+1,1:ny,1:nz) + fhalo(i,1:ny,1:nz) )
       end do

    else if (p .eq.2) then   ! convert f from xix to xsx

       ! Construct fhalo in y-direction
       fhalo(1:nx,1:ny,1:nz) = f
       fhalo(1:nx,ny+1,1:nz) = f(:,1,:)

       ! Compute i2s
       do j=1,ny
          i2s(:,j,:) = 1./2.*( fhalo(1:nx,j+1,1:nz) + fhalo(1:nx,j,1:nz) )
       end do
       
    else if (p .eq. 3) then  ! convert f from xxi to xxs 

       ! Initialize  zhalo
       zhalo(1:nz) = coord(:)
       zhalo(0) = -coord(1)
       zhalo(nz+1) = 2.*coord(nz)-coord(nz-1)

       ! Construct fhalo in z-direction. Use BCs.
       fhalo(1:nx,1:ny,1:nz) = f
       if (bc .eq. 'ns') then
          fhalo(1:nx,1:ny,nz+1) = 0
       else if (bc .eq. 'fs') then
          fhalo(1:nx,1:ny,nz+1) = 2.*fhalo(1:nx,1:ny,nz)-fhalo(1:nx,1:ny,nz-1)
       end if
       ! Compute i2s
       do k=1,nz
          i2s(:,:,k) = 1./(zhalo(k+1)-zhalo(k-1))* &
               ( fhalo(1:nx,1:ny,k)*(zhalo(k+1)-zhalo(k)) &
                 + fhalo(1:nx,1:ny,k+1)*(zhalo(k)-zhalo(k-1)) )
       end do
    end if
  end function i2s

  !================================!
  ! function s2i                   !
  !                                !
  ! Moves pth coordinate position  !
  ! of 3d field from scalar (s) to !
  ! interface (i).                 !
  ! Requires optional argument 'bc'!
  ! if moving along z-axis         !
  !================================!

  function s2i(p,coord,f,bc)

    implicit none
    
    ! In
    integer :: p                                          ! Coordinate rank
    real, dimension(:), allocatable, intent(in) :: coord  ! Coordinate vector   
    real, dimension(:,:,:), allocatable, intent(in) :: f  ! 3d  field
    character*(*), intent(in), optional :: bc ! ns (no-slip) or fs (free-slip)

    ! Private
    real, dimension(:,:,:), allocatable :: fhalo ! halo field
    real, dimension(:), allocatable :: zhalo
    integer :: n, nx, ny, nz, i, j, k  

    ! Out
    real, dimension(:,:,:), allocatable :: s2i ! translated field

    ! Get dimension sizes
    nx = size(f,dim=1)
    ny = size(f,dim=2)
    nz = size(f,dim=3)

    allocate( fhalo(0:nx,0:ny,0:nz), & 
         s2i(nx,ny,nz), zhalo(0:(nz+1)) )

    if (p .eq.1) then   ! convert f from sxx to ixx 

       ! Construct f-halo on x-direction
       fhalo(1:nx,1:ny,1:nz) = f
       fhalo(0,1:ny,1:nz) = f(nx,:,:)

       ! Compute s2i
       do i=1,nx
          s2i(i,:,:) = 1./2.*( fhalo(i-1,1:ny,1:nz) + fhalo(i,1:ny,1:nz) )
       end do

    else if (p .eq.2) then   ! convert f from xsx to xix

       ! Construct fhalo in y-direction
       fhalo(1:nx,1:ny,1:nz) = f
       fhalo(1:nx,0,1:nz) = f(:,ny,:)

       ! Compute s2i
       do j=1,ny
          s2i(:,j,:) = 1./2.*( fhalo(1:nx,j-1,1:nz) + fhalo(1:nx,j,1:nz) )
       end do
       
    else if (p .eq. 3) then  ! convert f from xxs to xxi

       ! Initialize  zhalo
       zhalo(1:nz) = coord(:)
       zhalo(0) = -coord(1)
       zhalo(nz+1) = 2.*coord(nz)-coord(nz-1)

       ! Construct fhalo in z-direction. Use BCs.
       fhalo(1:nx,1:ny,1:nz) = f
       if (bc .eq. 'ns') then
          fhalo(1:nx,1:ny,0) = -f(1:nx,1:ny,1)
       else if (bc .eq. 'fs') then
          fhalo(1:nx,1:ny,0) = 2.*f(1:nx,1:ny,1)-f(1:nx,1:ny,2)          
       end if

       ! Compute s2i
       do k=1,nz
          s2i(:,:,k) = 1./2.*(fhalo(1:nx,1:ny,k) + fhalo(1:nx,1:ny,k-1))
       end do
    end if
  end function s2i

  !================================!
  ! Function partialder_s2i        !
  !                                !
  ! Computes partial derivative    !
  ! of field with respect to ith   !
  ! coordinate. Assumes ith        !
  ! coordinate in input lives on   !
  ! scalar position, and on        !
  ! interface position in output.  !
  ! Requires optional argument 'bc'!
  ! if differentiating along       ! 
  ! z-axis                         ! 
  !================================!

  function partialder_s2i(p,coord,f,bc)

    implicit none

    ! In
    integer :: p                                          ! Coordinate rank
    real, dimension(:), allocatable, intent(in) :: coord  ! Coordinate vector   
    real, dimension(:,:,:), allocatable, intent(in) :: f  ! 3d field
    character, intent(in), optional :: bc ! ns (no-slip) or fs (free-slip)

    ! Private
    real, dimension(:,:,:), allocatable :: fhalo ! halo field
    real :: delta                                ! uniform grid spacing
    real, dimension(:), allocatable :: zhalo
    integer :: n, nx, ny, nz, i, j, k  

    ! Out
    real, dimension(:,:,:), allocatable :: partialder_s2i ! differentiated field

    ! Get dimension sizes
    nx = size(f,dim=1)
    ny = size(f,dim=2)
    nz = size(f,dim=3)

    allocate( fhalo(0:(nx+1),0:(ny+1),0:(nz+1)), & 
         partialder_s2i(nx,ny,nz), zhalo(0:(nz+1)) )

    if (p .eq.1) then   ! ppix

       delta = coord(2) - coord(1)

       ! Construct f-halo on x-direction
       fhalo(1:nx,1:ny,1:nz) = f
       fhalo(0,1:ny,1:nz) = f(nx,:,:)
       fhalo(nx+1,1:ny,1:nz) = f(1,:,:)

       ! Compute partialder_s2i
       do i=1,nx
          do j=1,ny
             do k=1,nz
                partialder_s2i(i,j,k) = 1./delta*( fhalo(i,j,k) -fhalo(i-1,j,k) )
             end do
          end do
       end do

    else if (p .eq. 2) then   ! ppiy

       delta = coord(2)-coord(1)

       ! Construct fhalo in y-direction
       fhalo(1:nx,1:ny,1:nz) = f
       fhalo(1:nx,0,1:nz) = f(:,ny,:)
       fhalo(1:nx,ny+1,1:nz) = f(:,1,:)
    
       ! Compute partialder_s2i
       do i=1,nx
          do j=1,ny
             do k=1,nz
                partialder_s2i(i,j,k) = 1./delta*( fhalo(i,j,k) -fhalo(i,j-1,k) )
             end do
          end do
       end do

    else if (p .eq. 3) then   ! ppiz

       ! Initialize  zhalo
       zhalo(1:nz) = coord(:)
       zhalo(0) = -coord(1)
       zhalo(nz+1) = 2.*coord(nz)-coord(nz-1)

       ! Construct fhalo in z-direction. Use BCs.
       fhalo(1:nx,1:ny,1:nz) = f
       if (bc .eq. 'ns') then
          fhalo(1:nx,1:ny,0) = -f(:,:,1)
          fhalo(1:nx,1:ny,nz+1) = -f(:,:,nz)
       else if (bc .eq. 'fs') then
          fhalo(1:nx,1:ny,0) = 2.*f(:,:,1)-f(:,:,2)
          fhalo(1:nx,1:ny,nz+1) = 2.*f(:,:,nz)-f(:,:,nz-1)
       end if
          
       ! Compute partialder_s2i
       do i=1,nx
          do j=1,ny
             do k=1,nz
                partialder_s2i(i,j,k) = 1./(zhalo(k)-zhalo(k-1))*( fhalo(i,j,k)-fhalo(i,j,k-1) )
             end do
          end do
       end do
    end if

    end function partialder_s2i

  !================================!
  ! Function partialder_i2s        !
  !                                !
  ! Computes partial derivative    !
  ! of field with respect to ith   !
  ! coordinate. Assumes ith        !
  ! coordinate in input lives on   !
  ! interface position, and on     !
  ! scalar position in output      !
  ! Requires optional argument 'bc'!
  ! if differentiating along       ! 
  ! z-axis                         ! 
  !================================!

  function partialder_i2s(p,coord,f,bc)

    implicit none

    ! In
    integer :: p                                          ! Coordinate rank
    real, dimension(:), allocatable, intent(in) :: coord  ! Coordinate vector   
    real, dimension(:,:,:), allocatable, intent(in) :: f  ! 3d field
    character, intent(in), optional :: bc ! ns (no-slip) or fs (free-slip)

    ! Private
    real, dimension(:,:,:), allocatable :: fhalo ! halo field
    real :: delta                                ! uniform grid spacing
    real, dimension(:), allocatable :: zhalo
    integer :: n, nx, ny, nz, i, j, k  

    ! Out
    real, dimension(:,:,:), allocatable :: partialder_i2s ! differentiated field

    ! Get dimension sizes
    nx = size(f,dim=1)
    ny = size(f,dim=2)
    nz = size(f,dim=3)

    allocate( fhalo(0:(nx+1),0:(ny+1),0:(nz+1)), & 
         partialder_i2s(nx,ny,nz), zhalo(0:(nz+1)) )

    if (p .eq.1) then   ! ppix

       delta = coord(2) - coord(1)

       ! Construct f-halo on x-direction
       fhalo(1:nx,1:ny,1:nz) = f
       fhalo(0,1:ny,1:nz) = f(nx,:,:)
       fhalo(nx+1,1:ny,1:nz) = f(1,:,:)

       ! Compute partialder_i2s
       do i=1,nx
          do j=1,ny
             do k=1,nz
                partialder_i2s(i,j,k) = 1./delta*( fhalo(i+1,j,k) -fhalo(i,j,k) )
             end do
          end do
       end do

    else if (p .eq. 2) then   ! ppiy

       delta = coord(2)-coord(1)

       ! Construct fhalo in y-direction
       fhalo(1:nx,1:ny,1:nz) = f
       fhalo(1:nx,0,1:nz) = f(:,ny,:)
       fhalo(1:nx,ny+1,1:nz) = f(:,1,:)
    
       ! Compute partialder_i2s
       do i=1,nx
          do j=1,ny
             do k=1,nz
                partialder_i2s(i,j,k) = 1./delta*( fhalo(i,j+1,k) -fhalo(i,j,k) )
             end do
          end do
       end do

    else if (p .eq. 3) then   ! ppiz

       ! Initialize  zhalo
       zhalo(1:nz) = coord(:)
       zhalo(0) = -coord(1)
       zhalo(nz+1) = 2.*coord(nz)-coord(nz-1)

       ! Construct fhalo in z-direction. Use BCs.
       fhalo(1:nx,1:ny,1:nz) = f
       if (bc .eq. 'ns') then 
          fhalo(1:nx,1:ny,nz+1) = 0   
       else if (bc .eq. 'fs') then
          fhalo(1:nx,1:ny,nz+1) = 2.*f(:,:,nz)-f(:,:,nz-1)
       end if
          
       ! Compute partialder_i2s
       do i=1,nx
          do j=1,ny
             do k=1,nz
                partialder_i2s(i,j,k) = 2./( zhalo(k+1) - zhalo(k-1) ) *( fhalo(i,j,k+1)-fhalo(i,j,k) )
             end do
          end do
       end do
    end if

    end function partialder_i2s


  !=================================!
  ! Function partialder2            !
  !                                 !
  ! Computes 2nd partial derivative !
  ! of field with respect to ith    !
  ! coordinate. Works for scalar or !
  ! interface z-positions using     !
  ! specification s_or_i. Assume    !
  ! Dirichlet BCs                   !
  !=================================!

  function partialder2(p,coord,f,s_or_i)

    implicit none

    ! In
    integer :: p                                          ! Coordinate rank
    real, dimension(:), allocatable, intent(in) :: coord  ! Coordinate vector   
    real, dimension(:,:,:), allocatable, intent(in) :: f  ! 3d field
    character(1) :: s_or_i                                ! scalar or interface

    ! Private
    real, dimension(:,:,:), allocatable :: fhalo ! halo field
    real :: delta                                ! uniform grid spacing
    real, dimension(:), allocatable :: zhalo
    real, dimension(:,:), allocatable :: A
    integer :: n, nx, ny, nz, i, j, k  

    ! Out
    real, dimension(:,:,:), allocatable :: partialder2 ! differentiated field

    ! Get dimension sizes
    nx = size(f,dim=1)
    ny = size(f,dim=2)
    nz = size(f,dim=3)

    allocate( fhalo(0:(nx+1),0:(ny+1),0:(nz+1)), & 
         partialder2(nx,ny,nz), zhalo(0:(nz+1)) )

    if (p .eq.1) then   ! ppix

       delta = coord(2) - coord(1)

       ! Construct f-halo on x-direction
       fhalo(1:nx,1:ny,1:nz) = f
       fhalo(0,1:ny,1:nz) = f(nx,:,:)
       fhalo(nx+1,1:ny,1:nz) = f(1,:,:)

       ! Compute partialder2
       do i=1,nx
          do j=1,ny
             do k=1,nz
                partialder2(i,j,k) = 1./(delta)**2*( fhalo(i+1,j,k) - &
                     2.*fhalo(i,j,k) + fhalo(i-1,j,k) )
             end do
          end do
       end do

    else if (p .eq. 2) then   ! ppiy

       delta = coord(2)-coord(1)

       ! Construct fhalo in y-direction
       fhalo(1:nx,1:ny,1:nz) = f
       fhalo(1:nx,0,1:nz) = f(:,ny,:)
       fhalo(1:nx,ny+1,1:nz) = f(:,1,:)
    
       ! Compute partialder2
       do i=1,nx
          do j=1,ny
             do k=1,nz
                partialder2(i,j,k) = 1./(delta)**2*( fhalo(i,j+1,k) - &
                     2.*fhalo(i,j,k) + fhalo(i,j-1,k) )
             end do
          end do
       end do

    else if (p .eq. 3) then   ! ppiz

       if ( s_or_i .eq. 's') then  ! scalar values

          ! Construct fhalo in z-direction with f = 0 on z-boundaries
          fhalo(1:nx,1:ny,1:nz) = f
          fhalo(1:nx,1:ny,0) = -f(:,:,1)
          fhalo(1:nx,1:ny,nz+1) = -f(:,:,nz)
    
          ! Use ddz2 matrix to compute 2nd derivative
          call ddz2matrix(coord,A,s_or_i,'d')
          do i=1,nx
             do j=1,ny
                partialder2(i,j,:)=matmul(A,fhalo(i,j,:))
             end do
          end do

       else  if (s_or_i .eq. 'i') then ! interface values
          
          ! Construct fhalo in z-direction with f = 0 on z-boundaries
          fhalo(1:nx,1:ny,1:nz) = f
          fhalo(1:nx,1:ny,0) = -f(:,:,2)  ! since f(:,:,1) = 0 for interface
          fhalo(1:nx,1:ny,nz+1) = 0       ! top boundary
    
          ! Use ddz2 matrix to compute 2nd derivative
          call ddz2matrix(coord,A,s_or_i,'d')
          do i=1,nx
             do j=1,ny
                partialder2(i,j,:)=matmul(A,fhalo(i,j,:))
             end do
          end do
       else 
         write (*,*) ' Variable s_or_i not set to allowed value. Exiting.'
         stop       
       end if
    
    end if

  end function partialder2

  !================================!
  ! Function laplacian2d           !
  !                                !
  ! Computes horizontal laplacian  !
  ! of  doubly-periodic field.     !
  ! nz is along for the ride,      !
  ! allows us to compute laplacian !
  ! for all levels at once.        !
  !================================!

  function laplacian2d(f,dx,dy)

      implicit none

      ! In
      real, dimension(:,:,:), intent(in) :: f ! field of interest
      real, intent(in) :: dx, dy

      ! Private
      integer :: nx, ny, nz, i, j
      real, dimension(:,:,:), allocatable :: fhalo ! extend f with 1 unit halo

      ! Out
      real, dimension(:,:,:), allocatable :: laplacian2d ! Laplacian of f

      ! Get array sizes
      nx = size(f,dim=1)
      ny = size(f,dim=2)
      nz = size(f,dim=3)

      allocate( fhalo(0:(nx+1),0:(ny+1),1:nz), laplacian2d(nx,ny,nz) )

      ! Construct fhalo
      fhalo(1:nx,1:ny,:) = f
      fhalo(0,1:ny,:) = f(nx,1:ny,:)
      fhalo(nx+1,1:ny,:) = f(1,1:ny,:)
      fhalo(1:nx,0,:) = f(1:nx,ny,:)
      fhalo(1:nx,ny+1,:) = f(1:nx,1,:)
    
      ! Compute laplacian2d
      do i=1,nx
         do j = 1,ny 
            laplacian2d(i,j,:) = 1/(dx)**2*(fhalo(i+1,j,:) -2.*fhalo(i,j,:) &
            + fhalo(i-1,j,:)) + 1/(dy)**2*(fhalo(i,j+1,:) -2.*fhalo(i,j,:) &
            + fhalo(i,j-1,:) ) 
         end do
      end do

   end function laplacian2d

  !============================= !
  ! Subroutine laplacian3d       !
  !                              !
  ! Computes discrete laplacian  !
  ! of 3d doubly-periodic field  !
  ! for scalar or interface      !
  ! levels and for variable BCs. !
  !============================= !

  subroutine laplacian3d(f,dx,dy,z,s_or_i,bc,deltaf)

      implicit none

      ! In
      real, dimension(:,:,:), intent(in) :: f ! field of interest
      real, dimension(:), allocatable, intent(in) :: z ! vertical grid
      real :: dx, dy
      character*(*), intent(in) :: s_or_i, bc

      ! Private
      integer :: nx, ny,nz, i,j,k
      real, dimension(:,:,:), allocatable :: fhalo ! extend f with 1 unit halo
      real, dimension(:), allocatable :: Af ! 
      real, dimension(:,:), allocatable :: A ! 
      real :: dz
      logical :: ziso

      ! Out
      real, dimension(:,:,:), allocatable, intent(out) :: deltaf ! Laplacian of f

      ! Get array size
      nx = size(f,dim=1)
      ny = size(f,dim=2)
      nz = size(f,dim=3)

      ! Determine if z grd is isotropic
      ziso = ( (z(nz)-z(nz-1)) .eq. (z(2)-z(1)) )
      if (ziso) then
         dz = z(2) - z(1)
      end if
      
      allocate( fhalo(0:(nx+1),0:(ny+1),1:nz), &
           Af(1:nz), A(nz,nz),&
           deltaf(nx,ny,nz) )

      ! Construct fhalo
      fhalo(1:nx,1:ny,1:nz) = f
      fhalo(0,1:ny,1:nz) = f(nx,1:ny,1:nz)
      fhalo(nx+1,1:ny,1:nz) = f(1,1:ny,1:nz)
      fhalo(1:nx,0,1:nz) = f(1:nx,ny,1:nz)
      fhalo(1:nx,ny+1,1:nz) = f(1:nx,1,1:nz)

      ! Compute delta f
      call ddz2matrix(z,A,s_or_i,bc)
      do i=1,nx
         do j = 1,ny 
            Af=matmul(A,fhalo(i,j,:))
            do k = 1, nz
               deltaf(i,j,k) = 1/(dx)**2*(fhalo(i+1,j,k) -2*fhalo(i,j,k) &
                    + fhalo(i-1,j,k)) + 1/(dy)**2*( fhalo(i,j+1,k) &
                    -2*fhalo(i,j,k) + fhalo(i,j-1,k) ) + Af(k)
            end do
         end do
      end do


   end subroutine laplacian3d


  !==========================!
  ! Subroutine ddz2matrix    !
  !                          !
  ! For given vertical grid, !
  ! computes matrix that     !
  ! implements (d/dz)^2 on   !
  ! quantity f. Switch bc    !
  ! specifies Dirichlet BCs  !
  ! with f(0)=f(H)=0 or      ! 
  ! Neumann with df/dz(0)=   !
  ! df/dz(H)=0               !
  !                          !
  !==========================!

   subroutine ddz2matrix(z,A,s_or_i,bc)

      implicit none
  
      ! In
      real, dimension(:), allocatable, intent(in) :: z ! vertical grid
      character(1) :: s_or_i  ! scalar or interface
      character*(*) :: bc      ! Boundary conditions (d,n, or fs)
      
      ! Private
      integer :: nz, k, n 
      real, dimension(:), allocatable :: zhalo, z2halo ! Extend z with halo

      ! Out
      real, dimension(:,:), allocatable, intent(out) :: A ! Desired matrix

      write(*,*) 'Calculating ddz2matrix for ',s_or_i, ' levels with ', bc, &
                 ' boundary conditions'

      ! Check boundary conditions
      if ( (bc .ne. 'd').and.(bc .ne. 'n').and.(bc .ne. 'fs')) then
         print *, 'Boundary conditions must be either "d", "n", or "fs". '
      stop
      end if

      ! Check C grid position
      if ( (s_or_i .ne. 's').and.(s_or_i .ne. 'i') ) then
         print *, 'C grid position must be either "s" or "i". '
      stop
      end if
      
      ! Get nz, allocate arrays
      nz = size(z,dim=1)
      allocate( z2halo(0:(nz+1)), zhalo(0:(nz+1)), A(1:nz,1:nz) )
      
      ! Initialize A, zhalo, z2halo
      A=0. 
      zhalo(1:nz) = z
      zhalo(0) = -z(1)
      zhalo(nz+1) = 2*z(nz)-z(nz-1)
      z2halo(1:nz) = zint(z)
      z2halo(0) = - z2halo(2)
      z2halo(nz+1) = 2*z(nz)-z2halo(nz)

      ! Do scalar calculation
      if (s_or_i .eq. 's') then
         do k = 1,nz
            if (k .ne. 1) then
               A(k,k-1) = 2*( (zhalo(k+1)-zhalo(k-1))*(zhalo(k)-zhalo(k-1)) )**(-1)
            end if
            A(k,k) = - 2/(zhalo(k+1)-zhalo(k-1))*((zhalo(k+1)-zhalo(k))**(-1)                         + (zhalo(k)-zhalo(k-1))**(-1))
            if (k .ne. nz) then
               A(k,k+1) = 2*( (zhalo(k+1)-zhalo(k-1))*(zhalo(k+1)-zhalo(k)) )**(-1)
            end if
         end do
         ! Overwrite some values to implement Dirichlet BCs 
         if (trim(bc) .eq. 'd') then
            A(nz,nz) = -2/(zhalo(nz+1)-zhalo(nz-1))*( &
                 2/(zhalo(nz+1)-zhalo(nz)) + 1./(zhalo(nz)-zhalo(nz-1)) )
            A(1,1) = -2/(zhalo(2)-zhalo(0))*( 1./(zhalo(2)-zhalo(1)) &
                 + 2/(zhalo(1)-zhalo(0)) )            
         ! Ditto for Neumann BCs as per my notes, 4/30/14
         else if (trim(bc) .eq. 'n') then
            A(1,1) = -A(1,2)
            A(nz,nz) = -A(nz,nz-1)
         end if

      ! Interface calculation  
      else if (s_or_i .eq. 'i') then
         do k= 1,nz
            if (k .ne. 1) then
               A(k,k-1)=( (zhalo(k)-zhalo(k-1))*(z2halo(k)-z2halo(k-1)) )**(-1)
            end if
            A(k,k) = - 1/( zhalo(k)-zhalo(k-1) )* &
                  ((z2halo(k+1)-z2halo(k))**(-1)+(z2halo(k)-z2halo(k-1))**(-1))
            if (k .ne. nz) then
               A(k,k+1)=( (zhalo(k)-zhalo(k-1))*(z2halo(k+1)-z2halo(k)) )**(-1)
            end if
         end do
         !  BCs
         if (trim(bc) .eq. 'd') then ! Dirichlet, top BCs implmntd automatically
            A(1,1) = 0. 
            A(2,1) = 0. 
         else if (trim(bc) .eq. 'n') then ! Neumann as per notes 
                                          ! 5/2-5/5, 5/15 2014
            A(1,2) = -A(1,1)          
            A(nz,nz) = 1./3.*A(nz,nz)
            A(nz,nz-1) = 2./3.*A(nz,nz-1)
         else if (trim(bc) .eq. 'fs') then
            A(1,1) = 0.   ! Since d2f/dz2 = 0 at bottom
            A(1,2) = 0.
            A(nz,nz-2) = 1./2.*A(nz,nz-1) ! As per notes 5/14/14
            A(nz,nz-1) = 1./2.*A(nz,nz)
            A(nz,nz) = -1./4.*A(nz,nz)
         end if ! bc
      end if ! s_or_i


   end subroutine ddz2matrix



end module calculus_mod
