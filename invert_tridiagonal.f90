!==================================!
! Das Atmosphaerische Modell (DAM) !
! Copyright (c) David M. Romps     !
! http://www.romps.org             !
!==================================!

subroutine invert_tridiagonal(nx,ny,n,tmpa,tmpb,tmpc,tmprhs)

!=======================================================!
! Solves A*x = b, where A is an n x n matrix,           !
!  and x and b are length-nz vectors.                   !
!                                                       !
! The nx and ny indices are just along for the ride.    !
! They allow an entire subdomain's worth of columns     !
! to be calculated at once.  The integer n will         !
! be nz from the domain.                                !
!                                                       !
! tmpa = lower diagonal of A, first element is not used !
! tmpb = main diagonal of A                             !
! tmpc = upper diagonal of A, last element is not used  !
! tmprhs = b                                            !
!                                                       !
! Adapted from code by Peter Blossey                    !
!=======================================================!

implicit none

! In
integer :: nx,ny,n

! In/Out
complex, dimension(1:nx,1:ny,n) :: tmpa, tmpb, tmpc, tmprhs

integer :: k, nm

nm = n-1

! invert tridiagonal system: no pivoting.
do k = 1,nm

   ! forward sweep: normalize by diagonal element.
   ! note that diagonal element is one after this step
   tmprhs(:,:,k) = tmprhs(:,:,k)/tmpb(:,:,k)
   tmpc(:,:,k)   = tmpc(:,:,k)/tmpb(:,:,k)

   ! forward sweep: eliminate lower diagonal element from next eqn.
   tmpb(:,:,k+1)   = tmpb(:,:,k+1)   - tmpa(:,:,k+1)*tmpc(:,:,k)
   tmprhs(:,:,k+1) = tmprhs(:,:,k+1) - tmpa(:,:,k+1)*tmprhs(:,:,k)

end do

tmprhs(:,:,n) = tmprhs(:,:,n)/tmpb(:,:,n)

do k = nm,1,-1
   ! backward sweep
   tmprhs(:,:,k) = tmprhs(:,:,k) - tmpc(:,:,k)*tmprhs(:,:,k+1)
end do



end subroutine invert_tridiagonal
