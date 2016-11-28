! 8/27/12
!
! Smooths a scalar field by averaging over a 5 point stencil around each grid 
! point with a user speficied weighting.  'Weight' is the weighting for the 
! central grid point, and 'weight_nb' is the weighting for the 4 neighbor 
! points.  weight_nb = (1-weight)/4  
! 
SUBROUTINE smooth(phi,weight)
  USE btype
  IMPLICIT NONE 
  REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: phi
  REAL(WP), INTENT(IN) :: weight
  REAL(WP), DIMENSION(:,:), ALLOCATABLE :: phi_s
  INTEGER(I4B) :: nx1,ny1,i,j
  REAL(WP) :: weight_nb
  ! 
  nx1=size(phi,1)
  ny1=size(phi,2)
  !
  allocate(phi_s(nx1,ny1))
  phi_s=0.0_wp
  !
  weight_nb=(1.0_wp-weight)/4.0_wp
  !
  forall (i=2:nx1-1,j=2:ny1-1)
     phi_s(i,j)=weight*phi(i,j)+weight_nb*phi(i-1,j)+weight_nb*phi(i+1,j)+&
          weight_nb*phi(i,j+1)+weight_nb*phi(i,j-1)
  end forall
  phi=phi_s
  deallocate(phi_s)
END SUBROUTINE smooth
  
