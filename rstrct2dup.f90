! Restriction function using half-weighting for finite volume quantities 
! at cell centers (p,viscosity,temp,etc.)
! Takes solution at one grid level and projects it to next 
! coarser grid level 
! 
! 7/20/2011
!
FUNCTION rstrct2dup(uf,ng,mynx,myny)
  USE btype
  IMPLICIT NONE 
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: uf
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
  INTEGER(I4B), INTENT(IN) :: ng
  REAL(WP), DIMENSION(mynx(ng-1),myny(ng-1)) :: rstrct2dup
  INTEGER(I4B) :: nx2,ny2,nx1,ny1
  REAL(WP) :: h
  nx2=size(uf,1)
  ny2=size(uf,2)
  nx1=size(rstrct2dp,1)
  ny1=size(rstrct2dp,2)
  rstrct2dp=0.0_wp	
  h=1.0_wp/4.0_wp
  h2=1.0_wp/8.0_wp

  rstrct2dup(2:nx1-1,2:ny1-1)=h*(uf(3:nx2-1:2,2:ny2-2:2)+&
       uf(3:nx2-1:2,3:ny2-1:2))+h2*(uf(2:nx2-2:2,2:ny2-2:2)+&
       uf(2:nx2-2:2,3:ny2-1:2)+uf(4:nx2:2,2:ny2-2:2)+uf(4:nx2:2,3:ny2-1:2))
  
  ! Boundary grid points
  rstrct2dup(1,2:ny1-1)=0.5_wp*(uf(1,2:ny2-2:2)+uf(1,3:ny2-1:2))
  rstrct2dup(nx1,2:ny1-1)=0.5_wp*(uf(nx2,2:ny2-2:2)+uf(nx2,3:ny2-1:2))
  rstrct2dup(2:nx1-1,1)=2.0_wp*(h*uf(3:nx2-1:2,1)+&
       h2*(uf(2:nx2-2:2,1)+uf(4:nx2:2,1)))
  rstrct2dup(2:nx1-1,ny1)=2.0_wp*(h*uf(3:nx2-1:2,ny2)+&
       h2*(uf(2:nx2-2:2,ny2)+uf(4:nx2:2,ny2)))
  rstrct2dup(1,1)=uf(1,1)
  rstrct2dup(1,ny1)=uf(1,ny2)
  rstrct2dup(nx1,ny1)=uf(nx2,ny2)
  rstrct2dup(nx1,1)=uf(nx2,1)
  ! Do neighbor exchange after this to get ghost points on coarser grid
END FUNCTION rstrct2dup
