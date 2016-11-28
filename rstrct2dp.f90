! Restriction function using half-weighting
! Takes solution at one grid level and projects it to next 
! coarser grid level 
!
! The variable "grid" determines what part of the staggered grid to use
! grid=0 is for cell centers
! grid=1 is for horizontal velocity, u 
! grid=2 is for vertical velocity, w
!
! 8/15/2011
!
FUNCTION rstrct2dp(uf,ng,mynx,myny,grid)
  USE btype
  IMPLICIT NONE 
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: uf
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
  INTEGER(I4B), INTENT(IN) :: ng,grid
  REAL(WP), DIMENSION(mynx(ng-1),myny(ng-1)) :: rstrct2dp
  INTEGER(I4B) :: nx2,ny2,nx1,ny1
  REAL(WP) :: h,h2
  nx2=size(uf,1)
  ny2=size(uf,2)
  nx1=size(rstrct2dp,1)
  ny1=size(rstrct2dp,2)
  rstrct2dp=0.0_wp	
  h=1.0_wp/4.0_wp
  h2=1.0_wp/8.0_wp

  if (grid==0) then
     rstrct2dp(2:nx1-1,2:ny1-1)=h*(uf(2:nx2-1:2,2:ny2-1:2)+&
          uf(3:nx2:2,2:ny2-1:2)+uf(2:nx2-1:2,3:ny2:2)+uf(3:nx2:2,3:ny2:2))
     ! Boundary grid points
     rstrct2dp(1,2:ny1-1)=0.5_wp*(uf(1,2:ny2-2:2)+uf(1,3:ny2-1:2))
     rstrct2dp(nx1,2:ny1-1)=0.5_wp*(uf(nx2,2:ny2-2:2)+uf(nx2,3:ny2-1:2))
     rstrct2dp(2:nx1-1,1)=0.5_wp*(uf(2:nx2-2:2,1)+uf(3:nx2-1:2,1))
     rstrct2dp(2:nx1-1,ny1)=0.5_wp*(uf(2:ny2-2:2,ny2)+uf(3:nx1-1:2,ny2))
     rstrct2dp(1,1)=uf(1,1)
     rstrct2dp(1,ny1)=uf(1,ny2)
     rstrct2dp(nx1,ny1)=uf(nx2,ny2)
     rstrct2dp(nx1,1)=uf(nx2,1)
     ! Do neighbor exchange after this to get ghost points on coarser grid
  elseif (grid==1) then
       rstrct2dp(2:nx1-1,2:ny1-1)=h*(uf(3:nx2-1:2,2:ny2-2:2)+&
            uf(3:nx2-1:2,3:ny2-1:2))+h2*(uf(2:nx2-2:2,2:ny2-2:2)+&
            uf(2:nx2-2:2,3:ny2-1:2)+uf(4:nx2:2,2:ny2-2:2)+&
            uf(4:nx2:2,3:ny2-1:2))
       ! Boundary grid points
       rstrct2dp(1,2:ny1-1)=0.5_wp*(uf(1,2:ny2-2:2)+uf(1,3:ny2-1:2))
       rstrct2dp(nx1,2:ny1-1)=0.5_wp*(uf(nx2,2:ny2-2:2)+uf(nx2,3:ny2-1:2))
       rstrct2dp(2:nx1-1,1)=2.0_wp*(h*uf(3:nx2-1:2,1)+&
            h2*(uf(2:nx2-2:2,1)+uf(4:nx2:2,1)))
       rstrct2dp(2:nx1-1,ny1)=2.0_wp*(h*uf(3:nx2-1:2,ny2)+&
            h2*(uf(2:nx2-2:2,ny2)+uf(4:nx2:2,ny2)))
       rstrct2dp(1,1)=uf(1,1)
       rstrct2dp(1,ny1)=uf(1,ny2)
       rstrct2dp(nx1,ny1)=uf(nx2,ny2)
       rstrct2dp(nx1,1)=uf(nx2,1)
       ! Do neighbor exchange after this to get ghost points on coarser grid
  elseif (grid==2) then
       rstrct2dp(2:nx1-1,2:ny1-1)=h*(uf(2:nx2-2:2,3:ny2-1:2)+&
            uf(3:nx2-1:2,3:ny2-1:2))+h2*(uf(2:nx2-2:2,2:ny2-2:2)+&
            uf(3:nx2-1:2,2:ny2-2:2)+uf(2:nx2-2:2,4:ny2:2)+&
            uf(3:nx2-1:2,4:ny2:2))
       ! Boundary grid points
       rstrct2dp(2:nx1-1,1)=0.5_wp*(uf(2:nx2-2:2,1)+uf(3:nx2-1:2,1))
       rstrct2dp(2:nx1-1,ny1)=0.5_wp*(uf(2:nx2-2:2,ny2)+uf(3:nx2-1:2,ny2))
       rstrct2dp(1,2:ny1-1)=2.0_wp*(h*uf(1,3:ny2-1:2)+&
            h2*(uf(1,2:ny2-2:2)+uf(1,4:ny2:2)))
       rstrct2dp(nx1,2:ny1-1)=2.0_wp*(h*uf(nx2,3:ny2-1:2)+&
            h2*(uf(nx2,2:ny2-2:2)+uf(nx2,4:ny2:2)))
       rstrct2dp(1,1)=uf(1,1)
       rstrct2dp(1,ny1)=uf(1,ny2)
       rstrct2dp(nx1,ny1)=uf(nx2,ny2)
       rstrct2dp(nx1,1)=uf(nx2,1)
       ! Do neighbor exchange after this to get ghost points on coarser grid
  endif
END FUNCTION rstrct2dp
