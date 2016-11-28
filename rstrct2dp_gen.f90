! New full weighting restriction operator for general grid (allows for uneven
! grid cell sizes) 
!
! The variable "grid" determines what part of the staggered grid to use
! grid=0 is for cell centers
! grid=1 is for horizontal velocity, u 
! grid=2 is for vertical velocity, w
!
! 9/20/16
!
FUNCTION rstrct2dp_gen(uf,ng,mynx,myny,grid,x1,y1,x2,y2)
  USE btype
  IMPLICIT NONE 
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: uf
  REAL(WP), DIMENSION(:), INTENT(IN) :: x1,y1,x2,y2
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
  INTEGER(I4B), INTENT(IN) :: ng,grid
  REAL(WP), DIMENSION(mynx(ng-1),myny(ng-1)) :: rstrct2dp_gen
  INTEGER(I4B) :: nx2,ny2,nx1,ny1,i,j
  REAL(WP) :: h,h2
  nx2=size(uf,1)
  ny2=size(uf,2)
  nx1=size(rstrct2dp_gen,1)
  ny1=size(rstrct2dp_gen,2)
  rstrct2dp_gen=0.0_wp	
  h=1.0_wp/4.0_wp
  h2=1.0_wp/8.0_wp

  if (grid==0) then
     do j=2,ny1-1
        do i=2,nx1-1
           rstrct2dp_gen(i,j)=(1.0_wp/4.0_wp)*(((y2(2*j-3)-y1(j-1))/&
                (y1(j)-y1(j-1)))*((x2(2*i-3)-x1(i-1))/&
                (x1(i)-x1(i-1)))*uf(2*i-3,2*j-3)+&
                ((y2(2*j-3)-y1(j-1))/(y1(j)-y1(j-1)))*&
                ((x2(2*i-2)-x1(i-1))/(x1(i)-x1(i-1)))*uf(2*i-2,2*j-3)+&
                ((y2(2*j-3)-y1(j-1))/(y1(j)-y1(j-1)))*&
                ((x1(i+1)-x2(2*i-1))/(x1(i+1)-x1(i)))*uf(2*i-1,2*j-3)+&
                ((y2(2*j-3)-y1(j-1))/(y1(j)-y1(j-1)))*&
                ((x1(i+1)-x2(2*i))/(x1(i+1)-x1(i)))*uf(2*i,2*j-3)+&
                ((y2(2*j-2)-y1(j-1))/(y1(j)-y1(j-1)))*&
                ((x2(2*i-3)-x1(i-1))/(x1(i)-x1(i-1)))*uf(2*i-3,2*j-2)+&
                ((y2(2*j-2)-y1(j-1))/(y1(j)-y1(j-1)))*&
                ((x2(2*i-2)-x1(i-1))/(x1(i)-x1(i-1)))*uf(2*i-2,2*j-2)+&
                ((y2(2*j-2)-y1(j-1))/(y1(j)-y1(j-1)))*&
                ((x1(i+1)-x2(2*i-1))/(x1(i+1)-x1(i)))*uf(2*i-1,2*j-2)+&
                ((y2(2*j-2)-y1(j-1))/(y1(j)-y1(j-1)))*&
                ((x1(i+1)-x2(2*i))/(x1(i+1)-x1(i)))*uf(2*i,2*j-2)+&
                ((y1(j+1)-y2(2*j-1))/(y1(j+1)-y1(j)))*&
                ((x2(2*i-3)-x1(i-1))/(x1(i)-x1(i-1)))*uf(2*i-3,2*j-1)+&
                ((y1(j+1)-y2(2*j-1))/(y1(j+1)-y1(j)))*&
                ((x2(2*i-2)-x1(i-1))/(x1(i)-x1(i-1)))*uf(2*i-2,2*j-1)+&
                ((y1(j+1)-y2(2*j-1))/(y1(j+1)-y1(j)))*&
                ((x1(i+1)-x2(2*i-1))/(x1(i+1)-x1(i)))*uf(2*i-1,2*j-1)+&
                ((y1(j+1)-y2(2*j-1))/(y1(j+1)-y1(j)))*&
                ((x1(i+1)-x2(2*i))/(x1(i+1)-x1(i)))*uf(2*i,2*j-1)+&
                ((-y2(2*j)+y1(j+1))/(y1(j+1)-y1(j)))*&
                ((x2(2*i-3)-x1(i-1))/(x1(i)-x1(i-1)))*uf(2*i-3,2*j)+&
                ((-y2(2*j)+y1(j+1))/(y1(j+1)-y1(j)))*&
                ((x2(2*i-2)-x1(i-1))/(x1(i)-x1(i-1)))*uf(2*i-2,2*j)+&
                ((-y2(2*j)+y1(j+1))/(y1(j+1)-y1(j)))*&
                ((x1(i+1)-x2(2*i-1))/(x1(i+1)-x1(i)))*uf(2*i-1,2*j)+&
                ((-y2(2*j)+y1(j+1))/(y1(j+1)-y1(j)))*&
                ((x1(i+1)-x2(2*i))/(x1(i+1)-x1(i)))*uf(2*i,2*j))
        enddo
     enddo
     ! Boundary grid points
     rstrct2dp_gen(1,2:ny1-1)=0.5_wp*(uf(1,2:ny2-2:2)+uf(1,3:ny2-1:2))
     rstrct2dp_gen(nx1,2:ny1-1)=0.5_wp*(uf(nx2,2:ny2-2:2)+uf(nx2,3:ny2-1:2))
     rstrct2dp_gen(2:nx1-1,1)=0.5_wp*(uf(2:nx2-2:2,1)+uf(3:nx2-1:2,1))
     rstrct2dp_gen(2:nx1-1,ny1)=0.5_wp*(uf(2:ny2-2:2,ny2)+uf(3:nx1-1:2,ny2))
     rstrct2dp_gen(1,1)=uf(1,1)
     rstrct2dp_gen(1,ny1)=uf(1,ny2)
     rstrct2dp_gen(nx1,ny1)=uf(nx2,ny2)
     rstrct2dp_gen(nx1,1)=uf(nx2,1)
     ! Do neighbor exchange after this to get ghost points on coarser grid
  elseif (grid==1) then
     do j=2,ny1-1
        do i=2,nx1-1
           rstrct2dp_gen(i,j)=(1.0_wp/4.0_wp)*&
                (((y2(2*j-3)-y1(j-1))/(y1(j)-y1(j-1)))*&
                ((x2(2*i-2)-x1(i-1))/(x1(i)-x1(i-1)))*uf(2*i-2,2*j-3)+&
                ((y2(2*j-2)-y1(j-1))/(y1(j)-y1(j-1)))*&
                ((x2(2*i-2)-x1(i-1))/(x1(i)-x1(i-1)))*uf(2*i-2,2*j-2)+&
                ((y1(j+1)-y2(2*j-1))/(y1(j+1)-y1(j)))*&
                ((x2(2*i-2)-x1(i-1))/(x1(i)-x1(i-1)))*uf(2*i-2,2*j-1)+&
                ((y1(j+1)-y2(2*j))/(y1(j+1)-y1(j)))*&
                ((x2(2*i-2)-x1(i-1))/(x1(i)-x1(i-1)))*uf(2*i-2,2*j)+&
                ((y2(2*j-3)-y1(j-1))/(y1(j)-y1(j-1)))*uf(2*i-1,2*j-3)+&
                ((y2(2*j-2)-y1(j-1))/(y1(j)-y1(j-1)))*uf(2*i-1,2*j-2)+&
                ((y1(j+1)-y2(2*j-1))/(y1(j+1)-y1(j)))*uf(2*i-1,2*j-1)+&
                ((y1(j+1)-y2(2*j))/(y1(j+1)-y1(j)))*uf(2*i-1,2*j)+&
                ((y2(2*j-3)-y1(j-1))/(y1(j)-y1(j-1)))*&
                ((x1(i+1)-x2(2*i))/(x1(i+1)-x1(i)))*uf(2*i,2*j-3)+&
                ((y2(2*j-2)-y1(j-1))/(y1(j)-y1(j-1)))*&
                ((x1(i+1)-x2(2*i))/(x1(i+1)-x1(i)))*uf(2*i,2*j-2)+&
                ((y1(j+1)-y2(2*j-1))/(y1(j+1)-y1(j)))*&
                ((x1(i+1)-x2(2*i))/(x1(i+1)-x1(i)))*uf(2*i,2*j-1)+&
                ((y1(j+1)-y2(2*j))/(y1(j+1)-y1(j)))*&
                ((x1(i+1)-x2(2*i))/(x1(i+1)-x1(i)))*uf(2*i,2*j))
        enddo
     enddo
     ! Boundary grid points
     rstrct2dp_gen(1,2:ny1-1)=0.5_wp*(uf(1,2:ny2-2:2)+uf(1,3:ny2-1:2))
     rstrct2dp_gen(nx1,2:ny1-1)=0.5_wp*(uf(nx2,2:ny2-2:2)+uf(nx2,3:ny2-1:2))
     rstrct2dp_gen(2:nx1-1,1)=2.0_wp*(h*uf(3:nx2-1:2,1)+&
          h2*(uf(2:nx2-2:2,1)+uf(4:nx2:2,1)))
     rstrct2dp_gen(2:nx1-1,ny1)=2.0_wp*(h*uf(3:nx2-1:2,ny2)+&
          h2*(uf(2:nx2-2:2,ny2)+uf(4:nx2:2,ny2)))
     rstrct2dp_gen(1,1)=uf(1,1)
     rstrct2dp_gen(1,ny1)=uf(1,ny2)
     rstrct2dp_gen(nx1,ny1)=uf(nx2,ny2)
     rstrct2dp_gen(nx1,1)=uf(nx2,1)
     ! Do neighbor exchange after this to get ghost points on coarser grid
  elseif (grid==2) then
     do j=2,ny1-1
        do i=2,nx1-1
           rstrct2dp_gen(i,j)=(1.0_wp/4.0_wp)*&
                (((y2(2*j-2)-y1(j-1))/(y1(j)-y1(j-1)))*&
                ((x2(2*i-3)-x1(i-1))/(x1(i)-x1(i-1)))*uf(2*i-3,2*j-2)+&
                ((y2(2*j-2)-y1(j-1))/(y1(j)-y1(j-1)))*&
                ((x2(2*i-2)-x1(i-1))/(x1(i)-x1(i-1)))*uf(2*i-2,2*j-2)+&
                ((y2(2*j-2)-y1(j-1))/(y1(j)-y1(j-1)))*&
                ((x1(i)-x2(2*i-1))/(x1(i)-x1(i-1)))*uf(2*i-1,2*j-2)+&
                ((y2(2*j-2)-y1(j-1))/(y1(j)-y1(j-1)))*&
                ((x1(i)-x2(2*i))/(x1(i)-x1(i-1)))*uf(2*i,2*j-2)+&
                ((x2(2*i-3)-x1(i-1))/(x1(i)-x1(i-1)))*uf(2*i-3,2*j-1)+&
                ((x2(2*i-2)-x1(i-1))/(x1(i)-x1(i-1)))*uf(2*i-2,2*j-1)+&
                ((x1(i+1)-x2(2*i-1))/(x1(i+1)-x1(i)))*uf(2*i-1,2*j-1)+&
                ((x1(i+1)-x2(2*i))/(x1(i+1)-x1(i)))*uf(2*i,2*j-1)+&
                ((y1(j+1)-y2(2*j))/(y1(j+1)-y1(j)))*&
                ((x2(2*i-3)-x1(i-1))/(x1(i)-x1(i-1)))*uf(2*i-3,2*j)+&
                ((y1(j+1)-y2(2*j))/(y1(j+1)-y1(j)))*&
                ((x2(2*i-2)-x1(i-1))/(x1(i)-x1(i-1)))*uf(2*i-2,2*j)+&
                ((y1(j+1)-y2(2*j))/(y1(j)-y1(j-1)))*&
                ((x1(i)-x2(2*i-1))/(x1(i)-x1(i-1)))*uf(2*i-1,2*j)+&
                ((y1(j+1)-y2(2*j))/(y1(j)-y1(j-1)))*&
                ((x1(i)-x2(2*i))/(x1(i)-x1(i-1)))*uf(2*i,2*j))
        enddo
     enddo
     ! Boundary grid points
     rstrct2dp_gen(2:nx1-1,1)=0.5_wp*(uf(2:nx2-2:2,1)+uf(3:nx2-1:2,1))
     rstrct2dp_gen(2:nx1-1,ny1)=0.5_wp*(uf(2:nx2-2:2,ny2)+uf(3:nx2-1:2,ny2))
     rstrct2dp_gen(1,2:ny1-1)=2.0_wp*(h*uf(1,3:ny2-1:2)+&
          h2*(uf(1,2:ny2-2:2)+uf(1,4:ny2:2)))
     rstrct2dp_gen(nx1,2:ny1-1)=2.0_wp*(h*uf(nx2,3:ny2-1:2)+&
          h2*(uf(nx2,2:ny2-2:2)+uf(nx2,4:ny2:2)))
     rstrct2dp_gen(1,1)=uf(1,1)
     rstrct2dp_gen(1,ny1)=uf(1,ny2)
     rstrct2dp_gen(nx1,ny1)=uf(nx2,ny2)
     rstrct2dp_gen(nx1,1)=uf(nx2,1)
     ! Do neighbor exchange after this to get ghost points on coarser grid
  endif
END FUNCTION rstrct2dp_gen
