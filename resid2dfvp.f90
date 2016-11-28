! Calculates residual using the solution u, and forcing function rhs
! Residual can only be calculated on interior points, so may require 
! neighbor exchange to update ghost points
!
! 11/8/2010, finite volume discretization for points in cell centers 
! (temp,press,visc,etc)
!
FUNCTION resid2dfvp(u,rhs,ng,dx,dy,dtime)
  USE btype
  IMPLICIT NONE
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,rhs
  REAL(WP), DIMENSION(:), INTENT(IN) :: dx,dy
  REAL(WP), DIMENSION(size(u,1),size(u,2)) :: resid2dfvp
  INTEGER(I4B), INTENT(IN) :: ng
  REAL(WP), INTENT(IN) :: dtime
  REAL(WP), DIMENSION(size(dx,1)-1,size(dy,1)-1) :: ap
  REAL(WP), DIMENSION(size(dx,1)-1) :: ae,aw
  REAL(WP), DIMENSION(size(dy,1)-1) :: an,as
  INTEGER(I4B) :: n,i,j,nx,ny,b,c
  REAL(WP) :: h,h2,h2i
  h=1.0_wp/(2**ng)
  h2=h*h
  h2i=1.0_wp/h2
  nx=size(u,1)
  ny=size(u,2)
  b=size(ae,1)
  c=size(as,1)

  aw(1:b)=h/dx(1:b)
  ae(1:b)=h/dx(2:b+1)
  as(1:c)=h/dy(1:c)
  an(1:c)=h/dy(2:c+1)

  ! Do calculation
  resid2dfvp=0.0_wp
  if (dtime==0.0) then
     ! All nodes except boundary points
     do j=2,ny-1
        do i=2,nx-1
           resid2dfvp(i,j)=-h2i*(-(ae(i-1)+aw(i-1)+as(j-1)+an(j-1))*u(i,j)&
                +aw(i-1)*u(i-1,j)+ae(i-1)*u(i+1,j)+as(j-1)*u(i,j-1)+&
                an(j-1)*u(i,j+1))+rhs(i,j)
        enddo
     enddo
  else
     do j=2,ny-1
        do i=2,nx-1
           resid2dfvp(i,j)=(h2/dtime+0.5_wp*(ae(i-1)+aw(i-1)+as(j-1)+an(j-1)))&
                *u(i,j)-0.5_wp*(aw(i-1)*u(i-1,j)+ae(i-1)*u(i+1,j)+&
                as(j-1)*u(i,j-1)+an(j-1)*u(i,j+1))-rhs(i,j)
        enddo
     enddo
  endif
  ! Boundary points are either bc's or ghost points, so don't need to 
  ! calculate residual for them
END FUNCTION resid2dfvp
