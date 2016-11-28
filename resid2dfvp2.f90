! Calculates residual using the solution u, and forcing function rhs
! Residual can only be calculated on interior points, so may require 
! neighbor exchange to update ghost points
!
! 11/8/2010, finite volume discretization for points in cell centers 
! (temp,press,visc,etc)
!
FUNCTION resid2dfvp2(u,rhs,ng,ae,aw,an,as,ap0,dtime,x,y,x_con,y_con,&
     deltax,deltay,ioff,joff)
  USE btype
  IMPLICIT NONE
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,rhs,ae,aw,an,as
  REAL(WP), DIMENSION(size(u,1),size(u,2)) :: resid2dfvp2
  REAL(WP), DIMENSION(:), INTENT(IN) :: x,y,deltax,deltay
  INTEGER(I4B), INTENT(IN) :: ng,ioff,joff
  REAL(WP), INTENT(IN) :: dtime,ap0,x_con,y_con
  INTEGER(I4B) :: n,i,j,nx,ny,b,c
  REAL(WP) :: h,h2,h2i
  h=1.0_wp/(2**ng)
  h2=h*h
  h2i=1.0_wp/h2
  nx=size(u,1)
  ny=size(u,2)

  ! Do calculation
  resid2dfvp2=0.0_wp
  if (dtime==0.0) then
     ! All nodes except boundary points
     do j=2,ny-1
        do i=2,nx-1
           if (x(i)/=x_con.or.y(j)/=y_con) then
              resid2dfvp2(i,j)=(-1.0_wp/(deltax(i-ioff)*deltay(j-joff)))*&
                   (-(ae(i-1,j-1)+aw(i-1,j-1)+as(i-1,j-1)+&
                   an(i-1,j-1)+ap0)*u(i,j)+aw(i-1,j-1)*u(i-1,j)+&
                   ae(i-1,j-1)*u(i+1,j)+as(i-1,j-1)*u(i,j-1)+&
                   an(i-1,j-1)*u(i,j+1))+rhs(i,j)
           else
              resid2dfvp2(i,j)=0.0_wp
           endif
        enddo
     enddo
  else
     do j=2,ny-1
        do i=2,nx-1
           if (x(i)/=x_con.or.y(j)/=y_con) then
              resid2dfvp2(i,j)=((deltax(i-ioff)*deltay(j-joff))/dtime+&
                   0.5_wp*(ae(i-1,j-1)+aw(i-1,j-1)+&
                   as(i-1,j-1)+an(i-1,j-1)))*u(i,j)-&
                   0.5_wp*(aw(i-1,j-1)*u(i-1,j)&
                   +ae(i-1,j-1)*u(i+1,j)+as(i-1,j-1)*u(i,j-1)+&
                   an(i-1,j-1)*u(i,j+1))-rhs(i,j)
           else
              resid2dfvp2(i,j)=0.0_wp
           endif
        enddo
     enddo
  endif
  ! Boundary points are either bc's or ghost points, so don't need to 
  ! calculate residual for them
END FUNCTION resid2dfvp2
