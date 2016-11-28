! General 2-D interpolation function, using bi-linear interpolation, for 
! interpolating velocity field to tracer location
SUBROUTINE interpgen_tracer(uc,wc,interpgen_u,interpgen_w,x1,y1,x_u,x,&
     y_w,y,nx,ny,sx,sxu,sy,syw)
  USE btype
  IMPLICIT NONE
  ! uc is the velocity field component, defined on the regular grid
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: uc,wc
  ! x1,y1 are the x,y tracer positions to interpolate velocity field to
  REAL(WP), INTENT(IN) :: x1,y1
  REAL(WP), DIMENSION(:), INTENT(IN) :: x_u,x,y_w,y
  INTEGER(I4B), INTENT(IN) :: nx,ny,sx,sy,sxu,syw
  INTEGER(I4B) :: nx1,ny1,nx2,ny2,nx1_u,nx2_u,ny1_w,ny2_w
  REAL(WP), INTENT(OUT) :: interpgen_u,interpgen_w

  ! Find the locations of neighboring grid points 
  nx1_u=floor(x1*ny+3-sxu)
  nx2_u=ceiling(x1*ny+3-sxu)
  nx1=floor((x1+(1/(2*ny)))*ny+2-sx)
  nx2=ceiling((x1+(1/(2*ny)))*ny+2-sx)
  ny1_w=floor(y1*ny+3-syw)
  ny2_w=ceiling(y1*ny+3-syw)
  ny1=floor((y1+(1/(2*ny)))*ny+2-sy)
  ny2=ceiling((y1+(1/(2*ny)))*ny+2-sy)

  ! Initialize as zero
  interpgen_u=0.0_wp
  interpgen_w=0.0_wp
  
  if (nx1_u==nx2_u) then
     interpgen_u=(1.0_wp/(y(ny2)-y(ny1)))*(uc(nx1_u,ny1)*&
          (y(ny2)-y1)+uc(nx1_u,ny2)*(y1-y(ny1)))
     interpgen_w=(1.0_wp/((x(nx2)-x(nx1))*&
          (y_w(ny2_w)-y_w(ny1_w))))*(uc(nx1,ny1_w)*(x(nx2)-x1)*&
          (y_w(ny2_w)-y1)+uc(nx2,ny1_w)*(x1-x(nx1))*(y_w(ny2_w)-y1)+&
          uc(nx1,ny2_w)*(x(nx2)-x1)*(y1-y_w(ny1_w))+&
          uc(nx2,ny2_w)*(x1-x(nx1))*(y1-y_w(ny1_w)))
  elseif (ny1==ny2) then
     interpgen_u=(1.0_wp/(x_u(nx2_u)-x_u(nx1_u)))*(uc(nx1_u,ny1)*&
          (x_u(nx2_u)-x1)+uc(nx2_u,ny1)*(x1-x_u(nx1_u)))
     interpgen_w=(1.0_wp/((x(nx2)-x(nx1))*&
          (y_w(ny2_w)-y_w(ny1_w))))*(uc(nx1,ny1_w)*(x(nx2)-x1)*&
          (y_w(ny2_w)-y1)+uc(nx2,ny1_w)*(x1-x(nx1))*(y_w(ny2_w)-y1)+&
          uc(nx1,ny2_w)*(x(nx2)-x1)*(y1-y_w(ny1_w))+&
          uc(nx2,ny2_w)*(x1-x(nx1))*(y1-y_w(ny1_w)))
  elseif (nx1==nx2) then
     interpgen_u=(1.0_wp/((x_u(nx2_u)-x_u(nx1_u))*&
          (y(ny2)-y(ny1))))*(uc(nx1_u,ny1)*(x_u(nx2_u)-x1)*&
          (y(ny2)-y1)+uc(nx2_u,ny1)*(x1-x_u(nx1_u))*(y(ny2)-y1)+&
          uc(nx1_u,ny2)*(x_u(nx2_u)-x1)*(y1-y(ny1))+&
          uc(nx2_u,ny2)*(x1-x_u(nx1_u))*(y1-y(ny1)))
     interpgen_w=(1.0_wp/(y_w(ny2_w)-y_w(ny1_w)))*(uc(nx1,ny1_w)*&
          (y_w(ny2_w)-y1)+uc(nx1,ny2_w)*(y1-y_w(ny1_w)))
  elseif (ny1_w==ny2_w) then
     interpgen_u=(1.0_wp/((x_u(nx2_u)-x_u(nx1_u))*&
          (y(ny2)-y(ny1))))*(uc(nx1_u,ny1)*(x_u(nx2_u)-x1)*&
          (y(ny2)-y1)+uc(nx2_u,ny1)*(x1-x_u(nx1_u))*(y(ny2)-y1)+&
          uc(nx1_u,ny2)*(x_u(nx2_u)-x1)*(y1-y(ny1))+&
          uc(nx2_u,ny2)*(x1-x_u(nx1_u))*(y1-y(ny1)))
     interpgen_w=(1.0_wp/(x(nx2)-x(nx1)))*(uc(nx1,ny1_w)*&
          (x(nx2)-x1)+uc(nx2,ny1_w)*(x1-x(nx1)))
  else
     interpgen_u=(1.0_wp/((x_u(nx2_u)-x_u(nx1_u))*&
          (y(ny2)-y(ny1))))*(uc(nx1_u,ny1)*(x_u(nx2_u)-x1)*&
          (y(ny2)-y1)+uc(nx2_u,ny1)*(x1-x_u(nx1_u))*(y(ny2)-y1)+&
          uc(nx1_u,ny2)*(x_u(nx2_u)-x1)*(y1-y(ny1))+&
          uc(nx2_u,ny2)*(x1-x_u(nx1_u))*(y1-y(ny1)))
     interpgen_w=(1.0_wp/((x(nx2)-x(nx1))*&
          (y_w(ny2_w)-y_w(ny1_w))))*(uc(nx1,ny1_w)*(x(nx2)-x1)*&
          (y_w(ny2_w)-y1)+uc(nx2,ny1_w)*(x1-x(nx1))*(y_w(ny2_w)-y1)+&
          uc(nx1,ny2_w)*(x(nx2)-x1)*(y1-y_w(ny1_w))+&
          uc(nx2,ny2_w)*(x1-x(nx1))*(y1-y_w(ny1_w)))
  endif
END SUBROUTINE interpgen_tracer
