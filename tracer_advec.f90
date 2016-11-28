! Runge Kutta scheme for advecting tracers 
SUBROUTINE tracer_advec(x1,y1,x0,y0,deltat,u,w,asp,x,x_u,y,y_w,nx,ny,sx,&
     sxu,sy,syw,periodic)
  USE btype
  IMPLICIT NONE
  REAL(WP), INTENT(OUT) :: x1,y1
  REAL(WP), INTENT(IN) :: x0,y0,deltat,asp
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
  REAL(WP), DIMENSION(:), INTENT(IN) :: x,x_u,y,y_w
  INTEGER(I4B), INTENT(IN) :: nx,ny,sx,sxu,sy,syw
  LOGICAL, INTENT(IN) :: periodic
  REAL(WP) :: k1_x,k1_y,k2_x,k2_y,k3_x,k3_y,k4_x,k4_y,x_new1,y_new1,&
       x_new2,y_new2,x_new3,y_new3

  ! Input: timestep size (deltat), current x,y position (x0,y0), 
  ! velocity field (u,w)
  call interpgen_tracer(u,w,k1_x,k1_y,x0,y0,x_u,x,y_w,y,nx,ny,sx,sxu,&
       sy,syw,periodic)
  x_new1=x0+(deltat/2.0_wp)*k1_x
  y_new1=y0+(deltat/2.0_wp)*k1_y
  call interpgen_tracer(u,w,k2_x,k2_y,x_new1,y_new1,x_u,x,y_w,y,nx,ny,sx,sxu,&
       sy,syw,periodic) 
  x_new2=x0+(deltat/2.0_wp)*k2_x
  y_new2=y0+(deltat/2.0_wp)*k2_y
  call interpgen_tracer(u,w,k3_x,k3_y,x_new2,y_new2,x_u,x,y_w,y,nx,ny,sx,sxu,&
       sy,syw,periodic)
  x_new3=x0+deltat*k3_x
  y_new3=y0+deltat*k3_y
  call interpgen_tracer(u,w,k4_x,k4_y,x_new3,y_new3,x_u,x,y_w,y,nx,ny,sx,sxu,&
       sy,syw,periodic)
  x1=x0+(deltat/6.0_wp)*(k1_x+2.0_wp*k2_x+2.0_wp*k3_x+k4_x)
  if (x1<0.0_wp) then
     x1=asp+x1
  elseif (x1>asp) then
     x1=x1-asp   
  endif
  y1=y0+(deltat/6.0_wp)*(k1_y+2.0_wp*k2_y+2.0_wp*k3_y+k4_y)

CONTAINS

  ! General 2-D interpolation function, using bi-linear interpolation, for 
  ! interpolating velocity field to tracer location
  SUBROUTINE interpgen_tracer(uc,wc,interpgen_u,interpgen_w,x1,y1,x_u,x,&
       y_w,y,nx,ny,sx,sxu,sy,syw,periodic)
    USE btype
    IMPLICIT NONE
    ! uc is the velocity field component, defined on the regular grid
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: uc,wc
    ! x1,y1 are the x,y tracer positions to interpolate velocity field to
    REAL(WP), INTENT(IN) :: x1,y1
    REAL(WP), DIMENSION(:), INTENT(IN) :: x_u,x,y_w,y
    INTEGER(I4B), INTENT(IN) :: nx,ny,sx,sy,sxu,syw
    LOGICAL, INTENT(IN) :: periodic
    INTEGER(I4B) :: nx1,ny1,nx2,ny2,nx1_u,nx2_u,ny1_w,ny2_w
    REAL(WP), INTENT(OUT) :: interpgen_u,interpgen_w

    ! Find the locations of neighboring grid points 
    if (periodic) then
       nx1_u=floor(x1*ny+4.0_wp-sxu)
       nx2_u=ceiling(x1*ny+4.0_wp-sxu)
    else
       nx1_u=floor(x1*ny+3.0_wp-sxu)
       nx2_u=ceiling(x1*ny+3.0_wp-sxu)
    endif
    nx1=floor((x1+(1.0_wp/(2.0_wp*ny)))*ny+2.0_wp-sx)
    nx2=ceiling((x1+(1.0_wp/(2.0_wp*ny)))*ny+2.0_wp-sx)
    ny1_w=floor(y1*ny+3.0_wp-syw)
    ny2_w=ceiling(y1*ny+3.0_wp-syw)
    ny1=floor((y1+(1.0_wp/(2.0_wp*ny)))*ny+2.0_wp-sy)
    ny2=ceiling((y1+(1.0_wp/(2.0_wp*ny)))*ny+2.0_wp-sy)

    ! Initialize as zero
    interpgen_u=0.0_wp
    interpgen_w=0.0_wp

    if (nx1_u==nx2_u.and.ny1==ny2) then
       interpgen_u=uc(nx1_u,ny1)
       interpgen_w=(1.0_wp/((x(nx2)-x(nx1))*&
            (y_w(ny2_w)-y_w(ny1_w))))*(wc(nx1,ny1_w)*(x(nx2)-x1)*&
            (y_w(ny2_w)-y1)+wc(nx2,ny1_w)*(x1-x(nx1))*(y_w(ny2_w)-y1)+&
            wc(nx1,ny2_w)*(x(nx2)-x1)*(y1-y_w(ny1_w))+&
            wc(nx2,ny2_w)*(x1-x(nx1))*(y1-y_w(ny1_w)))
    elseif (nx1_u==nx2_u.and.ny1_w==ny2_w) then
       interpgen_u=(1.0_wp/(y(ny2)-y(ny1)))*(uc(nx1_u,ny1)*&
            (y(ny2)-y1)+uc(nx1_u,ny2)*(y1-y(ny1)))
       interpgen_w=(1.0_wp/(x(nx2)-x(nx1)))*(wc(nx1,ny1_w)*&
            (x(nx2)-x1)+wc(nx2,ny1_w)*(x1-x(nx1)))
    elseif (ny1_w==ny2_w.and.nx1==nx2) then
       interpgen_u=(1.0_wp/((x_u(nx2_u)-x_u(nx1_u))*&
            (y(ny2)-y(ny1))))*(uc(nx1_u,ny1)*(x_u(nx2_u)-x1)*&
            (y(ny2)-y1)+uc(nx2_u,ny1)*(x1-x_u(nx1_u))*(y(ny2)-y1)+&
            uc(nx1_u,ny2)*(x_u(nx2_u)-x1)*(y1-y(ny1))+&
            uc(nx2_u,ny2)*(x1-x_u(nx1_u))*(y1-y(ny1)))
       interpgen_w=wc(nx1,ny1_w)
    elseif (nx1_u==nx2_u) then
       interpgen_u=(1.0_wp/(y(ny2)-y(ny1)))*(uc(nx1_u,ny1)*&
            (y(ny2)-y1)+uc(nx1_u,ny2)*(y1-y(ny1)))
       interpgen_w=(1.0_wp/((x(nx2)-x(nx1))*&
            (y_w(ny2_w)-y_w(ny1_w))))*(wc(nx1,ny1_w)*(x(nx2)-x1)*&
            (y_w(ny2_w)-y1)+wc(nx2,ny1_w)*(x1-x(nx1))*(y_w(ny2_w)-y1)+&
            wc(nx1,ny2_w)*(x(nx2)-x1)*(y1-y_w(ny1_w))+&
            wc(nx2,ny2_w)*(x1-x(nx1))*(y1-y_w(ny1_w)))
    elseif (ny1==ny2) then
       interpgen_u=(1.0_wp/(x_u(nx2_u)-x_u(nx1_u)))*(uc(nx1_u,ny1)*&
            (x_u(nx2_u)-x1)+uc(nx2_u,ny1)*(x1-x_u(nx1_u)))
       interpgen_w=(1.0_wp/((x(nx2)-x(nx1))*&
            (y_w(ny2_w)-y_w(ny1_w))))*(wc(nx1,ny1_w)*(x(nx2)-x1)*&
            (y_w(ny2_w)-y1)+wc(nx2,ny1_w)*(x1-x(nx1))*(y_w(ny2_w)-y1)+&
            wc(nx1,ny2_w)*(x(nx2)-x1)*(y1-y_w(ny1_w))+&
            wc(nx2,ny2_w)*(x1-x(nx1))*(y1-y_w(ny1_w)))
    elseif (nx1==nx2) then
       interpgen_u=(1.0_wp/((x_u(nx2_u)-x_u(nx1_u))*&
            (y(ny2)-y(ny1))))*(uc(nx1_u,ny1)*(x_u(nx2_u)-x1)*&
            (y(ny2)-y1)+uc(nx2_u,ny1)*(x1-x_u(nx1_u))*(y(ny2)-y1)+&
            uc(nx1_u,ny2)*(x_u(nx2_u)-x1)*(y1-y(ny1))+&
            uc(nx2_u,ny2)*(x1-x_u(nx1_u))*(y1-y(ny1)))
       interpgen_w=(1.0_wp/(y_w(ny2_w)-y_w(ny1_w)))*(wc(nx1,ny1_w)*&
            (y_w(ny2_w)-y1)+wc(nx1,ny2_w)*(y1-y_w(ny1_w)))
    elseif (ny1_w==ny2_w) then
       interpgen_u=(1.0_wp/((x_u(nx2_u)-x_u(nx1_u))*&
            (y(ny2)-y(ny1))))*(uc(nx1_u,ny1)*(x_u(nx2_u)-x1)*&
            (y(ny2)-y1)+uc(nx2_u,ny1)*(x1-x_u(nx1_u))*(y(ny2)-y1)+&
            uc(nx1_u,ny2)*(x_u(nx2_u)-x1)*(y1-y(ny1))+&
            uc(nx2_u,ny2)*(x1-x_u(nx1_u))*(y1-y(ny1)))
       interpgen_w=(1.0_wp/(x(nx2)-x(nx1)))*(wc(nx1,ny1_w)*&
            (x(nx2)-x1)+wc(nx2,ny1_w)*(x1-x(nx1)))
    else
       interpgen_u=(1.0_wp/((x_u(nx2_u)-x_u(nx1_u))*&
            (y(ny2)-y(ny1))))*(uc(nx1_u,ny1)*(x_u(nx2_u)-x1)*&
            (y(ny2)-y1)+uc(nx2_u,ny1)*(x1-x_u(nx1_u))*(y(ny2)-y1)+&
            uc(nx1_u,ny2)*(x_u(nx2_u)-x1)*(y1-y(ny1))+&
            uc(nx2_u,ny2)*(x1-x_u(nx1_u))*(y1-y(ny1)))
       interpgen_w=(1.0_wp/((x(nx2)-x(nx1))*&
            (y_w(ny2_w)-y_w(ny1_w))))*(wc(nx1,ny1_w)*(x(nx2)-x1)*&
            (y_w(ny2_w)-y1)+wc(nx2,ny1_w)*(x1-x(nx1))*(y_w(ny2_w)-y1)+&
            wc(nx1,ny2_w)*(x(nx2)-x1)*(y1-y_w(ny1_w))+&
            wc(nx2,ny2_w)*(x1-x(nx1))*(y1-y_w(ny1_w)))
    endif
  END SUBROUTINE interpgen_tracer

END SUBROUTINE tracer_advec
