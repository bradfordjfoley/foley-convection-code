PROGRAM tracer_test
  USE btype
  USE bfunc2d
  IMPLICIT NONE
  INTEGER(I4B) :: sx,sy,sxu,syw,nx,ny,i,j,nsteps
  REAL(WP), DIMENSION(:), ALLOCATABLE :: x,y,xu,yw
  REAL(WP), DIMENSION(:,:), ALLOCATABLE :: u,w
  REAL(WP) :: psi_max,umax,wmax,vmax_tot,deltat,deltax,time,x0,y0,x1,y1

  ! 
  sx=1; sy=1; sxu=2; syw=2
  nx=128; ny=128 
  allocate(x(nx+2),y(ny+2),xu(nx+1),yw(ny+1))
  do i=1,nx+2
     x(i)=(1.0_wp/ny)*(sx+i-2)-0.5_wp/ny
  enddo
  do i=1,ny+2
     y(i)=(1.0_wp/ny)*(sx+i-2)-0.5_wp/ny
  enddo
  do i=1,nx+1
     xu(i)=(1.0_wp/ny)*(sxu+i-3)
  enddo
  do i=1,ny+1
     yw(i)=(1.0_wp/ny)*(syw+i-3)
  enddo

  ! Calculate velocity field
  allocate(u(nx+1,ny+2),w(nx+2,ny+1))
  psi_max=250.0_wp/(pi**3.0_wp)
  do i=1,nx+1
     do j=1,ny+2
        u(i,j)=(psi_max*pi)*sin(pi*xu(i))*cos(pi*y(j))
     enddo
  enddo
  do i=1,nx+2
     do j=1,ny+1
        w(i,j)=-(psi_max*pi)*cos(pi*x(i))*sin(pi*yw(j))
     enddo
  enddo
  
  ! Timestep size
  umax=maxval(u)
  wmax=maxval(w)
  vmax_tot=max(umax,wmax)
  deltax=x(2)-x(1)
  deltat=(deltax/vmax_tot)*0.25_wp
  nsteps=nint(0.2203712_wp/deltat)
  
  ! Now set up and do tracer advection
  x0=0.5_wp; y0=0.0159221_wp
  time=0
  do i=1,nsteps
     call tracer_advec(x1,y1,x0,y0,deltat,u,w,x,xu,y,yw,nx,ny,sx,sxu,sy,syw)
     x0=x1; y0=y1; time=time+deltat
  enddo
  print*,x1,y1,time

  deallocate(x,y,xu,yw,u,w)
END PROGRAM tracer_test
