! Jacobi iteration solver; takes solution array and forcing 
! function as input, returns updated solution in w
! Needs to be paired with neighbor exchanges to update 
! ghost points for each subdomain after each iteration
! performs one jacobi iteration	
! 
! Set up for finite volume discretization, here with variable at cell center
! (like p,temp,visc,etc)
! Calculates coefficents from grid spacing on the fly, so accepts uneven grid
! spacing as long as control volumes are same size (maybe modify for fully
! uneven grid in the future)
!
! 11/8/2010
!
SUBROUTINE jacobi2dfvp2(u,rhs,ng,w,ae,aw,an,as,ap0,dtime,&
     x,y,x_con,y_con,value)
  USE btype
  IMPLICIT NONE
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,ae,aw,an,as
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: rhs
  INTEGER(I4B), INTENT(IN) :: ng
  REAL(WP), DIMENSION(:), INTENT(IN) :: x,y
  REAL(WP), DIMENSION(:,:), INTENT(OUT) :: w
  REAL(WP), INTENT(IN) :: dtime,ap0,x_con,y_con,value
  INTEGER(I4B) :: nx1,ny1,i,j,b,c
  REAL(WP) :: h,h2
  nx1=size(u,1)
  ny1=size(u,2)
  b=size(ae,1)
  c=size(as,1)
  h=1.0_wp/(2**ng)
  h2=h*h
  w=0.0_wp

  if (dtime==0.0_wp) then
     ! Works for steady-state problems (with ap0 set to 0) and for fully 
     ! implicit time differencing (with ap0 set to dx*dy/dt and 
     ! rhs = -(1/dt)*u_p)
     ! Unfortunately, using this for implicit time-dependent problems means 
     ! setting dtime=0 in the call to multigrid (or jacobi), which is confusing
     ! and SHOULD BE FIXED
     ! All nodes except boundary points
     do j=2,ny1-1
        do i=2,nx1-1
           if (x(i)/=x_con.or.y(j)/=y_con) then
              w(i,j)=(1.0_wp/(ae(i-1,j-1)+aw(i-1,j-1)+as(i-1,j-1)+&
                   an(i-1,j-1)+ap0))*(aw(i-1,j-1)*u(i-1,j)+&
                   ae(i-1,j-1)*u(i+1,j)+as(i-1,j-1)*u(i,j-1)+&
                   an(i-1,j-1)*u(i,j+1)-h2*rhs(i,j))
           else
              w(i,j)=value
           endif
        enddo
     enddo
  else
     ! Crank-Nicholson Time Integration
     do j=2,ny1-1
        do i=2,nx1-1
           if (x(i)/=x_con.or.y(j)/=y_con) then
              w(i,j)=0.5_wp/(h2/dtime+0.5_wp*(ae(i-1,j-1)+aw(i-1,j-1)+&
                   as(i-1,j-1)+an(i-1,j-1)))*(aw(i-1,j-1)*u(i-1,j)+&
                   ae(i-1,j-1)*u(i+1,j)+as(i-1,j-1)*u(i,j-1)+&
                   an(i-1,j-1)*u(i,j+1)+2.0_wp*rhs(i,j))
           else
              w(i,j)=value
           endif
        enddo
     enddo
  endif 

  ! Boundary conditions
  ! Copy all boundary grid points from old solution to new solution
  ! This will preserve boundary conditions, and will be written over 
  ! after communication step for ghost points
  w(1:nx1,1)=u(1:nx1,1)
  w(1:nx1,ny1)=u(1:nx1,ny1)
  w(1,1:ny1)=u(1,1:ny1)
  w(nx1,1:ny1)=u(nx1,1:ny1)
END SUBROUTINE jacobi2dfvp2
