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
SUBROUTINE jacobi2dfvp(u,rhs,ng,w,dx,dy,dtime)
  USE btype
  IMPLICIT NONE
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: u
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: rhs
  INTEGER(I4B), INTENT(IN) :: ng
  REAL(WP), DIMENSION(:), INTENT(IN) :: dx,dy
  REAL(WP), DIMENSION(:,:), INTENT(OUT) :: w
  REAL(WP), INTENT(IN) :: dtime
  REAL(WP), DIMENSION(size(dx,1)-1) :: ae,aw
  REAL(WP), DIMENSION(size(dy,1)-1) :: an,as
  INTEGER(I4B) :: nx1,ny1,i,j,b,c
  REAL(WP) :: h,h2
  nx1=size(u,1)
  ny1=size(u,2)
  b=size(ae,1)
  c=size(as,1)
  h=1.0_wp/(2**ng)
  h2=h*h
  w=0.0_wp

  aw(1:b)=h/dx(1:b)
  ae(1:b)=h/dx(2:b+1)
  as(1:c)=h/dy(1:c)
  an(1:c)=h/dy(2:c+1)  

  if (dtime==0.0) then
  ! All nodes except boundary points
     do j=2,ny1-1
        do i=2,nx1-1
           w(i,j)=(1.0_wp/(ae(i-1)+aw(i-1)+as(j-1)+an(j-1)))*&
                (aw(i-1)*u(i-1,j)+ae(i-1)*u(i+1,j)+as(j-1)*u(i,j-1)+&
                an(j-1)*u(i,j+1)-h2*rhs(i,j))
        enddo
     enddo
  else 
     do j=2,ny1-1
        do i=2,nx1-1
           w(i,j)=0.5_wp/(h2/dtime+0.5_wp*(ae(i-1)+aw(i-1)+as(j-1)+an(j-1)))*&
                (aw(i-1)*u(i-1,j)+ae(i-1)*u(i+1,j)+as(j-1)*u(i,j-1)+&
                an(j-1)*u(i,j+1)+2.0_wp*rhs(i,j))
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
END SUBROUTINE jacobi2dfvp
