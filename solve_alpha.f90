! Subroutine to solve for grain-growth/reduction (separately from advection)
!
! Uses a semi-implicit time discretization (Crank-Nicholson) with source 
! term linearization as in Patankar 4.2-5 for the implicit part  
!
SUBROUTINE solve_alpha(alpha,source,sink,u,w,dt,m,p,h)
  USE btype
  IMPLICIT NONE
  REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: alpha
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: source,sink,u,w
  REAL(WP), DIMENSION(:,:), ALLOCATABLE :: sc,s_p,alpha_n,ae,aw,an,as
  REAL(WP), INTENT(IN) :: dt,m,p,h
  INTEGER(I4B) :: i,j,nx,ny,nxa,nya
  REAL(WP) :: ap1
  !
  nx=size(source,1)
  ny=size(source,2)
  nxa=size(alpha,1)
  nya=size(alpha,2)
  !
  allocate(sc(nx,ny),s_p(nx,ny),alpha_n(nxa,nya),ae(nx,ny),aw(nx,ny),&
       an(nx,ny),as(nx,ny))
  sc=0.0_wp; s_p=0.0_wp; alpha_n=1.0_wp; ae=0.0_wp; aw=0.0_wp; an=0.0_wp
  as=0.0_wp
  ap1=1.0_wp/dt
  !
  ! Make linearized source term Sc and Sp (Patankar 4.2-5) from previous
  ! solution for alpha
  forall (i=1:nx,j=1:ny)
     sc(i,j)=(1.0_wp+m)*source(i,j)*alpha(i+1,j+1)**(-m) +& 
          (p-1.0_wp)*sink(i,j)*alpha(i+1,j+1)**p
     s_p(i,j)=-m*source(i,j)*alpha(i+1,j+1)**(-m-1.0_wp) - &
          p*sink(i,j)*alpha(i+1,j+1)**(p-1.0_wp)
     ae(i,j)=max(-u(i+1,j+1),0.0_wp)/h
     aw(i,j)=max(u(i,j+1),0.0_wp)/h
     an(i,j)=max(-w(i+1,j+1),0.0_wp)/h
     as(i,j)=max(w(i+1,j),0.0_wp)/h
  end forall
  ! Now solve for grain-growth/reduction
  forall (i=1:nx,j=1:ny)
     alpha_n(i+1,j+1)=(1.0_wp/ap1)*(-(ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+&
          (1.0_wp/h)*(u(i+1,j+1)-u(i,j+1))+(1.0_wp/h)*(w(i+1,j+1)-w(i+1,j))-&
          ap1-(s_p(i,j)/2.0_wp))*alpha(i+1,j+1)+ae(i,j)*alpha(i+2,j+1)+&
          aw(i,j)*alpha(i,j+1)+an(i,j)*alpha(i+1,j+2)+&
          as(i,j)*alpha(i+1,j)+0.5_wp*sc(i,j)+&
          0.5_wp*(source(i,j)*alpha(i+1,j+1)**(-m)-&
          sink(i,j)*alpha(i+1,j+1)**p))
  end forall
  !
  alpha=alpha_n
  deallocate(sc,s_p,alpha_n)
END SUBROUTINE solve_alpha
