! Subroutine to solve for grain-growth/reduction (separately from advection)
!
! Uses a semi-implicit time discretization (Crank-Nicholson) with source 
! term linearization as in Patankar 4.2-5 for the implicit part  
!
SUBROUTINE fineness(alpha,source,sink,dt,m,p)
  USE btype
  IMPLICIT NONE
  REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: alpha
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: source,sink
  REAL(WP), DIMENSION(:,:), ALLOCATABLE :: sc,s_p,alpha_n
  REAL(WP), INTENT(IN) :: dt,m,p
  INTEGER(I4B) :: i,j,nx,ny,nxa,nya
  !
  nx=size(source,1)
  ny=size(source,2)
  nxa=size(alpha,1)
  nya=size(alpha,2)
  !
  allocate(sc(nx,ny),s_p(nx,ny),alpha_n(nxa,nya))
  sc=0.0_wp; s_p=0.0_wp; alpha_n=1.0_wp
  !
  ! Make linearized source term Sc and Sp (Patankar 4.2-5) from previous
  ! solution for alpha
  forall (i=1:nx,j=1:ny)
     sc(i,j)=(1.0_wp+m)*source(i,j)*alpha(i+1,j+1)**(-m) +& 
          (p-1.0_wp)*sink(i,j)*alpha(i+1,j+1)**p
     s_p(i,j)=-m*source(i,j)*alpha(i+1,j+1)**(-m-1.0_wp) - &
          p*sink(i,j)*alpha(i+1,j+1)**(p-1.0_wp)
  end forall
  ! Now solve for grain-growth/reduction
  forall (i=1:nx,j=1:ny)
     alpha_n(i+1,j+1)=1.0_wp/(1.0_wp/dt - s_p(i,j)/2.0_wp)*&
          ((1.0_wp/dt)*alpha(i+1,j+1) + 0.5_wp*sc(i,j) + &
          0.5_wp*(source(i,j)*alpha(i+1,j+1)**(-m) - &
          sink(i,j)*alpha(i+1,j+1)**p))
  end forall
  !
  alpha=alpha_n
  deallocate(sc,s_p,alpha_n)
END SUBROUTINE fineness
