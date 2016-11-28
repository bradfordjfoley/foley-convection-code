! Subroutine to calculate strain-rate for use in solving fineness equation
! Developing as of 4/19/12
! Arrays du_dx,du_dy,dw_dx,dw_dy will be defined without ghost points 
! they'll be allocated in the main program
SUBROUTINE strainrate(u,w,dx,dy,deltax,deltay,du_dx,du_dy,dw_dx,dw_dy)
  USE btype
  IMPLICIT NONE
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
  REAL(WP), DIMENSION(:), INTENT(IN) :: dx,dy,deltax,deltay
  REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: du_dx,du_dy,dw_dx,dw_dy
  INTEGER :: nx_sr,ny_sr,i,j

  nx_sr=size(du_dx,1)
  ny_sr=size(du_dx,2)

  do j=1,ny_sr
     do i=1,nx_sr
        du_dx(i,j)=(u(i+1,j+1)-u(i,j+1))*(1.0_wp/deltax(i))
        dw_dy(i,j)=(w(i+1,j+1)-w(i+1,j))*(1.0_wp/deltay(j))
        du_dy(i,j)=0.5_wp*((u(i,j+2)-u(i,j))/(dy(j)+dy(j+1)) + &
             (u(i+1,j+2)-u(i+1,j))/(dy(j)+dy(j+1)))
        dw_dx(i,j)=0.5_wp*((w(i+2,j)-w(i,j))/(dx(i)+dx(i+1)) + &
             (w(i+2,j+1)-w(i,j+1))/(dx(i)+dx(i+1)))
     enddo
  enddo

END SUBROUTINE strainrate
