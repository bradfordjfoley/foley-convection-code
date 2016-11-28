! Subroutine to calculate strain-rate for use in solving fineness equation
! Developing as of 7/24/12
! Using higher order approximatation of derivatives (5 point stencil), though
! boundary grid points will use standard center difference 
! Formulation assumes constant grid spaing deltax and deltay
! Arrays du_dx,du_dy,dw_dx,dw_dy will be defined without ghost points 
! they'll be allocated in the main program
SUBROUTINE strainrate2(u,w,topbdy,botbdy,deltax,deltay,du_dx,du_dy,&
     dw_dx,dw_dy)
  USE btype
  IMPLICIT NONE
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
  LOGICAL, INTENT(IN) :: topbdy,botbdy
  REAL(WP), INTENT(IN) :: deltax,deltay
  REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: du_dx,du_dy,dw_dx,dw_dy
  INTEGER :: nx_sr,ny_sr,i,j

  nx_sr=size(du_dx,1)
  ny_sr=size(du_dx,2)

  ! ASSUMING SIDE BOUNDARIES ARE PERIODIC
  ! Do all grid points but top and bottom
  do j=2,ny_sr-1
     do i=1,nx_sr
        du_dx(i,j)=(-u(i+3,j+2)+27.0_wp*u(i+2,j+2)-27.0_wp*u(i+1,j+2)+&
             u(i,j+2))/(24.0_wp*deltax)
        dw_dy(i,j)=(-w(i+2,j+3)+27.0_wp*w(i+2,j+2)-27.0_wp*w(i+2,j+1)+&
             w(i+2,j))/(24.0_wp*deltay)
        du_dy(i,j)=0.5_wp*((-u(i+2,j+4)+8.0_wp*u(i+2,j+3)-8.0_wp*u(i+2,j+1)+&
             u(i+2,j))/(12.0_wp*deltay) + &
             (-u(i+1,j+4)+8.0_wp*u(i+1,j+3)-8.0_wp*u(i+1,j+1)+&
             u(i+1,j))/(12.0_wp*deltay))
        dw_dx(i,j)=0.5_wp*((-w(i+4,j+2)+8.0_wp*w(i+3,j+2)-8.0_wp*w(i+1,j+2)+&
             w(i,j+2))/(12.0_wp*deltax) + &
              (-w(i+4,j+1)+8.0_wp*w(i+3,j+1)-8.0_wp*w(i+1,j+1)+&
             w(i,j+1))/(12.0_wp*deltax))
     enddo
  enddo

  ! BOTTOM
  if (botbdy) then
     do i=1,nx_sr
        du_dx(i,1)=(-u(i+3,3)+27.0_wp*u(i+2,3)-27.0_wp*u(i+1,3)+&
             u(i,3))/(24.0_wp*deltax)
        dw_dy(i,1)=(w(i+2,3)-w(i+2,2))/deltay
        du_dy(i,1)=0.5_wp*((u(i+2,4)-u(i+2,2))/(2.0_wp*deltay) + &
             (u(i+1,4)-u(i+1,2))/(2.0_wp*deltay))
        dw_dx(i,1)=0.5_wp*((-w(i+4,3)+8.0_wp*w(i+3,3)-8.0_wp*w(i+1,3)+&
             w(i,3))/(12.0_wp*deltax) + &
             (-w(i+4,2)+8.0_wp*w(i+3,2)-8.0_wp*w(i+1,2)+&
             w(i,2))/(12.0_wp*deltax))
     enddo
  else
     do i=1,nx_sr
        du_dx(i,1)=(-u(i+3,3)+27.0_wp*u(i+2,3)-27.0_wp*u(i+1,3)+&
             u(i,3))/(24.0_wp*deltax)
        dw_dy(i,1)=(-w(i+2,4)+27.0_wp*w(i+2,3)-27.0_wp*w(i+2,2)+&
             w(i+2,1))/(24.0_wp*deltay)
        du_dy(i,1)=0.5_wp*((-u(i+2,5)+8.0_wp*u(i+2,4)-8.0_wp*u(i+2,2)+&
             u(i+2,1))/(12.0_wp*deltay) + &
             (-u(i+1,5)+8.0_wp*u(i+1,4)-8.0_wp*u(i+1,2)+&
             u(i+1,1))/(12.0_wp*deltay))
        dw_dx(i,1)=0.5_wp*((-w(i+4,3)+8.0_wp*w(i+3,3)-8.0_wp*w(i+1,3)+&
             w(i,3))/(12.0_wp*deltax) + &
             (-w(i+4,2)+8.0_wp*w(i+3,2)-8.0_wp*w(i+1,2)+&
             w(i,2))/(12.0_wp*deltax))
     enddo
  endif
  ! TOP
  if (topbdy) then
     do i=1,nx_sr
        du_dx(i,ny_sr)=(-u(i+3,ny_sr+2)+27.0_wp*u(i+2,ny_sr+2)-&
             27.0_wp*u(i+1,ny_sr+2)+u(i,ny_sr+2))/(24.0_wp*deltax)
        dw_dy(i,ny_sr)=(w(i+2,ny_sr+2)-w(i+2,ny_sr+1))/deltay
        du_dy(i,ny_sr)=0.5_wp*((u(i+2,ny_sr+3)-u(i+2,ny_sr+1))&
             /(2.0_wp*deltay) + (u(i+1,ny_sr+3)-u(i+1,ny_sr+1))/&
             (2.0_wp*deltay))
        dw_dx(i,ny_sr)=0.5_wp*((-w(i+4,ny_sr+2)+8.0_wp*w(i+3,ny_sr+2)&
             -8.0_wp*w(i+1,ny_sr+2)+w(i,ny_sr+2))/(12.0_wp*deltax) + &
             (-w(i+4,ny_sr+1)+8.0_wp*w(i+3,ny_sr+1)-8.0_wp*w(i+1,ny_sr+1)+&
             w(i,ny_sr+1))/(12.0_wp*deltax))
     enddo
  else
     do i=1,nx_sr
        du_dx(i,ny_sr)=(-u(i+3,ny_sr+2)+27.0_wp*u(i+2,ny_sr+2)-&
             27.0_wp*u(i+1,ny_sr+2)+u(i,ny_sr+2))/(24.0_wp*deltax)
        dw_dy(i,ny_sr)=(-w(i+2,ny_sr+3)+27.0_wp*w(i+2,ny_sr+2)-&
             27.0_wp*w(i+2,ny_sr+1)+w(i+2,ny_sr))/(24.0_wp*deltay)
        du_dy(i,ny_sr)=0.5_wp*((-u(i+2,ny_sr+4)+8.0_wp*u(i+2,ny_sr+3)-&
             8.0_wp*u(i+2,ny_sr+1)+u(i+2,ny_sr))/(12.0_wp*deltay) + &
             (-u(i+1,ny_sr+4)+8.0_wp*u(i+1,ny_sr+3)-8.0_wp*u(i+1,ny_sr+1)+&
             u(i+1,ny_sr))/(12.0_wp*deltay))
        dw_dx(i,j)=0.5_wp*((-w(i+4,ny_sr+2)+8.0_wp*w(i+3,ny_sr+2)-&
             8.0_wp*w(i+1,ny_sr+2)+w(i,ny_sr+2))/(12.0_wp*deltax) + &
              (-w(i+4,ny_sr+1)+8.0_wp*w(i+3,ny_sr+1)-8.0_wp*w(i+1,ny_sr+1)+&
             w(i,ny_sr+1))/(12.0_wp*deltax))
     enddo
  endif
END SUBROUTINE strainrate2
