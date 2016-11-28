! Jacobi iteration solver; takes solution array and forcing 
! function as input, returns updated solution in w
! Needs to be paired with neighbor exchanges to update 
! ghost points for each subdomain after each iteration
! performs one jacobi iteration	
SUBROUTINE jacobi2dp(u,rhs,ng,w,lhbdy,lhbc,rhbdy,rhbc,topbdy,&
     topbc,botbdy,botbc,lbv,rbv,tbv,bbv)
  USE btype
  IMPLICIT NONE
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: u 
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: rhs
  REAL(WP), DIMENSION(:,:), INTENT(OUT) :: w
  INTEGER(I4B), INTENT(IN) :: ng
  LOGICAL, INTENT(IN) :: lhbdy,rhbdy,topbdy,botbdy
  CHARACTER(LEN=3), INTENT(IN) :: lhbc,rhbc,topbc,botbc
  REAL(WP), INTENT(IN) :: lbv,rbv,tbv,bbv
  INTEGER(I4B) :: nx1,ny1
  REAL(WP) :: h,h2,p,p2,p3
  nx1=size(u,1)
  ny1=size(u,2)
  h=1.0_wp/(2**ng)
  h2=h*h
  p=1.0_wp/4.0
  p2=1.0_wp/5.0
  p3=1.0_wp/6.0
  w=0.0_wp
  ! Do only interior grid points of a process's computational domain
  w(3:nx1-2,3:ny1-2)=p*(u(4:nx1-1,3:ny1-2)+u(2:nx1-3,3:ny1-2)+&
       u(3:nx1-2,4:ny1-1)+u(3:nx1-2,2:ny1-3)-h2*rhs(3:nx1-2,3:ny1-2)) 
  ! Now do boundary grid points *********************************************
  ! LEFT ********************************************************************
  if (lhbdy) then
     if (lhbc=='con') then
        w(2,3:ny1-2)=p2*(u(3,3:ny1-2)+2.0_wp*lbv+u(2,4:ny1-1)+&
             u(2,2:ny1-3)-h2*rhs(2,3:ny1-2))
     endif
     if (lhbc=='per') then
        w(2,3:ny1-2)=p*(u(3,3:ny1-2)+u(1,3:ny1-2)+u(2,4:ny1-1)+&
             u(2,2:ny1-3)-h2*rhs(2,3:ny1-2)) 
     endif
  else
     w(2,3:ny1-2)=p*(u(3,3:ny1-2)+u(1,3:ny1-2)+u(2,4:ny1-1)+&
          u(2,2:ny1-3)-h2*rhs(2,3:ny1-2)) 
  endif
  ! RIGHT ********************************************************************
  if (rhbdy) then
     if (rhbc=='con') then
        w(nx1-1,3:ny1-2)=p2*(2.0_wp*rbv+u(nx1-2,3:ny1-2)+&
             u(nx1-1,4:ny1-1)+u(nx1-1,2:ny1-3)-h2*rhs(nx1-1,3:ny1-2))
     endif
     if (rhbc=='per') then
        w(nx1-1,3:ny1-2)=p*(u(nx1,3:ny1-2)+u(nx1-2,3:ny1-2)+&
             u(nx1-1,4:ny1-1)+u(nx1-1,2:ny1-3)-h2*rhs(nx1-1,3:ny1-2))
     endif
  else
     w(nx1-1,3:ny1-2)=p*(u(nx1,3:ny1-2)+u(nx1-2,3:ny1-2)+&
          u(nx1-1,4:ny1-1)+u(nx1-1,2:ny1-3)-h2*rhs(nx1-1,3:ny1-2))
  endif
  ! BOTTOM *******************************************************************
  if (botbdy) then
     if (botbc=='con') then
        w(3:nx1-2,2)=p2*(u(4:nx1-1,2)+u(2:nx1-3,2)+u(3:nx1-2,3)+&
             2.0_wp*bbv-h2*rhs(3:nx1-2,2))
     endif
     if (botbc=='per') then
         w(3:nx1-2,2)=p*(u(4:nx1-1,2)+u(2:nx1-3,2)+u(3:nx1-2,3)+&
             u(3:nx1-2,1)-h2*rhs(3:nx1-2,2))  
      endif
  else
     w(3:nx1-2,2)=p*(u(4:nx1-1,2)+u(2:nx1-3,2)+u(3:nx1-2,3)+&
          u(3:nx1-2,1)-h2*rhs(3:nx1-2,2))
  endif
  ! TOP **********************************************************************
  if (topbdy) then
     if (topbc=='con') then
        w(3:nx1-2,ny1-1)=p2*(u(4:nx1-1,ny1-1)+u(2:nx1-3,ny1-1)+&
             2.0_wp*tbv+u(3:nx1-2,ny1-2)-h2*rhs(3:nx1-2,ny1-1))
     endif
     if (topbc=='per') then
        w(3:nx1-2,ny1-1)=p*(u(4:nx1-1,ny1-1)+u(2:nx1-3,ny1-1)+&
             u(3:nx1-2,ny1)+u(3:nx1-2,ny1-2)-h2*rhs(3:nx1-2,ny1-1))
     endif
  else
     w(3:nx1-2,ny1-1)=p*(u(4:nx1-1,ny1-1)+u(2:nx1-3,ny1-1)+&
          u(3:nx1-2,ny1)+u(3:nx1-2,ny1-2)-h2*rhs(3:nx1-2,ny1-1))
  endif
  ! BOTTOM LEFT ************************************************************* 
  if (lhbdy.and.botbdy) then
     if (lhbc=='con'.and.botbc=='con') then
        w(2,2)=p3*(2.0_wp*lbv+u(3,2)+u(2,3)+2.0_wp*bbv-h2*rhs(2,2))
     endif
     if (lhbc=='per'.and.botbc=='con') then
        w(2,2)=p2*(u(1,2)+u(3,2)+u(2,3)+2.0_wp*bbv-h2*rhs(2,2))
     endif
  elseif (lhbdy.and.botbdy==.false.) then 
     if (lhbc=='con') then
        w(2,2)=p2*(2.0_wp*lbv+u(3,2)+u(2,3)+u(2,1)-h2*rhs(2,2))
     endif
     if (lhbc=='per') then
        w(2,2)=p*(u(1,2)+u(3,2)+u(2,3)+u(2,1)-h2*rhs(2,2))
     endif
  elseif (lhbdy==.false..and.botbdy) then 
     if (botbc=='con') then
        w(2,2)=p2*(u(1,2)+u(3,2)+u(2,3)+2.0_wp*bbv-h2*rhs(2,2))
     endif
     if (botbc=='per') then
        w(2,2)=p*(u(1,2)+u(3,2)+u(2,3)+u(2,1)-h2*rhs(2,2))
     endif
  elseif (lhbdy==.false..and.botbdy==.false.) then 
     w(2,2)=p*(u(1,2)+u(3,2)+u(2,3)+u(2,1)-h2*rhs(2,2))
  endif
! TOP LEFT ***************************************************************** 
  if (lhbdy.and.topbdy) then
     if (lhbc=='con'.and.topbc=='con') then
        w(2,ny1-1)=p3*(2.0_wp*lbv+u(3,ny1-1)+2.0_wp*tbv+&
             u(2,ny1-2)-h2*rhs(2,ny1-1))
     endif
     if (lhbc=='per'.and.topbc=='con') then
        w(2,ny1-1)=p2*(u(1,ny1-1)+u(3,ny1-1)+2.0_wp*tbv+&
             u(2,ny1-2)-h2*rhs(2,ny1-1))
     endif
  elseif (lhbdy.and.topbdy==.false.) then 
     if (lhbc=='con') then
        w(2,ny1-1)=p2*(2.0_wp*lbv+u(3,ny1-1)+u(2,ny1)+&
             u(2,ny1-2)-h2*rhs(2,ny1-1))
     endif
     if (lhbc=='per') then
        w(2,ny1-1)=p*(u(1,ny1-1)+u(3,ny1-1)+u(2,ny1)+&
             u(2,ny1-2)-h2*rhs(2,ny1-1))
     endif
  elseif (lhbdy==.false..and.topbdy) then 
     if (topbc=='con') then
        w(2,ny1-1)=p2*(u(1,ny1-1)+u(3,ny1-1)+u(2,ny1-2)+2.0_wp*tbv&
             -h2*rhs(2,ny1-1))
     endif
     if (topbc=='per') then
        w(2,ny1-1)=p*(u(1,ny1-1)+u(3,ny1-1)+u(2,ny1-2)+u(2,ny1)&
             -h2*rhs(2,ny1-1))
     endif
  elseif (lhbdy==.false..and.topbdy==.false.) then 
     w(2,ny1-1)=p*(u(1,ny1-1)+u(3,ny1-1)+u(2,ny1-2)+u(2,ny1)-h2*rhs(2,ny1-1))
  endif
! TOP RIGHT ***************************************************************** 
  if (rhbdy.and.topbdy) then
     if (rhbc=='con'.and.topbc=='con') then
        w(nx1-1,ny1-1)=p3*(u(nx1-2,ny1-1)+2.0_wp*rbv+&
             2.0_wp*tbv+u(nx1-1,ny1-2)-h2*rhs(nx1-1,ny1-1))
     endif
     if (rhbc=='per'.and.topbc=='con') then
        w(nx1-1,ny1-1)=p2*(u(nx1,ny1-1)+u(nx1-2,ny1-1)+&
             2.0_wp*tbv+u(nx1-1,ny1-2)-h2*rhs(nx1-1,ny1-1))
     endif
  elseif (rhbdy.and.topbdy==.false.) then 
     if (rhbc=='con') then
        w(nx1-1,ny1-1)=p2*(2.0_wp*rbv+u(nx1-2,ny1-1)+u(nx1-1,ny1)+&
             u(nx1-1,ny1-2)-h2*rhs(nx1-1,ny1-1))
     endif
     if (rhbc=='per') then
        w(nx1-1,ny1-1)=p*(u(nx1,ny1-1)+u(nx1-2,ny1-1)+u(nx1-1,ny1)+&
             u(nx1-1,ny1-2)-h2*rhs(nx1-1,ny1-1))
     endif
  elseif (rhbdy==.false..and.topbdy) then 
     if (topbc=='con') then
        w(nx1-1,ny1-1)=p2*(u(nx1,ny1-1)+u(nx1-2,ny1-1)+2.0_wp*tbv+&
             u(nx1-1,ny1-2)-h2*rhs(nx1-1,ny1-1))
     endif
     if (topbc=='per') then
        w(nx1-1,ny1-1)=p*(u(nx1,ny1-1)+u(nx1-2,ny1-1)+u(nx1-1,ny1)+&
             u(nx1-1,ny1-2)-h2*rhs(nx1-1,ny1-1))
     endif
  elseif (rhbdy==.false..and.topbdy==.false.) then 
     w(nx1-1,ny1-1)=p*(u(nx1,ny1-1)+u(nx1-2,ny1-1)+u(nx1-1,ny1)+&
          u(nx1-1,ny1-2)-h2*rhs(nx1-1,ny1-1))
  endif
! BOTTOM RIGHT ************************************************************* 
  if (rhbdy.and.botbdy) then
     if (rhbc=='con'.and.botbc=='con') then
        w(nx1-1,2)=p3*(u(nx1-2,2)+2.0_wp*rbv+2.0_wp*bbv+&
             u(nx1-1,3)-h2*rhs(nx1-1,2))
     endif
     if (rhbc=='per'.and.botbc=='con') then
        w(nx1-1,2)=p2*(u(nx1-2,2)+u(nx1,2)+2.0_wp*bbv+&
             u(nx1-1,3)-h2*rhs(nx1-1,2))
     endif
  elseif (rhbdy.and.botbdy==.false.) then
     if (rhbc=='con') then
        w(nx1-1,2)=p2*(u(nx1-2,2)+2.0_wp*rbv+u(nx1-1,1)+&
             u(nx1-1,3)-h2*rhs(nx1-1,2))  
     endif
     if (rhbc=='per') then
         w(nx1-1,2)=p*(u(nx1-2,2)+u(nx1,2)+u(nx1-1,1)+&
             u(nx1-1,3)-h2*rhs(nx1-1,2)) 
      endif
  elseif (rhbdy==.false..and.botbdy) then 
     if (botbc=='con') then
        w(nx1-1,2)=p2*(u(nx1-2,2)+u(nx1,2)+2.0_wp*bbv+&
             u(nx1-1,3)-h2*rhs(nx1-1,2))  
     endif
     if (botbc=='per') then
         w(nx1-1,2)=p*(u(nx1-2,2)+u(nx1,2)+u(nx1-1,1)+&
             u(nx1-1,3)-h2*rhs(nx1-1,2)) 
      endif
  elseif (rhbdy==.false..and.botbdy==.false.) then
     w(nx1-1,2)=p*(u(nx1-2,2)+u(nx1,2)+u(nx1-1,1)+&
          u(nx1-1,3)-h2*rhs(nx1-1,2)) 
  endif
END SUBROUTINE jacobi2dp
