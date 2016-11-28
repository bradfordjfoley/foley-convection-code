SUBROUTINE mkgridu(s,e,asp,ng,nglow,myn,x,dx,asp_nx)
  USE btype
  IMPLICIT NONE
  TYPE(ptr1d), INTENT(INOUT) :: x(:),dx(:)
  INTEGER(I4B), INTENT(IN) :: s,e,ng,nglow,asp_nx
  REAL(WP), INTENT(IN) :: asp
  INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: myn
  INTEGER(I4B) :: ng1,s1,i,nx1

  ng1=ng
  s1=s
  do 
     if (ng1<nglow) exit
     if (ng1==ng) then
        myn(ng1)=e-s+3
     else
        if (mod(myn(ng1+1),2)==0) then
           myn(ng1)=(myn(ng1+1))/2 + 1
        endif
        if (mod(myn(ng1+1),2)/=0) then
           myn(ng1)=(myn(ng1+1)-1)/2 +1
        endif
     endif
     ! Allocate X on current grid level
     allocate(x(ng1)%d(myn(ng1)))
     ! X on current grid level
     nx1=(2**ng1)*asp_nx
!!$     if (ng1==ng) then
!!$        do i=1,64
!!$           x(ng1)%d(i)=(0.4_wp/64)*(i-1)
!!$        enddo        
!!$        do i=65,192
!!$           x(ng1)%d(i)=(0.2_wp/128)*(i-65)+0.4_wp
!!$        enddo
!!$        do i=193,myn(ng1)
!!$           x(ng1)%d(i)=(0.4_wp/64)*(i-193)+0.6_wp
!!$        enddo
!!$        else
        do i=1,myn(ng1)
           x(ng1)%d(i)=(1.0_wp/(nx1/asp))*(s1+i-3)
        enddo
!!$     endif
     ! Allocate and determine dx on current grid level
     allocate(dx(ng1)%d(myn(ng1)-1))
     do i=1,myn(ng1)-1
        dx(ng1)%d(i)=x(ng1)%d(i+1)-x(ng1)%d(i) 
     enddo
     ! Next coarser grid
     ng1=ng1-1
     ! Adjust starting X and Y points to current grid level
     if (mod(s1,2)==0) then
        s1=s1/2 + 1
     elseif (mod(s1,2)/=0) then
        s1=(s1+1)/2
     endif
  enddo
END SUBROUTINE mkgridu
