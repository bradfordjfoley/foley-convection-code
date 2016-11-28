SUBROUTINE mkgridu(s,e,asp,ng,nglow,myn,x,dx)
  USE btype
  IMPLICIT NONE
  TYPE(ptr1d), INTENT(INOUT) :: x(:),dx(:)
  INTEGER(I4B), INTENT(IN) :: s,e,ng,nglow
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
     nx1=(2**ng1)*asp
     do i=1,myn(ng1)
        x(ng1)%d(i)=(1.0_wp/(nx1/asp))*(s1+i-3)
     enddo
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
