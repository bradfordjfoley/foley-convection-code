SUBROUTINE mkgrid(sx,ex,sy,ey,asp,ng,nglow,lhbdy,lhbc,lhbcw,lhbcc,rhbdy,rhbc,&
     rhbcw,rhbcc,topbdy,topbc,topbcu,topbcc,botbdy,botbc,botbcu,botbcc,mynx,&
     myny,x,xw,xc,y,yu,yc,dx,dxw,dxc,dy,dyu,dyc,asp_nx)
  USE btype
  IMPLICIT NONE
  TYPE(ptr1d), INTENT(INOUT) :: x(:),y(:),xw(:),yu(:),xc(:),yc(:),dx(:),&
       dxw(:),dxc(:),dy(:),dyu(:),dyc(:)
  INTEGER(I4B), INTENT(IN) :: sx,ex,sy,ey,ng,nglow,asp_nx
  REAL(WP), INTENT(IN) :: asp
  LOGICAL, INTENT(IN) :: lhbdy,rhbdy,topbdy,botbdy
  CHARACTER(LEN=3), INTENT(IN) :: rhbc,rhbcw,rhbcc,lhbc,lhbcw,lhbcc,topbc,&
       topbcu,topbcc,botbc,botbcu,botbcc
  INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: mynx,myny
  INTEGER(I4B) :: ng1,sx1,sy1,i,nx1,ny1

  ng1=ng
  sx1=sx
  sy1=sy
  do 
     if (ng1<nglow) exit
     if (ng1==ng) then
        mynx(ng1)=ex-sx+3
        myny(ng1)=ey-sy+3
     else
        if (mod(mynx(ng1+1),2)==0) then
           mynx(ng1)=(mynx(ng1+1))/2 + 1
        endif
        if (mod(mynx(ng1+1),2)/=0) then
           mynx(ng1)=(mynx(ng1+1)-1)/2 +1
        endif
        if (mod(myny(ng1+1),2)==0) then
           myny(ng1)=(myny(ng1+1))/2 + 1
        endif
        if (mod(myny(ng1+1),2)/=0) then
           myny(ng1)=(myny(ng1+1)-1)/2 +1
        endif
     endif
     ! Allocate X and Y on current grid level
     allocate(x(ng1)%d(mynx(ng1)),xw(ng1)%d(mynx(ng1)),xc(ng1)%d(mynx(ng1)))
     allocate(y(ng1)%d(myny(ng1)),yu(ng1)%d(myny(ng1)),yc(ng1)%d(myny(ng1)))
     ! X and Y on current grid level
     nx1=(2**ng1)*asp_nx
     ny1=2**ng1
!!$     if (ng1==ng) then
!!$        do i=1,65
!!$           x(ng1)%d(i)=(0.4_wp/64)*(i-1)-0.5_wp*(0.4_wp/64)
!!$        enddo
!!$        do i=1,128
!!$           x(ng1)%d(i+65)=(0.2_wp/128)*i-0.5_wp*(0.2_wp/128)+0.4_wp
!!$        enddo
!!$        do i=1,65
!!$           x(ng1)%d(i+128+65)=(0.4_wp/64)*i-0.5_wp*(0.4_wp/64)+0.6_wp
!!$        enddo
!!$        do i=1,65
!!$           y(ng1)%d(i)=(0.4_wp/64)*(i-1)-0.5_wp*(0.4_wp/64)
!!$        enddo
!!$        do i=1,128
!!$           y(ng1)%d(i+65)=(0.2_wp/128)*i-0.5_wp*(0.2_wp/128)+0.4_wp
!!$        enddo
!!$        do i=1,65
!!$           y(ng1)%d(i+128+65)=(0.4_wp/64)*i-0.5_wp*(0.4_wp/64)+0.6_wp
!!$        enddo
!!$     else
        do i=1,mynx(ng1)
           x(ng1)%d(i)=(1.0_wp/(nx1/asp))*(sx1+i-2)-0.5_wp/(nx1/asp)
        enddo
        do i=1,myny(ng1)
           y(ng1)%d(i)=(1.0_wp/ny1)*(sy1+i-2)-0.5_wp/ny1
        enddo
!!$     endif
     ! Make x-grid for w and y grid for u because u/w can have 
     ! different bc's from T
     xw(ng1)%d=x(ng1)%d
     xc(ng1)%d=x(ng1)%d
     yu(ng1)%d=y(ng1)%d
     yc(ng1)%d=y(ng1)%d
     ! Boundaries on X and Y 
     if (lhbdy.and.lhbc=='con') then
        x(ng1)%d(1)=0.0_wp
     endif
     if (lhbdy.and.lhbcw=='con') then
        xw(ng1)%d(1)=0.0_wp
     endif
     if (lhbdy.and.lhbcc=='con') then
        xc(ng1)%d(1)=0.0_wp
     endif
     if (rhbdy.and.rhbc=='con') then
        x(ng1)%d(mynx(ng1))=1.0_wp*asp
     endif
     if (rhbdy.and.rhbcw=='con') then
        xw(ng1)%d(mynx(ng1))=1.0_wp*asp
     endif
     if (rhbdy.and.rhbcc=='con') then
        xc(ng1)%d(mynx(ng1))=1.0_wp*asp
     endif
     if (botbdy.and.botbc=='con') then
        y(ng1)%d(1)=0.0_wp
     endif
     if (botbdy.and.botbcu=='con') then
        yu(ng1)%d(1)=0.0_wp
     endif
     if (botbdy.and.botbcc=='con') then
        yc(ng1)%d(1)=0.0_wp
     endif
     if (topbdy.and.topbc=='con') then
        y(ng1)%d(myny(ng1))=1.0_wp
     endif
     if (topbdy.and.topbcu=='con') then
        yu(ng1)%d(myny(ng1))=1.0_wp
     endif
     if (topbdy.and.topbcc=='con') then
        yc(ng1)%d(myny(ng1))=1.0_wp
     endif
     ! Allocate and determine dx,dy on current grid level
     allocate(dx(ng1)%d(mynx(ng1)-1),dy(ng1)%d(myny(ng1)-1))
     allocate(dxw(ng1)%d(mynx(ng1)-1),dyu(ng1)%d(myny(ng1)-1))
     allocate(dxc(ng1)%d(mynx(ng1)-1),dyc(ng1)%d(myny(ng1)-1))
     do i=1,mynx(ng1)-1
        dx(ng1)%d(i)=x(ng1)%d(i+1)-x(ng1)%d(i) 
        dxw(ng1)%d(i)=xw(ng1)%d(i+1)-xw(ng1)%d(i)
        dxc(ng1)%d(i)=xc(ng1)%d(i+1)-xc(ng1)%d(i)
     enddo
     do i=1,myny(ng1)-1
        dy(ng1)%d(i)=y(ng1)%d(i+1)-y(ng1)%d(i)
        dyu(ng1)%d(i)=yu(ng1)%d(i+1)-yu(ng1)%d(i)
        dyc(ng1)%d(i)=yc(ng1)%d(i+1)-yc(ng1)%d(i)
     enddo
     ! Next coarser grid
     ng1=ng1-1
     ! Adjust starting X and Y points to current grid level
     if (mod(sx1,2)==0) then
        sx1=sx1/2 + 1
     endif
     if (mod(sx1,2)/=0) then
        sx1=(sx1+1)/2
     endif
     if (mod(sy1,2)==0) then
        sy1=sy1/2 + 1
     endif
     if (mod(sy1,2)/=0) then
        sy1=(sy1+1)/2
     endif
  enddo
END SUBROUTINE mkgrid
