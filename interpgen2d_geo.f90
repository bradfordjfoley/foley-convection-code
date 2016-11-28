! General 2-D interpolation function, using geometric interpolation
FUNCTION interpgen2d_geo(uc,x1,y1,x2,y2)
  USE btype
  IMPLICIT NONE
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: uc
  REAL(WP), DIMENSION(:), INTENT(IN) :: x1,y1,x2,y2
  INTEGER(I4B) :: nx1,ny1,nx2,ny2,i,j,l,m
  REAL(WP), DIMENSION(size(x2,1),size(y2,1)) :: interpgen2d_geo
  REAL(WP) :: c1,c2,c3,c4

  nx1=size(uc,1)
  ny1=size(uc,2)
  nx2=size(interpgen2d_geo,1)
  ny2=size(interpgen2d_geo,2)

  ! Initialize as zeroes
  interpgen2d_geo=0.0_wp
  
  m=2
  do j=1,ny2
     if (y2(j).gt.y1(ny1).or.y2(j).lt.y2(1)) exit
     l=2
     do i=1,nx2
        if (x2(i).gt.x1(nx1).or.x2(i).lt.x1(1)) exit
        if (x2(i).gt.x1(l)) then
           l=l+1
        endif
        if (y2(j).gt.y1(m)) then
           m=m+1
        endif
        if (x2(i)==x1(l-1).and.y2(j)==y1(m-1)) then
           interpgen2d_geo(i,j)=uc(l-1,m-1)
        elseif (x2(i)==x1(l).and.y2(j)==y1(m-1)) then
           interpgen2d_geo(i,j)=uc(l,m-1)
        elseif (x2(i)==x1(l).and.y2(j)==y1(m)) then
           interpgen2d_geo(i,j)=uc(l,m)
        elseif (x2(i)==x1(l-1).and.y2(j)==y1(m)) then
           interpgen2d_geo(i,j)=uc(l-1,m)
        elseif (x2(i)==x1(l-1)) then
           c1=(y1(m)-y2(j))/(y1(m)-y1(m-1))
           c2=(y2(j)-y1(m-1))/(y1(m)-y1(m-1))
           interpgen2d_geo(i,j)=exp(c1*log(uc(l-1,m-1))+c2*log(uc(l-1,m)))
        elseif (x2(i)==x1(l)) then
           c1=(y1(m)-y2(j))/(y1(m)-y1(m-1))
           c2=(y2(j)-y1(m-1))/(y1(m)-y1(m-1))
           interpgen2d_geo(i,j)=exp(c1*log(uc(l,m-1))+c2*log(uc(l,m)))
        elseif (y2(j)==y1(m-1)) then
           c1=(x1(l)-x2(i))/(x1(l)-x1(l-1))
           c2=(x2(i)-x1(l-1))/(x1(l)-x1(l-1))
           interpgen2d_geo(i,j)=exp(c1*log(uc(l-1,m-1))+c2*log(uc(l,m-1)))
        elseif (y2(j)==y1(m)) then
           c1=(x1(l)-x2(i))/(x1(l)-x1(l-1))
           c2=(x2(i)-x1(l-1))/(x1(l)-x1(l-1))
           interpgen2d_geo(i,j)=exp(c1*log(uc(l-1,m))+c2*log(uc(l,m)))
        else   
           c1=((x1(l)-x2(i))*(y1(m)-y2(j)))/&
                ((x1(l)-x1(l-1))*(y1(m)-y1(m-1)))
           c2=((x2(i)-x1(l-1))*(y1(m)-y2(j)))/&
                ((x1(l)-x1(l-1))*(y1(m)-y1(m-1)))
           c3=((x2(i)-x1(l-1))*(y2(j)-y1(m-1)))/&
                ((x1(l)-x1(l-1))*(y1(m)-y1(m-1)))
           c4=((x1(l)-x2(i))*(y2(j)-y1(m-1)))/&
                ((x1(l)-x1(l-1))*(y1(m)-y1(m-1)))
           interpgen2d_geo(i,j)=exp(c1*log(uc(l-1,m-1))+c2*log(uc(l,m-1))+&
                c3*log(uc(l,m))+c4*log(uc(l-1,m)))
        endif
     enddo
  enddo
END FUNCTION interpgen2d_geo
