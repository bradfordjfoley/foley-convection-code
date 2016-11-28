! Subroutine to calculate finite volume discretization coefficients ae, aw, an,
! and as.  Mu is the viscosity (or "conductivity"), deltax and deltay are the 
! control volume sizes, and dx and dy are the grid spacing for the variable 
! that we will eventually solve for 
SUBROUTINE calc_ap(ae,aw,an,as,ng,nglow,mynx,myny,mynxu,mynyw,aeu,awu,anu,&
     asu,aew,aww,asw,anw,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
  USE mpi
  USE btype 
  USE bfunc2d
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: ng,nglow,nbrleft,nbrright,nbrtop,nbrbottom,comm2d
  TYPE(ptr2d), INTENT(IN) :: aeu(:),awu(:),anu(:),asu(:),aew(:),aww(:),&
       anw(:),asw(:)
  TYPE(ptr2d), INTENT(INOUT) :: ae(:),aw(:),an(:),as(:)
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny,mynxu,mynyw
  INTEGER(I4B) :: i,j,ngrid,status(MPI_STATUS_SIZE),ierr
  REAL(WP) :: h,h2
!
  h=1.0_wp/(2**ng)
  h2=h*h
  allocate(ae(ng)%a(mynx(ng)-2,myny(ng)-2),aw(ng)%a(mynx(ng)-2,myny(ng)-2),&
       an(ng)%a(mynx(ng)-2,myny(ng)-2),as(ng)%a(mynx(ng)-2,myny(ng)-2))
  ae(ng)%a=0.0_wp
  aw(ng)%a=0.0_wp
  an(ng)%a=0.0_wp
  as(ng)%a=0.0_wp
  ! Determine coefficients on finest grid level
  ae(ng)%a(1:mynxu(ng)-2,1:myny(ng)-2)=h2/&
       (aeu(ng)%a(1:mynxu(ng)-2,1:myny(ng)-2)&
       +awu(ng)%a(1:mynxu(ng)-2,1:myny(ng)-2)&
       +anu(ng)%a(1:mynxu(ng)-2,1:myny(ng)-2)&
       +asu(ng)%a(1:mynxu(ng)-2,1:myny(ng)-2))
  aw(ng)%a(2:mynx(ng)-2,1:myny(ng)-2)=ae(ng)%a(1:mynx(ng)-3,1:myny(ng)-2)
  ! Do a neighbor exchange where left column of aw gets right column of ae
  call MPI_SENDRECV(ae(ng)%a(mynx(ng)-2,1:myny(ng)-2),myny(ng)-2,&
       MPI_DOUBLE_PRECISION,nbrright,0,aw(ng)%a(1,1:myny(ng)-2),myny(ng)-2,&
       MPI_DOUBLE_PRECISION,nbrleft,0,comm2d,status,ierr)
  an(ng)%a(1:mynx(ng)-2,1:mynyw(ng)-2)=h2/&
       (aew(ng)%a(1:mynx(ng)-2,1:mynyw(ng)-2)&
       +aww(ng)%a(1:mynx(ng)-2,1:mynyw(ng)-2)&
       +anw(ng)%a(1:mynx(ng)-2,1:mynyw(ng)-2)&
       +asw(ng)%a(1:mynx(ng)-2,1:mynyw(ng)-2))
  as(ng)%a(1:mynx(ng)-2,2:myny(ng)-2)=an(ng)%a(1:mynx(ng)-2,1:myny(ng)-3)
  ! Do a neighbor exchange where bottom row of as gets top row of an
  call MPI_SENDRECV(an(ng)%a(1:mynx(ng)-2,myny(ng)-2),mynx(ng)-2,&
       MPI_DOUBLE_PRECISION,nbrtop,0,as(ng)%a(1:mynx(ng)-2,1),mynx(ng)-2,&
       MPI_DOUBLE_PRECISION,nbrbottom,0,comm2d,status,ierr)
  ngrid=ng
  do 
     if (ngrid<=nglow) exit
     ngrid=ngrid-1
     h=1.0_wp/(2**ngrid)
     h2=h*h
     allocate(ae(ngrid)%a(mynx(ngrid)-2,myny(ngrid)-2),&
          aw(ngrid)%a(mynx(ngrid)-2,myny(ngrid)-2),&
          an(ngrid)%a(mynx(ngrid)-2,myny(ngrid)-2),&
          as(ngrid)%a(mynx(ngrid)-2,myny(ngrid)-2))
     ae(ngrid)%a=0.0_wp
     aw(ngrid)%a=0.0_wp
     an(ngrid)%a=0.0_wp
     as(ngrid)%a=0.0_wp
     ! Determine coefficients on finest grid level
     ae(ngrid)%a(1:mynxu(ngrid)-2,1:myny(ngrid)-2)=h2/&
          (aeu(ngrid)%a(1:mynxu(ngrid)-2,1:myny(ngrid)-2)&
          +awu(ngrid)%a(1:mynxu(ngrid)-2,1:myny(ngrid)-2)&
          +anu(ngrid)%a(1:mynxu(ngrid)-2,1:myny(ngrid)-2)&
          +asu(ngrid)%a(1:mynxu(ngrid)-2,1:myny(ngrid)-2))
     aw(ngrid)%a(2:mynx(ngrid)-2,1:myny(ngrid)-2)=&
          ae(ngrid)%a(1:mynx(ngrid)-3,1:myny(ngrid)-2)
     ! Do a neighbor exchange where left edge of aw gets right edge of ae
     call MPI_SENDRECV(ae(ngrid)%a(mynx(ngrid)-2,1:myny(ngrid)-2),&
          myny(ngrid)-2,MPI_DOUBLE_PRECISION,nbrright,0,&
          aw(ngrid)%a(1,1:myny(ngrid)-2),myny(ngrid)-2,&
          MPI_DOUBLE_PRECISION,nbrleft,0,comm2d,status,ierr)    
     an(ngrid)%a(1:mynx(ngrid)-2,1:mynyw(ngrid)-2)=h2/&
          (aew(ngrid)%a(1:mynx(ngrid)-2,1:mynyw(ngrid)-2)&
          +aww(ngrid)%a(1:mynx(ngrid)-2,1:mynyw(ngrid)-2)&
          +anw(ngrid)%a(1:mynx(ngrid)-2,1:mynyw(ngrid)-2)&
          +asw(ngrid)%a(1:mynx(ngrid)-2,1:mynyw(ngrid)-2))
     as(ngrid)%a(1:mynx(ngrid)-2,2:myny(ngrid)-2)=&
          an(ngrid)%a(1:mynx(ngrid)-2,1:myny(ngrid)-3)
     ! Do a neighbor exchange where bottom row of as gets top row of an
     call MPI_SENDRECV(an(ngrid)%a(1:mynx(ngrid)-2,myny(ngrid)-2),&
          mynx(ngrid)-2,MPI_DOUBLE_PRECISION,nbrtop,0,&
          as(ngrid)%a(1:mynx(ngrid)-2,1),mynx(ngrid)-2,&
          MPI_DOUBLE_PRECISION,nbrbottom,0,comm2d,status,ierr)
  enddo
END SUBROUTINE calc_ap
  
