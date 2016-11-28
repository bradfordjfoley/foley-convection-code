! Subroutine to calculate finite volume discretization coefficients ae, aw, an,
! and as.  Mu is the viscosity (or "conductivity"), deltax and deltay are the 
! control volume sizes, and dx and dy are the grid spacing for the variable 
! that we will eventually solve for 
SUBROUTINE calc_ap(ae,aw,an,as,ng,nglow,mynx,myny,mynxu,mynyw,aeu,awu,anu,&
     asu,aew,aww,asw,anw,nbrleft,nbrright,nbrtop,nbrbottom,comm2d,&
     deltax,deltay)
  USE mpi
  USE btype 
  USE bfunc2d
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: ng,nglow,nbrleft,nbrright,nbrtop,nbrbottom,comm2d
  TYPE(ptr2d), INTENT(IN) :: aeu(:),awu(:),anu(:),asu(:),aew(:),aww(:),&
       anw(:),asw(:)
  TYPE(ptr2d), INTENT(INOUT) :: ae(:),aw(:),an(:),as(:)
  TYPE(ptr1d), INTENT(IN) :: deltax(:),deltay(:)
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny,mynxu,mynyw
  INTEGER(I4B) :: i,j,ngrid,status(MPI_STATUS_SIZE),ierr
  REAL(WP) :: h,h2
!
  allocate(ae(ng)%a(mynx(ng)-2,myny(ng)-2),aw(ng)%a(mynx(ng)-2,myny(ng)-2),&
       an(ng)%a(mynx(ng)-2,myny(ng)-2),as(ng)%a(mynx(ng)-2,myny(ng)-2))
  ae(ng)%a=0.0_wp
  aw(ng)%a=0.0_wp
  an(ng)%a=0.0_wp
  as(ng)%a=0.0_wp
  ! Determine coefficients on finest grid level
  do j=1,myny(ng)-2
     do i=1,mynxu(ng)-2
        ae(ng)%a(i,j)=(deltay(ng)%d(j)**2.0_wp)/&
             (aeu(ng)%a(i,j)+awu(ng)%a(i,j)+anu(ng)%a(i,j)+asu(ng)%a(i,j))
     enddo
  enddo
  aw(ng)%a(2:mynx(ng)-2,1:myny(ng)-2)=ae(ng)%a(1:mynx(ng)-3,1:myny(ng)-2)
  ! Do a neighbor exchange where left column of aw gets right column of ae
  call MPI_SENDRECV(ae(ng)%a(mynx(ng)-2,1:myny(ng)-2),myny(ng)-2,&
       MPI_DOUBLE_PRECISION,nbrright,0,aw(ng)%a(1,1:myny(ng)-2),myny(ng)-2,&
       MPI_DOUBLE_PRECISION,nbrleft,0,comm2d,status,ierr)
  do j=1,mynyw(ng)-2
     do i=1,mynx(ng)-2
        an(ng)%a(i,j)=(deltax(ng)%d(i)**2.0_wp)/&
             (aew(ng)%a(i,j)+aww(ng)%a(i,j)+anw(ng)%a(i,j)+asw(ng)%a(i,j))
     enddo
  enddo
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
     do j=1,myny(ngrid)-2
        do i=1,mynxu(ngrid)-2
           ae(ngrid)%a(i,j)=(deltay(ngrid)%d(j)**2.0_wp)/&
                (aeu(ngrid)%a(i,j)+awu(ngrid)%a(i,j)+&
                anu(ngrid)%a(i,j)+asu(ngrid)%a(i,j))
        enddo
     enddo
     aw(ngrid)%a(2:mynx(ngrid)-2,1:myny(ngrid)-2)=&
          ae(ngrid)%a(1:mynx(ngrid)-3,1:myny(ngrid)-2)
     ! Do a neighbor exchange where left edge of aw gets right edge of ae
     call MPI_SENDRECV(ae(ngrid)%a(mynx(ngrid)-2,1:myny(ngrid)-2),&
          myny(ngrid)-2,MPI_DOUBLE_PRECISION,nbrright,0,&
          aw(ngrid)%a(1,1:myny(ngrid)-2),myny(ngrid)-2,&
          MPI_DOUBLE_PRECISION,nbrleft,0,comm2d,status,ierr)  
     do j=1,mynyw(ngrid)-2
        do i=1,mynx(ngrid)-2
           an(ngrid)%a(i,j)=(deltax(ngrid)%d(i)**2.0_wp)/&
                (aew(ngrid)%a(i,j)+aww(ngrid)%a(i,j)+&
                anw(ngrid)%a(i,j)+asw(ngrid)%a(i,j))
        enddo
     enddo
     as(ngrid)%a(1:mynx(ngrid)-2,2:myny(ngrid)-2)=&
          an(ngrid)%a(1:mynx(ngrid)-2,1:myny(ngrid)-3)
     ! Do a neighbor exchange where bottom row of as gets top row of an
     call MPI_SENDRECV(an(ngrid)%a(1:mynx(ngrid)-2,myny(ngrid)-2),&
          mynx(ngrid)-2,MPI_DOUBLE_PRECISION,nbrtop,0,&
          as(ngrid)%a(1:mynx(ngrid)-2,1),mynx(ngrid)-2,&
          MPI_DOUBLE_PRECISION,nbrbottom,0,comm2d,status,ierr)
  enddo
END SUBROUTINE calc_ap
  
