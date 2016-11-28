!Module containing interface blocks for functions and subroutines 
!that themselves call functions or subroutines with interaces in bfunc2d
MODULE bfunc2d_main
  INTERFACE 
     SUBROUTINE mglin2d(u,rhs,dtime,ng,nglow,grid,ncycle,nrelax,mynx,myny,x,y,&
          ae,aw,an,as,ap0,myid,nbrleft,nbrright,nbrtop,nbrbottom,&
          tfsl,bfsl,lfsl,rfsl,x_con,y_con,value,comm2d)
       !USE mpi
       USE btype
       USE bfunc2d
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: u
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: rhs
       INTEGER(I4B), INTENT(IN) :: ncycle,myid,nbrleft,nbrright,nbrtop,&
            nbrbottom,ng,nglow,comm2d,grid,nrelax
       TYPE(ptr1d), INTENT(IN) :: x(:),y(:)
       TYPE(ptr2d), INTENT(IN) :: ae(:),aw(:),an(:),as(:)
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
       REAL(WP), INTENT(IN) :: dtime,ap0,x_con,y_con,value
       LOGICAL, INTENT(IN) :: tfsl,bfsl,lfsl,rfsl
     END SUBROUTINE mglin2d
  END INTERFACE
  INTERFACE 
     SUBROUTINE mg2d(u,rhs,dtime,ng,nglow,grid,ncycle,nrelax,mynx,myny,x,y,&
          ae,aw,an,as,ap0,myid,nbrleft,nbrright,nbrtop,nbrbottom,tfsl,bfsl,&
          lfsl,rfsl,x_con,y_con,value,comm2d,deltax,deltay,ioff,joff)
       !USE mpi
       USE btype
       USE bfunc2d
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: u
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: rhs
       INTEGER(I4B), INTENT(IN) :: ncycle,myid,nbrleft,nbrright,&
            nbrtop,nbrbottom,ng,nglow,comm2d,grid,nrelax,ioff,joff
       TYPE(ptr1d), INTENT(IN) :: x(:),y(:),deltax(:),deltay(:)
       TYPE(ptr2d), INTENT(IN) :: ae(:),aw(:),an(:),as(:)
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
       REAL(WP), INTENT(IN) :: dtime,ap0,x_con,y_con,value
       LOGICAL, INTENT(IN) :: tfsl,bfsl,lfsl,rfsl
     END SUBROUTINE mg2d
  END INTERFACE
  INTERFACE 
     SUBROUTINE calc_au(ng,nglow,ae,aw,an,as,mu,mynx,myny,mynxu,mynyw,&
     deltax,deltay,dy,x,y,xu,yw,myid,nbrleft,nbrright,nbrtop,&
     nbrbottom,comm2d,lhbdy,lhbc,rhbdy,rhbc,topbdy,topbc,botbdy,botbc)
       !USE mpi
       USE btype
       USE bfunc2d
       IMPLICIT NONE
       INTEGER(I4B), INTENT(IN) :: ng,nglow
       TYPE(ptr2d), INTENT(INOUT) :: ae(:),aw(:),an(:),as(:)
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: mu
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny,mynxu,mynyw
       TYPE(ptr1d), INTENT(IN) :: deltax(:),deltay(:),dy(:),x(:),y(:),&
            xu(:),yw(:)
       INTEGER(I4B), INTENT(IN) :: myid,nbrleft,nbrright,nbrtop,nbrbottom,&
            comm2d
       LOGICAL, INTENT(IN) :: lhbdy,rhbdy,topbdy,botbdy
       CHARACTER(LEN=3), INTENT(IN) :: lhbc,rhbc,topbc,botbc
     END SUBROUTINE calc_au
  END INTERFACE
  INTERFACE 
     SUBROUTINE calc_aw(ng,nglow,ae,aw,an,as,mu,mynx,myny,mynxu,mynyw,&
     deltax,deltay,dx,x,y,xu,yw,myid,nbrleft,nbrright,nbrtop,&
     nbrbottom,comm2d,lhbdy,lhbc,rhbdy,rhbc)
       !USE mpi
       USE btype
       USE bfunc2d
       IMPLICIT NONE
       INTEGER(I4B), INTENT(IN) :: ng,nglow
       TYPE(ptr2d), INTENT(INOUT) :: ae(:),aw(:),an(:),as(:)
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: mu
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny,mynxu,mynyw
       TYPE(ptr1d), INTENT(IN) :: deltax(:),deltay(:),dx(:),x(:),y(:),&
            xu(:),yw(:)
       INTEGER(I4B), INTENT(IN) :: myid,nbrleft,nbrright,nbrtop,nbrbottom,&
            comm2d
       LOGICAL, INTENT(IN) :: lhbdy,rhbdy
       CHARACTER(LEN=3), INTENT(IN) :: lhbc,rhbc
     END SUBROUTINE calc_aw
  END INTERFACE
  INTERFACE 
     SUBROUTINE calc_ap(ae,aw,an,as,ng,nglow,mynx,myny,mynxu,mynyw,aeu,awu,&
          anu,asu,aew,aww,asw,anw,nbrleft,nbrright,nbrtop,nbrbottom,comm2d,&
          deltax,deltay)
       !USE mpi
       USE btype
       USE bfunc2d
       IMPLICIT NONE
       INTEGER(I4B), INTENT(IN) :: ng,nglow,nbrleft,nbrright,nbrtop,&
            nbrbottom,comm2d
       TYPE(ptr2d), INTENT(IN) :: aeu(:),awu(:),anu(:),asu(:),aew(:),aww(:),&
            anw(:),asw(:)
       TYPE(ptr2d), INTENT(INOUT) :: ae(:),aw(:),an(:),as(:)
       TYPE(ptr1d), INTENT(IN) :: deltax(:),deltay(:)
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny,mynxu,mynyw
     END SUBROUTINE calc_ap
  END INTERFACE
  INTERFACE 
     SUBROUTINE mpdata(phi,u,w,dtime,ng,iord,mynx,myny,myid,nbrleft,&
          nbrright,nbrtop,nbrbottom,tfsl,bfsl,comm2d)
       USE btype
       USE bfunc2d
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: phi
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
       INTEGER(I4B), INTENT(IN) :: iord,myid,nbrleft,nbrright,nbrtop,&
            nbrbottom,ng,comm2d
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
       REAL(WP), INTENT(IN) :: dtime
       LOGICAL, INTENT(IN) :: tfsl,bfsl
     END SUBROUTINE mpdata
  END INTERFACE
  INTERFACE 
     SUBROUTINE mpdata_flux(phi,u,w,dtime,ng,iord,mynx,myny,myid,nbrleft,&
          nbrright,nbrtop,nbrbottom,tfsl,bfsl,lfsl,rfsl,comm2d,&
          tbdy,bbdy,tbc,bbc,utfsl,ubfsl,wlfsl,wrfsl)
       USE btype
       USE bfunc2d
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: phi
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
       INTEGER(I4B), INTENT(IN) :: iord,myid,nbrleft,nbrright,nbrtop,&
            nbrbottom,ng,comm2d
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
       REAL(WP), INTENT(IN) :: dtime
       LOGICAL, INTENT(IN) :: tfsl,bfsl,tbdy,bbdy,lfsl,rfsl,utfsl,ubfsl,&
            wlfsl,wrfsl
       CHARACTER(LEN=3), INTENT(IN) :: tbc,bbc
     END SUBROUTINE mpdata_flux
  END INTERFACE
  INTERFACE
     SUBROUTINE mpdata_flux_gen(phi,u,w,dtime,ng,iord,mynx,myny,myid,nbrleft,&
          nbrright,nbrtop,nbrbottom,tfsl,bfsl,lfsl,rfsl,comm2d,tbdy,bbdy,tbc,&
          bbc,utfsl,ubfsl,wlfsl,wrfsl,deltax,deltay,dx,dy,x,y,xu,yw)
       USE btype
       USE bfunc2d
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: phi
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
       REAL(WP), DIMENSION(:), INTENT(IN) :: deltax,deltay,dx,dy,x,y,xu,yw
       INTEGER(I4B), INTENT(IN) :: iord,myid,nbrleft,nbrright,nbrtop,&
            nbrbottom,ng,comm2d
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
       REAL(WP), INTENT(IN) :: dtime
       LOGICAL, INTENT(IN) :: tfsl,bfsl,tbdy,bbdy,lfsl,rfsl,utfsl,ubfsl,&
            wlfsl,wrfsl
       CHARACTER(LEN=3), INTENT(IN) :: tbc,bbc
       END SUBROUTINE mpdata_flux_gen
    END INTERFACE
  INTERFACE 
     SUBROUTINE diffusion(phi,lewis,dtime,Q,ng,nglow,dx,dy,mynx,myny,myid,&
     nbrleft,nbrright,nbrtop,nbrbottom,tfsl,bfsl,lfsl,rfsl,comm2d,x,y,&
     deltax,deltay)
       USE btype
       USE bfunc2d 
       IMPLICIT NONE 
       REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: phi
       TYPE(ptr1d), INTENT(IN) :: dx(:),dy(:),deltax(:),deltay(:)
       REAL(WP), INTENT(IN) :: lewis,dtime,Q 
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
       REAL(WP), DIMENSION(:), INTENT(IN) :: x,y
       INTEGER(I4B), INTENT(IN) :: nbrleft,nbrright,nbrtop,nbrbottom,&
            comm2d,ng,nglow,myid
       LOGICAL, INTENT(IN) :: tfsl,bfsl,lfsl,rfsl
     END SUBROUTINE diffusion
  END INTERFACE
  INTERFACE 
     SUBROUTINE advecdiff(phi,lewis,dtime,u,w,ng,nglow,x,y,dx,dy,mynx,&
          myny,myid,nbrleft,nbrright,nbrtop,nbrbottom,tfsl,bfsl,comm2d)
       USE btype
       USE bfunc2d 
       IMPLICIT NONE 
       REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: phi
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
       TYPE(ptr1d), INTENT(IN) :: x(:),y(:),dx(:),dy(:)
       REAL(WP), INTENT(IN) :: lewis,dtime 
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
       INTEGER(I4B), INTENT(IN) :: nbrleft,nbrright,nbrtop,nbrbottom,&
            comm2d,ng,nglow,myid
       LOGICAL, INTENT(IN) :: tfsl,bfsl
     END SUBROUTINE advecdiff
  END INTERFACE
END MODULE bfunc2d_main
