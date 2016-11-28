! Module containing interface blocks for multi-grid functions and subroutines
MODULE bfunc
  INTERFACE 
     SUBROUTINE decomp(s,e,dims,coords,n,bdy1,bdy2)
       USE btype
       IMPLICIT NONE
       INTEGER(I4B), INTENT(IN) :: n,dims,coords
       INTEGER(I4B), INTENT(OUT) :: s,e
       LOGICAL, INTENT(OUT) :: bdy1,bdy2
     END SUBROUTINE decomp
  END INTERFACE
  INTERFACE 
     SUBROUTINE decompu(s,e,dims,coords,n)
       USE btype
       IMPLICIT NONE
       INTEGER(I4B), INTENT(IN) :: n,dims,coords
       INTEGER(I4B), INTENT(OUT) :: s,e
     END SUBROUTINE decompu
  END INTERFACE
  INTERFACE 
     SUBROUTINE mkgrid(sx,ex,sy,ey,asp,ng,nglow,lhbdy,lhbc,rhbdy,rhbc,&
     topbdy,topbc,botbdy,botbc,mynx,myny,x,y,dx,dy)
       USE btype
       IMPLICIT NONE
       TYPE(ptr1d), INTENT(INOUT) :: x(:),y(:),dx(:),dy(:)
       INTEGER(I4B), INTENT(IN) :: sx,ex,sy,ey,ng,nglow
       REAL(WP), INTENT(IN) :: asp
       LOGICAL, INTENT(IN) :: lhbdy,rhbdy,topbdy,botbdy
       CHARACTER(LEN=3), INTENT(IN) :: rhbc,lhbc,topbc,botbc
       INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: mynx,myny
     END SUBROUTINE mkgrid
  END INTERFACE
  INTERFACE 
     SUBROUTINE mglin2d(u,rhs,dtime,ng,nglow,ncycle,mynx,myny,x,y,dx,dy,&
     myid,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
       USE mpi
       USE btype
       USE bfunc2d
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: u
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: rhs
       INTEGER(I4B), INTENT(IN) :: ncycle,myid,nbrleft,nbrright,&
            nbrtop,nbrbottom,ng,nglow,comm2d
       TYPE(ptr1d), INTENT(IN) :: x(:),y(:),dx(:),dy(:)
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
       REAL(WP), INTENT(IN) :: dtime
     END SUBROUTINE mglin2d
  END INTERFACE
END MODULE bfunc
