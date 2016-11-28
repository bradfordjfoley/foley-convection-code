! Module containing interface blocks for multi-grid functions and subroutines
MODULE bfunc2d
  INTERFACE
     SUBROUTINE nbrex2d(ue,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
          nbrbottom,comm2d)
       !USE mpi
       USE btype
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: ue
       INTEGER(I4B), INTENT(IN) :: myid,nbrleft,nbrright,comm2d,&
            nbrtop,nbrbottom
       LOGICAL, INTENT(IN) :: tfsl,bfsl,lfsl,rfsl
     END SUBROUTINE nbrex2d
  END INTERFACE
  INTERFACE
     SUBROUTINE nbrex2d_int(ue,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,&
          nbrright,nbrtop,nbrbottom,comm2d)
       !USE mpi
       USE btype
       IMPLICIT NONE
       INTEGER(I4B), DIMENSION(:,:), INTENT(INOUT) :: ue
       INTEGER(I4B), INTENT(IN) :: myid,nbrleft,nbrright,comm2d,&
            nbrtop,nbrbottom
       LOGICAL, INTENT(IN) :: tfsl,bfsl,lfsl,rfsl
     END SUBROUTINE nbrex2d_int
  END INTERFACE
  INTERFACE
     SUBROUTINE nbrex2d2(ue,tfsl,bfsl,myid,nbrleft,nbrright,nbrtop,&
          nbrbottom,comm2d)
       !USE mpi
       USE btype
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: ue
       INTEGER(I4B), INTENT(IN) :: myid,nbrleft,nbrright,comm2d,&
            nbrtop,nbrbottom
       LOGICAL, INTENT(IN) :: tfsl,bfsl
     END SUBROUTINE nbrex2d2
  END INTERFACE
  INTERFACE
     FUNCTION interpgen2d(uc,x1,y1,x2,y2)
       USE btype
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: uc
       REAL(WP), DIMENSION(:), INTENT(IN) :: x1,y1,x2,y2
       REAL(WP), DIMENSION(size(x2,1),size(y2,1)) :: interpgen2d
     END FUNCTION interpgen2d
  END INTERFACE
  INTERFACE
     FUNCTION interpgen2d_geo(uc,x1,y1,x2,y2)
       USE btype
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: uc
       REAL(WP), DIMENSION(:), INTENT(IN) :: x1,y1,x2,y2
       REAL(WP), DIMENSION(size(x2,1),size(y2,1)) :: interpgen2d_geo
     END FUNCTION interpgen2d_geo
  END INTERFACE
  INTERFACE
     SUBROUTINE jacobi2dp(u,rhs,ng,w,lhbdy,lhbc,rhbdy,rhbc,topbdy,&
          topbc,botbdy,botbc,lbv,rbv,tbv,bbv)
       USE btype
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: u
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: rhs
       INTEGER(I4B), INTENT(IN) :: ng
       REAL(wp), DIMENSION(:,:), INTENT(OUT) :: w
       LOGICAL, INTENT(IN) :: lhbdy,rhbdy,topbdy,botbdy
       CHARACTER(LEN=3), INTENT(IN) :: lhbc,rhbc,topbc,botbc
       REAL(WP), INTENT(IN) :: lbv,rbv,tbv,bbv
     END SUBROUTINE jacobi2dp
  END INTERFACE
  INTERFACE
     SUBROUTINE jacobi2dfvp(u,rhs,ng,w,dx,dy,dtime)
       USE btype
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: u
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: rhs
       INTEGER(I4B), INTENT(IN) :: ng
       REAL(WP), DIMENSION(:), INTENT(IN) :: dx,dy
       REAL(wp), DIMENSION(:,:), INTENT(OUT) :: w
       REAL(WP), INTENT(IN) :: dtime
     END SUBROUTINE jacobi2dfvp
  END INTERFACE
  INTERFACE
     SUBROUTINE jacobi2dfvp2(u,rhs,ng,w,ae,aw,an,as,ap0,dtime,&
          x,y,x_con,y_con,value,deltax,deltay,ioff,joff)
       USE btype
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,ae,aw,an,as
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: rhs
       INTEGER(I4B), INTENT(IN) :: ng,ioff,joff
       REAL(WP), DIMENSION(:,:), INTENT(OUT) :: w
       REAL(WP), INTENT(IN) :: dtime,ap0,x_con,y_con,value
       REAL(WP), DIMENSION(:), INTENT(IN) :: x,y,deltax,deltay
     END SUBROUTINE jacobi2dfvp2
  END INTERFACE
  INTERFACE     
     FUNCTION rstrct2dp(uf,ng,mynx,myny,grid)
       USE btype
       IMPLICIT NONE 
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: uf
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
       INTEGER(I4B), INTENT(IN) :: ng,grid
       REAL(WP), DIMENSION(mynx(ng-1),myny(ng-1)) :: rstrct2dp
     END FUNCTION rstrct2dp
  END INTERFACE
  INTERFACE     
     FUNCTION rstrct2dp_gen(uf,ng,mynx,myny,grid,x1,y1,x2,y2)
       USE btype
       IMPLICIT NONE 
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: uf
       REAL(WP), DIMENSION(:), INTENT(IN) :: x1,y1,x2,y2
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
       INTEGER(I4B), INTENT(IN) :: ng,grid
       REAL(WP), DIMENSION(mynx(ng-1),myny(ng-1)) :: rstrct2dp_gen
     END FUNCTION rstrct2dp_gen
  END INTERFACE
  INTERFACE     
     FUNCTION resid2dfvp(u,rhs,ng,dx,dy,dtime)
         USE btype
         IMPLICIT NONE
         REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,rhs
         REAL(WP), DIMENSION(:), INTENT(IN) :: dx,dy
         REAL(WP), DIMENSION(size(u,1),size(u,2)) :: resid2dfvp
         INTEGER(I4B), INTENT(IN) :: ng
         REAL(WP), INTENT(IN) :: dtime
     END FUNCTION resid2dfvp
  END INTERFACE
  INTERFACE     
     FUNCTION resid2dfvp2(u,rhs,ng,ae,aw,an,as,ap0,dtime,x,y,x_con,y_con,&
          deltax,deltay,ioff,joff)
         USE btype
         IMPLICIT NONE
         REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,rhs,ae,aw,an,as
         REAL(WP), DIMENSION(size(u,1),size(u,2)) :: resid2dfvp2
         INTEGER(I4B), INTENT(IN) :: ng,ioff,joff
         REAL(WP), INTENT(IN) :: dtime,ap0,x_con,y_con
         REAL(WP), DIMENSION(:), INTENT(IN) :: x,y,deltax,deltay
     END FUNCTION resid2dfvp2
  END INTERFACE
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
     SUBROUTINE mkgrid(sx,ex,sy,ey,asp,ng,nglow,lhbdy,lhbc,lhbcw,lhbcc,rhbdy,&
          rhbc,rhbcw,rhbcc,topbdy,topbc,topbcu,topbcc,botbdy,botbc,botbcu,&
          botbcc,mynx,myny,x,xw,xc,y,yu,yc,dx,dxw,dxc,dy,dyu,dyc,asp_nx)
       USE btype
       IMPLICIT NONE
       TYPE(ptr1d), INTENT(INOUT) :: x(:),y(:),xw(:),yu(:),xc(:),yc(:),dx(:),&
            dxw(:),dxc(:),dy(:),dyu(:),dyc(:)
       INTEGER(I4B), INTENT(IN) :: sx,ex,sy,ey,ng,nglow,asp_nx
       REAL(WP), INTENT(IN) :: asp
       LOGICAL, INTENT(IN) :: lhbdy,rhbdy,topbdy,botbdy
       CHARACTER(LEN=3), INTENT(IN) :: rhbc,rhbcw,rhbcc,lhbc,lhbcw,lhbcc,&
            topbc,topbcu,topbcc,botbc,botbcu,botbcc
       INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: mynx,myny
     END SUBROUTINE mkgrid
  END INTERFACE
  INTERFACE 
     SUBROUTINE mkgridu(s,e,asp,ng,nglow,myn,x,dx,asp_nx)
       USE btype
       IMPLICIT NONE
       TYPE(ptr1d), INTENT(INOUT) :: x(:),dx(:)
       INTEGER(I4B), INTENT(IN) :: s,e,ng,nglow,asp_nx
       REAL(WP), INTENT(IN) :: asp
       INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: myn
     END SUBROUTINE mkgridu
  END INTERFACE
  INTERFACE 
     SUBROUTINE strainrate(u,w,dx,dy,deltax,deltay,du_dx,du_dy,dw_dx,dw_dy)
       USE btype
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
       REAL(WP), DIMENSION(:), INTENT(IN) :: dx,dy,deltax,deltay
       REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: du_dx,du_dy,dw_dx,dw_dy
     END SUBROUTINE strainrate
  END INTERFACE
  INTERFACE 
     SUBROUTINE strainrate2(u,w,topbdy,botbdy,deltax,deltay,du_dx,du_dy,&
          dw_dx,dw_dy)
       USE btype
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
       LOGICAL, INTENT(IN) :: topbdy,botbdy
       REAL(WP), INTENT(IN) :: deltax,deltay
       REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: du_dx,du_dy,dw_dx,dw_dy
     END SUBROUTINE strainrate2
  END INTERFACE
  INTERFACE 
     SUBROUTINE fineness(alpha,source,sink,dt,m,p)
       USE btype
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: alpha
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: source,sink
       REAL(WP), INTENT(IN) :: dt,m,p
     END SUBROUTINE fineness
  END INTERFACE
  INTERFACE 
     SUBROUTINE solve_alpha(alpha,source,sink,u,w,dt,m,p,h)
       USE btype
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: alpha
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: source,sink,u,w
       REAL(WP), INTENT(IN) :: dt,m,p,h
     END SUBROUTINE solve_alpha
  END INTERFACE
  INTERFACE 
     SUBROUTINE smooth(phi,weight)
       USE btype
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: phi
       REAL(WP), INTENT(IN) :: weight
     END SUBROUTINE smooth
  END INTERFACE
  INTERFACE 
     SUBROUTINE tracer_advec(x1,y1,x0,y0,deltat,u,w,asp,x,x_u,y,y_w,nx,ny,sx,&
          sxu,sy,syw,periodic)
       USE btype
       IMPLICIT NONE
       REAL(WP), INTENT(OUT) :: x1,y1
       REAL(WP), INTENT(IN) :: x0,y0,deltat,asp
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
       REAL(WP), DIMENSION(:), INTENT(IN) :: x,x_u,y,y_w
       INTEGER(I4B), INTENT(IN) :: nx,ny,sx,sxu,sy,syw
       LOGICAL, INTENT(IN) :: periodic
     END SUBROUTINE tracer_advec
  END INTERFACE
  INTERFACE 
     SUBROUTINE nbrex2d_tracer_counts(tracers_out,tracers_in,myid,nbrleft,&
     nbrright,nbrtop,nbrbottom,nbrtopleft,nbrtopright,&
     nbrbotleft,nbrbotright,comm2d)
       !USE mpi
       USE btype
       IMPLICIT NONE
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: tracers_out
       INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: tracers_in
       INTEGER(I4B), INTENT(IN) :: myid,nbrleft,nbrright,comm2d,nbrtop,&
            nbrbottom,nbrtopleft,nbrtopright,nbrbotleft,nbrbotright
     END SUBROUTINE nbrex2d_tracer_counts
  END INTERFACE
  INTERFACE 
     SUBROUTINE nbrex2d_tracer_dat(tracers,oldtracers,send_count,recv_count,&
          tracers_out,nbr_to,nbr_from,nstart,comm2d)
       !USE mpi
       USE btype
       IMPLICIT NONE
       TYPE(ptr1d), INTENT(INOUT) :: tracers(:)
       TYPE(ptr1d), INTENT(IN) :: oldtracers(:)
       INTEGER(I4B), INTENT(IN) :: nbr_to,nbr_from,comm2d,nstart,send_count,&
            recv_count
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: tracers_out
     END SUBROUTINE nbrex2d_tracer_dat
  END INTERFACE
END MODULE bfunc2d
