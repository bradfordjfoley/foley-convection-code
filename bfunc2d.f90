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
          x,y,x_con,y_con,value)
       USE btype
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,ae,aw,an,as
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: rhs
       INTEGER(I4B), INTENT(IN) :: ng
       REAL(WP), DIMENSION(:,:), INTENT(OUT) :: w
       REAL(WP), INTENT(IN) :: dtime,ap0,x_con,y_con,value
       REAL(WP), DIMENSION(:), INTENT(IN) :: x,y
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
     FUNCTION resid2dfvp2(u,rhs,ng,ae,aw,an,as,ap0,dtime,x,y,x_con,y_con)
         USE btype
         IMPLICIT NONE
         REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,rhs,ae,aw,an,as
         REAL(WP), DIMENSION(size(u,1),size(u,2)) :: resid2dfvp2
         INTEGER(I4B), INTENT(IN) :: ng
         REAL(WP), INTENT(IN) :: dtime,ap0,x_con,y_con
         REAL(WP), DIMENSION(:), INTENT(IN) :: x,y
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
          botbcc,mynx,myny,x,xw,xc,y,yu,yc,dx,dxw,dxc,dy,dyu,dyc)
       USE btype
       IMPLICIT NONE
       TYPE(ptr1d), INTENT(INOUT) :: x(:),y(:),xw(:),yu(:),xc(:),yc(:),dx(:),&
            dxw(:),dxc(:),dy(:),dyu(:),dyc(:)
       INTEGER(I4B), INTENT(IN) :: sx,ex,sy,ey,ng,nglow
       REAL(WP), INTENT(IN) :: asp
       LOGICAL, INTENT(IN) :: lhbdy,rhbdy,topbdy,botbdy
       CHARACTER(LEN=3), INTENT(IN) :: rhbc,rhbcw,rhbcc,lhbc,lhbcw,lhbcc,&
            topbc,topbcu,topbcc,botbc,botbcu,botbcc
       INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: mynx,myny
     END SUBROUTINE mkgrid
  END INTERFACE
  INTERFACE 
     SUBROUTINE mkgridu(s,e,asp,ng,nglow,myn,x,dx)
       USE btype
       IMPLICIT NONE
       TYPE(ptr1d), INTENT(INOUT) :: x(:),dx(:)
       INTEGER(I4B), INTENT(IN) :: s,e,ng,nglow
       REAL(WP), INTENT(IN) :: asp
       INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: myn
     END SUBROUTINE mkgridu
  END INTERFACE
  INTERFACE 
     SUBROUTINE mglin2d(u,rhs,dtime,ng,nglow,grid,ncycle,nrelax,mynx,myny,x,y,&
          ae,aw,an,as,ap0,myid,nbrleft,nbrright,nbrtop,nbrbottom,&
          tfsl,bfsl,lfsl,rfsl,x_con,y_con,value,comm2d)
       !USE mpi
       USE btype
       !USE bfunc2d
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
          lfsl,rfsl,x_con,y_con,value,comm2d)
       !USE mpi
       USE btype
       !USE bfunc2d
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: u
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: rhs
       INTEGER(I4B), INTENT(IN) :: ncycle,myid,nbrleft,nbrright,&
            nbrtop,nbrbottom,ng,nglow,comm2d,grid,nrelax
       TYPE(ptr1d), INTENT(IN) :: x(:),y(:)
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
       !USE bfunc2d
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
       !USE bfunc2d
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
          anu,asu,aew,aww,asw,anw,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
       !USE mpi
       USE btype
       !USE bfunc2d
       IMPLICIT NONE
       INTEGER(I4B), INTENT(IN) :: ng,nglow,nbrleft,nbrright,nbrtop,&
            nbrbottom,comm2d
       TYPE(ptr2d), INTENT(IN) :: aeu(:),awu(:),anu(:),asu(:),aew(:),aww(:),&
            anw(:),asw(:)
       TYPE(ptr2d), INTENT(INOUT) :: ae(:),aw(:),an(:),as(:)
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny,mynxu,mynyw
     END SUBROUTINE calc_ap
  END INTERFACE
  INTERFACE 
     SUBROUTINE mpdata(phi,u,w,dtime,ng,iord,mynx,myny,myid,nbrleft,&
          nbrright,nbrtop,nbrbottom,tfsl,bfsl,comm2d)
       USE btype
       !USE bfunc2d
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
       !USE bfunc2d
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
     SUBROUTINE diffusion(phi,lewis,dtime,Q,ng,nglow,dx,dy,mynx,myny,myid,&
     nbrleft,nbrright,nbrtop,nbrbottom,tfsl,bfsl,lfsl,rfsl,comm2d,x,y)
       USE btype
       !USE bfunc2d 
       IMPLICIT NONE 
       REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: phi
       TYPE(ptr1d), INTENT(IN) :: dx(:),dy(:)
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
       !USE bfunc2d 
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
  INTERFACE 
     SUBROUTINE strainrate(u,w,dx,dy,deltax,deltay,du_dx,du_dy,dw_dx,dw_dy)
       USE btype
       IMPLICIT NONE
       REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
       REAL(WP), DIMENSION(:), INTENT(IN) :: dx,dy
       REAL(WP), INTENT(IN) :: deltax,deltay
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
