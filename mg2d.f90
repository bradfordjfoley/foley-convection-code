! 2d parallel multigrid solver that does specified number of v cycles at 
! a certain grid level 
!
! NOTE THIS IS NOT A FULL MULTIGRID SOLVER
!
SUBROUTINE mg2d(u,rhs,dtime,ng,nglow,grid,ncycle,nrelax,mynx,myny,x,y,ae,aw,&
     an,as,ap0,myid,nbrleft,nbrright,nbrtop,nbrbottom,tfsl,bfsl,&
     lfsl,rfsl,x_con,y_con,value,comm2d,deltax,deltay,ioff,joff)
  USE mpi
  USE btype
  USE bfunc2d
  IMPLICIT NONE
  REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: u
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: rhs
  INTEGER(I4B), INTENT(IN) :: ncycle,myid,nbrleft,nbrright,nbrtop,&
       nbrbottom,ng,nglow,comm2d,grid,nrelax,ioff,joff
  TYPE(ptr1d), INTENT(IN) :: x(:),y(:),deltax(:),deltay(:)
  TYPE(ptr2d), INTENT(IN) :: ae(:),aw(:),an(:),as(:)
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
  REAL(WP), INTENT(IN) :: dtime,ap0,x_con,y_con,value
  LOGICAL, INTENT(IN) :: tfsl,bfsl,lfsl,rfsl
  INTEGER(I4B) :: j
  REAL(WP) :: tol
  !
  do j=1,ncycle
     call mgp(ng,u,rhs)
  end do
  ! Recursive routine obtains the solution at each grid level 
CONTAINS
  RECURSIVE SUBROUTINE mgp(j,u,rhs)
    USE mpi
    USE btype
    USE bfunc2d
    IMPLICIT NONE
    REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: u
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: rhs
    INTEGER(I4B), INTENT(IN) :: j
    INTEGER(I4B) :: NPRE,NPOST
    INTEGER(I4B) :: jpost,jpre
    REAL(WP), DIMENSION(mynx(j-1),myny(j-1)) :: res2,v
    REAL(WP), DIMENSION(size(u,1),size(u,2)) :: w,res1
    NPRE=nrelax
    NPOST=nrelax
    if (j == nglow) then
       do jpre=1,10
          call jacobi2dfvp2(u,rhs,j,w,ae(j)%a,aw(j)%a,an(j)%a,as(j)%a,&
               ap0,dtime,x(j)%d,y(j)%d,x_con,y_con,value,&
               deltax(j)%d,deltay(j)%d,ioff,joff)
          call nbrex2d(w,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
               nbrbottom,comm2d)
          call jacobi2dfvp2(w,rhs,j,u,ae(j)%a,aw(j)%a,an(j)%a,as(j)%a,&
               ap0,dtime,x(j)%d,y(j)%d,x_con,y_con,value,&
               deltax(j)%d,deltay(j)%d,ioff,joff)
          call nbrex2d(u,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
               nbrbottom,comm2d)
       enddo
    else
       ! pre-relaxation
       do jpre=1,j*NPRE
          call jacobi2dfvp2(u,rhs,j,w,ae(j)%a,aw(j)%a,an(j)%a,as(j)%a,&
               ap0,dtime,x(j)%d,y(j)%d,x_con,y_con,value,&
               deltax(j)%d,deltay(j)%d,ioff,joff)
          call nbrex2d(w,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
               nbrbottom,comm2d)
          call jacobi2dfvp2(w,rhs,j,u,ae(j)%a,aw(j)%a,an(j)%a,as(j)%a,&
               ap0,dtime,x(j)%d,y(j)%d,x_con,y_con,value,&
               deltax(j)%d,deltay(j)%d,ioff,joff)
          call nbrex2d(u,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
               nbrbottom,comm2d)
       enddo
       ! calculate the residual
       res1=resid2dfvp2(u,rhs,j,ae(j)%a,aw(j)%a,an(j)%a,as(j)%a,ap0,dtime,&
            x(j)%d,y(j)%d,x_con,y_con,deltax(j)%d,deltay(j)%d,ioff,joff)   
       call nbrex2d(res1,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
            nbrbottom,comm2d)
       res2=rstrct2dp_gen(res1,j,mynx,myny,grid,x(j-1)%d,y(j-1)%d,&
            x(j)%d,y(j)%d)
!!$       res2=rstrct2dp(res1,j,mynx,myny,grid)
       call nbrex2d(res2,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
            nbrbottom,comm2d)
       v=0.0_wp
       ! Use residual as rhs for calculating the coarse grid correction
       call mgp(j-1,v,res2)
       ! Interpolation of coarse grid correction improves the current solution
       ! Don't need to communicate because we exchanged neighbor
       ! points after each jacobi iteration
       u=u+interpgen2d(v,x(j-1)%d,y(j-1)%d,x(j)%d,y(j)%d)
       ! No communication needed after interp 
       ! Post relaxtion to resolve short-wavelength features
       do jpost=1,j*NPOST
          call jacobi2dfvp2(u,rhs,j,w,ae(j)%a,aw(j)%a,an(j)%a,as(j)%a,&
               ap0,dtime,x(j)%d,y(j)%d,x_con,y_con,value,&
               deltax(j)%d,deltay(j)%d,ioff,joff)
          call nbrex2d(w,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
               nbrbottom,comm2d)
          call jacobi2dfvp2(w,rhs,j,u,ae(j)%a,aw(j)%a,an(j)%a,as(j)%a,&
               ap0,dtime,x(j)%d,y(j)%d,x_con,y_con,value,&
               deltax(j)%d,deltay(j)%d,ioff,joff)
          call nbrex2d(u,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
               nbrbottom,comm2d)
       end do
    end if
  END SUBROUTINE mgp
END SUBROUTINE mg2d
