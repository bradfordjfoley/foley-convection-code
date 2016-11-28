! 2d parallel multigrid solver
SUBROUTINE mglin2d(u,rhs,dtime,ng,nglow,grid,ncycle,nrelax,mynx,myny,x,y,&
     ae,aw,an,as,ap0,myid,nbrleft,nbrright,nbrtop,nbrbottom,tfsl,bfsl,&
     lfsl,rfsl,x_con,y_con,value,comm2d)
  USE mpi
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
  TYPE(ptr2d), ALLOCATABLE :: rho(:)
  TYPE(ptr1d), ALLOCATABLE :: x_lbdy(:),y_tbdy(:),x_rbdy(:),y_bbdy(:)
  INTEGER(I4B) :: ngrid,j,jcycle,i,nfinal
  REAL(WP) :: tol
  REAL(WP), DIMENSION(:,:), POINTER :: uj_1,uj
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: wj
!
  ngrid=ng
  allocate(rho(ng),x_rbdy(ng),y_tbdy(ng),x_lbdy(ng),y_bbdy(ng))
  allocate(rho(ngrid)%a(mynx(ngrid),myny(ngrid)))
  allocate(uj(mynx(ngrid),myny(ngrid)))
  if (tfsl .neqv. .true.) then
     allocate(y_tbdy(ng)%d(mynx(ng)))
     y_tbdy(ng)%d(1:mynx(ng))=u(1:mynx(ng),myny(ng))
  elseif (bfsl .neqv. .true.) then
     allocate(y_bbdy(ng)%d(mynx(ng)))
     y_bbdy(ng)%d(1:mynx(ng))=u(1:mynx(ng),1)
  elseif (lfsl .neqv. .true.) then
     allocate(x_lbdy(ng)%d(myny(ng)))
     x_lbdy(ng)%d(1:myny(ng))=u(1,1:myny(ng))
  elseif (rfsl .neqv. .true.) then
     allocate(x_rbdy(ng)%d(myny(ng)))
     x_rbdy(ng)%d(1:myny(ng))=u(mynx(ng),1:myny(ng))
  endif
  ! Set up the right hand side on each grid level as an array of pointers
  rho(ngrid)%a(1:mynx(ngrid),1:myny(ngrid))=rhs(1:mynx(ngrid),1:myny(ngrid))
  uj(1:mynx(ngrid),1:myny(ngrid))=u(1:mynx(ngrid),1:myny(ngrid))
  do 
     if (ngrid<=nglow) exit
     ngrid=ngrid-1
     uj_1=>uj
     allocate(rho(ngrid)%a(mynx(ngrid),myny(ngrid)),uj(mynx(ngrid),&
          myny(ngrid)))
     rho(ngrid)%a=0.0_wp
     ! Form right hand side on each grid level
     rho(ngrid)%a=rstrct2dp(rho(ngrid+1)%a,ngrid+1,mynx,myny,grid)
     call nbrex2d(rho(ngrid)%a,.false.,.false.,.false.,.false.,myid,&
          nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
     ! Restrict old solution to lowest grid level to preserve bc's
     uj=rstrct2dp(uj_1,ngrid+1,mynx,myny,grid)
     deallocate(uj_1)
     call nbrex2d(uj,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
          nbrbottom,comm2d)
     if (tfsl .neqv. .true.) then
        allocate(y_tbdy(ngrid)%d(mynx(ngrid)))
        y_tbdy(ngrid)%d(1:mynx(ngrid))=uj(1:mynx(ngrid),myny(ngrid))
     elseif (bfsl .neqv. .true.) then
        allocate(y_bbdy(ngrid)%d(mynx(ngrid)))
        y_bbdy(ngrid)%d(1:mynx(ngrid))=uj(1:mynx(ngrid),1)
     elseif (lfsl .neqv. .true.) then
        allocate(x_lbdy(ngrid)%d(myny(ngrid)))
        x_lbdy(ngrid)%d(1:myny(ngrid))=uj(1,1:myny(ngrid))
     elseif (rfsl .neqv. .true.) then
        allocate(x_rbdy(ngrid)%d(myny(ngrid)))
        x_rbdy(ngrid)%d(1:myny(ngrid))=uj(mynx(ngrid),1:myny(ngrid))
     endif
  enddo
  allocate(wj(mynx(nglow),myny(nglow)))
  wj=0.0_wp
  do i=1,5 !max(3,numprocs)
     call jacobi2dfvp2(uj,rho(nglow)%a,nglow,wj,ae(nglow)%a,aw(nglow)%a,&
          an(nglow)%a,as(nglow)%a,ap0,dtime,x(nglow)%d,y(nglow)%d,x_con,&
          y_con,value)
     call nbrex2d(wj,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
          nbrbottom,comm2d)
     call jacobi2dfvp2(wj,rho(nglow)%a,nglow,uj,ae(nglow)%a,aw(nglow)%a,&
          an(nglow)%a,as(nglow)%a,ap0,dtime,x(nglow)%d,y(nglow)%d,x_con,&
          y_con,value)
     call nbrex2d(uj,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
          nbrbottom,comm2d)
  enddo
  deallocate(wj)
  ! Now solve at each higher grid level
  do j=nglow+1,ng
     uj_1=>uj
     allocate(uj(mynx(j),myny(j)))
     ! Interpolate from j-1 level to current grid level
     ! Everyone does coarsest grid interpolation 
     uj=interpgen2d(uj_1,x(j-1)%d,y(j-1)%d,x(j)%d,y(j)%d)
     ! This is not general, but works for my decomposition and 
     ! having even number of processors in each direction
     ! 
     ! Use the restricted boundaries from the original input u to 
     ! enforce boundary conditions after interpolation
     if (tfsl .neqv. .true.) then
        uj(1:mynx(j),myny(j))=y_tbdy(j)%d(1:mynx(j))
     elseif (bfsl .neqv. .true.) then
        uj(1:mynx(j),1)=y_bbdy(j)%d(1:mynx(j))
     elseif (lfsl .neqv. .true.) then
        uj(1,1:myny(j))=x_lbdy(j)%d(1:myny(j))
     elseif (rfsl .neqv. .true.) then
        uj(mynx(j),1:myny(j))=x_rbdy(j)%d(1:myny(j))
     endif
     deallocate(uj_1)
     do jcycle=1,ncycle
        ! Call recursive subroutine to obtain solution at each grid level
        call mgp(j,uj,rho(j)%a)
     enddo
  enddo
!!$  ! Determine if solution has converged
!!$  allocate(wj(mynx(ng),myny(ng)))
!!$  tol=1e-6
!!$  nfinal=1000
!!$  i=0
!!$  do while (i.lt.nfinal)
!!$     i=i+1
!!$     call jacobi2dfvp(uj,rhs,ng,wj,dx(ng)%d,dy(ng)%d,dtime)
!!$     call nbrex2d(wj,myid,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
!!$     call jacobi2dfvp(wj,rhs,ng,uj,dx(ng)%d,dy(ng)%d,dtime)
!!$     call nbrex2d(uj,myid,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
!!$     if (maxval(abs(uj-wj)).lt.tol) exit
!!$  enddo
  u=uj
  do j=nglow,ng
     deallocate(rho(j)%a)
     if (tfsl .neqv. .true.) then
        deallocate(y_tbdy(j)%d)
     elseif (bfsl .neqv. .true.) then
        deallocate(y_bbdy(j)%d)
     elseif (lfsl .neqv. .true.) then
        deallocate(x_lbdy(j)%d)
     elseif (rfsl .neqv. .true.) then
        deallocate(x_rbdy(j)%d)
     endif
  enddo
  deallocate(rho)
  deallocate(uj)
  if (tfsl) then
     deallocate(y_tbdy)
  elseif (bfsl) then
     deallocate(y_bbdy)
  elseif (lfsl) then
     deallocate(x_lbdy)
  elseif (rfsl) then
     deallocate(x_rbdy)
  endif

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
               ap0,dtime,x(j)%d,y(j)%d,x_con,y_con,value)
          call nbrex2d(w,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
               nbrbottom,comm2d)
          call jacobi2dfvp2(w,rhs,j,u,ae(j)%a,aw(j)%a,an(j)%a,as(j)%a,&
               ap0,dtime,x(j)%d,y(j)%d,x_con,y_con,value)
          call nbrex2d(u,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
               nbrbottom,comm2d)
       enddo
    else
       ! pre-relaxation
       do jpre=1,j*NPRE
          call jacobi2dfvp2(u,rhs,j,w,ae(j)%a,aw(j)%a,an(j)%a,as(j)%a,&
               ap0,dtime,x(j)%d,y(j)%d,x_con,y_con,value)
          call nbrex2d(w,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
               nbrbottom,comm2d)
          call jacobi2dfvp2(w,rhs,j,u,ae(j)%a,aw(j)%a,an(j)%a,as(j)%a,&
               ap0,dtime,x(j)%d,y(j)%d,x_con,y_con,value)
          call nbrex2d(u,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
               nbrbottom,comm2d)
       enddo
       ! calculate the residual
       res1=resid2dfvp2(u,rhs,j,ae(j)%a,aw(j)%a,an(j)%a,as(j)%a,ap0,dtime,&
            x(j)%d,y(j)%d,x_con,y_con)   
       call nbrex2d(res1,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
            nbrbottom,comm2d)
       res2=rstrct2dp(res1,j,mynx,myny,grid)
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
               ap0,dtime,x(j)%d,y(j)%d,x_con,y_con,value)
          call nbrex2d(w,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
               nbrbottom,comm2d)
          call jacobi2dfvp2(w,rhs,j,u,ae(j)%a,aw(j)%a,an(j)%a,as(j)%a,&
               ap0,dtime,x(j)%d,y(j)%d,x_con,y_con,value)
          call nbrex2d(u,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
               nbrbottom,comm2d)
       end do
    end if
  END SUBROUTINE mgp
END SUBROUTINE mglin2d
