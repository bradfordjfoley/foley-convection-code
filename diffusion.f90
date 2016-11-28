! 2-D diffusion solver using center-differencing in space and crank-nicholson
! time discretization 
SUBROUTINE diffusion(phi,lewis,dtime,Q,ng,nglow,dx,dy,mynx,myny,myid,&
     nbrleft,nbrright,nbrtop,nbrbottom,tfsl,bfsl,lfsl,rfsl,comm2d,x,y,&
     deltax,deltay)
  USE btype
  USE bfunc2d 
  USE mpi
  IMPLICIT NONE 
  REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: phi
  TYPE(ptr1d), INTENT(IN) :: dx(:),dy(:),deltax(:),deltay(:)
  REAL(WP), INTENT(IN) :: lewis,dtime,Q 
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
  INTEGER(I4B), INTENT(IN) :: nbrleft,nbrright,nbrtop,nbrbottom,comm2d,ng,&
       nglow,myid
  LOGICAL, INTENT(IN) :: tfsl,bfsl,lfsl,rfsl
  REAL(WP), DIMENSION(:), INTENT(IN) :: x,y
  TYPE(ptr2d), ALLOCATABLE :: ae(:),aw(:),an(:),as(:) 
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: phi_n,rhs 
  INTEGER(I4B) :: i,j,n,nfinal,ngrid,ierr
  REAL(WP) :: tol,tolphi,errphi,tolphi_tot,errphi_tot
! 
  ! Allocate space for FV coefficients 
  allocate(ae(ng),aw(ng),an(ng),as(ng))
  ! Calculate finite volume coefficients
  call calc_ad(ae,aw,an,as,lewis,ng,ng,dx,dy,mynx,myny,deltax,deltay)
  ! Build rhs using solution previous time step 
!!$  h=1.0_wp/(2.0_wp**ng)
!!$  h2=h*h
  allocate(rhs(mynx(ng),myny(ng)))
  rhs=0.0_wp
  do j=2,myny(ng)-1
     do i=2,mynx(ng)-1
        rhs(i,j)=0.5_wp*(ae(ng)%a(i-1,j-1)*phi(i+1,j)+&
             aw(ng)%a(i-1,j-1)*phi(i-1,j)+an(ng)%a(i-1,j-1)*phi(i,j+1)+&
             as(ng)%a(i-1,j-1)*phi(i,j-1))+&
             ((deltax(ng)%d(i-1)*deltay(ng)%d(j-1))/dtime - &
             0.5_wp*(ae(ng)%a(i-1,j-1)+aw(ng)%a(i-1,j-1)+&
             an(ng)%a(i-1,j-1)+as(ng)%a(i-1,j-1)))*phi(i,j)+&
             (deltax(ng)%d(i-1)*deltay(ng)%d(j-1))*Q
     enddo
  enddo
  ! Now call multi-grid to solve for diffusion
!!$  call mg2d(phi,rhs,dtime,ng,nglow,0,2,mynx,myny,x,y,ae,aw,an,as,myid,&
!!$       nbrleft,nbrright,nbrtop,nbrbottom,tfsl,bfsl,comm2d)
  ! Check if solution has converged 
  allocate(phi_n(mynx(ng),myny(ng)))
  phi_n=0.0_wp
  n=0
  nfinal=1000
  tol=1e-10
  do while (n.lt.nfinal) 
     n=n+1
     call jacobi2dfvp2(phi,rhs,ng,phi_n,ae(ng)%a,aw(ng)%a,an(ng)%a,&
          as(ng)%a,0.0_wp,dtime,x,y,-100.0_wp,-100.0_wp,0.0_wp,&
          deltax(ng)%d,deltay(ng)%d,1,1)
     call nbrex2d(phi_n,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
          nbrbottom,comm2d)
     tolphi=maxval(phi_n)*tol
     if (tolphi==0.0_wp) then
        tolphi=1e-3;
     endif
     call MPI_ALLREDUCE(tolphi,tolphi_tot,1,MPI_DOUBLE_PRECISION,&
          MPI_MAX,comm2d,ierr)
     errphi=maxval(abs(phi_n-phi))
     call MPI_ALLREDUCE(errphi,errphi_tot,1,MPI_DOUBLE_PRECISION,&
          MPI_MAX,comm2d,ierr)     
     if (errphi_tot.lt.tolphi_tot) exit
     phi=phi_n
  enddo
  phi=phi_n
  ngrid=ng
  do 
     if (ngrid.lt.ng) exit
     deallocate(ae(ngrid)%a,aw(ngrid)%a,an(ngrid)%a,as(ngrid)%a)
     ngrid=ngrid-1
  enddo
  deallocate(ae,aw,an,as,phi_n,rhs)
!
CONTAINS
  SUBROUTINE calc_ad(ae,aw,an,as,le,ng,nglow,dx,dy,mynx,myny,deltax,deltay)
    USE btype
    USE bfunc2d
    IMPLICIT NONE 
    TYPE(ptr2d), INTENT(INOUT) :: ae(:),aw(:),an(:),as(:)
    ! Le is the Lewis number, thermal diffusivity over chemical diffusivity
    REAL(WP), INTENT(IN) :: le
    TYPE(ptr1d), INTENT(IN) :: dx(:),dy(:),deltax(:),deltay(:)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny 
    INTEGER(I4B), INTENT(IN) :: ng,nglow
    INTEGER(I4B) :: i,j,ngrid
    REAL(WP) :: h,h2
    !
    ! Allocate space for coefficient arrays on finest grid 
    allocate(ae(ng)%a(mynx(ng)-2,myny(ng)-2),aw(ng)%a(mynx(ng)-2,myny(ng)-2),&
         an(ng)%a(mynx(ng)-2,myny(ng)-2),as(ng)%a(mynx(ng)-2,myny(ng)-2))
    ae(ng)%a=0.0_wp
    aw(ng)%a=0.0_wp
    an(ng)%a=0.0_wp
    as(ng)%a=0.0_wp
    ! Calculate coefficients on finest grid level
    forall (i=1:mynx(ng)-2,j=1:myny(ng)-2)
       aw(ng)%a(i,j)=(deltay(ng)%d(j)/le)/dx(ng)%d(i)
       ae(ng)%a(i,j)=(deltay(ng)%d(j)/le)/dx(ng)%d(i+1)
       as(ng)%a(i,j)=(deltax(ng)%d(i)/le)/dy(ng)%d(j)
       an(ng)%a(i,j)=(deltax(ng)%d(i)/le)/dy(ng)%d(j+1)
    end forall
    ! Now loop through, allocate and calculate coefficients on all grid levels
    ngrid=ng
    do 
       if (ngrid<=nglow) exit
       ngrid=ngrid-1
       h=1.0_wp/(2**ngrid)
       allocate(ae(ngrid)%a(mynx(ngrid)-2,myny(ngrid)-2),&
            aw(ngrid)%a(mynx(ngrid)-2,myny(ngrid)-2),&
            an(ngrid)%a(mynx(ngrid)-2,myny(ngrid)-2),&
            as(ngrid)%a(mynx(ngrid)-2,myny(ngrid)-2))
       ae(ngrid)%a=0.0_wp
       aw(ngrid)%a=0.0_wp
       an(ngrid)%a=0.0_wp
       as(ngrid)%a=0.0_wp
       ! Calculate coefficients on next coarsest grid level
       forall (i=1:mynx(ngrid)-2,j=1:myny(ngrid)-2)
          aw(ngrid)%a(i,j)=(deltay(ng)%d(j)/le)/dx(ngrid)%d(i)
          ae(ngrid)%a(i,j)=(deltay(ng)%d(j)/le)/dx(ngrid)%d(i+1)
          as(ngrid)%a(i,j)=(deltax(ng)%d(i)/le)/dy(ngrid)%d(j)
          an(ngrid)%a(i,j)=(deltax(ng)%d(i)/le)/dy(ngrid)%d(j+1)
       end forall 
    enddo
    !
  END SUBROUTINE calc_ad
END SUBROUTINE diffusion
