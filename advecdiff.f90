! 2-D advection diffusion solver, following Patankar 1980
! Use power-law upwind scheme for advection, central difference for diffusion 
! and fully implicit time-stepping 
SUBROUTINE advecdiff(phi,lewis,dtime,u,w,ng,nglow,x,y,dx,dy,mynx,myny,myid,&
     nbrleft,nbrright,nbrtop,nbrbottom,tfsl,bfsl,comm2d)
  USE btype
  USE bfunc2d 
  IMPLICIT NONE 
  REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: phi
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
  TYPE(ptr1d), INTENT(IN) :: x(:),y(:),dx(:),dy(:)
  REAL(WP), INTENT(IN) :: lewis,dtime 
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
  INTEGER(I4B), INTENT(IN) :: nbrleft,nbrright,nbrtop,nbrbottom,comm2d,ng,&
       nglow,myid
  LOGICAL, INTENT(IN) :: tfsl,bfsl
  TYPE(ptr2d), ALLOCATABLE :: ae(:),aw(:),an(:),as(:) 
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: phi_n,rhs, phi_last 
  INTEGER(I4B) :: i,j,n,nfinal,ngrid
  REAL(WP) :: h,h2,tol,tolphi,errphi,ap0
! 
  ! Allocate variables needed for this subroutine 
  allocate(ae(ng),aw(ng),an(ng),as(ng),rhs(mynx(ng),myny(ng)))
  ! Calculate finite volume coefficients
  call calc_ad(ae,aw,an,as,lewis,ng,ng,dx,dy,mynx,myny,u,w)
  ! Build rhs using solution previous time step 
  h=1.0_wp/(2.0_wp**ng)
  h2=h*h
  ap0=h2/dtime
!!$  ! Convert temperature perturbation to full temperature
!!$  forall (i=1:mynx(ng),j=1:myny(ng))
!!$     phi(i,j)=phi(i,j) + (1.0_wp-y(ng)%d(j))
!!$  endforall
  ! Set up right hand side
  rhs=(-1.0_wp/dtime)*phi
  
  ! Now call multi-grid to solve for diffusion
!!$  call mg2d(phi,rhs,dtime,ng,nglow,0,2,mynx,myny,x,y,ae,aw,an,as,myid,&
!!$       nbrleft,nbrright,nbrtop,nbrbottom,tfsl,bfsl,comm2d)
  ! Check if solution has converged 

  ! Solve via Jacobi iterations (only takes a few since time step is small 
  ! and initial guess from old time step is already close to the solution at 
  ! next time step)
  allocate(phi_n(mynx(ng),myny(ng)),phi_last(mynx(ng),myny(ng)))
  phi_n=0.0_wp; phi_last=0.0_wp
  n=0
  nfinal=1000
  tol=1e-10
  do while (n.lt.nfinal) 
     n=n+1
     phi_last=phi
     call jacobi2dfvp2(phi,rhs,ng,phi_n,ae(ng)%a,aw(ng)%a,an(ng)%a,&
          as(ng)%a,ap0,0.0_wp)
     call nbrex2d(phi_n,tfsl,bfsl,myid,nbrleft,nbrright,nbrtop,&
          nbrbottom,comm2d)
     call jacobi2dfvp2(phi_n,rhs,ng,phi,ae(ng)%a,aw(ng)%a,an(ng)%a,&
          as(ng)%a,ap0,0.0_wp)
     call nbrex2d(phi,tfsl,bfsl,myid,nbrleft,nbrright,nbrtop,&
          nbrbottom,comm2d)
     tolphi=maxval(phi)*tol
     errphi=maxval(abs(phi_last-phi))
     if (errphi.lt.tolphi) exit
  enddo
  
!!$  ! Convert full temperature back into temperature perturbation
!!$  forall (i=1:mynx(ng),j=1:myny(ng))
!!$     phi(i,j)=phi(i,j) - (1.0_wp-y(ng)%d(j))
!!$  endforall

  ! Deallocate arrays used in advec/diff solver
  ngrid=ng
  do 
     if (ngrid.lt.ng) exit
     deallocate(ae(ngrid)%a,aw(ngrid)%a,an(ngrid)%a,as(ngrid)%a)
     ngrid=ngrid-1
  enddo
  deallocate(ae,aw,an,as,phi_n,phi_last,rhs)
!
CONTAINS
  SUBROUTINE calc_ad(ae,aw,an,as,le,ng,nglow,dx,dy,mynx,myny,u,w)
    USE btype
    USE bfunc2d
    IMPLICIT NONE 
    TYPE(ptr2d), INTENT(INOUT) :: ae(:),aw(:),an(:),as(:)
    ! Le is the Lewis number, thermal diffusivity over chemical diffusivity
    REAL(WP), INTENT(IN) :: le
    TYPE(ptr1d), INTENT(IN) :: dx(:),dy(:)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny 
    INTEGER(I4B), INTENT(IN) :: ng,nglow
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
    INTEGER(I4B) :: i,j,ngrid
    REAL(WP) :: h,h2,gamma
    !
    ! Allocate space for coefficient arrays on finest grid 
    allocate(ae(ng)%a(mynx(ng)-2,myny(ng)-2),aw(ng)%a(mynx(ng)-2,myny(ng)-2),&
         an(ng)%a(mynx(ng)-2,myny(ng)-2),as(ng)%a(mynx(ng)-2,myny(ng)-2))
    ae(ng)%a=0.0_wp
    aw(ng)%a=0.0_wp
    an(ng)%a=0.0_wp
    as(ng)%a=0.0_wp
    h=1.0_wp/(2**ng)
    gamma=(1.0_wp/le)*h
    ! Calculate coefficients on finest grid level
    forall (i=1:mynx(ng)-2,j=1:myny(ng)-2)
       aw(ng)%a(i,j)=(gamma/dx(ng)%d(i))*max(0.0_wp,(1.0_wp-&
            0.1_wp*abs((h*u(i,j+1))/(gamma/dx(ng)%d(i))))**5.0_wp)+&
            max(u(i,j+1)*h,0.0_wp)
       ae(ng)%a(i,j)=(gamma/dx(ng)%d(i+1))*max(0.0_wp,(1.0_wp-&
            0.1_wp*abs((h*u(i+1,j+1))/(gamma/dx(ng)%d(i+1))))**5.0_wp)+&
            max(-u(i+1,j+1)*h,0.0_wp)
       as(ng)%a(i,j)=(gamma/dy(ng)%d(j))*max(0.0_wp,(1.0_wp-&
            0.1_wp*abs((h*w(i+1,j))/(gamma/dy(ng)%d(j))))**5.0_wp)+&
            max(w(i+1,j)*h,0.0_wp)
       an(ng)%a(i,j)=(gamma/dy(ng)%d(j+1))*max(0.0_wp,(1.0_wp-&
            0.1_wp*abs((h*w(i+1,j+1))/(gamma/dy(ng)%d(j+1))))**5.0_wp)+&
            max(-w(i+1,j+1)*h,0.0_wp)
    end forall
    ! Now loop through, allocate and calculate coefficients on all grid levels
    ! Need to do restriction on u/w to calculate 
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
          aw(ngrid)%a(i,j)=(h/le)/dx(ngrid)%d(i)
          ae(ngrid)%a(i,j)=(h/le)/dx(ngrid)%d(i+1)
          as(ngrid)%a(i,j)=(h/le)/dy(ngrid)%d(j)
          an(ngrid)%a(i,j)=(h/le)/dy(ngrid)%d(j+1)
       end forall 
    enddo
    !
  END SUBROUTINE calc_ad
END SUBROUTINE advecdiff
