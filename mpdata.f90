! 2d advection solver using the MPDATA scheme of Smolarkiewicz (1984)
!
! This solves for a pure advection problem; for an advection-diffusion 
! problem one can employ operator splitting, and solve for the diffusion 
! seperately with a pure diffusion step
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
  INTEGER(I4B) :: i,j,iter
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: ae,aw,an,as,phi_n,u_t,w_t,u1,w1
  !
  ! Allocate arrays for FV coefficients
  allocate(ae(mynx(ng)-2,myny(ng)-2),aw(mynx(ng)-2,myny(ng)-2),&
       an(mynx(ng)-2,myny(ng)-2),as(mynx(ng)-2,myny(ng)-2),&
       phi_n(mynx(ng),myny(ng)),u_t(size(u,1),size(u,2)),&
       w_t(size(w,1),size(w,2)),u1(size(u,1),size(u,2)),&
       w1(size(w,1),size(w,2)))
  ae=0.0_wp
  aw=0.0_wp
  an=0.0_wp
  as=0.0_wp
  phi_n=0.0_wp
  u_t=0.0_wp
  w_t=0.0_wp
  u1=u
  w1=w
  iter=0
  do 
     ! Calculate the advective fluxes/coefficents for upwind advection 
     call calc_a(ae,aw,an,as,mynx,myny,ng,u1,w1)
     ! Now do upwind advection with current velocity field
     call upwind(phi,phi_n,dtime,ae,aw,an,as,ng,u1,w1)
     ! Neighbor exchange 
     call nbrex2d(phi_n,tfsl,bfsl,myid,nbrleft,nbrright,nbrtop,&
          nbrbottom,comm2d)  
     phi=phi_n
     iter=iter+1
     if (iter.ge.iord) exit
     ! If we have more corrective steps to do, then
     ! calcuate anti-diffusion velocities (u_t, w_t) 
     call calc_u_tilda(phi_n,u_t,w_t,dtime,ng,u1,w1)
     ! Neighbor exchange for u_t 
     call nbrex2d(u_t,.false.,.false.,myid,nbrleft,nbrright,nbrtop,&
          nbrbottom,comm2d)
     !Neighbor exchange for w_t 
     call nbrex2d(w_t,.false.,.false.,myid,nbrleft,nbrright,nbrtop,&
          nbrbottom,comm2d)
     u1=u_t
     w1=w_t
  end do
  deallocate(ae,aw,as,an,phi_n)
  !
CONTAINS 
  !
  SUBROUTINE calc_a(ae,aw,an,as,mynx,myny,ng,u,w)
    USE btype
    USE bfunc2d
    IMPLICIT NONE
    REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: ae,aw,an,as
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
    INTEGER(I4B), INTENT(IN) :: ng 
    INTEGER(I4B) :: i,j
    REAL(WP) :: h
    !
    h=1.0_wp/(2.0_wp**ng)
    ! Calculate coefficients (assume Dx = Dz) and constant control 
    ! volume size
    forall (i=1:mynx(ng)-2,j=1:myny(ng)-2)
       ae(i,j)=max(-u(i+1,j+1),0.0_wp)*h
       aw(i,j)=max(u(i,j+1),0.0_wp)*h
       an(i,j)=max(-w(i+1,j+1),0.0_wp)*h 
       as(i,j)=max(w(i+1,j),0.0_wp)*h
    endforall
  END SUBROUTINE calc_a
!
  SUBROUTINE upwind(phi,phi_n,dt,ae,aw,an,as,ng,u,w)
    USE btype
    USE bfunc2d
    IMPLICIT NONE
    REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: ae,aw,an,as
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w,phi
    REAL(WP), DIMENSION(:,:), INTENT(OUT) :: phi_n
    REAL(WP), INTENT(IN) :: dt
    INTEGER(I4B), INTENT(IN) :: ng 
    INTEGER(I4B) :: i,j,nx1,ny1
    REAL(WP) :: h,h2,ap1
    !
    h=1.0_wp/(2.0_wp**ng)
    h2=h*h
    ap1=h2/dt
    ! 
    nx1=size(phi,1)
    ny1=size(phi,2)
    ! Now do explicit time step upwind scheme
    do j=2,ny1-1
       do i=2,nx1-1
          phi_n(i,j)=(1.0_wp/ap1)*((ap1-(ae(i-1,j-1)+aw(i-1,j-1)+as(i-1,j-1)&
               +an(i-1,j-1))-(u(i,j)-u(i-1,j))*h-(w(i,j)-w(i,j-1))*h)*&
               phi(i,j)+aw(i-1,j-1)*phi(i-1,j)+ae(i-1,j-1)*phi(i+1,j)+&
               as(i-1,j-1)*phi(i,j-1)+an(i-1,j-1)*phi(i,j+1))
       enddo
    enddo
    ! Boundary conditions
    ! Copy all boundary grid points from old solution to new solution
    ! This will preserve boundary conditions, and will be written over 
    ! after communication step for ghost points
    phi_n(1:nx1,1)=phi(1:nx1,1)
    phi_n(1:nx1,ny1)=phi(1:nx1,ny1)
    phi_n(1,1:ny1)=phi(1,1:ny1)
    phi_n(nx1,1:ny1)=phi(nx1,1:ny1)
  END SUBROUTINE upwind
!
  SUBROUTINE calc_u_tilda(phi_n,u_t,w_t,dt,ng,u,w)
    USE btype
    USE bfunc2d
    IMPLICIT NONE
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: phi_n
    REAL(WP), DIMENSION(:,:), INTENT(OUT) :: u_t,w_t
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
    REAL(WP) :: dt
    INTEGER(I4B), INTENT(IN) :: ng 
    INTEGER(I4B) :: i,j,nxu,nyu,nxw,nyw
    REAL(WP) :: h,eps
    !
    h=1.0_wp/(2.0_wp**ng)
    nxu=size(u,1)
    nyu=size(u,2)
    nxw=size(w,1)
    nyw=size(w,2)
    eps=epsilon(0.0_wp)
    ! Horizontal velocity
    forall (i=2:nxu-1,j=2:nyu-1)
       u_t(i,j)=(abs(u(i,j))*h-dt*u(i,j)**2)*((phi_n(i+1,j)-phi_n(i,j))/&
            ((phi_n(i+1,j)+phi_n(i,j)+eps)*h))-0.5_wp*dt*u(i,j)*&
            0.25_wp*(w(i+1,j)+w(i,j)+w(i+1,j-1)+w(i,j-1))*&
            ((phi_n(i+1,j+1)+phi_n(i,j+1)-phi_n(i+1,j-1)-phi_n(i,j-1))/&
            ((phi_n(i+1,j+1)+phi_n(i,j+1)+phi_n(i+1,j-1)+phi_n(i,j-1)+eps)*h))
    end forall
    ! Boundary conditions 
    u_t(1:nxu,1)=u(1:nxu,1)
    u_t(1:nxu,nyu)=u(1:nxu,nyu)
    u_t(1,1:nyu)=u(1,1:nyu)
    u_t(nxu,1:nyu)=u(nxu,1:nyu)
    ! Vertical velocity
    forall (i=2:nxw-1,j=2:nyw-1)
       w_t(i,j)=(abs(w(i,j))*h-dt*w(i,j)**2)*((phi_n(i,j+1)-phi_n(i,j))/&
            ((phi_n(i,j+1)+phi_n(i,j)+eps)*h))-0.5_wp*dt*w(i,j)*&
            0.25_wp*(u(i,j+1)+u(i,j)+u(i-1,j+1)+u(i-1,j))*&
            ((phi_n(i+1,j+1)+phi_n(i+1,j)-phi_n(i-1,j+1)-phi_n(i-1,j))/&
            ((phi_n(i+1,j+1)+phi_n(i+1,j)+phi_n(i-1,j+1)+phi_n(i-1,j)+eps)*h))
    end forall
    ! Boundary conditions 
    w_t(1:nxw,1)=w(1:nxw,1)
    w_t(1:nxw,nyw)=w(1:nxw,nyw)
    w_t(1,1:nyw)=w(1,1:nyw)
    w_t(nxw,1:nyw)=w(nxw,1:nyw)
  END SUBROUTINE calc_u_tilda
!
END SUBROUTINE mpdata
