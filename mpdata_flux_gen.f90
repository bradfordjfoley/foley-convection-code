! 2d advection solver using the MPDATA scheme of Smolarkiewicz (1984)
!
! This solves for a pure advection problem; for an advection-diffusion 
! problem one can employ operator splitting, and solve for the diffusion 
! seperately with a pure diffusion step
SUBROUTINE mpdata_flux_gen(phi,u,w,dtime,ng,iord,mynx,myny,myid,nbrleft,&
     nbrright,nbrtop,nbrbottom,tfsl,bfsl,lfsl,rfsl,comm2d,tbdy,bbdy,tbc,bbc,&
     utfsl,ubfsl,wlfsl,wrfsl,deltax,deltay,dx,dy,x,y,xu,yw)
  USE btype
  USE bfunc2d
  IMPLICIT NONE
  REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: phi
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
  ! Deltax and deltay are control volume size, dx and dy are spacing between
  ! cell center grid points
  REAL(WP), DIMENSION(:), INTENT(IN) :: deltax,deltay,dx,dy,x,y,xu,yw
  INTEGER(I4B), INTENT(IN) :: iord,myid,nbrleft,nbrright,nbrtop,&
       nbrbottom,ng,comm2d
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
  REAL(WP), INTENT(IN) :: dtime
  LOGICAL, INTENT(IN) :: tfsl,bfsl,tbdy,bbdy,lfsl,rfsl,utfsl,ubfsl,wlfsl,wrfsl
  CHARACTER(LEN=3), INTENT(IN) :: tbc,bbc
  INTEGER(I4B) :: i,j,iter
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: ae,aw,an,as,phi_n,u_t,w_t,u1,w1,&
       beta_up_u,beta_down_u,beta_up_w,beta_down_w,u_prime,w_prime,phi_old
  !
  ! Allocate arrays for FV coefficients
  allocate(ae(mynx(ng)-2,myny(ng)-2),aw(mynx(ng)-2,myny(ng)-2),&
       an(mynx(ng)-2,myny(ng)-2),as(mynx(ng)-2,myny(ng)-2),&
       phi_n(mynx(ng),myny(ng)),u_t(size(u,1),size(u,2)),&
       w_t(size(w,1),size(w,2)),u1(size(u,1),size(u,2)),&
       w1(size(w,1),size(w,2)),beta_up_u(mynx(ng),myny(ng)),&
       beta_down_u(mynx(ng),myny(ng)),u_prime(size(u,1),size(u,2)),&
       w_prime(size(w,1),size(w,2)),beta_up_w(mynx(ng),myny(ng)),&
       beta_down_w(mynx(ng),myny(ng)),phi_old(mynx(ng),myny(ng)))
  ae=0.0_wp
  aw=0.0_wp
  an=0.0_wp
  as=0.0_wp
  phi_n=0.0_wp
  u_t=0.0_wp
  w_t=0.0_wp
  u1=u
  w1=w
  phi_old=phi
  beta_up_u=0.0_wp; beta_down_u=0.0_wp; beta_up_w=0.0_wp; beta_down_w=0.0_wp
  u_prime=0.0_wp; w_prime=0.0_wp
  iter=0
  do 
     ! Calculate the advective fluxes/coefficents for upwind advection 
     call calc_a(ae,aw,an,as,mynx,myny,ng,u1,w1,deltax,deltay)
     ! Now do upwind advection with current velocity field
     call upwind(phi,phi_n,dtime,ae,aw,an,as,ng,u1,w1,deltax,deltay)
     ! Neighbor exchange 
     call nbrex2d(phi_n,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
          nbrbottom,comm2d)  
     phi=phi_n
     iter=iter+1
     if (iter.ge.iord) exit
     ! If we have more corrective steps to do, then
     ! calcuate anti-diffusion velocities (u_t, w_t) 
     call calc_u_tilda(phi_n,u_t,w_t,dtime,ng,u1,w1,tbdy,bbdy,tbc,bbc,&
          dx,dy,x,y,xu,yw)
     ! Neighbor exchange for u_t 
     call nbrex2d(u_t,utfsl,ubfsl,.false.,.false.,myid,nbrleft,nbrright,&
          nbrtop,nbrbottom,comm2d)
     ! Neighbor exchange for w_t 
     call nbrex2d(w_t,.false.,.false.,wlfsl,wrfsl,myid,nbrleft,nbrright,&
          nbrtop,nbrbottom,comm2d)
     ! Calculate beta_up and beta_down
     call flux_lim(beta_up_u,beta_down_u,beta_up_w,beta_down_w,phi_n,phi_old,&
          u_t,w_t,dtime,ng,dx,dy)
     ! Neighbor exchange for beta_up and down
     call nbrex2d(beta_up_u,.false.,.false.,.false.,.false.,myid,&
          nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
     call nbrex2d(beta_down_u,.false.,.false.,.false.,.false.,myid,&
          nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
     call nbrex2d(beta_up_w,.false.,.false.,.false.,.false.,myid,&
          nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
     call nbrex2d(beta_down_w,.false.,.false.,.false.,.false.,myid,&
          nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
     ! Now calculate u_prime and w_prime
     call calc_uprime(u_prime,w_prime,beta_up_u,beta_down_u,beta_up_w,&
          beta_down_w,u_t,w_t)
     ! Neighbor exchange for u_prime and w_prime 
     call nbrex2d(u_prime,utfsl,ubfsl,.false.,.false.,myid,nbrleft,&
          nbrright,nbrtop,nbrbottom,comm2d) 
     call nbrex2d(w_prime,.false.,.false.,wlfsl,wrfsl,myid,nbrleft,&
          nbrright,nbrtop,nbrbottom,comm2d) 
     u1=u_prime
     w1=w_prime
  end do
  deallocate(ae,aw,as,an,phi_n,u_t,w_t,u1,w1,u_prime,w_prime,&
       beta_down_u,beta_up_u,beta_up_w,beta_down_w,phi_old)
  !
CONTAINS 
  !
  SUBROUTINE calc_a(ae,aw,an,as,mynx,myny,ng,u,w,deltax,deltay)
    USE btype
    USE bfunc2d
    IMPLICIT NONE
    REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: ae,aw,an,as
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
    REAL(WP), DIMENSION(:), INTENT(IN) :: deltax,deltay
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
    INTEGER(I4B), INTENT(IN) :: ng 
    INTEGER(I4B) :: i,j
    ! Calculate coefficients 
    ! Here deltax and deltay are control volume size
    forall (i=1:mynx(ng)-2,j=1:myny(ng)-2)
       ae(i,j)=max(-u(i+1,j+1),0.0_wp)*deltay(j)
       aw(i,j)=max(u(i,j+1),0.0_wp)*deltay(j)
       an(i,j)=max(-w(i+1,j+1),0.0_wp)*deltax(i) 
       as(i,j)=max(w(i+1,j),0.0_wp)*deltax(i)
    endforall
  END SUBROUTINE calc_a
!
  SUBROUTINE upwind(phi,phi_n,dt,ae,aw,an,as,ng,u,w,deltax,deltay)
    USE btype
    USE bfunc2d
    IMPLICIT NONE
    REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: ae,aw,an,as
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w,phi
    REAL(WP), DIMENSION(:,:), INTENT(OUT) :: phi_n
    REAL(WP), DIMENSION(:), INTENT(IN) :: deltax,deltay
    REAL(WP), INTENT(IN) :: dt
    INTEGER(I4B), INTENT(IN) :: ng 
    INTEGER(I4B) :: i,j,nx1,ny1
    ! 
    ! Here deltax and deltay are control volume size
    nx1=size(phi,1)
    ny1=size(phi,2)
    ! Now do explicit time step upwind scheme
    do j=2,ny1-1
       do i=2,nx1-1
          phi_n(i,j)=(dt/(deltax(i-1)*deltay(j-1)))*&
               ((((deltax(i-1)*deltay(j-1))/dt)-(ae(i-1,j-1)+aw(i-1,j-1)+&
               as(i-1,j-1)+an(i-1,j-1))-(u(i,j)-u(i-1,j))*deltay(j-1)-&
               (w(i,j)-w(i,j-1))*deltax(i-1))*phi(i,j)+&
               aw(i-1,j-1)*phi(i-1,j)+ae(i-1,j-1)*phi(i+1,j)+&
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
  SUBROUTINE calc_u_tilda(phi_n,u_t,w_t,dt,ng,u,w,tbdy,bbdy,tbc,bbc,dx,&
       dy,x,y,xu,yw)
    USE btype
    USE bfunc2d
    IMPLICIT NONE
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: phi_n
    REAL(WP), DIMENSION(:,:), INTENT(OUT) :: u_t,w_t
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
    REAL(WP), DIMENSION(:), INTENT(IN) :: dx,dy,x,y,xu,yw
    LOGICAL, INTENT(IN) :: tbdy,bbdy
    CHARACTER(LEN=3), INTENT(IN) :: tbc,bbc
    REAL(WP), INTENT(IN) :: dt
    INTEGER(I4B), INTENT(IN) :: ng 
    INTEGER(I4B) :: i,j,nxu,nyu,nxw,nyw
    REAL(WP) :: h,eps,w_avg,phi_deriv_e,phi_deriv_w,phi_avg,&
         u_avg,phi_deriv_n,phi_deriv_s
    !
    nxu=size(u,1)
    nyu=size(u,2)
    nxw=size(w,1)
    nyw=size(w,2)
    eps=epsilon(0.0_wp)
    ! Horizontal velocity
    ! dx & dy are spacing between cell center grid points (T,P,etc)
    do j=2,nyu-1
       do i=2,nxu-1
          w_avg=((yw(j)-y(j))/(yw(j)-yw(j-1)))*&
               ((x(i+1)-xu(i))/(x(i+1)-x(i)))*w(i,j-1)+&
               ((yw(j)-y(j))/(yw(j)-yw(j-1)))*&
               ((xu(i)-x(i))/(x(i+1)-x(i)))*w(i+1,j-1)+&
               ((y(j)-yw(j-1))/(yw(j)-yw(j-1)))*&
               ((x(i+1)-xu(i))/(x(i+1)-x(i)))*w(i,j)+&
               ((y(j)-yw(j-1))/(yw(j)-yw(j-1)))*&
               ((xu(i)-x(i))/(x(i+1)-x(i)))*w(i+1,j)
          phi_avg=((y(j+1)-y(j))/(y(j+1)-y(j-1)))*&
               ((x(i+1)-xu(i))/(x(i+1)-x(i)))*phi_n(i,j-1)+&
               ((y(j+1)-y(j))/(y(j+1)-y(j-1)))*&
               ((xu(i)-x(i))/(x(i+1)-x(i)))*phi_n(i+1,j-1)+&
               ((y(j)-y(j-1))/(y(j+1)-y(j-1)))*&
               ((x(i+1)-xu(i))/(x(i+1)-x(i)))*phi_n(i,j+1)+&
               ((y(j)-y(j-1))/(y(j+1)-y(j-1)))*&
               ((xu(i)-x(i))/(x(i+1)-x(i)))*phi_n(i+1,j+1)
          phi_deriv_e=(phi_n(i+1,j)-phi_n(i+1,j-1))/dy(j-1)+&
               ((phi_n(i+1,j+1)-phi_n(i+1,j))/dy(j)-&
               (phi_n(i+1,j)-phi_n(i+1,j-1))/dy(j-1))*&
               (dy(j-1)/(dy(j)+dy(j-1)))
          phi_deriv_w=(phi_n(i,j)-phi_n(i,j-1))/dy(j-1)+&
               ((phi_n(i,j+1)-phi_n(i,j))/dy(j)-&
               (phi_n(i,j)-phi_n(i,j-1))/dy(j-1))*&
               (dy(j-1)/(dy(j)+dy(j-1)))
          u_t(i,j)=0.5_wp*(abs(u(i,j))*dx(i)-dt*u(i,j)**2)*&
               ((phi_n(i+1,j)-phi_n(i,j))/dx(i))*(phi_n(i,j)+&
               (phi_n(i+1,j)-phi_n(i,j))*((xu(i)-x(i))/&
               (dx(i)))+eps)**(-1.0_wp)-0.5_wp*dt*u(i,j)*w_avg*&
               (phi_deriv_w+(xu(i)-x(i))*((phi_deriv_e-phi_deriv_w)/x(i)))*&
               (phi_n(i,j)+(xu(i)-x(i))*((phi_n(i+1,j)-phi_n(i,j))/&
               deltax(i))+eps)**(-1.0_wp)
       enddo
    enddo
!!$    forall (i=2:nxu-1,j=3:nyu-2)
!!$       u_t(i,j)=0.5_wp*(abs(u(i,j))*deltax(i)-dt*u(i,j)**2)*&
!!$            ((phi_n(i+1,j)-phi_n(i,j))/deltax(i))*(phi_n(i,j)+&
!!$            (phi_n(i+1,j)-phi_n(i,j))*((xu(i)-x(i))/(deltax(i)))+eps)**(-1.0_wp)-&
!!$            0.5_wp*dt*u(i,j)*(((yw(j)-y(j))/(yw(j)-yw(j-1)))*&
!!$            ((x(i+1)-xu(i))/(x(i+1)-x(i)))*w(i,j-1)+&
!!$            ((yw(j)-y(j))/(yw(j)-yw(j-1)))*((xu(i)-x(i))/(x(i+1)-x(i)))*w(i+1,j-1)+&
!!$            ((y(j)-yw(j-1))/(yw(j)-yw(j-1)))*((x(i+1)-xu(i))/(x(i+1)-x(i)))*w(i,j)+&
!!$            ((y(j)-yw(j-1))/(yw(j)-yw(j-1)))*((xu(i)-x(i))/(x(i+1)-x(i)))*w(i+1,j))*&
!!$            ((phi_n(i+1,j+1)+phi_n(i,j+1)-phi_n(i+1,j-1)-phi_n(i,j-1))/&
!!$            ((phi_n(i+1,j+1)+phi_n(i,j+1)+phi_n(i+1,j-1)+phi_n(i,j-1)+eps)*h))
!!$    end forall
!!$    ! Do top and bottom differently if tbc/bbc is constant value
!!$    if (tbdy.and.tbc=='con') then
!!$       j=nyu-1
!!$       do i=2,nxu-1 
!!$          u_t(i,j)=(abs(u(i,j))*h-dt*u(i,j)**2)*((phi_n(i+1,j)-phi_n(i,j))/&
!!$               ((phi_n(i+1,j)+phi_n(i,j)+eps)*h))-1.5_wp*dt*u(i,j)*&
!!$               0.25_wp*(w(i+1,j)+w(i,j)+w(i+1,j-1)+w(i,j-1))*&
!!$               ((2.0_wp*phi_n(i+1,j+1)+2.0_wp*phi_n(i,j+1)-phi_n(i,j)-&
!!$               phi_n(i+1,j)-phi_n(i+1,j-1)-phi_n(i,j-1))/&
!!$               ((2.0_wp*phi_n(i+1,j+1)+2.0_wp*phi_n(i,j+1)+&
!!$               phi_n(i+1,j-1)+phi_n(i,j-1)+eps)*h))
!!$       enddo
!!$    else
!!$       j=nyu-1
!!$       do i=2,nxu-1 
!!$          u_t(i,j)=(abs(u(i,j))*h-dt*u(i,j)**2)*((phi_n(i+1,j)-phi_n(i,j))/&
!!$               ((phi_n(i+1,j)+phi_n(i,j)+eps)*h))-0.5_wp*dt*u(i,j)*&
!!$               0.25_wp*(w(i+1,j)+w(i,j)+w(i+1,j-1)+w(i,j-1))*&
!!$               ((phi_n(i+1,j+1)+phi_n(i,j+1)-phi_n(i+1,j-1)-phi_n(i,j-1))/&
!!$               ((phi_n(i+1,j+1)+phi_n(i,j+1)+phi_n(i+1,j-1)+phi_n(i,j-1)+&
!!$               eps)*h))
!!$       enddo
!!$    endif
!!$    if (bbdy.and.bbc=='con') then
!!$       j=2
!!$       do i=2,nxu-1 
!!$          u_t(i,j)=(abs(u(i,j))*h-dt*u(i,j)**2)*((phi_n(i+1,j)-phi_n(i,j))/&
!!$               ((phi_n(i+1,j)+phi_n(i,j)+eps)*h))-1.5_wp*dt*u(i,j)*&
!!$               0.25_wp*(w(i+1,j)+w(i,j)+w(i+1,j-1)+w(i,j-1))*&
!!$               ((phi_n(i+1,j+1)+phi_n(i,j+1)+phi_n(i,j)+phi_n(i+1,j)-&
!!$               2.0_wp*phi_n(i+1,j-1)-2.0_wp*phi_n(i,j-1))/((phi_n(i+1,j+1)+&
!!$               phi_n(i,j+1)+2.0_wp*phi_n(i+1,j-1)+2.0_wp*phi_n(i,j-1)+eps)*h))
!!$       enddo
!!$    else
!!$       j=2
!!$       do i=2,nxu-1  
!!$          u_t(i,j)=(abs(u(i,j))*h-dt*u(i,j)**2)*((phi_n(i+1,j)-phi_n(i,j))/&
!!$               ((phi_n(i+1,j)+phi_n(i,j)+eps)*h))-0.5_wp*dt*u(i,j)*&
!!$               0.25_wp*(w(i+1,j)+w(i,j)+w(i+1,j-1)+w(i,j-1))*&
!!$               ((phi_n(i+1,j+1)+phi_n(i,j+1)-phi_n(i+1,j-1)-phi_n(i,j-1))/&
!!$               ((phi_n(i+1,j+1)+phi_n(i,j+1)+phi_n(i+1,j-1)+phi_n(i,j-1)+&
!!$               eps)*h))
!!$       enddo
!!$    endif
    ! Boundary conditions 
    u_t(1:nxu,1)=u(1:nxu,1)
    u_t(1:nxu,nyu)=u(1:nxu,nyu)
    u_t(1,1:nyu)=u(1,1:nyu)
    u_t(nxu,1:nyu)=u(nxu,1:nyu)
    ! Vertical velocity
    do j=2,nyw-1
       do i=2,nxw-1
          u_avg=((y(j+1)-yw(j))/(y(j+1)-y(j)))*&
               ((xu(i)-x(i))/(xu(i)-xu(i-1)))*u(i-1,j)+&
               ((y(j+1)-yw(j))/(y(j+1)-y(j)))*&
               ((x(i)-xu(i-1))/(xu(i)-xu(i-1)))*u(i,j)+&
               ((yw(j)-y(j))/(y(j+1)-y(j)))*&
               ((xu(i)-x(i))/(xu(i)-xu(i-1)))*u(i-1,j+1)+&
               ((yw(j)-y(j))/(y(j+1)-y(j)))*&
               ((x(i)-xu(i-1))/(xu(i)-xu(i-1)))*u(i,j+1)
          phi_avg=((y(j+1)-yw(j))/(y(j+1)-y(j)))*&
               ((x(i+1)-x(i))/(x(i+1)-x(i-1)))*phi_n(i-1,j)+&
               ((y(j+1)-yw(j))/(y(j+1)-y(j)))*&
               ((x(i)-x(i-1))/(x(i+1)-x(i-1)))*phi_n(i+1,j)+&
               ((yw(j)-y(j))/(y(j+1)-y(j)))*&
               ((x(i+1)-x(i))/(x(i+1)-x(i-1)))*phi_n(i-1,j+1)+&
               ((yw(j)-y(j))/(y(j+1)-y(j)))*&
               ((x(i)-x(i-1))/(x(i+1)-x(i)-1))*phi_n(i+1,j+1)
          phi_deriv_s=(phi_n(i,j)-phi_n(i-1,j))/dx(i-1)+&
               ((phi_n(i+1,j)-phi_n(i,j))/dx(i)-&
               (phi_n(i,j)-phi_n(i-1,j))/dx(i-1))*&
               (dx(i-1)/(dx(i)+dx(i-1)))
          phi_deriv_n=(phi_n(i,j+1)-phi_n(i-1,j+1))/dx(i-1)+&
               ((phi_n(i+1,j+1)-phi_n(i,j+1))/dx(i)-&
               (phi_n(i,j+1)-phi_n(i-1,j+1))/dx(i-1))*&
               (dx(i-1)/(dx(i)+dx(i-1)))
          w_t(i,j)=0.5_wp*(abs(w(i,j))*dy(j)-dt*w(i,j)**2)*&
               ((phi_n(i,j+1)-phi_n(i,j))/dy(i))*(phi_n(i,j)+&
               (phi_n(i,j+1)-phi_n(i,j))*((yw(j)-y(j))/&
               (dy(j)))+eps)**(-1.0_wp)-0.5_wp*dt*w(i,j)*u_avg*&
               (phi_deriv_s+(yw(j)-y(j))*((phi_deriv_n-phi_deriv_s)/dy(j)))*&
               (phi_n(i,j)+(yw(j)-y(j))*((phi_n(i,j+1)-phi_n(i,j))/&
               dy(j))+eps)**(-1.0_wp)
       enddo
    enddo
!!$    forall (i=2:nxw-1,j=2:nyw-1)
!!$       w_t(i,j)=(abs(w(i,j))*h-dt*w(i,j)**2)*((phi_n(i,j+1)-phi_n(i,j))/&
!!$            ((phi_n(i,j+1)+phi_n(i,j)+eps)*h))-0.5_wp*dt*w(i,j)*&
!!$            0.25_wp*(u(i,j+1)+u(i,j)+u(i-1,j+1)+u(i-1,j))*&
!!$            ((phi_n(i+1,j+1)+phi_n(i+1,j)-phi_n(i-1,j+1)-phi_n(i-1,j))/&
!!$            ((phi_n(i+1,j+1)+phi_n(i+1,j)+phi_n(i-1,j+1)+phi_n(i-1,j)+eps)*h))
!!$    end forall
    ! Boundary conditions 
    w_t(1:nxw,1)=w(1:nxw,1)
    w_t(1:nxw,nyw)=w(1:nxw,nyw)
    w_t(1,1:nyw)=w(1,1:nyw)
    w_t(nxw,1:nyw)=w(nxw,1:nyw)
  END SUBROUTINE calc_u_tilda
!
  SUBROUTINE flux_lim(beta_up_u,beta_down_u,beta_up_w,beta_down_w,phi_n,phi,&
       u_t,w_t,dt,ng,dx,dy)
    USE btype
    USE bfunc2d
    IMPLICIT NONE
    REAL(WP), DIMENSION(:,:), INTENT(OUT) :: beta_up_u,beta_down_u,beta_up_w,&
         beta_down_w
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: phi_n,phi,u_t,w_t 
    REAL(WP), DIMENSION(:), INTENT(IN) :: dx,dy
    REAL(WP),INTENT(IN) :: dt
    INTEGER(I4B), INTENT(IN) :: ng
    INTEGER(I4B) :: nx1,ny1
    REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: phi_n_max_u,phi_n_min_u,&
         phi_n_max_w,phi_n_min_w
    REAL(WP) :: h,eps
    !
    nx1=size(beta_up_u,1)
    ny1=size(beta_up_u,2)
    h=1.0_wp/(2.0_wp**ng)
    eps=epsilon(0.0_wp)
    allocate(phi_n_max_u(nx1,ny1),phi_n_min_u(nx1,ny1),phi_n_min_w(nx1,ny1),&
         phi_n_max_w(nx1,ny1))
    phi_n_max_u=0.0_wp; phi_n_min_u=0.0_wp; phi_n_max_w=0.0_wp 
    phi_n_min_w=0.0_wp 
    ! dx and dy are differences/spacing between cell centers

    ! Calculate phi_n_max and phi_n_min
    do j=2,ny1-1
       do i=2,nx1-1
          phi_n_max_u(i,j)=max(phi_n(i,j),phi_n(i-1,j),phi_n(i+1,j),&
               phi(i,j),phi(i-1,j),phi(i+1,j))
          phi_n_min_u(i,j)=min(phi_n(i,j),phi_n(i-1,j),phi_n(i+1,j),&
               phi(i,j),phi(i-1,j),phi(i+1,j))
          phi_n_max_w(i,j)=max(phi_n(i,j),phi_n(i,j-1),phi_n(i,j+1),&
               phi(i,j),phi(i,j-1),phi(i,j+1))
          phi_n_min_w(i,j)=min(phi_n(i,j),phi_n(i,j-1),phi_n(i,j+1),&
               phi(i,j),phi(i,j-1),phi(i,j+1))
          beta_up_u(i,j)=(phi_n_max_u(i,j)-phi_n(i,j))/((dt/dx(i))*&
               (max(u_t(i-1,j),0.0_wp)*phi_n(i-1,j)-min(u_t(i,j),0.0_wp)*&
               phi_n(i+1,j))+(dt/dy(j))*(max(w_t(i,j-1),0.0_wp)*phi_n(i,j-1)-&
               min(w_t(i,j),0.0_wp)*phi_n(i,j+1))+eps)
          beta_down_u(i,j)=(phi_n(i,j)-phi_n_min_u(i,j))/((dt/dx(i))*&
               (max(u_t(i,j),0.0_wp)*phi_n(i,j)-min(u_t(i-1,j),0.0_wp)*&
               phi_n(i,j))+(dt/dy(j))*(max(w_t(i,j),0.0_wp)*phi_n(i,j)-&
               min(w_t(i,j-1),0.0_wp)*phi_n(i,j))+eps)
          beta_up_w(i,j)=(phi_n_max_w(i,j)-phi_n(i,j))/((dt/dx(i))*&
               (max(u_t(i-1,j),0.0_wp)*phi_n(i-1,j)-min(u_t(i,j),0.0_wp)*&
               phi_n(i+1,j))+(dt/dy(j))*(max(w_t(i,j-1),0.0_wp)*phi_n(i,j-1)-&
               min(w_t(i,j),0.0_wp)*phi_n(i,j+1))+eps)
          beta_down_w(i,j)=(phi_n(i,j)-phi_n_min_w(i,j))/((dt/dx(i))*&
               (max(u_t(i,j),0.0_wp)*phi_n(i,j)-min(u_t(i-1,j),0.0_wp)*&
               phi_n(i,j))+(dt/dy(j))*(max(w_t(i,j),0.0_wp)*phi_n(i,j)-&
               min(w_t(i,j-1),0.0_wp)*phi_n(i,j))+eps)
       enddo
    enddo
    ! Deallocate phi_n_min,phi_n_max
    deallocate(phi_n_min_u,phi_n_max_u,phi_n_min_w,phi_n_max_w)
  END SUBROUTINE flux_lim
  !
  SUBROUTINE calc_uprime(u_prime,w_prime,beta_up_u,beta_down_u,&
       beta_up_w,beta_down_w,u_t,w_t)
    USE btype
    USE bfunc2d
    IMPLICIT NONE
    REAL(WP), DIMENSION(:,:), INTENT(OUT) :: u_prime,w_prime
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: beta_up_u,beta_down_u,beta_up_w,&
         beta_down_w,u_t,w_t
    INTEGER(I4B) :: nxu,nyu,nxw,nyw
    !
    nxu=size(u,1)
    nyu=size(u,2)
    nxw=size(w,1)
    nyw=size(w,2)
    ! U_prime
    forall (i=2:nxu-1,j=2:nyu-1)
       u_prime(i,j)=min(1.0_wp,beta_down_u(i,j),beta_up_u(i+1,j))*&
            max(u_t(i,j),0.0_wp)+min(1.0_wp,beta_up_u(i,j),&
            beta_down_u(i+1,j))*min(u_t(i,j),0.0_wp)
    end forall
    ! Boundary points
    u_prime(1:nxu,1)=u_t(1:nxu,1)
    u_prime(1:nxu,nyu)=u_t(1:nxu,nyu)
    u_prime(1,1:nyu)=u_t(1,1:nyu)
    u_prime(nxu,1:nyu)=u_t(nxu,1:nyu)   
    ! W_prime
    forall (i=2:nxw-1,j=2:nyw-1)
       w_prime(i,j)=min(1.0_wp,beta_down_w(i,j),beta_up_w(i,j+1))*&
            max(w_t(i,j),0.0_wp)+min(1.0_wp,beta_up_w(i,j),&
            beta_down_w(i,j+1))*min(w_t(i,j),0.0_wp)
    end forall
    ! Boundary points
    w_prime(1:nxw,1)=w_t(1:nxw,1)
    w_prime(1:nxw,nyw)=w_t(1:nxw,nyw)
    w_prime(1,1:nyw)=w_t(1,1:nyw)
    w_prime(nxw,1:nyw)=w_t(nxw,1:nyw) 
  END SUBROUTINE calc_uprime
  !
END SUBROUTINE mpdata_flux_gen
