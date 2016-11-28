! 2d advection solver using the MPDATA scheme of Smolarkiewicz (1984)
!
! This solves for a pure advection problem; for an advection-diffusion 
! problem one can employ operator splitting, and solve for the diffusion 
! seperately with a pure diffusion step
SUBROUTINE mpdata_flux(phi,u,w,dtime,ng,iord,mynx,myny,myid,nbrleft,&
     nbrright,nbrtop,nbrbottom,tfsl,bfsl,lfsl,rfsl,comm2d,tbdy,bbdy,tbc,bbc,&
     utfsl,ubfsl,wlfsl,wrfsl)
  USE btype
  USE bfunc2d
  IMPLICIT NONE
  REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: phi
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
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
     call calc_a(ae,aw,an,as,mynx,myny,ng,u1,w1)
     ! Now do upwind advection with current velocity field
     call upwind(phi,phi_n,dtime,ae,aw,an,as,ng,u1,w1)
     ! Neighbor exchange 
     call nbrex2d(phi_n,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,nbrtop,&
          nbrbottom,comm2d)  
     phi=phi_n
     iter=iter+1
     if (iter.ge.iord) exit
     ! If we have more corrective steps to do, then
     ! calcuate anti-diffusion velocities (u_t, w_t) 
     call calc_u_tilda(phi_n,u_t,w_t,dtime,ng,u1,w1,tbdy,bbdy,tbc,bbc)
     ! Neighbor exchange for u_t 
     call nbrex2d(u_t,utfsl,ubfsl,.false.,.false.,myid,nbrleft,nbrright,&
          nbrtop,nbrbottom,comm2d)
     ! Neighbor exchange for w_t 
     call nbrex2d(w_t,.false.,.false.,wlfsl,wrfsl,myid,nbrleft,nbrright,&
          nbrtop,nbrbottom,comm2d)
     ! Calculate beta_up and beta_down
     call flux_lim(beta_up_u,beta_down_u,beta_up_w,beta_down_w,phi_n,phi_old,&
          u_t,w_t,dtime,ng)
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
  SUBROUTINE calc_u_tilda(phi_n,u_t,w_t,dt,ng,u,w,tbdy,bbdy,tbc,bbc)
    USE btype
    USE bfunc2d
    IMPLICIT NONE
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: phi_n
    REAL(WP), DIMENSION(:,:), INTENT(OUT) :: u_t,w_t
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
    LOGICAL, INTENT(IN) :: tbdy,bbdy
    CHARACTER(LEN=3), INTENT(IN) :: tbc,bbc
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
    forall (i=2:nxu-1,j=3:nyu-2)
       u_t(i,j)=(abs(u(i,j))*h-dt*u(i,j)**2)*((phi_n(i+1,j)-phi_n(i,j))/&
            ((phi_n(i+1,j)+phi_n(i,j)+eps)*h))-0.5_wp*dt*u(i,j)*&
            0.25_wp*(w(i+1,j)+w(i,j)+w(i+1,j-1)+w(i,j-1))*&
            ((phi_n(i+1,j+1)+phi_n(i,j+1)-phi_n(i+1,j-1)-phi_n(i,j-1))/&
            ((phi_n(i+1,j+1)+phi_n(i,j+1)+phi_n(i+1,j-1)+phi_n(i,j-1)+eps)*h))
    end forall
    ! Do top and bottom differently if tbc/bbc is constant value
    if (tbdy.and.tbc=='con') then
       j=nyu-1
       do i=2,nxu-1 
          u_t(i,j)=(abs(u(i,j))*h-dt*u(i,j)**2)*((phi_n(i+1,j)-phi_n(i,j))/&
               ((phi_n(i+1,j)+phi_n(i,j)+eps)*h))-1.5_wp*dt*u(i,j)*&
               0.25_wp*(w(i+1,j)+w(i,j)+w(i+1,j-1)+w(i,j-1))*&
               ((2.0_wp*phi_n(i+1,j+1)+2.0_wp*phi_n(i,j+1)-phi_n(i,j)-&
               phi_n(i+1,j)-phi_n(i+1,j-1)-phi_n(i,j-1))/&
               ((2.0_wp*phi_n(i+1,j+1)+2.0_wp*phi_n(i,j+1)+&
               phi_n(i+1,j-1)+phi_n(i,j-1)+eps)*h))
       enddo
    else
       j=nyu-1
       do i=2,nxu-1 
          u_t(i,j)=(abs(u(i,j))*h-dt*u(i,j)**2)*((phi_n(i+1,j)-phi_n(i,j))/&
               ((phi_n(i+1,j)+phi_n(i,j)+eps)*h))-0.5_wp*dt*u(i,j)*&
               0.25_wp*(w(i+1,j)+w(i,j)+w(i+1,j-1)+w(i,j-1))*&
               ((phi_n(i+1,j+1)+phi_n(i,j+1)-phi_n(i+1,j-1)-phi_n(i,j-1))/&
               ((phi_n(i+1,j+1)+phi_n(i,j+1)+phi_n(i+1,j-1)+phi_n(i,j-1)+&
               eps)*h))
       enddo
    endif
    if (bbdy.and.bbc=='con') then
       j=2
       do i=2,nxu-1 
          u_t(i,j)=(abs(u(i,j))*h-dt*u(i,j)**2)*((phi_n(i+1,j)-phi_n(i,j))/&
               ((phi_n(i+1,j)+phi_n(i,j)+eps)*h))-1.5_wp*dt*u(i,j)*&
               0.25_wp*(w(i+1,j)+w(i,j)+w(i+1,j-1)+w(i,j-1))*&
               ((phi_n(i+1,j+1)+phi_n(i,j+1)+phi_n(i,j)+phi_n(i+1,j)-&
               2.0_wp*phi_n(i+1,j-1)-2.0_wp*phi_n(i,j-1))/((phi_n(i+1,j+1)+&
               phi_n(i,j+1)+2.0_wp*phi_n(i+1,j-1)+2.0_wp*phi_n(i,j-1)+eps)*h))
       enddo
    else
       j=2
       do i=2,nxu-1  
          u_t(i,j)=(abs(u(i,j))*h-dt*u(i,j)**2)*((phi_n(i+1,j)-phi_n(i,j))/&
               ((phi_n(i+1,j)+phi_n(i,j)+eps)*h))-0.5_wp*dt*u(i,j)*&
               0.25_wp*(w(i+1,j)+w(i,j)+w(i+1,j-1)+w(i,j-1))*&
               ((phi_n(i+1,j+1)+phi_n(i,j+1)-phi_n(i+1,j-1)-phi_n(i,j-1))/&
               ((phi_n(i+1,j+1)+phi_n(i,j+1)+phi_n(i+1,j-1)+phi_n(i,j-1)+&
               eps)*h))
       enddo
    endif
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
  SUBROUTINE flux_lim(beta_up_u,beta_down_u,beta_up_w,beta_down_w,phi_n,phi,&
       u_t,w_t,dt,ng)
    USE btype
    USE bfunc2d
    IMPLICIT NONE
    REAL(WP), DIMENSION(:,:), INTENT(OUT) :: beta_up_u,beta_down_u,beta_up_w,&
         beta_down_w
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: phi_n,phi,u_t,w_t 
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
          beta_up_u(i,j)=(phi_n_max_u(i,j)-phi_n(i,j))/((dt/h)*&
               (max(u_t(i-1,j),0.0_wp)*phi_n(i-1,j)-min(u_t(i,j),0.0_wp)*&
               phi_n(i+1,j))+(dt/h)*(max(w_t(i,j-1),0.0_wp)*phi_n(i,j-1)-&
               min(w_t(i,j),0.0_wp)*phi_n(i,j+1))+eps)
          beta_down_u(i,j)=(phi_n(i,j)-phi_n_min_u(i,j))/((dt/h)*&
               (max(u_t(i,j),0.0_wp)*phi_n(i,j)-min(u_t(i-1,j),0.0_wp)*&
               phi_n(i,j))+(dt/h)*(max(w_t(i,j),0.0_wp)*phi_n(i,j)-&
               min(w_t(i,j-1),0.0_wp)*phi_n(i,j))+eps)
          beta_up_w(i,j)=(phi_n_max_w(i,j)-phi_n(i,j))/((dt/h)*&
               (max(u_t(i-1,j),0.0_wp)*phi_n(i-1,j)-min(u_t(i,j),0.0_wp)*&
               phi_n(i+1,j))+(dt/h)*(max(w_t(i,j-1),0.0_wp)*phi_n(i,j-1)-&
               min(w_t(i,j),0.0_wp)*phi_n(i,j+1))+eps)
          beta_down_w(i,j)=(phi_n(i,j)-phi_n_min_w(i,j))/((dt/h)*&
               (max(u_t(i,j),0.0_wp)*phi_n(i,j)-min(u_t(i-1,j),0.0_wp)*&
               phi_n(i,j))+(dt/h)*(max(w_t(i,j),0.0_wp)*phi_n(i,j)-&
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
END SUBROUTINE mpdata_flux
