! Subroutine to calculate finite volume discretization coefficients ae, aw, an,
! and as.  Mu is the viscosity (or "conductivity"), deltax and deltay are the 
! control volume sizes, dy is the vertical grid spacing for the variable 
! that we will eventually solve for (horizontal grid spacing, dx, is same as
! deltax for calculating horizontal velocity, u, so dx is not needed)
SUBROUTINE calc_au(ng,nglow,ae,aw,an,as,mu,mynx,myny,mynxu,mynyw,&
     deltax,deltay,dy,x,y,xu,yw,myid,nbrleft,nbrright,nbrtop,&
     nbrbottom,comm2d,lhbdy,lhbc,rhbdy,rhbc,topbdy,topbc,botbdy,botbc)
  USE btype 
  USE bfunc2d
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: ng,nglow
  TYPE(ptr2d), INTENT(INOUT) :: ae(:),aw(:),an(:),as(:)
  REAL(WP), DIMENSION(:,:), INTENT(IN) :: mu
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny,mynxu,mynyw
  TYPE(ptr1d), INTENT(IN) :: deltax(:),deltay(:),dy(:),x(:),y(:),&
       xu(:),yw(:)
  INTEGER(I4B), INTENT(IN) :: myid,nbrleft,nbrright,nbrtop,nbrbottom,comm2d
  LOGICAL, INTENT(IN) :: lhbdy,rhbdy,topbdy,botbdy
  CHARACTER(LEN=3), INTENT(IN) :: lhbc,rhbc,topbc,botbc
  REAL(WP), ALLOCATABLE , DIMENSION(:,:) :: mu_c
  REAL(WP), DIMENSION(:,:), POINTER :: mu1,mu2
  INTEGER(I4B) :: i,j,ngrid

  ! Determine coefficients on finest grid level
  allocate(mu_c(mynxu(ng),mynyw(ng)))
  mu_c=0.0_wp
  ! Use interpolation function to get mu at cell corners
  mu_c=interpgen2d(mu,x(ng)%d,y(ng)%d,xu(ng)%d,yw(ng)%d)
  ! NOW CALCULATE COEFFICIENTS 
  allocate(ae(ng)%a(mynxu(ng)-2,myny(ng)-2),aw(ng)%a(mynxu(ng)-2,myny(ng)-2),&
       an(ng)%a(mynxu(ng)-2,myny(ng)-2),as(ng)%a(mynxu(ng)-2,myny(ng)-2))
  forall (i=1:mynxu(ng)-2,j=1:myny(ng)-2)
     ae(ng)%a(i,j)=(2.0_wp*mu(i+2,j+1)*deltay(ng)%d(j))/deltax(ng)%d(i+1)
     aw(ng)%a(i,j)=(2.0_wp*mu(i+1,j+1)*deltay(ng)%d(j))/deltax(ng)%d(i)
  end forall
  forall (i=1:mynxu(ng)-2,j=1:myny(ng)-2) 
     an(ng)%a(i,j)=(mu_c(i+1,j+1)*deltax(ng)%d(i))/dy(ng)%d(j+1)
     as(ng)%a(i,j)=(mu_c(i+1,j)*deltax(ng)%d(i))/dy(ng)%d(j)
  end forall
  deallocate(mu_c)

  ! Set boundary coefficients to zero if bc is free slip 
  if (lhbdy.and.lhbc=='fsl') then
     aw(ng)%a(1,:)=0.0_wp
  endif
  if (rhbdy.and.rhbc=='fsl') then
     ae(ng)%a(mynx(ng)-2,:)=0.0_wp
  endif
  if (topbdy.and.topbc=='fsl') then
     an(ng)%a(:,myny(ng)-2)=0.0_wp
  endif
  if (botbdy.and.botbc=='fsl') then
     as(ng)%a(:,1)=0.0_wp
  endif

  ! Now calculate coefficients on coarser grids
  allocate(mu1(mynx(ng),myny(ng)))
  mu1=mu
  ngrid=ng
  do 
     if (ngrid<=nglow) exit
     ngrid=ngrid-1
     mu2=>mu1
     allocate(mu_c(mynxu(ngrid),mynyw(ngrid)),mu1(mynx(ngrid),myny(ngrid)))
     mu_c=0.0_wp
     mu1=0.0_wp
     ! Restrict viscosity to coraser grid
     mu1=rstrct2dp(mu2,ngrid+1,mynx,myny,0)
     deallocate(mu2)
     call nbrex2d(mu1,.false.,.false.,.false.,.false.,myid,nbrleft,&
          nbrright,nbrtop,nbrbottom,comm2d)
     ! Use interpolation function to get mu at cell corners
     mu_c=interpgen2d(mu1,x(ngrid)%d,y(ngrid)%d,xu(ngrid)%d,yw(ngrid)%d)
     ! NOW CALCULATE COEFFICIENTS 
     allocate(ae(ngrid)%a(mynxu(ngrid)-2,myny(ngrid)-2),&
          aw(ngrid)%a(mynxu(ngrid)-2,myny(ngrid)-2),&
          an(ngrid)%a(mynxu(ngrid)-2,myny(ngrid)-2),&
          as(ngrid)%a(mynxu(ngrid)-2,myny(ngrid)-2))
     forall (i=1:mynxu(ngrid)-2,j=1:myny(ngrid)-2)
        ae(ngrid)%a(i,j)=(2.0_wp*mu1(i+2,j+1)*deltay(ngrid)%d(j))&
             /deltax(ngrid)%d(i+1)
        aw(ngrid)%a(i,j)=(2.0_wp*mu1(i+1,j+1)*deltay(ngrid)%d(j))&
             /deltax(ngrid)%d(i)
     end forall
     forall (i=1:mynxu(ngrid)-2,j=1:myny(ngrid)-2) 
        an(ngrid)%a(i,j)=(mu_c(i+1,j+1)*deltax(ngrid)%d(i))/dy(ngrid)%d(j+1)
        as(ngrid)%a(i,j)=(mu_c(i+1,j)*deltax(ngrid)%d(i))/dy(ngrid)%d(j)
     end forall

     ! Set boundary coefficients to zero if bc is free slip 
     if (lhbdy.and.lhbc=='fsl') then
        aw(ngrid)%a(1,:)=0.0_wp
     endif
     if (rhbdy.and.rhbc=='fsl') then
        ae(ngrid)%a(mynx(ngrid)-2,:)=0.0_wp
     endif
     if (topbdy.and.topbc=='fsl') then
        an(ngrid)%a(:,myny(ngrid)-2)=0.0_wp
     endif
     if (botbdy.and.botbc=='fsl') then
        as(ngrid)%a(:,1)=0.0_wp
     endif

     deallocate(mu_c)
  enddo
  deallocate(mu1)
     
END SUBROUTINE calc_au
  
