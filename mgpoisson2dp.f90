! Brad Foley
! 8/1/11 Time-dependent conduction added and works for mg with multiple 
! processors; since time steps are small and change from previous solution to 
! current solution is small, mg is "overkill" and just doing a few jacobi 
! iterations converges on the solution quicker than mg
!
! Changing relaxation routine to work on interior points, then use 
! sepearte routine to calculate boundary points so that boundary 
! conditions are flexible
! Finite volume descritization so grid points are at cell centers 
! (i.e. for temperature equation) 
! Decomposing domains to include all interior points including first 
! and last rows/columns, with a row/column on either side for ghost 
! points/boundary values 
!
PROGRAM mgpoisson2dp 
  USE mpi
  USE btype
  USE bfunc2d
  IMPLICIT NONE
  INTEGER, PARAMETER :: NX=128,NY=128
  INTEGER, PARAMETER :: nxc=2,nyc=2
  REAL(WP), PARAMETER :: asp=1.0
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: u,rhs,err,utrue,w
  REAL(WP), ALLOCATABLE, DIMENSION(:) :: ae,aw,an,as,f
  INTEGER, ALLOCATABLE, DIMENSION(:) :: mynx,myny
  REAL(WP), DIMENSION(:,:), POINTER :: uj_1,uj,wj
  REAL(WP) :: t1,t2,erra,h,h2,lbv,rbv,tbv,bbv,dtime,dt,time,bn,tol
  TYPE(ptr2d), ALLOCATABLE :: rho(:)
  TYPE(ptr1d), ALLOCATABLE :: dx(:),dy(:),x(:),y(:)
  INTEGER :: n,i,j,k,ierr,myid,numprocs,npx,npy,comm2d,ng2,n2,&
       nbrleft,nbrright,nbrtop,nbrbottom,nlocal,deficit,&
       dims(2),coords(2),ng,ngrid,jcycle,ncycles,ng1,sx,ex,sy,ey,&
       nx1,ny1,sx1,sy1,nfinal,t
  INTEGER :: nglow
  LOGICAL :: periods(2),lhbdy,rhbdy,topbdy,botbdy
  INTEGER :: status(MPI_STATUS_SIZE)
  CHARACTER(LEN=10) :: aname,uname,errname 
  CHARACTER(LEN=3) :: lhbc,rhbc,topbc,botbc
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
  ! Determine number of grid levels
  ng=nint(log(real(nx/asp))/log(2.0_wp))
  ! Determine lowest grid level 
  nglow=nint(log(real(nxc/asp))/log(2.0_wp))

  ! Set boundary conditions
  lhbc='per'
  rhbc='per'
  topbc='con'  
  botbc='con'
  lbv=0.0_wp
  rbv=0.0_wp
  tbv=0.0_wp
  bbv=0.0_wp

  ! Set number of processors in x direction (npx) & y direction (npy) 
  npx=2
  npy=2

  ! Stuff for decomposition of domain
  dims(1)=npx      ! If not 0, specifies n procs in x direction
  dims(2)=npy      ! If not 0, specifies n procs in y direction
  if (lhbc=='per') then
     periods(1)=.true.
  else 
     periods(1)=.false.
  endif
  if (topbc=='per') then
     periods(2)=.true.
  else
     periods(2)=.false.
  endif
  if (myid == 0) then
     print*,'number of grids=',ng
     print*,'bottom grid level=',nglow
     if (npx*npy /= numprocs) then
        print*,'Number of processors must equal npx*npy*npz!'
        print*,'Program terminated!'
        stop
     endif
  endif

  ! Set up domain in mpi (for communication)
  call MPI_DIMS_CREATE(numprocs,2,dims,ierr)
  call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periods,.true.,comm2d,ierr)
  ! get position in the domain and neighbors
  call MPI_COMM_RANK(comm2d,myid,ierr)
  call MPI_Cart_shift(comm2d,0,1,nbrleft,nbrright,ierr)
  call MPI_Cart_shift(comm2d,1,1,nbrbottom,nbrtop,ierr)

  ! Allocate starting, ending arrays
  allocate(mynx(ng),myny(ng))
  mynx=0
  myny=0

  ! Decompose domain in x and y directions, and determine if my 
  ! subdomain is a boundary
  call MPI_Cart_get(comm2d,2,dims,periods,coords,ierr)
  call decomp(sx,ex,dims(1),coords(1),nx,lhbc,rhbc,lhbdy,rhbdy)
  call decomp(sy,ey,dims(2),coords(2),ny,botbc,topbc,botbdy,topbdy)
  print*,myid,sx,ex,sy,ey 

  ! Set up X-Y grid; allocate and determine x,y on each grid level
  allocate(x(ng),y(ng),dx(ng),dy(ng))
  call mkgrid(sx,ex,sy,ey,asp,ng,nglow,lhbdy,lhbc,rhbdy,rhbc,&
     topbdy,topbc,botbdy,botbc,mynx,myny,x,y,dx,dy)

  ! Allocate 2d arrays
  allocate(rhs(mynx(ng),myny(ng)),u(mynx(ng),myny(ng)),w(mynx(ng),myny(ng)))
  u=0.0_wp
  rhs=0.0_wp
  w=0.0_wp

  ! Set boundary conditions on solution
  if (lhbdy.and.lhbc=='con') then
     u(1,1:myny(ng))=lbv
  endif
  if (rhbdy.and.rhbc=='con') then
     u(mynx(ng),1:myny(ng))=rbv
  endif
  if (botbdy.and.botbc=='con') then
     u(1:mynx(ng),1)=bbv
  endif
  if (topbdy.and.topbc=='con') then
     u(1:mynx(ng),myny(ng))=tbv
  endif

  forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
     u(i,j)=sin(pi*y(ng)%d(j))
  end forall 
  call nbrex2d(u,myid,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)

!!$  forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
!!$     rhs(i,j)=sin(pi*x(ng)%d(i))*sin(pi*y(ng)%d(j))
!!$  end forall 

  ! Solve unsteady conduction problem
  h=1.0_wp/(2**ng)
  h2=h*h
  ! Max timestep size
  dt=h2/4.0_wp
  ! Make timestep size smaller than max 
  dtime=dt*1e-2

  ! Determine finite volume coefficients to construct RHS
  ! from solution from previous timestep
  allocate(ae(mynx(ng)-2),aw(mynx(ng)-2),an(myny(ng)-2),as(myny(ng)-2))
  aw(1:mynx(ng)-2)=h/dx(ng)%d(1:mynx(ng)-2)
  ae(1:mynx(ng)-2)=h/dx(ng)%d(2:mynx(ng)-1)
  as(1:myny(ng)-2)=h/dy(ng)%d(1:myny(ng)-2)
  an(1:myny(ng)-2)=h/dy(ng)%d(2:myny(ng)-1)

  tol=1e-9
  nfinal=50000
  t=0
  time=0.0_wp
  ! number of multigrid cycles
  ncycles=2
  do while (time.lt.1e-3)
     ! Now construct RHS
     do j=2,myny(ng)-1
        do i=2,mynx(ng)-1
           rhs(i,j)=0.5_wp*(ae(i-1)*u(i+1,j) + aw(i-1)*u(i-1,j) + &
                an(j-1)*u(i,j+1) + as(j-1)*u(i,j-1)) + (h2/dtime - &
                0.5_wp*(ae(i-1)+aw(i-1)+an(j-1)+as(j-1)))*u(i,j)
        enddo
     enddo
     i=0
!!$     t1=MPI_WTIME()
!!$     do while (i.lt.nfinal) 
!!$        i=i+1
!!$        call jacobi2dfvp(u,rhs,ng,w,dx(ng)%d,dy(ng)%d,dtime)
!!$        call nbrex2d(w,myid,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
!!$        call jacobi2dfvp(w,rhs,ng,u,dx(ng)%d,dy(ng)%d,dtime)
!!$        call nbrex2d(u,myid,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
!!$        if (maxval(abs(u-w)).lt.tol) exit
!!$     enddo
!!$     t2=MPI_WTIME()
     ! Call multigrid solver
     t1=MPI_WTIME()
     call mglin2d(u,rhs,dtime,ng,nglow,ncycles,mynx,myny,x,y,dx,dy,&
          myid,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
     t2=MPI_WTIME()
     time = time + dtime
     t=t+1
  enddo

  print*,time,t

  ! Allocate space for analytical solution and error
  allocate(utrue(mynx(ng),myny(ng)),err(mynx(ng),myny(ng)))
  utrue=0.0_wp
  err=0.0_wp

!!$  ! Analytical solution for sin wave forcing
!!$  forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
!!$     utrue(i,j)=(-1.0_wp/(2.0_wp*(pi**2)))*sin(pi*x(ng)%d(i))*&
!!$          sin(pi*y(ng)%d(j))
!!$  end forall

!!$  ! Analytical solution for conduction across a layer
!!$  forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
!!$     utrue(i,j)=1.0_wp-y(ng)%d(j)
!!$  end forall 

  ! Analytical solution for cooling sine initial temperature
  forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
     utrue(i,j)=exp(-pi**2*time)*sin(pi*y(ng)%d(j))
  end forall 

!!$  print*
!!$  print*,u(1,1:myny(ng))
!!$  print*
!!$  print*,utrue(2,1:myny(ng))

  ! Calculate error in solution versus analytical solution
  forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
     err(i,j)=100.0_wp*(abs(u(i,j)-utrue(i,j)))/abs(utrue(i,j))
  end forall

  ! sum error
  erra=0.0_wp
  do i=2,mynx(ng)-1
     do j=2,myny(ng)-1
        erra=erra+err(i,j)
     enddo
  enddo
  print*,"average percent error for",myid,"is",erra/((mynx(ng)-2)*(myny(ng)-2))
  print*,myid,maxval(err)

  if (myid==0) then
     print*
     print*,"multi-grid time is",t2-t1
  endif

  call MPI_FINALIZE(ierr)
END PROGRAM mgpoisson2dp
