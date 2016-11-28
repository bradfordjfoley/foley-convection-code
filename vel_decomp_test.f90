! Brad Foley
! 8/11/11 Implementing simpler algorithm from Patankar to solve x 
! and z momentum equations
!
! Finite volume descritization so grid points are at cell centers for 
! pressure, temperature; u velocity is at horizontal cell edges, and w 
! velocity is at vertical cell edges.  Grid for x-mom equation is nx+1 by nz, 
! and grid for z-mom equation is nx by nz+1  
!
! Decomposing domains to include all interior points including first 
! and last rows/columns, with a row/column on either side for ghost 
! points/boundary values 
!
PROGRAM mgpoisson2dp 
  USE mpi
  USE btype
  USE bfunc2d
  IMPLICIT NONE
  INTEGER, PARAMETER :: NX=32,NY=16
  INTEGER, PARAMETER :: nxc=4,nyc=2
  REAL(WP), PARAMETER :: asp=2.0
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: u,rhs,err,utrue,w
  REAL(WP), ALLOCATABLE, DIMENSION(:) :: ae,aw,an,as,f
  INTEGER, ALLOCATABLE, DIMENSION(:) :: mynx,myny,mynxu,mynyw
  REAL(WP), DIMENSION(:,:), POINTER :: uj_1,uj,wj
  REAL(WP) :: t1,t2,erra,h,h2,lbvt,rbvt,tbvt,bbvt,dtime,dt,time,bn,tol,&
       mstep,mtime,tolu,erru,rbvu,lbvu,tbvu,bbvu,rbvw,lbvw,tbvw,bbvw
  TYPE(ptr2d), ALLOCATABLE :: rho(:)
  TYPE(ptr1d), ALLOCATABLE :: dx(:),dy(:),x(:),y(:),dxu(:),dyw(:),xu(:),yw(:)
  INTEGER :: n,i,j,k,ierr,myid,numprocs,npx,npy,comm2d,ng2,n2,&
       nbrleft,nbrright,nbrtop,nbrbottom,nlocal,deficit,&
       dims(2),coords(2),ng,ngrid,jcycle,ncycles,ng1,sx,ex,sy,ey,&
       nx1,ny1,sx1,sy1,nfinal,t,nglow,filenumber,sxu,exu,syw,eyw
  LOGICAL :: periods(2),lhbdy,rhbdy,topbdy,botbdy
  INTEGER :: status(MPI_STATUS_SIZE)
  CHARACTER(LEN=10) :: aname,uname,errname 
  CHARACTER(LEN=3) :: lhbct,rhbct,topbct,botbct,lhbcu,rhbcu,topbcu,botbcu,&
       lhbcw,rhbcw,topbcw,botbcw
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
  ! Determine number of grid levels
  ng=nint(log(real(nx/asp))/log(2.0_wp))
  ! Determine lowest grid level 
  nglow=nint(log(real(nxc/asp))/log(2.0_wp))

!!!!!!!! Set boundary conditions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Boundary conditions for temperature (t)
  lhbct='con'
  rhbct='con'
  topbct='con'  
  botbct='con'
  lbvt=0.0_wp
  rbvt=0.0_wp
  tbvt=0.0_wp
  bbvt=0.0_wp
  ! Boundary conditions for horizontal velocity (u)
  lhbcu='con'
  rhbcu='con'
  topbcu='con'  
  botbcu='con'
  lbvu=0.0_wp
  rbvu=0.0_wp
  tbvu=0.0_wp
  bbvu=0.0_wp  
  ! Boundary conditions for vertical velocity (w) 
  lhbcw='con'
  rhbcw='con'
  topbcw='con'  
  botbcw='con'
  lbvw=0.0_wp
  rbvw=0.0_wp
  tbvw=0.0_wp
  bbvw=0.0_wp  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Set number of processors in x direction (npx) & y direction (npy) 
  npx=1
  npy=2

  ! Stuff for decomposition of domain
  dims(1)=npx      ! If not 0, specifies n procs in x direction
  dims(2)=npy      ! If not 0, specifies n procs in y direction
  if (lhbct=='per') then
     periods(1)=.true.
  else 
     periods(1)=.false.
  endif
  if (topbct=='per') then
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
  allocate(mynx(ng),myny(ng),mynxu(ng),mynyw(ng))
  mynx=0
  myny=0
  mynxu=0
  mynyw=0

  ! Decompose domain in x and y directions, and determine if my 
  ! subdomain is a boundary
  call MPI_Cart_get(comm2d,2,dims,periods,coords,ierr)
  call decomp(sx,ex,dims(1),coords(1),nx,lhbdy,rhbdy)
  call decomp(sy,ey,dims(2),coords(2),ny,botbdy,topbdy)
  print*,myid,sx,ex,sy,ey 
  ! for periodic side bc's on u, u-grid goes from 2 to nx+1
  if (rhbcu=='per'.and.lhbcu=='per') then
     call decompu(sxu,exu,dims(1),coords(1),nx+1)
  endif
  ! for constant top and bottom bc's on w, w-grid goes from 2 to ny
  if (topbcw=='con'.and.botbcw=='con') then
     call decompu(syw,eyw,dims(2),coords(2),ny)
  endif
  ! for constant side bc's on u, u-grid goes from 2 to nx 
  if (rhbcu=='con'.and.lhbcu=='con') then
     call decompu(sxu,exu,dims(1),coords(1),nx)
  endif
  print*,myid,sxu,exu,syw,eyw

  ! Set up X-Y grid for cell centers and for u/w at cell boundaries  
  ! Allocate pointers for x,y,dx,dy,xu,yw,dxu,dwy
  allocate(x(ng),y(ng),dx(ng),dy(ng),xu(ng),yw(ng),dxu(ng),dyw(ng))
  ! Makes X-Y grid on each grid level for cell centers 
  call mkgrid(sx,ex,sy,ey,asp,ng,nglow,lhbdy,lhbct,rhbdy,rhbct,&
     topbdy,topbct,botbdy,botbct,mynx,myny,x,y,dx,dy)
  ! Make X grid for horizontal velocity 
  call mkgridu(sxu,exu,asp,ng,nglow,mynxu,xu,dxu)
  ! Make Y grid for vertical velocity
  call mkgridu(syw,eyw,1.0_wp,ng,nglow,mynyw,yw,dyw)

  ! Allocate 2d arrays
  allocate(rhs(mynx(ng),myny(ng)),u(mynx(ng),myny(ng)),w(mynx(ng),myny(ng)))
  u=0.0_wp
  rhs=0.0_wp
  w=0.0_wp

  ! Set boundary conditions on solution
  if (lhbdy.and.lhbct=='con') then
     u(1,1:myny(ng))=lbvt
  endif
  if (rhbdy.and.rhbct=='con') then
     u(mynx(ng),1:myny(ng))=rbvt
  endif
  if (botbdy.and.botbct=='con') then
     u(1:mynx(ng),1)=bbvt
  endif
  if (topbdy.and.topbct=='con') then
     u(1:mynx(ng),myny(ng))=tbvt
  endif

  forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
     rhs(i,j)=sin(pi*x(ng)%d(i))*sin(pi*y(ng)%d(j))
  end forall

  ncycles=2
  ! Call multi-grid to solve for u
  t1=MPI_WTIME()
  call mglin2d(u,rhs,0.0_wp,ng,nglow,0,ncycles,mynx,myny,x,y,dx,dy,&
       myid,nbrleft,nbrright,nbrtop,nbrbottom,.false.,.false.,comm2d)
  t2=MPI_WTIME()

  call nbrex2d(rhs,.true.,.true.,myid,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)

  if (myid==0) then 
     print*,rhs(1:mynx(ng),myny(ng))
     print*
     print*,rhs(1:mynx(ng),myny(ng)-1)
  endif

  ! Allocate space for analytical solution and error
  allocate(utrue(mynx(ng),myny(ng)),err(mynx(ng),myny(ng)))
  utrue=0.0_wp
  err=0.0_wp

  ! Analytical solution for sin wave forcing
  forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
     utrue(i,j)=(-1.0_wp/(2.0_wp*(pi**2)))*sin(pi*x(ng)%d(i))*&
          sin(pi*y(ng)%d(j))
  end forall

  ! Calculate error in solution versus analytical solution
  forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
     err(i,j)=100.0_wp*(abs(u(i,j)-utrue(i,j)))/abs(utrue(i,j))
  end forall

  ! sum error
  erra=0.0_wp
  do i=2,mynx(ng)-1
     do j=2,mynyw(ng)-1
        erra=erra+err(i,j)
     enddo
  enddo
  print*,"average percent error for",&
       myid,"is",erra/((mynx(ng)-2)*(myny(ng)-2))
  print*,myid,maxval(err)

  if (myid==0) then
     print*
     print*,"multi-grid time is",t2-t1
  endif

!!$  forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
!!$     u(i,j)=sin(pi*y(ng)%d(j))
!!$  end forall 
!!$  call nbrex2d(u,myid,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
!!$
!!$  ! Solve unsteady conduction problem
!!$  h=1.0_wp/(2**ng)
!!$  h2=h*h
!!$  ! Max timestep size
!!$  dt=h2/4.0_wp
!!$  ! Make timestep size smaller than max 
!!$  dtime=dt*1e-2
!!$  ! Time spacing and time counter for movie files
!!$  mstep=0.05_wp
!!$  mtime=mstep
!!$  filenumber=0
!!$
!!$  ! Determine finite volume coefficients to construct RHS
!!$  ! from solution from previous timestep
!!$  allocate(ae(mynx(ng)-2),aw(mynx(ng)-2),an(myny(ng)-2),as(myny(ng)-2))
!!$  aw(1:mynx(ng)-2)=h/dx(ng)%d(1:mynx(ng)-2)
!!$  ae(1:mynx(ng)-2)=h/dx(ng)%d(2:mynx(ng)-1)
!!$  as(1:myny(ng)-2)=h/dy(ng)%d(1:myny(ng)-2)
!!$  an(1:myny(ng)-2)=h/dy(ng)%d(2:myny(ng)-1)
!!$
!!$  ! Allocate space for analytical solution and error
!!$  allocate(utrue(mynx(ng),myny(ng)),err(mynx(ng),myny(ng)))
!!$  utrue=0.0_wp
!!$  err=0.0_wp
!!$
!!$  tol=1e-9
!!$  nfinal=50000
!!$  t=0
!!$  time=0.0_wp
!!$  ! number of multigrid cycles
!!$  ncycles=2

!!$  ! open file which contains timeseries info	
!!$  open(21,status='unknown',file='data.t')
!!$  do while (time.lt.1e-4)
!!$
!!$     ! Record whole solution, analytical solution, and percent error at a 
!!$     ! certain time spacing (mstep)
!!$     if (time.eq.0.0_wp) then
!!$        ! Analytical solution for cooling sine initial temperature
!!$        forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
!!$           utrue(i,j)=exp(-pi**2*time)*sin(pi*y(ng)%d(j))
!!$        end forall
!!$        call nbrex2d(utrue,myid,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
!!$        ! Calculate error in solution versus analytical solution
!!$        forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
!!$           err(i,j)=100.0_wp*(abs(u(i,j)-utrue(i,j)))/abs(utrue(i,j))
!!$        end forall
!!$        call nbrex2d(err,myid,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
!!$        write(aname,"('f000.',i1)") myid 
!!$        open(10,file=aname)
!!$        do i=1,mynx(ng)
!!$           do j=1,myny(ng)
!!$              write(10,'(5(1pe11.4,1x))') x(ng)%d(i),y(ng)%d(j),u(i,j),&
!!$                   utrue(i,j),err(i,j)
!!$           enddo
!!$        enddo
!!$        close(10)
!!$     endif      
!!$     if (time >= mtime ) then
!!$        mtime = time + mstep
!!$        filenumber = filenumber + 1
!!$        if (filenumber.lt.10) then
!!$           write(aname,"('f','00',i1,'.',i1)") filenumber,myid
!!$        elseif ((filenumber.ge.10).and.(filenumber.le.99)) then
!!$           write(aname,"('f','0',i2,'.',i1)") filenumber,myid
!!$        else
!!$           write(aname,"('f',i3,'.',i1)") filenumber,myid
!!$        endif
!!$        forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
!!$           utrue(i,j)=exp(-pi**2*time)*sin(pi*y(ng)%d(j))
!!$        end forall
!!$        call nbrex2d(utrue,myid,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
!!$        ! Calculate error in solution versus analytical solution
!!$        forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
!!$           err(i,j)=100.0_wp*(abs(u(i,j)-utrue(i,j)))/abs(utrue(i,j))
!!$        end forall
!!$        call nbrex2d(err,myid,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
!!$        open(10,file=aname)
!!$        do i=1,mynx(ng)
!!$           do j=1,myny(ng)
!!$              write(10,'(5(1pe11.4,1x))') x(ng)%d(i),y(ng)%d(j),u(i,j),&
!!$                   utrue(i,j),err(i,j)
!!$           enddo
!!$        enddo
!!$        close(10)
!!$     endif
!!$     
!!$     ! Every 20 timesteps, print out error information
!!$     if (mod(t,20)==0) then
!!$        ! Analytical solution for cooling sine initial temperature
!!$        forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
!!$           utrue(i,j)=exp(-pi**2*time)*sin(pi*y(ng)%d(j))
!!$        end forall
!!$        ! Calculate error in solution versus analytical solution
!!$        forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
!!$           err(i,j)=100.0_wp*(abs(u(i,j)-utrue(i,j)))/abs(utrue(i,j))
!!$        end forall
!!$
!!$        ! sum error
!!$        erra=0.0_wp
!!$        do i=2,mynx(ng)-1
!!$           do j=2,myny(ng)-1
!!$              erra=erra+err(i,j)
!!$           enddo
!!$        enddo
!!$        write(21,'(I1.2,4(1pe11.4,1x))') myid,&
!!$             time,erra/((mynx(ng)-2)*(myny(ng)-2)),t2-t1,tolu
!!$     endif
!!$
!!$     ! Construct RHS from solution at current timestep
!!$     do j=2,myny(ng)-1
!!$        do i=2,mynx(ng)-1
!!$           rhs(i,j)=0.5_wp*(ae(i-1)*u(i+1,j) + aw(i-1)*u(i-1,j) + &
!!$                an(j-1)*u(i,j+1) + as(j-1)*u(i,j-1)) + (h2/dtime - &
!!$                0.5_wp*(ae(i-1)+aw(i-1)+an(j-1)+as(j-1)))*u(i,j)
!!$        enddo
!!$     enddo
!!$     n=0
!!$     t1=MPI_WTIME()
!!$     do while (n.lt.nfinal) 
!!$        n=n+1
!!$        call jacobi2dfvp(u,rhs,ng,w,dx(ng)%d,dy(ng)%d,dtime)
!!$        call nbrex2d(w,myid,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
!!$        call jacobi2dfvp(w,rhs,ng,u,dx(ng)%d,dy(ng)%d,dtime)
!!$        call nbrex2d(u,myid,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
!!$        tolu=maxval(u)*tol
!!$        erru=maxval(abs(u-w))
!!$        print*,t,n,tolu,erru
!!$        if (erru.lt.tolu) exit
!!$     enddo
!!$     t2=MPI_WTIME()
!!$
!!$     ! Call multigrid to solve at next timestep
!!$     t1=MPI_WTIME()
!!$     call mglin2d(u,rhs,dtime,ng,nglow,ncycles,mynx,myny,x,y,dx,dy,&
!!$          myid,nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
!!$     t2=MPI_WTIME()
!!$
!!$     ! Update time and timestep
!!$     time = time + dtime
!!$     t=t+1
!!$  enddo
!!$  close(21)
  call MPI_FINALIZE(ierr)
END PROGRAM mgpoisson2dp
