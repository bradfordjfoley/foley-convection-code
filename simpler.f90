! 9/23/16: Code now accepts non-square control volumes where dx /= dy
! and irregular grids, as long as cells are always rectangular. Not 
! fully tested; will test first momentum solvers, then new mpdata 
! subroutine
!
! 8/25/14: Added tracer ratio method for advecting composition; only works 
! for having a step function like change in chemical buoyancy, e.g. "dense" 
! particles and "normal" particles.  
!
! 11/12/2013: Updated to accept reflecting side boundary conditions 
!
! 6/18/2012: I now set the coefficient to zero at the boundary where the bc is
! no flux (for u and w).  This fixes an error I had in the pressure equation 
! where the coefficents (ae and aw) at the top and bottom boundary were not 
! consistent with setting an/as to zero
!
! SIMPLER algorithm from Patankar for solving momentum equations on a 
! staggered grid.  With current finite volume discretization and 
! parallelization approach boundary conditions for u and w are as follows: 
!
! U top and bottom boundaries- no flux (fsl), constant value (con), or 
! periodic (per), though using periodic top/bot boundaries would be rare
!
! U left and right - con and per, though per will be by far most commonly used
!
! W top and bottom - con and per, though constant will be almost always used 
! (due to non-penetrating, free slip boundaries) 
!
! W left and right - con,per,fsl; Periodic boundaries will be most common, 
! but others could be useful
!
! For all other variables to calculate (temperature, composition, etc.) 
! con,per, and fsl conditions can be used on any boundary
PROGRAM simpler 
  USE mpi
  USE btype
  USE bfunc2d
  USE bfunc2d_main
  IMPLICIT NONE
  INTEGER :: nx,ny
  INTEGER :: nxc,nyc
  REAL(WP) :: asp
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: u,w,u_rhs,w_rhs,u_hat,w_hat,err,&
       utrue,mu,mu_c,temp,press,p_rhs,press_last,u_last,w_last,p_cor,&
       u_last_p,w_last_p,chem,alpha,u2,w2,press2,chem2,temp2,&
       du_dx,du_dy,dw_dx,dw_dy,deform,mu_t,source,sink,temp_tot,chem_tot,&
       alpha_tot,mu_tot,u_tot,w_tot,press_tot,errvel,alpha2,u_big,w_big,&
       tau,tau_xx,tau_zz,tau_xz,tau_tot,tauxx_tot,tauzz_tot,tauxz_tot,u_tmp,&
       tracers_tot,r_bound,l_bound,s_bound,n_bound,n_tr_weight_dense,&
       n_tr_weight_reg,errua,errwa,errpa,erru_all,errw_all,errp_all,&
       errvel_all,theta,dissipation,dissipation_tot
  REAL(WP), ALLOCATABLE, DIMENSION(:) :: x_tot,y_tot,xu_tot,yw_tot,temp_ha,&
       chem_ha,alpha_ha,mu_ha,xu_tmp,d_chem_layer
  INTEGER, ALLOCATABLE, DIMENSION(:) :: mynx,myny,mynxu,mynyw
  REAL(WP), DIMENSION(:,:), POINTER :: uj_1,uj,wj
  REAL(WP) :: t1,t2,erra,h,h2,lbvt,rbvt,tbvt,bbvt,dtime,dt,dta,dtd,time,bn,&
       mstep,mtime,rbvu,lbvu,tbvu,bbvu,rbvw,lbvw,tbvw,bbvw,Ra,xx1,xx2,&
       xx3,xx4,xx5,relax_p,relax,tol,tolu,tolw,tolp,erru,errw,errp,etae,&
       etat,etan,mv,t0,lewis,time_max,nu,v_rms,u_surf,rn,rbvc,lbvc,&
       tbvc,bbvc,Bu,tolu_tot,tolw_tot,tolp_tot,erru_tot,errw_tot,errp_tot,&
       p_scale,p_scale0,umin,umax,umin_tot,umax_tot,t3,dta_tot,&
       xx6,xx7,dtc,deltax,deltay,D,H_dam,Eh,m_alpha,p_alpha,t_off,&
       errv_tot,errv,maxdu,maxdu_tot,maxdw,maxdw_tot,tolv_tot,tolv,u0,&
       alpha_max,alpha_maxtot,timestep,press_max,press_max_tot,press_min,&
       press_min_tot,weight,ts_init,temp_cutoff,tau_avg,tauxx_avg,tauzz_avg,&
       tauxz_avg,temp_avg,alpha_avg,tau_xx_max,tauxx_maxtot,u_max,w_max,&
       max_vel,max_vel_tot,alpha_init,Q,d_crust,div_max,div_max_tot,cfluxavg,&
       avg_crustdepth,max_crustdepth,xnew,ynew,myxmin,myxmax,myymin,myymax,&
       myxmin_r,myxmax_r,myymin_r,myymax_r,t4,x_min_box,x_max_box,y_min_box,&
       y_max_box,trstep,trtime,mu_jump,x_zero,y_zero,diss_max,diss_maxtot,&
       deltax_tmp,dtd_tot,dtc_tot
  TYPE(ptr2d), ALLOCATABLE :: rho(:),aeu(:),awu(:),anu(:),asu(:),&
       aew(:),aww(:),anw(:),asw(:),ae_p(:),aw_p(:),an_p(:),as_p(:)
  TYPE(ptr1d), ALLOCATABLE :: dx(:),dy(:),x(:),y(:),dxu(:),dyw(:),xu(:),yw(:),&
       xw(:),yu(:),dxw(:),dyu(:),xc(:),yc(:),dxc(:),dyc(:),tracers(:),&
       oldtracers(:),tracers2(:)
  INTEGER :: n,i,j,k,l,ierr,myid,numprocs,npx,npy,comm2d,ng2,n2,&
       nbrleft,nbrright,nbrtop,nbrbottom,nlocal,deficit,&
       dims(2),coords(2),ng,ngrid,jcycle,ncycles,ng1,sx,ex,sy,ey,&
       nx1,ny1,sx1,sy1,nfinal,t,nglow,filenumber,sxu,exu,syw,eyw,iter,loop,&
       iord,iter_total,ibc,file_max,itime,ntdump,jkf,ntotal,iter_uw_loop,&
       p_cyc,uw_cyc,maxrow,maxroww,p_it,uw_it,dam_iord,ntracers,myntracers,&
       my_tracers_out_tot,nbrtopright,nbrtopleft,nbrbotright,nbrbotleft,&
       my_tracers_in_tot,myntracers_new,count,count_tr,count_tl,count_br,&
       count_bl,count_t,count_b,count_r,count_l,tracers_per_cell,tracers_ic,&
       i_start,i_end,j_start,j_end,filenumber_tr,ntracers_init,ni,asp_nx
  INTEGER, ALLOCATABLE, DIMENSION(:) :: sxall,syall,mynxall,mynyall,sxuall,&
       sywall,mynxuall,mynywall,eyall,eywall
  LOGICAL :: periods(2),lhbdy,rhbdy,topbdy,botbdy,utfsl,ubfsl,ttfsl,tbfsl,&
       ctfsl,cbfsl,remove_net_u,input_alpha,tlfsl,trfsl,clfsl,crfsl,wlfsl,&
       wrfsl,chem_layer_init,crust_production
  INTEGER :: status(MPI_STATUS_SIZE),seed_tracers
  INTEGER, DIMENSION(8) :: seed,my_tracers_in,my_tracers_out 
  INTEGER, DIMENSION(1000) :: tr_out_right,tr_out_left,tr_out_top,tr_out_bot,&
       tr_out_topright,tr_out_topleft,tr_out_botright,tr_out_botleft
  CHARACTER(LEN=100) :: init_f,init_u,init_w,init_tr
  CHARACTER(LEN=7) :: dname,uname,wname,hname,errname,tname,trname
  CHARACTER(LEN=3) :: lhbct,rhbct,topbct,botbct,lhbcu,rhbcu,topbcu,botbcu,&
       lhbcw,rhbcw,topbcw,botbcw,lhbcc,rhbcc,topbcc,botbcc
  NAMELIST/param/Ra,lewis,asp,nx,ny,nxc,nyc,relax,relax_p,tol,etae,etat,mv,&
       t0,lhbct,rhbct,topbct,botbct,lhbcu,rhbcu,topbcu,botbcu,lhbcw,rhbcw,&
       topbcw,botbcw,lbvt,rbvt,tbvt,bbvt,rbvu,lbvu,tbvu,bbvu,rbvw,lbvw,tbvw,&
       bbvw,npx,npy,iter_total,iord,ibc,init_f,init_u,init_w,time_max,&
       file_max,mstep,ntdump,botbcc,topbcc,lhbcc,rhbcc,rbvc,lbvc,tbvc,bbvc,Bu,&
       iter_uw_loop,p_cyc,uw_cyc,D,H_dam,Eh,m_alpha,p_alpha,t_off,timestep,&
       p_it,uw_it,dam_iord,weight,ts_init,temp_cutoff,remove_net_u,&
       input_alpha,alpha_init,Q,d_crust,chem_layer_init,crust_production,&
       tracers_per_cell,tracers_ic,init_tr,trstep,mu_jump,ntracers_init
  OPEN(4,status='old',file='conv.par',action='READ')
  READ(4,nml=param)
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
  ! Determine number of grid levels
  ng=nint(log(real(ny))/log(2.0_wp))
  ! Determine lowest grid level 
  nglow=nint(log(real(nyc))/log(2.0_wp))

  ! Total number of grid points 
  ntotal=nx*ny

  asp_nx=nx/ny

  ! Total number of tracers 
  if (tracers_ic==2) then
     ntracers=ntracers_init
  else
     ntracers=tracers_per_cell*ntotal
  endif

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
        print*,'Number of processors must equal npx*npy!'
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

  nbrtopright=-1
  nbrtopleft=-1
  nbrbotleft=-1
  nbrbotright=-1

  ! Each processor needs to get it's neighbor at the corners (e.g. top right, 
  ! etc.) 
  ! My top-right neighbor is nbrtop of my nbrright; so I send my nbrtop to my 
  ! nbrleft
  call MPI_SENDRECV(nbrtop,1,MPI_INTEGER,nbrleft,0,nbrtopright,1,MPI_INTEGER,&
       nbrright,0,comm2d,status,ierr)
  ! My top-left neighbor is nbrtop of my nbrleft; so I send my nbrtop to my 
  ! nbrright
  call MPI_SENDRECV(nbrtop,1,MPI_INTEGER,nbrright,1,nbrtopleft,1,MPI_INTEGER,&
       nbrleft,1,comm2d,status,ierr)
  ! My bot-right neighbor is nbrbot of my nbrright; so I send my nbrbot to my 
  ! nbrleft
  call MPI_SENDRECV(nbrbottom,1,MPI_INTEGER,nbrleft,2,nbrbotright,1,&
       MPI_INTEGER,nbrright,2,comm2d,status,ierr)
  ! My bot-left neighbor is nbrbot of my nbrleft; so I send my nbrbot to my 
  ! nbrright
  call MPI_SENDRECV(nbrbottom,1,MPI_INTEGER,nbrright,3,nbrbotleft,1,&
       MPI_INTEGER,nbrleft,3,comm2d,status,ierr)

  print*,myid,nbrtop,nbrbottom,nbrleft,nbrright,nbrtopleft,nbrtopright,&
       nbrbotleft,nbrbotright

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

  ! Set up X-Y grid for cell centers and for u/w at cell boundaries  
  ! Allocate pointers for x,y,dx,dy,xu,yw,dxu,dwy
  allocate(x(ng),y(ng),dx(ng),dy(ng),xu(ng),yw(ng),dxu(ng),dyw(ng),&
       xw(ng),yu(ng),dxw(ng),dyu(ng),xc(ng),yc(ng),dxc(ng),dyc(ng))
  ! Makes X-Y grid on each grid level for cell centers 
  call mkgrid(sx,ex,sy,ey,asp,ng,nglow,lhbdy,lhbct,lhbcw,lhbcc,rhbdy,rhbct,&
       rhbcw,rhbcc,topbdy,topbct,topbcu,topbcc,botbdy,botbct,botbcu,botbcc,&
       mynx,myny,x,xw,xc,y,yu,yc,dx,dxw,dxc,dy,dyu,dyc,asp_nx)

  ! Make X grid for horizontal velocity 
  call mkgridu(sxu,exu,asp,ng,nglow,mynxu,xu,dxu,asp_nx)
  ! Make Y grid for vertical velocity
  call mkgridu(syw,eyw,1.0_wp,ng,nglow,mynyw,yw,dyw,1)

  ! Allocate 2d arrays
  allocate(u_rhs(mynxu(ng),myny(ng)),w_rhs(mynx(ng),mynyw(ng)),&
       u(mynxu(ng),myny(ng)),w(mynx(ng),mynyw(ng)),u_hat(mynxu(ng),myny(ng)),& 
       w_hat(mynx(ng),mynyw(ng)),temp(mynx(ng),myny(ng)),&
       press(mynx(ng),myny(ng)),p_rhs(mynx(ng),myny(ng)),&
       press_last(mynx(ng),myny(ng)),u_last(mynxu(ng),myny(ng)),&
       w_last(mynx(ng),mynyw(ng)),p_cor(mynx(ng),myny(ng)),&
       w_last_p(mynx(ng),mynyw(ng)),u_last_p(mynxu(ng),myny(ng)),&
       mu(mynx(ng),myny(ng)),mu_c(mynxu(ng),mynyw(ng)),&
       du_dx(mynx(ng)-2,myny(ng)-2),dw_dy(mynx(ng)-2,myny(ng)-2),&
       du_dy(mynx(ng)-2,myny(ng)-2),dw_dx(mynx(ng)-2,myny(ng)-2),&
       alpha(mynx(ng),myny(ng)),mu_t(mynx(ng),myny(ng)),&
       errvel(mynx(ng),myny(ng)),tau(mynx(ng),myny(ng)),&
       tau_xx(mynx(ng),myny(ng)),tau_zz(mynx(ng),myny(ng)),&
       tau_xz(mynx(ng),myny(ng)),errua(mynxu(ng),myny(ng)),&
       errwa(mynx(ng),mynyw(ng)),errpa(mynx(ng),myny(ng)),&
       theta(mynx(ng),myny(ng)),dissipation(mynx(ng),myny(ng)),&
       d_chem_layer(mynx(ng)))
  ! Defining composition at same grid points as for vertical velocity 
  ! (only for when I am using tracers for composition)
  allocate(chem(mynx(ng),mynyw(ng)))

  !
  u=0.0_wp; w=0.0_wp; u_rhs=0.0_wp; w_rhs=0.0_wp; u_hat=0.0_wp
  w_hat=0.0_wp; temp=0.0_wp; p_rhs=0.0_wp; press=0.0_wp; press_last=0.0_wp
  u_last=0.0_wp; w_last=0.0_wp; u_last_p=0.0_wp; w_last_p=0.0_wp
  p_cor=0.0_wp; mu=0.0_wp; mu_c=0.0_wp; du_dx=0.0_wp; dw_dy=0.0_wp; 
  du_dy=0.0_wp; dw_dx=0.0_wp; chem=0.0_wp; alpha=1.0_wp; mu_t=0.0_wp
  dissipation=0.0_wp

  ! Allocate grid space for deformational work, grain growth source if 
  ! we are solving the damage equation
  if (D/=0.0_wp) then
     allocate(deform(mynx(ng)-2,myny(ng)-2),source(mynx(ng)-2,myny(ng)-2),&
          sink(mynx(ng)-2,myny(ng)-2))
     deform=0.0_wp
     source=0.0_wp
     sink=0.0_wp
  endif

  ! Set top/bottom or right/left fsl boundary conditions
  ! temp equation
  if (topbct=='fsl') then
     ttfsl=.true.
  else
     ttfsl=.false.
  endif
  if (botbct=='fsl') then
     tbfsl=.true.
  else
     tbfsl=.false.
  endif
  if (lhbct=='fsl') then
     tlfsl=.true.
  else
     tlfsl=.false.
  endif
  if (rhbct=='fsl') then
     trfsl=.true.
  else
     trfsl=.false.
  endif
  ! chem equation
  if (topbcc=='fsl') then
     ctfsl=.true.
  else
     ctfsl=.false.
  endif
  if (botbcc=='fsl') then
     cbfsl=.true.
  else
     cbfsl=.false.
  endif
  if (lhbcc=='fsl') then
     clfsl=.true.
  else
     clfsl=.false.
  endif
  if (rhbcc=='fsl') then
     crfsl=.true.
  else
     crfsl=.false.
  endif
  ! U equation (won't have no flux side boundaries)
  if (topbcu=='fsl') then
     utfsl=.true.
  else
     utfsl=.false.
  endif
  if (botbcu=='fsl') then
     ubfsl=.true.
  else
     ubfsl=.false.
  endif
  ! W equation (won't have no flux top boundaries)
  if (lhbcw=='fsl') then
     wlfsl=.true.
  else
     wlfsl=.false.
  endif
  if (rhbcw=='fsl') then
     wrfsl=.true.
  else
     wrfsl=.false.
  endif

  ! Gather all processors starting x and y positions onto root process
  ! for recombining data at output
  if (myid==0) then
     allocate(sxall(numprocs),syall(numprocs),mynxall(numprocs),&
          mynyall(numprocs),sxuall(numprocs),sywall(numprocs),&
          mynxuall(numprocs),mynywall(numprocs),x_tot(nx),y_tot(ny),&
          eyall(numprocs),eywall(numprocs))
  endif
  ! Every process will have the full xu and yw
  allocate(xu_tot(nx+1),yw_tot(ny+1))
  call MPI_GATHER(sx,1,MPI_INTEGER,sxall,1,MPI_INTEGER,0,comm2d,ierr)
  call MPI_GATHER(sy,1,MPI_INTEGER,syall,1,MPI_INTEGER,0,comm2d,ierr)
  call MPI_GATHER(ey,1,MPI_INTEGER,eyall,1,MPI_INTEGER,0,comm2d,ierr)
  call MPI_GATHER(sxu,1,MPI_INTEGER,sxuall,1,MPI_INTEGER,0,comm2d,ierr)
  call MPI_GATHER(syw,1,MPI_INTEGER,sywall,1,MPI_INTEGER,0,comm2d,ierr)
  call MPI_GATHER(eyw,1,MPI_INTEGER,eywall,1,MPI_INTEGER,0,comm2d,ierr)
  call MPI_GATHER(mynx(ng),1,MPI_INTEGER,mynxall,1,MPI_INTEGER,0,comm2d,ierr)
  call MPI_GATHER(mynxu(ng),1,MPI_INTEGER,mynxuall,1,MPI_INTEGER,0,comm2d,ierr)
  call MPI_GATHER(myny(ng),1,MPI_INTEGER,mynyall,1,MPI_INTEGER,0,comm2d,ierr)
  call MPI_GATHER(mynyw(ng),1,MPI_INTEGER,mynywall,1,MPI_INTEGER,0,comm2d,ierr)

  ! Root process also needs total grid, x, y, xu, yw
  if (myid==0) then
     do i=1,nx
        x_tot(i)=(1.0_wp/(nx/asp))*i-0.5_wp/(nx/asp)
     enddo
     do i=1,ny
        y_tot(i)=(1.0_wp/ny)*i-0.5_wp/ny
     enddo
  endif

  do i=1,nx+1
     xu_tot(i)=(1.0_wp/(nx/asp))*(i-1)
  enddo
  do i=1,ny+1
     yw_tot(i)=(1.0_wp/ny)*(i-1)
  enddo
  
  maxrow=ceiling(ny/real(npy))
  maxroww=ceiling((ny+1)/real(npy))

  ! Figure out the boundaries of my chunk of the domain
  myxmin=xu(ng)%d(1)
  myxmax=xu(ng)%d(mynxu(ng)-1)
  myymin=yw(ng)%d(1)
  myymax=yw(ng)%d(mynyw(ng)-1)
  if (rhbdy.and.rhbcu=='con') then
     myxmax=xu(ng)%d(mynxu(ng))
  endif
  if (topbdy.and.topbcw=='con') then
     myymax=yw(ng)%d(mynyw(ng))
  endif
  ! Now figure out boundaries of my sub-domain for use when reading in 
  ! tracer position data from file; this ensures that tracers with an x 
  ! position equal to the max x of the domain, or a y position of 1, still 
  ! get properly distributed 
  if (rhbdy.and.topbdy) then
     myxmax_r=myxmax+epsilon(1.0_wp)
     myymax_r=myymax+epsilon(1.0_wp)
     myxmin_r=myxmin
     myymin_r=myymin
  elseif (rhbdy) then
     myxmin_r=myxmin
     myxmax_r=myxmax+epsilon(1.0_wp)
     myymin_r=myymin
     myymax_r=myymax
  elseif (topbdy) then
     myxmin_r=myxmin
     myxmax_r=myxmax
     myymin_r=myymin
     myymax_r=myymax+epsilon(1.0_wp)
  else
     myxmin_r=myxmin
     myxmax_r=myxmax
     myymin_r=myymin
     myymax_r=myymax
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set up array of tracer position and composition/state data 
  ! element 1: composition (1 for dense, 0 for light)
  ! element 2: x coordinate, element 3: z coordinate
  if (ntracers>0) then
     if (tracers_ic==1) then
        myntracers=(mynx(ng)-2)*(myny(ng)-2)*tracers_per_cell
        ! Set up array of tracer position and composition/state data 
        ! element 1: composition (1 for dense, 0 for light)
        ! element 2: x coordinate, element 3: z coordinate
        allocate(tracers(myntracers))
        l=1
        ! tracers_ic=1 initializes tracers as a given number per cell (set by 
        ! tracers_per_cell), and puts then within each cell at random x,y 
        ! positions
        call system_clock(seed_tracers)
        seed=myid+seed_tracers*seed_tracers
        call random_seed(PUT=seed)
        call random_number(rn)
        call random_number(rn)
        do i=2,mynx(ng)-1
           do j=2,myny(ng)-1
              do k=1,tracers_per_cell
                 allocate(tracers(l)%d(3))
                 tracers(l)%d(1)=0.0_wp
                 call random_number(rn)
                 ! Initial x,y position for this single tracer test case
                 call random_number(rn)
                 tracers(l)%d(2)=x(ng)%d(i)+dxu(ng)%d(i-1)*(rn-0.5_wp)
                 call random_number(rn)
                 tracers(l)%d(3)=y(ng)%d(j)+dyw(ng)%d(j-1)*(rn-0.5_wp)
                 l=l+1
              enddo
           enddo
        enddo
        if (l-1.ne.myntracers) then
           print*,'Error initializng tracers!'  
           print*,'Mismatch between initialized tracers and &
                number of tracers allocated!'
           stop
        endif
     elseif (tracers_ic==2) then
        ! tracers_ic=2 reads in the tracer positions from a file
        ! like for temp,velocity, each processor reads in all of the tracer 
        ! position data, and picks out just the tracers in it's subdomain
        allocate(tracers2(ntracers))
        ! Read in tracer data
        open(1,file=init_tr)
        k=0
        do i=1,ntracers
           read(1,*) xx1,xx2,xx3,xx4
           allocate(tracers2(i)%d(3))
           tracers2(i)%d(1)=xx2
           tracers2(i)%d(2)=xx3
           tracers2(i)%d(3)=xx4
           if (xx3<myxmax_r.and.xx3>=myxmin_r.and.&
                xx4<myymax_r.and.xx4>=myymin_r) then
              ! Count number of tracers in my subdomain
              k=k+1
           endif
        enddo
        close(1)
        myntracers=k
        allocate(tracers(myntracers))
        l=1
        do i=1,ntracers
           if (tracers2(i)%d(2)<myxmax_r.and.tracers2(i)%d(2)>=myxmin_r.and.&
                tracers2(i)%d(3)<myymax_r.and.tracers2(i)%d(3)>=myymin_r) then
              allocate(tracers(l)%d(3))
              tracers(l)%d=tracers2(i)%d
              l=l+1
           endif
        enddo
        do i=1,ntracers
           deallocate(tracers2(i)%d)
        enddo
        deallocate(tracers2)
        if (l-1.ne.myntracers) then
           print*,'Error initializng tracers!'  
           print*,'Mismatch between initialized tracers and &
                number of tracers allocated!'
           stop
        endif
     elseif (tracers_ic==3) then
        ! Make a square with a high density of tracers in it
        ! Initially the root process will have all the tracers
        if (myid==0) then
           x_min_box=0.2_wp
           i_start=nint((x_min_box+(1.0_wp/(2.0_wp*ny)))*ny+2.0_wp-sx)
           x_max_box=0.3_wp
           i_end=nint((x_max_box+(1.0_wp/(2.0_wp*ny)))*ny+2.0_wp-sx)
           y_min_box=0.1_wp
           j_start=nint((y_min_box+(1.0_wp/(2.0_wp*ny)))*ny+2.0_wp-sy)
           y_max_box=0.2_wp
           j_end=nint((y_max_box+(1.0_wp/(2.0_wp*ny)))*ny+2.0_wp-sy)
           ntracers=(i_end-i_start+1)*(j_end-j_start+1)*tracers_per_cell
           myntracers=ntracers
           allocate(tracers(myntracers))
        else
           myntracers=0
           allocate(tracers(myntracers))
        endif
        call MPI_BCAST(ntracers,1,MPI_INTEGER,0,comm2d,ierr)
        if (myid==0) then
           l=1
           ! tracers_ic=1 initializes tracers as a given number per cell 
           ! (set by tracers_per_cell), and puts then within each cell 
           ! at random x,y positions
           do i=i_start,i_end
              do j=j_start,j_end
                 do k=1,tracers_per_cell
                    allocate(tracers(l)%d(3))
                    tracers(l)%d(1)=0.0_wp
                    !call system_clock(seed_tracers)
                    seed=l+myid+1
                    call random_seed(PUT=seed)
                    ! Initial x,y position for this single tracer test case
                    call random_number(rn)
                    tracers(l)%d(2)=x(ng)%d(i)+dxu(ng)%d(i)*(rn-0.5_wp)
                    call random_number(rn)
                    tracers(l)%d(3)=y(ng)%d(j)+dyw(ng)%d(j)*(rn-0.5_wp)
                    l=l+1
                 enddo
              enddo
           enddo
           if (l-1.ne.myntracers) then
              print*,'Error initializng tracers!'  
              print*,'Mismatch between initialized tracers and &
                   number of tracers allocated!'
              stop
           endif
        endif
     endif
  endif

  ! Initial condition 1 is a conductive state
  if (ibc==1) then
     do j=1,myny(ng)
        do i=1,mynx(ng)
           temp(i,j)=1.0_wp-y(ng)%d(j)
        enddo
     enddo
     if (botbdy) then
        seed=myid
        call random_seed(PUT=seed)
        do i=2,mynx(ng)-1
           call random_number(rn)
           temp(i,myny(ng)/2)=temp(i,myny(ng)/2) + 5E-3*(0.5_wp-rn) 
        enddo
     endif
     if (ntracers==0) then
        chem=0.0_wp
     endif
     press=0.0_wp
     u=0.0_wp
     w=0.0_wp
  ! Initial condition 2 is to restart from old results
  ! Have to put all data into dummy arrays, which will then 
  ! be parted out to each processor   
  elseif (ibc==2) then
     allocate(temp2(nx,ny),press2(nx,ny),chem2(nx,ny),u2(nx+1,ny),w2(nx,ny+1),&
          alpha2(nx,ny))
     ! Read in temperature field here
     if (ntracers>0) then
        open(1,file=init_f)
        do i=1,nx
           do j=1,ny
              read(1,*) xx1,xx2,xx3,xx4,xx5,xx6
              temp2(i,j) = xx3 
              alpha2(i,j) = xx4
              press2(i,j)= xx5
           enddo
        enddo
        close(1)
     else 
        open(1,file=init_f)
        do i=1,nx
           do j=1,ny
              read(1,*) xx1,xx2,xx3,xx4,xx5,xx6 !,xx7
              temp2(i,j) = xx3 
              alpha2(i,j) = xx4
              !press2(i,j)= 0.0_wp
              press2(i,j)= xx5
              !chem2(i,j) = xx7
              chem2(i,j) = 0.0_wp
           enddo
        enddo
        close(1)
     endif
     ! Read in horizontal velocity field (u) here
     open(1,file=init_u)
     do i=1,nx+1
        do j=1,ny
           read(1,*) xx1,xx2,xx3
           u2(i,j) = xx3
        enddo
     enddo
     close(1)
     ! Read in vertical velocity field (w) here
     open(1,file=init_w)
     do i=1,nx
        do j=1,ny+1
           read(1,*) xx1,xx2,xx3
           w2(i,j) = xx3 
        enddo
     enddo
     close(1)

     ! Now each processor takes it's piece of the data to work on
     temp(2:mynx(ng)-1,2:myny(ng)-1)=temp2(sx:ex,sy:ey)
     if (ntracers==0) then
        chem(2:mynx(ng)-1,2:myny(ng)-1)=chem2(sx:ex,sy:ey)
     endif
     press(2:mynx(ng)-1,2:myny(ng)-1)=press2(sx:ex,sy:ey)
     alpha(2:mynx(ng)-1,2:myny(ng)-1)=alpha2(sx:ex,sy:ey)
     u(2:mynxu(ng)-1,2:myny(ng)-1)=u2(sxu:exu,sy:ey)
     w(2:mynx(ng)-1,2:mynyw(ng)-1)=w2(sx:ex,syw:eyw)
     deallocate(temp2,chem2,press2,u2,w2,alpha2)
  ! Initial condition 3: hot mantle, cooling from above; no chemistry
  elseif (ibc==3) then
     do j=1,myny(ng)
        do i=1,mynx(ng)
           temp(i,j)=1.0_wp - (1.0_wp - ts_init)*y(ng)%d(j)
        enddo
     enddo
     seed=myid
     if (botbdy) then
        seed=myid
        call random_seed(PUT=seed)
        do i=2,mynx(ng)-1
           call random_number(rn)
           temp(i,myny(ng)/2)=temp(i,myny(ng)/2) + 5E-3*(0.5_wp-rn)
        enddo
     endif
     if (ntracers==0) then
        chem=0.0_wp
     endif
     press=0.0_wp
     u=0.0_wp
     w=0.0_wp 
  ! Initial condition 4 is a hot start, cooling from above
  elseif (ibc==4) then
     temp=1.0_wp
     do j=1,myny(ng)
        do i=1,mynx(ng)
           chem(i,j)=y(ng)%d(j)
        enddo
     enddo
     seed=myid
     if (botbdy) then
        seed=myid
        call random_seed(PUT=seed)
        do i=2,mynx(ng)-1
           call random_number(rn)
           chem(i,myny(ng)/2)=chem(i,myny(ng)/2) + 1E-3*(0.5_wp-rn)
        enddo
     endif
     press=0.0_wp
     u=0.0_wp
     w=0.0_wp 
  ! Initial condition 5 is an unstable linear chemical profile with small 
  ! unstable temp profile
  elseif (ibc==5) then
     do j=1,myny(ng)
        do i=1,mynx(ng)
           temp(i,j)=1.0_wp - (1.0_wp - ts_init)*y(ng)%d(j)
        enddo
     enddo
     seed=myid
     do j=1,myny(ng)
        do i=1,mynx(ng)
           chem(i,j)=y(ng)%d(j)
        enddo
     enddo
     if (botbdy) then
        call random_seed(PUT=seed)
        do i=2,mynx(ng)-1
           call random_number(rn)
           chem(i,myny(ng)/2)=chem(i,myny(ng)/2) + 5E-3*(0.5_wp-rn)
        enddo
     endif
     press=0.0_wp
     u=0.0_wp
     w=0.0_wp 
  ! Initial condition 6 is for shear calculation; linear profile with depth 
  ! in u. C,T = 0  
  elseif (ibc==6) then
     temp=0.0_wp
     chem=0.0_wp
     press=0.0_wp
     do j=1,myny(ng)
        do i=1,mynx(ng)
           u(i,j)=(tbvu-bbvu)*y(ng)%d(j) + bbvu 
        enddo
     enddo
     w=0.0_wp
  ! Initial condition for thermochem benchmark; t=0, u & w = 0  
  elseif (ibc==7) then
     temp=0.0_wp
     chem=0.0_wp
     press=0.0_wp
     u=0.0_wp
     w=0.0_wp
  endif

  ! If I want to specify the initial value of alpha everywhere
  if (input_alpha) then
     alpha=alpha_init
  endif

  !8/27/14: Making an initial chemical layer is now done with tracers when 
  !ntracers is greater than 0; otherwise it is done with the grid based scheme
  ! Set up a chemical layer at the top of the mantle
  if (chem_layer_init.and.ntracers==0) then
     do j=1,myny(ng)
        do i=1,mynx(ng)
           if (y(ng)%d(j).gt.d_crust) then
              chem(i,j)=1.0_wp
           else
              chem(i,j)=0.0_wp
           endif
        enddo
     enddo
  elseif (chem_layer_init.and.ntracers>0) then
     ! For RT problem in PvK thermo-chem benchmark paper
     d_chem_layer(1:mynx(ng))=d_crust+0.02_wp*cos(pi*x(ng)%d(1:mynx(ng))/asp)
     do i=1,myntracers
        ! Invert x-coordinate of tracer for nearest x-coordinate of grid
        ni=nint(tracers(i)%d(2)*(nx/asp)+2.5_wp-sx)
        if (tracers(i)%d(3)<d_chem_layer(ni)) then
           tracers(i)%d(1)=1.0_wp
        else
           tracers(i)%d(1)=0.0_wp
        endif
     enddo
     ! Now convert tracer positions/property (e.g. dense versus regular) 
     ! into compositional field
     ! call calc_chem_tr_ratio(tracers,x,y,chem,myntracers,nx,ny,sx,sy)
     call calc_chem_tr_ratio(tracers,chem,xw,y,myntracers,asp,nx,ny,sx,syw,&
          mynx(ng),mynyw(ng),nbrright,nbrleft,nbrtop,nbrbottom,comm2d)
  endif

  ! Set boundary conditions 
  call boundary(temp,lhbdy,rhbdy,topbdy,botbdy,lhbct,rhbct,&
         topbct,botbct,lbvt,rbvt,tbvt,bbvt)
  call boundary(chem,lhbdy,rhbdy,topbdy,botbdy,lhbcc,rhbcc,&
         topbcc,botbcc,lbvc,rbvc,tbvc,bbvc)
  call boundary(u,lhbdy,rhbdy,topbdy,botbdy,lhbcu,rhbcu,&
         topbcu,botbcu,lbvu,rbvu,tbvu,bbvu)
  call boundary(w,lhbdy,rhbdy,topbdy,botbdy,lhbcw,rhbcw,&
         topbcw,botbcw,lbvw,rbvw,tbvw,bbvw)

  call nbrex2d(w,.false.,.false.,wlfsl,wrfsl,myid,nbrleft,nbrright,&
       nbrtop,nbrbottom,comm2d)
  call nbrex2d(u,utfsl,ubfsl,.false.,.false.,myid,nbrleft,nbrright,&
       nbrtop,nbrbottom,comm2d)
  call nbrex2d(press,.false.,.false.,.false.,.false.,myid,nbrleft,&
       nbrright,nbrtop,nbrbottom,comm2d)
  call nbrex2d(temp,ttfsl,tbfsl,tlfsl,trfsl,myid,nbrleft,nbrright,&
       nbrtop,nbrbottom,comm2d)
  ! Tracer method with 0 diffusion
  call nbrex2d(chem,.false.,.false.,.false.,.false.,myid,nbrleft,nbrright,&
       nbrtop,nbrbottom,comm2d)
  call nbrex2d(alpha,.true.,.true.,.true.,.true.,myid,nbrleft,nbrright,&
       nbrtop,nbrbottom,comm2d)

  filenumber=0
  filenumber_tr=0
  trtime=0.0_wp
  mtime=0.0_wp
  open(21,status='unknown',file='data.t')
!!$  open(23,status='unknown',file='ilog')

!!!!!!!!!!!!!!!!!!!!TIME LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  time=0.0_wp
  do while (time.lt.time_max)

     ! Define viscosity 
!!$     ! With power law approximation
!!$     forall(i=1:mynx(ng),j=1:myny(ng))
!!$        mu(i,j) = exp((etae*(temp(i,j)+(1.0_wp-y(ng)%d(j))+etat)**mv)/&
!!$             ((etat**(mv+1))/mv+(temp(i,j)+(1.0_wp-y(ng)%d(j))+etat)**(mv+1))-&
!!$             (etae*(etat+t0)**mv)/((etat**(mv+1))/mv+(etat+t0)**(mv+1)))
!!$     end forall
!!$     if (mu_jump==1.0_wp) then
     do j=1,myny(ng)
        do i=1,mynx(ng)
           if (temp(i,j).lt.temp_cutoff) then
              mu(i,j) = exp(etae/(temp_cutoff+etat)-etae/(t0+etat))*&
                   alpha(i,j)**(-m_alpha)*(1.0_wp-chem(i,j)*&
                   (1.0_wp-mu_jump))
           else
              mu(i,j) = exp(etae/(temp(i,j)+etat)-etae/(t0+etat))*&
                   alpha(i,j)**(-m_alpha)*(1.0_wp-chem(i,j)*&
                   (1.0_wp-mu_jump))
           endif
        enddo
     enddo
     ! Only temperature dependent part of viscosity for use in damage eq
     do j=1,myny(ng)
        do i=1,mynx(ng)
           if (temp(i,j).lt.temp_cutoff) then
              mu_t(i,j) = exp(etae/(temp_cutoff+etat)-etae/(t0+etat))*&
                   (1.0_wp-chem(i,j)*(1.0_wp-mu_jump))
           else
              mu_t(i,j) = exp(etae/(temp(i,j)+etat)-etae/(t0+etat))*&
                   (1.0_wp-chem(i,j)*(1.0_wp-mu_jump))
           endif
        enddo
     enddo
!!$     else
!!$        ! Adding viscosity jump due to composition; this is for using the 
!!$        ! number of tracers in a cell to determine whether the viscosity jump
!!$        ! should be applied to that cell or not   
!!$        call calc_visc_tracer(tracers,mu,mu_t,alpha,temp,etae,etat,t0,&
!!$             temp_cutoff,m_alpha,mu_jump,myntracers,nx,ny,sxu,syw,asp,&
!!$             ntracers,ntotal)    
!!$     endif
     mu_c=interpgen2d(mu,x(ng)%d,y(ng)%d,xu(ng)%d,yw(ng)%d)

     ! CALCULATE COEFFICIENTS FOR SOLVING MOMENTUM EQUATIONS
     ! Coefficients for horizontal velocity, u
     allocate(aeu(ng),awu(ng),anu(ng),asu(ng))
     ! Call calc_a routine to calculate coefficents on each grid level
     call calc_au(ng,nglow,aeu,awu,anu,asu,mu,mynx,myny,mynxu,mynyw,&
          dxu,dyw,dyu,x,y,xu,yw,myid,nbrleft,nbrright,nbrtop,nbrbottom,&
          comm2d,lhbdy,lhbcu,rhbdy,rhbcu,topbdy,topbcu,botbdy,botbcu)

     ! Coefficients for vertical velocity, w
     allocate(aew(ng),aww(ng),anw(ng),asw(ng))
     ! Call calc_a routine to calculate coefficents on each grid level
     call calc_aw(ng,nglow,aew,aww,anw,asw,mu,mynx,myny,mynxu,mynyw,&
          dxu,dyw,dxw,x,y,xu,yw,myid,nbrleft,nbrright,nbrtop,nbrbottom,&
          comm2d,lhbdy,lhbcw,rhbdy,rhbcw)

     ! Now calculate coefficients for pressure equation
     allocate(ae_p(ng),aw_p(ng),an_p(ng),as_p(ng))
     call calc_ap(ae_p,aw_p,an_p,as_p,ng,nglow,mynx,myny,mynxu,mynyw,aeu,awu,&
          asu,anu,aew,aww,asw,anw,nbrleft,nbrright,nbrtop,nbrbottom,comm2d,&
          dxu,dyw)

!!$     ! Set relax and relax_p to 0.1 after 10 timesteps
!!$     if (itime.ge.10) then
!!$        relax=0.1_wp
!!$        relax_p=0.1_wp
!!$     endif

     ! SIMPLER algorithm for solving momentum equations
     do iter=1,iter_total

        press_last=press
        u_last_p=u
        w_last_p=w
        ! Set up rhs for w
        ! For tracer method with composition set on same grid as w
        forall (i=2:mynx(ng)-1,j=2:mynyw(ng)-1)
           w_rhs(i,j)=-Ra*((temp(i,j)+(temp(i,j+1)-temp(i,j))*&
                ((yw(ng)%d(j)-y(ng)%d(j))/(y(ng)%d(j+1)-y(ng)%d(j))))-&
                Bu*(chem(i,j)))-1.0_wp/(dy(ng)%d(j)*dxu(ng)%d(i-1))*&
                (mu_c(i,j)*(u(i,j+1)-u(i,j))&
                -mu_c(i-1,j)*(u(i-1,j+1)-u(i-1,j))) 
        end forall
        ! Calculate psuedovelocity, w_hat (eq 6.26 in Patankar) 
        forall (i=2:mynx(ng)-1,j=2:mynyw(ng)-1)
           w_hat(i,j)=(1.0_wp/(aew(ng)%a(i-1,j-1)+aww(ng)%a(i-1,j-1)+&
                anw(ng)%a(i-1,j-1)+asw(ng)%a(i-1,j-1)))*(aew(ng)%a(i-1,j-1)&
                *w(i+1,j)+aww(ng)%a(i-1,j-1)*w(i-1,j)+anw(ng)%a(i-1,j-1)*&
                w(i,j+1)+asw(ng)%a(i-1,j-1)*w(i,j-1)-&
                (dy(ng)%d(j)*dxu(ng)%d(i-1))*w_rhs(i,j))
        end forall
        ! BC's 
        call boundary(w_hat,lhbdy,rhbdy,topbdy,botbdy,lhbcw,rhbcw,&
         topbcw,botbcw,lbvw,rbvw,tbvw,bbvw)
        ! Neighbor exchange for w_hat
        call nbrex2d(w_hat,.false.,.false.,wlfsl,wrfsl,myid,nbrleft,&
             nbrright,nbrtop,nbrbottom,comm2d)
        ! Set up rhs for u
        forall (i=2:mynxu(ng)-1,j=2:myny(ng)-1)
           u_rhs(i,j)=-1.0_wp/(dyw(ng)%d(j-1)*dx(ng)%d(i))*&
                (mu_c(i,j)*(w(i+1,j)-w(i,j))&
                -mu_c(i,j-1)*(w(i+1,j-1)-w(i,j-1)))
        end forall
        ! Calculate psuedovelocity, u_hat (eq 6.26 in Patankar) 
        forall (i=2:mynxu(ng)-1,j=2:myny(ng)-1)
           u_hat(i,j)=1.0_wp/(aeu(ng)%a(i-1,j-1)+awu(ng)%a(i-1,j-1)+&
                anu(ng)%a(i-1,j-1)+asu(ng)%a(i-1,j-1))*(aeu(ng)%a(i-1,j-1)*&
                u(i+1,j)+awu(ng)%a(i-1,j-1)*u(i-1,j)+anu(ng)%a(i-1,j-1)*&
                u(i,j+1)+asu(ng)%a(i-1,j-1)*u(i,j-1)-&
                (dyw(ng)%d(j-1)*dx(ng)%d(i))*u_rhs(i,j))
        end forall
        ! BC's       
        call boundary(u_hat,lhbdy,rhbdy,topbdy,botbdy,lhbcu,rhbcu,&
         topbcu,botbcu,lbvu,rbvu,tbvu,bbvu)
        ! Neighbor exchange for u_hat
        call nbrex2d(u_hat,utfsl,ubfsl,.false.,.false.,myid,nbrleft,&
             nbrright,nbrtop,nbrbottom,comm2d)

        ! PRESSURE EQUATION
        ! Set up rhs using psuedovelocities
        forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
           p_rhs(i,j)=(1.0_wp/dxu(ng)%d(i-1))*(u_hat(i,j)-u_hat(i-1,j))+&
                (1.0_wp/dyw(ng)%d(j-1))*(w_hat(i,j)-w_hat(i,j-1))
        end forall
        ! Set arbitrary point in the pressure field to zero 
        ! (see Patankar 6.7-4 for explanation)
        if (myid==0) then
           x_zero=x(ng)%d(mynx(ng)-1)
           y_zero=y(ng)%d(myny(ng)-1)
        else
           x_zero=-100.0_wp
           y_zero=-100.0_wp
        endif
        ! Now call multigrid to solve pressure equation
        call mg2d(press,p_rhs,0.0_wp,ng,nglow,0,p_cyc,p_it,mynx,myny,x,y,&
             ae_p,aw_p,an_p,as_p,0.0_wp,myid,nbrleft,nbrright,nbrtop,&
             nbrbottom,.false.,.false.,.false.,.false.,x_zero,y_zero,&
             0.0_wp,comm2d,dxu,dyw,1,1)

        ! Under-relax pressure solution 
        press = relax_p*press + (1.0_wp-relax_p)*press_last

        ! Treat this pressure as p*, solve for u* and w* (as in Patankar, eqs. 
        ! 6.8-6.10)
        ! Iterate over u* and w* at this pressure since they are coupled 
        ! However, do not need to iterate to convergence since pressure field
        ! will be relcalculated at next overall iteration

        do loop=1,iter_uw_loop
           u_last=u
           w_last=w
           ! Set up rhs for w
           ! For tracer method with composition set on same grid as w
           forall (i=2:mynx(ng)-1,j=2:mynyw(ng)-1)
              w_rhs(i,j)=1.0_wp/(dy(ng)%d(j))*(press(i,j+1)-press(i,j))&
                   -Ra*((temp(i,j)+(temp(i,j+1)-temp(i,j))*&
                   ((yw(ng)%d(j)-y(ng)%d(j))/(y(ng)%d(j+1)-y(ng)%d(j))))-&
                   Bu*(chem(i,j)))-1.0_wp/(dy(ng)%d(j)*dxu(ng)%d(i-1))*&
                   (mu_c(i,j)*(u(i,j+1)-u(i,j))&
                   -mu_c(i-1,j)*(u(i-1,j+1)-u(i-1,j))) 
           end forall
           ! Call multigrid to solve for w* 
           call mg2d(w,w_rhs,0.0_wp,ng,nglow,2,uw_cyc,uw_it,mynx,mynyw,xw,yw,&
                aew,aww,anw,asw,0.0_wp,myid,nbrleft,nbrright,nbrtop,nbrbottom,&
                .false.,.false.,wlfsl,wrfsl,-100.0_wp,-100.0_wp,0.0_wp,comm2d,&
                dxu,dy,1,0)

           ! Now under-relax solution for w*
           w=relax*w + (1.0_wp-relax)*w_last
           ! Set up rhs for u
           forall (i=2:mynxu(ng)-1,j=2:myny(ng)-1)
              u_rhs(i,j)=(1.0_wp/dx(ng)%d(i))*(press(i+1,j)-press(i,j))&
                   -1.0_wp/(dyw(ng)%d(j-1)*dx(ng)%d(i))*&
                   (mu_c(i,j)*(w(i+1,j)-w(i,j))&
                   -mu_c(i,j-1)*(w(i+1,j-1)-w(i,j-1)))
           end forall
           ! Call multigrid to solve for u* 
           call mg2d(u,u_rhs,0.0_wp,ng,nglow,1,uw_cyc,uw_it,mynxu,myny,xu,yu,&
                aeu,awu,anu,asu,0.0_wp,myid,nbrleft,nbrright,nbrtop,nbrbottom,&
                utfsl,ubfsl,.false.,.false.,-100.0_wp,-100.0_wp,0.0_wp,comm2d,&
                dx,dyw,0,1)
           ! Now under-relax solution for u* 
           u=relax*u + (1.0_wp-relax)*u_last
           ! Now check if solutions for u and w have converged 
           tolu=tol*relax*(maxval(abs(u)))
           if (tolu.eq.0.0_wp) tolu=tol
           ! Combine tolerances of u from all processors
           call MPI_ALLREDUCE(tolu,tolu_tot,1,MPI_DOUBLE_PRECISION,&
                MPI_MAX,comm2d,ierr)
           erru=maxval(abs(u-u_last))
           call MPI_ALLREDUCE(erru,erru_tot,1,MPI_DOUBLE_PRECISION,&
                MPI_MAX,comm2d,ierr)
           tolw=tol*relax*(maxval(abs(w)))
           if (tolw.eq.0.0_wp) tolw=tol
           call MPI_ALLREDUCE(tolw,tolw_tot,1,MPI_DOUBLE_PRECISION,&
                MPI_MAX,comm2d,ierr)
           errw=maxval(abs(w-w_last))
           call MPI_ALLREDUCE(errw,errw_tot,1,MPI_DOUBLE_PRECISION,&
                MPI_MAX,comm2d,ierr)
           if (erru_tot<tolu_tot.and.errw_tot<tolw_tot) exit
        enddo

        ! Now solve for pressure correction 
        ! Set up rhs using u* and w*
        forall (i=2:mynx(ng)-1,j=2:myny(ng)-1)
           p_rhs(i,j)=(1.0_wp/dxu(ng)%d(i-1))*(u(i,j)-u(i-1,j))+&
                (1.0_wp/dyw(ng)%d(j-1))*(w(i,j)-w(i,j-1))
        end forall
        p_cor=0.0_wp
        ! Solve for pressure correction
        call mg2d(p_cor,p_rhs,0.0_wp,ng,nglow,0,uw_cyc,uw_it,mynx,myny,x,y,&
             ae_p,aw_p,an_p,as_p,0.0_wp,myid,nbrleft,nbrright,nbrtop,&
             nbrbottom,.false.,.false.,.false.,.false.,-100.0_wp,-100.0_wp,&
             0.0_wp,comm2d,dxu,dyw,1,1)
!!$        call mglin2d(p_cor,p_rhs,0.0_wp,ng,nglow,0,uw_cyc,uw_it,mynx,myny,x,y,&
!!$             ae_p,aw_p,an_p,as_p,0.0_wp,myid,nbrleft,nbrright,nbrtop,&
!!$             nbrbottom,.false.,.false.,.false.,.false.,-100.0_wp,-100.0_wp,&
!!$             0.0_wp,comm2d)
        ! Now correct velocities (equations 6.17 and 6.19 in Patankar)
        forall (i=2:mynxu(ng)-1,j=2:myny(ng)-1)
           u(i,j)=u(i,j)+(ae_p(ng)%a(i-1,j-1)/dyw(ng)%d(j-1))*&
                (p_cor(i,j)-p_cor(i+1,j))
        end forall
        call nbrex2d(u,utfsl,ubfsl,.false.,.false.,myid,nbrleft,&
             nbrright,nbrtop,nbrbottom,comm2d)
        forall (i=2:mynx(ng)-1,j=2:mynyw(ng)-1)
           w(i,j)=w(i,j)+(an_p(ng)%a(i-1,j-1)/dxu(ng)%d(i-1))*&
                (p_cor(i,j)-p_cor(i,j+1))
        end forall
        call nbrex2d(w,.false.,.false.,wlfsl,wrfsl,myid,nbrleft,&
             nbrright,nbrtop,nbrbottom,comm2d)

        if (remove_net_u) then
           ! Remove net horizontal velocity 
           call average_u(u,dyw(ng)%d,dxu(ng)%d,asp,u0,mynx,myny,ng)
           u=u-u0
        endif

        ! Check if velocities and pressure have converged 
        tolu=tol*relax*(maxval(abs(u)))
        if (tolu.eq.0.0_wp) tolu=tol*relax
        call MPI_ALLREDUCE(tolu,tolu_tot,1,MPI_DOUBLE_PRECISION,&
             MPI_MAX,comm2d,ierr)
        erru=maxval(abs(u-u_last_p))
        call MPI_ALLREDUCE(erru,erru_tot,1,MPI_DOUBLE_PRECISION,&
             MPI_MAX,comm2d,ierr)
        tolw=tol*relax*(maxval(abs(w)))
        if (tolw.eq.0.0_wp) tolw=tol*relax
        call MPI_ALLREDUCE(tolw,tolw_tot,1,MPI_DOUBLE_PRECISION,&
             MPI_MAX,comm2d,ierr)
        errw=maxval(abs(w-w_last_p))
        call MPI_ALLREDUCE(errw,errw_tot,1,MPI_DOUBLE_PRECISION,&
             MPI_MAX,comm2d,ierr)
        tolp=tol*relax_p*(maxval(abs(press)))
        if (tolp.eq.0.0) tolp=tol*relax_p
        call MPI_ALLREDUCE(tolp,tolp_tot,1,MPI_DOUBLE_PRECISION,& 
             MPI_MAX,comm2d,ierr)
        errp=maxval(abs(press-press_last))
        call MPI_ALLREDUCE(errp,errp_tot,1,MPI_DOUBLE_PRECISION,&
             MPI_MAX,comm2d,ierr)
!!$        errv=maxval(abs(dxu(ng)%d(1)*p_rhs))
!!$        call MPI_ALLREDUCE(errv,errv_tot,1,MPI_DOUBLE_PRECISION,&
!!$             MPI_MAX,comm2d,ierr)
!!$        tolv_tot=(1.0_wp/relax)*max(tolu_tot,tolw_tot)
        errv=maxval(abs(p_rhs))
        call MPI_ALLREDUCE(errv,errv_tot,1,MPI_DOUBLE_PRECISION,&
             MPI_MAX,comm2d,ierr)
        tolv=max(tolu_tot/maxval(dxu(ng)%d),tolw_tot/maxval(dyw(ng)%d))
        call MPI_ALLREDUCE(tolv,tolv_tot,1,MPI_DOUBLE_PRECISION,&
             MPI_MAX,comm2d,ierr)
        tolv_tot=tolv_tot*(1.0_wp/relax) 

!!$        if (myid==0) then
!!$           print *, 'iter # =', iter
!!$           print *, 'tol u =',tolu_tot,'err u =',erru_tot
!!$           print *, 'tol w =',tolw_tot,'err w =',errw_tot
!!$           print *, 'tol p =',tolp_tot,'err p =',errp_tot
!!$           print *, 'tol v =',tolv_tot,'err v=',errv_tot
!!$        endif

        if (erru_tot<tolu_tot.and.errw_tot<tolw_tot.and.errp_tot<tolp_tot&
             .and.errv_tot<tolv_tot) exit
        if (erru_tot<1e-8.and.errw_tot<1e-8.and.errp_tot<1e-6&
             .and.errv_tot<1e-6) exit
!!$        if (erru_tot<tolu_tot.and.errw_tot<tolw_tot.and.errp_tot<tolp_tot) exit
!!$        if (erru_tot<1e-9.and.errw_tot<1e-9.and.errp_tot<1e-7) exit
        if (iter.eq.iter_total) then
           print*, 'Solution for u and w has not converged'
           print*, 'Program terminated'
           STOP
        endif
     enddo

     do j=nglow,ng
        deallocate(aeu(j)%a,awu(j)%a,anu(j)%a,asu(j)%a,aew(j)%a,aww(j)%a,&
             anw(j)%a,asw(j)%a,ae_p(j)%a,aw_p(j)%a,an_p(j)%a,as_p(j)%a)
     enddo
     deallocate(aeu,awu,anu,asu,aew,aww,anw,asw,ae_p,aw_p,an_p,as_p)

     ! Now calculate the stress components, tau_xx, tau_zz, tau_xz 
     ! First calculate strain-rate 
     call strainrate(u,w,dxw(ng)%d,dyu(ng)%d,dxu(ng)%d,dyw(ng)%d,du_dx,du_dy,&
          dw_dx,dw_dy)
     ! Now determine stress components
     forall (i=1:mynx(ng)-2,j=1:myny(ng)-2)
        tau_xx(i+1,j+1)=2.0_wp*mu(i+1,j+1)*du_dx(i,j)
        tau_zz(i+1,j+1)=2.0_wp*mu(i+1,j+1)*dw_dy(i,j)
        tau_xz(i+1,j+1)=mu(i+1,j+1)*(du_dy(i,j)+dw_dx(i,j))
        tau(i+1,j+1)=sqrt(0.5_wp*tau_xx(i+1,j+1)**2.0_wp +&
             0.5_wp*tau_zz(i+1,j+1)**2.0_wp + tau_xz(i+1,j+1)**2.0_wp)
        dissipation(i+1,j+1)=tau_xx(i+1,j+1)*du_dx(i,j)+&
             tau_zz(i+1,j+1)*dw_dy(i,j)+tau_xz(i+1,j+1)*(du_dy(i,j)+dw_dx(i,j))
     endforall

!!!!!!!!!!!!!!!!!!!!!WRITE OUT SOLUTION FOR PLOTTING!!!!!!!!!!!!!!!!!!!!!!!!!! 

     if (time.eq.0.0_wp) then
        ! First use mpi_gather to recombine temperature, composition, 
        ! fineness, and viscosity fields onto root process
        ! Root process allocates space for combined solutions
        if (myid==0) then
           allocate(temp_tot(nx,ny),alpha_tot(nx,ny),mu_tot(nx,ny),&
                press_tot(nx,ny),u_tot(nx+1,ny),w_tot(nx,ny+1),&
                temp_ha(ny),chem_ha(ny),alpha_ha(ny),mu_ha(ny),&
                tau_tot(nx,ny),tauxx_tot(nx,ny),tauzz_tot(nx,ny),&
                tauxz_tot(nx,ny),dissipation_tot(nx,ny))
           allocate(chem_tot(nx,ny+1))
        endif

        ! Now use subroutines combine_root and combine_rest to recombine
        ! a decomposed array
        if (myid==0) then
           call combine_root(temp,temp_tot,maxrow,nx,mynx(ng),myny(ng),&
                numprocs,sxall,syall,mynxall,mynyall,comm2d)
           call combine_root(press,press_tot,maxrow,nx,mynx(ng),myny(ng),&
                numprocs,sxall,syall,mynxall,mynyall,comm2d)
           call combine_root(alpha,alpha_tot,maxrow,nx,mynx(ng),myny(ng),&
                numprocs,sxall,syall,mynxall,mynyall,comm2d)
           call combine_root(mu,mu_tot,maxrow,nx,mynx(ng),myny(ng),numprocs,&
                sxall,syall,mynxall,mynyall,comm2d)
           call combine_root(tau,tau_tot,maxrow,nx,mynx(ng),myny(ng),&
                numprocs,sxall,syall,mynxall,mynyall,comm2d)
           call combine_root(tau_xx,tauxx_tot,maxrow,nx,mynx(ng),myny(ng),&
                numprocs,sxall,syall,mynxall,mynyall,comm2d)
           call combine_root(tau_zz,tauzz_tot,maxrow,nx,mynx(ng),myny(ng),&
                numprocs,sxall,syall,mynxall,mynyall,comm2d)
           call combine_root(tau_xz,tauxz_tot,maxrow,nx,mynx(ng),myny(ng),&
                numprocs,sxall,syall,mynxall,mynyall,comm2d)
           call combine_root(dissipation,dissipation_tot,maxrow,nx,mynx(ng),&
                myny(ng),numprocs,sxall,syall,mynxall,mynyall,comm2d)
        else 
           call combine_rest(temp,maxrow,mynx(ng),myny(ng),comm2d)
           call combine_rest(press,maxrow,mynx(ng),myny(ng),comm2d)
           call combine_rest(alpha,maxrow,mynx(ng),myny(ng),comm2d)
           call combine_rest(mu,maxrow,mynx(ng),myny(ng),comm2d)
           call combine_rest(tau,maxrow,mynx(ng),myny(ng),comm2d)
           call combine_rest(tau_xx,maxrow,mynx(ng),myny(ng),comm2d)
           call combine_rest(tau_zz,maxrow,mynx(ng),myny(ng),comm2d)
           call combine_rest(tau_xz,maxrow,mynx(ng),myny(ng),comm2d)
           call combine_rest(dissipation,maxrow,mynx(ng),myny(ng),comm2d)
        endif

        ! Now recombine u 
        if (myid==0) then
           call combine_root(u,u_tot,maxrow,nx+1,mynxu(ng),myny(ng),numprocs,&
                sxuall,syall,mynxuall,mynyall,comm2d)
        else
           call combine_rest(u,maxrow,mynxu(ng),myny(ng),comm2d) 
        endif
        if (myid==0) then
           ! Make sure boundary conditions are enforced on u
           call boundary(u_tot,.true.,.true.,.false.,.false.,lhbcu,rhbcu,&
                topbcu,botbcu,lbvu,rbvu,tbvu,bbvu)
           if (lhbcu=='per') then
              u_tot(1,1:ny)=u_tot(nx+1,1:ny)
           endif
        endif
 
        ! Now recombine w
        if (myid==0) then
           call combine_root(w,w_tot,maxroww,nx,mynx(ng),mynyw(ng),&
                numprocs,sxall,sywall,mynxall,mynywall,comm2d)
        else
           call combine_rest(w,maxroww,mynx(ng),mynyw(ng),comm2d) 
        endif
        if (myid==0) then
           ! Make sure boundary conditions are enforced on w
           call boundary(w_tot,.false.,.false.,.true.,.true.,lhbcw,rhbcw,&
                topbcw,botbcw,lbvw,rbvw,tbvw,bbvw)
        endif

        ! Now recombine chem for tracer method where grid is same as w
        if (myid==0) then
           call combine_root(chem,chem_tot,maxroww,nx,mynx(ng),mynyw(ng),&
                numprocs,sxall,sywall,mynxall,mynywall,comm2d)
           call combine_root_bdy(chem,chem_tot,mynx(ng),mynyw(ng),numprocs,&
                sxall,sywall,eywall,mynxall,topbdy,botbdy,comm2d)
        else
           call combine_rest(chem,maxroww,mynx(ng),mynyw(ng),comm2d) 
           call combine_rest_bdy(chem,topbdy,botbdy,mynx(ng),mynyw(ng),comm2d)
        endif

        ! Now root process prints out combined files
        if (myid==0) then
           write(dname,"('f000')")
           open(10,file=dname)
           do i=1,nx
              do j=1,ny
                 write(10,'(7(1pe11.4,1x))') x_tot(i),y_tot(j),temp_tot(i,j),&
                      alpha_tot(i,j),press_tot(i,j),mu_tot(i,j),&
                      dissipation_tot(i,j)
              enddo
           enddo
           close(10)
           write(tname,"('t000')")
           open(10,file=tname)
           do i=1,nx
              do j=1,ny
                 write(10,'(6(1pe11.4,1x))') x_tot(i),y_tot(j),tau_tot(i,j),&
                      tauxx_tot(i,j),tauzz_tot(i,j),tauxz_tot(i,j)
              enddo
           enddo
           close(10)
           write(uname,"('u000')") 
           open(10,file=uname)
           do i=1,nx+1
              do j=1,ny
                 write(10,'(3(1pe11.4,1x))') xu_tot(i),y_tot(j),u_tot(i,j)
              enddo
           enddo
           close(10)
           write(wname,"('w000')") 
           open(10,file=wname)
           do i=1,nx
              do j=1,ny+1
                 write(10,'(4(1pe11.4,1x))') x_tot(i),yw_tot(j),w_tot(i,j),&
                      chem_tot(i,j)
              enddo
           enddo
           close(10)
        endif
        ! Write out tracer position info
        if (ntracers>0) then
           if (myid==0) then 
              allocate(tracers_tot(ntracers,3))
              ! Use mpi_gather to put all tracer data on root process
              call gather_tracers_root(tracers,tracers_tot,myntracers,&
                   ntracers,comm2d)
           else
              ! Use mpi_gather to put all tracer data on root process
              call gather_tracers_rest(tracers,myntracers,comm2d)
           endif
           if (myid==0) then
              write(trname,"('tr000')")
              open(10,file=trname)
              do i=1,ntracers
                 write(10,'(4(1pe11.4,1x))') time+t_off,tracers_tot(i,1),&
                   tracers_tot(i,2),tracers_tot(i,3)
              enddo
              close(10)
              deallocate(tracers_tot)
           endif
        endif
        
        ! Compute horizontal averages 
        if (myid==0) then
           call ha_root(temp,temp_ha,dxu(ng)%d,mynx(ng),myny(ng),sxall,syall,&
                mynyall,npx,asp,numprocs,maxrow,comm2d)
           call ha_root_chem(chem,chem_ha,dxu(ng)%d,mynx(ng),myny(ng),sxall,&
                syall,mynyall,npx,asp,numprocs,comm2d)
           call ha_root(alpha,alpha_ha,dxu(ng)%d,mynx(ng),myny(ng),sxall,syall,&
                mynyall,npx,asp,numprocs,maxrow,comm2d)
           call ha_root(mu,mu_ha,dxu(ng)%d,mynx(ng),myny(ng),sxall,syall,&
           mynyall,npx,asp,numprocs,maxrow,comm2d)
        else
           call ha_rest(temp,dxu(ng)%d,mynx(ng),myny(ng),maxrow,comm2d)
           call ha_rest_chem(chem,dxu(ng)%d,mynx(ng),myny(ng),comm2d)
           call ha_rest(alpha,dxu(ng)%d,mynx(ng),myny(ng),maxrow,comm2d)
           call ha_rest(mu,dxu(ng)%d,mynx(ng),myny(ng),maxrow,comm2d)
        endif
        
        ! Now write out ha file
        if (myid==0) then
           write(hname,"('ha000')")
           open(10,file=hname)
              do j=1,ny
                 write(10,'(6(1pe11.4,1x))') time+t_off,y_tot(j),temp_ha(j),&
                      chem_ha(j),alpha_ha(j),mu_ha(j)
              enddo
           close(10)
        endif

        mtime=time+mstep
        trtime = time + trstep
        if (myid==0) then
           deallocate(temp_tot,chem_tot,alpha_tot,mu_tot,press_tot,&
                u_tot,w_tot,temp_ha,chem_ha,alpha_ha,mu_ha,tau_tot,&
                tauxx_tot,tauzz_tot,tauxz_tot,dissipation_tot)
        endif
     endif

     if (time >= mtime) then
        mtime = time + mstep
        filenumber = filenumber + 1
        if (filenumber.gt.file_max) then
           print *, "Too many files for movie, Abort!!!"
           STOP
        endif

        ! First use mpi_gather to recombine temperature, composition, 
        ! fineness, and viscosity fields onto root process
        ! Root process allocates space for combined solutions
        if (myid==0) then
           allocate(temp_tot(nx,ny),alpha_tot(nx,ny),&
                mu_tot(nx,ny),press_tot(nx,ny),u_tot(nx+1,ny),w_tot(nx,ny+1),&
                temp_ha(ny),chem_ha(ny),alpha_ha(ny),mu_ha(ny),&
                tau_tot(nx,ny),tauxx_tot(nx,ny),tauzz_tot(nx,ny),&
                tauxz_tot(nx,ny),dissipation_tot(nx,ny))
           allocate(chem_tot(nx,ny+1))
        endif

        ! Now use subroutines combine_root and combine_rest to recombine
        ! a decomposed array
        if (myid==0) then
           call combine_root(temp,temp_tot,maxrow,nx,mynx(ng),myny(ng),&
                numprocs,sxall,syall,mynxall,mynyall,comm2d)
           call combine_root(chem,chem_tot,maxrow,nx,mynx(ng),myny(ng),&
                numprocs,sxall,syall,mynxall,mynyall,comm2d)
           call combine_root(press,press_tot,maxrow,nx,mynx(ng),myny(ng),&
                numprocs,sxall,syall,mynxall,mynyall,comm2d)
           call combine_root(alpha,alpha_tot,maxrow,nx,mynx(ng),myny(ng),&
                numprocs,sxall,syall,mynxall,mynyall,comm2d)
           call combine_root(mu,mu_tot,maxrow,nx,mynx(ng),myny(ng),numprocs,&
                sxall,syall,mynxall,mynyall,comm2d)
           call combine_root(tau,tau_tot,maxrow,nx,mynx(ng),myny(ng),&
                numprocs,sxall,syall,mynxall,mynyall,comm2d)
           call combine_root(tau_xx,tauxx_tot,maxrow,nx,mynx(ng),myny(ng),&
                numprocs,sxall,syall,mynxall,mynyall,comm2d)
           call combine_root(tau_zz,tauzz_tot,maxrow,nx,mynx(ng),myny(ng),&
                numprocs,sxall,syall,mynxall,mynyall,comm2d)
           call combine_root(tau_xz,tauxz_tot,maxrow,nx,mynx(ng),myny(ng),&
                numprocs,sxall,syall,mynxall,mynyall,comm2d)
           call combine_root(dissipation,dissipation_tot,maxrow,nx,mynx(ng),&
                myny(ng),numprocs,sxall,syall,mynxall,mynyall,comm2d)
        else 
           call combine_rest(temp,maxrow,mynx(ng),myny(ng),comm2d)
           call combine_rest(chem,maxrow,mynx(ng),myny(ng),comm2d)
           call combine_rest(press,maxrow,mynx(ng),myny(ng),comm2d)
           call combine_rest(alpha,maxrow,mynx(ng),myny(ng),comm2d)
           call combine_rest(mu,maxrow,mynx(ng),myny(ng),comm2d)
           call combine_rest(tau,maxrow,mynx(ng),myny(ng),comm2d)
           call combine_rest(tau_xx,maxrow,mynx(ng),myny(ng),comm2d)
           call combine_rest(tau_zz,maxrow,mynx(ng),myny(ng),comm2d)
           call combine_rest(tau_xz,maxrow,mynx(ng),myny(ng),comm2d)
           call combine_rest(dissipation,maxrow,mynx(ng),myny(ng),comm2d)
        endif

        ! Now recombine u 
        if (myid==0) then
           call combine_root(u,u_tot,maxrow,nx+1,mynxu(ng),myny(ng),numprocs,&
                sxuall,syall,mynxuall,mynyall,comm2d)
        else
           call combine_rest(u,maxrow,mynxu(ng),myny(ng),comm2d) 
        endif
        if (myid==0) then
           ! Make sure boundary conditions are enforced on u
           call boundary(u_tot,.true.,.true.,.false.,.false.,lhbcu,rhbcu,&
                topbcu,botbcu,lbvu,rbvu,tbvu,bbvu)
           if (lhbcu=='per') then
              u_tot(1,1:ny)=u_tot(nx+1,1:ny)
           endif
        endif
           
        ! Now recombine w
        if (myid==0) then
           call combine_root(w,w_tot,maxroww,nx,mynx(ng),mynyw(ng),numprocs,&
                sxall,sywall,mynxall,mynywall,comm2d)
        else
           call combine_rest(w,maxroww,mynx(ng),mynyw(ng),comm2d) 
        endif
        if (myid==0) then
           ! Make sure boundary conditions are enforced on w
           call boundary(w_tot,.false.,.false.,.true.,.true.,lhbcw,rhbcw,&
                topbcw,botbcw,lbvw,rbvw,tbvw,bbvw)
        endif

        ! Now recombine chem for tracer method where grid is same as w
        if (myid==0) then
           call combine_root(chem,chem_tot,maxroww,nx,mynx(ng),mynyw(ng),&
                numprocs,sxall,sywall,mynxall,mynywall,comm2d)
           call combine_root_bdy(chem,chem_tot,mynx(ng),mynyw(ng),numprocs,&
                sxall,sywall,eywall,mynxall,topbdy,botbdy,comm2d)
        else
           call combine_rest(chem,maxroww,mynx(ng),mynyw(ng),comm2d) 
           call combine_rest_bdy(chem,topbdy,botbdy,mynx(ng),mynyw(ng),comm2d)
        endif

        ! Root process writes output
        if (myid==0) then
           if (filenumber.lt.10) then
              write(dname,"('f','00',i1)") filenumber
              write(uname,"('u','00',i1)") filenumber
              write(wname,"('w','00',i1)") filenumber
              write(hname,"('ha','00',i1)") filenumber
              write(tname,"('t','00',i1)") filenumber
           elseif ((filenumber.ge.10).and.(filenumber.le.99)) then
              write(dname,"('f','0',i2)") filenumber
              write(uname,"('u','0',i2)") filenumber
              write(wname,"('w','0',i2)") filenumber
              write(hname,"('ha','0',i2)") filenumber
              write(tname,"('t','0',i2)") filenumber
           else
              write(dname,"('f',i3)") filenumber
              write(uname,"('u',i3)") filenumber
              write(wname,"('w',i3)") filenumber
              write(hname,"('ha',i3)") filenumber
              write(tname,"('t',i3)") filenumber
           endif
           open(10,file=dname)
           do i=1,nx
              do j=1,ny
                 write(10,'(7(1pe11.4,1x))') x_tot(i),y_tot(j),temp_tot(i,j),&
                      alpha_tot(i,j),press_tot(i,j),mu_tot(i,j),&
                      dissipation_tot(i,j)
              enddo
           enddo
           close(10)
           open(10,file=tname)
           do i=1,nx
              do j=1,ny
                 write(10,'(6(1pe11.4,1x))') x_tot(i),y_tot(j),tau_tot(i,j),&
                      tauxx_tot(i,j),tauzz_tot(i,j),tauxz_tot(i,j)
              enddo
           enddo
           close(10)
           open(10,file=uname)
           do i=1,nx+1
              do j=1,ny
                 write(10,'(3(1pe11.4,1x))') xu_tot(i),y_tot(j),u_tot(i,j)
              enddo
           enddo
           close(10)
           open(10,file=wname)
           do i=1,nx
              do j=1,ny+1
                 write(10,'(4(1pe11.4,1x))') x_tot(i),yw_tot(j),w_tot(i,j),&
                      chem_tot(i,j)
              enddo
           enddo
           close(10)
        endif

        ! Compute horizontal averages 
        if (myid==0) then
           call ha_root(temp,temp_ha,dxu(ng)%d,mynx(ng),myny(ng),sxall,syall,&
                mynyall,npx,asp,numprocs,maxrow,comm2d)
           call ha_root_chem(chem,chem_ha,dxu(ng)%d,mynx(ng),myny(ng),sxall,&
                syall,mynyall,npx,asp,numprocs,comm2d)
           call ha_root(alpha,alpha_ha,dxu(ng)%d,mynx(ng),myny(ng),sxall,syall,&
                mynyall,npx,asp,numprocs,maxrow,comm2d)
           call ha_root(mu,mu_ha,dxu(ng)%d,mynx(ng),myny(ng),sxall,syall,&
                mynyall,npx,asp,numprocs,maxrow,comm2d)
        else
           call ha_rest(temp,dxu(ng)%d,mynx(ng),myny(ng),maxrow,comm2d)
           call ha_rest_chem(chem,dxu(ng)%d,mynx(ng),myny(ng),comm2d)
           call ha_rest(alpha,dxu(ng)%d,mynx(ng),myny(ng),maxrow,comm2d)
           call ha_rest(mu,dxu(ng)%d,mynx(ng),myny(ng),maxrow,comm2d)
        endif
        
        ! Now write out ha file
        if (myid==0) then
           open(10,file=hname)
              do j=1,ny
                 write(10,'(6(1pe11.4,1x))') time+t_off,y_tot(j),temp_ha(j),&
                      chem_ha(j),alpha_ha(j),mu_ha(j)
              enddo
           close(10)
        endif

        ! Deallocate output arrays
        if (myid==0) then
           deallocate(temp_tot,chem_tot,alpha_tot,mu_tot,press_tot,&
                u_tot,w_tot,temp_ha,chem_ha,alpha_ha,mu_ha,tau_tot,&
                tauxx_tot,tauzz_tot,tauxz_tot,dissipation_tot)
        endif
     endif
     ! So I can write out tracer position files on a different time step 
     if (time >= trtime) then
        trtime = time + trstep
        filenumber_tr = filenumber_tr + nint(trstep/mstep)
        if (filenumber_tr.gt.file_max) then
           print *, "Too many files for movie, Abort!!!"
           STOP
        endif
        ! Write out tracer positions, if the number of tracers is greater 
        ! than 0
        if (ntracers>0) then
           ! Print out tracer data half as frequently to save disk space
           if (myid==0) then 
              allocate(tracers_tot(ntracers,3))
              ! Use mpi_gather to put all tracer data on root process
              call gather_tracers_root(tracers,tracers_tot,myntracers,&
                   ntracers,comm2d)
           else
              ! Use mpi_gather to put all tracer data on root process
              call gather_tracers_rest(tracers,myntracers,comm2d)
           endif
           if (myid==0) then
              if (filenumber_tr.lt.10) then
                 write(trname,"('tr','00',i1)") filenumber_tr
              elseif ((filenumber_tr.ge.10).and.(filenumber_tr.le.99)) then
                 write(trname,"('tr','0',i2)") filenumber_tr
              else
                 write(trname,"('tr',i3)") filenumber_tr
              endif
              open(10,file=trname)
              do i=1,ntracers
                 write(10,'(4(1pe11.4,1x))') time+t_off,tracers_tot(i,1),&
                   tracers_tot(i,2),tracers_tot(i,3)
              enddo
              close(10)
              deallocate(tracers_tot)
           endif
        endif
     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!WRITE OUT TIME SERIES INFO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (mod(itime,ntdump) == 0.or.time.eq.0.0_wp) then
        call nusselt(temp(:,myny(ng)-1),dy(ng)%d(myny(ng)-1),dxu(ng)%d,asp,&
             nu,mynx,ng,topbdy)
        call usurf(u(:,myny(ng)-1),dxu(ng)%d,asp,u_surf,mynxu,ng,topbdy)
        if (topbdy) then
           umax=maxval(u(:,myny(ng)-1))
        else 
           umax=-100000000000.0_wp
        endif
        call MPI_REDUCE(umax,umax_tot,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
             comm2d,ierr)
        if (topbdy) then
           umin=minval(u(:,myny(ng)-1))
        else
           umin=100000000000.0_wp
        endif
        call MPI_REDUCE(umin,umin_tot,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,&
             comm2d,ierr)
        call rms(u,w,dyw(ng)%d,dxu(ng)%d,asp,v_rms,mynx,myny,ng)
        ! Find the global maximum of u, w (and which is bigger of those 
        ! two)
        u_max=maxval(abs(u))
        w_max=maxval(abs(w))
        max_vel=max(u_max,w_max)
        call MPI_REDUCE(max_vel,max_vel_tot,1,MPI_DOUBLE_PRECISION,&
             MPI_MAX,0,comm2d,ierr)        
        ! Max fineness
        alpha_max=maxval(alpha)
        call MPI_REDUCE(alpha_max,alpha_maxtot,1,MPI_DOUBLE_PRECISION,&
             MPI_MAX,0,comm2d,ierr)   
        ! Max dissipation
        diss_max=maxval(abs(dissipation))
        call MPI_REDUCE(diss_max,diss_maxtot,1,MPI_DOUBLE_PRECISION,&
             MPI_MAX,0,comm2d,ierr) 
        ! Compute global average of tau, rms average of stress components 
        call average(tau,1.0_wp,dyw(ng)%d,dxu(ng)%d,asp,tau_avg,mynx,myny,ng)
        call average(tau_xx,2.0_wp,dyw(ng)%d,dxu(ng)%d,asp,tauxx_avg,mynx,myny,ng)
        !call average(tau_zz,2.0_wp,deltay,deltax,asp,tauzz_avg,mynx,myny,ng)
        call average(tau_xz,2.0_wp,dyw(ng)%d,dxu(ng)%d,asp,tauxz_avg,mynx,myny,ng)
        ! Find global max of tau_xx
        tau_xx_max=maxval(abs(tau_xx))
        call MPI_REDUCE(tau_xx_max,tauxx_maxtot,1,MPI_DOUBLE_PRECISION,&
             MPI_MAX,0,comm2d,ierr)   
        ! Compute global averages of tempertaure, alpha
        call average(temp,1.0_wp,dyw(ng)%d,dxu(ng)%d,asp,temp_avg,mynx,myny,ng)
        call average(alpha,1.0_wp,dyw(ng)%d,dxu(ng)%d,asp,alpha_avg,mynx,myny,ng)
        ! Calculate average chemical flux in cold downwellings 
        if (crust_production.or.chem_layer_init) then
           call chem_flux_avg(chem,temp,cfluxavg,dxu(ng)%d,dyw(ng)%d,&
                mynx,myny,ng)
           ! Calculate average and max depth of c=0.5 contour
           call crust_depth(chem,avg_crustdepth,max_crustdepth,mynx,&
                mynyw,ng)
        else
           cfluxavg=0.0_wp; avg_crustdepth=0.0_wp; max_crustdepth=0.0_wp
        endif
        if (myid==0) then
           write(21,'(16(1pe11.4,1x))') time+t_off,nu,temp_avg,v_rms,&
                abs(umax_tot-umin_tot)/2.0_wp,max_vel_tot,alpha_maxtot,&
                alpha_avg,tau_avg,tauxx_avg,tauxx_maxtot,tauxz_avg,cfluxavg,&
                avg_crustdepth,1.0_wp-max_crustdepth,diss_maxtot
        endif
        itime = itime + 1
     else
        itime = itime + 1
     endif

     dtd=minval(dxu(ng)%d)*minval(dyw(ng)%d)/2.0_wp
     call MPI_ALLREDUCE(dtd,dtd_tot,1,MPI_DOUBLE_PRECISION,&
          MPI_MIN,comm2d,ierr)
     dtc=(minval(dxu(ng)%d)*minval(dyw(ng)%d)*lewis)/2.0_wp
     call MPI_ALLREDUCE(dtc,dtc_tot,1,MPI_DOUBLE_PRECISION,&
          MPI_MIN,comm2d,ierr)
     ! This is conservative/restrictive since highest velocity may not
     ! coencide with smallest grid cell size
     dta=min(minval(dxu(ng)%d)/maxval(abs(u)),&
          minval(dyw(ng)%d)/maxval(abs(w)))
     call MPI_ALLREDUCE(dta,dta_tot,1,MPI_DOUBLE_PRECISION,&
             MPI_MIN,comm2d,ierr)
!!$     dt=min(dtd_tot,dtc_tot,dta_tot)*timestep
     dt=dta_tot*timestep

!!!!!!!!!!!!!!!!!!!!!!Tracer Advection!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (ntracers>0) then
        my_tracers_out=0
        my_tracers_in=0
        ! Set the arrays of tracer numbers leaving my domain to -1; -1 will 
        ! signify the end of the list of tracers leaving my domain in a 
        ! particular direction
        tr_out_left=-1; tr_out_right=-1; tr_out_bot=-1; tr_out_top=-1
        tr_out_topright=-1; tr_out_topleft=-1; tr_out_botright=-1; 
        tr_out_botleft=-1
        
        ! Processors without tracers will just allocate tracers(0), which is
        ! like not allocating it anyway (they also will skip the calculations
        ! because loops would go from 1 to 0, and thus not happen
        !
        
        ! If side boundaries are periodic, each processor needs and extra 
        ! set of ghost points on it's left hand side in order to deal with 
        ! tracers that may cross out of the domain to the left
        !
        ! Also need a new xu that includes the x-coordinate of this extra 
        ! ghost point
        if (rhbcu=='per'.and.lhbcu=='per') then
           allocate(u_tmp(mynxu(ng)+1,myny(ng)))
           allocate(xu_tmp(mynxu(ng)+1))
           do i=1,mynxu(ng)
              xu_tmp(i+1)=xu(ng)%d(i)
              do j=1,myny(ng)
                 u_tmp(i+1,j)=u(i,j)
              enddo
           enddo
           ! Now fill extra left side ghost points by doing a neighbor 
           ! exchange; everyone passes to the right and receives from 
           ! the left
           call MPI_SENDRECV(u(mynxu(ng)-2,1:myny(ng)),myny(ng),&
                MPI_DOUBLE_PRECISION,nbrright,0,u_tmp(1,1:myny(ng)),&
                myny(ng),MPI_DOUBLE_PRECISION,nbrleft,0,comm2d,status,ierr)
           call MPI_SENDRECV(dxu(mynxu(ng)-2),1,MPI_DOUBLE_PRECISION,nbrright,&
                1,deltax_tmp,1,MPI_DOUBLE_PRECISION,nbrleft,1,comm2d,&
                status,ierr)
           xu_tmp(1)=xu_tmp(2)-deltax_tmp
        endif
        
        ! Allocate space for old tracers (tracer positions after this 
        ! upcoming advection step) because I need this to remake the array
        ! of tracer data
        allocate(oldtracers(myntracers))
        do i=1,myntracers
           allocate(oldtracers(i)%d(3))
        enddo
        count_tr=1; count_tl=1; count_br=1; count_bl=1
        count_t=1; count_l=1; count_r=1; count_b=1

        do i=1,myntracers
           if (rhbcu=='per'.and.lhbcu=='per') then
              call tracer_advec(xnew,ynew,tracers(i)%d(2),tracers(i)%d(3),dt,&
                   u_tmp,w,asp,x(ng)%d,xu_tmp,y(ng)%d,yw(ng)%d,nx,ny,sx,sxu,&
                   sy,syw,.true.)
           else
              call tracer_advec(xnew,ynew,tracers(i)%d(2),tracers(i)%d(3),dt,&
                   u,w,asp,x(ng)%d,xu(ng)%d,y(ng)%d,yw(ng)%d,nx,ny,sx,sxu,&
                   sy,syw,.false.)
           endif
           tracers(i)%d(2)=xnew
           tracers(i)%d(3)=ynew 
           oldtracers(i)%d=tracers(i)%d

           ! Now figure out how many tracers left my domain, 
           ! and where they went
           !
           ! Set up number of tracers out as an array of integers:
           ! 1: out_top, 2: out_topright, 3: out_right, 4: out_botright
           ! 5: out_bot, 6: out_botleft, 7: out_left, 8: out_topleft
           if (xnew>myxmax.and.ynew>myymax.and.xnew<xu_tot(nx)) then
              ! last statement ensures that tracers leaving the left edge
              ! of the domain (to wrap around to the right side) aren't 
              ! erroniously sent to the right neighbor (should go to left 
              ! neighbor in periodic domain)
              !
              ! Don't need to worry about this for top and bottom because
              ! those boundaries are never periodic
              !
              my_tracers_out(2)=my_tracers_out(2)+1
              tr_out_topright(count_tr)=i
              count_tr=count_tr+1
           elseif (xnew>myxmax.and.ynew>myymax.and.xnew>=xu_tot(nx)) then
              my_tracers_out(8)=my_tracers_out(8)+1
              tr_out_topleft(count_tl)=i
              count_tl=count_tl+1
           elseif (xnew>myxmax.and.ynew<myymin.and.xnew<xu_tot(nx)) then
              my_tracers_out(4)=my_tracers_out(4)+1
              tr_out_botright(count_br)=i
              count_br=count_br+1
           elseif (xnew>myxmax.and.ynew<myymin.and.xnew>=xu_tot(nx)) then
              my_tracers_out(6)=my_tracers_out(6)+1
              tr_out_botleft(count_bl)=i
              count_bl=count_bl+1
           elseif (xnew<myxmin.and.ynew<myymin.and.xnew>xu_tot(2)) then
              my_tracers_out(6)=my_tracers_out(6)+1
              tr_out_botleft(count_bl)=i
              count_bl=count_bl+1
           elseif (xnew<myxmin.and.ynew<myymin.and.xnew<=xu_tot(2)) then
              my_tracers_out(4)=my_tracers_out(4)+1
              tr_out_botright(count_br)=i
              count_br=count_br+1
           elseif (xnew<myxmin.and.ynew>myymax.and.xnew>xu_tot(2)) then
              my_tracers_out(8)=my_tracers_out(8)+1
              tr_out_topleft(count_tl)=i
              count_tl=count_tl+1
           elseif (xnew<myxmin.and.ynew>myymax.and.xnew<=xu_tot(2)) then
              my_tracers_out(2)=my_tracers_out(2)+1
              tr_out_topright(count_tr)=i
              count_tr=count_tr+1
           elseif (xnew>myxmax.and.xnew<xu_tot(nx)) then
              my_tracers_out(3)=my_tracers_out(3)+1
              tr_out_right(count_r)=i
              count_r=count_r+1
           elseif (xnew>myxmax.and.xnew>=xu_tot(nx)) then
              my_tracers_out(7)=my_tracers_out(7)+1
              tr_out_left(count_l)=i
              count_l=count_l+1
           elseif (xnew<myxmin.and.xnew>xu_tot(2)) then
              my_tracers_out(7)=my_tracers_out(7)+1
              tr_out_left(count_l)=i
              count_l=count_l+1
           elseif (xnew<myxmin.and.xnew<=xu_tot(2)) then
              my_tracers_out(3)=my_tracers_out(3)+1
              tr_out_right(count_r)=i
              count_r=count_r+1
           elseif (ynew>myymax) then
              my_tracers_out(1)=my_tracers_out(1)+1
              tr_out_top(count_t)=i
              count_t=count_t+1
           elseif (ynew<myymin) then
              my_tracers_out(5)=my_tracers_out(5)+1
              tr_out_bot(count_b)=i
              count_b=count_b+1
           endif
        enddo
        
        if (rhbcu=='per'.and.lhbcu=='per') then
           deallocate(u_tmp)
           deallocate(xu_tmp)
        endif
        do i=1,myntracers
           deallocate(tracers(i)%d)
        enddo
        deallocate(tracers)
        ! Total number of tracers leaving my domain
        my_tracers_out_tot=sum(my_tracers_out)
       
        ! Now need to know incoming tracers from every side of my domain,
        ! and the total number of incoming tracers; subroutine takes care of 
        ! that, returns my_tracers_in, interger arry of incoming tracers
        ! in same format as my_tracers_out
        call nbrex2d_tracer_counts(my_tracers_out,my_tracers_in,myid,nbrleft,&
             nbrright,nbrtop,nbrbottom,nbrtopleft,nbrtopright,nbrbotleft,&
             nbrbotright,comm2d)
        my_tracers_in_tot=sum(my_tracers_in)
        
        ! New number of tracers
        myntracers_new=myntracers-my_tracers_out_tot+my_tracers_in_tot
        
        ! Reallocate tracers array
        allocate(tracers(myntracers_new))
        ! Copy in the oldtracers that stayed in my domain
        count=1
        do i=1,size(oldtracers,1)
           if (oldtracers(i)%d(2)<=myxmax.and.oldtracers(i)%d(2)>=myxmin&
                .and.oldtracers(i)%d(3)<=myymax.and.&
                oldtracers(i)%d(3)>=myymin) then
              allocate(tracers(count)%d(3))
              tracers(count)%d=oldtracers(i)%d
              count=count+1
           endif
        enddo
        
        ! Last thing to do, send tracers that left my domain into their new 
        ! domain, and receive incoming tracers
        !
        ! Send up, receive from below 
        if (my_tracers_out(1)>0.or.my_tracers_in(5)>0) then
           call nbrex2d_tracer_dat(tracers,oldtracers,my_tracers_out(1),&
                my_tracers_in(5),tr_out_top,nbrtop,nbrbottom,count,comm2d)
           count=count+my_tracers_in(5)
        endif
        ! Send down, receive from above
        if (my_tracers_out(5)>0.or.my_tracers_in(1)>0) then
           call nbrex2d_tracer_dat(tracers,oldtracers,my_tracers_out(5),&
                my_tracers_in(1),tr_out_bot,nbrbottom,nbrtop,count,comm2d)
           count=count+my_tracers_in(1)
        endif
        ! Send right, receive from left
        if (my_tracers_out(3)>0.or.my_tracers_in(7)>0) then
           call nbrex2d_tracer_dat(tracers,oldtracers,my_tracers_out(3),&
                my_tracers_in(7),tr_out_right,nbrright,nbrleft,count,comm2d)
           count=count+my_tracers_in(7)
        endif
        ! Send left, receive from right
        if (my_tracers_out(7)>0.or.my_tracers_in(3)>0) then
           call nbrex2d_tracer_dat(tracers,oldtracers,my_tracers_out(7),&
                   my_tracers_in(3),tr_out_left,nbrleft,nbrright,count,comm2d)
           count=count+my_tracers_in(3)
        endif
        ! Send top-right, receive from bottom-left
        if (my_tracers_out(2)>0.or.my_tracers_in(6)>0) then
           call nbrex2d_tracer_dat(tracers,oldtracers,my_tracers_out(2),&
                my_tracers_in(6),tr_out_topright,nbrtopright,&
                nbrbotleft,count,comm2d)
           count=count+my_tracers_in(6)
        endif
        ! Send bottom-left, receive from top-right
        if (my_tracers_out(6)>0.or.my_tracers_in(2)>0) then
           call nbrex2d_tracer_dat(tracers,oldtracers,my_tracers_out(6),&
                my_tracers_in(2),tr_out_botleft,nbrbotleft,&
                nbrtopright,count,comm2d)
           count=count+my_tracers_in(2)
        endif
        ! Send top-left, receive from bottom-right
        if (my_tracers_out(8)>0.or.my_tracers_in(4)>0) then
           call nbrex2d_tracer_dat(tracers,oldtracers,my_tracers_out(8),&
                my_tracers_in(4),tr_out_topleft,nbrtopleft,&
                nbrbotright,count,comm2d)
           count=count+my_tracers_in(4)
        endif
        ! Send bottom-right, receive from top-left
        if (my_tracers_out(4)>0.or.my_tracers_in(8)>0) then
           call nbrex2d_tracer_dat(tracers,oldtracers,my_tracers_out(4),&
                my_tracers_in(8),tr_out_botright,nbrbotright,&
                nbrtopleft,count,comm2d)
           count=count+my_tracers_in(8)
        endif

        ! Deallocate old tracers
        do i=1,myntracers
           deallocate(oldtracers(i)%d)
        enddo
        deallocate(oldtracers)
        myntracers=myntracers_new

        ! Now covert tracer distribution into composition at cell centers
        ! Need to find the tracers within +/- deltax (deltay) of each 
        ! gridpoint, use the shape function around each grid point to determine
        ! the weighted number of tracers
        call calc_chem_tr_ratio(tracers,chem,xw,y,myntracers,asp,nx,ny,sx,syw,&
             mynx(ng),mynyw(ng),nbrright,nbrleft,nbrtop,nbrbottom,comm2d)
        call nbrex2d(chem,.false.,.false.,.false.,.false.,myid,nbrleft,&
             nbrright,nbrtop,nbrbottom,comm2d)
     endif

!!!!!!!!!!!!!!!!!!!!!Crustal production for grid based scheme!!!!!!!!!!!!!!!!!!
     if (crust_production.and.ntracers==0) then
        ! First determine maximum divergence at surface 
        if (topbdy) then
           div_max=maxval(du_dx(1:mynx(ng)-2,myny(ng)-2))
        endif
        call MPI_ALLREDUCE(div_max,div_max_tot,1,MPI_DOUBLE_PRECISION,&
             MPI_MAX,comm2d,ierr)
        ! Now make crust at given thickness (set by d_crust) in region of 
        ! peak surface divergence (corresponding to spreading centers)
        do i=1,mynx(ng)-2
           if (du_dx(i,myny(ng)-2) >= div_max_tot*0.9_wp) then
              do j=1,myny(ng)-2
                 if (y(ng)%d(j+1) >= d_crust) then
                    chem(i+1,j+1)=1.0_wp
                 endif
              enddo
           endif
        enddo
     endif

!!!!!!!!!!!!!!!!!Damage Equation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (D/=0.0_wp) then
        ! If using the first order strain rate calculation use below, and 
        ! comment out all below
!!$       call strainrate(u,w,dxw(ng)%d,dyu(ng)%d,deltax,deltay,du_dx,du_dy,&
!!$             dw_dx,dw_dy)

        ! Calculate non-dimensional psi (correctly)
        ! Note that psi is not actually the deformational work, because 
        ! viscosity was factored out
        forall (i=1:mynx(ng)-2,j=1:myny(ng)-2)
           deform(i,j)=2.0_wp*(du_dx(i,j)**2.0_wp + dw_dy(i,j)**2.0_wp + &
                du_dy(i,j)*dw_dx(i,j)) + du_dy(i,j)**2.0_wp + &
                dw_dx(i,j)**2.0_wp
        end forall
        
        ! Now formulate the fineness source term
        forall (i=1:mynx(ng)-2,j=1:myny(ng)-2)
           source(i,j)=D*mu_t(i+1,j+1)*deform(i,j)
        end forall     

        ! Grain growth (i.e. fineness sink) term
        do j=1,myny(ng)-2
           do i=1,mynx(ng)-2
              if (temp(i+1,j+1).lt.temp_cutoff) then
                 sink(i,j)=H_dam*exp(-(Eh/(temp_cutoff+etat) - Eh/(t0+etat)))
              else
                 sink(i,j)=H_dam*exp(-(Eh/(temp(i+1,j+1)+etat) - Eh/(t0+etat)))
              endif
           enddo
        enddo

        ! Solve for grain-growth/reduction using operator splitting and 
        ! MPDATA for advection 
        call fineness(alpha,source,sink,dt,m_alpha,p_alpha)
        ! Neighbor exchange to get ghost points
        call nbrex2d(alpha,.true.,.true.,.true.,.true.,myid,nbrleft,&
             nbrright,nbrtop,nbrbottom,comm2d)
        ! Advect fineness using MPDATA
!!$        call mpdata_flux_gen(alpha,u,w,dt,ng,dam_iord,mynx,myny,myid,nbrleft,&
!!$             nbrright,nbrtop,nbrbottom,ctfsl,cbfsl,clfsl,crfsl,comm2d,&
!!$             topbdy,botbdy,'fsl','fsl',utfsl,ubfsl,wlfsl,wrfsl,dxu(ng)%d,&
!!$             dyw(ng)%d,dx(ng)%d,dy(ng)%d,x(ng)%d,y(ng)%d,xu(ng)%d,yw(ng)%d)
        call mpdata_flux(alpha,u,w,dt,ng,dam_iord,mynx,myny,myid,nbrleft,&
             nbrright,nbrtop,nbrbottom,ctfsl,cbfsl,clfsl,crfsl,comm2d,&
             topbdy,botbdy,'fsl','fsl',utfsl,ubfsl,wlfsl,wrfsl)
        call nbrex2d(alpha,.true.,.true.,.true.,.true.,myid,nbrleft,&
             nbrright,nbrtop,nbrbottom,comm2d)        
     endif

!!!!!!!! Advection diffusion equation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Will solve via operator splitting method, where I do a pure advection 
     ! step, then use the output of the advection step as input for a pure 
     ! diffusion step.  Advection step uses MPDATA from Smolarkiewicz; 
     ! diffusion uses center difference with crank-nicholson time discretizatio
     ! Solve for diffusion

     ! First do temperature equation 
     ! Lewis number is set to 1 because it can only be 1 for temp equation
!!$     call mpdata_flux_gen(temp,u,w,dt,ng,iord,mynx,myny,myid,nbrleft,&
!!$          nbrright,nbrtop,nbrbottom,ttfsl,tbfsl,tlfsl,trfsl,comm2d,topbdy,&
!!$          botbdy,topbct,botbct,utfsl,ubfsl,wlfsl,wrfsl,dxu(ng)%d,&
!!$             dyw(ng)%d,dx(ng)%d,dy(ng)%d,x(ng)%d,y(ng)%d,xu(ng)%d,yw(ng)%d)
!!$     call mpdata_flux(temp,u,w,dt,ng,iord,mynx,myny,myid,nbrleft,&
!!$          nbrright,nbrtop,nbrbottom,ttfsl,tbfsl,tlfsl,trfsl,comm2d,topbdy,&
!!$          botbdy,topbct,botbct,utfsl,ubfsl,wlfsl,wrfsl)
!!$     call diffusion(temp,1.0_wp,dt,Q,ng,nglow,dx,dy,mynx,myny,myid,nbrleft,&
!!$          nbrright,nbrtop,nbrbottom,ttfsl,tbfsl,tlfsl,trfsl,comm2d,&
!!$          x(ng)%d,y(ng)%d,dxu,dyw)
     
! 8/27/14: Now code will only enter this routine if ntracers is 0; i.e. if we 
! are using a grid based scheme for composition.  Should also look into how 
! to do diffusion with tracers     
     if (Bu.ne.0.0_wp.and.ntracers==0) then
        ! Now do chemistry and use the specified lewis number
!!$        call mpdata_flux_gen(chem,u,w,dt,ng,iord,mynx,myny,myid,nbrleft,nbrright,&
!!$             nbrtop,nbrbottom,ctfsl,cbfsl,clfsl,crfsl,comm2d,topbdy,botbdy,&
!!$             topbcc,botbcc,utfsl,ubfsl,wlfsl,wrfsl,dxu(ng)%d,&
!!$             dyw(ng)%d,dx(ng)%d,dy(ng)%d,x(ng)%d,y(ng)%d,xu(ng)%d,yw(ng)%d)
        call mpdata_flux(chem,u,w,dt,ng,iord,mynx,myny,myid,nbrleft,nbrright,&
             nbrtop,nbrbottom,ctfsl,cbfsl,clfsl,crfsl,comm2d,topbdy,botbdy,&
             topbcc,botbcc,utfsl,ubfsl,wlfsl,wrfsl)
        call diffusion(chem,lewis,dt,0.0_wp,ng,nglow,dxc,dyc,mynx,myny,myid,&
             nbrleft,nbrright,nbrtop,nbrbottom,ctfsl,cbfsl,clfsl,crfsl,comm2d,&
             x(ng)%d,y(ng)%d,dxu,dyw)
     endif

!!$! Power-law upwind scheme for advection, fully implicit time integration for
!!$! both diffusion and advection, following chapter 5 of Patankar
!!$call advecdiff(temp,lewis,dt,u,w,ng,nglow,x,y,dx,dy,mynx,myny,myid,&
!!$     nbrleft,nbrright,nbrtop,nbrbottom,ttfsl,tbfsl,comm2d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     ! Update time
     time=time+dt

  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF TIME LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  close(21) 
!!$  close(23)

  call MPI_FINALIZE(ierr)
!
  CONTAINS
!
    SUBROUTINE chem_flux_avg(chem,temp,chem_flux_avg_out,dx,dy,mynx,myny,ng)
      USE btype
      USE bfunc2d
      USE mpi
      IMPLICIT NONE
      REAL(WP), DIMENSION(:,:), INTENT(IN) :: chem,temp
      REAL(WP), DIMENSION(:), INTENT(IN) :: dx,dy
      REAL(WP), INTENT(OUT) :: chem_flux_avg_out
      INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
      INTEGER(I4B), INTENT(IN) :: ng
      REAL(WP) :: sum1,sum2,sum1_tot,sum2_tot,c_max,chem_avg
      INTEGER(I4B) :: i,j
      c_max=0.9
      sum1=0.0_wp
      sum2=0.0_wp
      do i=2,mynx(ng)-1
         do j=2,myny(ng)-1
            chem_avg=0.5_wp*(chem(i,j-1)+chem(i,j))
            if (chem_avg>c_max.and.temp(i,j+1)<0.5_wp) then
               sum1=sum1+0.5_wp*(chem(i,j)*w(i,j)+chem(i,j-1)*w(i,j-1))*&
                    dx(i-1)*dy(j-1)
               sum2=sum2+dx(i-1)*dy(j-1)
            endif
         enddo
      enddo
      call MPI_REDUCE(sum1,sum1_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
           comm2d,ierr)
      call MPI_REDUCE(sum2,sum2_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
           comm2d,ierr)
      if (myid==0) then
         if (sum1_tot==0.0_wp) then
            chem_flux_avg_out=0.0_wp
         else
            chem_flux_avg_out=sum1_tot/sum2_tot
         endif
      endif
    END SUBROUTINE chem_flux_avg
!
    SUBROUTINE crust_depth(chem,crust_depth_avg,crust_depth_max,mynx,&
         myny,ng)
      USE btype
      USE bfunc2d
      USE mpi
      IMPLICIT NONE
      REAL(WP), DIMENSION(:,:), INTENT(IN) :: chem
      INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
      INTEGER(I4B), INTENT(IN) :: ng
      REAL(WP), INTENT(OUT) :: crust_depth_avg,crust_depth_max
      REAL(WP), DIMENSION(:), ALLOCATABLE :: crust_z,crust_z_tot
      INTEGER(I4B) :: i,j
      REAL(WP) :: ne,nw,nn,ns,sum,depth_chem_contour,c_max
      !
      c_max=0.9
      allocate(crust_z(mynx(ng)-2))
      crust_z=0.0_wp
      do i=2,mynx(ng)-1
         do j=2,myny(ng)-1
            if (chem(i,j).ge.c_max) then
               ne=chem(i+1,j)
               nw=chem(i-1,j)
               nn=chem(i,j+1)
               ns=chem(i,j-1)
               if (ne<c_max.or.nw<c_max.or.ns<c_max.or.nn<c_max) then
                  if (i-1==1) then
                     crust_z(i-1)=yw(ng)%d(j)
                     exit
                  elseif (yw(ng)%d(j)>crust_z(i-2)-0.1_wp) then
                     crust_z(i-1)=yw(ng)%d(j)
                     exit
                  endif
               endif
            elseif (chem(i,j)<c_max.and.j==myny(ng)-1) then
               crust_z(i-1)=1.0_wp
            endif
         enddo
      enddo
      ! Now need to use mpi_gather to get deepest occurence of c=c_max 
      ! contour at each x position
      ! 
      ! Putting max depth of c=c_max contour for all x on root processor
      if (myid==0) then
         allocate(crust_z_tot(nx))
      endif
      sum=0.0_wp
      do i=1,nx
         if (i>=sx.and.i<=ex) then
            depth_chem_contour=crust_z(i-sx+1)
         else
            depth_chem_contour=1.0_wp
         endif
         call MPI_REDUCE(depth_chem_contour,crust_z_tot(i),1,&
              MPI_DOUBLE_PRECISION,MPI_MIN,0,comm2d,ierr)
         if (myid==0) then
            sum=sum+(1.0_wp-crust_z_tot(i))
         endif
      enddo
      ! Now average crust depth, and pick out max crust depth
      if (myid==0) then
         crust_depth_avg=sum/nx
         crust_depth_max=minval(crust_z_tot)
         deallocate(crust_z_tot)
      endif
      deallocate(crust_z)
    END SUBROUTINE crust_depth
!
    SUBROUTINE nusselt(temp,dy,dx,asp,nuss,mynx,ng,tbdy) 
      USE btype
      USE bfunc2d
      USE mpi
      IMPLICIT NONE 
      REAL(WP), DIMENSION(:), INTENT(IN) :: temp,dx
      REAL(WP), INTENT(IN) :: dy,asp
      REAL(WP), INTENT(OUT) :: nuss
      INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx
      INTEGER(I4B), INTENT(IN) :: ng
      LOGICAL, INTENT(IN) :: tbdy
      INTEGER(I4B) :: i
      REAL(WP) :: sum,sum_tot
      !
      sum=0.0_wp
      if (tbdy) then
         do i=2,mynx(ng)-1
            sum=sum+((temp(i)-tbvt)/dy)*dx(i-1)
         enddo
      endif
      call MPI_REDUCE(sum,sum_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
           comm2d,ierr)
      if (myid==0) then
         nuss=sum_tot/asp
      endif
    END SUBROUTINE nusselt
!
    SUBROUTINE usurf(u,dx,asp,uavg,mynx,ng,tbdy) 
      USE btype
      USE bfunc2d
      USE mpi
      IMPLICIT NONE 
      REAL(WP), DIMENSION(:), INTENT(IN) :: u,dx
      REAL(WP), INTENT(IN) :: asp
      REAL(WP), INTENT(OUT) :: uavg
      INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx
      INTEGER(I4B), INTENT(IN) :: ng
      LOGICAL, INTENT(IN) :: tbdy
      INTEGER(I4B) :: i
      REAL(WP) :: sum,sum_tot
      !
      sum=0.0_wp
      if (tbdy) then
         do i=1,mynx(ng)-2
            sum=sum+0.5_wp*(u(i)**2.0_wp+u(i+1)**2.0_wp)*dx(i)
         enddo
      endif
      call MPI_REDUCE(sum,sum_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
           comm2d,ierr)
      if (myid==0) then
         uavg=sqrt(sum_tot/asp)
      endif
    END SUBROUTINE usurf
!
    SUBROUTINE rms(u,w,dy,dx,asp,v_rms,mynx,myny,ng) 
      USE btype
      USE bfunc2d
      USE mpi
      IMPLICIT NONE 
      REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,w
      REAL(WP), DIMENSION(:), INTENT(IN) :: dy,dx
      REAL(WP), INTENT(IN) :: asp
      REAL(WP), INTENT(OUT) :: v_rms
      INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
      INTEGER(I4B), INTENT(IN) :: ng
      INTEGER(I4B) :: i,j
      REAL(WP) :: sum,sum_tot
      !
      sum=0.0_wp
      do i=1,mynx(ng)-2
         do j=1,myny(ng)-2
            sum=sum+0.5_wp*(u(i,j+1)**2.0_wp+u(i+1,j+1)**2.0_wp)*dx(i)*dy(j)+&
                0.5_wp*(w(i+1,j)**2.0_wp+w(i+1,j+1)**2.0_wp)*dx(i)*dy(j) 
         enddo
      enddo
      call MPI_REDUCE(sum,sum_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm2d,ierr)
      if (myid==0) then
         v_rms=sqrt(sum_tot/asp)
      endif
    END SUBROUTINE rms
!
    SUBROUTINE average_u(phi,dy,dx,asp,avg,mynx,myny,ng) 
      USE btype
      USE bfunc2d
      USE mpi
      IMPLICIT NONE 
      REAL(WP), DIMENSION(:,:), INTENT(IN) :: phi
      REAL(WP), DIMENSION(:), INTENT(IN) :: dy,dx
      REAL(WP), INTENT(IN) :: asp
      REAL(WP), INTENT(OUT) :: avg
      INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
      INTEGER(I4B), INTENT(IN) :: ng
      INTEGER(I4B) :: i,j
      REAL(WP) :: sum,sum_tot
      !
      sum=0.0_wp
      do i=1,mynx(ng)-2
         do j=1,myny(ng)-2
            sum=sum+0.5_wp*(phi(i,j+1)+phi(i+1,j+1))*dx(i)*dy(j) 
         enddo
      enddo
      call MPI_ALLREDUCE(sum,sum_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
           comm2d,ierr)
      avg=sum_tot/asp
    END SUBROUTINE average_u
!    
    SUBROUTINE average(phi,power,dy,dx,asp,avg,mynx,myny,ng) 
      USE btype
      USE bfunc2d
      USE mpi
      IMPLICIT NONE 
      REAL(WP), DIMENSION(:,:), INTENT(IN) :: phi
      REAL(WP), DIMENSION(:), INTENT(IN) :: dy,dx
      REAL(WP), INTENT(IN) :: asp,power
      REAL(WP), INTENT(OUT) :: avg
      INTEGER(I4B), DIMENSION(:), INTENT(IN) :: mynx,myny
      INTEGER(I4B), INTENT(IN) :: ng
      INTEGER(I4B) :: i,j
      REAL(WP) :: sum,sum_tot
      !
      sum=0.0_wp
      do i=2,mynx(ng)-1
         do j=2,myny(ng)-1
            sum=sum+phi(i,j)**power*dx(i-1)*dy(j-1) 
         enddo
      enddo
      call MPI_ALLREDUCE(sum,sum_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
           comm2d,ierr)
      avg=(sum_tot/asp)**(1.0_wp/power)
    END SUBROUTINE average
!   
    SUBROUTINE boundary(phi,lhbdy,rhbdy,topbdy,botbdy,lhbc,rhbc,&
         topbc,botbc,lbv,rbv,tbv,bbv)
      USE btype
      USE bfunc2d
      IMPLICIT NONE
      REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: phi
      LOGICAL, INTENT(IN) :: lhbdy,rhbdy,topbdy,botbdy
      CHARACTER(len=3), INTENT(IN) :: lhbc,rhbc,topbc,botbc
      REAL(WP), INTENT(IN) :: lbv,rbv,tbv,bbv
      INTEGER(I4B) :: nx,ny
      !
      nx=size(phi,1)
      ny=size(phi,2)
      if (lhbdy) then 
         if (lhbc=='con') then
            phi(1,:)=lbv
         endif
      endif
      if (rhbdy) then 
         if (rhbc=='con') then
            phi(nx,:)=rbv
         endif
      endif
      if (topbdy) then 
         if (topbc=='con') then
            phi(:,ny)=tbv
         endif
      endif
      if (botbdy) then 
         if (botbc=='con') then
            phi(:,1)=bbv
         endif
      endif
    END SUBROUTINE boundary
!
    SUBROUTINE combine_root(u,w,max,nx1,mynx,myny,numprocs,sxall,syall,&
         mynxall,mynyall,comm)
      USE btype
      USE mpi 
      IMPLICIT NONE 
      REAL(WP), DIMENSION(:,:), INTENT(IN) :: u
      REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: w
      INTEGER, INTENT(IN) :: nx1,max,mynx,myny,numprocs,comm
      INTEGER, DIMENSION(:), INTENT(IN) :: sxall,syall,mynxall,mynyall
      INTEGER :: i,j,recv(numprocs),displ(numprocs),send
      do i=1,max
         do j=1,numprocs
            displ(j)=nx1*(syall(j)-1) + (sxall(j)-1) + (i-1)*nx1
            if (i>mynyall(j)-2) then
               recv(j)=0
            else
               recv(j)=mynxall(j)-2
            endif
         enddo
         if (i>myny-2) then
            send=0
         else
            send=mynx-2
         endif
         call MPI_GATHERV(u(2:mynx-1,i+1),send,MPI_DOUBLE_PRECISION,&
              w,recv,displ,MPI_DOUBLE_PRECISION,0,comm,ierr)
      enddo
    END SUBROUTINE combine_root
!
    SUBROUTINE combine_rest(u,max,mynx,myny,comm)
      USE btype
      USE mpi 
      IMPLICIT NONE 
      REAL(WP), DIMENSION(:,:), INTENT(IN) :: u
      INTEGER, INTENT(IN) :: max,mynx,myny,comm
      INTEGER :: i,send
      do i=1,max
         if (i>myny-2) then
            send=0
         else
            send=mynx-2
         endif
         call MPI_GATHERV(u(2:mynx-1,i+1),send,MPI_DOUBLE_PRECISION,0,0,0,&
              MPI_DOUBLE_PRECISION,0,comm,ierr)
      enddo
    END SUBROUTINE combine_rest
!
    SUBROUTINE combine_root_bdy(u,w,mynx,myny,numprocs,sxall,syall,&
         eyall,mynxall,topbdy,botbdy,comm)
      USE btype
      USE mpi 
      IMPLICIT NONE 
      REAL(WP), DIMENSION(:,:), INTENT(IN) :: u
      REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: w
      INTEGER, INTENT(IN) :: mynx,myny,numprocs,comm
      INTEGER, DIMENSION(:), INTENT(IN) :: sxall,syall,mynxall,eyall
      LOGICAL, INTENT(IN) :: topbdy,botbdy
      REAL(WP), DIMENSION(nx) :: top,bottom
      INTEGER :: i,j,recv_top(numprocs),recv_bot(numprocs),&
           displ_top(numprocs),displ_bot(numprocs),send
      do j=1,numprocs
         if (eyall(j)==ny) then
            displ_top(j)=sxall(j)-1
            recv_top(j)=mynxall(j)-2
         else
            displ_top(j)=0
            recv_top(j)=0
         endif
         if (syall(j)-1==1) then
            displ_bot(j)=sxall(j)-1
            recv_bot(j)=mynxall(j)-2
         else
            displ_bot(j)=0
            recv_bot(j)=0
         endif
      enddo
      if (topbdy) then
         send=mynx-2
      else
         send=0
      endif
      call MPI_GATHERV(u(2:mynx-1,myny),send,MPI_DOUBLE_PRECISION,&
           top,recv_top,displ_top,MPI_DOUBLE_PRECISION,0,comm,ierr)
      if (botbdy) then
         send=mynx-2
      else
         send=0
      endif
      call MPI_GATHERV(u(2:mynx-1,1),send,MPI_DOUBLE_PRECISION,&
           bottom,recv_bot,displ_bot,MPI_DOUBLE_PRECISION,0,comm,ierr)
      w(1:nx,1)=bottom
      w(1:nx,ny+1)=top
    END SUBROUTINE combine_root_bdy
!
    SUBROUTINE combine_rest_bdy(u,topbdy,botbdy,mynx,myny,comm)
      USE btype
      USE mpi 
      IMPLICIT NONE 
      REAL(WP), DIMENSION(:,:), INTENT(IN) :: u
      INTEGER, INTENT(IN) :: mynx,myny,comm
      LOGICAL, INTENT(IN) :: topbdy,botbdy
      INTEGER :: i,send
      if (topbdy) then
         send=mynx-2
      else
         send=0
      endif
      ! May not work because I've set the send counts to 0, but still have
      ! u(2:mynx-1,myny) as the send buffer
      call MPI_GATHERV(u(2:mynx-1,myny),send,MPI_DOUBLE_PRECISION,0,0,0,&
           MPI_DOUBLE_PRECISION,0,comm,ierr)
      if (botbdy) then
         send=mynx-2
      else
         send=0
      endif
      call MPI_GATHERV(u(2:mynx-1,1),send,MPI_DOUBLE_PRECISION,0,0,0,&
           MPI_DOUBLE_PRECISION,0,comm,ierr)
    END SUBROUTINE combine_rest_bdy
!
    SUBROUTINE ha_root(phi,ha,dx,mynx,myny,sxall,syall,mynyall,npx,asp,&
         numprocs,maxrow,comm)
      USE btype
      USE mpi
      IMPLICIT NONE
      REAL(WP), DIMENSION(:,:), INTENT(IN) :: phi
      REAL(WP), DIMENSION(:), INTENT(INOUT) :: ha
      REAL(WP), DIMENSION(:), INTENT(IN) :: dx
      INTEGER, DIMENSION(:), INTENT(IN) :: sxall,syall,mynyall
      INTEGER, INTENT(IN) :: mynx,myny,npx,numprocs,maxrow,comm
      REAL(WP), INTENT(IN) :: asp
      REAL(WP) :: sum,sum2
      REAL(WP), DIMENSION(npx,ny) :: hap
      INTEGER :: i,j,recv(numprocs),displ(numprocs),send
      !
      ! Integrate over each row 
      do j=1,maxrow
         sum=0.0_wp
         do i=2,mynx-1
            sum = sum + phi(i,j+1)*dx(i-1)
         enddo
         do i=1,numprocs
            displ(i)=npx*(syall(i)-1) + nint((sxall(i)-1)/real(nx/npx)) &
                 + (j-1)*npx
            if (j>mynyall(i)-2) then
               recv(i)=0
            else
               recv(i)=1
            endif
         enddo
         if (j>myny-2) then
            send=0
         else
            send=1
         endif
         call MPI_GATHERV(sum,send,MPI_DOUBLE_PRECISION,hap,recv,displ,&
              MPI_DOUBLE_PRECISION,0,comm,ierr)
      enddo
      do j=1,ny
         sum2=0.0_wp
         do i=1,npx
            sum2 = sum2 + hap(i,j)
         enddo
         ha(j)=sum2/asp
      enddo
    END SUBROUTINE ha_root
!
    SUBROUTINE ha_rest(phi,dx,mynx,myny,maxrow,comm)
      USE btype
      USE mpi
      IMPLICIT NONE
      REAL(WP), DIMENSION(:,:), INTENT(IN) :: phi
      INTEGER, INTENT(IN) :: mynx,myny,maxrow,comm
      REAL(WP), DIMENSION(:), INTENT(IN) :: dx
      REAL(WP) :: sum
      INTEGER :: i,send
      !
      ! Integrate over each row 
      do j=1,maxrow
         sum=0.0_wp
         do i=2,mynx-1
            sum = sum + phi(i,j+1)*dx(i-1)
         enddo
         if (j>myny-2) then
            send=0
         else
            send=1
         endif
         call MPI_GATHERV(sum,send,MPI_DOUBLE_PRECISION,0,0,0,&
              MPI_DOUBLE_PRECISION,0,comm,ierr)
      enddo
    END SUBROUTINE ha_rest
!
    SUBROUTINE ha_root_chem(phi,ha,dx,mynx,myny,sxall,syall,mynyall,&
         npx,asp,numprocs,comm)
      USE btype
      USE mpi
      IMPLICIT NONE
      REAL(WP), DIMENSION(:,:), INTENT(IN) :: phi
      REAL(WP), DIMENSION(:), INTENT(INOUT) :: ha
      REAL(WP), DIMENSION(:), INTENT(IN) :: dx
      INTEGER, DIMENSION(:), INTENT(IN) :: sxall,syall,mynyall
      INTEGER, INTENT(IN) :: mynx,myny,npx,numprocs,comm
      REAL(WP), INTENT(IN) :: asp
      REAL(WP) :: sum,sum2
      REAL(WP), DIMENSION(npx,ny) :: hap
      INTEGER :: i,j,recv(numprocs),displ(numprocs),send
      !
      ! Integrate over each row 
      ! Myny should be for the grid defined at cell centers
      do j=1,myny-2
         sum=0.0_wp
         do i=2,mynx-1
            sum = sum + 0.5_wp*(phi(i,j+1)+phi(i,j))*dx(i-1)
         enddo
         do i=1,numprocs
            displ(i)=npx*(syall(i)-1) + nint((sxall(i)-1)/real(nx/npx)) &
                 + (j-1)*npx 
            recv(i)=1
         enddo
         send=1
         call MPI_GATHERV(sum,send,MPI_DOUBLE_PRECISION,hap,recv,displ,&
              MPI_DOUBLE_PRECISION,0,comm,ierr)
      enddo
      do j=1,ny
         sum2=0.0_wp
         do i=1,npx
            sum2 = sum2 + hap(i,j)
         enddo
         ha(j)=sum2/asp
      enddo
    END SUBROUTINE ha_root_chem
!
    SUBROUTINE ha_rest_chem(phi,dx,mynx,myny,comm) 
      !myny should be for cell centers
      USE btype
      USE mpi
      IMPLICIT NONE
      REAL(WP), DIMENSION(:,:), INTENT(IN) :: phi
      INTEGER, INTENT(IN) :: mynx,myny,comm
      REAL(WP), DIMENSION(:), INTENT(IN) :: dx
      REAL(WP) :: sum
      INTEGER :: i,send
      !
      ! Integrate over each row 
      do j=2,myny-1
         sum=0.0_wp
         do i=2,mynx-1
            sum = sum + 0.5_wp*(phi(i,j)+phi(i,j-1))*dx(i-1)
         enddo
         send=1
         call MPI_GATHERV(sum,send,MPI_DOUBLE_PRECISION,0,0,0,&
              MPI_DOUBLE_PRECISION,0,comm,ierr)
      enddo
    END SUBROUTINE ha_rest_chem
!
    SUBROUTINE gather_tracers_root(tracers,tracers_tot,send_count,&
         ntracers_tot,comm2d)
      USE btype
      USE mpi
      IMPLICIT NONE 
      TYPE(ptr1d), INTENT(IN) :: tracers(:) 
      REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: tracers_tot
      INTEGER(I4B), INTENT(IN) :: send_count,ntracers_tot,comm2d
      REAL(WP), DIMENSION(send_count*3) :: send_dat
      REAL(WP), DIMENSION(ntracers*3) :: recv_dat
      INTEGER(I4B) :: ierr,i,j,k,myntracers_tot(numprocs),stride(numprocs),&
           recv_count(numprocs)
      ! Put tracer data into a dummy array so it can be sent to root
      k=1
      do i=1,send_count
         do j=1,3
            send_dat(k)=tracers(i)%d(j)
            k=k+1
         enddo
      enddo
      ! All processes send their send_count (i.e. number of tracers in their 
      ! sub-domain to the root
      call MPI_GATHER(send_count,1,MPI_INTEGER,myntracers_tot,1,MPI_INTEGER,&
           0,comm2d,ierr)
      recv_count=myntracers_tot*3
      ! Now determine the stride for each processor's incoming tracer data
      stride(1)=0
      do i=2,numprocs
         stride(i)=stride(i-1)+myntracers_tot(i-1)*3
      enddo
      ! gather to root process
      call MPI_GATHERV(send_dat,send_count*3,MPI_DOUBLE_PRECISION,recv_dat,&
           recv_count,stride,MPI_DOUBLE_PRECISION,0,comm2d,ierr)
      k=1
      do i=1,ntracers
         do j=1,3
            tracers_tot(i,j)=recv_dat(k)
            k=k+1
         enddo
      enddo
    END SUBROUTINE gather_tracers_root
!
    SUBROUTINE gather_tracers_rest(tracers,send_count,comm2d)
      USE btype
      USE mpi
      IMPLICIT NONE 
      TYPE(ptr1d), INTENT(IN) :: tracers(:) 
      INTEGER(I4B), INTENT(IN) :: send_count,comm2d
      REAL(WP), DIMENSION(send_count*3) :: send_dat
      INTEGER(I4B) :: ierr,i,j,k
      ! Put tracer data into a dummy array so it can be sent to root
      call MPI_GATHER(send_count,1,MPI_INTEGER,0,0,MPI_INTEGER,&
           0,comm2d,ierr)
      k=1
      do i=1,send_count
         do j=1,3
            send_dat(k)=tracers(i)%d(j)
            k=k+1
         enddo
      enddo
      ! gather to root process
      call MPI_GATHERV(send_dat,send_count*3,MPI_DOUBLE_PRECISION,0,0,0,&
           MPI_DOUBLE_PRECISION,0,comm2d,ierr)
    END SUBROUTINE gather_tracers_rest
! 
    SUBROUTINE calc_chem_tr_ratio(tracers,chem,x,y,myntracers,asp,nx,ny,sx,sy,&
         mynx1,myny1,nbrright,nbrleft,nbrtop,nbrbottom,comm2d)
      USE btype
      USE mpi
      IMPLICIT NONE
      TYPE(ptr1d), INTENT(IN) :: tracers(:),x(:),y(:)
      REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: chem
      REAL(WP), INTENT(IN) :: asp
      INTEGER(I4B), INTENT(IN) :: myntracers,nx,ny,sx,sy,nbrright,nbrleft,&
           nbrtop,nbrbottom,comm2d,mynx1,myny1
      REAL(WP), DIMENSION(mynx1,myny1) :: n_tr_weight_dense,&
           n_tr_weight_reg
      REAL(WP), DIMENSION(mynx1) :: partial_dense_bot,partial_dense_top,&
           partial_reg_bot,partial_reg_top
      REAL(WP), DIMENSION(myny1) :: partial_dense_right,&
           partial_dense_left,partial_reg_right,partial_reg_left
      REAL(WP) :: partial_dense_tr,partial_dense_tl,partial_dense_bl,&
      partial_dense_br,partial_reg_tr,partial_reg_tl,partial_reg_bl,&
      partial_reg_br
      INTEGER :: i,j,n,m,nx1,nx2,ny1,ny2 !,fix_i(100),fix_j(100) 
      !
      ! 8/29/14: Making this routine more efficient by looping through the 
      ! tracers array only once; make n_tr_weight_ arrays, same size as chem,
      ! and calculate the n_tr_weight_(i,j) on the fly as I loop through all 
      ! the tracers in my domain.  In other words, loop through my tracers 
      ! and calculate the contribution each tracer makes to the n_tr_weight_
      ! at the grid points nearest this tracer
      !
      ! Will need to sort out grid points on the boundary between processors, 
      ! because I don't have "ghost tracers." I'll need to calculate a partial 
      ! sum on each grid point and exchange that between neighbors
      n_tr_weight_dense=0.0_wp
      n_tr_weight_reg=0.0_wp
      partial_dense_right=0.0_wp; partial_dense_left=0.0_wp 
      partial_dense_top=0.0_wp; partial_dense_bot=0.0_wp
      partial_reg_right=0.0_wp; partial_reg_left=0.0_wp 
      partial_reg_top=0.0_wp; partial_reg_bot=0.0_wp 
      partial_dense_tr=0.0_wp; partial_dense_tl=0.0_wp; partial_dense_bl=0.0_wp
      partial_dense_br=0.0_wp
      partial_reg_tr=0.0_wp; partial_reg_tl=0.0_wp; partial_reg_bl=0.0_wp
      partial_reg_br=0.0_wp
      do n=1,myntracers
         ! First figure out the four nearest grid points for chem by 
         ! determining i of the grid points to the left and right, and j of 
         ! the grid points above and below 
         nx1=floor((tracers(n)%d(2)+(1.0_wp*asp/(2.0_wp*nx)))*(nx/asp)+2.0_wp-sx)
         nx2=ceiling((tracers(n)%d(2)+(1.0_wp*asp/(2.0_wp*nx)))*(nx/asp)+2.0_wp-sx)
         ny1=floor(tracers(n)%d(3)*ny+3.0_wp-sy)
         ny2=ceiling(tracers(n)%d(3)*ny+3.0_wp-sy)

         ! Do all grid points; ghost points (e.g. points along a 
         ! processor's boundary will end up calculating a partial sum of 
         ! tracers shared with neighboring processors.  This way I can just 
         ! exchange these partial sums between processors and get the final, 
         ! weighted sums of dense and regular tracers
         if (tracers(n)%d(1)==1.0_wp.and.nx1/=nx2.and.ny1/=ny2) then
            n_tr_weight_dense(nx1,ny1)=n_tr_weight_dense(nx1,ny1)+&
                 max(0.0_wp,1.0_wp-(abs(tracers(n)%d(2)-x(ng)%d(nx1)))/&
                 (x(ng)%d(nx2)-x(ng)%d(nx1)))*max(0.0_wp,1.0_wp-&
                 (abs(tracers(n)%d(3)-y(ng)%d(ny1)))/(y(ng)%d(ny2)-&
                 y(ng)%d(ny1)))
            n_tr_weight_dense(nx2,ny1)=n_tr_weight_dense(nx2,ny1)+&
                 max(0.0_wp,1.0_wp-(abs(tracers(n)%d(2)-x(ng)%d(nx2)))/&
                 (x(ng)%d(nx2)-x(ng)%d(nx1)))*max(0.0_wp,1.0_wp-&
                 (abs(tracers(n)%d(3)-y(ng)%d(ny1)))/(y(ng)%d(ny2)-&
                 y(ng)%d(ny1)))
            n_tr_weight_dense(nx1,ny2)=n_tr_weight_dense(nx1,ny2)+&
                 max(0.0_wp,1.0_wp-(abs(tracers(n)%d(2)-x(ng)%d(nx1)))/&
                 (x(ng)%d(nx2)-x(ng)%d(nx1)))*max(0.0_wp,1.0_wp-&
                 (abs(tracers(n)%d(3)-y(ng)%d(ny2)))/(y(ng)%d(ny2)-&
                 y(ng)%d(ny1)))
            n_tr_weight_dense(nx2,ny2)=n_tr_weight_dense(nx2,ny2)+&
                 max(0.0_wp,1.0_wp-(abs(tracers(n)%d(2)-x(ng)%d(nx2)))/&
                 (x(ng)%d(nx2)-x(ng)%d(nx1)))*max(0.0_wp,1.0_wp-&
                 (abs(tracers(n)%d(3)-y(ng)%d(ny2)))/(y(ng)%d(ny2)-&
                 y(ng)%d(ny1)))
         elseif (tracers(n)%d(1)==1.0_wp.and.nx1/=nx2.and.ny1==ny2) then
            n_tr_weight_dense(nx1,ny1)=n_tr_weight_dense(nx1,ny1)+&
                 max(0.0_wp,1.0_wp-(abs(tracers(n)%d(2)-x(ng)%d(nx1)))/&
                 (x(ng)%d(nx2)-x(ng)%d(nx1)))
            n_tr_weight_dense(nx2,ny1)=n_tr_weight_dense(nx2,ny1)+&
                 max(0.0_wp,1.0_wp-(abs(tracers(n)%d(2)-x(ng)%d(nx2)))/&
                 (x(ng)%d(nx2)-x(ng)%d(nx1)))
         elseif (tracers(n)%d(1)==1.0_wp.and.ny1/=ny2.and.nx1==nx2) then
            n_tr_weight_dense(nx1,ny1)=n_tr_weight_dense(nx1,ny1)+&
                 max(0.0_wp,1.0_wp-(abs(tracers(n)%d(3)-y(ng)%d(ny1)))&
                 /(y(ng)%d(ny2)-y(ng)%d(ny1)))
            n_tr_weight_dense(nx1,ny2)=n_tr_weight_dense(nx1,ny2)+&
                 max(0.0_wp,1.0_wp-(abs(tracers(n)%d(3)-y(ng)%d(ny2)))&
                 /(y(ng)%d(ny2)-y(ng)%d(ny1)))
         elseif (tracers(n)%d(1)==1.0_wp.and.ny1==ny2.and.nx1==nx2) then
            n_tr_weight_dense(nx1,ny1)=n_tr_weight_dense(nx1,ny1)+1.0_wp
         elseif (tracers(n)%d(1)==0.0_wp.and.nx1/=nx2.and.ny1/=ny2) then
            n_tr_weight_reg(nx1,ny1)=n_tr_weight_reg(nx1,ny1)+&
                 max(0.0_wp,1.0_wp-(abs(tracers(n)%d(2)-x(ng)%d(nx1)))/&
                 (x(ng)%d(nx2)-x(ng)%d(nx1)))*max(0.0_wp,1.0_wp-&
                 (abs(tracers(n)%d(3)-y(ng)%d(ny1)))/(y(ng)%d(ny2)-&
                 y(ng)%d(ny1)))
            n_tr_weight_reg(nx2,ny1)=n_tr_weight_reg(nx2,ny1)+&
                 max(0.0_wp,1.0_wp-(abs(tracers(n)%d(2)-x(ng)%d(nx2)))/&
                 (x(ng)%d(nx2)-x(ng)%d(nx1)))*max(0.0_wp,1.0_wp-&
                 (abs(tracers(n)%d(3)-y(ng)%d(ny1)))/(y(ng)%d(ny2)-&
                 y(ng)%d(ny1)))
            n_tr_weight_reg(nx1,ny2)=n_tr_weight_reg(nx1,ny2)+&
                 max(0.0_wp,1.0_wp-(abs(tracers(n)%d(2)-x(ng)%d(nx1)))/&
                 (x(ng)%d(nx2)-x(ng)%d(nx1)))*max(0.0_wp,1.0_wp-&
                 (abs(tracers(n)%d(3)-y(ng)%d(ny2)))/(y(ng)%d(ny2)-&
                 y(ng)%d(ny1)))
            n_tr_weight_reg(nx2,ny2)=n_tr_weight_reg(nx2,ny2)+&
                 max(0.0_wp,1.0_wp-(abs(tracers(n)%d(2)-x(ng)%d(nx2)))/&
                 (x(ng)%d(nx2)-x(ng)%d(nx1)))*max(0.0_wp,1.0_wp-&
                 (abs(tracers(n)%d(3)-y(ng)%d(ny2)))/(y(ng)%d(ny2)-&
                 y(ng)%d(ny1)))
         elseif (tracers(n)%d(1)==0.0_wp.and.nx1/=nx2.and.ny1==ny2) then
            n_tr_weight_reg(nx1,ny1)=n_tr_weight_reg(nx1,ny1)+&
                 max(0.0_wp,1.0_wp-(abs(tracers(n)%d(2)-x(ng)%d(nx1)))/&
                 (x(ng)%d(nx2)-x(ng)%d(nx1)))
            n_tr_weight_reg(nx2,ny1)=n_tr_weight_reg(nx2,ny1)+&
                 max(0.0_wp,1.0_wp-(abs(tracers(n)%d(2)-x(ng)%d(nx2)))/&
                 (x(ng)%d(nx2)-x(ng)%d(nx1)))
         elseif (tracers(n)%d(1)==0.0_wp.and.ny1/=ny2.and.nx1==nx2) then
            n_tr_weight_reg(nx1,ny1)=n_tr_weight_reg(nx1,ny1)+&
                 max(0.0_wp,1.0_wp-(abs(tracers(n)%d(3)-y(ng)%d(ny1)))&
                 /(y(ng)%d(ny2)-y(ng)%d(ny1)))
            n_tr_weight_reg(nx1,ny2)=n_tr_weight_reg(nx1,ny2)+&
                 max(0.0_wp,1.0_wp-(abs(tracers(n)%d(3)-y(ng)%d(ny2)))&
                 /(y(ng)%d(ny2)-y(ng)%d(ny1)))
         elseif (tracers(n)%d(1)==0.0_wp.and.ny1==ny2.and.nx1==nx2) then
            n_tr_weight_reg(nx1,ny1)=n_tr_weight_reg(nx1,ny1)+1.0_wp
         endif
      enddo
      ! Now I need to exchange partial sums that come from grid points 
      ! on the boundary with neighboring processors.  First do left and right
      ! exchanges, send and receiving the entire column of ghost points.  
      ! Ghost points at the top of my sub-domain will have tracer sums of 0, 
      ! because they lie completely outside my sub-domain (except when the 
      ! sub-domain includes the top boundary of the whole domain)
      ! Ghost points at the bottom of my subdomain will have the partial sum
      ! for all tracers in the cell above it, and the cells on the top right 
      ! and top left, unless we are at a sub-domain left or right boundary
      !
      ! Send the partial sum on my left boundary to my neighbor on the left, 
      ! receive partial sum from my neighbor on the right and put it in 
      ! partial__right
      call MPI_SENDRECV(n_tr_weight_dense(1,1:myny1),myny1,&
           MPI_DOUBLE_PRECISION,nbrleft,0,partial_dense_right(1:myny1),&
           myny1,MPI_DOUBLE_PRECISION,nbrright,0,comm2d,status,ierr)
      call MPI_SENDRECV(n_tr_weight_reg(1,1:myny1),myny1,&
           MPI_DOUBLE_PRECISION,nbrleft,1,partial_reg_right(1:myny1),&
           myny1,MPI_DOUBLE_PRECISION,nbrright,1,comm2d,status,ierr)
      ! Now add the recieved partial sum to n_tr_weight
      n_tr_weight_dense(mynx1-1,1:myny1)=&
           n_tr_weight_dense(mynx1-1,1:myny1)+&
           partial_dense_right(1:myny1)
      n_tr_weight_reg(mynx1-1,1:myny1)=&
           n_tr_weight_reg(mynx1-1,1:myny1)+&
           partial_reg_right(1:myny1)
      ! Send right and receive from left
      call MPI_SENDRECV(n_tr_weight_dense(mynx1,1:myny1),myny1,&
           MPI_DOUBLE_PRECISION,nbrright,2,partial_dense_left(1:myny1),&
           myny1,MPI_DOUBLE_PRECISION,nbrleft,2,comm2d,status,ierr)
      call MPI_SENDRECV(n_tr_weight_reg(mynx1,1:myny1),myny1,&
           MPI_DOUBLE_PRECISION,nbrright,3,partial_reg_left(1:myny1),&
           myny1,MPI_DOUBLE_PRECISION,nbrleft,3,comm2d,status,ierr)
      ! Now add the recieved partial sum to n_tr_weight
      n_tr_weight_dense(2,1:myny1)=n_tr_weight_dense(2,1:myny1)+&
           partial_dense_left(1:myny1)
      n_tr_weight_reg(2,1:myny1)=n_tr_weight_reg(2,1:myny1)+&
           partial_reg_left(1:myny1)
      !
      ! Now send my bottom row of ghost points to the process below me (and 
      ! receive ghost points from process above) for i=2,mynx-1; i.e. don't 
      ! include left and right ghost points, they aren't needed. 
      ! Only have to send below and receive from above, to complete the tracer
      ! sum on my topmost real grid point, because the top real grid point sits
      ! on the sub-domain boundary, and needs info from the process above. 
      ! The bottommost real grid point sits on cell length within my domain, 
      ! and therefore already sees all the tracers it's supposed to
      !
      ! Send below and receive from above
      call MPI_SENDRECV(n_tr_weight_dense(2:mynx1-1,1),mynx1-2,&
           MPI_DOUBLE_PRECISION,nbrbottom,6,partial_dense_top(2:mynx1-1),&
           mynx1-2,MPI_DOUBLE_PRECISION,nbrtop,6,comm2d,status,ierr)
      call MPI_SENDRECV(n_tr_weight_reg(2:mynx1-1,1),mynx1-2,&
           MPI_DOUBLE_PRECISION,nbrbottom,7,partial_reg_top(2:mynx1-1),&
           mynx1-2,MPI_DOUBLE_PRECISION,nbrtop,7,comm2d,status,ierr)
      ! Now add the recieved partial sum to n_tr_weight
      n_tr_weight_dense(2:mynx1-1,myny1-1)=&
           n_tr_weight_dense(2:mynx1-1,myny1-1)+&
           partial_dense_top(2:mynx1-1)
      n_tr_weight_reg(2:mynx1-1,myny1-1)=&
           n_tr_weight_reg(2:mynx1-1,myny1-1)+&
           partial_reg_top(2:mynx1-1)   
      ! Don't need to exchange with neighbors on the corners now, because the
      ! left-right exchanges coupled with the above-below exchange gives 
      ! complete tracer sums all around
      
      ! Now convert weighted tracer counts into compositional field
      do j=1,myny1
         ! We are including the bottom and top ghost points, because this 
         ! ensures that processes on the boundary calculate composition at 
         ! the top and bottom boundaries.  Everywhere else, the ghost point 
         ! value of composition (either top or bottom) will be wrong, but a 
         ! neighbor exchange will fix that after calling this subroutine
         do i=2,mynx1-1
            if (n_tr_weight_reg(i,j)/=0.0_wp.or.&
                 n_tr_weight_dense(i,j)/=0.0_wp) then
               chem(i,j)=n_tr_weight_dense(i,j)/(n_tr_weight_reg(i,j)+&
                    n_tr_weight_dense(i,j))
            endif
         enddo
      enddo
    END SUBROUTINE calc_chem_tr_ratio
! 
    SUBROUTINE calc_visc_tracer(tracers,mu,mu_t,alpha,temp,etae,etat,t0,&
         temp_cutoff,m_alpha,mu_jump,myntracers,nx,ny,sx,sy,asp,ntracers,&
         ntotal)
      USE btype
      USE mpi
      IMPLICIT NONE
      TYPE(ptr1d), INTENT(IN) :: tracers(:)
      REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: mu,mu_t
      REAL(WP), DIMENSION(:,:), INTENT(IN) :: alpha,temp
      INTEGER(I4B), INTENT(IN) :: myntracers,nx,ny,sx,sy,ntracers,ntotal
      REAL(WP), INTENT(IN) :: mu_jump,asp,etae,etat,t0,temp_cutoff,m_alpha
      INTEGER(I4B) :: i,j,n,nx1,ny1,mynx1,myny1,tracer_min,&
           visc_scale_factor,k
      INTEGER(I4B), DIMENSION(size(mu,1),size(mu,2)) :: trace_count
      REAL(WP) :: tracer_avg
      REAL(WP), DIMENSION(size(mu,1),size(mu,2)) :: visc_scale_factor2,smooth2
      trace_count=0
      mynx1=size(mu,1)
      myny1=size(mu,2)
      ! Loop through all tracers and find the ceiling 
      do n=1,myntracers
         ! First figure out the four nearest grid points for chem by 
         ! determining i of the grid points to the left and right, and j of 
         ! the grid points above and below 
         nx1=ceiling(tracers(n)%d(2)*(nx/asp)+3.0_wp-sx)
         ny1=ceiling(tracers(n)%d(3)*ny+3.0_wp-sy)
         trace_count(nx1,ny1)=trace_count(nx1,ny1)+1
      enddo
      call nbrex2d_int(trace_count,.false.,.false.,.false.,.false.,myid,&
           nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
      tracer_avg=real(ntracers)/real(ntotal)
      do j=1,myny1
         do i=1,mynx1
            visc_scale_factor=min(nint(trace_count(i,j)/tracer_avg),10)
            visc_scale_factor2(i,j)=(real(visc_scale_factor)/10.0_wp)*mu_jump
         enddo
      enddo
      ! Now smooth out visc scale factor using two passes of a weighted 
      ! relaxation routine
      ! Do 2 smoothing passes
      do k=1,1
         do j=2,myny1-1
            do i=2,mynx1-1
               smooth2(i,j)=0.2_wp*visc_scale_factor2(i,j)+&
                    0.2_wp*visc_scale_factor2(i-1,j)+&
                    0.2_wp*visc_scale_factor2(i+1,j)+&
                    0.2_wp*visc_scale_factor2(i,j-1)+&
                    0.2_wp*visc_scale_factor2(i,j+1)
!!$               smooth2(i,j)=0.6_wp*visc_scale_factor2(i,j)+&
!!$                    0.1_wp*visc_scale_factor2(i-1,j)+&
!!$                    0.1_wp*visc_scale_factor2(i+1,j)+&
!!$                    0.1_wp*visc_scale_factor2(i,j-1)+&
!!$                    0.1_wp*visc_scale_factor2(i,j+1)
            enddo
         enddo
         call nbrex2d(smooth2,.true.,.true.,.true.,.true.,myid,&
              nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
         do j=2,myny1-1
            do i=2,mynx1-1
               visc_scale_factor2(i,j)=0.2_wp*smooth2(i,j)+&
                    0.2_wp*smooth2(i-1,j)+0.2_wp*smooth2(i+1,j)+&
                    0.2_wp*smooth2(i,j-1)+0.2_wp*smooth2(i,j+1)
!!$               visc_scale_factor2(i,j)=0.6_wp*smooth2(i,j)+&
!!$                    0.1_wp*smooth2(i-1,j)+0.1_wp*smooth2(i+1,j)+&
!!$                    0.1_wp*smooth2(i,j-1)+0.1_wp*smooth2(i,j+1)
            enddo
         enddo
         call nbrex2d(visc_scale_factor2,.true.,.true.,.true.,.true.,myid,&
              nbrleft,nbrright,nbrtop,nbrbottom,comm2d)
      enddo
      do j=1,myny1
         do i=1,mynx1
            if (temp(i,j).lt.temp_cutoff) then
               mu(i,j) = max(visc_scale_factor2(i,j),1.0_wp)*&
                    exp(etae/(temp_cutoff+etat)-&
                    etae/(t0+etat))*alpha(i,j)**(-m_alpha)
               mu_t(i,j) = max(visc_scale_factor2(i,j),1.0_wp)*&
                    exp(etae/(temp_cutoff+etat)-&
                    etae/(t0+etat))
            else
               mu(i,j) = max(visc_scale_factor2(i,j),1.0_wp)*&
                    exp(etae/(temp(i,j)+etat)-etae/(t0+etat))*&
                    alpha(i,j)**(-m_alpha)
               mu_t(i,j) = max(visc_scale_factor2(i,j),1.0_wp)*&
                    exp(etae/(temp(i,j)+etat)-etae/(t0+etat))
            endif
         enddo
      enddo
    END SUBROUTINE calc_visc_tracer
  END PROGRAM simpler


! Old way I combined partial sums to grid points on the boundaries of a 
! sub-domain.  This was for composition defined at cell centers
!!$      ! Now I need to exchange partial sums that come from grid points 
!!$      ! on the boundary with neighboring processors.  
!!$      ! First do top, bottom, right, and left (easiest cases)
!!$      ! from i=3 to mynx-2 and j=3 to myny-2
!!$      !
!!$      ! Send the partial sum on my left boundary to my neighbor on the left, 
!!$      ! receive partial sum from my neighbor on the right and put it in 
!!$      ! partial__right
!!$      call MPI_SENDRECV(n_tr_weight_dense(1,2:myny1-1),myny1-2,&
!!$           MPI_DOUBLE_PRECISION,nbrleft,0,partial_dense_right(2:myny1-1),&
!!$           myny1-2,MPI_DOUBLE_PRECISION,nbrright,0,comm2d,status,ierr)
!!$      call MPI_SENDRECV(n_tr_weight_reg(1,2:myny1-1),myny1-2,&
!!$           MPI_DOUBLE_PRECISION,nbrleft,1,partial_reg_right(2:myny1-1),&
!!$           myny1-2,MPI_DOUBLE_PRECISION,nbrright,1,comm2d,status,ierr)
!!$      ! Now add the recieved partial sum to n_tr_weight
!!$      n_tr_weight_dense(mynx1-1,2:myny1-1)=&
!!$           n_tr_weight_dense(mynx1-1,2:myny1-1)+&
!!$           partial_dense_right(2:myny1-1)
!!$      n_tr_weight_reg(mynx1-1,2:myny1-1)=&
!!$           n_tr_weight_reg(mynx1-1,2:myny1-1)+&
!!$           partial_reg_right(2:myny1-1)
!!$      ! Send right and receive from left
!!$      call MPI_SENDRECV(n_tr_weight_dense(mynx1,2:myny1-1),myny1-2,&
!!$           MPI_DOUBLE_PRECISION,nbrright,2,partial_dense_left(2:myny1-1),&
!!$           myny1-2,MPI_DOUBLE_PRECISION,nbrleft,2,comm2d,status,ierr)
!!$      call MPI_SENDRECV(n_tr_weight_reg(mynx1,2:myny1-1),myny1-2,&
!!$           MPI_DOUBLE_PRECISION,nbrright,3,partial_reg_left(2:myny1-1),&
!!$           myny1-2,MPI_DOUBLE_PRECISION,nbrleft,3,comm2d,status,ierr)
!!$      ! Now add the recieved partial sum to n_tr_weight
!!$      n_tr_weight_dense(2,2:myny1-1)=n_tr_weight_dense(2,2:myny1-1)+&
!!$           partial_dense_left(2:myny1-1)
!!$      n_tr_weight_reg(2,2:myny1-1)=n_tr_weight_reg(2,2:myny1-1)+&
!!$           partial_reg_left(2:myny1-1)
!!$      ! Send above and receive from below
!!$      call MPI_SENDRECV(n_tr_weight_dense(2:mynx1-1,myny1),mynx1-2,&
!!$           MPI_DOUBLE_PRECISION,nbrtop,4,partial_dense_bot(2:mynx1-1),&
!!$           mynx1-2,MPI_DOUBLE_PRECISION,nbrbottom,4,comm2d,status,ierr)
!!$      call MPI_SENDRECV(n_tr_weight_reg(2:mynx1-1,myny1),mynx1-2,&
!!$           MPI_DOUBLE_PRECISION,nbrtop,5,partial_reg_bot(2:mynx1-1),&
!!$           mynx1-2,MPI_DOUBLE_PRECISION,nbrbottom,5,comm2d,status,ierr)
!!$      ! Now add the recieved partial sum to n_tr_weight
!!$      n_tr_weight_dense(2:mynx1-1,2)=n_tr_weight_dense(2:mynx1-1,2)+&
!!$           partial_dense_bot(2:mynx1-1)
!!$      n_tr_weight_reg(2:mynx1-1,2)=n_tr_weight_reg(2:mynx1-1,2)+&
!!$           partial_reg_bot(2:mynx1-1)
!!$      ! Send below and receive from above
!!$      call MPI_SENDRECV(n_tr_weight_dense(2:mynx1-1,1),mynx1-2,&
!!$           MPI_DOUBLE_PRECISION,nbrbottom,6,partial_dense_top(2:mynx1-1),&
!!$           mynx1-2,MPI_DOUBLE_PRECISION,nbrtop,6,comm2d,status,ierr)
!!$      call MPI_SENDRECV(n_tr_weight_reg(2:mynx1-1,1),mynx1-2,&
!!$           MPI_DOUBLE_PRECISION,nbrbottom,7,partial_reg_top(2:mynx1-1),&
!!$           mynx1-2,MPI_DOUBLE_PRECISION,nbrtop,7,comm2d,status,ierr)
!!$      ! Now add the recieved partial sum to n_tr_weight
!!$      n_tr_weight_dense(2:mynx1-1,myny1-1)=&
!!$           n_tr_weight_dense(2:mynx1-1,myny1-1)+&
!!$           partial_dense_top(2:mynx1-1)
!!$      n_tr_weight_reg(2:mynx1-1,myny1-1)=&
!!$           n_tr_weight_reg(2:mynx1-1,myny1-1)+&
!!$           partial_reg_top(2:mynx1-1)   
!!$      ! Now get partial sums from neighbors on the corners; only needed for 
!!$      ! corner point of each sub-domain
!!$      ! Send to top-right and receive from bottom-left
!!$      call MPI_SENDRECV(n_tr_weight_dense(mynx1,myny1),1,&
!!$           MPI_DOUBLE_PRECISION,nbrtopright,8,partial_dense_bl,1,&
!!$           MPI_DOUBLE_PRECISION,nbrbotleft,8,comm2d,status,ierr)
!!$      call MPI_SENDRECV(n_tr_weight_reg(mynx1,myny1),1,&
!!$           MPI_DOUBLE_PRECISION,nbrtopright,9,partial_reg_bl,1,&
!!$           MPI_DOUBLE_PRECISION,nbrbotleft,9,comm2d,status,ierr)
!!$      n_tr_weight_dense(2,2)=n_tr_weight_dense(2,2)+partial_dense_bl
!!$      n_tr_weight_reg(2,2)=n_tr_weight_reg(2,2)+partial_reg_bl
!!$      ! Send to top-left and receive from bottom-right
!!$      call MPI_SENDRECV(n_tr_weight_dense(1,myny1),1,&
!!$           MPI_DOUBLE_PRECISION,nbrtopleft,10,partial_dense_br,1,&
!!$           MPI_DOUBLE_PRECISION,nbrbotright,10,comm2d,status,ierr)
!!$      call MPI_SENDRECV(n_tr_weight_reg(1,myny1),1,&
!!$           MPI_DOUBLE_PRECISION,nbrtopleft,11,partial_reg_br,1,&
!!$           MPI_DOUBLE_PRECISION,nbrbotright,11,comm2d,status,ierr)
!!$      n_tr_weight_dense(mynx1-1,2)=n_tr_weight_dense(mynx1-1,2)+&
!!$           partial_dense_br
!!$      n_tr_weight_reg(mynx1-1,2)=n_tr_weight_reg(mynx1-1,2)+&
!!$           partial_reg_br
!!$      ! Send to bottom-left and receive from top-right
!!$      call MPI_SENDRECV(n_tr_weight_dense(1,1),1,&
!!$           MPI_DOUBLE_PRECISION,nbrbotleft,12,partial_dense_tr,1,&
!!$           MPI_DOUBLE_PRECISION,nbrtopright,12,comm2d,status,ierr)
!!$      call MPI_SENDRECV(n_tr_weight_reg(1,1),1,&
!!$           MPI_DOUBLE_PRECISION,nbrbotleft,13,partial_reg_tr,1,&
!!$           MPI_DOUBLE_PRECISION,nbrtopright,13,comm2d,status,ierr)
!!$      n_tr_weight_dense(mynx1-1,myny1-1)=&
!!$           n_tr_weight_dense(mynx1-1,myny1-1)+partial_dense_tr
!!$      n_tr_weight_reg(mynx1-1,myny1-1)=&
!!$           n_tr_weight_reg(mynx1-1,myny1-1)+partial_reg_tr
!!$      ! Send to bottom-right and receive from top-left
!!$      call MPI_SENDRECV(n_tr_weight_dense(mynx1,1),1,&
!!$           MPI_DOUBLE_PRECISION,nbrbotleft,14,partial_dense_tl,1,&
!!$           MPI_DOUBLE_PRECISION,nbrtopright,14,comm2d,status,ierr)
!!$      call MPI_SENDRECV(n_tr_weight_reg(mynx1,1),1,&
!!$           MPI_DOUBLE_PRECISION,nbrbotleft,15,partial_reg_tl,1,&
!!$           MPI_DOUBLE_PRECISION,nbrtopright,15,comm2d,status,ierr)
!!$      n_tr_weight_dense(2,myny1-1)=n_tr_weight_dense(2,myny1-1)+&
!!$           partial_dense_tl
!!$      n_tr_weight_reg(2,myny1-1)=n_tr_weight_reg(2,myny1-1)+&
!!$           partial_reg_tl



! Old inefficient way I did the tracer ratio method for converting tracer 
! distribution into composition at cell centers
!!$      m=0
!!$      do i=2,mynx(ng)-1
!!$         do j=2,myny(ng)-1
!!$            r_bound=x(ng)%d(i)+deltax
!!$            l_bound=x(ng)%d(i)-deltax
!!$            n_bound=y(ng)%d(j)+deltay
!!$            s_bound=y(ng)%d(j)-deltay
!!$            n_tr_weight_dense=0.0_wp
!!$            n_tr_weight_reg=0.0_wp
!!$            do n=1,myntracers
!!$               if (tracers(n)%d(2)<=r_bound.and.tracers(n)%d(2)>=l_bound&
!!$                    .and.tracers(n)%d(3)<=n_bound.and.tracers(n)%d(3)>=&
!!$                    s_bound.and.tracers(n)%d(1)==1.0_wp) then
!!$                  n_tr_weight_dense=n_tr_weight_dense+max(0.0_wp,&
!!$                       1.0_wp-(abs(tracers(n)%d(2)-x(ng)%d(i)))/deltax)*&
!!$                       max(0.0_wp,1.0_wp-(abs(tracers(n)%d(3)-y(ng)%d(j)))&
!!$                       /deltay)
!!$               elseif (tracers(n)%d(2)<=r_bound.and.tracers(n)%d(2)>=l_bound&
!!$                    .and.tracers(n)%d(3)<=n_bound.and.tracers(n)%d(3)>=&
!!$                    s_bound.and.tracers(n)%d(1)==0.0_wp) then
!!$                  n_tr_weight_reg=n_tr_weight_reg+max(0.0_wp,&
!!$                       1.0_wp-(abs(tracers(n)%d(2)-x(ng)%d(i)))/deltax)*&
!!$                       max(0.0_wp,1.0_wp-(abs(tracers(n)%d(3)-y(ng)%d(j)))&
!!$                       /deltay)
!!$               endif
!!$            enddo
!!$            if (n_tr_weight_reg/=0.0_wp.or.n_tr_weight_dense/=0.0_wp) then
!!$               chem(i,j)=n_tr_weight_dense/&
!!$                    (n_tr_weight_reg+n_tr_weight_dense)
!!$            endif
!!$            if (n_tr_weight_reg==0.0_wp.and.n_tr_weight_dense==0.0_wp) then
!!$               m=m+1
!!$               fix_i(m)=i
!!$               fix_j(m)=j
!!$            endif
!!$         enddo
!!$      enddo
!!$      do i=1,m
!!$         chem(fix_i(i),fix_j(i))=0.25_wp*(chem(fix_i(i)-1,fix_j(i))+&
!!$              chem(fix_i(i)+1,fix_j(i))+chem(fix_i(i),fix_j(i)-1)+&
!!$              chem(fix_i(i),fix_j(i)+1))
!!$      enddo
