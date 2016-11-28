! Exchanges information with neighbors so each process gets updated
! ghost points, including corners
SUBROUTINE nbrex2d(ue,tfsl,bfsl,lfsl,rfsl,myid,nbrleft,nbrright,&
     nbrtop,nbrbottom,comm2d)
  USE mpi
  USE btype
  IMPLICIT NONE
  REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: ue
  INTEGER(I4B), INTENT(IN) :: myid,nbrleft,nbrright,comm2d,nbrtop,nbrbottom
  INTEGER(I4B) :: nx1,ny1,status(MPI_STATUS_SIZE),ierr
  LOGICAL, INTENT(IN) :: tfsl,bfsl,lfsl,rfsl
  nx1=size(ue,1)
  ny1=size(ue,2)
  call MPI_SENDRECV(ue(2:nx1-1,ny1-1),nx1-2,&
       MPI_DOUBLE_PRECISION,nbrtop,0,ue(2:nx1-1,1),&
       nx1-2,MPI_DOUBLE_PRECISION,nbrbottom,0,&
       comm2d,status,ierr)
  call MPI_SENDRECV(ue(2:nx1-1,2),nx1-2,&
       MPI_DOUBLE_PRECISION,nbrbottom,1,ue(2:nx1-1,ny1),&
       nx1-2,MPI_DOUBLE_PRECISION,nbrtop,1,&
       comm2d,status,ierr)
  ! Now can send entire columns in x direction 
  ! (including ghost points in y)
  call MPI_SENDRECV(ue(nx1-1,1:ny1),ny1,&
       MPI_DOUBLE_PRECISION,nbrright,2,ue(1,1:ny1),&
       ny1,MPI_DOUBLE_PRECISION,nbrleft,2,&
       comm2d,status,ierr)
  call MPI_SENDRECV(ue(2,1:ny1),ny1,&
       MPI_DOUBLE_PRECISION,nbrleft,3,ue(nx1,1:ny1),&
       ny1,MPI_DOUBLE_PRECISION,nbrright,3,&
       comm2d,status,ierr)
  if (tfsl.and.nbrtop.lt.0) then
     ue(1:nx1,ny1)=ue(1:nx1,ny1-1)
  endif
  if (bfsl.and.nbrbottom.lt.0) then
     ue(1:nx1,1)=ue(1:nx1,2)
  endif
  if (lfsl.and.nbrleft.lt.0) then
     ue(1,1:ny1)=ue(2,1:ny1)
  endif
  if (rfsl.and.nbrright.lt.0) then
     ue(nx1,1:ny1)=ue(nx1-1,1:ny1)
  endif
     
END SUBROUTINE nbrex2d
