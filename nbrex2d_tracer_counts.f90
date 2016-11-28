! Exchanges the number of tracers leaving/entering a processor's domain between
! all neighboring processors
SUBROUTINE nbrex2d_tracer_counts(tracers_out,tracers_in,myid,nbrleft,&
     nbrright,nbrtop,nbrbottom,nbrtopleft,nbrtopright,&
     nbrbotleft,nbrbotright,comm2d)
  USE mpi
  USE btype
  IMPLICIT NONE
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: tracers_out
  INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: tracers_in
  INTEGER(I4B), INTENT(IN) :: myid,nbrleft,nbrright,comm2d,nbrtop,nbrbottom,&
       nbrtopleft,nbrtopright,nbrbotleft,nbrbotright
  INTEGER(I4B) :: status(MPI_STATUS_SIZE),ierr
  !
  ! Exchange number of tracers leaving/entering my domain between top and
  ! bottom neighbors
  call MPI_SENDRECV(tracers_out(1),1,MPI_INTEGER,nbrtop,0,tracers_in(5),1,&
       MPI_INTEGER,nbrbottom,0,comm2d,status,ierr)
  call MPI_SENDRECV(tracers_out(5),1,MPI_INTEGER,nbrbottom,1,tracers_in(1),1,&
       MPI_INTEGER,nbrtop,1,comm2d,status,ierr)
  ! Left and right neighbors
  call MPI_SENDRECV(tracers_out(3),1,MPI_INTEGER,nbrright,2,tracers_in(7),1,&
       MPI_INTEGER,nbrleft,2,comm2d,status,ierr)
  call MPI_SENDRECV(tracers_out(7),1,MPI_INTEGER,nbrleft,3,tracers_in(3),1,&
       MPI_INTEGER,nbrright,3,comm2d,status,ierr)
  ! Top-right and bottom-left neighbors
  call MPI_SENDRECV(tracers_out(2),1,MPI_INTEGER,nbrtopright,4,tracers_in(6),&
       1,MPI_INTEGER,nbrbotleft,4,comm2d,status,ierr)
  call MPI_SENDRECV(tracers_out(6),1,MPI_INTEGER,nbrbotleft,5,tracers_in(2),&
       1,MPI_INTEGER,nbrtopright,5,comm2d,status,ierr)
  ! Top-left and bottom-right neighbors
  call MPI_SENDRECV(tracers_out(8),1,MPI_INTEGER,nbrtopleft,6,tracers_in(4),&
       1,MPI_INTEGER,nbrbotright,6,comm2d,status,ierr)
  call MPI_SENDRECV(tracers_out(4),1,MPI_INTEGER,nbrbotright,7,tracers_in(8),&
       1,MPI_INTEGER,nbrtopleft,7,comm2d,status,ierr)
END SUBROUTINE nbrex2d_tracer_counts
