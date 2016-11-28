! Subroutine for exchanging tracer data in a given direction (direction is 
! set by the choice of nbr_to and nbr_from; they should be a pair such that all
! processors are sending and receiving in the same direction, e.g. everyone
! sends up and receives from below)
! 
! OUT is the array of pointers, tracers, which contains the updated tracer 
! position data for a processor's sub-domain; this routine will place the 
! tracer data from tracers coming into "my" domain into the tracers 
! data structure
!
! INput is the array of pointers containing the old tracer data, oldtracers
!
! Other input needed: send_count, the number of tracers being sent in the 
! specified direction (e.g. my_tracers_out(1) in the main program, 
! for the example of sending up and receiving from below);
! recv_count, number of incoming tracers (e.g. my_tracers_in(5));
! tracers_out, array that gives the 'i' coordinate for tracer data leaving my 
! sub-domain (e.g. tr_out_top in our example);
! nbr_to & nbr_from; the nbr to send to, and receive from, respectively
! nstart, the index to start adding the new tracers into the updated tracer 
! data
SUBROUTINE nbrex2d_tracer_dat(tracers,oldtracers,send_count,recv_count,&
     tracers_out,nbr_to,nbr_from,nstart,comm2d)
  USE mpi
  USE btype
  IMPLICIT NONE
  TYPE(ptr1d), INTENT(INOUT) :: tracers(:)
  TYPE(ptr1d), INTENT(IN) :: oldtracers(:)
  INTEGER(I4B), INTENT(IN) :: nbr_to,nbr_from,comm2d,nstart,send_count,&
       recv_count
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: tracers_out
  INTEGER(I4B) :: status1(MPI_STATUS_SIZE),ierr,i
  REAL(WP), DIMENSION(send_count,3) :: send_dat
  REAL(WP), DIMENSION(recv_count,3) :: recv_dat
  !
  ! Will get number of tracers being sent and received in the direction 
  ! dictated by input
  !
  ! First pick out the tracers being sent and put in a dummy "send arry"
  do i=1,send_count
     send_dat(i,:)=oldtracers(tracers_out(i))%d
  enddo
  if (recv_count==0) then
     call MPI_SENDRECV(send_dat,send_count*3,MPI_DOUBLE_PRECISION,nbr_to,0,&
          recv_dat,recv_count*3,MPI_DOUBLE_PRECISION,MPI_PROC_NULL,0,comm2d,&
          status1,ierr)
  elseif (send_count==0) then
     call MPI_SENDRECV(send_dat,send_count*3,MPI_DOUBLE_PRECISION,&
          MPI_PROC_NULL,0,recv_dat,recv_count*3,MPI_DOUBLE_PRECISION,nbr_from,&
          0,comm2d,status1,ierr)
  else
     call MPI_SENDRECV(send_dat,send_count*3,MPI_DOUBLE_PRECISION,&
          nbr_to,0,recv_dat,recv_count*3,MPI_DOUBLE_PRECISION,nbr_from,&
          0,comm2d,status1,ierr)
  endif
  !
  ! Now put the recieved data into the full tracers data structure
  do i=1,recv_count
     allocate(tracers(i+nstart-1)%d(3))
     tracers(i+nstart-1)%d=recv_dat(i,:)
  enddo

END SUBROUTINE nbrex2d_tracer_dat
