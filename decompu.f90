! Do decomposition, following mpe_dcomp algorithm from using mpi
! Ensures 'fair' decomposition for processors that don't evenly divide grid
!
! Decomposition for staggered grid (u in x direction, w in z direction)
SUBROUTINE decompu(s,e,dims,coords,n)
  USE btype
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: n,dims,coords
  INTEGER(I4B), INTENT(OUT) :: s,e
  INTEGER(I4B) :: nlocal,deficit

  ! decomposition for pressure nodes
  nlocal=(n-1)/dims
  s=coords*nlocal+2
  deficit=mod(n-1,dims)
  s=s+min(coords,deficit)
  if (coords < deficit) then
     nlocal = nlocal + 1
  endif
  e=s+nlocal-1 
  if (e > n.or. coords == dims) then 
     e=n
  endif
END SUBROUTINE decompu
