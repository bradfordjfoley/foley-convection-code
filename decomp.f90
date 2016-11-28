! Do decomposition, following mpe_dcomp algorithm from using mpi
! Ensures 'fair' decomposition for processors that don't evenly divide grid

SUBROUTINE decomp(s,e,dims,coords,n,bdy1,bdy2)
  USE btype
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: n,dims,coords
  INTEGER(I4B), INTENT(OUT) :: s,e
  LOGICAL, INTENT(OUT) :: bdy1,bdy2
  INTEGER(I4B) :: nlocal,deficit

  ! decomposition for pressure nodes
  nlocal=n/dims
  s=coords*nlocal+1
  deficit=mod(n,dims)
  s=s+min(coords,deficit)
  if (coords < deficit) then
     nlocal = nlocal + 1
  endif
  e=s+nlocal-1 
  if (e > n.or. coords == dims) then 
     e=n
  endif

  ! Now determine if my subdomain has a boundary
  if (s == 1) then
     bdy1 = .true.
  else 
     bdy1 = .false.
  endif
  if (e == n) then
     bdy2 = .true.
  else
     bdy2 = .false.
  endif
END SUBROUTINE decomp
