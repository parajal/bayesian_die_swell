! ddot replacement for SPARSKIT

function distdot( n, x, ix, y, iy )

  integer :: n, ix, iy
  double precision :: distdot, x(*), y(*), ddot

  distdot = ddot(n,x,ix,y,iy)

end function distdot
