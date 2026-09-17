
! Copyright (C) 2021-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Nodal numbering of elements in gmsh relative to the structured tfem nodal
! numbering for high order elements.
!
! The routines optionally produce the index arrays, that are defined by
!
!       tfem_node = gmsh_index ( gmsh_node )
!
! and
!
!       gmsh_node = gmsh_inv_index ( tfem_node )
!
! where tfem_node is according to the structured nodal numbering as given in
! Section 4.2 of the tfem userguide. The gmsh numbering is hierarchic and can
! be found in Section 9.2 of the Gmsh Reference Manual.
!
! The index arrays can be used as follows. Suppose we have nodal values that
! are stored in the arrays u_gmsh and u_tfem according to the gmsh
! and structured tfem numbering, respectively. Then we can use the following
! statements for transfering the data:
! a. Transfering from gmsh -> tfem
!
!       u_tfem = u_gmsh ( gmsh_inv_index )
!       u_tfem ( gmsh_index ) = u_gmsh
!
! b. Transfering from tfem -> gmsh
!
!       u_gmsh ( gmsh_inv_index ) = u_tfem
!       u_gmsh = u_tfem ( gmsh_index )
!

module meshgen_gmsh_index_m

  implicit none

contains


! line element

  subroutine gmsh_index_line ( p, gmsh_index, gmsh_inv_index )

!   polynomial order (p>=1)
    integer, intent(in) :: p

!   index array for reordering nodes according to:
!       tfem_node = gmsh_index ( gmsh_node )
    integer, dimension(:), intent(out), optional :: gmsh_index

!   index array for reordering nodes according to:
!       gmsh_node = gmsh_inv_index ( tfem_node )
    integer, dimension(:), intent(out), optional :: gmsh_inv_index

    integer :: i

    if ( present(gmsh_index) ) then
      gmsh_index = [ 1, p+1, ( i, i = 2, p ) ]
!      print *, gmsh_index
    end if

    if ( present(gmsh_inv_index) ) then
      gmsh_inv_index = [ 1, ( i, i = 3, p+1 ), 2 ]
!      print *, gmsh_inv_index
    end if

  end subroutine gmsh_index_line


! triangle element

  subroutine gmsh_index_triangle ( p, gmsh_index, gmsh_inv_index )

!   polynomial order (p>=1)
    integer, intent(in) :: p

!   index array for reordering nodes according to:
!       tfem_node = gmsh_index ( gmsh_node )
    integer, dimension(:), intent(out), optional :: gmsh_index

!   index array for reordering nodes according to:
!       gmsh_node = gmsh_inv_index ( tfem_node )
    integer, dimension(:), intent(out), optional :: gmsh_inv_index

    integer :: mat(p+1,p+1), nr, r, s, i, j, k, b, e, nne

!   number of edge rings
    nr = (p+2)/3

!   set nodal point counter to zero
    k = 0

!   loop over all edge rings

    do r = 1, nr

      s = p+3-2*r ! end vertex index in matrix

!     three vertex points
      mat(r,r) = k+1
      mat(s,r) = k+2
      mat(r,s) = k+3
      k = k + 3

!     internal nodes on three edges
      nne = p+2-3*r   ! number of edge nodes
      b = r + 1       ! begin edge node index in matrix
      e = p+2-2*r     ! end edge node index in matrix

!     edge 1
      mat(b:e,r) = [ (k + i, i = 1, nne) ]
      k = k + nne

!     edge 2
      do i = 1, nne
        mat(e-i+1,b+i-1) = k + i
      end do
      k = k + nne

!     edge 3
      mat(r,e:b:-1) = [ (k + i, i = 1, nne) ]
      k = k + nne

    end do

!   center node
    if ( mod(p,3) == 0 ) mat(nr+1,nr+1) = k + 1

!    do i = p+1,1,-1
!      print *, mat(1:p+2-i,i)
!    end do

    if ( present(gmsh_index) ) then

!     fill gmsh_index

      k = 0
      do j = 1, p+1
        do i = 1, p+2-j
          k = k + 1
          gmsh_index(mat(i,j)) = k
        end do
      end do

!      print *, gmsh_index

     end if

    if ( present(gmsh_inv_index) ) then

!     fill gmsh_inv_index

      k = 0
      do j = 1, p+1
        do i = 1, p+2-j
          k = k + 1
          gmsh_inv_index(k) = mat(i,j)
        end do
      end do

!      print *, gmsh_inv_index

     end if

  end subroutine gmsh_index_triangle


! quadrilateral element

  subroutine gmsh_index_quadrilateral ( p, gmsh_index, gmsh_inv_index )

!   polynomial order (p>=1)
    integer, intent(in) :: p

!   index array for reordering nodes according to:
!       tfem_node = gmsh_index ( gmsh_node )
    integer, dimension(:), intent(out), optional :: gmsh_index

!   index array for reordering nodes according to:
!       gmsh_node = gmsh_inv_index ( tfem_node )
    integer, dimension(:), intent(out), optional :: gmsh_inv_index

    integer :: mat(p+1,p+1), nr, r, s, i, j, k, b, e, nne

!   number of edge rings
    nr = (p+1)/2

!   set nodal point counter to zero
    k = 0

!   loop over all edge rings

    do r = 1, nr

      s = p+2-r ! end vertex index in matrix

!     four vertex points
      mat(r,r) = k+1
      mat(s,r) = k+2
      mat(s,s) = k+3
      mat(r,s) = k+4
      k = k + 4

!     internal nodes on four edges
      nne = p+1-2*r   ! number of edge nodes
      b = r + 1       ! begin edge node index in matrix
      e = p+1-r       ! end edge node index in matrix

!     edge 1
      mat(b:e,r) = [ (k + i, i = 1, nne) ]
      k = k + nne

!     edge 2
      mat(s,b:e) = [ (k + i, i = 1, nne) ]
      k = k + nne

!     edge 3
      mat(e:b:-1,s) = [ (k + i, i = 1, nne) ]
      k = k + nne

!     edge 4
      mat(r,e:b:-1) = [ (k + i, i = 1, nne) ]
      k = k + nne

    end do

!   center node
    if ( mod(p,2) == 0 ) mat(nr+1,nr+1) = k + 1

!    do i = p+1,1,-1
!      print *, mat(:,i)
!    end do

    if ( present(gmsh_index) ) then

!     fill gmsh_index

      k = 0
      do j = 1, p+1
        do i = 1, p+1
          k = k + 1
          gmsh_index(mat(i,j)) = k
        end do
      end do

!     print *, gmsh_index

    end if

    if ( present(gmsh_inv_index) ) then

!     fill gmsh_inv_index

      k = 0
      do j = 1, p+1
        do i = 1, p+1
          k = k + 1
          gmsh_inv_index(k) = mat(i,j)
        end do
      end do

!     print *, gmsh_inv_index

    end if

  end subroutine gmsh_index_quadrilateral


! hexahedron element

  subroutine gmsh_index_hexahedron ( p, gmsh_index, gmsh_inv_index )

!   polynomial order (p>=1)
    integer, intent(in) :: p

!   index array for reordering nodes according to:
!       tfem_node = gmsh_index ( gmsh_node )
    integer, dimension(:), intent(out), optional :: gmsh_index

!   index array for reordering nodes according to:
!       gmsh_node = gmsh_inv_index ( tfem_node )
    integer, dimension(:), intent(out), optional :: gmsh_inv_index

    write(*,'(/a/)') &
      'Error in gmsh_index_hexahedron: not yet implemented.'
    stop

  end subroutine gmsh_index_hexahedron

end module meshgen_gmsh_index_m
