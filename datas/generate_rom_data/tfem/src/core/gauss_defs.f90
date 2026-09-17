
! Copyright (C) 2007-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Definitions for Gauss integration

module gauss_defs_m

  use kind_defs_m

  implicit none


! type definition of a (gauss) integration rule

  type gauss_t

!   Global shape of the element:
!     line
!     quadrilateral
!     triangle
!     hexahedron
!     tetrahedron
!     prism
!     pyramid
    character (len=13) :: globalshape = ''

!   integration type:
!     0: standard Gauss-Legendre (using analytical expressions)
!     1: Gauss-Legendre-Lobatto
!     2: mapped Gauss (from quad->triangle, hexahedron->pyramid)
!     3: standard Gauss-Legendre (using tabulated or numerical results)
    integer :: inttype = 0

!   intrule: integration rule
!     inttype=0:
!        - triangle/tetrahedron/prism:
!            number of Gauss points (prism: base triangle only)
!        - line/quadrilateral/hexahedron:
!            number of Gauss points in one dimension
!     inttype=1:
!        - line/quadrilateral/hexahedron:
!            number of Gauss points in one dimension
!     inttype=2:
!        - triangle: number of Gauss points in one dimension
!        - pyramid: number of Gauss points in one dimension (quad base only)
!     inttype=3:
!        - triangle/tetrahedron/prism/pyramid:
!            maximum polynomial order p for exact integration
!            (prism: base triangle only)
!        - line/quadrilateral/hexahedron:
!            number of Gauss points in one dimension
    integer :: intrule = 0

!   intrule2: secondary integration rule
!     inttype=0:
!        - prism: number of Gauss points in height direction
!     inttype=2:
!        - pyramid: number of Gauss points in height direction
!     inttype=3:
!        - prism: number of Gauss points in height direction
    integer :: intrule2 = 0

!   number of subdomains in one dimension for a composite/generalized rule
    integer :: nsubint = 1

  end type gauss_t

end module gauss_defs_m
