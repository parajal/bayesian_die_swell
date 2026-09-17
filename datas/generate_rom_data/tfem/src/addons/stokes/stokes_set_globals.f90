
! Copyright (C) 2006-2024 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! This module provides the set_globals_stokes_xxx and (un)set_stokes_elem for
! the Stokes equation in 2D and 3D such that they can be used by other
! element routines


module stokes_set_globals_m

  use stokes_elements_generic_m, only: set_globals_stokes_vp, &
          set_globals_stokes_vp_boun, set_globals_stokes_l_boun, &
          set_globals_stokes_object, &
          set_stokes_elem, unset_stokes_elem, set_globals_stokes_tALE, &
          set_stokes_shape_function_global

  implicit none

end module stokes_set_globals_m

