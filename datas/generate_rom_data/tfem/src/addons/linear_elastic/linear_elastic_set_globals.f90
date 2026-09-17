
! Copyright (C) 2020-2020 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! This module provides the set_globals_linear_elastic_xxx and
! (un)set_linear_elastic_elem for the linear elastic equation in 2D and 3D
! such that they can be used by other element routines


module linear_elastic_set_globals_m

  use linear_elastic_elements_generic_m, only: set_globals_linear_elastic_up, &
      set_globals_linear_elastic_up_boun, set_globals_linear_elastic_l_boun, &
      set_linear_elastic_elem, unset_linear_elastic_elem

  implicit none

end module linear_elastic_set_globals_m

