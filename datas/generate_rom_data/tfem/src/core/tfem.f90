
! Copyright (C) 2004-2010 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! A module for TFEM that `uses' all commonly used core modules in a main
! program in one go.
! Note that a solver is missing. This means that a solver must always
! be explicitly used.

module tfem_m

  use kind_defs_m
  use mesh_m
  use meshgen_m
  use problem_m
  use system_m
  use element_defs_m
  use vector_m
  use elvector_m
  use eltree_m
  use postprocessing_m

  implicit none

end module tfem_m
