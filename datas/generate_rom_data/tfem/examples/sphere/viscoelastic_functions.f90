module viscoelastic_functions_m

  use kind_defs_m

  implicit none

contains


! inflow boundary conditions for the conformation tensor in DG

  function cinflow ( n, x, un )

!   number of components in the function
    integer, intent(in) :: n

!   coordinates of the point
    real(dp), intent(in), dimension(:) :: x

!   normal velocity (should be negative)

    real(dp), intent(in) :: un

!   the inflow value of the c tensor at the point x
    real(dp), dimension(n) :: cinflow

    cinflow = 0

  end function cinflow

end module viscoelastic_functions_m
