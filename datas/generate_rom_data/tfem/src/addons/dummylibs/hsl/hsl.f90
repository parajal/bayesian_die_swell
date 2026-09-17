subroutine hsl_error(routine)
  character(len=*) :: routine
  write(*,'(/a,a/)') ' HSL library not installed. Routine missing: ', routine
  stop
end subroutine hsl_error

subroutine ma41id
  call hsl_error('ma41id')
end

subroutine ma41ad
  call hsl_error('ma41ad')
end

subroutine ma57id
  call hsl_error('ma57id')
end

subroutine ma57ad
  call hsl_error('ma57ad')
end

subroutine ma57bd
  call hsl_error('ma57bd')
end

subroutine ma57cd
  call hsl_error('ma57cd')
end

subroutine ma57dd
  call hsl_error('ma57dd')
end

subroutine ea23ad
  call hsl_error('ea23ad')
end
