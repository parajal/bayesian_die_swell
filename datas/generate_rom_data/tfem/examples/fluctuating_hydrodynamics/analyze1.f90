! Program to analyze the data in out1

program analyze1

use kind_defs_m
implicit none

integer :: i, ios, n
real(dp) :: s(3,3), d(3)

open(unit=10,file='out1')

n= 0
s = 0
do
  read(unit=10,fmt=*,iostat=ios) d
  if ( ios /= 0 ) exit
  n = n + 1
  do i = 1, 3
    s(i,:) = s(i,:) + d(i)*d
  end do
end do
s = s/n
print *, 'estimated diffusion coefficients:'
print *, s/2
print *, 'estimated standard deviation of diagonal:'
do i = 1, 3
  d(i) = s(i,i)
end do
print *, sqrt(2._dp*(n-1))/n*d/2

end program analyze1
