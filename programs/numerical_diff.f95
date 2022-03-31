program nmdf
  use tools
  ! ----------------------------------------------------------------
  ! Newton's differentiation:
  !   Calculates the differentiation of y at a given x using newton's
  !   interpolation polynomial and the difference table calculated
  !   by the subroutine newinter()
  ! inputs:
  !   points    -- number of points for interpolation
  !   x()       -- the arrays of x's
  !   y()       -- the arrays of corresponding y's (f(x)'s)
  !   x_point   -- the value of x at which f'(x) is to be calculated
  ! output:
  !   res       -- the result of the derivation
  !   poly      -- the used interpolation polynomial
  ! Note:
  !   - used the subroutine newinter()
  !   - used forward difference table and modified the formula a bit
  !     to calculate the backward differentiation
  !   - used module "tools"
  ! env details:
  !   compiler -- gfortran v10.2.0
  !   OS       -- macOS 12.0 x86_64F
  !   Kernel   -- Darwin 21.0.0
  ! ----------------------------------------------------------------
  implicit none
  integer, parameter :: points = 8
  real(dp), dimension(points) :: x, y
  real(dp) :: x_point, res
  integer :: i
  character(len=10) :: poly

  ! given data points
  data(x(i), i=1, 8)/1.21, 1.31, 1.41, 1.51, 1.61, 1.71, 1.81, 1.91/
  data(y(i), i=1, 8)/4.43, 6.84, 10.06, 12.67, 15.73, 19.87, 21.86, 24.61/

  x_point = 1.61 ! the calculation point
  res = 0.0

  call neum_diff(points, x, y, x_point, poly, res)

  ! result printing
  print *, ''
  print *, 'result:'
  write (*, form_result2c) poly, '-- interpolation polynomial'
  write (*, form_result2r) x_point, '-- given point'
  write (*, form_result2r) res, '-- value of derivative'
end program nmdf
