program mc_pi
  use tools

  implicit none
  integer :: sq_points, cr_points
  real(dp) :: result, dev

  sq_points = 10000 !number of points inside the square

  ! use grphmonte_pi() subroutine
  ! to get those data points in a file and use the grphing tools in graph/
  call monte_pi(sq_points, cr_points, result)

  dev = abs(pi - result)

  print *, ''
  print *, 'result:'
  write (*, form_result2rp) result, '-- approximate value of pi'
  write (*, form_result2rp) pi, '-- acutal value of pi'
  write (*, form_result2rp) dev, '-- deviation from actual value'
  write (*, form_result2i) cr_points, '-- points inside circle'
end program mc_pi
