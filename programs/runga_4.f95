program rk_4
  use tools

  implicit none
  real(dp) :: x0, x_point, y0, h

  x0 = 0.0; y0 = 1.0; h = 0.05
  x_point = 0.25 !the value of x at which y is to be calculated

  write (*, form_header3) 'itn', 'x', 'y'
  write (*, line3) '-------', '-------', '-------'

  call rungkut_4(x0, y0, x_point, h)

  print *, ''
  print *, 'result:'
  write (*, form_result2r) x0, '-- given value of x'
  write (*, form_result2r) y0, '-- the required value of y'
end program rk_4
