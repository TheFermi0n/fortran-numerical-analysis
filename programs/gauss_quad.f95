program gauss_quad
  use tools

  implicit none
  real(dp) :: res
  integer :: a, b, k

  res = 0.0
  a = -2; b = 4
  k = 5

  call g_quad(k, a, b, res)

  print *, ''
  print *, 'result:'
  write (*, form_result2r) res, '-- the required result of the integration'
end program gauss_quad
