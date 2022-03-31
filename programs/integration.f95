program int_mc
  use tools

  implicit none
  real(dp) :: res
  integer :: a, b, sample_size

  a = 1; b = 2 ! given lower and upper limit respectively
  sample_size = 10000 ! the size of the sampaling.

  res = 0.0

  print *, ''
  print *, 'result:'
  call monte_int(a, b, sample_size, res)
  write (*, form_result2r) res, '-- value of the integration'

end program int_mc
