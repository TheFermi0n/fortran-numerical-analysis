program jacobi
  use tools

  implicit none
  integer, parameter :: dim = 3 ! matirx dimension
  real(dp), parameter :: eps = 0.5e-4
  real(dp), dimension(dim, dim) :: a
  real(dp), dimension(dim, 2) :: x
  real(dp), dimension(dim) :: b

  real(dp) :: error
  integer :: i, itn

  ! the A,B matrix
  ! data(a(1, i), i=1, 4)/10.0, -2.0, -1.0, -1.0/
  ! data(a(2, i), i=1, 4)/-2.0, 10.0, -1.0, -1.0/
  ! data(a(3, i), i=1, 4)/-1.0, -1.0, 10.0, -2.0/
  ! data(a(4, i), i=1, 4)/-1.0, -1.0, -2.0, 10.0/
  ! data(b(i), i=1, 4)/3.0, 15.0, 27.0, -9.0/

  ! data(a(1, i), i=1, 3)/9.0, 2.0, 4.0/
  ! data(a(2, i), i=1, 3)/1.0, 10.0, 4.0/
  ! data(a(3, i), i=1, 3)/2.0, -4.0, 10.0/
  ! data(b(i), i=1, 3)/20, 6, -15/

  ! data(a(1, i), i=1, 3)/4.0, -1.0, -1.0/
  ! data(a(2, i), i=1, 3)/-2.0, 6.0, 1.0/
  ! data(a(3, i), i=1, 3)/-1.0, 1.0, 7.0/
  ! data(b(i), i=1, 3)/3.0, 9.0, -6.0/

  data(a(1, i), i=1, 3)/6.0, 1.0, 1.0/
  data(a(2, i), i=1, 3)/1.0, 4.0, -1.0/
  data(a(3, i), i=1, 3)/1.0, -1.0, 5.0/
  data(b(i), i=1, 3)/20.0, 6.0, 7.0/

  call dominance(dim, a) ! to check if it is strictly dominant

  print *, ''
  print *, 'table:'
  write (*, form_header6) 'itn', 'x1', 'x2', 'x3', 'x4', 'error'
  write (*, line6) '-------', '-------', '-------', '-------', '-------', '---------'

  x(:, 1) = 1.0
  itn = 0

  call jac(dim, a, x, b, error, itn)

  do while (error >= eps)
    x(:, 1) = x(:, 2)
    itn = itn + 1
    call jac(dim, a, x, b, error, itn)
  end do

  print *, ''
  print *, 'result:'
  write (*, form_result1rj) (x(i, 2), '-- value of x', i, i=1, dim)
end program jacobi
