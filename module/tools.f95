module tools
  use types
  ! ----------------------------------------------------------------
  ! tools:
  !   This module contains all the functions and subroutines
  ! author:
  !   Exam-roll -- 20419PHY083
  ! env details:
  !   compiler -- gfortran v10.2.0
  !   OS       -- macOS 12.0 x86_64
  !   Kernel   -- Darwin 21.0.0
  ! Note:
  !   - depends on the module "types"
  ! ----------------------------------------------------------------
  implicit none

contains
  function fn(x) result(z)
    ! ----------------------------------------------------------------
    ! Single variable function:
    !   function with one variable that needs to be calculated
    ! input:
    !   x -- input of the function
    ! output:
    !   z -- output of the function f(x)
    ! ----------------------------------------------------------------
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: z

    z = 2.0*x**3.0 - 3.0*x**2.0 + 4.0*x - 5.0
    ! z = 0.5*exp(-((1 + x)/2.0)**2.0)
    ! z = 1/(x + 2)
  end function fn

  function f(x, y) result(z)
    ! ----------------------------------------------------------------
    ! Two variable function:
    !   function with two variables that needs to be calculated
    ! inputs:
    !   x,y -- inputs of the function
    ! output:
    !   z   -- output of the function f(x,y)
    ! ----------------------------------------------------------------
    implicit none
    real(dp), intent(in) :: x, y
    real(dp) :: z

    ! z = exp(x) + y
    z = -x*y
    ! z = sqrt(1 - x**2.0)
  end function f

  subroutine oela(x0, y0, x_point, h)
    ! ----------------------------------------------------------------
    ! Euler method:
    !   Finds the value of y at a given x by using a picewise linear
    !   approximation of the solution
    ! inputs:
    !   h         -- interval size (keep it small for better result)
    !   x0        -- initial value of x (also used as output)
    !   y0        -- initial value of y (also used as output)
    !   x_point   -- the value of x at which y (f(x)) is to be calculated
    ! internals:
    !   x_new     -- new value of x's calculated
    !   y_new     -- new value of y's calculated
    !   k1        -- slope
    ! outputs:
    !   y0, x0    -- values of the required y at given x
    ! ----------------------------------------------------------------
    implicit none
    real(dp), intent(in) :: x_point, h
    real(dp) :: x0, y0
    real(dp) :: k1, x_new, y_new
    integer :: itn

    x_new = x0
    y_new = y0
    itn = 0

    do while (x_new <= x_point)
      ! displaying the values
      write (*, form_table3) itn, x_new, y_new
      x0 = x_new
      y0 = y_new

      ! valus of slope
      k1 = f(x_new, y_new)
      ! calculating the y
      y_new = y_new + h*k1

      ! preparing for the next value
      x_new = x_new + h
      itn = itn + 1
    end do
  end subroutine oela

  subroutine mod_oela(x0, y0, x_point, h)
    ! ----------------------------------------------------------------
    ! Modified Euler method:
    !   Finds the value of y at a given x also known as Heun's method
    ! inputs:
    !   h         -- interval size (keep it small for better result)
    !   x0        -- initial value of x (also used as output)
    !   y0        -- initial value of y (also used as output)
    !   x_point   -- the value of x at which y (f(x)) is to be calculated
    ! internals:
    !   x_new     -- new value of x's calculated
    !   y_new     -- new value of y's calculated
    !   k1, k2    -- slope
    ! outputs:
    !   y0, x0    -- values of the required y at given x
    ! ----------------------------------------------------------------
    implicit none
    real(dp), intent(in) :: h, x_point
    real(dp) :: x0, y0
    real(dp) :: k1, k2, k, x_new, y_new
    integer :: itn

    x_new = x0
    y_new = y0
    itn = 0

    do while (x_new <= x_point)
      ! displaying the values
      write (*, form_table3) itn, x_new, y_new
      x0 = x_new
      y0 = y_new

      ! values of slope
      k1 = f(x_new, y_new)
      k2 = f(x_new + h, y_new + h*k1)

      ! weighted average of the slopes
      k = (k1 + k2)/2.0

      ! calculating the y
      y_new = y_new + h*k

      ! preparing for the next value
      x_new = x_new + h
      itn = itn + 1
    end do
  end subroutine mod_oela

  subroutine g_quad(k, a, b, int_val)
    ! ----------------------------------------------------------------
    ! Gaussian quadrature method:
    !   Integrates a given function by fitting a polynomial at unequal
    !   intervals whose abscissa and weights are calculated using the
    !   zeros of the legendre polynomial and Newton's method.
    ! inputs:
    !   k         -- number of points needed
    !   x()       -- the arrays of x's
    !   y()       -- the arrays of corresponding y's (f(x)'s)
    !   x_point   -- the value of x at which f'(x) is to be calculated
    ! internal:
    !   h         -- class height
    !   loc       -- determines the position of x_point on the table
    !                and gets the respective index for calculation
    !   y         -- value of the integral before scaling
    ! output:
    !   int_val   -- the result of the derivation after scaling
    ! Note:
    !   - used function fn() because only 1 variable
    !   - used subroutine absc_wght() for abscissa and weights
    ! ----------------------------------------------------------------
    implicit none
    integer, intent(in) :: k
    integer, intent(in) :: a, b
    real(dp), intent(out) :: int_val

    real(dp), dimension(k) :: x, w
    real(dp) :: fun, y, p, q
    integer :: i

    call absc_wght(k, x, w) ! calculating the weighting factors

    ! printing the above weighting factors
    write (*, form_header2) 'abscissa', 'weights'
    write (*, line2) '--------', '--------'
    do i = 1, k
      write (*, form_table2) x(i), w(i)
    end do

    ! main integration
    y = 0.0
    p = (a + b)/2.0; q = (b - a)/2.0

    ! loop to calculate the values for each weighted factor
    do i = 1, k
      fun = fn(p + q*x(i))
      y = y + w(i)*fun
    end do

    int_val = q*y ! final result
  end subroutine g_quad

  subroutine absc_wght(points, abscissa, weights)
    ! ----------------------------------------------------------------
    ! Gaussian quadrature Weighting Factors:
    !   calculate the abscissa and weights  for the gaussian quadrature
    !   method {g_quad()} using the zeros of the legendre polynomial
    ! inputs:
    !   points    -- number of data points
    !   abscissa  -- the abscissa
    !   weights   -- the weights
    ! internals:
    !   p0        -- Pn coefficient of the legendre polynomial
    !   p1        -- Pn-1 coefficient of the legendre polynomial
    !   p2        -- Pn-2 coefficient of the legendre polynomial
    !   m         -- finds the mid-point
    !   x         -- initial guess of the abscissa which is further refined
    !               using the Newton-rhapson method
    !   p0_dash   -- the derivative of the legendre polynomial
    ! Note:
    !   - Newton's method is used to refine the value of abscissa (x)
    !   - since the weighted factors are symmetric we can break it into half
    !     and calculate only one half and then set the other half
    ! ----------------------------------------------------------------
    implicit none
    integer, intent(in) :: points
    real(dp), dimension(points), intent(out) :: abscissa, weights

    real(dp), parameter :: eps = 1.0e-15
    integer :: i, j, m
    real(dp) :: p0, p1, p2, p0_dash, x, x_old

    if (mod(points, 2) == 0) then ! finding the mid-point
      m = points/2
    else
      m = (points + 1)/2
    end if

    ! loop to calculate the zeros of legendre polynomial
    do i = 1, m
      x = cos(pi*(i - 0.25)/(points + 0.5)) ! initial guess of x

      do while (abs(x - x_old) > eps)
        p0 = 1.0
        p1 = x

        ! loop to calculate the legendre polynomial
        do j = 1, points
          p2 = p1
          p1 = p0
          p0 = ((2.0*j - 1.0)*x*p1 - (j - 1.0)*p2)/j ! legendre polynomial
        end do

        p0_dash = points*(x*p0 - p1)/(x**2.0 - 1.0) ! legendre polynomial derivative

        x_old = x
        x = x_old - p0/p0_dash ! newton's method for refining x

        ! arranging the abscissa
        abscissa(i) = -x ! the -abscissa's
        abscissa(points + 1 - i) = x ! the +abscissa's

        ! arranging the weights
        weights(i) = (2.0*(1 - x**2.0))/(points*p1)**2.0 ! the weights for -abscissa
        weights(points + 1 - i) = weights(i) ! the weights for +abscissa
      end do
    end do
  end subroutine absc_wght

  subroutine rungkut_2(x0, y0, x_point, h)
    ! ----------------------------------------------------------------
    ! Runge-Kutte two point method:
    !   Finds the value of y at a given x also known as Polygon method
    ! inputs:
    !   h         -- interval size (keep it small for better result)
    !   x0        -- initial value of x (also used as output)
    !   y0        -- initial value of y (also used as output)
    !   x_point   -- the value of x at which y (f(x)) is to be calculated
    ! internals:
    !   x_new     -- new value of x's calculated
    !   y_new     -- new value of y's calculated
    !   k1, k2    -- slopes
    ! outputs:
    !   y0, x0    -- values of the required y at given x
    ! ----------------------------------------------------------------
    implicit none
    real(dp), intent(in) :: h, x_point
    real(dp) :: x0, y0
    real(dp) :: k1, k2, x_new, y_new
    integer :: itn

    x_new = x0
    y_new = y0
    itn = 0

    do while (x_new <= x_point)
      ! displaying the values
      write (*, form_table3) itn, x_new, y_new
      x0 = x_new
      y0 = y_new
      ! value of slopes
      k1 = f(x_new, y_new)
      k2 = f(x_new + h/2.0, y_new + k1*h/2.0)
      ! calculating the y
      y_new = y_new + h*k2
      ! preparing for next values
      x_new = x_new + h
      itn = itn + 1
    end do
  end subroutine rungkut_2

  subroutine rungkut_4(x0, y0, x_point, h)
    ! ----------------------------------------------------------------
    ! Runge-Kutte two point method:
    !   Finds the value of y at a given x also known as Polygon method
    ! inputs:
    !   h         -- interval size (keep it small for better result)
    !   x0        -- initial value of x (also used as output)
    !   y0        -- initial value of y (also used as output)
    !   x_point   -- the value of x at which y (f(x)) is to be calculated
    ! internals:
    !   x_new     -- new value of x's calculated
    !   y_new     -- new value of y's calculated
    !   k1,k2,k3,k4  -- slopes
    ! outputs:
    !   y0, x0    -- values of the required y at given x
    ! ----------------------------------------------------------------
    implicit none
    real(dp), intent(in) :: x_point, h
    real(dp) :: x0, y0
    real(dp) :: k1, k2, k3, k4, k, x_new, y_new
    integer :: itn

    x_new = x0
    y_new = y0
    itn = 0

    do while (x_new <= x_point)
      ! displaying the values
      write (*, form_table3) itn, x_new, y_new
      x0 = x_new
      y0 = y_new

      ! values of slopes
      k1 = f(x_new, y_new)
      k2 = f(x_new + h/2.0, y_new + k1*h/2.0)
      k3 = f(x_new + h/2.0, y_new + k2*h/2.0)
      k4 = f(x_new + h, y_new + k3*h)

      ! weighted average of the slopes
      k = (k1 + 2.0*(k2 + k3) + k4)/6.0

      ! calculating the y
      y_new = y_new + k*h

      ! preparingg for the next value
      x_new = x_new + h
      itn = itn + 1
    end do
  end subroutine rungkut_4

  subroutine neum_diff(n, x, y, x_point, poly, deriv)
    ! ----------------------------------------------------------------
    ! Newton's differentiation:
    !   Calculates the differentiation of y at a given x using newton's
    !   interpolation polynomial and the difference table calculated
    !   by the subroutine newinter()
    ! inputs:
    !   n         -- number of points for interpolation
    !   x()       -- the arrays of x's
    !   y()       -- the arrays of corresponding y's (f(x)'s)
    !   x_point   -- the value of x at which f'(x) is to be calculated
    ! internal:
    !   m         -- the mid-point of the table of n points
    !   h         -- class height
    !   loc       -- determines the position of x_point on the table
    !               and gets the respective index for calculation
    ! output:
    !   deriv     -- the result of the derivation
    !   poly      -- the used interpolation polynomial
    ! Note:
    !   - used the subroutine newinter
    !   - used forward difference table and modified the formula a bit
    !     to calculate the backward differentiation.
    ! ----------------------------------------------------------------
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(n), intent(in) :: x, y
    real(dp), intent(in) :: x_point
    real(dp), intent(out) :: deriv
    character(len=8), intent(out) :: poly

    real(dp), dimension(n, n) :: del_f
    real(dp) :: h, sum_1, sum_2, x_sum
    integer :: a, b, m, loc

    ! condition to get the mid-point of the data points
    if (mod(n, 2) == 0) then
      m = n/2
    else
      m = (n + 1)/2
    end if

    loc = minloc(x, mask=(x >= x_point), dim=1) ! position of an x in the array

    ! check if x_point is outside the given data points
    if ((loc == 0) .or. (loc == n + 1)) then
      stop 'your {x_point} is outside the given data points'
    end if

    call newinter(n, x, y, del_f) ! calculating the difference table

    ! conditon for forward and backward differentiation
    if (x_point <= x(m)) then ! forward case
      poly = 'forward'
      if (loc == 1) then
        a = loc
        b = a + 1
      else
        a = loc
        b = a - 1
      end if

      ! calculation of the derivative
      h = abs(x(a) - x(b))
      x_sum = ((x_point - x(b)) + (x_point - x(a)))
      sum_1 = del_f(b, 1)/h
      sum_2 = del_f(b, 2)/(2.0*h*h)

    else if (x_point > x(m)) then ! backward case
      poly = 'backward'
      if (loc == n) then
        a = loc
        b = a - 1
      else
        a = loc
        b = a + 1
      end if

      ! calculation of the derivative
      h = abs(x(a) - x(b))
      x_sum = ((x_point - x(b)) + (x_point - x(a)))
      sum_1 = del_f(b, 1)/h
      sum_2 = del_f(b - 1, 2)/(2.0*h*h)
    end if

    deriv = sum_1 + sum_2*x_sum ! final result
  end subroutine neum_diff

  subroutine newinter(n, x, y, del_f)
    ! ----------------------------------------------------------------
    ! Newton's interpolation table:
    !   Calculates the forward differences and gives us the F.D. table
    ! inputs:
    !   n         -- number of points for interpolation
    !   x()       -- the arrays of x's (used just for tabulation)
    !   y()       -- the arrays of corresponding y's (f(x)'s)
    ! output:
    !   del_f(,)  -- the forward difference coefficients
    ! Note:
    !   - This same table is used for backward difference because it is
    !     same as forward difference, just the position is different.
    ! ----------------------------------------------------------------
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(n), intent(in) :: y, x
    real(dp), dimension(n, n), intent(out) :: del_f

    real(dp), dimension(n, n) :: del_y
    integer :: i, j, p

    ! loop to calculate the differences
    do j = 1, n
      do i = 1, n - j
        if (j == 1) then ! calculating the 1st difference column
          del_y(i, 1) = y(i + 1) - y(i)
        else ! calculate the other coloumns
          del_y(i, j) = del_y(i + 1, j - 1) - del_y(i, j - 1)
        end if
      end do
    end do

    del_f = del_y

    ! printing the difference table
    print *, 'the difference table:'
    write (*, form_header3n) 'x', 'y', "del_f's"
    write (*, line3n) '------', '------', '-------------------------'
    do i = 1, n
      write (*, form_tableip) x(i), y(i), (del_f(i, p), p=1, n - i)
    end do
  end subroutine newinter

  subroutine jac(n, a, x, b, err, itn)
    ! ----------------------------------------------------------------
    ! Jacobi method:
    !   calculate the single value of x,y and z using the iterative
    !   method of simultaneous displacememnts
    ! inputs:
    !   n         -- number of data points (dimension of the matrix)
    !   a         -- matrix of the coefficients A
    !   b         -- matrix of the constants B
    !   x(:,1)    -- initial x1,x2,x3
    !   itn       -- iteration number
    ! internals:
    !   sum_ax()  -- sum of the off diagonal terms of A
    !   aberror   -- difference of the values of old x's with new x's
    ! output:
    !   err       -- the eucledian norm of "error"
    !   x(:,2)    -- final x1,x2,x3
    ! Note:
    !   - use the subroutine dominance() to check if the solution will
    !     converge or not
    ! ----------------------------------------------------------------
    implicit none
    integer, intent(in) :: n, itn
    real(dp), dimension(n), intent(in) :: b
    real(dp), dimension(n, n), intent(in) :: a
    real(dp), intent(out) :: err
    real(dp), dimension(n, 2) :: x

    real(dp), dimension(n) :: sum_ax, aberror
    integer :: p, q

    !loop to calculate x1,x2,x3
    do p = 1, n
      sum_ax(p) = 0.0
      do q = 1, n
        if (q /= p) then
          sum_ax(p) = sum_ax(p) + a(p, q)*x(q, 1)
        end if
      end do
      x(p, 2) = (b(p) - sum_ax(p))/a(p, p)
    end do

    ! error calculation
    aberror = x(:, 1) - x(:, 2)
    err = norm2(aberror) ! taking the eucledian norm

    ! printing the values
    write (*, form_table6) itn, (x(p, 2), p=1, n), err
  end subroutine jac

  subroutine dominance(n, a)
    ! ----------------------------------------------------------------
    ! Diagonally dominance checker:
    !   checks if the matrix is strictly diagonally dominant
    ! inputs:
    !   n         -- dimension of the matrix
    !   a         -- matrix of the coefficients A
    ! internals:
    !   row_sum   -- sum of the off diagonal row elements
    ! Note:
    !   - used for the jacobi method equation solver.
    !   - if dominance test failed, need to make the matrix dominant
    ! ----------------------------------------------------------------
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(n, n), intent(in) :: a

    real(dp) :: row_sum
    integer :: p, q

    ! loop to calculate the sum of off diagonal terms
    do p = 1, n
      row_sum = 0.0
      do q = 1, n
        if (q /= p) then
          row_sum = row_sum + abs(a(p, q)) ! off-diag element sum
        end if
      end do

      ! condition for strictly diagonally dominant
      if (abs(a(p, p)) < row_sum) then
        stop "-- !! the matirx is NOT strictly diagonally dominant !! --"
      end if
    end do
  end subroutine dominance

  subroutine monte_pi(sq_points, points_inside, approx_pi)
    ! ----------------------------------------------------------------
    ! Monte Carlo Pi:
    !   Calculates the approximate value of pi using Monte Carlo method
    ! inputs:
    !   sq_points     -- total points inside square including the circle
    !   points_inside -- total points inside the circle
    ! internals:
    !   rand_x, rand_y -- uniform random number in the range [0,1)
    ! ----------------------------------------------------------------
    implicit none
    integer, intent(in) :: sq_points
    integer, intent(out) :: points_inside
    real(dp), intent(out) :: approx_pi

    integer :: i, cr_points
    real(dp) :: x, y, rand_x, rand_y, cr_radi

    call random_seed() ! to maintain the state random_number generator

    do i = 1, sq_points
      ! generating uniform random numbers in the range [0,1)
      call random_number(rand_x)
      call random_number(rand_y)

      ! scaling in the range [a,b]: x = a + (b-a)*x
      x = -1 + 2*rand_x
      y = -1 + 2*rand_y

      cr_radi = x**2.0 + y**2.0 ! radius of the circle

      ! condition to check if the hits were perfect
      if (cr_radi <= 1.0) then
        cr_points = cr_points + 1
      end if
    end do

    ! final result
    points_inside = cr_points
    approx_pi = 4.0*(float(cr_points)/float(sq_points))
  end subroutine monte_pi

  subroutine monte_int(a, b, s_size, approx_int)
    ! ----------------------------------------------------------------
    ! Monte Carlo Integration:
    !   Calculates the approximate value of integral using MC method
    ! inputs:
    !   a,b           -- lower and upper limit
    !   s_size        -- size of the sample
    ! internals:
    !   rand_x, rand_y -- uniform random number in the range [0,1)
    !   sum_f         -- sum of the function from 1 to sample size
    ! output:
    !   approx_int    -- the approximate value of the integral
    ! ----------------------------------------------------------------
    implicit none
    integer, intent(in) :: a, b, s_size
    real(dp), intent(out) :: approx_int

    real(dp) :: x, y, sum_f, rand_x, rand_y
    integer :: i

    call random_seed() ! to maintain the state random_number generator

    ! loop to calculate the integral
    do i = 1, s_size
      ! generating uniform random numbers in the range [0,1)
      call random_number(rand_x)
      call random_number(rand_y)

      ! scaling in the range [a,b]: x = a + (b-a)*x
      x = a + rand_x*(b - a)
      y = a + rand_y*(b - a)

      sum_f = sum_f + f(x, y) ! summing all the values of the given function
    end do

    approx_int = float(b - a)*(sum_f/float(s_size)) ! final result of integral
  end subroutine monte_int

  subroutine garphmonte_pi(sq_points, points_inside, approx_pi)
    ! ----------------------------------------------------------------
    ! Monte Carlo Pi:
    !   Calculates the approximate value of pi using Monte Carlo method
    ! inputs:
    !   sq_points     -- total points inside square including the circle
    !   points_inside -- total points inside the circle
    ! internals:
    !   rand_x, rand_y -- uniform random number in the range [0,1)
    ! Note:
    !   - this subroutine can be used to get the points in a output file
    !     and can be plotted using a graphing tool
    ! ----------------------------------------------------------------
    implicit none
    integer, intent(in) :: sq_points
    integer, intent(out) :: points_inside
    real(dp), intent(out) :: approx_pi

    integer :: i, cr_points
    real(dp) :: cr_radi
    real(dp) :: x, y, rand_x, rand_y

    call random_seed() ! to maintain the state random_number generator

    open (1, file='values_inside.dat')
    open (2, file='values_outside.dat')

    do i = 1, sq_points
      ! generating uniform random numbers in the range [0,1)
      call random_number(rand_x)
      call random_number(rand_y)

      ! scaling in the range [a,b]: x = a + (b-a)*x
      x = -1 + 2*rand_x
      y = -1 + 2*rand_y

      cr_radi = x**2.0 + y**2.0 ! radius of the circle

      ! condition to check if the hits were perfect
      if (cr_radi <= 1.0) then
        cr_points = cr_points + 1
        write (1, "(2f15.5)") x, y
      else
        write (2, "(2f15.5)") x, y
      end if
    end do

    close (1); close (2)

    ! final result
    points_inside = cr_points
    approx_pi = 4.0*(float(cr_points)/float(sq_points))
  end subroutine garphmonte_pi

end module tools
