module types
  ! ----------------------------------------------------------------
  ! cosmetics:
  !   This module contains all types and formats
  ! author:
  !   Exam-roll -- 20419PHY083
  ! env details:
  !   compiler -- gfortran v10.2.0
  !   OS       -- macOS 12.0 x86_64
  !   Kernel   -- Darwin 21.0.0
  ! Note:
  !   - used for the module "tools"
  ! ----------------------------------------------------------------
  implicit none

  ! various parameters used across my programs
  integer, parameter :: dp = selected_real_kind(15, 307) ! double precession
  real(dp), parameter :: pi = 4.d0*atan(1.d0) ! value of PI

  ! format specifier section
  ! the "----" lines
  character(len=15), parameter :: line2 = "(1x,2a10)"
  character(len=20), parameter :: line3 = "(a8,2(2x,a8))"
  character(len=25), parameter :: line3n = "(2(2x,a8),10x,a25)"
  character(len=25), parameter :: line5 = "(a8,3(2x,a8),5x,a8)"
  character(len=25), parameter :: line6 = "(a8,4(2x,a8),5x,a8)"
  ! the table headers
  character(len=15), parameter :: form_header2 = "(3x,a8,2x,a8)"
  character(len=20), parameter :: form_header3 = "(/,a5,2(5x,a5))"
  character(len=25), parameter :: form_header3n = "(/,3x,a5,4x,a5,20x,a8)"
  character(len=25), parameter :: form_header5 = "(/,a6,3(5x,a5),8x,a6)"
  character(len=25), parameter :: form_header6 = "(/,a6,4(5x,a5),8x,a6)"
  ! the table datas
  character(len=10), parameter :: form_table2 = "(2f10.4)"
  character(len=20), parameter :: form_table3 = "(i5,3x,2f10.4)"
  character(len=30), parameter :: form_table5 = "(i5,3x,3F10.4,2(3x,f10.5))"
  character(len=30), parameter :: form_table6 = "(i5,3x,4F10.4,3x,f10.5)"
  character(len=30), parameter :: form_tableip = "(f10.2, 10f10.4)"
  ! the final results
  character(len=15), parameter :: form_result1rj = "(f10.4,2x,a,i1)"
  character(len=15), parameter :: form_result2r = "(f10.4,2x,a)"
  character(len=15), parameter :: form_result2rp = "(f15.8,2x,a)"
  character(len=15), parameter :: form_result2i = "(i15,2x,a)"
  character(len=15), parameter :: form_result2c = "(3x,a10,a28)"

end module types
