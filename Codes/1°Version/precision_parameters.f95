!#######################################################################
!Parameters related to numerical precision
!#######################################################################
module precision_parameters 
 
   implicit none 
   integer, parameter :: i8b = selected_int_kind(20)                    !Larger integer
   integer, parameter :: i4b = selected_int_kind(9)                     !4-byte integer
   integer, parameter :: i2b = selected_int_kind(4)                     !2-byte integer
   integer, parameter :: i1b = selected_int_kind(2)                     !1-byte integer
   integer,parameter :: sp = kind(0.0)                                  !Single precision
   integer,parameter :: dp = kind(0.d0)                                 !Double precision
   real(dp), parameter :: small = 1.e-30_dp                             !Small number
   real(dp), parameter :: eps = epsilon(1._dp)
   real(dp), parameter :: big=huge(1._dp)*eps               

   !Eventually remove this
   real(dp) :: epsilon_left_array(100),epsilon_right_array(100),&
      lambda_left_array(99),lambda_right_array(99)
   integer :: number_epsilon_left_intervals,&
      number_epsilon_right_intervals
   
end module precision_parameters 
