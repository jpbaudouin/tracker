!#######################################################################
!################## DEFINITION OF THE CORE PARAMETERS ##################
!#######################################################################


!25-10-2017 Creation
!12-12-2017 1st Commit -> pi, dtr


! CONTAINS :
!   The modules core with different parameters




module core
  ! core subroutine and parameters

  implicit none
    
  real, parameter :: pi = 4. * atan(1.)
  real, parameter :: dtr = 4. * atan(1.)/180.0 !pi/180




end module




