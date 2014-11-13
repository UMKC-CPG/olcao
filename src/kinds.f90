!#include "../config.h"
!****************************************************************************
!
! Machine independent declaration of number variable precision
!
!****************************************************************************
module O_Kinds

   ! Make sure that there are not accidental variable declarations.
   implicit none

!   integer, parameter :: single  = selected_real_kind(6,37)
!   integer, parameter :: double  = selected_real_kind(15,307)
   integer, parameter :: single  = kind(0)
   integer, parameter :: double  = kind(0.0d0)

end module O_Kinds
