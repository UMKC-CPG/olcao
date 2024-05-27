!#include "../config.h"
!****************************************************************************
!
! Physics and control constants used throughout all the programs
!
!****************************************************************************
module O_Constants

   ! Import the precision variables
   use O_Kinds

   ! Make sure that there are not accidental variable declarations.
   implicit none

   !Start defining the constants
   real (kind=double) :: pi
   real (kind=double) :: hartree
   real (kind=double) :: hPlanck
   real (kind=double) :: hBarPlanck
   real (kind=double) :: boltzConst
   real (kind=double) :: eCharge
   real (kind=double) :: eMass
   real (kind=double) :: bohrRad
   real (kind=double) :: auTime
   real (kind=double) :: lightFactor
   real (kind=double) :: lightSpeed
   real (kind=double) :: fineStructure
   real (kind=double) :: smallThresh, bigThresh
   integer :: dim3
   integer :: lAngMomCount
   integer :: num_lAngMomTerms

   ! All values are from NIST unless otherwise noted.  The exponentials are not
   !   included here and must be accounted for in the particular calculation
   !   that is being done in the program.
   parameter (pi       = 3.141592653589793_double) ! Pi not from NIST.
   parameter (hartree  = 27.211386245988_double)   ! Hartree from NIST. (eV)
   parameter (hPlanck  = 6.62607015_double)    ! Planck's constant (Js) (e-34)
   parameter (hBarPlanck = 1.054571817_double) ! Reduced (Js) (e-34)
   parameter (boltzConst = 8.617333262_double) ! Boltzmann's const (eV/K) (e-5)
   parameter (eCharge  = 1.602176634_double)   ! Elementary Charge (C) (e-19)
   parameter (eMass    = 9.1093837015_double)  ! Electron mass (kg) (e-31)
   parameter (bohrRad  = 0.529177210903_double)! Bohr Radius (m) (e-10)
   parameter (auTime   = 2.4188843265857_double)! Atomic unit time (s) (e-17)
   parameter (lightFactor = 8.9875517873681764_double) ! (1e-16 * c^2) (m/s)^2
   parameter (lightSpeed = 299792458.0_double) ! m/s
   parameter (fineStructure = 7.2973525693_double) ! Fine structure (e-3).
         !   NOTE: This is also equal to (1/lightSpeed = 1/c) in atomic units.
   parameter (smallThresh=1.0E-8_double) ! Threshold value to determine if some
         !   number is big enough or negligable.
   parameter (bigThresh=1.0E20_double) ! Threshold value to initialize a search
         !   for some number.  The answer must clearly be less than this value.
   parameter (dim3=3) ! The three real dimension.
   parameter (lAngMomCount=4) ! Allow program to work for S, P, D, and F
         !   orbitals. Note that this is a count of the *number* of different
         !   l angular momentum kinds, not the actual value of angular
         !   momentum. So, if lAngMomCount = 4 then the maximum l angular
         !   momentum value is actually 3. s=0 -> lAngMomCount=1, etc.
   parameter (num_lAngMomTerms=(2*(lAngMomCount-1)+2) &
         & * ((lAngMomCount-1)+1) / 2) ! Total number of different angular
         !   momentum terms 1s+3p+5d+7f+...

end module O_Constants
