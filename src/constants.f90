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
   real (kind=double) :: hPlank
   real (kind=double) :: boltzConst
   real (kind=double) :: eCharge
   real (kind=double) :: eMass
   real (kind=double) :: bohrRad
   real (kind=double) :: auTime
   real (kind=double) :: lightFactor
   real (kind=double) :: fineStructure
   real (kind=double) :: smallThresh, bigThresh
   integer :: dim3
   integer :: maxOrbitals

   ! All values are from NIST unless otherwise noted.  The exponentials are not
   !   included here and must be accounted for in the particular calculation
   !   that is being done in the program.
   parameter (pi       = 3.141592653589793_double) ! Pi not from NIST.
   parameter (hartree  = 27.21138386_double)   ! Hartree from NIST.
   parameter (hPlank   = 6.6260693_double)     ! Plank's constant (Js) (e-34)
   parameter (boltzConst = 8.617343_double)    ! Boltzmann's const (eV/K) (e-5)
   parameter (eCharge  = 1.60217653_double)    ! Elementary Charge (C) (e-19)
   parameter (eMass    = 9.1093826_double)     ! Electron mass (kg) (e-31)
   parameter (bohrRad  = 0.5291772180_double)  ! Bohr Radius (m) (e-10)
   parameter (auTime   = 2.418884326505_double)! Atomic unit of time (s) (e-17)
   parameter (lightFactor = 8.9875517873681764_double) ! 1/(1e-16 * c^2)(s/m)^2
   parameter (fineStructure = 7.2973525376_double) ! Fine structure (e-3).
   parameter (smallThresh=1.0E-8_double) ! Threshold value to determine if some
         !   number is big enough or negligable.
   parameter (bigThresh=1.0E20_double) ! Threshold value to initialize a search
         !   for some number.  The answer must clearly be less than this value.
   parameter (dim3=3) ! The three real dimension.
   parameter (maxOrbitals=4) ! Allow program to work for S, P, D, and F orbitals

end module O_Constants
