module MixPotentialSubs

public

contains

! Compute the new exchange correlation potential given the current guess and
!   the potential from the previous iteration.
subroutine mixPotential

   ! Import necessary subroutines.
   use O_Kinds
   use PotData
   use EnergyData
   use AtomData
   use O_RadialGrid
   use ExecutionData

   implicit none

   ! Define local variables.
   real (kind=double) :: mixingComp
   real (kind=double), allocatable, dimension(:) :: newPotDn
   real (kind=double), allocatable, dimension(:) :: newPotUp

   ! Allocate space to hold the temporary new potential.
   allocate (newPotDn(numRadialPoints))
   allocate (newPotUp(numRadialPoints))

   ! Compute the compliment of the mixing factor.
   mixingComp = 1.0_double - mixingFactor

   ! Compute the new potentials as a function of the current guess and the
   !   old one from the previous iteration.
   newPotDn(:) = mixingFactor * hartreePotDn(:) + mixingComp * hartreePotDn2(:)
   newPotUp(:) = mixingFactor * hartreePotUp(:) + mixingComp * hartreePotUp2(:)

   ! Copy the new potentials over the current guess.
   hartreePotDn2(:) = newPotDn(:)
   hartreePotUp2(:) = newPotUp(:)

   ! Deallocate the space for temporarily holding the new potentials.
   deallocate (newPotDn)
   deallocate (newPotUp)

end subroutine mixPotential

end module MixPotentialSubs
