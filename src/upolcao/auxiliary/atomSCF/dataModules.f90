module ExecutionData

   ! Import the necessary modules.
   use O_Kinds

   integer :: maxIteration
   integer :: doSpinPol
   integer :: doRelativistic
   integer :: iteration    ! Range = 1..maxIteration
   integer :: converged    ! 0=not converged;   1=converged
   integer :: solverChoice ! 1=diagonalization; 2=direct integration
   integer :: exchCorrCode ! 0=xAlpha only; 1=wigner; 2=hedinLundqvist;
         ! 3=gunnarsonLundqvistWilkins; 4=vonBarthHedin; 5=ceperleyAlder;
         ! 6=vonBarthHedinGC

   real (kind=double) :: tolerance    ! Convergence tolerance
   real (kind=double) :: mixingFactor
   real (kind=double) :: oldMaxDeltaPot ! The maximum potential change from the
         ! previous iteration.

end module ExecutionData



module AtomData

   ! Import the necessary modules.
   use O_Kinds

   character*2 :: elementName
   integer :: numCoreOrb
   integer :: numValeOrb
   integer :: numOrbitals
   integer :: maxQNl  ! Maximum l quantum number.
   real (kind=double) :: nuclearAlpha
   real (kind=double) :: atomicNumber
   real (kind=double) :: shellCharge
   real (kind=double) :: shellRadius
   real (kind=double) :: coreCharge
   real (kind=double) :: valeCharge
   real (kind=double) :: electronCharge
   real (kind=double) :: ionicCharge
   integer, allocatable, dimension (:) :: orbitalQNn ! Dimension=numorbitals
   integer, allocatable, dimension (:) :: orbitalQNl ! Dimension=numOrbitals
   real (kind=double), allocatable, dimension (:) :: orbitalQNs !Dim=numOrbitals
   real (kind=double), allocatable, dimension (:) :: orbitalCharge ! As above

end module AtomData



module PotData

   ! Import the necessary modules.
   use O_Kinds

   ! Ionic potential as a function of r.
   real (kind=double), allocatable, dimension (:,:) :: ionicPotDn ! Old viod
   real (kind=double), allocatable, dimension (:,:) :: ionicPotUp ! Old viou

   ! Hartree potential as a function of r.
   real (kind=double), allocatable, dimension (:) :: hartreePotDn ! Old vod
   real (kind=double), allocatable, dimension (:) :: hartreePotUp ! Old vou

   ! Hartree potential as a function of r copied for some reason.
   real (kind=double), allocatable, dimension (:) :: hartreePotDn2 ! Old vid
   real (kind=double), allocatable, dimension (:) :: hartreePotUp2 ! Old viu

   ! Exchange-Correlation and Correlation potential as a function of r.
   real (kind=double), allocatable, dimension (:) :: vxcRadial
   real (kind=double), allocatable, dimension (:) ::  vcRadial

end module PotData



module EnergyData

   ! Import the necessary modules.
   use O_Kinds

   ! Energy eigen values.
   real (kind=double), allocatable, dimension (:) :: eigenValues

   ! Exchange-Correlation and Correlation energy as a function of r.
   real (kind=double), allocatable, dimension (:) :: excRadial
   real (kind=double), allocatable, dimension (:) ::  ecRadial

   ! Other contributions to total energy
   real (kind=double) :: hartreeE
   real (kind=double), allocatable, dimension(:) :: orbitalKE
   real (kind=double), allocatable, dimension(:) :: orbitalPE
   real (kind=double), dimension(10) :: totalEnergy
   
end module EnergyData



module WaveFnData

   ! Import the necessary modules.
   use O_Kinds

   ! Wave function solved from method 2 for individual l,s pairs.
   real (kind=double), allocatable, dimension (:) :: waveFn
   real (kind=double), allocatable, dimension (:) :: dr1WaveFn

! Not used.  Strange.
!   ! Wave function solved in method 2 accumulated for l.
!   real (kind=double), allocatable, dimension (:,:) :: waveFnQNl

end module WaveFnData



module ChargeDensityData

   ! Import the necessary modules.
   use O_Kinds

   ! Charge density arrays in various forms.
   real (kind=double), allocatable, dimension (:) :: rhoDn
   real (kind=double), allocatable, dimension (:) :: rhoUp
   real (kind=double), allocatable, dimension (:) :: rhoTotal
   real (kind=double), allocatable, dimension (:) :: dr1RhoTotal ! 1st deriv.
   real (kind=double), allocatable, dimension (:) :: dr2RhoTotal ! 2nd deriv.
   real (kind=double), allocatable, dimension (:) :: rhoOverR ! Total / r
   real (kind=double), allocatable, dimension (:) :: dr1RhoOverR ! 1st deriv.
   real (kind=double), allocatable, dimension (:) :: dr2RhoOverR ! 2nd deriv.

   ! Core charge density.
   real (kind=double), allocatable, dimension (:) :: rhoCore

   ! Charge density renormalization factor for each iteration.
   real (kind=double) :: renormFactor

   ! Cumulative integrated charge density.
   real (kind=double), allocatable, dimension (:) :: cumulInt1 ! rhoTotal
   real (kind=double), allocatable, dimension (:) :: cumulInt2 ! rhoTotal/r

end module ChargeDensityData



module TempOrbitalData

   ! Import necessary definition module
    use O_Kinds

   integer, allocatable, dimension (:) :: tempValeQNn
   integer, allocatable, dimension (:) :: tempValeQNl
   real (kind=double), allocatable, dimension (:) :: tempValeOrbChargeDn
   real (kind=double), allocatable, dimension (:) :: tempValeOrbChargeUp

end module TempOrbitalData
