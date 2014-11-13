program bond


   ! Use necessary modules.
   use O_CommandLine
   use O_TimeStamps
   use O_Bond
   use O_Input
   use O_KPoints
   use O_Populate
   use O_SecularEquation
   use O_PSCFBandHDF5
   use O_Potential


   ! Make sure that no funny variables are defined.
   implicit none


   ! Initialize the logging labels.
   call initOperationLabels


   ! Parse the command line parameters
   call parseBondCommandLine


   ! Open the bond files that will be written to.
   open (unit=10,file='fort.10',status='new',form='formatted')
   if (spin == 2) then
      open (unit=11,file='fort.11',status='new',form='formatted')
   endif


   ! Read in the input to initialize all the key data structure variables.
   call parseInput

   ! The bond calculation must be done without any smearing to determine
   !   the number of electrons for each atom.  So we set the thermal smearing
   !   factor to zero.
   thermalSigma = 0.0_double

   ! Find specific computational parameters not EXPLICITLY given in the input
   !   file.  These values can, however, be easily determined from the input
   !   file.
   call getImplicitInfo


   ! Access the HDF5 data stored by band.
   call accessPSCFBandHDF5

   ! Allocate space to store the energy eigen values, and then read them in.
   allocate (energyEigenValues (numStates,numKPoints,spin))
   call readEnergyEigenValuesBand

   ! Populate the electron states to find the highest occupied state (Fermi
   !   energ for metals).
   call populateStates

   ! Shift the energy eigen values according to the highest occupied state and
   !   convert them from au to eV.
   call shiftEnergyEigenValues(occupiedEnergy)
!   call convertEnergyEigenValuesToeV

   ! Call the bond subroutine to compute the bond order and effective charge.
   call computeBond

   ! Deallocate unnecesary matrices
   deallocate (energyEigenValues)

   ! Close access to the band HDF5 data.
   call closeAccessPSCFBandHDF5

   ! Open a file to signal completion of the program.
   open (unit=2,file='fort.2',status='new')

end program bond



subroutine getImplicitInfo

   ! Import necessary modules.
   use O_ExchangeCorrelation
   use O_AtomicSites
   use O_AtomicTypes
   use O_PotSites
   use O_PotTypes
   use O_Lattice
   use O_KPoints
   use O_Potential
   use O_Populate
   use O_TimeStamps

   implicit none

   call timeStampStart(2)

   ! Subroutines need to be called in this order due to data dependencies.
   call makeSampleVectors

   call getAtomicTypeImplicitInfo
   call getAtomicSiteImplicitInfo
   call getPotSiteImplicitInfo
   call getPotTypeImplicitInfo

   call getRecipCellVectors
   call convertKPointsToXYZ

   call initPotStructures

   call timeStampEnd(2)

end subroutine getImplicitInfo
