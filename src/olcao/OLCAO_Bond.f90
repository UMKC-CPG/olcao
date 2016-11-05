program bond


   ! Use necessary modules.
   use O_TimeStamps
   use O_Potential,       only: spin
   use O_KPoints,         only: numKPoints
   use O_Bond,            only: computeBond
   use O_Bond3C,          only: computeBond3C
   use O_CommandLine,     only: parseBondCommandLine
   use O_Populate,        only: occupiedEnergy, populateStates
   use O_Input,           only: doBond3C, thermalSigma, numStates, parseInput
   use O_PSCFBandHDF5,    only: accessPSCFBandHDF5, closeAccessPSCFBandHDF5
   use O_SecularEquation, only: energyEigenValues, readEnergyEigenValuesBand, &
         & shiftEnergyEigenValues


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
   if (doBond3C == 1) then
      open (unit=12,file='fort.12',status='new',form='formatted')
      if (spin == 2) then
         open (unit=13,file='fort.13',status='new',form='formatted')
      endif
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
   call accessPSCFBandHDF5(numStates)

   ! Allocate space to store the energy eigen values, and then read them in.
   allocate (energyEigenValues (numStates,numKPoints,spin))
   call readEnergyEigenValuesBand(numStates)

   ! Populate the electron states to find the highest occupied state (Fermi
   !   energ for metals).
   call populateStates

   ! Shift the energy eigen values according to the highest occupied state.
   call shiftEnergyEigenValues(occupiedEnergy,numStates)

   ! Call the bond subroutine to compute the bond order and effective charge.
   call computeBond

   ! Compute the three center bond order if requested.
   if (doBond3C == 1) then
      call computeBond3C
   endif

   ! Deallocate unnecesary matrices
   deallocate (energyEigenValues)

   ! Close access to the band HDF5 data.
   call closeAccessPSCFBandHDF5

   ! Close the BOND files that were written to.
   close (10)
   if (spin == 2) then
      close (11)
   endif
   if (doBond3C == 1) then
      close (12)
      if (spin == 2) then
         close (13)
      endif
   endif

   ! Open a file to signal completion of the program.
   open (unit=2,file='fort.2',status='new')

end program bond



subroutine getImplicitInfo

   ! Import necessary modules.
   use O_TimeStamps
   use O_ExchangeCorrelation, only: makeSampleVectors
   use O_AtomicSites,         only: getAtomicSiteImplicitInfo
   use O_AtomicTypes,         only: getAtomicTypeImplicitInfo
   use O_PotSites,            only: getPotSiteImplicitInfo
   use O_PotTypes,            only: getPotTypeImplicitInfo
   use O_Lattice,             only: getRecipCellVectors
   use O_KPoints,             only: convertKPointsToXYZ
   use O_Potential,           only: initPotStructures

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
