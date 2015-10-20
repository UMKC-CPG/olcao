program DOS


   ! Use necessary modules.
   use O_TimeStamps
   use O_Potential,       only: spin
   use O_DOS,             only: computeDOS
   use O_KPoints,         only: numKPoints
   use O_CommandLine,     only: parseDOSCommandLine
   use O_Input,           only: numStates, parseInput
   use O_Populate,        only: occupiedEnergy, populateStates
   use O_PSCFBandHDF5,    only: accessPSCFBandHDF5, closeAccessPSCFBandHDF5
   use O_SecularEquation, only: energyEigenValues, readEnergyEigenValuesBand, &
         & shiftEnergyEigenValues


   ! Make sure that no funny variables are defined.
   implicit none

   ! Initialize the logging labels.
   call initOperationLabels


   ! Parse the command line parameters
   call parseDOSCommandLine


   ! Open the DOS files that will be written to.  If a spin polarized
   !   calculation is being done, then 60, 70, 80 hold spin up and 61, 71, 81
   !   hold spin down.  80,81=Localization Index
   open (unit=60,file='fort.60',status='new',form='formatted') ! TDOS
   open (unit=70,file='fort.70',status='new',form='formatted') ! PDOS
   open (unit=80,file='fort.80',status='new',form='formatted') ! LI
   if (spin == 2) then
      open (unit=61,file='fort.61',status='new',form='formatted') ! Spin TDOS
      open (unit=71,file='fort.71',status='new',form='formatted') ! Spin PDOS
      open (unit=81,file='fort.81',status='new',form='formatted') ! Spin LI
   endif


   ! Read in the input to initialize all the key data structure variables.
   call parseInput


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

   ! Shift the energy eigen values according to the highest occupied state and
   !   convert them from au to eV.
   call shiftEnergyEigenValues(occupiedEnergy,numStates)
!   call convertEnergyEigenValuesToeV

   ! Call the DOS subroutine to compute the total and partial density of states
   !   as well as the localization index.
   call computeDOS

   ! Deallocate unnecesary matrices
   deallocate (energyEigenValues)

   ! Close access to the band HDF5 data.
   call closeAccessPSCFBandHDF5

   ! Open a file to signal completion of the program.
   open (unit=2,file='fort.2',status='new')

end program DOS



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
