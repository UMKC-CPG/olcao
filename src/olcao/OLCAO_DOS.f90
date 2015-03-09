program DOS


   ! Use necessary modules.
   use O_CommandLine
   use O_TimeStamps
   use O_DOS
   use O_Input
   use O_KPoints
   use O_Populate
   use O_SecularEquation
   use O_PSCFBandHDF5
   use O_Potential


   ! Make sure that no funny variables are defined.
   implicit none
   
   type(commandLineParameters) :: clp ! from O_CommandLine
   type(inputData) :: inDat ! from O_Input

   ! Initialize the logging labels.
   call initOperationLabels


   ! Parse the command line parameters
   call parseDOSCommandLine(clp)


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
   call parseInput(inDat,clp)


   ! Find specific computational parameters not EXPLICITLY given in the input
   !   file.  These values can, however, be easily determined from the input
   !   file.
   call getImplicitInfo


   ! Access the HDF5 data stored by band.
   call accessPSCFBandHDF5(inDat%numStates,clp)

   ! Allocate space to store the energy eigen values, and then read them in.
   allocate (energyEigenValues (inDat%numStates,numKPoints,spin))
   call readEnergyEigenValuesBand(inDat%numStates)

   ! Populate the electron states to find the highest occupied state (Fermi
   !   energ for metals).
   call populateStates(inDat,clp)

   ! Shift the energy eigen values according to the highest occupied state and
   !   convert them from au to eV.
   call shiftEnergyEigenValues(occupiedEnergy,inDat%numStates)
!   call convertEnergyEigenValuesToeV

   ! Call the DOS subroutine to compute the total and partial density of states
   !   as well as the localization index.
   call computeDOS(inDat,clp)

   ! Deallocate unnecesary matrices
   deallocate (energyEigenValues)

   ! Close access to the band HDF5 data.
   call closeAccessPSCFBandHDF5(clp)

   ! Open a file to signal completion of the program.
   open (unit=2,file='fort.2',status='new')

end program DOS



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
