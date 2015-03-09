program bond


   ! Use necessary modules.
   use O_CommandLine
   use O_TimeStamps
   use O_Bond
   use O_Bond3C
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
   call parseBondCommandLine(clp)


   ! Open the bond files that will be written to.
   open (unit=10,file='fort.10',status='new',form='formatted')
   if (spin == 2) then
      open (unit=11,file='fort.11',status='new',form='formatted')
   endif
   if (inDat%doBond3C == 1) then
      open (unit=12,file='fort.12',status='new',form='formatted')
      if (spin == 2) then
         open (unit=13,file='fort.13',status='new',form='formatted')
      endif
   endif

   ! Read in the input to initialize all the key data structure variables.
   call parseInput(inDat,clp)

   ! The bond calculation must be done without any smearing to determine
   !   the number of electrons for each atom.  So we set the thermal smearing
   !   factor to zero.
   inDat%thermalSigma = 0.0_double

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

   ! Call the bond subroutine to compute the bond order and effective charge.
   call computeBond(inDat,clp)

   ! Compute the three center bond order if requested.
   if (inDat%doBond3C == 1) then
      call computeBond3C(inDat,clp)
   endif

   ! Deallocate unnecesary matrices
   deallocate (energyEigenValues)

   ! Close access to the band HDF5 data.
   call closeAccessPSCFBandHDF5(clp)

   ! Close the BOND files that were written to.
   close (10)
   if (spin == 2) then
      close (11)
   endif
   if (inDat%doBond3C == 1) then
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
