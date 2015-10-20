module O_Optc

   contains

subroutine computeOptc


   ! Use necessary modules.
   use O_TimeStamps
   use O_Potential,       only: spin
   use O_PSCFIntgHDF5,    only: getMOMEStatus
   use O_OptcPrint,       only: printOptcResults
   use O_CommandLine,     only: parseOptcCommandLine, stateSet
   use O_KPoints,         only: numKPoints, computePhaseFactors
   use O_Populate,        only: occupiedEnergy, populateStates
   use O_Lattice,         only: initializeLattice, initializeFindVec
   use O_Input,           only: numStates, lastInitStatePACS, parseInput
   use O_PSCFBandHDF5,    only: accessPSCFBandHDF5, closeAccessPSCFBandHDF5
   use O_OptcTransitions, only: transCounter, energyDiff, transitionProb, &
         & getEnergyStatistics, computeTransitions
   use O_SecularEquation, only: energyEigenValues, readEnergyEigenValuesBand, &
         & appendExcitedEnergyEigenValuesBand, shiftEnergyEigenValues


   ! Make sure that no funny variables are defined.
   implicit none


   ! Define local variables.
   integer :: doneMOME


   ! Initialize the logging labels.
   call initOperationLabels


   ! Parse the command line parameters
   call parseOptcCommandLine


   ! Open the optc files that will be written to.
   if (stateSet /= 1) then
      open (unit=40,file='fort.40',status='new',form='formatted')
   endif
   open (unit=50,file='fort.50',status='new',form='formatted')
   if (spin == 2) then
      if (stateSet /= 1) then
         open (unit=41,file='fort.41',status='new',form='formatted')
      endif
      open (unit=51,file='fort.51',status='new',form='formatted')
   endif


   ! Read in the input to initialize all the key data structure variables.
   call parseInput

   ! Find specific computational parameters not EXPLICITLY given in the input
   !   file.  These values can, however, be easily determined from the input
   !   file.
   call getImplicitInfo

   ! Create real-space superlattice out of the given cell.
   call initializeLattice(0)

   ! Initialize data structures to find arbitrary lattice points easily.
   call initializeFindVec

   ! Compute the kpoint phase factors.
   call computePhaseFactors

   ! Determine if it is necessary to compute the momentum matrix elements or
   !   if they have already been computed.
   call getMOMEStatus (doneMOME)

   ! Compute the momentum matrix elements if necessary.
   if (doneMOME == 0) then
      call addOnMOME
   endif

   ! Access the HDF5 data stored by band.
   call accessPSCFBandHDF5(numStates)

   ! Allocate space to store the energy eigen values, and then read them in.
   allocate (energyEigenValues (numStates,numKPoints,spin))
   call readEnergyEigenValuesBand(numStates)

   ! Populate the electron states to find the highest occupied state (Fermi
   !   energ for metals).
   call populateStates

   if (stateSet == 1) then ! Doing PACS calculation and we need to modify
      ! the energy eigen values and their position. Now that the occupied
      ! energy is known we can append the unoccupied excited states.
      call appendExcitedEnergyEigenValuesBand(lastInitStatePACS+1,&
            & numStates)
   endif

   ! Shift the energy eigen values according to the highest occupied state.
   call shiftEnergyEigenValues(occupiedEnergy,numStates)

   ! Compute some statistics and variables concerning the energy values.
   call getEnergyStatistics

   ! Compute the transition pairs and energy values of those transitions.
   call computeTransitions

   ! Print the output if necessary.
   if (stateSet /= 2) then  ! Not doing a Sigma(E) calculation.
      call printOptcResults
   endif

   ! Deallocate unused matrices
   deallocate (energyEigenValues)
   deallocate (transCounter)
   if (stateSet /= 2) then  ! Not doing a sigma(E) calculation.
      deallocate (energyDiff)
      deallocate (transitionProb)
   endif

   ! Close access to the band HDF5 data.
   call closeAccessPSCFBandHDF5

   ! Open a file to signal completion of the program.
   open (unit=2,file='fort.2',status='new')

end subroutine computeOptc


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


subroutine addOnMOME

   ! Import the necessary HDF modules
   use HDF5
   use O_PSCFIntgHDF5
   use O_IntegralsPSCF
   use O_CommandLine
   use O_Lattice
   use O_Basis

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Set the doMOME command line parameter to 1 (even though it was not
   !   given on the command line).  This will allow us to use the intgAndOrMom
   !   subroutine.
   doMOME = 1

   ! Prepare the necessary data structures for computing the MOME.
   call renormalizeBasis

   ! Open up the integral HDF5 file so that we can write the MOME to it.
   call openMOMEPSCFIntgHDF5

   ! Compute momentum matrix integrals and add them to the integral HDF5 file.
   call intgAndOrMom(0,doMOME)

   ! Record that the momentum matrix elements were calculated.
   call setMOMEStatus(doMOME)

   ! Close the integral HDF5 file so that the rest of optc can proceed as if
   !   this never happened (except now of course the file has the momentum
   !   matrix elements).
   call closeMOMEPSCFIntgHDF5

end subroutine addOnMOME

end module O_Optc
