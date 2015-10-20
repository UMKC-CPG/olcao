module O_Main

contains

subroutine mainSCF (totalEnergy, fromExternal)

   ! Import necessary modules.
   use O_Kinds
   use O_SetupHDF5,       only: accessSetupHDF5, closeSetupHDF5
   use O_MainHDF5,        only: initMainHDF5, closeMainHDF5
   use O_CommandLine,     only: doDOS, doBond, parseMainCommandLine
   use O_Input,           only: numStates, iterFlagTDOS, parseInput
   use O_Potential,       only: converged, currIteration, lastIteration, spin, &
                              & initPotCoeffs, cleanUpPotential
   use O_PotentialUpdate, only: makeSCFPot
   use O_AtomicSites,     only: valeDim, cleanUpAtomSites
   use O_AtomicTypes,     only: cleanUpAtomTypes
   use O_PotSites,        only: cleanUpPotSites
   use O_PotTypes,        only: cleanUpPotTypes
   use O_SecularEquation, only: secularEqnAllKP, cleanUpSecularEqn
   use O_ValeCharge,      only: makeValenceRho
   use O_KPoints,         only: cleanUpKPoints
   use O_Populate,        only: populateStates
   use O_LAPACKParameters, only: setBlockSize
   use O_ExchangeCorrelation, only: cleanUpExchCorr
   use O_DOS,  only: computeIterationTDOS, printIterationTDOS, computeDOS
   use O_Bond, only: computeBond
   use O_TimeStamps, only: initOperationLabels

   ! Import the HDF5 module
   use HDF5

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the parameters that are passed to this subroutine.
   real (kind=double) :: totalEnergy
   integer :: fromExternal

   ! Define local variables.
   integer :: i
   integer :: OLCAOkill
   integer :: hdferr

   ! Open (almost) all the text files that will be written to in this program.
   open (unit=7,file='fort.7',status='unknown',form='formatted')
   open (unit=8,file='fort.8',status='old',form='formatted')
   if (spin == 2) then
      open (unit=13,file='fort.13',status='unknown',form='formatted')
   endif
   open (unit=14,file='fort.14',status='unknown',form='formatted')

   ! Initialize tau timer and start it.
!   call TAU_PROFILE_INIT()
!   call TAU_PROFILE_TIMER(profiler, 'setup')
!   call TAU_PROFILE_START(profiler)
!   call TAU_PROFILE_SET_NODE(1)


   ! Initialize the logging labels.
   call initOperationLabels


   ! Parse the command line parameters
   call parseMainCommandLine


   ! Read in the input to initialize all the key data structure variables.
   call parseInput


   ! Find specific computational parameters not EXPLICITLY given in the input
   !   file.  These values can, however, be easily determined from the input
   !   file.
   call getImplicitInfo


   ! Prepare the HDF5 files for main.
   call initMainHDF5(numStates)


   ! Access the setup hdf5 files.
   call accessSetupHDF5


   ! Note the columns that record statistics for each iteration.  The first is
   !   the iteration
   if (spin == 1) then
      write (7,*)"Iter. Occ. Energy   El. Error    Convergence     Total Energy"
   else
      write (7,*)"Iter. Occ. Energy   El. Error    Convergence     Total Energy      Mag. Mom."
   endif
   call flush (7)

   ! Note the columns that record energy statistics for each iteration.
   write (14,*) "Iter.  Kinetic   ElecStat   ExchCorr   Total"
   call flush (14)


   ! Obtain the initial potential coefficients.
   call initPotCoeffs


   ! Set up LAPACK machine parameters.
   call setBlockSize(valeDim)


   ! Begin the self consistance iterations

   do while (.true.)

      ! Solve the schrodinger equation
      do i = 1, spin
         call secularEqnAllKP(i,numStates)
      enddo


      ! Find the Fermi level, and the number of occupied bands, and the
      !   number of electrons occupying each of those bands.
      call populateStates


      ! Compute the TDOS if it was requested for each iteration.
      if (iterFlagTDOS == 1) then
         call computeIterationTDOS
      endif


      ! Calculate the valence charge density
      call makeValenceRho


      ! Compute the new self consistant potential
      call makeSCFPot(totalEnergy)


      ! Determine if the computation is complete
      if ((converged == 1) .or. (currIteration > lastIteration)) exit

      ! Check if the job was requested to stop.
      open (unit=666,file="OLCAOkill",status="OLD",IOSTAT=OLCAOkill)
      if (OLCAOkill == 0) then
         exit
      endif
   enddo



   ! Print the accumulated TDOS in a useful format if requested in the input.
   if (iterFlagTDOS == 1) then
      call printIterationTDOS
   endif



   ! Compute the final band structure if convergence was achieved or if the
   !   last iteration was done.  (I.e. we are out of the SCF loop.)
   do i = 1, spin
      call secularEqnAllKP(i, numStates)
   enddo

   ! Find the Fermi level, the number of occupied bands, and the number of
   !   electrons occupying each of those bands.
   if ((doDOS .eq. 1) .or. (doBond .eq. 1)) then
      call populateStates
   endif

   ! Compute the dos if it was requested.  For spin polarized calculations the
   !   spin up is in 60, 70, 80 and spin down is in 61, 71, 81.
   if (doDOS == 1) then
      open (unit=60,file='fort.60',status='new',form='formatted') ! TDOS
      open (unit=70,file='fort.70',status='new',form='formatted') ! PDOS
      open (unit=80,file='fort.80',status='new',form='formatted') ! LI
      if (spin == 2) then
         open (unit=61,file='fort.61',status='new',form='formatted') !Spin TDOS
         open (unit=71,file='fort.71',status='new',form='formatted') !Spin PDOS
         open (unit=81,file='fort.81',status='new',form='formatted') !Spin LI
      endif
      call computeDos
   endif

   ! Compute the bond order if it was requested.
   if (doBond == 1) then
      open (unit=10,file='fort.10',status='new',form='formatted')
      if (spin == 2) then
         open (unit=11,file='fort.11',status='new',form='formatted')
      endif
      call computeBond
   endif

   ! Close the HDF objects that were used.
   call closeSetupHDF5
   call closeMainHDF5

   ! Close the HDF5 interface.
   call h5close_f (hdferr)
   if (hdferr /= 0) stop 'Failed to close the HDF5 interface.'

   ! End the tau timer
!   call TAU_PROFILE_STOP(profiler)

   ! Close the output files
   close (7)
   close (8)
   close (13)
   close (20)

   ! Deallocate all the other as of yet un-deallocated arrays.
   call cleanUpAtomTypes
   call cleanUpAtomSites
   call cleanUpPotTypes
   call cleanUpPotSites
   call cleanUpKPoints
   call cleanUpExchCorr
   call cleanUpPotential
   call cleanUpSecularEqn

   ! Open file to signal completion of the program to the calling olcao script.
   open (unit=2,file='fort.2',status='unknown')

end subroutine mainSCF


subroutine getImplicitInfo

   ! Import necessary modules.
   use O_TimeStamps
   use O_ExchangeCorrelation, only: makeSampleVectors
   use O_AtomicTypes,         only: getAtomicTypeImplicitInfo
   use O_AtomicSites,         only: getAtomicSiteImplicitInfo
   use O_PotSites,            only: getPotSiteImplicitInfo
   use O_PotTypes,            only: getPotTypeImplicitInfo
   use O_Lattice,             only: getRecipCellVectors
   use O_KPoints,             only: convertKPointsToXYZ
   use O_Potential,           only: initPotStructures
   use O_Populate,            only: initCoreStateStructures

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

   call initCoreStateStructures

   call timeStampEnd(2)

end subroutine getImplicitInfo


end module O_Main
