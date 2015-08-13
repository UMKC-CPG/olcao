module O_Band


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define module data.
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: coreValeOL
#else
   real (kind=double), allocatable, dimension (:,:)      :: coreValeOLGamma
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine band

   ! Import the necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_AtomicSites,      only: valeDim
   use O_LAPACKParameters, only: setBlockSize
   use O_Input,            only: numStates, parseInput
   use O_CommandLine,      only: parseBandCommandLine, doSYBD
   use O_PSCFBandHDF5,     only: initPSCFBandHDF5, closePSCFBandHDF5
   use O_Lattice,          only: initializeLattice, initializeFindVec
   use O_KPoints,          only: makePathKPoints, computePhaseFactors

   ! Make sure that there are no accidental variable declarations.
   implicit none


   ! Initialize the logging labels.
   call initOperationLabels


   ! Parse the command line parameters
   call parseBandCommandLine


   ! Read in the input to initialize all the key data structure variables.
   call parseInput


   ! Find specific computational parameters not EXPLICITLY given in the input
   !   file.  These values can, however, be easily determined from the input
   !   file.
   call getImplicitInfo


   ! Set up LAPACK machine parameters.
   call setBlockSize(valeDim)



   ! In the case of a symmetric band structure calculation we must modify the
   !   kpoints according to the lattice type information given in the
   !   olcao.dat input file.
   if (doSYBD == 1) then
      call makePathKPoints
   endif


   ! Create real-space super lattice out of the primitive lattice.  These
   !   "supercells" must be big enough so as to include all the points within
   !   sphere bounded by the negligability limit.  Points outside the sphere
   !   are considered negligable and are ignored.
   call initializeLattice (0)


   ! Setup the necessary data structures so that we can easily find the lattice
   !   vector that is closest to any other arbitrary vector.
   call initializeFindVec


   ! Compute the kpoint phase factors.
   call computePhaseFactors
   call flush(20)


   ! Prepare the HDF5 files for the post SCF band structure calculation.
   call initPSCFBandHDF5(numStates)


   ! Compute the solid state wave function.
   call computeBands


   ! Close the HDF5 band structure file.
   call closePSCFBandHDF5


   ! Open a file to signal completion of the program.
   open (unit=2,file='fort.2',status='new')


end subroutine band



subroutine computeBands

   ! Use necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Potential,     only: spin
   use O_CommandLine,   only: doSYBD
   use O_Input,         only: numStates
   use O_KPoints,       only: numKPoints
   use O_IntegralsPSCF, only: getIntgResults
   use O_AtomicSites,   only: coreDim, valeDim
   use O_PSCFBandHDF5,  only: valeValeBand, valeValeBand_did, saveCoreValeOL
#ifndef GAMMA
   use O_SecularEquation, only: valeVale, valeValeOL, energyEigenValues, &
         & preserveValeValeOL, restoreValeValeOL, secularEqnOneKP
#else
   use O_SecularEquation, only: valeValeGamma, valeValeOLGamma, &
         & energyEigenValues, preserveValeValeOL, restoreValeValeOL, &
         & secularEqnOneKP
#endif

   ! Define local variables.
   integer :: i,j

   ! Record the date and time we start.
   call timeStampStart(15)


   ! Allocate space to hold the energyEigenValues, coreValeOL, and valeVale
   !   matrices.
   allocate (energyEigenValues(numStates,numKPoints,spin))
#ifndef GAMMA
   allocate (valeVale(valeDim,valeDim,1,spin))
   allocate (valeValeOL(valeDim,valeDim,1,spin))
   allocate (coreValeOL(coreDim,valeDim,1))
#else
   allocate (valeValeGamma(valeDim,valeDim,spin))
   allocate (valeValeOLGamma(valeDim,valeDim,spin))
   allocate (coreValeOLGamma(coreDim,valeDim))
#endif

!write (20,*) "numKPoints=",numKPoints
!write (20,*) "valeDim=",valeDim
!call flush (20)

   ! Loop over all kpoints and compute the solid state wave function and
   !   energy eigen values for each. But first, ...
   ! Record the fact that we are starting the k-point loop so that when
   !   someone looks at the output file as the job is running they will
   !   know that all the little dots represent a count of the number of
   !   k-points. Then they can figure out the progress and progress rate.
   write (20,*) "Beginning k-point loop."
   if (numKPoints > 1) write (20,*) "Expecting ",numKPoints," iterations."
   call flush (20)

   do i = 1, numKPoints
!write (20,*) "current kpoints = ",i
!call flush (20)

      ! Read the overlap integral results and apply kpoints effects.
#ifndef GAMMA
      call getIntgResults(valeValeOL(:,:,:,1),coreValeOL,i,1,&
            & valeValeBand_did(i),valeValeBand,doSYBD,1)
#else
      call getIntgResults(valeValeOLGamma(:,:,1),coreValeOLGamma,1,&
            & valeValeBand_did(i),valeValeBand,doSYBD,1)
#endif
!write (20,*) "Got here 0"
!call flush (20)

      ! Write the coreValeOL for later use by other programs (if needed).
#ifndef GAMMA
      call saveCoreValeOL (coreValeOL,i)
#else
      call saveCoreValeOL (coreValeOLGamma,i)
#endif

      ! Preserve the valeValeOL (if needed).
      call preserveValeValeOL

      do j = 1, spin
#ifndef GAMMA
         ! Read the hamiltonian integral results and apply kpoint effects.
         call getIntgResults(valeVale(:,:,:,j),coreValeOL,i,2,&
               & valeValeBand_did(i),valeValeBand,doSYBD,j)
#else
         ! Read the hamiltonian integral results and apply kpoint effects.
         call getIntgResults(valeValeGamma(:,:,j),coreValeOLGamma,2,&
               & valeValeBand_did(i),valeValeBand,doSYBD,j)
#endif

         ! Solve the wave equation for this kpoint.
         call secularEqnOneKP(j,i,numStates,doSYBD)

         ! Restore the valeValeOL (if needed).
         if (j == 1) then
            call restoreValeValeOL
         endif
      enddo

      ! Record that this kpoint has been finished.
      if (mod(i,10) .eq. 0) then
         write (20,ADVANCE="NO",FMT="(a1)") "|"
      else
         write (20,ADVANCE="NO",FMT="(a1)") "."
      endif
      if (mod(i,50) .eq. 0) then
         write (20,*) " ",i
      endif
      call flush (20)

   enddo

   ! Print the band results if necessary.
   if (doSYBD == 1) then
      call printSYBD
   endif

   ! Record the date and time we end.
   call timeStampEnd(15)

end subroutine computeBands


subroutine printSYBD

   ! Use necessary modules.
   use O_Potential,       only: spin
   use O_Constants,       only: hartree
   use O_Input,           only: numStates
   use O_KPoints,         only: numPathKP, pathKPointMag
   use O_Populate,        only: occupiedEnergy, populateStates
   use O_SecularEquation, only: energyEigenValues, shiftEnergyEigenValues

   ! Define local variables.
   integer :: i,j
   character*7 :: filename

   ! Obtain the fermi energy with the populate energy levels subroutine.
   call populateStates

   ! Adjust the energyEigenValues down by the fermi level determined above.
   call shiftEnergyEigenValues(occupiedEnergy,numStates)

   do i = 1, spin

      ! Create the filename.
      write (filename,fmt="(a5,i2)") "fort.",30+i

      ! Open the output file.
      open (unit=30+i,file=filename,status='new',form='formatted')

      ! Record the number of kpoints and the number of states
      write (30+i,*) numPathKP, numStates

      ! Record each kpoint's energy values
      do j = 1, numPathKP

         ! Record the distance of this kpoint from the beginning of the path.
         write (30+i,*) pathKPointMag(j)

         ! Record the energy values for this KPoint
         write (30+i,fmt="(10f12.5)") energyEigenValues(:numStates,j,i)*hartree
      enddo

      ! Close the output file.
      close (30+i)

   enddo

end subroutine printSYBD


subroutine getImplicitInfo

   ! Import necessary modules.
   use O_ExchangeCorrelation, only: makeSampleVectors
   use O_AtomicSites,         only: getAtomicSiteImplicitInfo
   use O_AtomicTypes,         only: getAtomicTypeImplicitInfo
   use O_PotSites,            only: getPotSiteImplicitInfo
   use O_PotTypes,            only: getPotTypeImplicitInfo
   use O_Lattice,             only: getRecipCellVectors
   use O_KPoints,             only: convertKPointsToXYZ
   use O_Potential,           only: initPotStructures
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


end module O_Band
