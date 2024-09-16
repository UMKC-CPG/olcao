module O_OLCAO

contains

subroutine OLCAO

   ! Import necessary modules. Note that 
   use O_SCFHDF5
   use O_CommandLine
   use O_TimeStamps, only: initOperationLabels
   use O_Input, only: parseInput
   use O_LocalEnv

   ! Initialize the logging labels.
   call initOperationLabels

   ! Parse the command line parameters
   call parseCommandLine

   ! Read in the input to initialize all the key data structure variables.
   call parseInput

   ! Find specific computational parameters not EXPLICITLY given in the input
   !   file, but which can, however, be easily determined from the input file.
   call getImplicitInfo

   ! Now, we are ready to do *either* SCF or Post SCF work.
   if (doSCF == 1) then
      call setupSCF ! Preparation for SCF cycle and wave function calculation.

      call mainSCF ! The actual SCF cycle and wave function calculation.

      if (doDOS_SCF == 1) then
         call dos(1)  ! Passing inSCF == 1.
      endif

      if (doBond_SCF >= 1) then
         call bond(1, doBond_SCF)  ! Passing inSCF == 1
      endif

!      if (doDIMO_SCF == 1) then
!         call dimo 
!      endif

      if (doOPTC_SCF >= 1) then
         call optc(1, doOPTC_SCF)  ! Passing inSCF == 1
      endif

      if (doField_SCF == 1) then
         !call field(1, doField_SCF)
      endif

      if (doPSCF == 1) then
         !call reset
      endif

      call cleanUpSCF

      call closeHDF5_SCF
   endif

   if (doPSCF == 1) then
      call intgPSCF ! Preparations for wave function calculation.

      call bandPSCF(doSYBD_PSCF) ! Wave function calculation.

!      if (doDIMO_PSCF == 1) then
!         call dimo
!      endif

      if (doDOS_PSCF == 1) then
         call dos(0)
      endif

      if (doBond_PSCF >= 1) then
         call bond(0, doBond_PSCF)
      endif

      if (doOptc_PSCF >= 1) then
         call optc(0, do_OPTC_PSCF)
      endif

      call closeHDF5_PSCF
   endif

   if (doLoEn == 1) then
      call analyzeLocalEnv
   endif

end subroutine OLCAO


subroutine setupSCF
   
   ! Import the precision variables
   use O_Kinds

   ! Import the necessary modules.
   use O_SCFHDF5, only: initHDF5_SCF
   use O_CommandLine, only: doDIMO_SCF, doOPTC_SCF
   use O_Input, only: numStates
   use O_Lattice, only: initializeLattice, initializeFindVec, &
         & cleanUpLattice
   use O_KPoints, only: numKPoints, computePhaseFactors, cleanUpKPoints
   use O_Basis, only: renormalizeBasis, cleanUpBasis
   use O_ExchangeCorrelation, only: maxNumRayPoints, getECMeshParameters, &
         & makeECMeshAndOverlap, cleanUpExchCorr
   use O_ElectroStatics, only: makeElectrostatics
   use O_GaussianRelations, only: makeAlphaDist, makeAlphaNucDist, &
         & makeAlphaPotDist, cleanUpGaussRelations
   use O_Integrals,   only: allocateIntegralsSCF, gaussOverlapOL, &
         & gaussOverlapKE, gaussOverlapMV, gaussOverlapNP, &
         & elecPotGaussOverlap, cleanUpIntegralsSCF, &
         & secondCleanUpIntegralsSCF
   use O_Integrals3Terms,  only: allocateIntegralsSCF3Terms, &
         & gaussOverlapDM, gaussOverlapMM, cleanUpIntegralsSCF3Terms
   use O_AtomicSites, only: coreDim, valeDim, cleanUpAtomSites
   use O_AtomicTypes, only: cleanUpRadialFns, cleanUpAtomTypes
   use O_PotSites, only: cleanUpPotSites
   use O_PotTypes, only: cleanUpPotTypes
   use O_Potential, only: rel, cleanUpPotential
   use O_CoreCharge, only: makeCoreRho

   ! Import the HDF5 module.
   use HDF5

   ! Import the necessary MPI files
!   use MPI

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define global mpi parameters
!   integer :: myWorldPid
!   integer :: numWorldProcs

   ! Define tau parameters
!   integer profiler(2) / 0, 0 /
!   save profiler

   ! Define error detectors for hdf, tau, and mpi.
   integer :: hdferr
!   integer :: tauerr
!   integer :: mpierr


   ! Initialize the MPI interface
!   call MPI_INIT (mpierr)
!   call MPI_COMM_RANK (MPI_COMM_WORLD,myWorldPid,mpierr)
!   call MPI_COMM_SIZE (MPI_COMM_WORLD,numWorldProcs,mpierr)

   ! Initialize tau timer and start it.

!   call TAU_PROFILE_INIT()
!   call TAU_PROFILE_TIMER(profiler, 'setup')
!   call TAU_PROFILE_START(profiler)
!   call TAU_PROFILE_SET_NODE(0)

   ! Read in the input to initialize all the key data structure variables.
   call parseInput


   ! Find specific computational parameters not EXPLICITLY given in the input
   !   file.  These values can, however, be easily determined from the input
   !   file.
   call getImplicitInfo


   ! Create real-space and reciprocal-space super lattices out of the primitive
   !   lattice.  These "supercells" must be big enough so as to include all the
   !   points within a sphere bounded by the negligability limit.  Points
   !   outside the sphere are considered negligable.
   call initializeLattice (1)


   ! Setup the necessary data structures so that we can easily find the lattice
   !   vector that is closest to any other arbitrary vector.
   call initializeFindVec


   ! Compute the kpoint phase factors.
   call computePhaseFactors


   ! Renormalize the basis functions
   call renormalizeBasis


   ! Determine the parameters for the exchange correlation mesh.
   call getECMeshParameters


   ! Now, the dimensions of the system are known.  Therefore we can
   !   initialize the HDF5 file structure format, and datasets.
   call initHDF5_SCF (maxNumRayPoints, numStates)


   ! Construct the exchange correlation overlap matrix, and sampling field.
   call makeECMeshAndOverlap


   ! Create the alpha distance matrices.
   call makeAlphaDist


   ! Allocate space to be used for each of the single matrix integrals.
   call allocateIntegralsSCF(coreDim,valeDim,numKPoints)


   ! Calculate the matrix elements of the overlap between all LCAO Bloch
   !   wave functions.
   call gaussOverlapOL


   ! Calculate the matrix elements of the kinetic energy between all LCAO Bloch
   !   basis functions.
   call gaussOverlapKE


   ! Calculate the matrix elements of the mass velocity between all LCAO Bloch
   !   basis functions if needed for the scalar relativistic calculation.
   if (rel == 1) then
      call gaussOverlapMV
   endif


   ! Create the alpha distance matrix with nuclear alpha factor
   call makeAlphaNucDist


   ! Calculate the matrix elements of the overlap between all LCAO Bloch
   !   wave functions and the nuclear potentials.
   call gaussOverlapNP


   ! Create the alpha distance matrix with potential alpha factor
   call makeAlphaPotDist


   ! Calculate the matrix elements of the overlap between all LCAO Bloch
   !   wave functions and the potential site potential alphas.
   call elecPotGaussOverlap


   ! Now that all the single matrices are done being made we can deallocate
   !   the data structures that were used in all the above subroutines but are
   !   not necessary now.
   call cleanUpIntegralsSCF

   ! If any supplementary three-term matrices (xyz) were requested for
   !   computing other properties, then allocate the space that is necessary
   !   for the integral calculation and do the integral.
   if ((doDIMO_SCF == 1) .or. (doOPTC_SCF >= 1)) then

      ! Consider in the future an option to do the XYZ independently to
      !   conserve memory if that becomes a problem.
      call allocateIntegralsSCF3Terms(coreDim,valeDim,numKPoints)

      ! Do the relevant integrals.
      if (doDIMO_SCF == 1) then
         call gaussOverlapDM
      endif
      if (doOPTC_SCF >= 1) then
         call gaussOverlapMM
      endif

      ! Clean up from the 3Term integrals.
      call cleanUpIntegralsSCF3Terms
   endif


   call secondCleanupIntegralsSCF ! Overlap matrix parts used for ortho.
   call cleanUpBasis
   call cleanUpGaussRelations


   ! Construct a vector describing the core charge density since it will not
   !   change throughout the SCF iterations because of core orthogonalization.
   if (coreDim /= 0) then
      call makeCoreRho
   endif


   ! Deallocate the data component of the atomType data structure that holds
   !   the radial functions since they are no longer needed.
   call cleanUpRadialFns


   ! Construct matrix operators and integral vectors that are used to later
   !   determine the electrostatic potential.
   call makeElectrostatics


   !! Close all the parts of the setup HDF5 file.
   !call closeSetupHDF5


   !! Close the HDF5 interface.
   !call h5close_f (hdferr)
   !if (hdferr /= 0) stop 'Failed to close the HDF5 interface.'

   !
   !! Deallocate all the other as of yet un-deallocated arrays.
   !call cleanUpAtomTypes
   !call cleanUpAtomSites
   !call cleanUpPotTypes
   !call cleanUpPotSites
   !call cleanUpKPoints
   !call cleanUpExchCorr
   !call cleanUpLattice
   !call cleanUpPotential

   !! End the tau timer
!  ! call TAU_PROFILE_STOP(profiler)

   !! End the MPI interface
!  ! call MPI_FINALIZE (mpierr)

   !! Close the output file
   !close (20)

   !! Open a file to signal completion of the program.
   !open (unit=2,file='fort.2',status='unknown')

end subroutine setupSCF


subroutine getImplicitInfo

   ! Import necessary modules.
   use O_TimeStamps
   use O_ExchangeCorrelation, only: makeSampleVectors
   use O_AtomicTypes, only: getAtomicTypeImplicitInfo
   use O_AtomicSites, only: getAtomicSiteImplicitInfo
   use O_PotSites, only: getPotSiteImplicitInfo
   use O_PotTypes, only: getPotTypeImplicitInfo
   use O_Lattice, only: getRecipCellVectors
   use O_KPoints, only: convertKPointsToXYZ
   use O_Potential, only: initPotStructures
   use O_Populate, only: initCoreStateStructures
   use O_Input, only: getDipoleMomentCenter

   implicit none

   call timeStampStart(2)

   ! If a dipole moment calculation is being done, then compute the
   !   coordinates of the moment center "C" in Cartesian coordinates.
   call getDipoleMomentCenter

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


subroutine mainSCF

   ! Import necessary modules.
   use HDF5
   use O_Kinds
   use O_CommandLine, only: doDOS_SCF, doBond_SCF, doDIMO_SCF, doForce_SCF
   use O_Input, only: numStates, iterFlagTDOS
   use O_Potential, only: converged, currIteration, lastIteration, &
         & spin, rel, initPotCoeffs, cleanUpPotential
   use O_PotentialUpdate, only: makeSCFPot
   use O_AtomicSites, only: valeDim, cleanUpAtomSites
   use O_AtomicTypes, only: cleanUpAtomTypes
   use O_PotSites, only: cleanUpPotSites
   use O_PotTypes, only: cleanUpPotTypes
   use O_SecularEquation, only: secularEqnAllKP, cleanUpSecularEqn, &
         & shiftEnergyEigenValues
   use O_ValeCharge, only: makeValenceRho
   use O_KPoints, only: cleanUpKPoints
   use O_Populate, only: occupiedEnergy, populateStates
   use O_LAPACKParameters, only: setBlockSize
   use O_ExchangeCorrelation, only: cleanUpExchCorr
   use O_DOS, only: computeIterationTDOS, printIterationTDOS, computeDOS
   use O_Bond, only: computeBond
   use O_Force, only: computeForceIntg

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define local variables.
   integer :: i
   integer :: OLCAOkill
   integer :: hdferr
   real (kind=double) :: totalEnergy

   ! Open (almost) all the text files that will be written to in this program.
   open (unit=7,file='fort.7',status='unknown',form='formatted')
   open (unit=8,file='fort.8',status='old',form='formatted')
   if (spin == 2) then
      open (unit=13,file='fort.13',status='unknown',form='formatted')
   endif
   open (unit=14,file='fort.14',status='unknown',form='formatted')

   ! Note the columns that record statistics for each iteration.  The first is
   !   the iteration
   if (spin == 1) then
      write (7,*) "Iter. Occ. Energy   El. Error    ", &
            & "Convergence     Total Energy"
   else
      write (7,*) "Iter. Occ. Energy   El. Error    ", &
            & "Convergence     Total Energy      Mag. Mom."
   endif
   call flush (7)

   ! Note the columns that record energy statistics for each iteration.
   if (rel == 0) then
      write (14,*) "Iter. Kinetic        ElecStat       ", &
            & "ExchCorr       Total"
   else
      write (14,*) "Iter. Kinetic        MassVel        ", &
            & "ElecStat       ExchCorr       Total"
   endif
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

!   ! Shift the energy eigen values according to the highest occupied state.
!   call shiftEnergyEigenValues(occupiedEnergy,numStates)
!
!   ! Find the Fermi level, the number of occupied bands, and the number of
!   !   electrons occupying each of those bands.
!   if ((doDOS == 1) .or. (doBond == 1) .or. (doDIMO == 1) .or. &
!         & (doForce == 1)) then
!      call populateStates
!   endif
!
!   ! Compute the dos if it was requested.  For spin polarized calculations the
!   !   spin up is in 60, 70, 80 and spin down is in 61, 71, 81.
!   if (doDOS_SCF == 1) then
!      open (unit=60,file='fort.60',status='new',form='formatted') ! TDOS
!      open (unit=70,file='fort.70',status='new',form='formatted') ! PDOS
!      open (unit=80,file='fort.80',status='new',form='formatted') ! LI
!      if (spin == 2) then
!         open (unit=61,file='fort.61',status='new',form='formatted') !Spin TDOS
!         open (unit=71,file='fort.71',status='new',form='formatted') !Spin PDOS
!         open (unit=81,file='fort.81',status='new',form='formatted') !Spin LI
!      endif
!      call computeDos
!   endif
!
!   ! Compute the bond order if it was requested.
!   if (doBond_SCF == 1) then
!      open (unit=10,file='fort.10',status='new',form='formatted')
!      if (spin == 2) then
!         open (unit=11,file='fort.11',status='new',form='formatted')
!      endif
!      call computeBond
!   endif
!
!   ! Use the dipole moment matrix to compute the dipole.
!   if (doDIMO_SCF == 1) then
!      open (unit=74,file='fort.74',status='new',form='formatted')
!      if (spin == 2) then
!         open (unit=75,file='fort.75',status='new',form='formatted')
!      endif
!      call makeValenceRho
!   endif
!
!   ! Compute forces between atoms.
!   if (doForce_SCF == 1) then
!      open (unit=98,file='fort.98',status='new',form='formatted')
!      if (spin == 2) then
!         open (unit=99,file='fort.99',status='new',form='formatted')
!      endif
!      call computeForceIntg(totalEnergy)
!      call makeValenceRho
!   endif

!   ! Close the HDF objects that were used.
!   call closeHDF5_SCF
!
!   ! Close the HDF5 interface.
!   call h5close_f (hdferr)
!   if (hdferr /= 0) stop 'Failed to close the HDF5 interface.'

end subroutine mainSCF


subroutine intgPSCF

   ! Import the necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Potential, only: initPotCoeffs
   use O_Basis,     only: renormalizeBasis
   use O_Integrals, only: allPSCFIntgCombo
   use O_PSCFHDF5,  only: initHDF5_PSCF
   use O_Lattice,   only: initializeLattice, initializeFindVec

   ! Make sure that there are not accidental variable declarations.
   implicit none


   ! Open the potential file that will be read from in this program.
   open (unit=8,file='fort.8',status='old',form='formatted')


   ! Create real-space super lattice out of the primitive lattice.  These
   !   "supercells" must be big enough so as to include all the points within
   !   sphere bounded by the negligability limit.  Points outside the sphere
   !   are considered negligable and are ignored.
   call initializeLattice (1)


   ! Setup the necessary data structures so that we can easily find the lattice
   !   vector that is closest to any other arbitrary vector.
   call initializeFindVec


   ! Renormalize the basis functions.
   call renormalizeBasis


   ! Prepare the HDF5 files for the post SCF calculations.
   call initHDF5_PSCF


   ! Read the provided potential coefficients.
   call initPotCoeffs

   ! Compute the integrals and the momentum matrix elements (MME or MOME) (if
   !   requested).  The MME request is given on the command line and the
   !   subroutine gets that flag from that module.  The integrals are also
   !   optional and are controlled by the intgAndOrMom parameter.  Obviously,
   !   in this program the whole point is to compute the integrals.  However,
   !   in other programs (optc) the integrals do not need to be computed and
   !   only the MME need to be computed.  In that case the same subroutine can
   !   be used, except that only the MME are computed, not the overlap and
   !   hamiltonian.
   call allPSCFIntgCombo


end subroutine intgPSCF


subroutine bandPSCF

      ! Import the necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Input,            only: numStates
   use O_Potential,        only: spin
   use O_AtomicSites,      only: valeDim
   use O_LAPACKParameters, only: setBlockSize
   use O_Input,            only: numStates, parseInput
   use O_CommandLine,      only: doSYBD_PSCF
   use O_Lattice,          only: initializeLattice, initializeFindVec
   use O_KPoints,          only: makePathKPoints, numKPoints, computePhaseFactors
!   use O_PSCFBandHDF5,  only: valeValeBand, valeValeBand_did, saveCoreValeOL
#ifndef GAMMA
   use O_SecularEquation, only: valeVale, valeValeOL, energyEigenValues, &
         & preserveValeValeOL, restoreValeValeOL, secularEqnOneKP
#else
   use O_SecularEquation, only: valeValeGamma, valeValeOLGamma, &
         & energyEigenValues, preserveValeValeOL, restoreValeValeOL, &
         & secularEqnOneKP
#endif
   use HDF5 ! Needed to define the kludgeInt hid_t integer.

   ! Define local variables.
   integer :: i,j
   integer (hid_t) :: kludgeInt ! Place holder integer to avoid accessing
         ! unallocated memory.

   ! Make sure that there are no accidental variable declarations.
   implicit none


   ! Set up LAPACK machine parameters.
   call setBlockSize(valeDim)


   ! In the case of a symmetric band structure calculation we must modify the
   !   kpoints according to the lattice type information given in the
   !   olcao.dat input file.
   if (doSYBD == 1) then
      call makePathKPoints
   endif


!   ! Create real-space super lattice out of the primitive lattice.  These
!   !   "supercells" must be big enough so as to include all the points within
!   !   sphere bounded by the negligability limit.  Points outside the sphere
!   !   are considered negligable and are ignored.
!   call initializeLattice (0)
!
!
!   ! Setup the necessary data structures so that we can easily find the lattice
!   !   vector that is closest to any other arbitrary vector.
!   call initializeFindVec
!
!
!   ! Compute the kpoint phase factors.
!   call computePhaseFactors
!   call flush(20)
!
!
!   ! Prepare the HDF5 files for the post SCF band structure calculation.
!   call initPSCFBandHDF5(numStates)


   ! Compute the solid state wave function for each kpoint.

   ! Allocate space to hold the energyEigenValues, coreValeOL, and valeVale
   !   matrices.
   allocate (energyEigenValues(numStates,numKPoints,spin))
#ifndef GAMMA
   allocate (valeVale(valeDim,valeDim,1,spin))
   allocate (valeValeOL(valeDim,valeDim,1,spin))
#else
   allocate (valeValeGamma(valeDim,valeDim,spin))
   allocate (valeValeOLGamma(valeDim,valeDim,spin))
#endif


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

      ! Read the overlap integral results and apply kpoints effects. Note
      !   that the doSYBD if statement is a bit of a kludge because we do
      !   not allocate the valeValeBand_did array in the doSYBD==1 case
      !   so it is improper to try to access and pass an unallocated data
      !   value. Hence we send a temporary integer instead. (This is a
      !   placeholder until the partial SYBD code is implemented.)
#ifndef GAMMA
      if (doSYBD == 0) then
         call getIntgResults(valeValeOL(:,:,:,1),coreValeOL,i,1,&
               & valeValeBand_did(i),valeValeBand,1,1)
      else
         call getIntgResults(valeValeOL(:,:,:,1),coreValeOL,i,1,&
               & kludgeInt,valeValeBand,0,1)
      endif
#else
      if (doSYBD == 0) then
         call getIntgResults(valeValeOLGamma(:,:,1),coreValeOLGamma,1,&
               & valeValeBand_did(i),valeValeBand,1,1)
      else
         call getIntgResults(valeValeOLGamma(:,:,1),coreValeOLGamma,1,&
               & kludgeInt,valeValeBand,0,1)
      endif
#endif

      ! Write the coreValeOL for later use by other programs (if needed).
#ifndef GAMMA
      call saveCoreValeOL (coreValeOL,i)
#else
      call saveCoreValeOL (coreValeOLGamma,i)
#endif

      ! Preserve the valeValeOL (if needed). It will be destroyed by the
      !   upcoming call to secularEqnOneKP where the ZHEGV routine is called.
      !   This should be preserved in the case that a spin polarized
      !   calculation is being done.
      call preserveValeValeOL

      do j = 1, spin
#ifndef GAMMA
         ! Read the hamiltonian integral results and apply kpoint effects. See
         !   note above about the doSYBD 'if' statement.
         if (doSYBD == 0) then
            call getIntgResults(valeVale(:,:,:,j),coreValeOL,i,2,&
                  & valeValeBand_did(i),valeValeBand,0,j)
         else
            call getIntgResults(valeVale(:,:,:,j),coreValeOL,i,2,&
                  & kludgeInt,valeValeBand,1,j)
         endif
#else
         ! Read the hamiltonian integral results and apply kpoint effects.
         if (doSYBD == 0) then
            call getIntgResults(valeValeGamma(:,:,j),coreValeOLGamma,2,&
                  & valeValeBand_did(i),valeValeBand,0,j)
         else
            call getIntgResults(valeValeGamma(:,:,j),coreValeOLGamma,2,&
                  & kludgeInt,valeValeBand,1,j)
         endif
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
   call computeBands

end subroutine bandPSCF


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
   use HDF5 ! Needed to define the kludgeInt hid_t integer.

   ! Define local variables.
   integer :: i,j
   integer (hid_t) :: kludgeInt ! Place holder integer to avoid accessing
         ! unallocated memory.

   ! Record the date and time we start.
   call timeStampStart(15)


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
         write (30+i,fmt="(11f14.5)") energyEigenValues(:numStates,j,i)*hartree
      enddo

      ! Close the output file.
      close (30+i)

   enddo

end subroutine printSYBD


subroutine dos(inSCF)

   ! Use necessary modules.
   use O_TimeStamps
   use O_Potential,       only: spin
   use O_DOS,             only: computeDOS
   use O_KPoints,         only: numKPoints
   use O_Input,           only: numStates, parseInput
   use O_Populate,        only: occupiedEnergy, populateStates
!   use O_PSCFBandHDF5,    only: accessPSCFBandHDF5, closeAccessPSCFBandHDF5
   use O_SecularEquation, only: energyEigenValues, &
         & shiftEnergyEigenValues
!   use O_SecularEquation, only: energyEigenValues, readEnergyEigenValuesBand, &
!         & shiftEnergyEigenValues


   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare passed parameters.
   integer, intent(in) :: inSCF

   ! Open the DOS files that will be written to.  If a spin polarized
   !   calculation is being done, then 60, 70, 80 hold spin up and 61, 71, 81
   !   hold spin down.  60,61=TDOS; 70,71= PDOS; 80,81=Localization Index
   open (unit=60,file='fort.60',status='new',form='formatted')
   open (unit=70,file='fort.70',status='new',form='formatted')
   open (unit=80,file='fort.80',status='new',form='formatted')
   if (spin == 2) then
      open (unit=61,file='fort.61',status='new',form='formatted')
      open (unit=71,file='fort.71',status='new',form='formatted')
      open (unit=81,file='fort.81',status='new',form='formatted')
   endif


   ! Populate the electron states to find the highest occupied state (Fermi
   !   energy for metals).
   call populateStates

   ! Shift the energy eigen values according to the highest occupied state.
   call shiftEnergyEigenValues(occupiedEnergy,numStates)

   ! Call the DOS subroutine to compute the total and partial density of states
   !   as well as the localization index.
   call computeDOS(inSCF)

   ! Close the output files.
   close(60)
   close(70)
   close(80)
   if (spin == 2) then
      close(61)
      close(71)
      close(81)
   endif

end subroutine dos


subroutine bond (inSCF, doBond)

   ! Use necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Potential,       only: spin
   use O_Bond,            only: computeBond
   use O_Bond3C,          only: computeBond3C
   use O_Populate,        only: occupiedEnergy, populateStates
   use O_Input,           only: thermalSigma, numStates, parseInput
!   use O_PSCFBandHDF5,    only: accessPSCFBandHDF5, closeAccessPSCFBandHDF5
   use O_SecularEquation, only: shiftEnergyEigenValues


   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare passed parameters.
   integer, intent(in) :: inSCF
   integer, intent(in) :: doBond

   ! Declare local variables.
   real(kind=double) :: thermalSigmaTemp


   ! Open the bond files that will be written to.
   open (unit=10,file='fort.10',status='new',form='formatted')
   if (spin == 2) then
      open (unit=11,file='fort.11',status='new',form='formatted')
   endif
   if (doBond == 2) then
      open (unit=12,file='fort.12',status='new',form='formatted')
      if (spin == 2) then
         open (unit=13,file='fort.13',status='new',form='formatted')
      endif
   endif


   ! The bond calculation must be done without any smearing to determine
   !   the number of electrons for each atom.  So we set the thermal smearing
   !   factor to zero after saving its input value for later.
   thermalSigmaTemp = thermalSigma
   thermalSigma = 0.0_double

   ! Populate the electron states to find the highest occupied state (Fermi
   !   energy for metals).
   call populateStates

   ! Shift the energy eigen values according to the highest occupied state.
   call shiftEnergyEigenValues(occupiedEnergy,numStates)

   ! Call the bond subroutine to compute the bond order and effective charge.
   call computeBond(inSCF)

   ! Compute the three center bond order if requested.
   if (doBond == 2) then
      call computeBond3C(inSCF)
   endif

   ! Restore the thermal sigma to the input vale.
   thermalSigma = thermalSigmaTemp

   ! Close the BOND files that were written to.
   close (10)
   if (spin == 2) then
      close (11)
   endif
   if (doBond == 2) then ! Close the 3-center files.
      close (12)
      if (spin == 2) then
         close (13)
      endif
   endif

end subroutine bond


!subroutine field(inSCF)
!
!! The goal of this program is to produce plottable data of the solid state
!   !   wave function, data derived form it (e.g. the charge density), and/or
!   !   ancillary data such as the potential function. It permits plotting of
!   !   specific wave function states according to energy level by numerically
!   !   evaluating the analytical wave functions, wave functions squared, etc.
!   !   on a defined 3D uniform mesh. The data may be stored in either plain
!   !   text OpenDX format or in XDMF+HDF5 format for Paraview.
!
!   ! Import the necessary modules.
!   use O_Kinds
!   use O_TimeStamps
!   use O_ElementData,     only: initElementData
!   use O_Populate,        only: populateStates
!   use O_Potential,       only: spin, initPotCoeffs
!   use O_Field,           only: computeFieldMesh, cleanUpField
!   use O_SecularEquation, only: energyEigenValues
!   use O_Lattice,         only: initialize3DMesh
!
!   ! Make sure that there are not accidental variable declarations.
!   implicit none
!
!
!   ! Open the potential file that will be read from in this program.
!   open (unit=8,file='fort.8',status='old',form='formatted')
!
!   ! Initialize element data from periodic table of the elements.
!   call initElementData
!
!
!   ! Populate the electron states to find the highest occupied state (Fermi
!   !   energ for metals).
!   call populateStates
!
!
!   ! Initialize certain parameters for constructing and traversing the 3D mesh.
!   call initialize3DMesh
!
!
!   ! Compute the requested field values for each mesh point and store in HDF5.
!   call computeFieldMesh
!
!
!   ! Clean up any left over arrays that need to be deallocated.
!   call cleanUpField
!
!
!end subroutine field


subroutine optc(inSCF,doOPTC)

   ! Use necessary modules.
   use O_TimeStamps
   use O_Potential,       only: spin
   use O_OptcPrint,       only: printOptcResults
   use O_KPoints,         only: numKPoints, computePhaseFactors
   use O_Populate,        only: occupiedEnergy, populateStates
   use O_Lattice,         only: initializeLattice, initializeFindVec
   use O_Input,           only: numStates, lastInitStatePACS, parseInput,&
         & detailCodePOPTC
!   use O_PSCFBandHDF5,    only: accessPSCFBandHDF5, closeAccessPSCFBandHDF5
   use O_OptcTransitions, only: transCounter, energyDiff, transitionProb, &
         & transitionProbPOPTC, getEnergyStatistics, computeTransitions
   use O_SecularEquation, only: energyEigenValues, &
         & shiftEnergyEigenValues
!   use O_SecularEquation, only: energyEigenValues, readEnergyEigenValuesBand,&
!         & appendExcitedEValsBand, shiftEnergyEigenValues


   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: inSCF
   integer, intent(in) :: doOPTC


   ! Define local variables.
   integer :: doneMOME


   ! Open the optc files that will be written to. 40,41=optical conductivity;
   !   50,51 = imaginary dielectric function (usually). Note, for XANES/ELNES
   !   calculations, we do not need the optical conductivity. Also note, for
   !   Sigma(E) calculations the 50,51 units are the Sigma(E) itself, not the
   !   imaginary dielectric function (epsilon2).
   if (doOPTC /= 2) then
      open (unit=40,file='fort.40',status='new',form='formatted')
   endif
   open (unit=50,file='fort.50',status='new',form='formatted')
   if (spin == 2) then
      if (doOPTC /= 2) then
         open (unit=41,file='fort.41',status='new',form='formatted')
      endif
      open (unit=51,file='fort.51',status='new',form='formatted')
   endif


   ! Populate the electron states to find the highest occupied state (Fermi
   !   energy for metals).
   call populateStates

!   if (doOPTC == 2) then ! Doing PACS calculation and we need to modify
!      ! the energy eigen values and their position. Now that the occupied
!      ! energy is known we can append the unoccupied excited states.
!      call appendExcitedEValsBand(lastInitStatePACS+1,numStates)
!   endif

   ! Shift the energy eigen values according to the highest occupied state.
   call shiftEnergyEigenValues(occupiedEnergy,numStates)

   ! Compute some statistics and variables concerning the energy values.
   call getEnergyStatistics(doOPTC)

   ! Compute the transition pairs and energy values of those transitions.
   call computeTransitions(inSCF,doOPTC)

   ! Print the output if necessary.
   if (doOPTC /= 3) then  ! Not doing a Sigma(E) calculation.
      call printOptcResults(doOPTC) ! Internally distinguishes optc, poptc.
   endif

   ! Deallocate unused matrices
   deallocate (energyEigenValues)
   deallocate (transCounter)
   if (doOPTC /= 3) then  ! Not doing a sigma(E) calculation.
      deallocate (energyDiff)
      if (detailCodePOPTC == 0) then
         deallocate (transitionProb)
      else
         deallocate (transitionProbPOPTC)
      endif
   endif

end subroutine optc


subroutine cleanUpSCF

   ! Use necessary modules.
   use O_CommandLine, only: doDIMO_SCF, doForce_SCF
   use O_Potential, only:  cleanUpPotential
   use O_AtomicSites, only: cleanUpAtomSites
   use O_AtomicTypes, only: cleanUpAtomTypes
   use O_PotSites, only: cleanUpPotSites
   use O_PotTypes, only: cleanUpPotTypes
   use O_SecularEquation, only: cleanUpSecularEqn
   use O_KPoints, only: cleanUpKPoints
   use O_ExchangeCorrelation, only: cleanUpExchCorr

   implicit none

   ! Close the output files
   close (7)
   close (8)
   close (13)
   close (14)
   close (20)

   ! Deallocate all the other as of yet un-deallocated arrays.
   call cleanUpAtomTypes
   call cleanUpAtomSites
   call cleanUpPotTypes
   call cleanUpPotSites
   call cleanUpKPoints
   call cleanUpExchCorr
   call cleanUpPotential

   ! Because we already deallocated in makeValence Rho we only deallocate
   !   if necessary.
   if ((doDIMO_SCF == 0) .and. (doForce_SCF == 0)) then
      call cleanUpSecularEqn
   endif

   ! Open file to signal completion of the program to the calling olcao script.
   open (unit=2,file='fort.2',status='unknown')

end subroutine cleanUpSCF


end module O_OLCAO
