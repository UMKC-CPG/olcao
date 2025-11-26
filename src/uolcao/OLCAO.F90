module O_OLCAO

contains

subroutine OLCAO

   ! Import necessary modules. Note that 
   use O_SCFHDF5
   use O_PSCFHDF5
   use O_CommandLine
   use O_TimeStamps, only: initOperationLabels
   use O_ElementData,     only: initElementData

   ! Initialize the logging labels.
   call initOperationLabels

   ! Initialize element data from periodic table of the elements.
   call initElementData

   ! Parse the command line parameters
   call parseCommandLine

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

      if (doDIMO_SCF == 1) then
         call dimo(1)  ! Passing inSCF == 1
      endif

      if (doOPTC_SCF >= 1) then
         call optc(1, doOPTC_SCF)  ! Passing inSCF == 1
      endif

      if (doField_SCF == 1) then
         call field(1)  ! Passing inSCF == 1
      endif

      if (doMTOP_SCF == 1) then
         call mtop(1)  ! Passing inSCF == 1
      endif

      if (doPSCF == 1) then
         !call reset
      endif

      call cleanUpSCF

      call closeHDF5_SCF
   endif

   if (doPSCF == 1) then
      call intgPSCF ! Preparations for wave function calculation.

      call bandPSCF ! Wave function calculation.

      if (doDOS_PSCF == 1) then
         call dos(0) ! Passing inSCF == 0
      endif

      if (doBond_PSCF >= 1) then
         call bond(0, doBond_PSCF) ! Passing inSCF == 0
      endif

      if (doDIMO_PSCF == 1) then
         call dimo(0) ! Passing inSCF == 0
      endif

      if (doOptc_PSCF >= 1) then
         call optc(0, doOPTC_PSCF) ! Passing inSCF == 0
      endif

      if (doField_PSCF == 1) then
         call field(0)  ! Passing inSCF == 0
      endif

      if (doMTOP_PSCF == 1) then
         call mtop(0)  ! Passing inSCF == 0
      endif

      call cleanUpPSCF

      call closeHDF5_PSCF
   endif

   if (doLoEn == 1) then
      call loen(0)
   endif

end subroutine OLCAO


subroutine setupSCF
   
   ! Import the precision variables
   use O_Kinds

   ! Import the necessary modules.
   use O_SCFHDF5, only: initHDF5_SCF
   use O_SCFIntegralsHDF5, only: atomOverlap_did, atomOverlapCV_did, &
         & atomKEOverlap_did, atomMVOverlap_did, atomNPOverlap_did, &
         & atomDMOverlap_did, atomMMOverlap_did, atomKOverlap_did, &
         & atomKOverlapPlusG_did, atomPotOverlap_did, atomOverlap_aid, &
         & atomKEOverlap_aid, atomMVOverlap_aid, atomNPOverlap_aid, &
         & atomDMOverlap_aid, atomMMOverlap_aid, atomKOverlap_aid, &
         & atomKOverlapPlusG_aid, atomPotTermOL_aid, &
         & numComponents, fullCVDims, packedVVDims
      use O_CommandLine, only: doDIMO_SCF, doOPTC_SCF, doField_SCF, doMTOP_SCF
   use O_Input, only: parseInput, numStates
   use O_Lattice, only: initializeLattice, initializeFindVec, &
         & cleanUpLattice, recipVectors
   use O_KPoints, only: numKPoints, initializeKPoints
   use O_Basis, only: renormalizeBasis, cleanUpBasis
   use O_ExchangeCorrelation, only: maxNumRayPoints, getECMeshParameters, &
         & makeECMeshAndOverlap, cleanUpExchCorr
   use O_ElectroStatics, only: makeElectrostatics
   use O_GaussianRelations, only: makeAlphaDist, makeAlphaNucDist, &
         & makeAlphaPotDist, cleanUpGaussRelations
   use O_Integrals,   only: allocateIntegralsSCF, gaussOverlapOL, &
         & gaussOverlapKE, gaussOverlapMV, gaussOverlapNP, &
         & elecPotGaussOverlap, cleanUpIntegrals, &
         & secondCleanUpIntegrals
   use O_Integrals3Terms ! Use all so we can exclude gaussKOverlap
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

   ! Define local variables.
   real (kind=double), dimension(3,3) :: zeroVectors

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

   ! Initialize local variables.
   zeroVectors(:,:) = 0.0_double

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
   call parseInput(1) ! inSCF=1


   ! Find specific computational parameters not EXPLICITLY given in the input
   !   file, but which can, however, be easily determined from the input file.
   call getImplicitInfo


   ! Create real-space and reciprocal-space super lattices out of the primitive
   !   lattice.  These "supercells" must be big enough so as to include all the
   !   points within a sphere bounded by the negligability limit.  Points
   !   outside the sphere are considered negligable.
   call initializeLattice (1)


   ! Setup the necessary data structures so that we can easily find the lattice
   !   vector that is closest to any other arbitrary vector.
   call initializeFindVec


   ! Compute the desired set of kpoints.
   call initializeKPoints (1) ! inSCF == 1


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


   ! Calculate the matrix elements of the overlap between all LCAO basis fns.
   call gaussOverlapOL(numComponents,fullCVDims,packedVVDims,atomOverlap_did,&
         & atomOverlapCV_did,atomOverlap_aid)


   ! Calculate the matrix elements of the kinetic energy between all LCAO
   !   basis functions.
   call gaussOverlapKE(packedVVDims,atomKEOverlap_did,atomKEOverlap_aid)


   ! Calculate the matrix elements of the mass velocity between all LCAO
   !   basis functions if needed for the scalar relativistic calculation.
   if (rel == 1) then
      call gaussOverlapMV(packedVVDims,atomMVOverlap_did,atomMVOverlap_aid)
   endif


   ! Create the alpha distance matrix with nuclear alpha factor
   call makeAlphaNucDist


   ! Calculate the matrix elements of the overlap between all LCAO
   !   wave functions and the nuclear potentials.
   call gaussOverlapNP(packedVVDims,atomNPOverlap_did,atomNPOverlap_aid)


   ! Create the alpha distance matrix with potential alpha factor
   call makeAlphaPotDist


   ! Calculate the matrix elements of the overlap between all LCAO
   !   wave functions and the potential site potential alphas.
   call elecPotGaussOverlap(packedVVDims,atomPotOverlap_did,atomPotTermOL_aid)


   ! Now that all the single matrices are done being made we can deallocate
   !   the data structures that were used in all the above subroutines but are
   !   not necessary now.
   call cleanUpIntegrals

   ! If any supplementary three-term matrices (xyz) were requested for
   !   computing other properties, then allocate the space that is necessary
   !   for the integral calculation and do the integral.
   if ((doDIMO_SCF == 1) .or. (doOPTC_SCF >= 1) .or. (doMTOP_SCF == 1)) then

      ! Consider in the future an option to do the XYZ independently to
      !   conserve memory if that becomes a problem.
      call allocateIntegrals3Terms(coreDim,valeDim,numKPoints)

      ! Do the relevant integrals.
      if (doDIMO_SCF == 1) then
         call gaussOverlapDM(packedVVDims,atomDMOverlap_did,atomDMOverlap_aid)
      endif
      if (doOPTC_SCF >= 1) then
         call gaussOverlapMM(packedVVDims,atomMMOverlap_did,atomMMOverlap_aid)
      endif
#ifndef GAMMA
      if (doMTOP_SCF == 1) then
         call gaussKOverlap(packedVVDims,atomKOverlap_did,atomKOverlap_aid,&
               & zeroVectors)
         call gaussKOverlap(packedVVDims,atomKOverlapPlusG_did,&
               & atomKOverlapPlusG_aid,recipVectors)
      endif
#endif

      ! Clean up from the 3Term integrals.
      call cleanUpIntegrals3Terms
   endif


   call secondCleanupIntegrals ! Overlap matrix parts used for ortho.
   call cleanUpGaussRelations
   if (doField_SCF /= 1) then
      call cleanUpBasis
   endif


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
   use O_SecularEquation, only: secularEqnSCF, cleanUpSecularEqn, &
         & shiftEnergyEigenValues
   use O_ValeCharge, only: makeValenceRho
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
         call secularEqnSCF(i,numStates)
      enddo


      ! Find the Fermi level, and the number of occupied bands, and the
      !   number of electrons occupying each of those bands.
      call populateStates


      ! Compute the TDOS if it was requested for each iteration.
      if (iterFlagTDOS == 1) then
         call computeIterationTDOS
      endif


      ! Calculate the valence charge density
      call makeValenceRho(1) ! inSCF == 1


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
      call secularEqnSCF(i, numStates)
   enddo

end subroutine mainSCF


subroutine intgPSCF

   ! Import the necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Input, only: parseInput, numStates
   use O_CommandLine, only: doDIMO_PSCF, doOPTC_PSCF, doSYBD_PSCF, &
         & doMTOP_PSCF, doField_PSCF
   use O_Potential, only: initPotCoeffs, spin
   use O_Basis, only: renormalizeBasis, cleanUpBasis
   use O_Integrals, only: allocateIntegralsPSCF, gaussOverlapOL,&
         & gaussOverlapHamPSCF, cleanUpIntegrals, secondCleanUpIntegrals, &
         & cleanUpHamIntegrals, reallocateIntegralsPSCF
   use O_Integrals3Terms ! Use all so we can exclude gaussKOverlap
   use O_PSCFHDF5, only: initHDF5_PSCF
   use O_PSCFIntegralsHDF5, only: atomOverlapPSCF_did, atomOverlapCV_PSCF_did,&
         & atomHamOverlapPSCF_did, atomDMOverlapPSCF_did,&
         & atomMMOverlapPSCF_did, atomKOverlapPSCF_did,&
         & atomKOverlapPlusGPSCF_did,atomOverlapPSCF_aid,&
         & atomHamOverlapPSCF_aid,atomDMOverlapPSCF_aid,atomMMOverlapPSCF_aid,&
         & atomKOverlapPSCF_aid, atomKOverlapPlusGPSCF_aid,&
         & numComponents,fullCVDimsPSCF,packedVVDimsPSCF
   use O_Lattice, only: initializeLattice, initializeFindVec, recipVectors
   use O_KPoints, only: numKPoints, makePathKPoints, initializeKPoints
   use O_GaussianRelations, only: makeAlphaDist, makeAlphaNucDist,&
         & makeAlphaPotDist, cleanUpGaussRelations
   use O_AtomicSites, only: coreDim, valeDim

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables.
   real (kind=double), dimension(3,3) :: zeroVectors

   ! Initialize local variables.
   zeroVectors(:,:) = 0.0_double

   ! Open the potential file that will be read from in this program.
   open (unit=8,file='fort.8',status='old',form='formatted')


   ! Read in the input to initialize all the key data structure variables.
   call parseInput(0) ! inSCF=0


   ! Find specific computational parameters not EXPLICITLY given in the input
   !   file, but which can, however, be easily determined from the input file.
   call getImplicitInfo


   ! Create real-space super lattice out of the primitive lattice.  These
   !   "supercells" must be big enough so as to include all the points within
   !   sphere bounded by the negligability limit.  Points outside the sphere
   !   are considered negligable and are ignored.
   call initializeLattice (1)


   ! Setup the necessary data structures so that we can easily find the lattice
   !   vector that is closest to any other arbitrary vector.
   call initializeFindVec


   ! Initialize the kpoints according to any input control parameters.
   call initializeKPoints(0) ! inSCF == 0


   ! Renormalize the basis functions.
   call renormalizeBasis


   ! Prepare the HDF5 files for the post SCF calculations.
   call initHDF5_PSCF(numStates)


   ! Create the alpha distance matrices.
   call makeAlphaDist


   ! Allocate space for the integral results.
   call allocateIntegralsPSCF(coreDim,ValeDim,numKPoints)


   ! Read the provided potential coefficients.
   call initPotCoeffs


   ! Calculate the matrix elements of the overlap between all LCAO basis fns.
   call gaussOverlapOL(numComponents,fullCVDimsPSCF,packedVVDimsPSCF,&
         & atomOverlapPSCF_did,atomOverlapCV_PSCF_did,atomOverlapPSCF_aid)


   ! Create the alpha distance matrix with nuclear alpha factor
   call makeAlphaNucDist


   ! Create the alpha distance matrix with potential alpha factor
   call makeAlphaPotDist


   ! Reallocate integrals from the overlap to the hamiltonian.
   call reallocateIntegralsPSCF(coreDim,valeDim,numKPoints,spin)


   ! Compute the hamiltonian matrix elements.
   call gaussOverlapHamPSCF(atomHamOverlapPSCF_did,atomHamOverlapPSCF_aid)


   ! Now that all the single matrices are done being made we can deallocate
   !   the data structures that were used in all the above subroutines but are
   !   not necessary now.
   call cleanUpHamIntegrals


   ! If any supplementary three-term matrices (xyz) were requested for
   !   computing other properties, then allocate the space that is necessary
   !   for the integral calculation and do the integral.
   if ((doDIMO_PSCF == 1) .or. (doOPTC_PSCF >= 1) .or. (doMTOP_PSCF == 1)) then

      ! Consider in the future an option to do the XYZ independently to
      !   conserve memory if that becomes a problem.
      call allocateIntegrals3Terms(coreDim,valeDim,numKPoints)

      ! Do the relevant integrals.
      if (doDIMO_PSCF == 1) then
         call gaussOverlapDM(packedVVDimsPSCF,atomDMOverlapPSCF_did,&
               & atomDMOverlapPSCF_aid)
      endif
      if (doOPTC_PSCF >= 1) then
         call gaussOverlapMM(packedVVDimsPSCF,atomMMOverlapPSCF_did,&
               & atomMMOverlapPSCF_aid)
      endif
#ifndef GAMMA
      if (doMTOP_PSCF == 1) then
         call gaussKOverlap(packedVVDimsPSCF,atomKOverlapPSCF_did,&
               & atomKOverlapPSCF_aid,zeroVectors)
         call gaussKOverlap(packedVVDimsPSCF,&
               & atomKOverlapPlusGPSCF_did,atomKOverlapPlusGPSCF_aid,&
               & recipVectors)
      endif
#endif

      ! Clean up from the 3Term integrals.
      call cleanUpIntegrals3Terms
   endif


   call secondCleanupIntegrals ! Overlap matrix parts used for ortho.
   call cleanUpGaussRelations
   if (doField_PSCF /= 1) then
      call cleanUpBasis
   endif


end subroutine intgPSCF


subroutine bandPSCF

      ! Import the necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Input,            only: numStates
   use O_Potential,        only: spin
   use O_AtomicSites,      only: valeDim
   use O_LAPACKParameters, only: setBlockSize
   use O_Input,            only: numStates
   use O_CommandLine,      only: doSYBD_PSCF, doMTOP_PSCF
   use O_Lattice,          only: initializeLattice, initializeFindVec
   use O_KPoints,          only: makePathKPoints, numKPoints, &
         & computePhaseFactors
   use O_PSCFIntegralsHDF5, only: numComponents, atomOverlapPSCF_did,&
         & atomHamOverlapPSCF_did
   use O_PSCFEigVecHDF5, only: eigenVectorsPSCF_did, eigenVectorsPSCF_aid
   use O_PSCFEigValHDF5, only: eigenValuesPSCF_did
   use O_SecularEquation, only: secularEqnPSCF, cleanUpSecularEqn

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define local variables.
   integer :: i


   ! Set up LAPACK machine parameters.
   call setBlockSize(valeDim)

   ! Compute the solid state wave function for each kpoint.

   ! Loop over all kpoints and compute the solid state wave function and
   !   energy eigen values for each. But first, ...
   ! Record the fact that we are starting the k-point loop so that when
   !   someone looks at the output file as the job is running they will
   !   know that all the little dots represent a count of the number of
   !   k-points. Then they can figure out the progress and progress rate.
   write (20,*) "Beginning k-point loop."
   write (20,*) "Expecting ",numKPoints," iterations."
   call flush (20)

   do i = 1, spin
      call secularEqnPSCF(i,numStates,numComponents,atomOverlapPSCF_did,&
            & atomHamOverlapPSCF_did,eigenValuesPSCF_did,&
            & eigenVectorsPSCF_did,eigenVectorsPSCF_aid)
   enddo

!   call cleanUpSecularEqn

   ! Print the band results if necessary.
   if (doSYBD_PSCF == 1) then
      call printSYBD
   endif

end subroutine bandPSCF


subroutine printSYBD

   ! Use necessary modules.
   use O_Potential,       only: spin
   use O_Constants,       only: hartree
   use O_Input,           only: numStates
   use O_KPoints,         only: numPathKP, pathKPointMag
   use O_Populate,        only: occupiedEnergy, populateStates
   use O_SecularEquation, only: energyEigenValues, shiftEnergyEigenValues
   use O_AtomicSites
   use O_AtomicTypes
   use O_PotTypes

   ! Define local variables.
   integer :: i,j,k,l
   integer :: basisFnCount
   integer :: currAtomType
   character*7 :: filename

   ! Obtain the fermi energy with the populate energy levels subroutine.
   call populateStates

   ! Adjust the energyEigenValues down by the fermi level determined above.
   call shiftEnergyEigenValues(occupiedEnergy)

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
         write (30+i,fmt="(11f14.5)") &
               & energyEigenValues(:numStates,j,i)*hartree
      enddo

      ! Close the output file.
      close (30+i)

      ! Create the filename.
      write (filename,fmt="(a5,i2)") "fort.",32+i

      ! Open a file that will contain the valence dimension meta-data.
      open(unit=32+i,file=filename,status='new',form='formatted')

      ! Print a metadata string for each n,l,m basis function.
      basisFnCount = 0
      do j = 1, numAtomSites
         currAtomType = atomSites(j)%atomTypeAssn
         do k = 1, atomTypes(currAtomType)%numValeRadialFns
            do l = -atomTypes(currAtomType)%valeQN_lList(k), &
                  & atomTypes(currAtomType)%valeQN_lList(k)
               basisFnCount = basisFnCount + 1
               write (32+i,fmt="(i5)",advance="NO") basisFnCount
               write (32+i,fmt="(i5)",advance="NO") j
               write (32+i,fmt="(a,a3)",advance="NO") " ", &
                     & atomTypes(currAtomType)%elementName
               write (32+i,fmt="(1e10.3)",advance="NO") &
                     & potTypes(currAtomType)%nucCharge
               write (32+i,fmt="(i5)",advance="NO") &
                     & atomTypes(currAtomType)%elementID
               write (32+i,fmt="(i5)",advance="NO") &
                     & atomTypes(currAtomType)%speciesID
               write (32+i,fmt="(i5)",advance="NO") &
                     & atomTypes(currAtomType)%typeID
               write (32+i,fmt="(i5)",advance="NO") &
                     & atomTypes(currAtomType)%valeQN_nList(k)
               write (32+i,fmt="(i5)",advance="NO") &
                     & atomTypes(currAtomType)%valeQN_lList(k)
               write (32+i,fmt="(i5)",advance="NO") l
               write (32+i,fmt="(3e15.6)") &
                     & atomSites(j)%cartPos(1:3)
               
            enddo
         enddo
      enddo

      ! Close the meta data file.
      close(32+i)
   enddo

end subroutine printSYBD


subroutine dos(inSCF)

   ! Use necessary modules.
   use O_TimeStamps
   use O_Potential,       only: spin
   use O_DOS,             only: computeDOS
   use O_KPoints,         only: numKPoints
   use O_Input,           only: numStates
   use O_Populate,        only: occupiedEnergy, populateStates
   use O_SecularEquation, only: energyEigenValues, &
         & shiftEnergyEigenValues


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
   call shiftEnergyEigenValues(occupiedEnergy)

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
   use O_Input,           only: thermalSigma, numStates
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
   call shiftEnergyEigenValues(occupiedEnergy)

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


subroutine dimo(inSCF)

   ! Use necessary modules.
   use O_TimeStamps
   use O_Potential,       only: spin
   use O_KPoints,         only: numKPoints
   use O_Input,           only: numStates
   use O_Populate,        only: occupiedEnergy, populateStates
   use O_SecularEquation, only: energyEigenValues, &
         & shiftEnergyEigenValues
   use O_ValeCharge, only: makeValenceRho


   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare passed parameters.
   integer, intent(in) :: inSCF

   ! Open the DOS files that will be written to.  If a spin polarized
   !   calculation is being done, then 60, 70, 80 hold spin up and 61, 71, 81
   !   hold spin down.  60,61=TDOS; 70,71= PDOS; 80,81=Localization Index
   open (unit=74,file='fort.74',status='new',form='formatted')
   if (spin == 2) then
      open (unit=75,file='fort.75',status='new',form='formatted')
   endif


   ! Populate the electron states to find the highest occupied state (Fermi
   !   energy for metals).
   call populateStates

   ! Shift the energy eigen values according to the highest occupied state.
   call shiftEnergyEigenValues(occupiedEnergy)

   ! Call the DOS subroutine to compute the total and partial density of states
   !   as well as the localization index.
   call makeValenceRho(inSCF)

   ! Close the output files.
   close(74)
   if (spin == 2) then
      close(75)
   endif

end subroutine dimo


subroutine field(inSCF)

! The goal of this program is to produce plottable data of the solid state
!   wave function, data derived form it (e.g. the charge density), and/or
!   ancillary data such as the potential function. It permits plotting of
!   specific wave function states according to energy level by numerically
!   evaluating the analytical wave functions, wave functions squared, etc.
!   on a defined 3D uniform mesh. The data may be stored in either plain
!   text OpenDX format or in XDMF+HDF5 format for Paraview.

   ! Import the necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Populate,        only: populateStates
   use O_Potential,       only: spin, initPotCoeffs
   use O_Field,           only: computeFieldMesh, cleanUpField
!   use O_SCFFieldHDF5,    only: wav_did, rho_did, pot_did, triggerAxis, &
!         & abcDimsChunk, fileFieldChunk_dsid
!   use O_PSCFFieldHDF5,    only: wavPSCF_did, rhoPSCF_did, potPSCF_did, &
!         & triggerAxisPSCF, abcDimsChunkPSCF, fileFieldChunkPSCF_dsid
   use O_SecularEquation, only: energyEigenValues
   use O_Lattice,         only: initialize3DMesh
   use O_FieldHDF5,       only: prepFieldHDF5, closeFieldHDF5

   ! Make sure that there are not accidental variable declarations.
   implicit none


   ! Define passed parameters.
   integer, intent(in) :: inSCF


   ! Initialize (or access) the HDF field file.
   call prepFieldHDF5(inSCF)


   ! Open the potential file that will be read from in this program.
   open (unit=8,file='fort.8',status='old',form='formatted')


   ! Populate the electron states to find the highest occupied state (Fermi
   !   energ for metals).
   call populateStates


   ! Initialize certain parameters for constructing and traversing the 3D mesh.
   call initialize3DMesh


   ! Compute requested field values for each mesh point and store in HDF5.
   call computeFieldMesh(inSCF)


   ! Clean up any left over arrays that need to be deallocated.
   call cleanUpField


   ! Close the field HDF5 file.
   call closeFieldHDF5


end subroutine field


subroutine optc(inSCF,doOPTC)

   ! Use necessary modules.
   use O_TimeStamps
   use O_Potential,       only: spin
   use O_OptcPrint,       only: printOptcResults
   use O_KPoints,         only: numKPoints
   use O_Populate,        only: occupiedEnergy, populateStates
   use O_Lattice,         only: initializeLattice, initializeFindVec
   use O_Input,           only: numStates, lastInitStatePACS, &
         & detailCodePOPTC
   use O_OptcTransitions, only: transCounter, energyDiff, transitionProb, &
         & transitionProbPOPTC, getEnergyStatistics, computeTransitions
   use O_SecularEquation, only: energyEigenValues, &
         & shiftEnergyEigenValues


   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: inSCF
   integer, intent(in) :: doOPTC


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
   call shiftEnergyEigenValues(occupiedEnergy)

   ! Compute some statistics and variables concerning the energy values.
   call getEnergyStatistics(doOPTC)

   ! Compute the transition pairs and energy values of those transitions.
   call computeTransitions(inSCF,doOPTC)

   ! Print the output if necessary.
   if (doOPTC /= 3) then  ! Not doing a Sigma(E) calculation.
      call printOptcResults(doOPTC) ! Internally distinguishes optc, poptc.
   endif

   ! Deallocate unused matrices
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


subroutine mtop(inSCF)

   ! Use necessary modules.
   use HDF5
   use O_TimeStamps
   use O_Potential,       only: spin
   use O_Populate,        only: occupiedEnergy, populateStates
   use O_SecularEquation, only: shiftEnergyEigenValues
   use O_MTOP ! Use all to avoid def GAMMA issues.
      
   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare passed parameters.
   integer, intent(in) :: inSCF

   ! Declare local variables.
   integer :: i
   real(kind=double), dimension(3,2) :: xyzP

   ! Open the MTOP files that will be written to.  If a spin polarized
   !   calculation is being done, then 180 holds spin up and 181 holds
   !   spin down.
   open (unit=180,file='fort.180',status='new',form='formatted')
   if (spin == 2) then
      open (unit=181,file='fort.181',status='new',form='formatted')
   endif


   ! Populate the electron states to find the highest occupied state (Fermi
   !   energy for metals).
   call populateStates

   ! Shift the energy eigen values according to the highest occupied state.
   call shiftEnergyEigenValues(occupiedEnergy)
#ifndef GAMMA
   call computeMTOPPolarization(inSCF,xyzP)

   !! Print out the MTOP result.
   !do i = 1, spin
   !   write(20,*) 'Polarization (a.u.):', xyzP(:,i)

   !   ! Print to stdout as requested
   !   if (i == 1) then
   !      write(20,fmt="(a)") "Spin Up"
   !      write(20,'(A,1X,F20.12)') 'P(x) =', xyzP(1,i)
   !      write(20,'(A,1X,F20.12)') 'P(y) =', xyzP(2,i)
   !      write(20,'(A,1X,F20.12)') 'P(z) =', xyzP(3,i)
   !   else
   !      write(20,fmt="(a)") "Spin Dn"
   !      write(20,'(A,1X,F20.12)') 'P(x) =', xyzP(1,i)
   !      write(20,'(A,1X,F20.12)') 'P(y) =', xyzP(2,i)
   !      write(20,'(A,1X,F20.12)') 'P(z) =', xyzP(3,i)
   !   endif

   !   ! Also write to fort.180 so the run script can copy it
   !   write(179+i,fmt="(3e16.8)") xyzP(:,i)
   !   flush(179+i)
   !enddo
#endif

   ! Close the output files.
   close(180)
   if (spin == 2) then
      close(181)
   endif

end subroutine mtop


subroutine loen(inSCF)

   use O_Kinds

   ! Import necessary module subroutines and data.
   use O_LocalEnv
   use O_Input,       only: parseInput, loenCode
   use O_Lattice,     only: initializeLattice, initializeFindVec,&
                          & cleanUpLattice
   use O_TimeStamps,  only: initOperationLabels

   ! Import the HDF5 module.
   use HDF5

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: inSCF


   ! Read in the input to initialize all the key data structure variables.
   call parseInput(inSCF)


   ! Find specific computational parameters not EXPLICITLY given in the input
   !   file, but which can, however, be easily determined from the input file.
   call getImplicitInfo


   ! Create real-space and reciprocal-space super lattices out of the primitive
   !   lattice.  These "supercells" must be big enough so as to include all the
   !   points within a sphere bounded by the negligability limit.  Points
   !   outside the sphere are considered negligable.
   call initializeLattice (1)


   ! Setup the necessary data structures so that we can easily find the lattice
   !   vector that is closest to any other arbitrary vector.
   call initializeFindVec


   ! Compute the requested local environment metric.
   if (loenCode == 1) then
      call bispec ! Use the bispectrum component method.
   endif

   ! Close the output file
   close (20)

   ! Open a file to signal completion of the program.
   open (unit=2,file='fort.2',status='unknown')
   

end subroutine loen


subroutine cleanUpSCF

   ! Use necessary modules.
   use O_CommandLine, only: doDIMO_SCF, doForce_SCF, doPSCF
   use O_Lattice, only: cleanUpLattice
   use O_Potential, only:  cleanUpPotential
   use O_AtomicSites, only: cleanUpAtomSites
   use O_AtomicTypes, only: cleanUpAtomTypes
   use O_PotSites, only: cleanUpPotSites
   use O_PotTypes, only: cleanUpPotTypes
   use O_SecularEquation, only: cleanUpSecularEqn
   use O_KPoints, only: cleanUpKPoints
   use O_ExchangeCorrelation, only: cleanUpExchCorr

   implicit none

   ! Close the output files.
   close (7) ! Iteration data
   close (8) ! SCF Potential
   close (13) ! Magnetic moments
   close (14) ! Energy data per iteration
   if(doPSCF < 0) then
      close (20) ! Primary output
   endif

   ! Deallocate all the other as of yet un-deallocated arrays.
   call cleanUpAtomTypes
   call cleanUpAtomSites
   call cleanUpPotTypes
   call cleanUpPotSites
   call cleanUpKPoints
   call cleanUpExchCorr
   call cleanUpPotential
   call cleanUpLattice

   ! Because we already deallocated in makeValence Rho we only deallocate
   !   if necessary.
   if ((doDIMO_SCF < 0) .and. (doForce_SCF < 0)) then
      call cleanUpSecularEqn
   endif

   ! Open file to signal completion of the program to the calling olcao script.
   open (unit=2,file='fort.2',status='unknown')

end subroutine cleanUpSCF


subroutine cleanUpPSCF

   ! Use necessary modules.
   use O_CommandLine, only: doSYBD_PSCF, doMTOP_PSCF
   use O_Lattice, only: cleanUpLattice
   use O_Potential, only:  cleanUpPotential
   use O_AtomicSites, only: cleanUpAtomSites
   use O_AtomicTypes, only: cleanUpAtomTypes
   use O_PotSites, only: cleanUpPotSites
   use O_PotTypes, only: cleanUpPotTypes
   use O_SecularEquation, only: cleanUpSecularEqn
   use O_KPoints, only: cleanUpKPoints

   implicit none

   ! Close any opened files.
   close (8)
   close (20)

   ! Deallocate all the other as of yet un-deallocated arrays.
   call cleanUpAtomTypes
   call cleanUpAtomSites
   call cleanUpPotTypes
   call cleanUpPotSites
   call cleanUpKPoints
   call cleanUpPotential
   call cleanUpLattice

   ! Deallocate if necessary.
   if (doSYBD_PSCF < 0) then
      call cleanUpSecularEqn
   endif

   ! Open file to signal completion of the program to the calling olcao script.
   open (unit=2,file='fort.2',status='unknown')

end subroutine cleanUpPSCF


end module O_OLCAO
