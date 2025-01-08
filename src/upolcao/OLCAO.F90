module O_OLCAO

contains

subroutine OLCAO

   ! Import necessary modules.
   use O_MPI
   use O_SCFHDF5
   use O_PSCFHDF5
   use O_CommandLine
   use O_TimeStamps, only: initOperationLabels
   use O_LocalEnv

   ! Initialize the MPI interface and process variables.
   call initMPI
!write(20+mpiRank,*) "Hello from process ", mpiRank, "out of ", mpiSize
!call flush(20+mpiRank)

   ! Only rank 0 needs to initialize the logging labels.
   call initOperationLabels

   ! All ranks parse the command line parameters
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

      call cleanUpPSCF

      call closeHDF5_PSCF
   endif

   if (doLoEn == 1) then
      call analyzeLocalEnv
   endif

   ! Finalize the MPI environment
   call closeMPI

end subroutine OLCAO


subroutine setupSCF
   
   ! Import the precision variables
   use O_Kinds

   ! Import the necessary modules.
   use O_SCFHDF5, only: initHDF5_SCF
   use O_SCFIntegralsHDF5, only: atomOverlap_did, atomOverlapCV_did, &
         & atomKEOverlap_did, atomMVOverlap_did, atomNPOverlap_did, &
         & atomDMOverlap_did, atomMMOverlap_did, atomPotOverlap_did, &
         & atomOverlap_aid, atomKEOverlap_aid, atomMVOverlap_aid, &
         & atomNPOverlap_aid, atomDMOverlap_aid, atomMMOverlap_aid, &
         & atomPotTermOL_aid, numComponents, fullCVDims, packedVVDims
      use O_CommandLine, only: doDIMO_SCF, doOPTC_SCF, doField_SCF
   use O_Input, only: parseInput, numStates
   use O_Lattice, only: initializeLattice, initializeFindVec, &
         & cleanUpLattice
   use O_KPoints, only: numKPoints, computePhaseFactors
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
   use O_Integrals3Terms,  only: allocateIntegrals3Terms, &
         & gaussOverlapDM, gaussOverlapMM, cleanUpIntegrals3Terms
   use O_AtomicSites, only: coreDim, valeDim, cleanUpAtomSites
   use O_AtomicTypes, only: cleanUpRadialFns, cleanUpAtomTypes
   use O_PotSites, only: cleanUpPotSites
   use O_PotTypes, only: cleanUpPotTypes
   use O_Potential, only: rel, cleanUpPotential
   use O_CoreCharge, only: makeCoreRho

   ! Import the HDF5 module.
   use HDF5

   ! Import the necessary OLCAO MPI data.
   use O_MPI

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! All ranks read in the input and initialize all key data structure
   !   variables.
   call parseInput(1) ! inSCF=1


   ! All ranks find specific computational parameters not EXPLICITLY given in
   !   the input file, but which can, however, be easily determined from the
   !   input file.
   call getImplicitInfo


   ! All ranks create real-space and reciprocal-space super lattices out of
   !   the primitive lattice.  These "super lattices" must be big enough so
   !   as to include all the points within a sphere bounded by the
   !   negligability limit.  Points outside the sphere are considered
   !   negligable.
   call initializeLattice (1)


   ! All ranks setup the necessary data structures so that we can easily find
   !   the lattice vector that is closest to any other arbitrary vector.
   call initializeFindVec


   ! All ranks compute the kpoint phase factors. FIX so that we parallelize
   !   over kpoints. Ranks compute phase factors of only the kpoints that
   !   are relevant to them.
   call computePhaseFactors


   ! All ranks renormalize the basis functions.
   call renormalizeBasis


   ! All ranks determine the parameters for the exchange correlation mesh.
   !   FIX so that all ranks work together to determine the parameters for
   !   the exchange correlation mesh.
   call getECMeshParameters

   ! Now, the dimensions of the system are known.  Therefore, we can
   !   initialize the HDF5 file structure format, and datasets.
   call initHDF5_SCF (maxNumRayPoints, numStates)


   ! All ranks construct the exchange correlation overlap matrix, and
   !   sampling field. FIX so that all ranks work together to compute this.
   !   Note that only rank 0 will actually save results to disk.
   call makeECMeshAndOverlap


   ! All ranks create the alpha distance matrices.
   call makeAlphaDist


   ! All ranks allocate space to be used for each of the single matrix
   !   integrals. FIX so that each rank only allocates the space that it
   !   actually needs.
   call allocateIntegralsSCF(coreDim,valeDim,numKPoints)


   ! All ranks calculate the matrix elements of the overlap between all LCAO
   !   basis fns. FIX so that each rank only computes the portion of this
   !   matrix that they are responsible for. Note that only rank 0 will
   !   actually save results to disk.
   call gaussOverlapOL(numComponents,fullCVDims,packedVVDims,atomOverlap_did,&
         & atomOverlapCV_did,atomOverlap_aid)


   ! All ranks calculate the matrix elements of the kinetic energy between
   !   all LCAO basis functions. FIX so that each rank only computes the
   !   portion of this matrix that they are reponsible for. Note that only
   !   rank 0 will actually save results to disk.
   call gaussOverlapKE(packedVVDims,atomKEOverlap_did,atomKEOverlap_aid)


   ! All ranks calculate the matrix elements of the mass velocity between all
   !   LCAO basis functions if needed for the scalar relativistic calculation.
   !   FIX so that each rank only computes the portion of this matrix that
   !   they are reponsible for. Note that only rank 0 will actually save
   !   results to disk.
   if (rel == 1) then
      call gaussOverlapMV(packedVVDims,atomMVOverlap_did,atomMVOverlap_aid)
   endif


   ! All ranks create the alpha distance matrix with nuclear alpha factor.
   !   FIX if possible so that each rank only computes the portion of this
   !   data structure that it actually needs.
   call makeAlphaNucDist


   ! All ranks calculate the matrix elements of the overlap between all LCAO
   !   wave functions and the nuclear potentials. FIX so that each rank only
   !   computes the portion of this matrix that they are reponsible for.
   !   Note that only rank 0 will actually save results to disk.
   call gaussOverlapNP(packedVVDims,atomNPOverlap_did,atomNPOverlap_aid)


   ! All ranks create the alpha distance matrix with potential alpha factor.
   !   FIX if possible so that each rank only computes the portion of this
   !   data structure that it actually needs.
   call makeAlphaPotDist


   ! All ranks calculate the matrix elements of the overlap between all LCAO
   !   wave functions and the potential site potential alphas. FIX so that
   !   each rank only computes the portion of this matrix that they are
   !   reponsible for. Note that only rank 0 will actually save results to
   !   disk.
   call elecPotGaussOverlap(packedVVDims,atomPotOverlap_did,atomPotTermOL_aid)


   ! Now that all the single matrices are done being made all ranks can
   !   deallocate the data structures that were used in all the above
   !   subroutines but that are not necessary now.
   call cleanUpIntegrals


   ! If any supplementary three-term matrices (xyz) were requested for
   !   computing other properties, then allocate the space that is necessary
   !   for the integral calculation and do the integral.
   if ((doDIMO_SCF == 1) .or. (doOPTC_SCF >= 1)) then

      ! All ranks allocate memory for the three-term integrals. Consider in
      !   the future an option to do the XYZ independently to conserve memory
      !   if that becomes a problem. First though, FIX so that each rank only
      !   allocates memory for the portion of the matrix that it is solving.
      call allocateIntegrals3Terms(coreDim,valeDim,numKPoints)

      ! All ranks do the relevant integrals. FIX so that each rank only
      !   computes the portion of this matrix that they are reponsible for.
      !   Note that only rank 0 will actually save results to disk.
      if (doDIMO_SCF == 1) then
         call gaussOverlapDM(packedVVDims,atomDMOverlap_did,atomDMOverlap_aid)
      endif
      if (doOPTC_SCF >= 1) then
         call gaussOverlapMM(packedVVDims,atomMMOverlap_did,atomMMOverlap_aid)
      endif

      ! All ranks clean up from the 3Term integrals.
      call cleanUpIntegrals3Terms
   endif


   ! All ranks do further clean up.
   call secondCleanupIntegrals ! Overlap matrix parts used for ortho.
   call cleanUpGaussRelations
   if (doField_SCF /= 1) then
      call cleanUpBasis
   endif


   ! All ranks construct a vector describing the core charge density since it
   !   will not change throughout the SCF iterations because of core
   !   orthogonalization. FIX so that each rank computes a portion of the
   !   vector. Note that only rank 0 saves the results to disk.
   if (coreDim /= 0) then
      call makeCoreRho
   endif


   ! All ranks deallocate the data component of the atomType data structure
   !   that holds the radial functions since they are no longer needed.
   call cleanUpRadialFns


   ! All ranks construct matrix operators and integral vectors that are used
   !   to later determine the electrostatic potential. FIX so that each rank
   !   only computes a portion of the operators and integral vectors. Note
   !   that only rank 0 saves results to disk.
   call makeElectrostatics

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
   use MPI_F08
   use O_MPI
   use O_Kinds
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
   use O_Populate, only: populateStates
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
   real (kind=double) :: totalEnergy

   ! Use rank 0 to open (almost) all the text files that will be written to.
   if (mpiRank == 0) then
      open (unit=7,file='fort.7',status='unknown',form='formatted')
      open (unit=8,file='fort.8',status='old',form='formatted')
      if (spin == 2) then
         open (unit=13,file='fort.13',status='unknown',form='formatted')
      endif
      open (unit=14,file='fort.14',status='unknown',form='formatted')
   endif

   ! Use rank 0 to label the columns that record statistics for each
   !   iteration.  Column 1 is the iteration number, column 2 is the electron
   !   occupation energy, column 3 is the electron fitting error when mapping
   !   the charge density into the auxiliary functional form, column 4 is
   !   the convergence criterion, column 5 is the total energy, and if the
   !   calculation is spin polarized, then the last column (6) is the magnetic
   !   moment.
   if (mpiRank == 0) then
      if (spin == 1) then
         write (7,*) "Iter. Occ. Energy   El. Error    ", &
               & "Convergence     Total Energy"
      else
         write (7,*) "Iter. Occ. Energy   El. Error    ", &
               & "Convergence     Total Energy      Mag. Mom."
      endif
      call flush (7)
   endif

   ! Use rank 0 to label the columns that record energy statistics for each
   !   iteration. Column 1 is the iteration number, column 2 is the kinetic
   !   energy, if the calculation is relativistic, then column 3 is the
   !   mass velocity and otherwise it is the electrostatic energy,
   !   the next column (4 or 5) is the exchange correlation energy, the
   !   next column (5 or 6) is the total energy.
   if (mpiRank == 0) then
      if (rel == 0) then
         write (14,*) "Iter. Kinetic        ElecStat       ", &
               & "ExchCorr       Total"
      else
         write (14,*) "Iter. Kinetic        MassVel        ", &
               & "ElecStat       ExchCorr       Total"
      endif
      call flush (14)
   endif

   ! All ranks obtain the initial potential coefficients.
   call initPotCoeffs


   ! All ranks set up LAPACK machine parameters. FIX so that all ranks set
   !   up the ScaLAPACK machine parameters.
   call setBlockSize(valeDim)


   ! All ranks begin the self consistance iterations.
   do while (.true.)

      ! All ranks solve the Kohn-Sham equation. FIX so that all ranks
      !   work together to solve the Kohn-Sham equation.
      do i = 1, spin
         call secularEqnSCF(i,numStates)
      enddo


      ! All ranks find the Fermi level, the number of occupied bands, and the
      !   number of electrons occupying each of those bands.
      call populateStates


      ! All ranks compute the TDOS if it was requested for each iteration.
      !   FIX so that all ranks work together to compute the TDOS.
      if (iterFlagTDOS == 1) then
         call computeIterationTDOS
      endif


      ! All ranks calculate the valence charge density.
      call makeValenceRho(1) ! inSCF == 1


      ! All ranks compute the new self consistant potential.
      call makeSCFPot(totalEnergy)


      ! Allranks determine if the computation is complete.
      if ((converged == 1) .or. (currIteration > lastIteration)) exit

      ! All ranks check if the job was requested to stop.
      if (mpiRank == 0) then
         open (unit=666,file="OLCAOkill",status="OLD",IOSTAT=OLCAOkill)
      endif
      call MPI_BCAST(OLCAOkill,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (OLCAOkill == 0) then
         exit
      endif
   enddo

   ! Only rank 0 prints the accumulated TDOS in a useful format if requested
   !   in the input.
   if (mpiRank == 0) then
      if (iterFlagTDOS == 1) then
         call printIterationTDOS
      endif
   endif

   ! All ranks compute the final band structure if convergence was achieved
   !   or if the last iteration was done.  (I.e. we are out of the SCF loop.)
   !   FIX so that all ranks work together to do this (as above).
   do i = 1, spin
      call secularEqnSCF(i, numStates)
   enddo

end subroutine mainSCF


subroutine intgPSCF

   ! Import the necessary modules.
   use O_Kinds
   use O_MPI
   use O_Input, only: parseInput, numStates
   use O_CommandLine, only: doDIMO_PSCF, doOPTC_PSCF, doSYBD_PSCF, doField_PSCF
   use O_Potential, only: initPotCoeffs, spin
   use O_Basis, only: renormalizeBasis, cleanUpBasis
   use O_Integrals, only: allocateIntegralsPSCF, gaussOverlapOL,&
         & gaussOverlapHamPSCF, cleanUpIntegrals, secondCleanUpIntegrals, &
         & cleanUpHamIntegrals, reallocateIntegralsPSCF
   use O_Integrals3Terms, only: gaussOverlapDM, gaussOverlapMM, &
         & allocateIntegrals3Terms, cleanUpIntegrals3Terms
   use O_PSCFHDF5, only: initHDF5_PSCF
   use O_PSCFIntegralsHDF5, only: atomOverlapPSCF_did, atomOverlapCV_PSCF_did,&
         & atomHamOverlapPSCF_did, atomDMOverlapPSCF_did,&
         & atomMMOverlapPSCF_did, atomOverlapSYBD_PSCF_did,&
         & atomOverlapCV_SYBD_PSCF_did, atomHamOverlapSYBD_PSCF_did,&
         & atomOverlapPSCF_aid,atomHamOverlapPSCF_aid,&
         & atomDMOverlapPSCF_aid,atomMMOverlapPSCF_aid,&
         & atomOverlapSYBD_PSCF_aid,atomHamOverlapSYBD_PSCF_aid,&
         & numComponents,fullCVDimsPSCF,packedVVDimsPSCF
   use O_Lattice, only: initializeLattice, initializeFindVec
   use O_KPoints, only: numKPoints, makePathKPoints, computePhaseFactors
   use O_GaussianRelations, only: makeAlphaDist, makeAlphaNucDist,&
         & makeAlphaPotDist, cleanUpGaussRelations
   use O_AtomicSites, only: coreDim, valeDim

   ! Make sure that there are not accidental variable declarations.
   implicit none


   ! Only rank 0 opens the potential file that will be read from in this
   !   program.
   if (mpiRank == 0) then
      open (unit=8,file='fort.8',status='old',form='formatted')
   endif


   ! All ranks read in the input to initialize all the key data structure
   !   variables.
   call parseInput(0) ! inSCF=0


   ! All ranks find specific computational parameters not EXPLICITLY given in
   !   the input file, but which can, however, be easily determined from the
   !   input file.
   call getImplicitInfo


   ! In the case of a symmetric band structure calculation all ranks must
   !   modify the kpoints according to the lattice type information given in
   !   the olcao.dat input file.
   if (doSYBD_PSCF == 1) then
      call makePathKPoints
   endif

   ! All ranks create real-space super lattice out of the primitive lattice.
   !   These "supercells" must be big enough so as to include all the points
   !   within sphere bounded by the negligability limit.  Points outside the
   !   sphere are considered negligable and are ignored.
   call initializeLattice (1)


   ! All ranks setup the necessary data structures so that we can easily find
   !   the lattice vector that is closest to any other arbitrary vector.
   call initializeFindVec


   ! All ranks compute the kpoint phase factors. FIX so that each rank only
   !   computes the phase factors for the kpoints it needs. (I.e., permit
   !   parallelization over kpoints.)
   call computePhaseFactors


   ! All ranks renormalize the basis functions.
   call renormalizeBasis


   ! Only rank 0 prepares the HDF5 files for the post SCF calculations.
   call initHDF5_PSCF(numStates)


   ! All ranks create the alpha distance matrices.
   call makeAlphaDist


   ! All ranks allocate space for the integral results. FIX so that each
   !   rank only allocates space for the portion of the integral matrices
   !   that they compute.
   call allocateIntegralsPSCF(coreDim,ValeDim,numKPoints)


   ! All ranks read the provided potential coefficients.
   call initPotCoeffs


   ! All ranks calculate the matrix elements of the overlap between all LCAO
   !   basis fns. FIX so that all ranks work together to compute the matrix.
   !   Note that only rank 0 will actually record results to disk.
   if (doSYBD_PSCF < 0) then
      call gaussOverlapOL(numComponents,fullCVDimsPSCF,packedVVDimsPSCF,&
            & atomOverlapPSCF_did,atomOverlapCV_PSCF_did,atomOverlapPSCF_aid)
   else
      call gaussOverlapOL(numComponents,fullCVDimsPSCF,packedVVDimsPSCF,&
            & atomOverlapSYBD_PSCF_did,atomOverlapCV_SYBD_PSCF_did,&
            & atomOverlapSYBD_PSCF_aid)
   endif


   ! All ranks create the alpha distance matrix with nuclear alpha factor.
   call makeAlphaNucDist


   ! All ranks create the alpha distance matrix with potential alpha factor.
   call makeAlphaPotDist


   ! All ranks reallocate integrals from the overlap to the hamiltonian. FIX
   !   so that each rank only allocates the space that it needs.
   call reallocateIntegralsPSCF(coreDim,valeDim,numKPoints,spin)


   ! All ranks compute the hamiltonian matrix elements. FIX so that all ranks
   !   work tegether to compute the matrix. Note that only rank 0 records
   !   the results to disk.
   if (doSYBD_PSCF < 0) then
      call gaussOverlapHamPSCF(atomHamOverlapPSCF_did,atomHamOverlapPSCF_aid)
   else
      call gaussOverlapHamPSCF(atomHamOverlapSYBD_PSCF_did,&
            & atomHamOverlapSYBD_PSCF_aid)
   endif


   ! Now that all the single matrices are done being made each rank can
   !   deallocate the data structures that were used in all the above
   !   subroutines but are not necessary now.
   call cleanUpHamIntegrals


   ! If any supplementary three-term matrices (xyz) were requested for
   !   computing other properties, then allocate the space that is necessary
   !   for the integral calculation and do the integral.
   if ((doDIMO_PSCF == 1) .or. (doOPTC_PSCF >= 1)) then

      ! All ranks allocate memory for the 3 term integrals. FIX so that each
      !   rank only allocates the space that it needs. Consider in the future
      !   an option to do the XYZ independently to conserve memory if that
      !   becomes a problem.
      call allocateIntegrals3Terms(coreDim,valeDim,numKPoints)

      ! All ranks do the relevant integrals. FIX so all ranks work together
      !   to compute the integrals. Note that only rank 0 actually records
      !   the results to disk.
      if (doDIMO_PSCF == 1) then
         call gaussOverlapDM(packedVVDimsPSCF,atomDMOverlapPSCF_did,&
               & atomDMOverlapPSCF_aid)
      endif
      if (doOPTC_PSCF >= 1) then
         call gaussOverlapMM(packedVVDimsPSCF,atomMMOverlapPSCF_did,&
               & atomMMOverlapPSCF_aid)
      endif

      ! All ranks clean up from the 3Term integrals.
      call cleanUpIntegrals3Terms
   endif

   ! All ranks clean up.
   call secondCleanupIntegrals ! Overlap matrix parts used for ortho.
   call cleanUpGaussRelations
   if (doField_PSCF /= 1) then
      call cleanUpBasis
   endif

end subroutine intgPSCF


subroutine bandPSCF

      ! Import the necessary modules.
   use O_Kinds
   use O_MPI
   use O_Input,            only: numStates
   use O_Potential,        only: spin
   use O_AtomicSites,      only: valeDim
   use O_LAPACKParameters, only: setBlockSize
   use O_Input,            only: numStates
   use O_CommandLine,      only: doSYBD_PSCF
   use O_Lattice,          only: initializeLattice, initializeFindVec
   use O_KPoints,          only: makePathKPoints, numKPoints, &
         & computePhaseFactors
   use O_PSCFIntegralsHDF5, only: numComponents, atomOverlapPSCF_did,&
         & atomHamOverlapPSCF_did, atomOverlapSYBD_PSCF_did,&
         & atomHamOverlapSYBD_PSCF_did
   use O_PSCFEigVecHDF5, only: eigenVectorsPSCF_did,eigenVectorsSYBD_PSCF_did,&
         & eigenVectorsPSCF_aid, eigenVectorsSYBD_PSCF_aid
   use O_PSCFEigValHDF5, only: eigenValuesPSCF_did,eigenValuesSYBD_PSCF_did
   use O_SecularEquation, only: secularEqnPSCF, cleanUpSecularEqn

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define local variables.
   integer :: i


   ! All ranks set up LAPACK machine parameters. FIX so all ranks set up
   !   the ScaLAPACK machine parameters.
   call setBlockSize(valeDim)

   ! Compute the solid state wave function for each kpoint.

   ! Loop over all kpoints and compute the solid state wave function and
   !   energy eigen values for each. But first, ...
   ! Rank 0 records the fact that we are starting the k-point loop so that
   !   when someone looks at the output file as the job is running they will
   !   know that all the little dots represent a count of the number of
   !   k-points. Then they can figure out the progress and progress rate.
   if (mpiRank == 0) then
      write (20,*) "Beginning k-point loop."
      write (20,*) "Expecting ",numKPoints," iterations."
      call flush (20)
   endif

   ! All ranks compute the secular equation. FIX so that all ranks work
   !   together to compute the secular equation.
   do i = 1, spin
      if (doSYBD_PSCF < 0) then  ! No SYBD
         call secularEqnPSCF(i,numStates,numComponents,atomOverlapPSCF_did,&
               & atomHamOverlapPSCF_did,eigenValuesPSCF_did,&
               & eigenVectorsPSCF_did,eigenVectorsPSCF_aid)
      else
         call secularEqnPSCF(i,numStates,numComponents,&
               & atomOverlapSYBD_PSCF_did,&
               & atomHamOverlapSYBD_PSCF_did,eigenValuesSYBD_PSCF_did,&
               & eigenVectorsSYBD_PSCF_did,eigenVectorsSYBD_PSCF_aid)
      endif
   enddo

   ! Only rank 0 prints the band results if necessary. FIX so that all ranks
   !   contribute their portion of the solution to the SYBD result for
   !   printing. Note, a side effect of printSYBD is to populat and shift
   !   the energy eigenvalues.
   if (doSYBD_PSCF == 1) then
      call printSYBD
   endif

end subroutine bandPSCF


subroutine printSYBD

   ! Use necessary modules.
   use O_MPI
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

   if (mpiRank == 0) then
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
   endif

end subroutine printSYBD


subroutine dos(inSCF)

   ! Use necessary modules.
   use O_MPI
   use O_Potential,       only: spin
   use O_DOS,             only: computeDOS
   use O_Populate,        only: occupiedEnergy, populateStates
   use O_SecularEquation, only: shiftEnergyEigenValues


   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare passed parameters.
   integer, intent(in) :: inSCF

   ! Rank 0 opens the DOS files that will be written to.  If a spin polarized
   !   calculation is being done, then 60, 70, 80 hold spin up and 61, 71, 81
   !   hold spin down.  60,61=TDOS; 70,71= PDOS; 80,81=Localization Index.
   if (mpiRank == 0) then
      open (unit=60,file='fort.60',status='new',form='formatted')
      open (unit=70,file='fort.70',status='new',form='formatted')
      open (unit=80,file='fort.80',status='new',form='formatted')
      if (spin == 2) then
         open (unit=61,file='fort.61',status='new',form='formatted')
         open (unit=71,file='fort.71',status='new',form='formatted')
         open (unit=81,file='fort.81',status='new',form='formatted')
      endif
   endif


   ! All ranks populate the electron states to find the highest occupied
   !   state (Fermi energy for metals).
   call populateStates

   ! All ranks shift the energy eigen values according to the highest
   !   occupied state.
   call shiftEnergyEigenValues(occupiedEnergy)

   ! All ranks call the DOS subroutine to compute the total and partial
   !   density of states as well as the localization index. FIX so all ranks
   !   work together to compute the DOS, PDOS, and LI.
   call computeDOS(inSCF)

   ! Rank 0 closes the output files.
   if (mpiRank == 0) then
      close(60)
      close(70)
      close(80)
      if (spin == 2) then
         close(61)
         close(71)
         close(81)
      endif
   endif

end subroutine dos


subroutine bond(inSCF, doBond)

   ! Use necessary modules.
   use O_Kinds
   use O_MPI
   use O_Potential,       only: spin
   use O_Bond,            only: computeBond
   use O_Bond3C,          only: computeBond3C
   use O_Populate,        only: occupiedEnergy, populateStates
   use O_Input,           only: thermalSigma
   use O_SecularEquation, only: shiftEnergyEigenValues


   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare passed parameters.
   integer, intent(in) :: inSCF
   integer, intent(in) :: doBond

   ! Declare local variables.
   real(kind=double) :: thermalSigmaTemp


   ! Rank 0 opens the bond files that will be written to.
   if (mpiRank == 0) then
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
   endif


   ! The bond calculation must be done without any smearing to determine
   !   the number of electrons for each atom.  So all ranks set the thermal
   !   smearing factor to zero after saving its input value for later.
   thermalSigmaTemp = thermalSigma
   thermalSigma = 0.0_double

   ! All ranks populate the electron states to find the highest occupied
   !   state (Fermi energy for metals).
   call populateStates

   ! All ranks shift the energy eigen values according to the highest
   !   occupied state.
   call shiftEnergyEigenValues(occupiedEnergy)

   ! All ranks call the bond subroutine to compute the bond order and
   !   effective charge. FIX so all ranks work together.
   call computeBond(inSCF)

   ! All ranks compute the three center bond order if requested.
   if (doBond == 2) then
      call computeBond3C(inSCF)
   endif

   ! All ranks restore the thermal sigma to the input vale.
   thermalSigma = thermalSigmaTemp

   ! Rank 0 closes the BOND files that were written to.
   if (mpiRank == 0) then
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
   endif

end subroutine bond


subroutine dimo(inSCF)

   ! Use necessary modules.
   use O_MPI
   use O_Potential,       only: spin
   use O_Populate,        only: occupiedEnergy, populateStates
   use O_SecularEquation, only: shiftEnergyEigenValues
   use O_ValeCharge, only: makeValenceRho


   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare passed parameters.
   integer, intent(in) :: inSCF

   ! Rank 0 opens the DOS files that will be written to.  If a spin polarized
   !   calculation is being done, then 60, 70, 80 hold spin up and 61, 71, 81
   !   hold spin down.  60,61=TDOS; 70,71= PDOS; 80,81=Localization Index
   if (mpiRank == 0) then
      open (unit=74,file='fort.74',status='new',form='formatted')
      if (spin == 2) then
         open (unit=75,file='fort.75',status='new',form='formatted')
      endif
   endif

   ! All ranks populate the electron states to find the highest occupied
   !   state (Fermi energy for metals).
   call populateStates

   ! All ranks shift the energy eigen values according to the highest
   !   occupied state.
   call shiftEnergyEigenValues(occupiedEnergy)

   ! All ranks compute the valence charge density matrix. FIX so each rank
   !   only computes the portion of the valence charge density matrix that
   !   they are responsible for.
   call makeValenceRho(inSCF)

   ! Rank 0 closes the output files.
   if (mpiRank == 0) then
      close(74)
      if (spin == 2) then
         close(75)
      endif
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
   use O_MPI
   use O_ElementData,     only: initElementData
   use O_Populate,        only: populateStates
   use O_Potential,       only: initPotCoeffs
   use O_Field,           only: computeFieldMesh, cleanUpField
!   use O_SCFFieldHDF5,    only: wav_did, rho_did, pot_did, triggerAxis, &
!         & abcDimsChunk, fileFieldChunk_dsid
!   use O_PSCFFieldHDF5,    only: wavPSCF_did, rhoPSCF_did, potPSCF_did, &
!         & triggerAxisPSCF, abcDimsChunkPSCF, fileFieldChunkPSCF_dsid
   use O_Lattice,         only: initialize3DMesh
   use O_FieldHDF5,       only: prepFieldHDF5, closeFieldHDF5

   ! Make sure that there are not accidental variable declarations.
   implicit none


   ! Define passed parameters.
   integer, intent(in) :: inSCF


   ! Rank 0 initializes (or accesses) the HDF field file. There are no
   !   allocations done in here.
   if (mpiRank == 0) then
      call prepFieldHDF5(inSCF)
   endif


   ! Rank 0 opens the potential file that will be read from in this program.
   if (mpiRank == 0) then
      open (unit=8,file='fort.8',status='old',form='formatted')
   endif


   ! All ranks initialize element data from periodic table of the elements.
   call initElementData


   ! All ranks populate the electron states to find the highest occupied
   !   state (Fermi energy for metals).
   call populateStates


   ! All ranks initialize certain parameters for constructing and traversing
   !   the 3D mesh. FIX so that each rank computes the parameters as relevant
   !   for the portion of the 3D mesh that they are responsible for.
   call initialize3DMesh


   ! All ranks compute requested field values for each mesh point and rank 0
   !   stores the results in an HDF5+XDMF format. FIX so each rank only
   !   computes the portion of the field that they are responsible for.
   call computeFieldMesh(inSCF)


   ! All ranks clean up any left over arrays that need to be deallocated.
   call cleanUpField


   ! Rank 0 closes the field HDF5 file. There are no deallocations done here.
   if (mpiRank == 0) then
      call closeFieldHDF5
   endif

end subroutine field


subroutine optc(inSCF,doOPTC)

   ! Use necessary modules.
   use O_MPI
   use O_Potential,       only: spin
   use O_OptcPrint,       only: printOptcResults
   use O_Populate,        only: occupiedEnergy, populateStates
   use O_Lattice,         only: initializeLattice, initializeFindVec
   use O_Input,           only: detailCodePOPTC!, numStates, lastInitStatePACS
   use O_OptcTransitions, only: transCounter, energyDiff, transitionProb, &
         & transitionProbPOPTC, getEnergyStatistics, computeTransitions
   use O_SecularEquation, only: shiftEnergyEigenValues


   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: inSCF
   integer, intent(in) :: doOPTC


   ! Rank 0 opens the optc files that will be written to. 40,41=optical
   !   conductivity; 50,51 = imaginary dielectric function (usually). Note,
   !   for XANES/ELNES calculations, we do not need the optical conductivity.
   !   Also note, for Sigma(E) calculations the 50,51 units are the Sigma(E)
   !   itself, not the imaginary dielectric function (epsilon2).
   if (mpiRank == 0) then
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
   endif


   ! All ranks populate the electron states to find the highest occupied
   !   state (Fermi energy for metals).
   call populateStates

   ! FIX
!   if (doOPTC == 2) then ! Doing PACS calculation and we need to modify
!      ! the energy eigen values and their position. Now that the occupied
!      ! energy is known we can append the unoccupied excited states.
!      call appendExcitedEValsBand(lastInitStatePACS+1,numStates)
!   endif

   ! All ranks shift the energy eigen values according to the highest
   !   occupied state.
   call shiftEnergyEigenValues(occupiedEnergy)

   ! All ranks compute some statistics and variables concerning the energy
   !   values. FIX so each rank obtains parameters for only the portion of the
   !   work that they are responsible for.
   call getEnergyStatistics(doOPTC)

   ! All ranks compute the transition pairs and energy values of those
   !   transitions. FIX so each rank only computes the transitions that they
   !   are responsible for.
   call computeTransitions(inSCF,doOPTC)

   ! Rank 0 prints the output if necessary.
   if (mpiRank == 0) then
      if (doOPTC /= 3) then  ! Not doing a Sigma(E) calculation.
         call printOptcResults(doOPTC) ! Internally distinguishes optc, poptc.
      endif
   endif

   ! All ranks deallocate unused matrices
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
   use O_MPI
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

   ! Rank 0 closes the output files.
   if (mpiRank == 0) then
      close (7) ! Iteration data
      close (8) ! SCF Potential
      close (13) ! Magnetic moments
      close (14) ! Energy data per iteration
      if(doPSCF < 0) then
         close (20) ! Primary output
      endif
   endif

   ! All ranks deallocate all the other as of yet un-deallocated arrays.
   call cleanUpAtomTypes
   call cleanUpAtomSites
   call cleanUpPotTypes
   call cleanUpPotSites
   call cleanUpKPoints
   call cleanUpExchCorr
   call cleanUpPotential
   call cleanUpLattice

   ! Because we already deallocated in makeValenceRho all ranks only
   !   deallocate if necessary.
   if ((doDIMO_SCF < 0) .and. (doForce_SCF < 0)) then
      call cleanUpSecularEqn
   endif

   ! Rank 0 opens a file to signal completion of the program to the calling
   !   upolcao script.
   if (mpiRank == 0) then
      open (unit=2,file='fort.2',status='unknown')
   endif

end subroutine cleanUpSCF


subroutine cleanUpPSCF

   ! Use necessary modules.
   use O_MPI
   use O_CommandLine, only: doSYBD_PSCF
   use O_Lattice, only: cleanUpLattice
   use O_Potential, only:  cleanUpPotential
   use O_AtomicSites, only: cleanUpAtomSites
   use O_AtomicTypes, only: cleanUpAtomTypes, cleanupRadialFns
   use O_PotSites, only: cleanUpPotSites
   use O_PotTypes, only: cleanUpPotTypes
   use O_SecularEquation, only: cleanUpSecularEqn
   use O_KPoints, only: cleanUpKPoints

   implicit none

   ! Rank 0 closes any opened files.
   if (mpiRank == 0) then
      close (8)
      close (20)
   endif

   ! All ranks deallocate all the other as of yet un-deallocated arrays.
   call cleanUpRadialFns
   call cleanUpAtomTypes
   call cleanUpAtomSites
   call cleanUpPotTypes
   call cleanUpPotSites
   call cleanUpKPoints
   call cleanUpPotential
   call cleanUpLattice

   ! All ranks deallocate if necessary.
   if (doSYBD_PSCF < 0) then
      call cleanUpSecularEqn
   endif

   ! Rank 0 opens a file to signal completion of the program to the calling
   !   upolcao script.
   if (mpiRank == 0) then
      open (unit=2,file='fort.2',status='unknown')
   endif

end subroutine cleanUpPSCF


end module O_OLCAO
