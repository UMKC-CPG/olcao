module O_Wave

   ! Import necessary modules.
   use O_Kinds

   ! Define the module data used by multiple subroutines.
   integer :: numCols
   character*13, allocatable, dimension(:) :: colLabel ! Labels for the
         ! columns of data in the profile output files.
real (kind=double), allocatable, dimension (:,:) :: accumCharge  ! Accumulated
         ! charge for each spin,kPoint combination.
real (kind=double) :: accumChargeKP
   real (kind=double), allocatable, dimension (:,:,:) :: profile ! Index1=a,b,c
         ! Index2 = As identified in subroutine initEnv.  The important thing
         !          is that in the non spin polarized case the last Index2 is
         !          always the potential while in the spin polarized case the
         !          last two Index2 values are the spin-up and spin-down
         !          potentials.  This is important because they have different
         !          units that the other charge columns.
         ! Index3 = 1..max(numMeshPoints(a,b,c))
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: accumWaveFnCoeffs
#else
   real    (kind=double), allocatable, dimension (:,:) :: &
         & accumWaveFnCoeffsGamma
#endif

   ! Begin listing module subroutines.

contains

subroutine computeWaveFnMesh

   ! The goal of this subroutine

   ! Import the necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Potential,    only: spin, potCoeffs
   use O_Basis,        only: initializeAtomSite
   use O_Populate,     only: electronPopulation
   use O_Constants,    only: smallThresh, hartree
   use O_PotTypes,     only: maxNumPotAlphas, potTypes
   use O_AtomicSites,  only: valeDim, numAtomSites, atomSites
   use O_PSCFBandHDF5, only: eigenVectorsBand_did, valeStatesBand
   use O_Kpoints,      only: numKPoints, kPointWeight, phaseFactor
   use O_Input,        only: numStates, styleWAVE, eminWAVE, emaxWAVE, doRho, &
         & numElectrons
   use O_OpenDX,       only: printODXFieldHead, printODXFieldTail, &
         & printODXAtomPos, printODXLattice
   use O_AtomicTypes,  only: numAtomTypes, atomTypes, maxNumAtomAlphas, &
         & maxNumStates, maxNumValeRadialFns
   use O_Lattice,      only: logBasisFnThresh, numCellsReal, cellSizesReal, &
         & cellDimsReal, numMeshPoints, realVectors, realFractStrideLength, &
         & findLatticeVector
#ifndef GAMMA
   use O_MatrixSubs,      only: readMatrix
   use O_SecularEquation, only: valeVale, energyEigenValues
#else
   use O_MatrixSubs,      only: readMatrixGamma
   use O_SecularEquation, only: valeValeGamma, energyEigenValues
#endif

   ! Make sure that no variables are accidentally defined.
   implicit none

   ! Define local variables.
   integer :: a,b,c
   integer :: i,j,k,l,m,n  ! Loop index variables.
   integer :: dim1
   integer :: skipKP
   integer :: currentPointCount
   integer :: energyLevelCounter
   integer :: minStateIndex
   integer :: maxStateIndex
   real (kind=double) :: cutoff
   real (kind=double) :: currentPopulation
   real (kind=double) :: shiftedSepSqrd  ! Separation (r) between atom & mesh.
   real (kind=double), dimension(3) :: shiftedVec ! Vector between atom & mesh.
   real (kind=double), dimension(3) :: latticeVector
   real (kind=double), dimension(3) :: shiftedAtomPos
   real (kind=double), allocatable, dimension (:) :: negligLimit
   real (kind=double), allocatable, dimension (:,:,:) :: &
         & structuredElectronPopulation
#ifndef GAMMA
   real (kind=double), allocatable, dimension (:,:) :: tempRealValeVale
   real (kind=double), allocatable, dimension (:,:) :: tempImagValeVale
#endif


   ! Data for each mesh point.  For spin non-polarized calculations, we have:
   !   Index1 = Total, Index2 = Neutral, Index3 = Potential.  For spin
   !   polarized calculations we have:  Index1 = Spin up, Index2 = spin down,
   !   Index3 = Neutral, Index4 = Neutral.  These values will be appropriately
   !   combined to form the data listed in the initEnv subroutine.
   real (kind=double), dimension (spin) :: potFnEval
   real (kind=double), allocatable, dimension (:,:) :: currNumElec
!complex (kind=double), allocatable, dimension (:) :: tempPointValue
   real (kind=double), allocatable, dimension (:) :: currentPointValue
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:) :: waveFnEval
#else
   real (kind=double), allocatable, dimension (:) :: waveFnEvalGamma
#endif

   ! Atom specific variables that change with each atom loop iteration.
   integer                               :: currentCumulAlphaSum
   integer                               :: currentNumPotAlphas
   integer                               :: currentValeStateIndex
   integer                               :: currentMaxValeQN_l
   integer                               :: currentNumValeRadialFns
   integer,              dimension (2)   :: currentAtomType
   integer,              dimension (2)   :: currentElements
   integer,              dimension (2)   :: currentNumAlphas
   integer,              dimension (2)   :: currentNumCoreStates
   integer,              dimension (2)   :: currentNumValeStates
   integer,              dimension (2)   :: currentNumTotalStates
   integer,              dimension (4)   :: currentNumOrbAlphas
   integer, allocatable, dimension (:)   :: currentValeQN_lList
   integer, allocatable, dimension (:,:) :: currentlmIndex
   integer, allocatable, dimension (:,:) :: currentlmAlphaIndex
   real (kind=double)                    :: atomMeshSepSqrd
   real (kind=double)                    :: maxLatticeRadius
   real (kind=double), dimension (3,2)   :: currentPosition
   real (kind=double), allocatable, dimension (:)      :: currentPotAlphas
   real (kind=double), allocatable, dimension (:,:)    :: currentPotCoeffs
   real (kind=double), allocatable, dimension (:,:)    :: currentAlphas
   real (kind=double), allocatable, dimension (:,:,:)  :: currentBasisFns
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:,:) :: modifiedBasisFns
#else
   real (kind=double), allocatable, dimension (:,:,:)  :: modifiedBasisFnsGamma
#endif

   ! Atom specific variables that will change with each real space cell loop.
   integer                                           :: orbitalCount
   integer                                           :: currentOrbType
   integer                                           :: lastContribAlphaIndex
   integer                                           :: lastContribPotAlphaIndex
   integer                                           :: tempNumAlphas
   real (kind=double)                                :: x,y,z,rSqrd
   real (kind=double), dimension (16)                :: angularFactor
   real (kind=double), allocatable, dimension (:)    :: potAlphaDist
   real (kind=double), allocatable, dimension (:)    :: alphaDist
   real (kind=double), allocatable, dimension (:)    :: expPotAlphaDist
   real (kind=double), allocatable, dimension (:)    :: expAlphaDist
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: atomicOrbital
#else
   real (kind=double), allocatable, dimension (:,:)  :: atomicOrbitalGamma
#endif

   ! Log the beginning of the wave function evaluation.
   call timeStampStart(25)

   ! Initialize the module environment
   call initEnv

   ! Allocate space for locally defined allocatable arrays
#ifndef GAMMA
   allocate (waveFnEval            (numKPoints,numCols-spin))
   allocate (modifiedBasisFns      (maxNumAtomAlphas,maxNumStates,&
                                  & numCols-spin,numKPoints))
   allocate (atomicOrbital         (maxNumStates,numKPoints,numCols-spin))
#else
   allocate (waveFnEvalGamma       (numCols-spin))
   allocate (modifiedBasisFnsGamma (maxNumAtomAlphas,maxNumStates,numCols-spin))
   allocate (atomicOrbitalGamma    (maxNumStates,numCols-spin))
#endif
   allocate (profile               (3,numCols,maxVal(numMeshPoints(:))))
   allocate (structuredElectronPopulation (numStates,numKPoints,spin))
allocate (accumCharge           (spin,numKPoints))
!allocate (tempPointValue        (numCols-spin))
   allocate (currentPointValue     (numCols))
   allocate (currentBasisFns       (maxNumAtomAlphas,maxNumStates,2))
   allocate (currentPotCoeffs      (maxNumPotAlphas,spin))
   allocate (currentPotAlphas      (maxNumPotAlphas))
   allocate (currentAlphas         (maxNumAtomAlphas,2))
   allocate (currentlmAlphaIndex   (maxNumAtomAlphas,2))
   allocate (currentlmIndex        (maxNumStates,2))
   allocate (currentValeQN_lList   (maxNumValeRadialFns))
   allocate (potAlphaDist          (maxNumPotAlphas))
   allocate (alphaDist             (maxNumAtomAlphas))
   allocate (expPotAlphaDist       (maxNumPotAlphas))
   allocate (expAlphaDist          (maxNumAtomAlphas))

allocate (currNumElec (numKPoints,spin))

   ! If we will create an OpenDX file, then we will print the header for the
   !   field data, the lattice information, and the atomic positions now.
   if ((styleWAVE == 1) .or. (styleWAVE == 2)) then
!write (20,*) "numCols =",numCols
      call printODXFieldHead (numCols)
      call printODXAtomPos
      call printODXLattice
   endif

   ! Initialize the profile data structure.
   profile (:,:,:) = 0.0_double

   ! Fill a matrix of electron populations from the electron population that
   !   was computed in populateLevels.  Note that electronPopulation is a one
   !   dimensional array that has some order, but is not sorted in the way
   !   that the energy eigen values were sorted.  Please read the comments in
   !   the populateLevels subroutine to understand the order.
   !   (You can also probably get it from the loop order here ;)
   energyLevelCounter=0
   do i = 1, numKPoints
      do j = 1, spin
         do k = 1, numStates
            energyLevelCounter = energyLevelCounter + 1
            structuredElectronPopulation (k,i,j) = electronPopulation(&
                  & energyLevelCounter)
         enddo
      enddo
   enddo

   ! Define whether the packed arrays have two rows (complex) or one (real).
#ifndef GAMMA
   dim1 = 2
#else
   dim1 = 1
#endif


   ! Allocate space to read the wave functions and accumulate the coeffs.  Also
   !   initialize the coefficient accumulators.
#ifndef GAMMA
!   allocate (valeVale(valeDim,valeDim,1,spin))
   allocate (valeVale(valeDim,numStates,1,spin))
   allocate (tempRealValeVale(valeDim,numStates))
   allocate (tempImagValeVale(valeDim,numStates))
   allocate (accumWaveFnCoeffs(valeDim,numCols-spin,numKPoints)) ! No coeffs
         ! are needed for the total or spin-up, spin-down potential.
   accumWaveFnCoeffs(:,:,:) = cmplx(0.0_double,0.0_double)
#else
!   allocate (valeValeGamma(valeDim,valeDim,spin))
   allocate (valeValeGamma(valeDim,numStates,spin))
   allocate (accumWaveFnCoeffsGamma(valeDim,numCols-spin)) ! No coeffs are
         ! for the total or spin-up, spin-down potential.  (Only 1 kp.)
   accumWaveFnCoeffsGamma(:,:) = 0.0_double
#endif


   ! For each atomic type in the system compute the negligability limit.  If
   !   the square of the separation between an atomic site and a mesh point is
   !   greater than this computed value then the contribution of this atom is
   !   considered to be negligable.  This negligability limit is computed based
   !   on the basis function negligability limit for gaussian integration
   !   which is derived from the olcao.dat input value of the BASISFUNCTION
   !   CUTOFF (typically 0.1e-15).  The difference is that this calculation
   !   will use a cutoff that is significantly reduced (0.5 * exponent value)
   !   (typcially 0.1e-7.5).
   allocate (negligLimit(numAtomTypes))  ! Many values may be the same.
   cutoff = logBasisFnThresh/2.0_double
!write (20,*) "cutoff=",cutoff
!call flush (20)
   do i = 1, numAtomTypes
      negligLimit(i) = cutoff / atomTypes(i)%alphas(1)
!write (20,*) "negligLimit(i),i=",negligLimit(i),i
!call flush (20)
   enddo

accumCharge(:,:) = 0.0_double
accumChargeKP = 0.0_double
   do i = 1, numKPoints

      ! Skip any kpoints with a negligable contribution for each state.
      skipKP = 0
      do j = 1, numStates
         if (sum(abs(structuredElectronPopulation(j,i,:)))>smallThresh) then
            skipKP = 1
            exit
         endif
      enddo
      if (skipKP == 0) then
!write (20,*) "Skipping kpoint ",i
!call flush (20)
         cycle
      endif

      ! Read the computed wave function into the valeVale matrix.  I know
      !   that the name is misleading, but so far this seems to work out
      !   the best in terms of not having too many names or too many
      !   allocate deallocate calls.
!write (20,*) "valeStatesBand = ",valeStatesBand
!write (20,*) "valeDim, numStates=",valeDim,numStates
      do j = 1, spin
#ifndef GAMMA
         call readMatrix (eigenVectorsBand_did(:,i,j),&
               & valeVale(:,:numStates,1,j),&
               & tempRealValeVale(:,:numStates),&
               & tempImagValeVale(:,:numStates),&
               & valeStatesBand,valeDim,numStates)
#else
         call readMatrixGamma (eigenVectorsBand_did(1,i,j),&
               & valeValeGamma(:,:numStates,j),&
               & valeStatesBand,valeDim,numStates)
#endif
      enddo

      ! Accumulate the wave function coefficients including either the kpoint
      !   weight factor or the charge occupancy weight factor depending on the
      !   value of doRho.
      do j = 1, spin

         ! Find the minimum and maximum state indices to include in the wave
         !   function summation for this kpoint and spin.
         minStateIndex = 0
         do k = 1, numStates
            if (energyEigenValues(k,i,j)*hartree >= eminWAVE) then
               minStateIndex = k
               exit
            endif
         enddo

         ! Assume that the last state will be the max state index.  This will
         !   be fixed to a lower state if one is found that exceeds the
         !   requested highest energy state.
         maxStateIndex = numStates
         do k = 2, numStates
            if (energyEigenValues(k,i,j)*hartree > emaxWAVE) then
               maxStateIndex = k-1
               exit
            endif
         enddo

currNumElec(i,j) = sum(structuredElectronPopulation(&
      & minStateIndex:maxStateIndex,i,j))
!write (20,*) "i,j,currNumElec(i,j) = ",i,j,currNumElec(i,j)
!write (20,*) "minStateIndex=",minStateIndex
!write (20,*) "maxStateIndex=",maxStateIndex
#ifndef GAMMA
         accumCharge(j,i) = 0.0_double
         do k = minStateIndex, maxStateIndex
!tempPointValue(1)=cmplx(0.0_double,0.0_double)
!do l = 1, valeDim
!!do m = 1, valeDim
!tempPointValue(1) = tempPointValue(1) + &
!   & conjg(valeVale(l,k,1,j))*valeVale(l,k,1,j)
!!enddo
!enddo
!write (20,*) "k,valeVale state squared=",k,tempPointValue(1)

            if (doRho == 0) then
!write (20,*) "kPointWeight(i),spin=",kPointWeight(i),spin
!call flush (20)
               currentPopulation = kPointWeight(i)/real(spin,double)
            else
!write (20,*) "k,i,j structEPop(k,i,j)=",&
! & k,i,j,structuredElectronPopulation(k,i,j)
!call flush (20)
               currentPopulation = structuredElectronPopulation(k,i,j)
               accumCharge(j,i) = accumCharge(j,i) + currentPopulation
            endif
!write (20,*) "before aWFC=",accumWaveFnCoeffs(:,j,i)
!write (20,*) "vV=",valeVale(:,k,1,j)
            accumWaveFnCoeffs(:,j,i) = accumWaveFnCoeffs(:,j,i) + &
                  & currentPopulation * valeVale(:,k,1,j)
         enddo
!write (20,*) "accumCharge=",accumCharge(j,i)
accumChargeKP = accumChargeKP + accumCharge(j,i)
!write (20,*) "before divide aWFC=",accumWaveFnCoeffs(:,j,i)
!tempPointValue(1)=cmplx(0.0_double,0.0_double)
!do k = 1, valeDim
!!do l = 1, valeDim
!tempPointValue(1) = tempPointValue(1) + &
!   & conjg(accumWaveFnCoeffs(k,j,i))*accumWaveFnCoeffs(k,j,i)
!!enddo
!enddo
!write (20,*) "accum state squared=",tempPointValue(1)

         if (doRho == 1) then
            accumWaveFnCoeffs(:,j,i) = accumWaveFnCoeffs(:,j,i) / &
                  & sqrt(dot_product(accumWaveFnCoeffs(:,j,i),&
                  & accumWaveFnCoeffs(:,j,i)))
         endif
!currNumElec(i,j) = sqrt(dot_product(accumWaveFnCoeffs(:,j,i),&
!& accumWaveFnCoeffs(:,j,i)))

!write (20,*) "after divide aWFC=",accumWaveFnCoeffs(:,j,i)
!tempPointValue(1)=cmplx(0.0_double,0.0_double)
!do k = 1, valeDim
!!do l = 1, valeDim
!tempPointValue(1) = tempPointValue(1) + &
!   & conjg(accumWaveFnCoeffs(k,j,i))*accumWaveFnCoeffs(k,j,i)
!!enddo
!enddo
!write (20,*) "accum state squared after divide=",tempPointValue(1)
#else
         do k = minStateIndex, maxStateIndex

            if (doRho == 0) then
               currentPopulation = kPointWeight(i)/real(spin,double)
            else
               currentPopulation = structuredElectronPopulation(k,i,j)
            endif

            accumWaveFnCoeffsGamma(:,j) = accumWaveFnCoeffsGamma(:,j) + &
                  & currentPopulation * valeValeGamma(:,k,j)
         enddo
#endif
      enddo
   enddo

!write (20,*) "accumChargeKP=",accumChargeKP

   ! Save special coefficients to represent the neutral non-interacting
   !   electronic distribution for each atom.
   call makeNeutralCoeffs

   ! There is a lot of calculation that needs to be done here.  We need to
   !   consider each mesh point in turn.  We calculate the accumulated
   !   contributions of each atomic orbital from each atom of each
   !   replicated cell in periodic space.

   ! The orbital contributions depend on the relative x,y,z location of the
   !   mesh point to the atomic site as well as r (the magnitude of the
   !   separation distance).  Recall that the orbitals have an angular
   !   component applied to them from the spherical harmonics.

   ! Initialize the count of the number of mesh points computed.
   currentPointCount = 0

   do c = 1, numMeshPoints(3)
   do b = 1, numMeshPoints(2)
   do a = 1, numMeshPoints(1)

      ! Initialize the evaluation of this mesh point.
      potFnEval(1:spin) = 0.0_double
#ifndef GAMMA
      waveFnEval(1:numKPoints,1:numCols-spin) = cmplx(0.0_double,0.0_double)
#else
      waveFnEvalGamma(1:numCols-spin) = 0.0_double
#endif


      ! Increment the counter of the number of mesh points computed so far.
      currentPointCount = currentPointCount + 1

      ! Compute the location in x,y,z coordinates of the current mesh point by
      !   multiplying the number of fractional strides taken by the lattice
      !   parameters.
      do i = 1, 3
         currentPosition(i,1) = &
               & ((a-1)*realFractStrideLength(1)) * realVectors(i,1) + &
               & ((b-1)*realFractStrideLength(2)) * realVectors(i,2) + &
               & ((c-1)*realFractStrideLength(3)) * realVectors(i,3)
      enddo
!write (20,*) "a,b,c=",a,b,c
!write (20,*) "currentMeshPos=",currentPosition(:,1)
!call flush (20)

      ! Initiate a loop over the number of atoms that will contribute to the
      !   waveFnSqrd at the current mesh point defined by a,b,c loop indices.
      do i = 1, numAtomSites

         ! Obtain key information about this atom.
         call initializeAtomSite(i,2,currentAtomType,currentElements,&
            & currentNumTotalStates,currentNumCoreStates,currentNumValeStates,&
            & currentNumAlphas,currentlmIndex,currentlmAlphaIndex,&
            & currentPosition,currentAlphas,currentBasisFns)
!write (20,*) "Starting atom site i=",i
!write (20,*) "currentAtomPos=",currentPosition(:,2)
!call flush (20)

         ! Obtain further key information about this atom.
         currentValeStateIndex  = atomSites(i)%cumulValeStates
         currentMaxValeQN_l     = atomTypes(currentAtomType(2))%maxValeQN_l
         currentValeQN_lList(:) = atomTypes(currentAtomType(2))%valeQN_lList(:)
         currentNumOrbAlphas(:) = &
               & atomTypes(currentAtomType(2))%numOrbAlphas(:)
         currentNumValeRadialFns = &
               atomTypes(currentAtomType(2))%numValeRadialFns

         ! Obtain information about this potential site.  Carefully note that
         !   I assumed the current atom type = the current pot type.  Bad. Bad.
         currentCumulAlphaSum = potTypes(currentAtomType(2))%cumulAlphaSum
         currentNumPotAlphas = potTypes(currentAtomType(2))%numAlphas
         currentPotAlphas(1:currentNumPotAlphas) = &
               & potTypes(currentAtomType(2))%alphas(1:currentNumPotAlphas)
         do j = 1, spin
            currentPotCoeffs(1:currentNumPotAlphas,j) = &
                  & potCoeffs(currentCumulAlphaSum+1:&
                  & currentCumulAlphaSum+currentNumPotAlphas,j)
         enddo

         ! Find the lattice point closest to the difference between the atom
         !   site and the mesh point site.
         call findLatticeVector((currentPosition(:,1)-currentPosition(:,2)),&
               & latticeVector)

         ! Determine the square of the minimum separation between the mesh
         !   point and the atomic site.
         atomMeshSepSqrd = sum((currentPosition(:,1) - currentPosition(:,2) - &
               & latticeVector(:))**2)

!write (20,*) "atomMeshSepSqrd,negligLimit=",atomMeshSepSqrd, &
!& negligLimit(currentAtomType(2))
!call flush (20)

         ! Compare the current separation to the maximum separation for this
         !   atomic type for the contribution to be considered non-negligable.
         if (atomMeshSepSqrd > negligLimit(currentAtomType(2))) cycle

         ! At this point we know that at least one atom from the set of
         !   lattice replicated atoms will contribute (possibly the atom in the
         !   original cell, but not necessarily).  We will now multiply the
         !   basis function coefficients for this atom by the accumulated wave
         !   function coefficients.
#ifndef GAMMA
         do j = 1, numKPoints
            do k = 1, numCols-spin  ! -spin because pot is not treated this way
               do l = 1,currentNumValeStates(2)
                  modifiedBasisFns(:,l,k,j) = currentBasisFns(:,l,2) * &
                        & accumWaveFnCoeffs(currentValeStateIndex+l,k,j)
               enddo
            enddo
         enddo
#else
         do k = 1, numCols-spin  ! -spin because pot is not treated this way
            do l = 1,currentNumValeStates(2)
               modifiedBasisFnsGamma(:,l,k) = currentBasisFns(:,l,2) * &
                     & accumWaveFnCoeffsGamma(currentValeStateIndex+l,k)
            enddo
         enddo
#endif

         ! Determine the maximum radius beyond which no lattice point will be
         !   considered to contribute to the wave function evaluation for this
         !   atom-mesh point pair.  (This is the law of cosines 
         !   c^2 = a^2 + b^2 + 2ab*cos(g) with g = gamma = angle between a and
         !   b.  Here atomMeshSepSqrd = a^2,negligLimit = b^2,g=0.)
         maxLatticeRadius = atomMeshSepSqrd + negligLimit(currentAtomType(2))+&
               & 2.0_double * sqrt(atomMeshSepSqrd * &
               & negligLimit(currentAtomType(2)))

         ! The maxLatticRadius will always be smaller than the one used for
         !   integration so we don't bother checking to see if we have enough
         !   lattice points to complete the calculation.

         ! Initiate a loop over the number of replicated real space cells.
         do j = 1, numCellsReal

            ! Exit the loop when we have exceeded the necessary number of
            !   lattice points based on distance.
!write (20,*) "j,cellSizesReal(j),maxLatticeRadius",j,cellSizesReal(j),&
!& maxLatticeRadius
!call flush (20)
            if (cellSizesReal(j) > maxLatticeRadius) exit

            ! Obtain the position of the atom shifted by the current lattice.
            shiftedAtomPos(:) = currentPosition(:,2) + latticeVector(:) + &
                  & cellDimsReal(:,j)

            ! Obtain the relative separation x,y,z vector.  This is the vector
            !   between the atomic site in some periodic cell and the mesh
            !   point in the origin cell.
            shiftedVec(:) =  shiftedAtomPos(:) - currentPosition(:,1)

            ! Obtain the squared seperation magnitude between the mesh point
            !   and the shifted position of the atom.
            shiftedSepSqrd = dot_product(shiftedVec(:),shiftedVec(:))

!write (20,*) "j,cellSizesReal(j),shiftedSepSqrd,negligLimit",&
!& j,cellSizesReal(j),shiftedSepSqrd,negligLimit(currentAtomType(2))
!call flush (20)

            ! Determine if this shifted atom position puts it outside of the
            !   above determined negligability limit for this atom pair.
            if (shiftedSepSqrd > negligLimit(currentAtomType(2))) cycle

            ! Now we start a loop over the alphas of this atom.  For each alpha
            !   the "alpha dist." is tested against the negligability distance.
            !   Once an alpha is found that has negligable contribution we can
            !   stop and we only need to consider the alphas up to that point.
            ! Note that alphaDist(k) is initially computed as alpha*r^2 and
            !   cutoff is calculated above as -log(threshold).  So, if the
            !   alphaDist(k) is *greater* than the cutoff, that means that
            !   the contribution is *less* than the minimum necessary to be
            !   non-negligable.  (Tricky!)
            ! Assume that the last alpha will be the last alpha to contribute.
            lastContribAlphaIndex = currentNumAlphas(2)
            do k = 1, currentNumAlphas(2)
               alphaDist(k) = currentAlphas(k,2)*shiftedSepSqrd
!write (20,*) "alphaDist(k),k,cutoff,currentAlphas(k,2)=",alphaDist(k),k,&
!& cutoff,currentAlphas(k,2)
!call flush(20)
               if (alphaDist(k) > cutoff) then
                  lastContribAlphaIndex = k
                  exit
               endif
            enddo
!write (20,*) "lastContribAlphaIndex=",lastContribAlphaIndex
!call flush (20)

            lastContribPotAlphaIndex = currentNumPotAlphas
            do k = 1, currentNumPotAlphas
               potAlphaDist(k) = currentPotAlphas(k)*shiftedSepSqrd

               if (potAlphaDist(k) > cutoff) then
                  lastContribPotAlphaIndex = k
                  exit
               endif
            enddo

            ! Compute the exponential for each of the alphas needed.  This will
            !   complete all the distance sensitive parts of the calculation.
            !   The rest is simply applying the appropriate coefficients and
            !   angular factors to evaulate the contribution of each type of
            !   orbital (e.g. s, px, py, pz, ...).
            expAlphaDist(1:lastContribAlphaIndex) = &
                  & exp(-alphaDist(1:lastContribAlphaIndex))
!write (20,*) "expAlphaDist(1),currentAlphas(1,2)=",&
!&expAlphaDist(1),currentAlphas(1,2)
!call flush (20)
            expPotAlphaDist(1:lastContribPotAlphaIndex) = &
                  & exp(-potAlphaDist(1:lastContribPotAlphaIndex))

            ! Note that the general expression for a Gaussian type orbital
            !   (GTO) is A * r^l * exp(-ar^2).  In OLCAO an atomic orbital is
            !   defined as a sum of N GTOs * Ylm(theta,phi).  Each of the GTOs
            !   has a different A and a (alpha) labeled Ai, and ai.  When
            !   expressed in cartesian coordinates the Ylm spherical harmonics
            !   have a 1/r^l term that will cancel the r^l term in the GTO.
            !   The coefficients (Ai) stored in the currentBasisFns(:,:,2)
            !   array for the current atom are actually the product of the
            !   regular GTO coefficients and the Ylm angular normalization
            !   coefficients.


            ! Compute the angular factors for all the necessary orbital types.
            !   Note that the order of the angular factors computed here must
            !   match the order of the coefficients computation in the
            !   renormalizeBasis subroutine in the basis.f90 file.
            angularFactor(1) = 1.0_double
            if (currentMaxValeQN_l > 0) then
               x=shiftedVec(1)
               y=shiftedVec(2)
               z=shiftedVec(3)

               angularFactor(2) = x
               angularFactor(3) = y
               angularFactor(4) = z
               if (currentMaxValeQN_l > 1) then
                  rSqrd = x*x + y*y + z*z

                  angularFactor(5) = x*y
                  angularFactor(6) = x*z
                  angularFactor(7) = y*z
                  angularFactor(8) = x**2 - y**2
                  angularFactor(9) = 3.0_double*z**2 - rSqrd
                  if (currentMaxValeQN_l > 2) then
                     angularFactor(10) = x*y*z
                     angularFactor(11) = z*(x**2 - y**2)
                     angularFactor(12) = x*(x**2 - 3.0_double*y**2)
                     angularFactor(13) = y*(y**2 - 3.0_double*x**2)
                     angularFactor(14) = z*(5.0_double*z**2 - 3.0_double*rSqrd)
                     angularFactor(15) = x*(5.0_double*z**2 - rSqrd)
                     angularFactor(16) = y*(5.0_double*z**2 - rSqrd)
                  endif
               endif
            endif


            ! For each orbital of this atom, at this position, accumulate the
            !   Gaussian coefficients of the contributing alphas.
            orbitalCount=0
            do k = 1, currentNumValeRadialFns

               ! Get the current orbital type (l quantum number)
               currentOrbType = currentValeQN_lList(k)

               ! Get a temporary number of alphas needed for the present k
               !   indexed valence radial function (derived from the QN_l of
               !   the radial function).
               tempNumAlphas = currentNumOrbAlphas(currentOrbType + 1)

               ! Use all the alphas for this radial function or only up to the
               !   last contributing one.
               tempNumAlphas = min(lastContribAlphaIndex,tempNumAlphas)

               ! Loop over the possible QN_m values for the current QN_l type
               !   of orbital.
               do l = 1, 2 * currentOrbType + 1

                  ! Increment the orbital counter to the next state.
                  orbitalCount = orbitalCount + 1

                  do m = 1, numKPoints

                     do n = 1, numCols-spin ! -spin b/c pot is done differently
                        ! Accumulate the product of the Gaussian coefficients
                        !   for this atomic orbital with the contributing
                        !   exponential factors (computed above).  This
                        !   accumulation will be the atomic orbital radial
                        !   function evaluated at the current mesh point for
                        !   each kpoint.
                        ! Then, multiply the atomic orbital radial function by
                        !   the appropriate angular coefficient to make the
                        !   actual atomic orbital evaluated at the current
                        !   mesh point.
#ifndef GAMMA
                        atomicOrbital(orbitalCount,m,n) = sum( &
                              & modifiedBasisFns(1:tempNumAlphas, &
                              & orbitalCount,n,m) * &
                              & expAlphaDist(1:tempNumAlphas)) * &
                              & abs(angularFactor(currentlmIndex( &
                              & orbitalCount,2)))
#else
                        atomicOrbitalGamma(orbitalCount,n) = sum( &
                              & modifiedBasisFnsGamma(1:tempNumAlphas,&
                              & orbitalCount,n) * &
                              & expAlphaDist(1:tempNumAlphas)) * &
                              & abs(angularFactor(currentlmIndex( &
                              & orbitalCount,2)))
#endif
                     enddo ! n   numCols-spin
                  enddo ! m   numKPoints

!write (20,*) "j oC aO(oC)",j,orbitalCount,&
!& atomicOrbital(orbitalCount,1)
!write (20,*) "currlmI aF(currlmI)",&
!& currentlmIndex(orbitalCount,2),&
!& angularFactor(currentlmIndex(orbitalCount,2))
!call flush (20)
               enddo ! l   Num of QN_m for this QN_l
            enddo ! k   numValeRadialWaveFns

            ! Accumulate the contribution of this replicated atom's orbitals
            !  onto the current mesh point.
#ifndef GAMMA
            do k = 1, numCols-spin ! -spin b/c pot is not treated this way
               do l = 1, numKPoints
                  waveFnEval(l,k) = waveFnEval(l,k) + &
                        & sum(atomicOrbital(1:orbitalCount,l,k)) * &
                        & phaseFactor(l,j)
               enddo
            enddo
!if (j==2) then
!write (20,*) "wFE11,wFE21=",waveFnEval(1,1),waveFnEval(2,1)
!endif
#else
            do k = 1, numCols-spin ! -spin b/c pot is not treated this way
               waveFnEvalGamma(k) = waveFnEvalGamma(k) + &
                     & sum(atomicOrbitalGamma(1:orbitalCount,k))
            enddo
#endif

            ! Accumulate the contribution of this replicated potential site's
            !   Gaussian set for the current mesh point.
            do k = 1, spin
               potFnEval(k) = potFnEval(k) + &
                     & sum(currentPotCoeffs(1:lastContribPotAlphaIndex,k) * &
                     & expPotAlphaDist(1:lastContribPotAlphaIndex))
            enddo

         enddo ! j=numCellsReal
!write (20,*) "Did j-1 cells, j-1=",j-1
      enddo ! i=numAtomSites

!write (20,*) "waveFnEval",waveFnEval(:)

      ! Obtain the square of the evaluated function.
#ifndef GAMMA
!write (20,*) "waveFnEval(:,1)=",waveFnEval(:,1)
!write (20,*) "currNumElec(:,1)=",currNumElec(:,1)
      if (spin == 1) then
         waveFnEval(:,1) = conjg(waveFnEval(:,1)) * waveFnEval(:,1)
         currentPointValue(1) = sum(waveFnEval(:,1) * currNumElec(:,1))
         waveFnEval(:,2) = conjg(waveFnEval(:,2)) * waveFnEval(:,2)
         currentPointValue(2) = sum(waveFnEval(:,2) * &
               & 0.5_double * kPointWeight(:))
      else
         waveFnEval(:,1) = conjg(waveFnEval(:,1)) * waveFnEval(:,1)
         currentPointValue(1) = sum(waveFnEval(:,1) * currNumElec(:,1))
         waveFnEval(:,2) = conjg(waveFnEval(:,2)) * waveFnEval(:,2)
         currentPointValue(2) = sum(waveFnEval(:,2) * currNumElec(:,2))
         waveFnEval(:,3) = conjg(waveFnEval(:,3)) * waveFnEval(:,3)
         currentPointValue(3) = sum(waveFnEval(:,3) * &
               & numElectrons * 0.5_double * kPointWeight(:))

         ! Compute the total charge as the sum of the up and down.
         currentPointValue(1) = currentPointValue(1) + currentPointValue(2)

         ! Compute the spin difference charge as the difference between the
         !   up and down.  Note that this computation *must* be done *after*
         !   the computation of the total charge.
         currentPointValue(2) = currentPointValue(1) - &
               & 2.0_double * currentPointValue(2)
      endif

      ! Compute (the charge density for this point) - (the charge density for
      !   this point assuming non-interacting neutral atoms).
      currentPointValue(spin+1) = currentPointValue(1) - &
            & currentPointValue(spin+1)
         
!write (20,*) "waveFnEval(:,1)=",waveFnEval(:,1)
!write (20,*) "cPV(1)=",currentPointValue(1)

#else
      currentPointValue(1:numCols-spin) = &
            & waveFnEvalGamma(1:numCols-spin) * waveFnEvalGamma(1:numCols-spin)
#endif

      ! Record the potential function evaluation.  Always the last or last two
      !   indices.
      currentPointValue(numCols-spin+1:numCols) = potFnEval(1:spin)

      ! If openDX Files are being created, then print to them.
      if ((styleWAVE == 1) .or. (styleWAVE == 2)) then
         do i = 1, numCols
            write (57+i,ADVANCE="NO",fmt="(1x,e13.5)") currentPointValue(i)

            ! Move to the next line if we are at the end of this one.
            if (mod(currentPointCount,5) .eq. 0) then
               write (57+i,*)
            endif
         enddo
      endif

      ! Accumulate the data for the profiles.
      do i = 1,numCols
         profile(1,i,a) = profile(1,i,a) + currentPointValue(i)  ! 1=a axis
         profile(2,i,b) = profile(2,i,b) + currentPointValue(i)  ! 2=b axis
         profile(3,i,c) = profile(3,i,c) + currentPointValue(i)  ! 3=c axis
      enddo

      ! Record that this loop has finished.
      if (mod(currentPointCount,10) .eq. 0) then
            write (20,ADVANCE="NO",FMT="(a1)") "|"
      else
         write (20,ADVANCE="NO",FMT="(a1)") "."
      endif
      if (mod(currentPointCount,50) .eq. 0) then
         write (20,*) " ",currentPointCount
      endif
      call flush (20)

   enddo
   enddo
   enddo
!write (20,*) "accumCharge(1,:)=",accumCharge(1,:)

   ! Finalize printing of the OpenDX field data.
   if ((styleWAVE == 1) .or. (styleWAVE == 2)) then

      ! Write a newline to finish out any uneven (incomplete) lines if needed.
      if (mod(currentPointCount,5) /= 0) then
         do i = 1,numCols
            write (57+i,*)
         enddo
      endif

      ! Print the tail for the field data.
      call printODXFieldTail (numCols)
   endif


   ! Print the profile data.
   call printProfileData


   ! Deallocate unnecessary arrays and matrices.
   deallocate (structuredElectronPopulation)
   deallocate (currentBasisFns)
   deallocate (currentAlphas)
   deallocate (currentlmAlphaIndex)
   deallocate (currentlmIndex)
   deallocate (currentValeQN_lList)
   deallocate (accumCharge)
#ifndef GAMMA
   deallocate (waveFnEval)
   deallocate (atomicOrbital)
   deallocate (modifiedBasisFns)
   deallocate (valeVale)
   deallocate (accumWaveFnCoeffs)
   deallocate (tempRealValeVale)
   deallocate (tempImagValeVale)
#else
   deallocate (waveFnEvalGamma)
   deallocate (atomicOrbitalGamma)
   deallocate (modifiedBasisFnsGamma)
   deallocate (valeValeGamma)
   deallocate (accumWaveFnCoeffsGamma)
#endif

   ! Log the end of the wave function evaluation.
   call timeStampEnd(25)

end subroutine computeWaveFnMesh


! This subroutine will set some default (hard coded) values for some array
!   sizes and other values.  (Basically, this will just be an easy place to
!   change values if the program is modified or extended in the future.)
subroutine initEnv

   ! Use necessary modules
   use O_Potential, only: spin

   ! Make sure that no variables are accidentally defined.
   implicit none

   numCols = spin*2+1

   allocate (colLabel(numCols))
   if (spin == 1) then
      colLabel(1) = '         rhoV'   ! Valence charge density (rho).
      colLabel(2) = '       rhoV-N'   ! Valence rho - Neutral atom rho.
      colLabel(3) = '          pot'   ! Electronic potential.
   else
      colLabel(1) = '   rhoV_up+dn'   ! Valence rho, spin up + spin down.
      colLabel(2) = '   rhoV_up-dn'   ! Valence rho, spin up - spin down.
      colLabel(3) = ' rhoV_up+dn-N'   ! Valence rho, up+dn - Neutral atom rho.
      colLabel(4) = '       pot-up'   ! Electronic potential for spin up.
      colLabel(5) = '       pot-dn'   ! Electronic potential for spin down.
   endif

end subroutine initEnv


subroutine printProfileData

   ! Use necessary modules
   use O_Kinds
   use O_Potential, only: spin
   use O_Constants, only: hartree, bohrRad
   use O_Lattice,   only: realMag, realFractCrossArea, realFractStrideLength, &
         & realPlaneAngles, numMeshPoints

   ! Make sure that no variables are accidentally defined.
   implicit none

   ! Declare local variables.
   integer :: i,j
   integer :: index1,index2
   real (kind=double), dimension(numCols) :: integral

   ! Open the profile data files.
   open (unit=30,file='fort.30',form='formatted')
   open (unit=31,file='fort.31',form='formatted')
   open (unit=32,file='fort.32',form='formatted')

   ! Print the header.
   write (30,fmt="(6a13)") "aPos",colLabel(:)
   write (31,fmt="(6a13)") "bPos",colLabel(:)
   write (32,fmt="(6a13)") "cPos",colLabel(:)

   ! Adjust the charge profiles to average out the cross sectional area effect
   !   and the fact that the profiles are simple accumulations of the other
   !   axes.  For the case of the potential we simply obtain the average over
   !   the plane.
   do i = 1, 3
!write (20,*) "i,realFractCrossArea(i)=",i,realFractCrossArea(i)
!write (20,*) "i,realFractStrideLength(i)=",i,realFractStrideLength(i)
!call flush (20)
      index1 = mod(i,3)+1
      index2 = mod(i+1,3)+1
      profile(i,1:numCols-spin,1:numMeshPoints(i)) = &
            & profile(i,1:numCols-spin,1:numMeshPoints(i)) * &
            & realFractCrossArea(i)
      profile(i,numCols-spin+1:numCols,1:numMeshPoints(i)) = &
            & profile(i,numCols-spin+1:numCols,1:numMeshPoints(i)) * hartree / &
            & numMeshPoints(index1) / numMeshPoints(index1)
   enddo

   ! Print the profiles.
   do i = 1,3
!write (20,*) "i,realMag(i)=",i,realMag(i)
!call flush (20)
      do j = 1, numMeshPoints(i)
         write (29+i,fmt="(6e13.4)") (j-1) * realFractStrideLength(i) * &
               & realMag(i) * bohrRad, profile(i,1:numCols-spin,j) / bohrRad, &
               & profile(i,numCols-spin+1:numCols,j)
      enddo
   enddo

   ! Close the profile data files.
   close (30)
   close (31)
   close (32)

   ! Integrate the charge and print the result to the primary output file.  For
   !   the case of the potential we simply obtain the average along each axis.
   do i = 1, 3
!write (20,*) "fSL(:)=",realFractStrideLength(:)
!write (20,*) "lM(:)=",realMag(:)
!write (20,*) "s(pA(:))=",sin(realPlaneAngles(:))
!call flush (20)
      integral(:) = 0.0_double

      ! Accumulate across all mesh points.
      do j = 1, numMeshPoints(i)
         integral(1:numCols-spin) = integral(1:numCols-spin) + &
               & profile(i,1:numCols-spin,j)
         integral(numCols-spin+1:numCols) = integral(numCols-spin+1:numCols) + &
               & profile(i,numCols-spin+1:numCols,j)
      enddo

      ! Factor in the spacing for charge, and average the potential.
      integral(1:numCols-spin) = integral(1:numCols-spin) * &
            & realFractStrideLength(i) * realMag(i) * sin(realPlaneAngles(i))
      integral(numCols-spin+1:numCols) = integral(numCols-spin+1:numCols) / &
            numMeshPoints(i)

      ! Print the results.
      write (20,fmt="(6a13)") colLabel(:),"Integrated"
      write (20,fmt="(5e13.4)") integral(:)
   enddo

end subroutine printProfileData


subroutine makeNeutralCoeffs

   ! Import any necessary modules.
   use O_Kinds
   use O_Constants,   only: maxOrbitals
   use O_ElementData, only: coreCharge, valeCharge
   use O_KPoints,     only: numKPoints
   use O_AtomicSites, only: valeDim, numAtomSites, atomSites
   use O_AtomicTypes, only: atomTypes
   use O_PotTypes,    only: potTypes
   use O_Potential,   only: spin

   ! Make sure no implicit variables are created.
   implicit none

   ! Define local variables.
   integer :: i,j,k
   integer :: currAtomType
   integer :: currQN_l
   integer :: totalStateCount
   real (kind=double) :: electronsPerQN_m
   real (kind=double), dimension(maxOrbitals) :: localValeCharge
   real (kind=double), allocatable, dimension(:) :: neutralCoeffs

   ! Allocate temp space to hold the neutral coefficients before copying them
   !   to the (possibly complex values) waveFnCoeffs.
   allocate (neutralCoeffs(valeDim))

   ! Initialize the count of the number of states in the system.  This will
   !   reach to valeDim at the end.
   totalStateCount = 0

   ! For each orbital of each atom we need to provide an electron occupation
   !   number that represents the distribution of electrons in a neutral non-
   !   interacting atom.  This is a crude approximation and not at all as
   !   rigorous as it should be.  I suppose that in the future this could be
   !   improved.
   ! Recall that each orbital is spin degenerate in the approximation to the
   !   neutral atom that we are using.
   do i = 1, numAtomSites

      ! Get the current atom type.
      currAtomType = atomSites(i)%atomTypeAssn

      ! Get a copy of the number of valence electrons in each s,p,d,f type
      !   atomic orbitals (in the neutral condition).  (We identify the index
      !   number by looking at the Z value of this atomic site's nuclear
      !   potential.  Yeah, I know...indirect.)
      localValeCharge(:) = valeCharge(:,int(potTypes(currAtomType)%nucCharge))

      ! In the case where this atom has no core radial basis functions we can
      !   deduce two possible reasons.  (1)  This atom actually never will have
      !   any core basis functions and all it's electrons are always considered
      !   to be in the valence (e.g. H and He).  (2)  This atom is a target
      !   atom in a XANES/ELNES calculation and it typically does have some
      !   number of core radial basis functions.
      ! For case #2 we must add the coreCharge from the database to the
      !   valeCharge from the database to make sure that we count all the
      !   electrons and populate all the atomic orbital basis functions
      !   properly.
      ! If there are zero core radial functions we can simply add the core
      !   charge to the valence charge to solve the problem with case #2, and
      !   this will not affect case #1 because the core charge values there are
      !   always zero.
      if (atomTypes(currAtomType)%numCoreRadialFns == 0) then
         localValeCharge(:) = localValeCharge(:) + &
               & coreCharge(:,int(potTypes(currAtomType)%nucCharge))
      endif
!write (20,*) "localVC=",localValeCharge(:)
!write (20,*) "numValeRadFns=",atomTypes(currAtomType)%numValeRadialFns

      do j = 1, atomTypes(currAtomType)%numValeRadialFns

         ! Get the current orbital type (l quantum number)
         currQN_l = atomTypes(currAtomType)%valeQN_lList(j)
!write (20,*) "j,currQN_l",j,currQN_l

         ! Determine the number of electrons to be evenly spread in all QN_m
         !   orbitals of this QN_l from the remaining available electrons for
         !   this QN_l.

         ! First check if there are any electrons left for these orbitals.  If
         !   not, then we will have 0 electrons per QN_m orbital.  If there are
         !   some then we must determine if it is enough to fill the orbital or
         !   if it will only partially fill it.  In either case we will
         !   subtract the electrons that are about to be spread from the local
         !   valence charge list of all electrons.
         ! Note that the +1 is because the currQN_l starts at 0 but arrays are
         !   indexed starting from 1.  This will be applied many times in the
         !   following sections.  Watch out because we must also compute
         !   2*QN_l+1 which has a different meaning.  (Number of QN_m for this
         !   QN_l.)
         if (localValeCharge(currQN_l+1) >= 0) then

            ! We have some electrons to spread.

            ! Check if we can fill this QN_l orbital.  (2 electrons per QN_m)
            if (localValeCharge(currQN_l+1) - 2*(2*currQN_l+1) >= 0) then

               ! Every QN_m gets 2 electrons.
               electronsPerQN_m = 2.0_double
!write (20,*) "if elecPerQN_m=",electronsPerQN_m

               ! Remove 2 electrons from the available set for each QN_m.
               localValeCharge(currQN_l+1) = localValeCharge(currQN_l+1) - &
                     & 2.0_double * (2*currQN_l+1) ! -2 for s, -6 for p, ...
!write (20,*) "if localVC=",localValeCharge(:)
            else

               ! Every QN_m gets x / (# QN_m orbitals) where x is the number of
               !   electrons remaining for this QN_l.
               electronsPerQN_m = localValeCharge(currQN_l+1) / &
                     & real((2*currQN_l+1),double)
!write (20,*) "else elecPerQN_m=",electronsPerQN_m

               ! We have used up all the electrons for this QN_l.
               localValeCharge(currQN_l+1) = 0.0_double
!write (20,*) "else localVC=",localValeCharge(:)
            endif
         else
            ! There are no electrons to spread.  (Either they were used up, or
            !   there were none to begin with.)
            electronsPerQN_m = 0.0_double
         endif

         ! Spread the electrons evenly to all QN_m orbitals of this QN_l.
         do k = 1, 2 * currQN_l + 1

            ! Increment the state counter.
            totalStateCount = totalStateCount + 1

            neutralCoeffs(totalStateCount) = electronsPerQN_m
         enddo
      enddo
   enddo

!write (20,*) "neutCoeffs=",neutralCoeffs(:)

#ifndef GAMMA
   ! Save the neutral coefficients in index 2 or 3 depending on whether or not
   !   the calculation is spin non-polarized or spin polarized.  Note that the
   !   non-neutral coefficients must also be a unit vector for evaluation
   !   on the mesh and then later scaled by the charge that this kpoint
   !   represents.
   ! Note that at this point we cannot get the neutral charge over a specific
   !   energy range because we do not have the energy values for each of the
   !   neutral atom electron states available here.  That would require a
   !   separate calculation for a neutral atom.  All we have available here is
   !   the total charge over all occupied atomic states.
   do i = 1, numKPoints
      accumWaveFnCoeffs(:,spin+1,i) = cmplx(neutralCoeffs(:),0.0_double) / &
            & sqrt(dot_product(neutralCoeffs(:),neutralCoeffs(:)))
   enddo
#else
   accumWaveFnCoeffsGamma(:,spin+1) = neutralCoeffs(:) / &
            & sqrt(dot_product(neutralCoeffs(:),neutralCoeffs(:)))
#endif

   ! All done.
   deallocate (neutralCoeffs)

end subroutine makeNeutralCoeffs


subroutine cleanUpWave

   ! Make sure that no variables are accidentally defined.
   implicit none

   ! Deallocate everything that was not previously deallocated.
   deallocate (colLabel)
   deallocate (profile)

end subroutine cleanUpWave

end module O_Wave
