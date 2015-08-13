module O_IntegralsSCF

   ! Import necessary modules.
   use O_Kinds
   use O_Constants
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define matrices for the interaction integrals in terms of the core and
   !   valence sections and their interaction.  It is necessary to retain the
   !   overlap valeVale and coreVale separately for the purpose of
   !   orthogonalization of the nuclear pot, electronic pot, and KE matrices.
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:)   :: coreCore
   complex (kind=double), allocatable, dimension (:,:)     :: valeCore
   complex (kind=double), allocatable, dimension (:,:,:)   :: coreVale
   complex (kind=double), allocatable, dimension (:,:,:)   :: coreValeOL
   complex (kind=double), allocatable, dimension (:,:,:,:) :: valeVale
   complex (kind=double), allocatable, dimension (:,:,:,:) :: valeValeOL
#else
   real    (kind=double), allocatable, dimension (:,:)     :: coreCoreGamma
   real    (kind=double), allocatable, dimension (:,:)     :: valeCoreGamma
   real    (kind=double), allocatable, dimension (:,:)     :: coreValeGamma
   real    (kind=double), allocatable, dimension (:,:)     :: coreValeOLGamma
   real    (kind=double), allocatable, dimension (:,:,:)   :: valeValeGamma
   real    (kind=double), allocatable, dimension (:,:,:)   :: valeValeOLGamma
#endif

   ! Define the parameters that are necessary to control the computation of
   !   each term in the total electronic potential.
   integer :: elecPotInteraction
   integer :: currPotTypeNumber
   integer :: currPotNumber
   integer :: currPotElement
   integer :: currAlphaNumber
   integer :: currMultiplicity
   real (kind=double) :: currPotAlpha
   integer (hsize_t), allocatable, dimension (:,:,:) :: anyElecPotInteraction
   integer :: hsize_t_bitsize ! Used to access any bit in the above array.
   integer :: matrixIndex
   integer :: cellIndex
 

   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine allocateIntegralsSCF(coreDim,valeDim,numKPoints)

   implicit none

   ! Define passed dummy parameters.
   integer :: coreDim
   integer :: valeDim
   integer :: numKPoints

#ifndef GAMMA
   allocate (coreCore (coreDim,coreDim,numKPoints))
   allocate (valeVale (valeDim,valeDim,numKPoints,1))
#else
   allocate (coreCoreGamma (coreDim,coreDim))
   allocate (valeValeGamma (valeDim,valeDim,1))
#endif

end subroutine allocateIntegralsSCF

! Standard two center overlap integral.
subroutine gaussOverlapOL

   ! Import necessary modules.
   use O_Kinds
   use O_Constants, only: dim3
   use O_KPoints, only: numKPoints
   use O_GaussianRelations, only: alphaDist
   use O_AtomicSites, only: valeDim, coreDim, numAtomSites, atomSites
   use O_AtomicTypes, only: maxNumAtomAlphas, maxNumStates, atomTypes
   use O_Lattice, only: numCellsReal, cellSizesReal, cellDimsReal, &
         & findLatticeVector
   use O_GaussianIntegrals, only: overlapInteg
   use O_Basis, only: initializeAtomSite
   use O_IntgSaving
   use O_TimeStamps

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables for logging and loop control
   integer :: i,j,k,l,m ! Loop index variables

   ! Atom specific variables that change with each atom pair loop iteration.
   integer,              dimension (2)    :: currentAtomType
   integer,              dimension (2)    :: currentElements
   integer,              dimension (2)    :: currentNumAlphas
   integer,              dimension (2)    :: currentNumCoreStates
   integer,              dimension (2)    :: currentNumValeStates
   integer,              dimension (2)    :: currentNumTotalStates
   integer, allocatable, dimension (:,:)  :: currentlmIndex
   integer, allocatable, dimension (:,:)  :: currentlmAlphaIndex
   real (kind=double), dimension (dim3,2) :: currentPosition
   real (kind=double), allocatable, dimension (:,:)   :: currentAlphas
   real (kind=double), allocatable, dimension (:,:,:) :: currentBasisFns

   ! Alpha loop variables.  Once an overlap has been determined to exist for
   !   two atoms (including lattice shifting) a loop is initiated that goes
   !   over the alphas of those two atoms.  These variables are used there.
   integer :: currentMode ! A simple flag indicating weather (1) the alphas from
         ! atom 1 are being compared to the smallest untested alpha of atom 2,
         ! or (2) the other way around.  Case (0) is just the initial state.
   integer :: maxAlpha1Used ! Simply the largest alpha used from the current
         ! atom pair.
   integer, dimension (2) :: alphaIndex ! This tracks which alphas of the
         ! current atom pair are being tested for overlap.
   real (kind=double), dimension (16,16) :: oneAlphaPair ! The overlap from
         ! one alpha pair.
   real (kind=double), allocatable, dimension (:,:,:) :: pairXBasisFn2 ! The
         ! above overlap times the basis function from atom 2.  This is
         ! accumulated with each iteration of the alpha loop.
   logical :: contrib  ! At least one alpha pair contributes.  For the nuc and
         ! elec calculations this also requires that the pot term contribute.


   ! Variables and data structures that change or are accumulated with each
   !   iteration of the lattice loop.
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: currentPair
#else
   real (kind=double), allocatable, dimension (:,:) :: currentPairGamma
#endif
   real (kind=double), allocatable, dimension (:,:) :: pairXBasisFn12


   ! Local position and direction vectors and radii
   real (kind=double), dimension (dim3) :: latticeVector ! Vector to lattice
         ! point closest to the difference between the unit cell positions for
         ! atom 1 and atom 2.
   real (kind=double), dimension (dim3) :: shiftedAtomPos ! The position of
         ! atom 2 shifted to each relevant lattice point.
   real (kind=double) :: atomSiteSepSqrd ! The square of the minimum distance
         ! seperating atom 1 and atom 2 according to their unit cell positions
         ! shifted by the lattice point closest to their difference.
   real (kind=double) :: shiftedAtomSiteSep ! The seperation distance between
         ! atom 1 and the shifted position of atom 2.
   real (kind=double) :: currentNegligLimit ! The distance beyond which all
         ! alpha pairs are considered to have no overlap.
   real (kind=double) :: maxLatticeRadius ! Maximum radius beyond which no
         ! lattice points will be considered for integration.

   ! Define variables for gauss integrals
   integer :: l1l2Switch
   integer, dimension(16) :: powerOfTwo = (/0,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3/)

   ! Record the beginning of this phase of the setup calculation.
   call timeStampStart (8)

   ! Allocate space for locally defined allocatable arrays
   allocate (currentBasisFns     (maxNumAtomAlphas,maxNumStates,2))
   allocate (currentAlphas       (maxNumAtomAlphas,2))
   allocate (currentlmAlphaIndex (maxNumAtomAlphas,2))
   allocate (currentlmIndex      (maxNumStates,2))
   allocate (pairXBasisFn2       (16,maxNumAtomAlphas,maxNumStates))
   allocate (pairXBasisFn12      (maxNumStates,maxNumStates))
#ifndef GAMMA
   allocate (coreValeOL          (coreDim,valeDim,numKPoints))
   allocate (currentPair         (maxNumStates,maxNumStates,numKPoints))
#else
   allocate (coreValeOLGamma     (coreDim,valeDim))
   allocate (currentPairGamma    (maxNumStates,maxNumStates))
#endif


   ! Initialize key matrices
#ifndef GAMMA
   coreValeOL    (:,:,:)      = 0.0_double
   valeVale      (:,:,:,:)    = 0.0_double
   coreCore      (:,:,:)      = 0.0_double
#else
   coreValeOLGamma    (:,:)   = 0.0_double
   valeValeGamma      (:,:,:) = 0.0_double
   coreCoreGamma      (:,:)   = 0.0_double
#endif

   ! Begin atom-atom overlap loops.
   do i = 1, numAtomSites

      ! Obtain local copies of key data from larger global data structures for
      !   the first looped atom.
      call initializeAtomSite(i,1,currentAtomType,currentElements,&
            & currentNumTotalStates,currentNumCoreStates,currentNumValeStates,&
            & currentNumAlphas,currentlmIndex,currentlmAlphaIndex,&
            & currentPosition,currentAlphas,currentBasisFns)

      ! Begin a loop over the other atoms in the system
      do j = i, numAtomSites

         ! Obtain local copies of key data from larger global data structures
         !   for the second looped atom.
         call initializeAtomSite(j,2,currentAtomType,currentElements,&
            & currentNumTotalStates,currentNumCoreStates,currentNumValeStates,&
            & currentNumAlphas,currentlmIndex,currentlmAlphaIndex,&
            & currentPosition,currentAlphas,currentBasisFns)


         ! At this point all the data that defines the nature of the two atoms
         !   in this pair have been copied and defined.

         ! Find the lattice point closest to the difference between the two
         !   atom sites.
         call findLatticeVector((currentPosition(:,1)-currentPosition(:,2)),&
               & latticeVector)

         ! Determine the square of the minimum seperation distance between the
         !   two atoms.
         atomSiteSepSqrd = sum((currentPosition(:,1) - currentPosition(:,2) - &
               & latticeVector(:))**2)

         ! Obtain the maximum distance from either atom where the overlap
         !   is considered to be non-negligable.
         currentNegligLimit = alphaDist(1,1,currentElements(1),&
               & currentElements(2))

         ! Determine if there are no alpha terms for this atom pair that fall
         !   within the current negligability limit.  Cycle if there are none.
         if (atomSiteSepSqrd > currentNegligLimit) cycle

         ! Determine the maximum radius beyond which no lattice point will be
         !   considered to contribute to the overlap integral for this atom
         !   pair.  (This is the law of cosines c^2 = a^2 + b^2 + 2ab*cos(g)
         !   with g = gamma = angle between a and b.  Here atomSiteSepSqrd=a^2,
         !   currentNegligLimit = b^2,g=0.)
         maxLatticeRadius = atomSiteSepSqrd + currentNegligLimit + 2.0_double*&
                          & sqrt(atomSiteSepSqrd * currentNegligLimit)

         ! Check if we have to check more lattice points than were determined
         !   to be needed by comparing the maxlatticeRadius to the maximum
         !   radius of the lattice points defined earlier.
         if (maxLatticeRadius > cellSizesReal(numCellsReal)) then
            write (20,*) 'More lattice points needed for this atom overlap pair'
            write (20,*) 'maxLatticeRadius requested=',maxLatticeRadius
            write (20,*) 'max available=',cellSizesReal(numCellsReal)
            stop
         endif

         ! Initialize the matrix for this atom pair that will record the
         !   overlap integrals with bloch vector (kpoint) phase factors.
#ifndef GAMMA
         currentPair(:,:,:) = 0.0_double
#else
         currentPairGamma(:,:) = 0.0_double
#endif

         ! Begin a loop over all the lattice points to shift the position of
         !   atom number 2 to all the replicated cells.
         do k = 1, numCellsReal

            ! Exit the loop when we have exceeded the necessary number of
            !   lattice points based on distance.
            if (cellSizesReal(k) > maxLatticeRadius) exit

            ! Obtain the position of atom #2 shifted by the current lattice.
            shiftedAtomPos(:) = currentPosition(:,2) + latticeVector(:) + &
                  & cellDimsReal(:,k)

            ! Obtain the square of the magnitude of the seperation vector
            !   between atom 1 and the shifted position of atom 2.
            shiftedAtomSiteSep = sum ((currentPosition(:,1) - &
                  & shiftedAtomPos(:))**2)

            ! Determine if this shifted atom position puts it outside of the
            !   above determined negligability limit for this atom pair.
            if (shiftedAtomSiteSep > currentNegligLimit) cycle


            ! Initialize a matrix to hold the product of the overlap integrals
            ! times the atom2 basis functions.
            pairXBasisFn2(:,:currentNumAlphas(1),&
                  &:currentNumTotalStates(2)) = 0.0_double

            ! Initialize a variable to track the largest atomic alpha used
            !   from atom 1.
            maxAlpha1Used = 0

            ! Now we start loops over the alphas of each atomic type.  For each
            !   alpha pair the overlap is tested to determine if it is
            !   negligable or not.  The alphas for atom 1 are tested in turn
            !   with the first alpha of atom 2.  After all the atom 1 alphas
            !   have been tested or the test fails we test all the alphas
            !   of atom 2 with the first alpha of atom 1 until all the alphas
            !   of atom 2 have been used or the test fails.  Then we increment
            !   the loop counter by 1.  The second iteration of the loop
            !   obviously works just like the first except that the only
            !   difference is that instead of starting with the first alphas
            !   of atom 1 and atom 2, we always start with the second alphas
            !   of both atoms.

            ! The loop must be up to the smaller number of alphas.
            do l = 1, min(currentNumAlphas(1),currentNumAlphas(2))

               ! Set the mode to zero to say that we are incrementing through
               !   neither the atom 1 alpha array nor the atom 2 alpha array.
               currentMode = 0

               ! Initialize the matrix index values for the two alphas.
               alphaIndex(:) = l

               ! Check if this alpha pair has any non-negligable contribution.
               if (alphaDist(l,l,currentElements(1),currentElements(2)) < &
                     & shiftedAtomSiteSep) exit

               ! Start looping over atomic alphas looking for a negligable
               !   contribution for each alpha pair.
               do while (.true.)

                  ! Check if we are going through atom 1 alphas.
                  if (currentMode == 1) then

                     ! Go to the next atom 1 alpha
                     alphaIndex(1) = alphaIndex(1) + 1

                     ! Check if there are no alphas left to do for atom 1
                     if (alphaIndex(1) > currentNumAlphas(1)) then

                        ! Switch mode to go through atom 2 alphas on the
                        !   current diagonal alpha of atom 1.
                        currentMode = 2
                        alphaIndex(1) = l
                        cycle
                     endif

                  ! Check if we are going through atom 2 alphas.
                  elseif (currentMode == 2) then

                     ! Go to the next atom 2 alpha
                     alphaIndex(2) = alphaIndex(2) + 1

                     ! Check if there are no alphas left to do for atom 2
                     if (alphaIndex(2) > currentNumAlphas(2)) exit
                  endif

                  ! Check if this atom alpha pair has any non negligable
                  !   overlap.
                  if (alphaDist(alphaIndex(1),alphaIndex(2),currentElements(1),&
                        & currentElements(2)) < shiftedAtomSiteSep) then

                     ! Switch mode to go through atom 2 alphas on the current
                     !   diagonal alpha of atom 1.
                     if (currentMode == 1) then
                        currentMode = 2
                        alphaIndex(1) = l
                        cycle
                     else
                        exit
                     endif
                  endif

                  ! At this point a sufficient overlap between atomic Gaussians
                  !   has been shown to exist.
                  contrib = .true.

                  ! Calculate the opcode to do the correct set of integrals
                  !   for the current alpha pair
                  l1l2Switch = ishft(1,&
                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
                        &+ ishft(16,&
                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))

                  ! We can proceed with the next step of the calculation.
                  !   This is the actual integral.             
                  call overlapInteg (currentAlphas(alphaIndex(1),1),&
                        & currentAlphas(alphaIndex(2),2), &
                        & currentPosition(:,1), shiftedAtomPos(:),&
                        & l1l2Switch, oneAlphaPair)

                  ! We can proceed with the next step of the calculation.
                  !   This is the actual integral.
                  !call overlapInteg (oneAlphaPair,&
                  !       & currentlmAlphaIndex (alphaIndex(1),1),&
                  !       & currentlmAlphaIndex (alphaIndex(2),2),&
                  !       & currentAlphas(alphaIndex(1),1),&
                  !       & currentAlphas(alphaIndex(2),2),&
                  !       & currentPosition(:,1),shiftedAtomPos(:))

                  ! Collect the results of the overlap of the current alpha
                  !   times the basis functions of atom 2.
                  if (contrib .eqv. .true.) then
                     do m = 1, currentNumTotalStates(2)
                        pairXBasisFn2(:currentlmAlphaIndex(alphaIndex(1),1),&
                              & alphaIndex(1),m) = &
                              & pairXBasisFn2(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),alphaIndex(1),m) + &
                              & oneAlphaPair(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),currentlmIndex(m,2)) * &
                              & currentBasisFns(alphaIndex(2),m,2)
                     enddo

                     ! Update the maximum alpha used from atom 1.
                     maxAlpha1Used = max(alphaIndex(1),maxAlpha1Used)
                  endif

                  ! Switch mode from the initial state of no searching to the
                  !   state of searching along alphas from atom 1.
                  if (currentMode == 0) then
                     currentMode = 1
                  endif
               enddo
            enddo  ! min number of alphas between the two atoms l

            ! At this point all the alpha loops are complete and we can form a
            !   product with the atom 1 basis function to give the overlap
            !   integral in a complete basis representation.
            call multWithBasisFn1 (currentBasisFns,pairXBasisFn2, &
                  & pairXBasisFn12,currentlmIndex,currentNumTotalStates, &
                  & maxAlpha1Used)

            ! Collect this atom 1, atom 2 basis function overlap matrix for all
            !   bloch vectors (kpoints) with phase factors appropriate to the
            !   current atom 2 lattice vector.  NOTE that the k index is over
            !   the number of cells in the superlattice.
#ifndef GAMMA
            call applyPhaseFactors (currentPair,pairXBasisFn12(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),k,0,0)
#else
            call applyPhaseFactorsGamma (currentPairGamma,pairXBasisFn12(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),k,0)
#endif
         enddo !(k superlattice)

         ! At this point all the lattice sums for the current atom pair are
         !   complete.

         ! So now we can arrange the data from this atom into a set of three
         !   large matrices.  A valence-valence matrix, a core-valence matrix,
         !   and a core-core matrix.

#ifndef GAMMA
         ! First we must make a correction for the atom 2 lattice origin shift.
         call kPointLatticeOriginShift (currentNumTotalStates,currentPair,&
               & latticeVector,numKPoints,0)
         call saveCurrentPair(i,j,numKPoints,currentPair,valeVale(:,:,:,1),&
               & coreValeOL,coreCore)
#else
         call saveCurrentPairGamma(i,j,currentPairGamma,&
               & valeValeGamma(:,:,1),coreValeOLGamma,coreCoreGamma)
#endif
      enddo ! (Atom loop #2)
   enddo    ! (Atom loop #1)

   ! At this point all atom pairs have been computed.  Now, we compute the
   !   valeVale with the effects of the core orthogonalization.  We can also
   !   deallocate a bunch of arrays and matrices that are no longer needed.
   deallocate (currentBasisFns)
   deallocate (currentAlphas)
   deallocate (currentlmAlphaIndex)
   deallocate (currentlmIndex)
   deallocate (pairXBasisFn2)
   deallocate (pairXBasisFn12)
#ifndef GAMMA
   deallocate (currentPair)
#else
   deallocate (currentPairGamma)
#endif

   ! Perform orthogonalization and save the results to disk.
   call orthoOL

   ! Record the completion of this gaussian integration set.
   call timeStampEnd (8)

end subroutine gaussOverlapOL

! Two center Kinetic Energy overlap integrals.
subroutine gaussOverlapKE

   ! Import necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Constants, only: dim3
   use O_KPoints, only: numKPoints
   use O_GaussianRelations, only: alphaDist
   use O_AtomicSites, only: valeDim, coreDim, numAtomSites, atomSites
   use O_AtomicTypes, only: maxNumAtomAlphas, maxNumStates, atomTypes
   use O_Lattice, only: numCellsReal, cellSizesReal, cellDimsReal, &
         & findLatticeVector
   use O_GaussianIntegrals, only: KEInteg
   use O_Basis, only: initializeAtomSite
   use O_IntgSaving

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables for logging and loop control
   integer :: i,j,k,l,m ! Loop index variables

   ! Atom specific variables that change with each atom pair loop iteration.
   integer,              dimension (2)    :: currentAtomType
   integer,              dimension (2)    :: currentElements
   integer,              dimension (2)    :: currentNumAlphas
   integer,              dimension (2)    :: currentNumCoreStates
   integer,              dimension (2)    :: currentNumValeStates
   integer,              dimension (2)    :: currentNumTotalStates
   integer, allocatable, dimension (:,:)  :: currentlmIndex
   integer, allocatable, dimension (:,:)  :: currentlmAlphaIndex
   real (kind=double), dimension (dim3,2) :: currentPosition
   real (kind=double), allocatable, dimension (:,:)   :: currentAlphas
   real (kind=double), allocatable, dimension (:,:,:) :: currentBasisFns

   ! Alpha loop variables.  Once an overlap has been determined to exist for
   !   two atoms (including lattice shifting) a loop is initiated that goes
   !   over the alphas of those two atoms.  These variables are used there.
   integer :: currentMode ! A simple flag indicating weather (1) the alphas from
         ! atom 1 are being compared to the smallest untested alpha of atom 2,
         ! or (2) the other way around.  Case (0) is just the initial state.
   integer :: maxAlpha1Used ! Simply the largest alpha used from the current
         ! atom pair.
   integer, dimension (2) :: alphaIndex ! This tracks which alphas of the
         ! current atom pair are being tested for overlap.
   real (kind=double), dimension (16,16) :: oneAlphaPair ! The overlap from
         ! one alpha pair.
   real (kind=double), allocatable, dimension (:,:,:) :: pairXBasisFn2 ! The
         ! above overlap times the basis function from atom 2.  This is
         ! accumulated with each iteration of the alpha loop.
   logical :: contrib  ! At least one alpha pair contributes.  For the nuc and
         ! elec calculations this also requires that the pot term contribute.


   ! Variables and data structures that change or are accumulated with each
   !   iteration of the lattice loop.
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: currentPair
#else
   real (kind=double), allocatable, dimension (:,:) :: currentPairGamma
#endif
   real (kind=double), allocatable, dimension (:,:) :: pairXBasisFn12


   ! Local position and direction vectors and radii
   real (kind=double), dimension (dim3) :: latticeVector ! Vector to lattice
         ! point closest to the difference between the unit cell positions for
         ! atom 1 and atom 2.
   real (kind=double), dimension (dim3) :: shiftedAtomPos ! The position of
         ! atom 2 shifted to each relevant lattice point.
   real (kind=double) :: atomSiteSepSqrd ! The square of the minimum distance
         ! seperating atom 1 and atom 2 according to their unit cell positions
         ! shifted by the lattice point closest to their difference.
   real (kind=double) :: shiftedAtomSiteSep ! The seperation distance between
         ! atom 1 and the shifted position of atom 2.
   real (kind=double) :: currentNegligLimit ! The distance beyond which all
         ! alpha pairs are considered to have no overlap.
   real (kind=double) :: maxLatticeRadius ! Maximum radius beyond which no
         ! lattice points will be considered for integration.

   ! Define variables for gauss integrals
   integer :: l1l2Switch
   integer, dimension(16) :: powerOfTwo = (/0,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3/)

   ! Record the beginning of this phase of the setup calculation.
   call timeStampStart (9)

   ! Allocate space for locally defined allocatable arrays
   allocate (currentBasisFns     (maxNumAtomAlphas,maxNumStates,2))
   allocate (currentAlphas       (maxNumAtomAlphas,2))
   allocate (currentlmAlphaIndex (maxNumAtomAlphas,2))
   allocate (currentlmIndex      (maxNumStates,2))
   allocate (pairXBasisFn2       (16,maxNumAtomAlphas,maxNumStates))
   allocate (pairXBasisFn12      (maxNumStates,maxNumStates))
#ifndef GAMMA
   allocate (coreVale            (coreDim,valeDim,numKPoints))
   allocate (currentPair         (maxNumStates,maxNumStates,numKPoints))
#else
   allocate (coreValeGamma       (coreDim,valeDim))
   allocate (currentPairGamma    (maxNumStates,maxNumStates))
#endif


   ! Initialize key matrices
#ifndef GAMMA
   coreVale      (:,:,:)      = 0.0_double
   valeVale      (:,:,:,:)    = 0.0_double
   coreCore      (:,:,:)      = 0.0_double
#else
   coreValeGamma      (:,:)   = 0.0_double
   valeValeGamma      (:,:,:) = 0.0_double
   coreCoreGamma      (:,:)   = 0.0_double
#endif

   ! Begin atom-atom overlap loops.
   do i = 1, numAtomSites

      ! Obtain local copies of key data from larger global data structures for
      !   the first looped atom.
      call initializeAtomSite(i,1,currentAtomType,currentElements,&
            & currentNumTotalStates,currentNumCoreStates,currentNumValeStates,&
            & currentNumAlphas,currentlmIndex,currentlmAlphaIndex,&
            & currentPosition,currentAlphas,currentBasisFns)

      ! Begin a loop over the other atoms in the system
      do j = i, numAtomSites

         ! Obtain local copies of key data from larger global data structures
         !   for the second looped atom.
         call initializeAtomSite(j,2,currentAtomType,currentElements,&
            & currentNumTotalStates,currentNumCoreStates,currentNumValeStates,&
            & currentNumAlphas,currentlmIndex,currentlmAlphaIndex,&
            & currentPosition,currentAlphas,currentBasisFns)


         ! At this point all the data that defines the nature of the two atoms
         !   in this pair have been copied and defined.

         ! Find the lattice point closest to the difference between the two
         !   atom sites.
         call findLatticeVector((currentPosition(:,1)-currentPosition(:,2)),&
               & latticeVector)

         ! Determine the square of the minimum seperation distance between the
         !   two atoms.
         atomSiteSepSqrd = sum((currentPosition(:,1) - currentPosition(:,2) - &
               & latticeVector(:))**2)

         ! Obtain the maximum distance from either atom where the overlap
         !   is considered to be non-negligable.
         currentNegligLimit = alphaDist(1,1,currentElements(1),&
               & currentElements(2))

         ! Determine if there are no alpha terms for this atom pair that fall
         !   within the current negligability limit.  Cycle if there are none.
         if (atomSiteSepSqrd > currentNegligLimit) cycle

         ! Determine the maximum radius beyond which no lattice point will be
         !   considered to contribute to the overlap integral for this atom
         !   pair.  (This is the law of cosines c^2 = a^2 + b^2 + 2ab*cos(g)
         !   with g = gamma = angle between a and b.  Here atomSiteSepSqrd=a^2,
         !   currentNegligLimit = b^2,g=0.)
         maxLatticeRadius = atomSiteSepSqrd + currentNegligLimit + 2.0_double*&
                          & sqrt(atomSiteSepSqrd * currentNegligLimit)

         ! Check if we have to check more lattice points than were determined
         !   to be needed by comparing the maxlatticeRadius to the maximum
         !   radius of the lattice points defined earlier.
         if (maxLatticeRadius > cellSizesReal(numCellsReal)) then
            write (20,*) 'More lattice points needed for this atom overlap pair'
            write (20,*) 'maxLatticeRadius requested=',maxLatticeRadius
            write (20,*) 'max available=',cellSizesReal(numCellsReal)
            stop
         endif

         ! Initialize the matrix for this atom pair that will record the
         !   overlap integrals with bloch vector (kpoint) phase factors.
#ifndef GAMMA
         currentPair(:,:,:) = 0.0_double
#else
         currentPairGamma(:,:) = 0.0_double
#endif

         ! Begin a loop over all the lattice points to shift the position of
         !   atom number 2 to all the replicated cells.
         do k = 1, numCellsReal

            ! Exit the loop when we have exceeded the necessary number of
            !   lattice points based on distance.
            if (cellSizesReal(k) > maxLatticeRadius) exit

            ! Obtain the position of atom #2 shifted by the current lattice.
            shiftedAtomPos(:) = currentPosition(:,2) + latticeVector(:) + &
                  & cellDimsReal(:,k)

            ! Obtain the square of the magnitude of the seperation vector
            !   between atom 1 and the shifted position of atom 2.
            shiftedAtomSiteSep = sum ((currentPosition(:,1) - &
                  & shiftedAtomPos(:))**2)

            ! Determine if this shifted atom position puts it outside of the
            !   above determined negligability limit for this atom pair.
            if (shiftedAtomSiteSep > currentNegligLimit) cycle


            ! Initialize a matrix to hold the product of the overlap integrals
            ! times the atom2 basis functions.
            pairXBasisFn2(:,:currentNumAlphas(1),&
                  &:currentNumTotalStates(2)) = 0.0_double

            ! Initialize a variable to track the largest atomic alpha used
            !   from atom 1.
            maxAlpha1Used = 0

            ! Now we start loops over the alphas of each atomic type.  For each
            !   alpha pair the overlap is tested to determine if it is
            !   negligable or not.  The alphas for atom 1 are tested in turn
            !   with the first alpha of atom 2.  After all the atom 1 alphas
            !   have been tested or the test fails we test all the alphas
            !   of atom 2 with the first alpha of atom 1 until all the alphas
            !   of atom 2 have been used or the test fails.  Then we increment
            !   the loop counter by 1.  The second iteration of the loop
            !   obviously works just like the first except that the only
            !   difference is that instead of starting with the first alphas
            !   of atom 1 and atom 2, we always start with the second alphas
            !   of both atoms.

            ! The loop must be up to the smaller number of alphas.
            do l = 1, min(currentNumAlphas(1),currentNumAlphas(2))

               ! Set the mode to zero to say that we are incrementing through
               !   neither the atom 1 alpha array nor the atom 2 alpha array.
               currentMode = 0

               ! Initialize the matrix index values for the two alphas.
               alphaIndex(:) = l

               ! Check if this alpha pair has any non-negligable contribution.
               if (alphaDist(l,l,currentElements(1),currentElements(2)) < &
                     & shiftedAtomSiteSep) exit

               ! Start looping over atomic alphas looking for a negligable
               !   contribution for each alpha pair.
               do while (.true.)

                  ! Check if we are going through atom 1 alphas.
                  if (currentMode == 1) then

                     ! Go to the next atom 1 alpha
                     alphaIndex(1) = alphaIndex(1) + 1

                     ! Check if there are no alphas left to do for atom 1
                     if (alphaIndex(1) > currentNumAlphas(1)) then

                        ! Switch mode to go through atom 2 alphas on the
                        !   current diagonal alpha of atom 1.
                        currentMode = 2
                        alphaIndex(1) = l
                        cycle
                     endif

                  ! Check if we are going through atom 2 alphas.
                  elseif (currentMode == 2) then

                     ! Go to the next atom 2 alpha
                     alphaIndex(2) = alphaIndex(2) + 1

                     ! Check if there are no alphas left to do for atom 2
                     if (alphaIndex(2) > currentNumAlphas(2)) exit
                  endif

                  ! Check if this atom alpha pair has any non negligable
                  !   overlap.
                  if (alphaDist(alphaIndex(1),alphaIndex(2),currentElements(1),&
                        & currentElements(2)) < shiftedAtomSiteSep) then

                     ! Switch mode to go through atom 2 alphas on the current
                     !   diagonal alpha of atom 1.
                     if (currentMode == 1) then
                        currentMode = 2
                        alphaIndex(1) = l
                        cycle
                     else
                        exit
                     endif
                  endif

                  ! At this point a sufficient overlap between atomic Gaussians
                  !   has been shown to exist.
                  contrib = .true.

	      ! Calculate the opcode to do the correct set of integrals
              ! for the current alpha pair
               l1l2Switch = ishft(1,&
                 &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
                 &+ ishft(16,&
                 &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))

              ! We can proceed with the next step of the calculation. This
              ! is the actual integral.
               call KEInteg (currentAlphas(alphaIndex(1),1),&
                 & currentAlphas(alphaIndex(2),2), &
                 & currentPosition(:,1), shiftedAtomPos(:),&
                 & l1l2Switch, oneAlphaPair)

                  ! We can proceed with the next step of the calculation. This
                  !   is the actual integral.
                  !call KEInteg (oneAlphaPair,&
                  !      & currentlmAlphaIndex (alphaIndex(1),1),&
                  !      & currentlmAlphaIndex (alphaIndex(2),2),&
                  !      & currentAlphas(alphaIndex(1),1),&
                  !      & currentAlphas(alphaIndex(2),2),&
                  !      & currentPosition(:,1),shiftedAtomPos(:))

                  ! Collect the results of the overlap of the current alpha
                  !   times the basis functions of atom 2.
                  if (contrib .eqv. .true.) then
                     do m = 1, currentNumTotalStates(2)
                        pairXBasisFn2(:currentlmAlphaIndex(alphaIndex(1),1),&
                              & alphaIndex(1),m) = &
                              & pairXBasisFn2(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),alphaIndex(1),m) + &
                              & oneAlphaPair(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),currentlmIndex(m,2)) * &
                              & currentBasisFns(alphaIndex(2),m,2)
                     enddo

                     ! Update the maximum alpha used from atom 1.
                     maxAlpha1Used = max(alphaIndex(1),maxAlpha1Used)
                  endif

                  ! Switch mode from the initial state of no searching to the
                  !   state of searching along alphas from atom 1.
                  if (currentMode == 0) then
                     currentMode = 1
                  endif
               enddo
            enddo  ! min number of alphas between the two atoms l

            ! At this point all the alpha loops are complete and we can form a
            !   product with the atom 1 basis function to give the overlap
            !   integral in a complete basis representation.
            call multWithBasisFn1 (currentBasisFns,pairXBasisFn2, &
                  & pairXBasisFn12,currentlmIndex,currentNumTotalStates, &
                  & maxAlpha1Used)

            ! Collect this atom 1, atom 2 basis function overlap matrix for all
            !   bloch vectors (kpoints) with phase factors appropriate to the
            !   current atom 2 lattice vector.  NOTE that the k index is over
            !   the number of cells in the superlattice.
#ifndef GAMMA
            call applyPhaseFactors (currentPair,pairXBasisFn12(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),k,0,0)
#else
            call applyPhaseFactorsGamma (currentPairGamma,pairXBasisFn12(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),k,0)
#endif
         enddo !(k superlattice)

         ! At this point all the lattice sums for the current atom pair are
         !   complete.

         ! So now we can arrange the data from this atom into a set of three
         !   large matrices.  A valence-valence matrix, a core-valence matrix,
         !   and a core-core matrix.

#ifndef GAMMA
         ! First we must make a correction for the atom 2 lattice origin shift.
         call kPointLatticeOriginShift (currentNumTotalStates,currentPair,&
               & latticeVector,numKPoints,0)
         call saveCurrentPair(i,j,numKPoints,currentPair,valeVale(:,:,:,1),&
               & coreVale,coreCore)
#else
         call saveCurrentPairGamma(i,j,currentPairGamma,&
               & valeValeGamma(:,:,1),coreValeGamma,coreCoreGamma)
#endif
      enddo ! (Atom loop #2)
   enddo    ! (Atom loop #1)

   ! At this point all atom pairs have been computed.  Now, we compute the
   !   valeVale with the effects of the core orthogonalization.  We can also
   !   deallocate a bunch of arrays and matrices that are no longer needed.
   deallocate (currentBasisFns)
   deallocate (currentAlphas)
   deallocate (currentlmAlphaIndex)
   deallocate (currentlmIndex)
   deallocate (pairXBasisFn2)
   deallocate (pairXBasisFn12)
#ifndef GAMMA
   deallocate (currentPair)
#else
   deallocate (currentPairGamma)
#endif

   ! Perform orthogonalization and save the results to disk.  The 2 is an
   !   operation code signifying that a non-overlap orthogonalization should be
   !   done, and specifically that the result is for kinetic energy and that
   !   it should be written to the KE portion of the hdf5 file.
   call ortho(2)

   ! Record the completion of this gaussian integration set.
   call timeStampEnd (9)
   

end subroutine gaussOverlapKE


! Nuclear potential three center with 1/r overlap integrals.
subroutine gaussOverlapNP

   ! Import necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Constants, only: dim3
   use O_KPoints, only: numKPoints
   use O_GaussianRelations, only: alphaDist
   use O_AtomicSites, only: valeDim, coreDim, numAtomSites, atomSites
   use O_AtomicTypes, only: maxNumAtomAlphas, maxNumStates, atomTypes
   use O_Lattice, only: numCellsReal, cellSizesReal, cellDimsReal, &
         & findLatticeVector
   use O_Basis, only: initializeAtomSite
   use O_IntgSaving

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables for logging and loop control
   integer :: i,j,k,l,m ! Loop index variables

   ! Atom specific variables that change with each atom pair loop iteration.
   integer,              dimension (2)    :: currentAtomType
   integer,              dimension (2)    :: currentElements
   integer,              dimension (2)    :: currentNumAlphas
   integer,              dimension (2)    :: currentNumCoreStates
   integer,              dimension (2)    :: currentNumValeStates
   integer,              dimension (2)    :: currentNumTotalStates
   integer, allocatable, dimension (:,:)  :: currentlmIndex
   integer, allocatable, dimension (:,:)  :: currentlmAlphaIndex
   real (kind=double), dimension (dim3,2) :: currentPosition
   real (kind=double), allocatable, dimension (:,:)   :: currentAlphas
   real (kind=double), allocatable, dimension (:,:,:) :: currentBasisFns

   ! Alpha loop variables.  Once an overlap has been determined to exist for
   !   two atoms (including lattice shifting) a loop is initiated that goes
   !   over the alphas of those two atoms.  These variables are used there.
   integer :: currentMode ! A simple flag indicating weather (1) the alphas from
         ! atom 1 are being compared to the smallest untested alpha of atom 2,
         ! or (2) the other way around.  Case (0) is just the initial state.
   integer :: maxAlpha1Used ! Simply the largest alpha used from the current
         ! atom pair.
   integer, dimension (2) :: alphaIndex ! This tracks which alphas of the
         ! current atom pair are being tested for overlap.
   real (kind=double), dimension (16,16) :: oneAlphaPair ! The overlap from
         ! one alpha pair.
   real (kind=double), allocatable, dimension (:,:,:) :: pairXBasisFn2 ! The
         ! above overlap times the basis function from atom 2.  This is
         ! accumulated with each iteration of the alpha loop.
   logical :: contrib  ! At least one alpha pair contributes.  For the nuc and
         ! elec calculations this also requires that the pot term contribute.


   ! Variables and data structures that change or are accumulated with each
   !   iteration of the lattice loop.
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: currentPair
#else
   real (kind=double), allocatable, dimension (:,:) :: currentPairGamma
#endif
   real (kind=double), allocatable, dimension (:,:) :: pairXBasisFn12


   ! Local position and direction vectors and radii
   real (kind=double), dimension (dim3) :: latticeVector ! Vector to lattice
         ! point closest to the difference between the unit cell positions for
         ! atom 1 and atom 2.
   real (kind=double), dimension (dim3) :: shiftedAtomPos ! The position of
         ! atom 2 shifted to each relevant lattice point.
   real (kind=double) :: atomSiteSepSqrd ! The square of the minimum distance
         ! seperating atom 1 and atom 2 according to their unit cell positions
         ! shifted by the lattice point closest to their difference.
   real (kind=double) :: shiftedAtomSiteSep ! The seperation distance between
         ! atom 1 and the shifted position of atom 2.
   real (kind=double) :: currentNegligLimit ! The distance beyond which all
         ! alpha pairs are considered to have no overlap.
   real (kind=double) :: maxLatticeRadius ! Maximum radius beyond which no
         ! lattice points will be considered for integration.

   ! Record the beginning of this phase of the setup calculation.
   call timeStampStart (10)

   ! Allocate space for locally defined allocatable arrays
   allocate (currentBasisFns     (maxNumAtomAlphas,maxNumStates,2))
   allocate (currentAlphas       (maxNumAtomAlphas,2))
   allocate (currentlmAlphaIndex (maxNumAtomAlphas,2))
   allocate (currentlmIndex      (maxNumStates,2))
   allocate (pairXBasisFn2       (16,maxNumAtomAlphas,maxNumStates))
   allocate (pairXBasisFn12      (maxNumStates,maxNumStates))
#ifndef GAMMA
   allocate (coreVale            (coreDim,valeDim,numKPoints))
   allocate (currentPair         (maxNumStates,maxNumStates,numKPoints))
#else
   allocate (coreValeGamma       (coreDim,valeDim))
   allocate (currentPairGamma    (maxNumStates,maxNumStates))
#endif


   ! Initialize key matrices
#ifndef GAMMA
   coreVale      (:,:,:)      = 0.0_double
   valeVale      (:,:,:,:)    = 0.0_double
   coreCore      (:,:,:)      = 0.0_double
#else
   coreValeGamma      (:,:)   = 0.0_double
   valeValeGamma      (:,:,:) = 0.0_double
   coreCoreGamma      (:,:)   = 0.0_double
#endif

   ! Begin atom-atom overlap loops.
   do i = 1, numAtomSites

      ! Obtain local copies of key data from larger global data structures for
      !   the first looped atom.
      call initializeAtomSite(i,1,currentAtomType,currentElements,&
            & currentNumTotalStates,currentNumCoreStates,currentNumValeStates,&
            & currentNumAlphas,currentlmIndex,currentlmAlphaIndex,&
            & currentPosition,currentAlphas,currentBasisFns)

      ! Begin a loop over the other atoms in the system
      do j = i, numAtomSites

         ! Obtain local copies of key data from larger global data structures
         !   for the second looped atom.
         call initializeAtomSite(j,2,currentAtomType,currentElements,&
            & currentNumTotalStates,currentNumCoreStates,currentNumValeStates,&
            & currentNumAlphas,currentlmIndex,currentlmAlphaIndex,&
            & currentPosition,currentAlphas,currentBasisFns)


         ! At this point all the data that defines the nature of the two atoms
         !   in this pair have been copied and defined.

         ! Find the lattice point closest to the difference between the two
         !   atom sites.
         call findLatticeVector((currentPosition(:,1)-currentPosition(:,2)),&
               & latticeVector)

         ! Determine the square of the minimum seperation distance between the
         !   two atoms.
         atomSiteSepSqrd = sum((currentPosition(:,1) - currentPosition(:,2) - &
               & latticeVector(:))**2)

         ! Obtain the maximum distance from either atom where the overlap
         !   is considered to be non-negligable.
         currentNegligLimit = alphaDist(1,1,currentElements(1),&
               & currentElements(2))

         ! Determine if there are no alpha terms for this atom pair that fall
         !   within the current negligability limit.  Cycle if there are none.
         if (atomSiteSepSqrd > currentNegligLimit) cycle

         ! Determine the maximum radius beyond which no lattice point will be
         !   considered to contribute to the overlap integral for this atom
         !   pair.  (This is the law of cosines c^2 = a^2 + b^2 + 2ab*cos(g)
         !   with g = gamma = angle between a and b.  Here atomSiteSepSqrd=a^2,
         !   currentNegligLimit = b^2,g=0.)
         maxLatticeRadius = atomSiteSepSqrd + currentNegligLimit + 2.0_double*&
                          & sqrt(atomSiteSepSqrd * currentNegligLimit)

         ! Check if we have to check more lattice points than were determined
         !   to be needed by comparing the maxlatticeRadius to the maximum
         !   radius of the lattice points defined earlier.
         if (maxLatticeRadius > cellSizesReal(numCellsReal)) then
            write (20,*) 'More lattice points needed for this atom overlap pair'
            write (20,*) 'maxLatticeRadius requested=',maxLatticeRadius
            write (20,*) 'max available=',cellSizesReal(numCellsReal)
            stop
         endif

         ! Initialize the matrix for this atom pair that will record the
         !   overlap integrals with bloch vector (kpoint) phase factors.
#ifndef GAMMA
         currentPair(:,:,:) = 0.0_double
#else
         currentPairGamma(:,:) = 0.0_double
#endif

         ! Begin a loop over all the lattice points to shift the position of
         !   atom number 2 to all the replicated cells.
         do k = 1, numCellsReal

            ! Exit the loop when we have exceeded the necessary number of
            !   lattice points based on distance.
            if (cellSizesReal(k) > maxLatticeRadius) exit

            ! Obtain the position of atom #2 shifted by the current lattice.
            shiftedAtomPos(:) = currentPosition(:,2) + latticeVector(:) + &
                  & cellDimsReal(:,k)

            ! Obtain the square of the magnitude of the seperation vector
            !   between atom 1 and the shifted position of atom 2.
            shiftedAtomSiteSep = sum ((currentPosition(:,1) - &
                  & shiftedAtomPos(:))**2)

            ! Determine if this shifted atom position puts it outside of the
            !   above determined negligability limit for this atom pair.
            if (shiftedAtomSiteSep > currentNegligLimit) cycle


            ! Initialize a matrix to hold the product of the overlap integrals
            ! times the atom2 basis functions.
            pairXBasisFn2(:,:currentNumAlphas(1),&
                  &:currentNumTotalStates(2)) = 0.0_double

            ! Initialize a variable to track the largest atomic alpha used
            !   from atom 1.
            maxAlpha1Used = 0

            ! Now we start loops over the alphas of each atomic type.  For each
            !   alpha pair the overlap is tested to determine if it is
            !   negligable or not.  The alphas for atom 1 are tested in turn
            !   with the first alpha of atom 2.  After all the atom 1 alphas
            !   have been tested or the test fails we test all the alphas
            !   of atom 2 with the first alpha of atom 1 until all the alphas
            !   of atom 2 have been used or the test fails.  Then we increment
            !   the loop counter by 1.  The second iteration of the loop
            !   obviously works just like the first except that the only
            !   difference is that instead of starting with the first alphas
            !   of atom 1 and atom 2, we always start with the second alphas
            !   of both atoms.

            ! The loop must be up to the smaller number of alphas.
            do l = 1, min(currentNumAlphas(1),currentNumAlphas(2))

               ! Set the mode to zero to say that we are incrementing through
               !   neither the atom 1 alpha array nor the atom 2 alpha array.
               currentMode = 0

               ! Initialize the matrix index values for the two alphas.
               alphaIndex(:) = l

               ! Check if this alpha pair has any non-negligable contribution.
               if (alphaDist(l,l,currentElements(1),currentElements(2)) < &
                     & shiftedAtomSiteSep) exit

               ! Start looping over atomic alphas looking for a negligable
               !   contribution for each alpha pair.
               do while (.true.)

                  ! Check if we are going through atom 1 alphas.
                  if (currentMode == 1) then

                     ! Go to the next atom 1 alpha
                     alphaIndex(1) = alphaIndex(1) + 1

                     ! Check if there are no alphas left to do for atom 1
                     if (alphaIndex(1) > currentNumAlphas(1)) then

                        ! Switch mode to go through atom 2 alphas on the
                        !   current diagonal alpha of atom 1.
                        currentMode = 2
                        alphaIndex(1) = l
                        cycle
                     endif

                  ! Check if we are going through atom 2 alphas.
                  elseif (currentMode == 2) then

                     ! Go to the next atom 2 alpha
                     alphaIndex(2) = alphaIndex(2) + 1

                     ! Check if there are no alphas left to do for atom 2
                     if (alphaIndex(2) > currentNumAlphas(2)) exit
                  endif

                  ! Check if this atom alpha pair has any non negligable
                  !   overlap.
                  if (alphaDist(alphaIndex(1),alphaIndex(2),currentElements(1),&
                        & currentElements(2)) < shiftedAtomSiteSep) then

                     ! Switch mode to go through atom 2 alphas on the current
                     !   diagonal alpha of atom 1.
                     if (currentMode == 1) then
                        currentMode = 2
                        alphaIndex(1) = l
                        cycle
                     else
                        exit
                     endif
                  endif

                  ! At this point a sufficient overlap between atomic Gaussians
                  !   has been shown to exist.
                  contrib = .true.

                  ! We can proceed with the next step of the calculation. In
                  !   this case there will be a further loop over the nuclear
                  !   potential sites.
                  call nuclearPE (contrib,alphaIndex,currentElements,&
                        & currentlmAlphaIndex,shiftedAtomSiteSep,&
                        & currentPosition,shiftedAtomPos,oneAlphaPair,&
                        & currentAlphas)

                  ! Collect the results of the overlap of the current alpha
                  !   times the basis functions of atom 2.
                  if (contrib .eqv. .true.) then
                     do m = 1, currentNumTotalStates(2)
                        pairXBasisFn2(:currentlmAlphaIndex(alphaIndex(1),1),&
                              & alphaIndex(1),m) = &
                              & pairXBasisFn2(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),alphaIndex(1),m) + &
                              & oneAlphaPair(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),currentlmIndex(m,2)) * &
                              & currentBasisFns(alphaIndex(2),m,2)
                     enddo

                     ! Update the maximum alpha used from atom 1.
                     maxAlpha1Used = max(alphaIndex(1),maxAlpha1Used)
                  endif

                  ! Switch mode from the initial state of no searching to the
                  !   state of searching along alphas from atom 1.
                  if (currentMode == 0) then
                     currentMode = 1
                  endif
               enddo
            enddo  ! min number of alphas between the two atoms l

            ! At this point all the alpha loops are complete and we can form a
            !   product with the atom 1 basis function to give the overlap
            !   integral in a complete basis representation.
            call multWithBasisFn1 (currentBasisFns,pairXBasisFn2, &
                  & pairXBasisFn12,currentlmIndex,currentNumTotalStates, &
                  & maxAlpha1Used)

            ! Collect this atom 1, atom 2 basis function overlap matrix for all
            !   bloch vectors (kpoints) with phase factors appropriate to the
            !   current atom 2 lattice vector.  NOTE that the k index is over
            !   the number of cells in the superlattice.
#ifndef GAMMA
            call applyPhaseFactors (currentPair,pairXBasisFn12(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),k,0,0)
#else
            call applyPhaseFactorsGamma (currentPairGamma,pairXBasisFn12(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),k,0)
#endif
         enddo !(k superlattice)

         ! At this point all the lattice sums for the current atom pair are
         !   complete.

         ! So now we can arrange the data from this atom into a set of three
         !   large matrices.  A valence-valence matrix, a core-valence matrix,
         !   and a core-core matrix.

#ifndef GAMMA
         ! First we must make a correction for the atom 2 lattice origin shift.
         call kPointLatticeOriginShift (currentNumTotalStates,currentPair,&
               & latticeVector,numKPoints,0)
         call saveCurrentPair(i,j,numKPoints,currentPair,valeVale(:,:,:,1),&
               & coreVale,coreCore)
#else
         call saveCurrentPairGamma(i,j,currentPairGamma,&
               & valeValeGamma(:,:,1),coreValeGamma,coreCoreGamma)
#endif
      enddo ! (Atom loop #2)
   enddo    ! (Atom loop #1)

   ! At this point all atom pairs have been computed.  Now, we compute the
   !   valeVale with the effects of the core orthogonalization.  We can also
   !   deallocate a bunch of arrays and matrices that are no longer needed.
   deallocate (currentBasisFns)
   deallocate (currentAlphas)
   deallocate (currentlmAlphaIndex)
   deallocate (currentlmIndex)
   deallocate (pairXBasisFn2)
   deallocate (pairXBasisFn12)
#ifndef GAMMA
   deallocate (currentPair)
#else
   deallocate (currentPairGamma)
#endif

   ! Perform orthogonalization and save the results to disk.  The 3 is an
   !   operation code signifying that a non-overlap orthogonalization should be
   !   done, and specifically that the result is for the nuclear potential and
   !   that it should be written to the NP portion of the hdf5 file.
   call ortho (3)

   ! Record the completion of this gaussian integration set.
   call timeStampEnd (10)
   

end subroutine gaussOverlapNP

! Electronic Potential three center overlap integrals.
subroutine gaussOverlapEP

   ! Import necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Constants, only: dim3
   use O_KPoints, only: numKPoints
   use O_GaussianRelations, only: alphaDist
   use O_Lattice, only: numCellsReal, cellSizesReal, cellDimsReal, &
         & findLatticeVector
   use O_AtomicSites, only: valeDim, coreDim, numAtomSites, atomSites
   use O_AtomicTypes, only: maxNumAtomAlphas, maxNumStates, atomTypes
   use O_Basis, only: initializeAtomSite
   use O_IntgSaving
   use O_GaussianIntegrals

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables for logging and loop control
   integer :: i,j,k,l,m ! Loop index variables

   ! Atom specific variables that change with each atom pair loop iteration.
   integer,              dimension (2)    :: currentAtomType
   integer,              dimension (2)    :: currentElements
   integer,              dimension (2)    :: currentNumAlphas
   integer,              dimension (2)    :: currentNumCoreStates
   integer,              dimension (2)    :: currentNumValeStates
   integer,              dimension (2)    :: currentNumTotalStates
   integer, allocatable, dimension (:,:)  :: currentlmIndex
   integer, allocatable, dimension (:,:)  :: currentlmAlphaIndex
   real (kind=double), dimension (dim3,2) :: currentPosition
   real (kind=double), allocatable, dimension (:,:)   :: currentAlphas
   real (kind=double), allocatable, dimension (:,:,:) :: currentBasisFns

   ! Alpha loop variables.  Once an overlap has been determined to exist for
   !   two atoms (including lattice shifting) a loop is initiated that goes
   !   over the alphas of those two atoms.  These variables are used there.
   integer :: currentMode ! A simple flag indicating weather (1) the alphas from
         ! atom 1 are being compared to the smallest untested alpha of atom 2,
         ! or (2) the other way around.  Case (0) is just the initial state.
   integer :: maxAlpha1Used ! Simply the largest alpha used from the current
         ! atom pair.
   integer, dimension (2) :: alphaIndex ! This tracks which alphas of the
         ! current atom pair are being tested for overlap.
   real (kind=double), dimension (16,16) :: oneAlphaPair ! The overlap from
         ! one alpha pair.
   real (kind=double), allocatable, dimension (:,:,:) :: pairXBasisFn2 ! The
         ! above overlap times the basis function from atom 2.  This is
         ! accumulated with each iteration of the alpha loop.
   logical :: contrib  ! At least one alpha pair contributes.  For the nuc and
         ! elec calculations this also requires that the pot term contribute.


   ! Variables and data structures that change or are accumulated with each
   !   iteration of the lattice loop.
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: currentPair
#else
   real (kind=double), allocatable, dimension (:,:) :: currentPairGamma
#endif
   real (kind=double), allocatable, dimension (:,:) :: pairXBasisFn12


   ! Local position and direction vectors and radii
   real (kind=double), dimension (dim3) :: latticeVector ! Vector to lattice
         ! point closest to the difference between the unit cell positions for
         ! atom 1 and atom 2.
   real (kind=double), dimension (dim3) :: shiftedAtomPos ! The position of
         ! atom 2 shifted to each relevant lattice point.
   real (kind=double) :: atomSiteSepSqrd ! The square of the minimum distance
         ! seperating atom 1 and atom 2 according to their unit cell positions
         ! shifted by the lattice point closest to their difference.
   real (kind=double) :: shiftedAtomSiteSep ! The seperation distance between
         ! atom 1 and the shifted position of atom 2.
   real (kind=double) :: currentNegligLimit ! The distance beyond which all
         ! alpha pairs are considered to have no overlap.
   real (kind=double) :: maxLatticeRadius ! Maximum radius beyond which no
         ! lattice points will be considered for integration.

   ! Define variables for gauss integrals
   integer :: l1l2Switch
   integer, dimension(16) :: powerOfTwo = (/0,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3/)

   ! Allocate space for locally defined allocatable arrays
   allocate (currentBasisFns     (maxNumAtomAlphas,maxNumStates,2))
   allocate (currentAlphas       (maxNumAtomAlphas,2))
   allocate (currentlmAlphaIndex (maxNumAtomAlphas,2))
   allocate (currentlmIndex      (maxNumStates,2))
   allocate (pairXBasisFn2       (16,maxNumAtomAlphas,maxNumStates))
   allocate (pairXBasisFn12      (maxNumStates,maxNumStates))
#ifndef GAMMA
   allocate (coreVale            (coreDim,valeDim,numKPoints))
   allocate (currentPair         (maxNumStates,maxNumStates,numKPoints))
#else
   allocate (coreValeGamma       (coreDim,valeDim))
   allocate (currentPairGamma    (maxNumStates,maxNumStates))
#endif


   ! Initialize key matrices
#ifndef GAMMA
   coreVale      (:,:,:)      = 0.0_double
   valeVale      (:,:,:,:)    = 0.0_double
   coreCore      (:,:,:)      = 0.0_double
#else
   coreValeGamma      (:,:)   = 0.0_double
   valeValeGamma      (:,:,:) = 0.0_double
   coreCoreGamma      (:,:)   = 0.0_double
#endif

   ! Begin atom-atom overlap loops.
   do i = 1, numAtomSites

      ! Obtain local copies of key data from larger global data structures for
      !   the first looped atom.
      call initializeAtomSite(i,1,currentAtomType,currentElements,&
            & currentNumTotalStates,currentNumCoreStates,currentNumValeStates,&
            & currentNumAlphas,currentlmIndex,currentlmAlphaIndex,&
            & currentPosition,currentAlphas,currentBasisFns)

      ! Begin a loop over the other atoms in the system
      do j = i, numAtomSites

         ! Obtain local copies of key data from larger global data structures
         !   for the second looped atom.
         call initializeAtomSite(j,2,currentAtomType,currentElements,&
            & currentNumTotalStates,currentNumCoreStates,currentNumValeStates,&
            & currentNumAlphas,currentlmIndex,currentlmAlphaIndex,&
            & currentPosition,currentAlphas,currentBasisFns)


         ! At this point all the data that defines the nature of the two atoms
         !   in this pair have been copied and defined.

         ! Find the lattice point closest to the difference between the two
         !   atom sites.
         call findLatticeVector((currentPosition(:,1)-currentPosition(:,2)),&
               & latticeVector)

         ! Determine the square of the minimum seperation distance between the
         !   two atoms.
         atomSiteSepSqrd = sum((currentPosition(:,1) - currentPosition(:,2) - &
               & latticeVector(:))**2)

         ! Obtain the maximum distance from either atom where the overlap
         !   is considered to be non-negligable.
         currentNegligLimit = alphaDist(1,1,currentElements(1),&
               & currentElements(2))

         ! Determine if there are no alpha terms for this atom pair that fall
         !   within the current negligability limit.  Cycle if there are none.
         if (atomSiteSepSqrd > currentNegligLimit) cycle

         ! Determine the maximum radius beyond which no lattice point will be
         !   considered to contribute to the overlap integral for this atom
         !   pair.  (This is the law of cosines c^2 = a^2 + b^2 + 2ab*cos(g)
         !   with g = gamma = angle between a and b.  Here atomSiteSepSqrd=a^2,
         !   currentNegligLimit = b^2,g=0.)
         maxLatticeRadius = atomSiteSepSqrd + currentNegligLimit + 2.0_double*&
                          & sqrt(atomSiteSepSqrd * currentNegligLimit)

         ! Check if we have to check more lattice points than were determined
         !   to be needed by comparing the maxlatticeRadius to the maximum
         !   radius of the lattice points defined earlier.
         if (maxLatticeRadius > cellSizesReal(numCellsReal)) then
            write (20,*) 'More lattice points needed for this atom overlap pair'
            write (20,*) 'maxLatticeRadius requested=',maxLatticeRadius
            write (20,*) 'max available=',cellSizesReal(numCellsReal)
            stop
         endif

         ! Initialize the matrix for this atom pair that will record the
         !   overlap integrals with bloch vector (kpoint) phase factors.
#ifndef GAMMA
         currentPair(:,:,:) = 0.0_double
#else
         currentPairGamma(:,:) = 0.0_double
#endif

         ! Begin a loop over all the lattice points to shift the position of
         !   atom number 2 to all the replicated cells.
         do k = 1, numCellsReal

            ! Exit the loop when we have exceeded the necessary number of
            !   lattice points based on distance.
            if (cellSizesReal(k) > maxLatticeRadius) exit

            ! For the electronic potential we have one more test to prevent
            !   other potential terms of the current potential type from
            !   passing this point.
            if (checkElecPotInteraction(i,j,k) == 1) then
               cycle
            else
               ! Assume now that this alpha will not have any interactions
               !   with this cell/atom/atom set.
               elecPotInteraction = 0
            endif

            ! Obtain the position of atom #2 shifted by the current lattice.
            shiftedAtomPos(:) = currentPosition(:,2) + latticeVector(:) + &
                  & cellDimsReal(:,k)

            ! Obtain the square of the magnitude of the seperation vector
            !   between atom 1 and the shifted position of atom 2.
            shiftedAtomSiteSep = sum ((currentPosition(:,1) - &
                  & shiftedAtomPos(:))**2)

            ! Determine if this shifted atom position puts it outside of the
            !   above determined negligability limit for this atom pair.
            if (shiftedAtomSiteSep > currentNegligLimit) cycle


            ! Initialize a matrix to hold the product of the overlap integrals
            ! times the atom2 basis functions.
            pairXBasisFn2(:,:currentNumAlphas(1),&
                  &:currentNumTotalStates(2)) = 0.0_double

            ! Initialize a variable to track the largest atomic alpha used
            !   from atom 1.
            maxAlpha1Used = 0

            ! Now we start loops over the alphas of each atomic type.  For each
            !   alpha pair the overlap is tested to determine if it is
            !   negligable or not.  The alphas for atom 1 are tested in turn
            !   with the first alpha of atom 2.  After all the atom 1 alphas
            !   have been tested or the test fails we test all the alphas
            !   of atom 2 with the first alpha of atom 1 until all the alphas
            !   of atom 2 have been used or the test fails.  Then we increment
            !   the loop counter by 1.  The second iteration of the loop
            !   obviously works just like the first except that the only
            !   difference is that instead of starting with the first alphas
            !   of atom 1 and atom 2, we always start with the second alphas
            !   of both atoms.

            ! The loop must be up to the smaller number of alphas.
            do l = 1, min(currentNumAlphas(1),currentNumAlphas(2))

               ! Set the mode to zero to say that we are incrementing through
               !   neither the atom 1 alpha array nor the atom 2 alpha array.
               currentMode = 0

               ! Initialize the matrix index values for the two alphas.
               alphaIndex(:) = l

               ! Check if this alpha pair has any non-negligable contribution.
               if (alphaDist(l,l,currentElements(1),currentElements(2)) < &
                     & shiftedAtomSiteSep) exit

               ! Start looping over atomic alphas looking for a negligable
               !   contribution for each alpha pair.
               do while (.true.)

                  ! Check if we are going through atom 1 alphas.
                  if (currentMode == 1) then

                     ! Go to the next atom 1 alpha
                     alphaIndex(1) = alphaIndex(1) + 1

                     ! Check if there are no alphas left to do for atom 1
                     if (alphaIndex(1) > currentNumAlphas(1)) then

                        ! Switch mode to go through atom 2 alphas on the
                        !   current diagonal alpha of atom 1.
                        currentMode = 2
                        alphaIndex(1) = l
                        cycle
                     endif

                  ! Check if we are going through atom 2 alphas.
                  elseif (currentMode == 2) then

                     ! Go to the next atom 2 alpha
                     alphaIndex(2) = alphaIndex(2) + 1

                     ! Check if there are no alphas left to do for atom 2
                     if (alphaIndex(2) > currentNumAlphas(2)) exit
                  endif

                  ! Check if this atom alpha pair has any non negligable
                  !   overlap.
                  if (alphaDist(alphaIndex(1),alphaIndex(2),currentElements(1),&
                        & currentElements(2)) < shiftedAtomSiteSep) then

                     ! Switch mode to go through atom 2 alphas on the current
                     !   diagonal alpha of atom 1.
                     if (currentMode == 1) then
                        currentMode = 2
                        alphaIndex(1) = l
                        cycle
                     else
                        exit
                     endif
                  endif

                  ! At this point a sufficient overlap between atomic Gaussians
                  !   has been shown to exist.
                  contrib = .true.

                  ! We can proceed with the next step of the calculation. In
                  !   this case a further loop over all the sites for the
                  !   current term of the electronic potential will be computed.
                  call electronicPE (contrib,alphaIndex,currentElements,&
                        & currentlmAlphaIndex,shiftedAtomSiteSep,&
                        & currentPosition,shiftedAtomPos,oneAlphaPair,&
                        & currentAlphas)

                  ! Collect the results of the overlap of the current alpha
                  !   times the basis functions of atom 2.
                  if (contrib .eqv. .true.) then
                     do m = 1, currentNumTotalStates(2)
                        pairXBasisFn2(:currentlmAlphaIndex(alphaIndex(1),1),&
                              & alphaIndex(1),m) = &
                              & pairXBasisFn2(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),alphaIndex(1),m) + &
                              & oneAlphaPair(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),currentlmIndex(m,2)) * &
                              & currentBasisFns(alphaIndex(2),m,2)
                     enddo

                     ! Update the maximum alpha used from atom 1.
                     maxAlpha1Used = max(alphaIndex(1),maxAlpha1Used)
                  endif

                  ! Switch mode from the initial state of no searching to the
                  !   state of searching along alphas from atom 1.
                  if (currentMode == 0) then
                     currentMode = 1
                  endif
               enddo
            enddo  ! min number of alphas between the two atoms l

            ! At this point all the alpha loops are complete and we can form a
            !   product with the atom 1 basis function to give the overlap
            !   integral in a complete basis representation.
            call multWithBasisFn1 (currentBasisFns,pairXBasisFn2, &
                  & pairXBasisFn12,currentlmIndex,currentNumTotalStates, &
                  & maxAlpha1Used)

            ! Collect this atom 1, atom 2 basis function overlap matrix for all
            !   bloch vectors (kpoints) with phase factors appropriate to the
            !   current atom 2 lattice vector.  NOTE that the k index is over
            !   the number of cells in the superlattice.
#ifndef GAMMA
            call applyPhaseFactors (currentPair,pairXBasisFn12(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),k,0,0)
#else
            call applyPhaseFactorsGamma (currentPairGamma,pairXBasisFn12(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),k,0)
#endif
            ! If there were no interactions for this cell then mark the
            !   appropriate bit to make sure that the next atomic alpha does
            !   not try to compute anything for this cell-atom pair.
            if (elecPotInteraction == 0) then
               anyElecPotInteraction(i,j,matrixIndex) = &
                     & ibclr(anyElecPotInteraction(i,j,matrixIndex),cellIndex)
            endif

         enddo !(k superlattice)

         ! At this point all the lattice sums for the current atom pair are
         !   complete.

         ! So now we can arrange the data from this atom into a set of three
         !   large matrices.  A valence-valence matrix, a core-valence matrix,
         !   and a core-core matrix.

#ifndef GAMMA
         ! First we must make a correction for the atom 2 lattice origin shift.
         call kPointLatticeOriginShift (currentNumTotalStates,currentPair,&
               & latticeVector,numKPoints,0)
         call saveCurrentPair(i,j,numKPoints,currentPair,valeVale(:,:,:,1),&
               & coreVale,coreCore)
#else
         call saveCurrentPairGamma(i,j,currentPairGamma,&
               & valeValeGamma(:,:,1),coreValeGamma,coreCoreGamma)
#endif

      enddo ! (Atom loop #2)
   enddo    ! (Atom loop #1)

   ! At this point all atom pairs have been computed.  Now, we compute the
   !   valeVale with the effects of the core orthogonalization.  We can also
   !   deallocate a bunch of arrays and matrices that are no longer needed.
   deallocate (currentBasisFns)
   deallocate (currentAlphas)
   deallocate (currentlmAlphaIndex)
   deallocate (currentlmIndex)
   deallocate (pairXBasisFn2)
   deallocate (pairXBasisFn12)
#ifndef GAMMA
   deallocate (currentPair)
#else
   deallocate (currentPairGamma)
#endif

   ! Perform orthogonalization and save the results to disk.  The 4 is an
   !   operation code signifying that a non-overlap orthogonalization should be
   !   done, and specifically that the result is for the electronic potential
   !   and that it should be written to the EP portion of the hdf5 file.
   call ortho(4)

end subroutine gaussOverlapEP


subroutine elecPotGaussOverlap

   ! Import the necessary modules
   use O_Kinds
   use O_TimeStamps
   use HDF5
   use O_Lattice, only: numCellsReal
   use O_AtomicSites, only: numAtomSites
   use O_PotTypes, only: potTypes
   use O_PotSites, only: potSites, numPotSites

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables that are extracted from the passed data structures
   integer :: numAlphas

   ! Define local variables for logging and loop control
   integer :: i,j ! Loop index variables
   integer :: currentIterCount
   integer (hsize_t) :: matrixSize ! Used to define the size of the
         ! anyElecPotInteraction matrix and the number of bits in hsize_t.

   ! Make a time stamp.
   call timeStampStart (11)

   ! Initialize the number of bits in the hsize_t variable.
   hsize_t_bitsize = bit_size(matrixSize)

   ! Allocate space to track whether or not an atom pair has any potential
   !   interactions.  The last index must be sufficiently large so that the
   !   sum of the bits from each numAtom X numAtom block is >= the number of
   !   real lattice cells.  (We need one bit for each cell for each atom pair.)
   if (mod(numCellsReal,bit_size(matrixSize)) == 0) then
      allocate (anyElecPotInteraction (numAtomSites,numAtomSites,numCellsReal/ &
            & bit_size(matrixSize)))
   else
      allocate (anyElecPotInteraction (numAtomSites,numAtomSites,numCellsReal/ &
            & bit_size(matrixSize)+1))
   endif

   ! Initialize the counter for the total number of iterations of the loops.
   currentIterCount = 0

   do i = 1, numPotSites

      ! Check if this potential site is the first to have some particular
      !   potential type.  If it is the first site to have a particular
      !   potential type then its "firstPotType" flag is 1.  Other sites with
      !   the same type have a "firstPotType" value equal to 0.
      if (potSites(i)%firstPotType == 0) cycle

      currPotTypeNumber = potSites(i)%potTypeAssn
      numAlphas         = potTypes(currPotTypeNumber)%numAlphas

      ! Initialize the anyElecPotInteraction matrix for this set of potential
      !   alphas.  The assumption is that every atom pair will have some
      !   potential interaction in every cell.  When a cell/atom/atom set is
      !   found to not have an interaction, then that bit is set to zero.
      !   (Note that the value of -1 will have every bit set to TRUE (1).)
      !   (I hope that there is no machine dependence to this.  ^_^)
      anyElecPotInteraction(:,:,:) = -1

      do j = 1, numAlphas

         ! Increment the counter of the number of total iterations.
         currentIterCount = currentIterCount + 1

         ! Record the current parameters for this iteration.
         currPotAlpha   = potTypes(currPotTypeNumber)%alphas(j)
         currPotNumber  = i
         currPotElement = potTypes(currPotTypeNumber)%elementID
         currAlphaNumber = j
         currMultiplicity = potTypes(currPotTypeNumber)%multiplicity

         ! Calculate the overlap for the current alpha of the current type
         !   with all the alphas of every pair of atoms.
         call gaussOverlapEP

         ! Record that this loop has finished
         if (mod(currentIterCount,10) .eq. 0) then
            write (20,ADVANCE="NO",FMT="(a1)") "|"
         else
            write (20,ADVANCE="NO",FMT="(a1)") "."
         endif
         if (mod(currentIterCount,50) .eq. 0) then
            write (20,*) " ",currentIterCount
         endif
         call flush (20)
      enddo
   enddo

   ! Deallocate matrices that are no longer used.
   deallocate (anyElecPotInteraction)


   ! Make a finishing time stamp.
   call timeStampEnd (11)

end subroutine elecPotGaussOverlap


function checkElecPotInteraction (i,j,k)

   ! Make sure no funny variables are used.
   implicit none

   ! Define passed dummy variables.
   integer, intent (in) :: i ! Index for atom loop 1.
   integer, intent (in) :: j ! Index for atom loop 2.
   integer, intent (in) :: k ! Index for current superlattice cell.

   ! Define the function return variable.
   integer :: checkElecPotInteraction

   ! Determine the matrix index for this cell.  The matrix
   !   index is the third index of the anyElecPotInteraction matrix.  It
   !   is used basically to increase the number of bits for each
   !   atom/atom pair where each bit represents one of the real space
   !   superlattice cells.  The +1 on the end is to make the index
   !   count start from 1 instead of 0.  k-1 must be used because
   !   when k=64 we want the matrix index to still be 1.  Only when
   !   k=65 should the matrix index become 2.
   matrixIndex = (k-1)/hsize_t_bitsize+1

   ! Determine the cell index for this cell.  Note that the result
   !   of this mod call is a number from 0 to 63 (for 8 byte
   !   integers).  This number will be the bit number for the current
   !   matrix index atom/atom set.  Note also that the bit numbers
   !   are indexed from 0 to 63 (for 8 byte integers.)  The k-1 is
   !   used for the same reason as above.  The real space cells (k)
   !   count from 1 to numCellsReal, but the cell index must start
   !   from 0 when k=1 and be equal to 63 when k=64.
   cellIndex = mod(k-1,hsize_t_bitsize)

   ! Only proceed with this pair of atoms if there was an interaction
   !   for the previous alpha (or this is the first alpha of the set).
   if (btest(anyElecPotInteraction(i,j,matrixIndex),cellIndex) &
         & .eqv. .false.) then
      checkElecPotInteraction = 1
   else
      checkElecPotInteraction = 0
   endif

end function checkElecPotInteraction


subroutine nuclearPE(contrib,alphaIndex,currentElements,currentlmAlphaIndex,&
      & shiftedAtomSiteSep,currentPosition,shiftedAtomPos,oneAlphaPair,&
      & currentAlphas)

   ! Use necessary modules
   use O_Kinds
   use O_Constants, only: dim3, smallThresh
   use O_PotSites, only: numPotSites
   use O_GaussianRelations, only: alphaDist, alphaCenter, alphaNucDist
   use O_Lattice, only: numCellsReal, cellSizesReal, cellDimsReal, &
         & findLatticeVector
   use O_Basis, only: initializePotSite
   use O_GaussianIntegrals, only: nucPotInteg

   ! Make sure no funny variables are used.
   implicit none

   ! Define passed dummy variables.
   logical, intent (out) :: contrib
   integer, dimension(2),    intent (in) :: alphaIndex
   integer, dimension(2),    intent (in) :: currentElements
   integer, dimension (:,:), intent (in) :: currentlmAlphaIndex
   real (kind=double),                     intent (in)  :: shiftedAtomSiteSep
   real (kind=double), dimension (dim3,2), intent (in)  :: currentPosition
   real (kind=double), dimension (dim3),   intent (in)  :: shiftedAtomPos
   real (kind=double), dimension (16,16),  intent (out) :: oneAlphaPair
   real (kind=double), dimension (:,:),    intent (in)  :: currentAlphas

   ! Declare local varibales.
   integer :: m,n
   integer :: potType ! The type number of the potential at the current site.
   integer :: potElement ! The element ID number of the potential sites.
   real (kind=double) :: zFactor  ! Charge on the current nucleus.
   real (kind=double) :: nucAlpha ! Alpha exponential factor for this nucleus.
   real (kind=double) :: threeAlphaDist ! A value like the alphaDist that
         ! can be pre-calculated since part of it will not change with the
         ! actual distance, but only with the choice of alphas (which are known
         ! at input time.)
   real (kind=double), dimension (dim3) :: shiftedPotPos ! The position of
         ! the nuclear pot shifted to each relevant lattice point.
   real (kind=double), dimension (dim3) :: overlapCenter ! The position of the
         ! center of the overlap between two alphas.  This is used to help get
         ! the distance to the third center (nuclear potential).
   real (kind=double), dimension (dim3) :: centerOriginVect ! The seperation
         ! vector between the center of the two alpha overlaps, and the origin
         ! to be used for the nuclear potential site.
   real (kind=double), dimension (dim3) :: potPosition ! Original cell
         ! position of the current potential site.
   real (kind=double), dimension (dim3) :: latticeVector2 ! Vector to lattice
         ! point closest to the difference between the unit cell positions for
         ! the center of the alpha1 alpha2 overlap and alpha3 origin.
   real (kind=double) :: centerOriginSep ! The sum of squares of the above
         ! vector.
   real (kind=double) :: shiftedCenterOriginSep ! The above value shifted
         ! according to the current lattice point loop.
   real (kind=double) :: maxLatticeRadius2 ! Maximum radius beyond which no
         ! lattice points will be considered for integration.  This is used
         ! for the nuclear alpha triple.
   real (kind=double), dimension (16,16) :: nucPotAlphaOverlap

   integer :: l1l2Switch
   integer, dimension(16) :: powerOfTwo = (/0,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3/)

   ! At this point a sufficient overlap has been found for the current alpha
   !   pair so we start looking through nuclear potentials.

   ! First we find the center of the overlap between the two alphas.
   overlapCenter = alphaCenter(alphaIndex(1),alphaIndex(2),&
         & currentElements(1),currentElements(2)) * &
         & (currentPosition(:,1) - shiftedAtomPos(:)) + shiftedAtomPos(:)

   ! Initialize the result matrix to zero before starting the nuc pot loop.
   nucPotAlphaOverlap(:currentlmAlphaIndex(alphaIndex(1),1), &
         & :currentlmAlphaIndex(alphaIndex(2),2)) = 0.0_double

   ! Mark a flag as false.  This flag is set true if any nuclear potential
   !   causes a contribution.
   contrib=.false.

   ! Initiate a loop over each potential site.
   do m = 1, numPotSites

      ! Initialize the current nuclear potential
      call initializePotSite (m, potType, potElement, zFactor, nucAlpha, &
            & potPosition)

      ! If the zFactor is sufficiently small we consider the effect of the
      !   overlap to be negligable and we cycle to the next one.
      if (zFactor < smallThresh) cycle

      ! Determine the maximum distance beyond which the overlap of the three
      !   gaussians is considered negligable.
      threeAlphaDist = (1 - shiftedAtomSiteSep / &
            & alphaDist(alphaIndex(1),alphaIndex(2),&
            & currentElements(1),currentElements(2))) * &
            & alphaNucDist(alphaIndex(1),alphaIndex(2), &
            & potElement,currentElements(1),currentElements(2))

      ! Locate the origin for the potential lattice sum.
      call findLatticeVector((overlapCenter(:)-potPosition(:)),latticeVector2)

      ! Find the seperation vector and distance between the minimum overlap
      !   center, and the origin.
      centerOriginVect(:) = overlapCenter(:) - potPosition(:) - latticeVector2
      centerOriginSep = sum(centerOriginVect(:)**2)

      ! Check if the largest potential is negligable or not.
      if (centerOriginSep > threeAlphaDist) cycle

      ! At this point it must be the case that at least one contribution will
      !   be calculated.
      contrib = .true.

      ! First, find the cut-off radius for the potential lattice summation by
      !   the triangle inequality.
      maxLatticeRadius2 = centerOriginSep +  threeAlphaDist + &
            & 2.0_double * sqrt (centerOriginSep*threeAlphaDist)

      ! Begin a loop over lattice points for the nuclear pot.
      do n = 1, numCellsReal

         ! Check if this lattice point extends beyond the range of negligability
         if (cellSizesReal(n) > maxLatticeRadius2) exit

         ! Get the overlap center shifted by the lattice point.
         shiftedCenterOriginSep = sum((centerOriginVect(:) - &
               & cellDimsReal(:,n))**2)

         ! Check if the shifted seperation between the center and the origin
         !   extends beyond the negligability limit.  If so, cycle to next cell.
         if (shiftedCenterOriginSep > threeAlphaDist) cycle

         ! Get the seperation shifted by the lattice point.
         shiftedPotPos(:) = potPosition(:) + latticeVector2(:) + &
               & cellDimsReal(:,n)

         ! Calculate the opcode to do the correct set of integrals
         ! for the current alpha pair
         l1l2Switch = ishft(1,&
           &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
           &+ ishft(16,&
           &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))

         call nucPotInteg (currentAlphas(alphaIndex(1),1),&
               & currentAlphas(alphaIndex(2),2),&
               & nucAlpha, currentPosition(:,1),&
               & shiftedAtomPos(:), shiftedPotPos(:),&
               & l1l2Switch, oneAlphaPair)

         !call nucPotInteg (oneAlphaPair,&
         !      & currentlmAlphaIndex (alphaIndex(1),1),&
         !      & currentlmAlphaIndex (alphaIndex(2),2),&
         !      & currentAlphas(alphaIndex(1),1),&
         !      & currentAlphas(alphaIndex(2),2),&
         !      & nucAlpha,currentPosition(:,1),&
         !      & shiftedAtomPos(:),shiftedPotPos(:))

         ! Accumulate the results returned for this alpha set.
         nucPotAlphaOverlap(:currentlmAlphaIndex (alphaIndex(1),1), &
               & :currentlmAlphaIndex (alphaIndex(2),2)) = &
               & nucPotAlphaOverlap(:currentlmAlphaIndex (alphaIndex(1),1),&
               & :currentlmAlphaIndex (alphaIndex(2),2)) - &
               & oneAlphaPair(:currentlmAlphaIndex (alphaIndex(1),1),&
               & :currentlmAlphaIndex (alphaIndex(2),2)) * zFactor

      enddo ! (n numCells)
   enddo ! (m numPots)

   ! Store the accumulated result in the matrix that is used in the common
   !   subroutine for SCF integrals.
   oneAlphaPair(:,:) = nucPotAlphaOverlap(:,:)

end subroutine nuclearPE


subroutine electronicPE(contrib,alphaIndex,currentElements,currentlmAlphaIndex,&
      & shiftedAtomSiteSep,currentPosition,shiftedAtomPos,oneAlphaPair,&
      & currentAlphas)

   ! Use necessary modules
   use O_Kinds
   use O_Constants, only: dim3
   use O_PotSites, only: potSites
   use O_GaussianRelations, only: alphaDist, alphaCenter, alphaPotDist
   use O_Lattice, only: numCellsReal, cellSizesReal, cellDimsReal, &
         & findLatticeVector
   use O_GaussianIntegrals, only: threeCentInteg

   ! Make sure no funny variables are used.
   implicit none

   ! Define passed dummy variables.
   logical, intent (out) :: contrib
   integer, dimension(2),    intent (in) :: alphaIndex
   integer, dimension(2),    intent (in) :: currentElements
   integer, dimension (:,:), intent (in) :: currentlmAlphaIndex
   real (kind=double),                     intent (in)  :: shiftedAtomSiteSep
   real (kind=double), dimension (dim3,2), intent (in)  :: currentPosition
   real (kind=double), dimension (dim3),   intent (in)  :: shiftedAtomPos
   real (kind=double), dimension (16,16),  intent (out) :: oneAlphaPair
   real (kind=double), dimension (:,:),    intent (in)  :: currentAlphas
   ! Define variables for gauss integrals
   integer :: l1l2Switch
   integer, dimension(16) :: powerOfTwo = (/0,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3/)

   ! Declare local varibales.
   integer :: m,n
   real (kind=double) :: threeAlphaDist ! A value like the alphaDist that
         ! can be pre-calculated since part of it will not change with the
         ! actual distance, but only with the choice of alphas (which are known
         ! at input time.)
   real (kind=double), dimension (dim3) :: shiftedPotPos ! The position of
         ! the nuclear pot shifted to each relevant lattice point.
   real (kind=double), dimension (dim3) :: overlapCenter ! The position of the
         ! center of the overlap between two alphas.  This is used to help get
         ! the distance to the third center (nuclear potential).
   real (kind=double), dimension (dim3) :: centerOriginVect ! The seperation
         ! vector between the center of the two alpha overlaps, and the origin
         ! to be used for the nuclear potential site.
   real (kind=double), dimension (dim3) :: potPosition ! Original cell
         ! position of the current potential site.
   real (kind=double), dimension (dim3) :: latticeVector2 ! Vector to lattice
         ! point closest to the difference between the unit cell positions for
         ! the center of the alpha1 alpha2 overlap and alpha3 origin.
   real (kind=double) :: centerOriginSep ! The sum of squares of the above
         ! vector.
   real (kind=double) :: shiftedCenterOriginSep ! The above value shifted
         ! according to the current lattice point loop.
   real (kind=double) :: maxLatticeRadius2 ! Maximum radius beyond which no
         ! lattice points will be considered for integration.  This is used
         ! for the nuclear alpha triple.
   real (kind=double), dimension (16,16) :: elecPotAlphaOverlap

   ! At this point a sufficient overlap has been found for the current alpha
   !   pair so we start looking through equivalent potential sites.

   ! First we find the center of the overlap between the two alphas.
   overlapCenter = alphaCenter(alphaIndex(1),alphaIndex(2),&
         & currentElements(1),currentElements(2)) * &
         & (currentPosition(:,1) - shiftedAtomPos(:)) + shiftedAtomPos(:)

   ! Determine the maximum distance beyond which the overlap of the three
   !   gaussians is considered negligable.
   threeAlphaDist = (1 - shiftedAtomSiteSep / &
         & alphaDist(alphaIndex(1),alphaIndex(2),&
         & currentElements(1),currentElements(2))) * &
         & alphaPotDist(alphaIndex(1),alphaIndex(2), &
         & currAlphaNumber,currentElements(1),currentElements(2),&
         & currPotElement)

   ! Initialize the result matrix to zero before starting the atomic potential
   !   loop.
   elecPotAlphaOverlap(:currentlmAlphaIndex(alphaIndex(1),1), &
         & :currentlmAlphaIndex(alphaIndex(2),2)) = 0.0_double

   ! Mark a flag as false.  This flag is set true if any atomic potential
   !   causes a contribution.
   contrib=.false.

   ! Initiate a loop over each equivalent potential site.
   do m = 0, currMultiplicity - 1

      ! Get the position of the current potential site.
      potPosition(:) = potSites(currPotNumber+m)%cartPos(:)

      ! Locate the origin for the potential lattice sum.
      call findLatticeVector((overlapCenter(:)-potPosition(:)),latticeVector2)

      ! Find the seperation vector and distance between the minimum overlap
      !   center, and the origin.
      centerOriginVect(:) = overlapCenter(:)-potPosition(:)-latticeVector2(:)
      centerOriginSep = sum(centerOriginVect(:)**2)

      ! Check if the current pot is negligable or not in even the closest case.
      if (centerOriginSep > threeAlphaDist) cycle

      ! At this point it must be the case that at least one contribution will
      !   be calculated.
      contrib=.true.
      elecPotInteraction=1

      ! First, find the cut-off radius for the potential lattice
      !   summation by the triangle inequality.
      maxLatticeRadius2 = centerOriginSep + threeAlphaDist + &
            & 2.0_double * sqrt (centerOriginSep*threeAlphaDist)

      ! Begin a loop over lattice points for the atom pot.
      do n = 1, numCellsReal

         ! Check if this lattice point extends beyond the negligability range.
         if (cellSizesReal(n) > maxLatticeRadius2) exit

         ! Get the overlap center shifted by the lattice point.
         shiftedCenterOriginSep = sum((centerOriginVect(:) - &
               & cellDimsReal(:,n))**2)

         ! Check if the shifted seperation between center and origin extends
         !   beyond the negligability limit.  If so, cycle to the next cell.
         if (shiftedCenterOriginSep > threeAlphaDist) cycle

         ! Get the seperation shifted by the lattice point.
         shiftedPotPos(:) = potPosition(:) + latticeVector2(:) + &
               & cellDimsReal(:,n)

         !call threeCentInteg (currentAlphas(alphaIndex(1),1),&
         !& currentAlphas(alphaIndex(2),2),currPotAlpha,&
         !& currentPosition(:,1), shiftedAtomPos(:), &
         !& shiftedPotPos(:), oneAlphaPair, l1l2Switch)

         ! Calculate the opcode to do the correct set of integrals
         ! for the current alpha pair
         l1l2Switch = ishft(1,&
               &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
               &+ ishft(16,&                     
               &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))

         call threeCentInteg (currentAlphas(alphaIndex(1),1),&
               & currentAlphas(alphaIndex(2),2),currPotAlpha,&
               & currentPosition(:,1), shiftedAtomPos(:), &
               & shiftedPotPos(:), l1l2Switch, oneAlphaPair)

         !call threeCentInteg (oneAlphaPair,&
         !      & currentlmAlphaIndex (alphaIndex(1),1),&
         !      & currentlmAlphaIndex (alphaIndex(2),2),&
         !      & currentAlphas(alphaIndex(1),1),&
         !      & currentAlphas(alphaIndex(2),2),currPotAlpha,&
         !      & currentPosition(:,1),shiftedAtomPos(:),shiftedPotPos(:))

         ! Accumulate the results returned for this alpha set.
         elecPotAlphaOverlap(:currentlmAlphaIndex(alphaIndex(1),1), &
               & :currentlmAlphaIndex (alphaIndex(2),2)) = &
               & elecPotAlphaOverlap(:currentlmAlphaIndex(alphaIndex(1),1), &
               & :currentlmAlphaIndex(alphaIndex(2),2)) + &
               & oneAlphaPair(:currentlmAlphaIndex(alphaIndex(1),1), &
               & :currentlmAlphaIndex(alphaIndex(2),2))

      enddo ! (n numCells)
   enddo ! (m multiplicity)

   ! Store the accumulated result in the matrix that is used in the common
   !   subroutine for SCF integrals.
   oneAlphaPair(:,:) = elecPotAlphaOverlap(:,:)

end subroutine electronicPE


subroutine orthoOL

   ! Use necessary modules.
   use O_Kinds
   use O_AtomicSites, only: coreDim, valeDim
   use O_KPoints, only: numKPoints
   use O_SetupIntegralsHDF5, only: atomOverlap_did, atomDims
   use O_Orthogonalization

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define local variables.
   integer :: i,j,k
   integer :: hdferr
   integer :: currIndex
   real (kind=double), allocatable, dimension (:,:) :: packedValeVale

   ! Allocate space for the valeCore matrix in a pre-transposed format for
   !   more efficient cache utilization during the matrix multiplication.  Also
   !   allocate space to store the packed valeVale matrix.
#ifndef GAMMA
   allocate (valeCore (coreDim,valeDim))
   allocate (packedValeVale (2,valeDim*(valeDim+1)/2))
#else
   allocate (valeCoreGamma  (coreDim,valeDim))
   allocate (packedValeVale (1,valeDim*(valeDim+1)/2))
#endif

   do i = 1, numKPoints

      ! It is only necessary to perform orthogonalization when the core
      !   dimension is non-zero.
      if (coreDim /= 0) then

#ifndef GAMMA
         ! Form a product of (valeCore)(coreValeOL) for the overlap
         !   integral (valeVale).
         call valeCoreCoreValeOL (valeDim,coreDim,valeVale(:,:,i,1),&
               & coreValeOL(:,:,i))

         ! Form product of (coreValeOL)(coreCore) in temp matrix (valeCore).
         call coreValeCoreCore (valeDim,coreDim,valeCore(:,:),&
               & coreValeOL(:,:,i),coreCore(:,:,i))

         ! Finally compute the product of the above
         !   (valeCore)(coreCore) with coreValeOL.
         call makeValeVale (valeDim,coreDim,valeDim,valeCore(:,:),&
               & coreValeOL(:,:,i),valeVale(:,:,i,1),packedValeVale,1,0)
#else
         ! Form a product of (valeCore)(coreValeOL) for the overlap
         !   integral (valeVale).
         call valeCoreCoreValeOLGamma (valeDim,coreDim,&
               & valeValeGamma(:,:,1),coreValeOLGamma)

         ! Form product of (coreValeOL)(coreCore) in temp matrix (valeCore).
         call coreValeCoreCoreGamma (valeDim,coreDim,valeCoreGamma,&
               & coreValeOLGamma,coreCoreGamma)

         ! Finally compute the product of the above
         !   (valeCore)(coreCore) with coreValeOL.
         call makeValeValeGamma (valeDim,coreDim,valeDim,valeCoreGamma,&
               & coreValeOLGamma,valeValeGamma(:,:,1),packedValeVale,1,0)
#endif
      else

         ! Initialize the index counter for packing.
         currIndex = 0

         ! Pack the valeVale matrix.
#ifndef GAMMA
         do j = 1, valeDim
            do k = 1, j
               currIndex = currIndex + 1
               packedValeVale(1,currIndex) = &
                     & real(valeVale(k,j,i,1),double)
               packedValeVale(2,currIndex) = aimag(valeVale(k,j,i,1))
            enddo
         enddo
#else
         do j = 1, valeDim
            do k = 1, j
               currIndex = currIndex + 1
               packedValeVale(1,currIndex) = valeValeGamma(k,j,1)
            enddo
         enddo
#endif
      endif

      ! Write the overlap valeVale data onto disk in HDF5 format.
      call h5dwrite_f(atomOverlap_did(i),H5T_NATIVE_DOUBLE,&
            & packedValeVale(:,:),atomDims,hdferr)
      if (hdferr /= 0) stop 'Can not write atom overlap.'
   enddo

   ! Deallocate remaining unnecessary matrices
#ifndef GAMMA
   deallocate (valeCore)
   deallocate (packedValeVale)
#else
   deallocate (valeCoreGamma)
   deallocate (packedValeVale)
#endif

end subroutine orthoOL


subroutine ortho (opCode)

   ! Use necessary modules.
   use O_Kinds
   use O_PotTypes, only: potTypes
   use O_KPoints, only: numKPoints
   use O_AtomicSites, only: coreDim, valeDim
   use O_SetupIntegralsHDF5, only: atomKEOverlap_did, atomNucOverlap_did, &
         & atomPotOverlap_did, atomDims
   use O_Orthogonalization

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed dummy arguments.
   integer :: opCode

   ! Define local variables.
   integer :: i,j,k
   integer :: hdferr
   integer :: currIndex
   real (kind=double), allocatable, dimension (:,:) :: packedValeVale

   ! Orthogonalizing against the overlap matrix is unnecessary when the core
   !   dimension is zero.  However, we must still allocate some matrices for
   !   use later.
   if (coreDim /= 0) then
#ifndef GAMMA
      do i = 1, numKPoints
        ! Form product of (valeCoreOL)(coreVale) and (valeCore)(coreValeOL).
        !   Subtract both from the kinetic energy (valeVale).
        call valeCoreCoreVale (valeDim,coreDim,valeVale(:,:,i,1),&
              & coreVale(:,:,i),coreValeOL(:,:,i))
      enddo
#else
      do i = 1, numKPoints
        ! Form a product of (valeCoreOL)(coreVale) and (valeCoreKE)
        !   (coreValeOL).  Subtract both from the kinetic energy (valeVale2).
        call valeCoreCoreValeGamma (valeDim,coreDim,valeValeGamma(:,:,1),&
              & coreValeGamma,coreValeOLGamma)
      enddo
#endif
   endif

   ! Allocate space to finish orthogonalization and pack the valeVale matrix.
#ifndef GAMMA
   deallocate (coreVale)
   allocate (valeCore(coreDim,valeDim)) ! Pre-transposed format.
   allocate (packedValeVale(2,valeDim*(valeDim+1)/2))
#else
   deallocate (coreValeGamma)
   allocate (valeCoreGamma(coreDim,valeDim)) ! Pre-transposed format.
   allocate (packedValeVale(1,valeDim*(valeDim+1)/2))
#endif

   do i = 1, numKPoints

      ! Orthogonalization is only necessary if the core dimension is non-zero.
      if (coreDim /= 0) then
#ifndef GAMMA
         ! Form product of (coreValeOL)(coreCore) in temp matrix (valeCore).
         call coreValeCoreCore (valeDim,coreDim,valeCore(:,:),&
               & coreValeOL(:,:,i),coreCore(:,:,i))

         ! Finally compute the product of the above
         !   (valeCore)(coreCore) with coreValeOL.
         call makeValeVale (valeDim,coreDim,valeDim,valeCore(:,:),&
               & coreValeOL(:,:,i),valeVale(:,:,i,1),packedValeVale(:,:),1,0)
#else
         ! Form product of (coreValeOL)(coreCore) in temp matrix (valeCore).
         call coreValeCoreCoreGamma (valeDim,coreDim,valeCoreGamma,&
               & coreValeOLGamma,coreCoreGamma)

         ! Finally compute the product of the above
         !   (valeCore)(coreCore) with coreValeOL.
         call makeValeValeGamma (valeDim,coreDim,valeDim,valeCoreGamma,&
               & coreValeOLGamma,valeValeGamma(:,:,1),&
               & packedValeVale(:,:),1,0)
#endif
      else

         ! Initialize the index for packing the valeVale matrix.
         currIndex = 0
#ifndef GAMMA
         do j = 1, valeDim
            do k = 1, j
               currIndex = currIndex + 1
               packedValeVale(1,currIndex) = &
                     & real(valeVale(k,j,i,1),double)
               packedValeVale(2,currIndex) = aimag(valeVale(k,j,i,1))
            enddo
         enddo
#else
         do j = 1, valeDim
            do k = 1, j
               currIndex = currIndex + 1
               packedValeVale(1,currIndex) = valeValeGamma(k,j,1)
            enddo
         enddo
#endif
      endif

      ! Write the valeVale term onto disk in HDF5 format.
      select case (opCode)
      case (2)
         call h5dwrite_f(atomKEOverlap_did(i),H5T_NATIVE_DOUBLE,&
               & packedValeVale(:,:),atomDims,hdferr)
         if (hdferr /= 0) stop 'Failed to write kinetic energy vale vale'
      case (3)
         call h5dwrite_f(atomNucOverlap_did(i),H5T_NATIVE_DOUBLE,&
               & packedValeVale(:,:),atomDims,hdferr)
         if (hdferr /= 0) stop 'Failed to write nuclear potential vale vale'
      case (4)
         call h5dwrite_f(atomPotOverlap_did(i,&
               & potTypes(currPotTypeNumber)%cumulAlphaSum+currAlphaNumber),&
               & H5T_NATIVE_DOUBLE,packedValeVale(:,:),atomDims,hdferr)
         if (hdferr /= 0) stop 'Failed to write electronic potential vale vale'
      end select
   enddo


   ! Deallocate matrices that are no longer necessary.
#ifndef GAMMA
   deallocate (valeCore)
   deallocate (packedValeVale)
#else
   deallocate (valeCoreGamma)
   deallocate (packedValeVale)
#endif

end subroutine ortho


subroutine cleanUpIntegralsSCF

   implicit none

#ifndef GAMMA
   deallocate (coreCore)
   deallocate (valeVale)
   deallocate (coreValeOL)
#else
   deallocate (coreCoreGamma)
   deallocate (valeValeGamma)
   deallocate (coreValeOLGamma)
#endif

end subroutine cleanUpIntegralsSCF


end module O_IntegralsSCF
