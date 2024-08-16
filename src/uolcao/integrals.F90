module O_Integrals

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
   !   orthogonalization of the nuclear pot, electronic pot, MV, and KE
   !   matrices.
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
   integer :: someIntegralDone

   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine allocateIntegralsSCF(coreDim,valeDim,numKPoints)

   implicit none

   ! Define passed dummy parameters.
   integer, intent(in) :: coreDim
   integer, intent(in) :: valeDim
   integer, intent(in) :: numKPoints

   ! Note that the "1" in the valeVale allocations is so that some other
   !   subroutines can be generalized to take valeVale variables with
   !   dimension=4 and where the fourth dimension can be a 1 or a 3.
#ifndef GAMMA
   allocate (coreCore (coreDim,coreDim,numKPoints))
   allocate (valeVale (valeDim,valeDim,numKPoints,1))
#else
   allocate (coreCoreGamma (coreDim,coreDim))
   allocate (valeValeGamma (valeDim,valeDim,1))
#endif

   ! Assume that no integrals will be performed.
   someIntegralDone = 0

end subroutine allocateIntegralsSCF



!subroutine allocateIntegralsPSCF(coreDim,valeDim,numKPoints)
!
!   implicit none
!
!   ! Define passed dummy parameters.
!   integer, intent(in) :: coreDim
!   integer, intent(in) :: valeDim
!   integer, intent(in) :: numKPoints
!
!#ifndef GAMMA
!   allocate (coreCore (coreDim,coreDim,numKPoints))
!   allocate (valeVale (valeDim,valeDim,numKPoints,1))
!#else
!   allocate (coreCoreGamma (coreDim,coreDim))
!   allocate (valeValeGamma (valeDim,valeDim,1))
!#endif
!
!end subroutine allocateIntegralsPSCF


! Standard two center overlap integral.
subroutine gaussOverlapOL

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
   use O_GaussianIntegrals, only: overlap2CIntg
   use O_SCFIntegralsHDF5, only: atomOverlap_aid
   use O_Basis, only: initializeAtomSite
   use O_IntgSaving

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables for logging and loop control
   integer :: i,j,k,l,m ! Loop index variables
   integer :: hdf5Status
   integer :: hdferr
   integer(hsize_t), dimension (1) :: attribIntDims ! Attribute dataspace dim

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
         ! from atom 1 are being compared to the smallest untested alpha of
         !  atom 2, or (2) the other way around.  Case (0) is just the initial
         !  state.
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
   call timeStampStart(8)

   ! Determine if this calculation has already been completed by a previous
   !   OLCAO execution.
   hdf5Status = 0
   attribIntDims(1) = 1
   call h5aread_f(atomOverlap_aid,H5T_NATIVE_INTEGER,hdf5Status,&
         & attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to read atom overlap status.'
   if (hdf5Status == 1) then
      write(20,*) "Two-center overlap already exists. Skipping."
      call timeStampEnd(8)
      call h5aclose_f(atomOverlap_aid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atom overlap status.'
      return
   endif

   ! At this point, we know that an integral will be performed. Note it.
   someIntegralDone = 1

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
                  !   overlap. If it does not, then cycle or exit as needed.
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

                  ! If we are in the initial state of "no searching" then
                  !   switch the mode from that state to the state of
                  !   searching along alphas from atom 1 while keeping the
                  !   atom 2 alpha fixed.
                  if (currentMode == 0) then
                     currentMode = 1
                  endif

                  ! At this point a sufficient overlap between atomic Gaussians
                  !   has been shown to exist. FIX: This appears to be a
                  !   useless term. Similarly, the .eqv. test below becomes
                  !   useless.
                  contrib = .true.

                  ! Calculate the opcode to do the correct set of integrals
                  !   for the current alpha pair
                  l1l2Switch = ishft(1,&
                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
                        &+ ishft(16,&
                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))

                  ! We can proceed with the next step of the calculation.
                  !   This is the actual integral.             
                  call overlap2CIntg (currentAlphas(alphaIndex(1),1),&
                        & currentAlphas(alphaIndex(2),2), &
                        & currentPosition(:,1), shiftedAtomPos(:),&
                        & l1l2Switch, oneAlphaPair)

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
               enddo
            enddo  ! min number of alphas between the two atoms l

            ! As long as at least *some* alpha from atom 1 was used then the
            !   pairXBasisFn2 will be non-zero and we should process it.
            if (maxAlpha1Used > 0) then

               ! At this point all the alpha loops are complete and we can form
               !   a product with the atom 1 basis function to give the overlap
               !   integral in a complete basis representation.
               call multWithBasisFn1 (currentBasisFns,pairXBasisFn2, &
                     & pairXBasisFn12,currentlmIndex,currentNumTotalStates, &
                     & maxAlpha1Used)

               ! Collect this atom 1, atom 2 basis function overlap matrix for
               !   all bloch vectors (kpoints) with phase factors appropriate
               !   to the current atom 2 lattice vector.  NOTE that the k
               !   index is over the number of cells in the superlattice.
#ifndef GAMMA
               call applyPhaseFactors (currentPair,pairXBasisFn12(1:&
                     & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                     & currentNumTotalStates(1),currentNumTotalStates(2),k,0,0)
#else
               call applyPhaseFactorsGamma (currentPairGamma,pairXBasisFn12(1:&
                     & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                     & currentNumTotalStates(1),currentNumTotalStates(2),0)
#endif
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
               & coreValeOL,coreCore,0)
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
   use O_GaussianIntegrals, only: kinetic2CIntg
   use O_SCFIntegralsHDF5, only: atomKEOverlap_aid
   use O_Basis, only: initializeAtomSite
   use O_IntgSaving

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables for logging and loop control
   integer :: i,j,k,l,m ! Loop index variables
   integer :: hdf5Status
   integer :: hdferr
   integer(hsize_t), dimension (1) :: attribIntDims ! Attribute dataspace dim

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
   integer :: currentMode ! A simple flag indicating weather (1) the alphas
         ! from atom 1 are being compared to the smallest untested alpha of
         ! atom 2, or (2) the other way around.  Case (0) is just the initial
         ! state.
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

   ! Determine if this calculation has already been completed by a previous
   !   OLCAO execution.
   hdf5Status = 0
   attribIntDims(1) = 1
   call h5aread_f(atomKEOverlap_aid,H5T_NATIVE_INTEGER,hdf5Status,&
         & attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to read atom KE overlap status.'
   if (hdf5Status == 1) then
      write(20,*) "Two-center KE overlap already exists. Skipping."
      call timeStampEnd(9)
      call h5aclose_f(atomKEOverlap_aid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atom KE overlap status.'
      return
   endif

   ! At this point, we know that an integral will be performed. Note it.
   someIntegralDone = 1

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
                  !   has been shown to exist. FIX: This appears to be a
                  !   useless term. Similarly, the .eqv. test below becomes
                  !   useless.
                  contrib = .true.

                  ! Calculate the opcode to do the correct set of integrals
                  ! for the current alpha pair
                  l1l2Switch = ishft(1,&
                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
                        &+ ishft(16,&
                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))

                  ! We can proceed with the next step of the calculation. This
                  ! is the actual integral.
                  call kinetic2CIntg (currentAlphas(alphaIndex(1),1),&
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
                  & currentNumTotalStates(1),currentNumTotalStates(2),0)
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
               & coreVale,coreCore,0)
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


! Two center mass velocity overlap integrals.
subroutine gaussOverlapMV

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
   use O_GaussianIntegrals, only: massVel2CIntg
   use O_SCFIntegralsHDF5, only: atomMVOverlap_aid
   use O_Basis, only: initializeAtomSite
   use O_IntgSaving

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables for logging and loop control
   integer :: i,j,k,l,m ! Loop index variables
   integer :: hdf5Status
   integer :: hdferr
   integer(hsize_t), dimension (1) :: attribIntDims ! Attribute dataspace dim

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
   call timeStampStart (30)

   ! Determine if this calculation has already been completed by a previous
   !   OLCAO execution.
   hdf5Status = 0
   attribIntDims(1) = 1
   call h5aread_f(atomMVOverlap_aid,H5T_NATIVE_INTEGER,hdf5Status,&
         & attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to read atom MV overlap status.'
   if (hdf5Status == 1) then
      write(20,*) "Two-center MV overlap already exists. Skipping."
      call timeStampEnd(30)
      call h5aclose_f(atomMVOverlap_aid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atom MV overlap status.'
      return
   endif

   ! At this point, we know that an integral will be performed. Note it.
   someIntegralDone = 1

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

                  ! Check if we are going through atom 1 alphas. NOTE: when
                  !   this while-loop starts currentMode==0 so that both the
                  !   if and elseif conditions will be false.
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
                  !   overlap. If it does not, then cycle or exit as needed.
                  if (alphaDist(alphaIndex(1),alphaIndex(2),currentElements(1),&
                        & currentElements(2)) < shiftedAtomSiteSep) then

                     ! If we are (were) iterating through atom 1 alphas, then
                     !   switch mode to go through atom 2 alphas on the current
                     !   diagonal alpha of atom 1. Otherwise, if we must have
                     !   already been going through the atom 2 alphas and so
                     !   we can exit.
                     if (currentMode == 1) then
                        currentMode = 2
                        alphaIndex(1) = l
                        cycle
                     else
                        exit
                     endif
                  endif

                  ! If we are in the initial state of "no searching" then
                  !   switch the mode from that state to the state of
                  !   searching along alphas from atom 1 while keeping the
                  !   atom 2 alpha fixed.
                  if (currentMode == 0) then
                     currentMode = 1
                  endif

                  ! At this point a sufficient overlap between atomic Gaussians
                  !   has been shown to exist. FIX: This appears to be a
                  !   useless term. Similarly, the .eqv. test below becomes
                  !   useless.
                  contrib = .true.

                  ! Calculate the opcode to do the correct set of integrals
                  ! for the current alpha pair
                  l1l2Switch = ishft(1,&
                        & (powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
                        & + ishft(16,&
                        & (powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))

                  ! We can proceed with the next step of the calculation. This
                  ! is the actual integral.
                  call massVel2CIntg (currentAlphas(alphaIndex(1),1),&
                        & currentAlphas(alphaIndex(2),2),&
                        & currentPosition(:,1), shiftedAtomPos(:),&
                        & l1l2Switch, oneAlphaPair)

                  ! Collect the results of the overlap of the current alpha
                  !   times the basis functions of atom 2. Note that a minus
                  !   sign is used in the accumulation because the mass
                  !   velocity term has an overall negative contribution to
                  !   the Hamiltonian.
                  if (contrib .eqv. .true.) then
                     do m = 1, currentNumTotalStates(2)
                        pairXBasisFn2(:currentlmAlphaIndex(alphaIndex(1),1),&
                              & alphaIndex(1),m) = &
                              & pairXBasisFn2(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),alphaIndex(1),m) - &
                              & oneAlphaPair(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),currentlmIndex(m,2)) * &
                              & currentBasisFns(alphaIndex(2),m,2)
                     enddo

                     ! Update the maximum alpha used from atom 1.
                     maxAlpha1Used = max(alphaIndex(1),maxAlpha1Used)
                  endif
               enddo
            enddo  ! min number of alphas between the two atoms l

            ! As long as at least *some* alpha from atom 1 was used then the
            !   pairXBasisFn2 will be non-zero and we should process it.
            if (maxAlpha1Used > 0) then

               ! At this point all the alpha loops are complete and we can
               !   form a product with the atom 1 basis function to give
               !   the overlap integral in a complete basis representation.
               call multWithBasisFn1 (currentBasisFns,pairXBasisFn2, &
                     & pairXBasisFn12,currentlmIndex,currentNumTotalStates, &
                     & maxAlpha1Used)

            ! Collect this atom 1, atom 2 basis function overlap matrix for
            !   all bloch vectors (kpoints) with phase factors appropriate to
            !   the current atom 2 lattice vector.  NOTE that the k index is
            !   over the number of cells in the superlattice.
#ifndef GAMMA
               call applyPhaseFactors (currentPair,pairXBasisFn12(1:&
                     & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                     & currentNumTotalStates(1),currentNumTotalStates(2),k,0,0)
#else
               call applyPhaseFactorsGamma (currentPairGamma,pairXBasisFn12(1:&
                     & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                     & currentNumTotalStates(1),currentNumTotalStates(2),0)
#endif
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
               & coreVale,coreCore,0)
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

   ! Perform orthogonalization and save the results to disk.  The 5 is an
   !   operation code signifying that a non-overlap orthogonalization should be
   !   done, and specifically that the result is for mass velocity and that
   !   it should be written to the MV portion of the hdf5 file.
   call ortho(5)

   ! Record the completion of this gaussian integration set.
   call timeStampEnd (30)

end subroutine gaussOverlapMV


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
   use O_SCFIntegralsHDF5, only: atomNucOverlap_aid
   use O_Basis, only: initializeAtomSite
   use O_IntgSaving

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables for logging and loop control
   integer :: i,j,k,l,m ! Loop index variables
   integer :: hdf5Status
   integer :: hdferr
   integer(hsize_t), dimension (1) :: attribIntDims ! Attribute dataspace dim

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

   ! Determine if this calculation has already been completed by a previous
   !   OLCAO execution.
   hdf5Status = 0
   attribIntDims(1) = 1
   call h5aread_f(atomNucOverlap_aid,H5T_NATIVE_INTEGER,hdf5Status,&
         & attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to read atom nuclear overlap status.'
   if (hdf5Status == 1) then
      write(20,*) "Three-center nuclear overlap already exists. Skipping."
      call timeStampEnd(10)
      call h5aclose_f(atomNucOverlap_aid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atom Nuc overlap status.'
      return
   endif

   ! At this point, we know that an integral will be performed. Note it.
   someIntegralDone = 1

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
                  !   has been shown to exist. FIX: Setting to .true. here is
                  !   useless because it is reset to .false. in nuclearPE.
                  !   Unlike in earlier subroutines, the .eqv. test below is
                  !   not useless.
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
                  & currentNumTotalStates(1),currentNumTotalStates(2),0)
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
               & coreVale,coreCore,0)
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
#ifndef GAMMA
   use O_KPoints, only: numKPoints
#endif
   use O_GaussianRelations, only: alphaDist
   use O_Lattice, only: numCellsReal, cellSizesReal, cellDimsReal, &
         & findLatticeVector
   use O_AtomicSites, only: valeDim, coreDim, numAtomSites
   use O_AtomicTypes, only: maxNumAtomAlphas, maxNumStates
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
                  !   has been shown to exist. FIX: Setting to .true. here is
                  !   useless because it is reset to .false. in electronicPE.
                  !   Unlike in earlier subroutines, the .eqv. test below is
                  !   not useless.
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
                  & currentNumTotalStates(1),currentNumTotalStates(2),0)
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
               & coreVale,coreCore,0)
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
   use O_SCFIntegralsHDF5, only: atomPotTermOL_aid

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables that are extracted from the passed data structures
   integer :: numAlphas

   ! Define local variables for logging and loop control
   integer :: i,j ! Loop index variables
   integer :: currentIterCount
   integer (hsize_t) :: matrixSize ! Used to define the size of the
         ! anyElecPotInteraction matrix and the number of bits in hsize_t.
   integer :: hdf5Status
   integer :: hdferr
   integer(hsize_t), dimension (1) :: attribIntDims ! Attribute dataspace dim

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

         ! Determine if this calculation has already been completed by a
         !   previous OLCAO execution.
         hdf5Status = 0
         attribIntDims(1) = 1
         call h5aread_f(atomPotTermOL_aid(currentIterCount),&
               & H5T_NATIVE_INTEGER,hdf5Status,attribIntDims,hdferr)
         if (hdferr /= 0) stop 'Failed to read atom pot term OL status.'
         if (hdf5Status == 1) then
            write(20,*) &
                  & "Three-center pot term OL already exists. Skipping: ", &
                  & currentIterCount
            call h5aclose_f(atomPotTermOL_aid(currentIterCount),hdferr)
            if (hdferr /= 0) stop 'Failed to close atom pot term OL status.'
            cycle
         endif

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
   use O_GaussianIntegrals, only: nuclear3CIntg

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

         call nuclear3CIntg (currentAlphas(alphaIndex(1),1),&
               & currentAlphas(alphaIndex(2),2),&
               & nucAlpha,currentPosition(:,1),&
               & shiftedAtomPos(:), shiftedPotPos(:),&
               & l1l2Switch, oneAlphaPair)

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
   use O_GaussianIntegrals, only: electron3CIntg

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

         ! Calculate the opcode to do the correct set of integrals
         ! for the current alpha pair
         l1l2Switch = ishft(1,&
               &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
               &+ ishft(16,&                     
               &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))

         call electron3CIntg (currentAlphas(alphaIndex(1),1),&
               & currentAlphas(alphaIndex(2),2),currPotAlpha,&
               & currentPosition(:,1), shiftedAtomPos(:), &
               & shiftedPotPos(:), l1l2Switch, oneAlphaPair)

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


!subroutine pscfHamiltonian
!
!   ! Import necessary modules.
!   use O_Kinds
!   use O_TimeStamps
!   use O_Constants, only: dim3
!   use O_CommandLine, only: doDIMO_PSCF, doOPTC_PSCF
!   use O_KPoints, only: numKPoints
!   use O_GaussianRelations, only: alphaDist
!   use O_AtomicSites, only: valeDim, coreDim, numAtomSites, atomSites
!   use O_AtomicTypes, only: maxNumAtomAlphas, maxNumStates, atomTypes
!   use O_Lattice, only: numCellsReal, cellSizesReal, cellDimsReal, &
!         & findLatticeVector
!   use O_GaussianIntegrals, only: overlap2CIntg
!   use O_PSCFIntegralsHDF5, only: atomOverlap_aid, atomHamOverlap_aid, &
!         & atomDMOverlap_gid, atomMMOverlap_gid
!   use O_Basis, only: initializeAtomSite
!   use O_IntgSaving
!
!   ! Make sure that there are not accidental variable declarations.
!   implicit none
!
!   ! Define local variables for logging and loop control
!   integer :: i,j,k,l,m ! Loop index variables
!   integer :: hdf5Status
!   integer :: hdferr
!   integer(hsize_t), dimension (1) :: attribIntDims ! Attribute dataspace dim
!
!   ! Atom specific variables that change with each atom pair loop iteration.
!   integer,              dimension (2)    :: currentAtomType
!   integer,              dimension (2)    :: currentElements
!   integer,              dimension (2)    :: currentNumAlphas
!   integer,              dimension (2)    :: currentNumCoreStates
!   integer,              dimension (2)    :: currentNumValeStates
!   integer,              dimension (2)    :: currentNumTotalStates
!   integer, allocatable, dimension (:,:)  :: currentlmIndex
!   integer, allocatable, dimension (:,:)  :: currentlmAlphaIndex
!   real (kind=double), dimension (dim3,2) :: currentPosition
!   real (kind=double), allocatable, dimension (:,:)   :: currentAlphas
!   real (kind=double), allocatable, dimension (:,:,:) :: currentBasisFns
!
!   ! Alpha loop variables.  Once an overlap has been determined to exist for
!   !   two atoms (including lattice shifting) a loop is initiated that goes
!   !   over the alphas of those two atoms.  These variables are used there.
!   integer :: currentMode ! A simple flag indicating weather (1) the alphas
!         ! from atom 1 are being compared to the smallest untested alpha of
!         ! atom 2, or (2) the other way around.  Case (0) is just the initial
!         ! state.
!   integer :: maxAlpha1Used ! Simply the largest alpha used from the current
!         ! atom pair.
!   integer, dimension (2) :: alphaIndex ! This tracks which alphas of the
!         ! current atom pair are being tested for overlap.
!   real (kind=double), dimension (16,16) :: oneAlphaPair ! The overlap from
!         ! one alpha pair.
!   real (kind=double), allocatable, dimension (:,:,:) :: pairXBasisFn2 ! The
!         ! above overlap times the basis function from atom 2.  This is
!         ! accumulated with each iteration of the alpha loop.
!   logical :: contrib  ! At least one alpha pair contributes.  For the nuc and
!         ! elec calculations this also requires that the pot term contribute.
!
!   ! Variables and data structures that change or are accumulated with each
!   !   iteration of the lattice loop.
!#ifndef GAMMA
!   complex (kind=double), allocatable, dimension (:,:,:) :: currentPair
!#else
!   real (kind=double), allocatable, dimension (:,:) :: currentPairGamma
!#endif
!   real (kind=double), allocatable, dimension (:,:) :: pairXBasisFn12
!
!
!   ! Local position and direction vectors and radii
!   real (kind=double), dimension (dim3) :: latticeVector ! Vector to lattice
!         ! point closest to the difference between the unit cell positions for
!         ! atom 1 and atom 2.
!   real (kind=double), dimension (dim3) :: shiftedAtomPos ! The position of
!         ! atom 2 shifted to each relevant lattice point.
!   real (kind=double) :: atomSiteSepSqrd ! The square of the minimum distance
!         ! seperating atom 1 and atom 2 according to their unit cell positions
!         ! shifted by the lattice point closest to their difference.
!   real (kind=double) :: shiftedAtomSiteSep ! The seperation distance between
!         ! atom 1 and the shifted position of atom 2.
!   real (kind=double) :: currentNegligLimit ! The distance beyond which all
!         ! alpha pairs are considered to have no overlap.
!   real (kind=double) :: maxLatticeRadius ! Maximum radius beyond which no
!         ! lattice points will be considered for integration.
!
!   ! Define variables for gauss integrals
!   integer :: l1l2Switch
!   integer, dimension(16) :: powerOfTwo = (/0,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3/)
!   integer :: ovlpDone ! Flag: the overlap was already computed.
!   integer :: hamDone ! Flag: the hamiltonian was already computed.
!
!   ! Record the beginning of this phase of the setup calculation.
!   call timeStampStart(20)
!
!   ! We need to understand which integrals were requested and which ones
!   !   have already been computed.
!   attribIntDims(1) = 1
!
!   ! Check on the overlap.
!   ovlpDone = 0
!   call h5aread_f(atomOverlap_aid,H5T_NATIVE_INTEGER,ovlpDone,&
!         & attribIntDims,hdferr)
!   if (hdferr /= 0) stop 'Failed to read atom overlap status.'
!   if (ovlpDone == 1) then
!      write(20,*) "Two-center overlap already exists. Skipping."
!      call h5aclose_f(atomOverlap_aid,hdferr)
!      if (hdferr /= 0) stop 'Failed to close atom overlap status.'
!   endif
!
!   ! Check on the hamiltonian.
!   hamDone = 0
!   call h5aread_f(atomHamOverlap_aid,H5T_NATIVE_INTEGER,hamDone,&
!         & attribIntDims,hdferr)
!   if (hdferr /= 0) stop 'Failed to read atom hamiltonian overlap status.'
!   if (hamDone == 1) then
!      write(20,*) "Two- and three-center hamiltonian already exists. Skipping."
!      call h5aclose_f(atomOverlap_aid,hdferr)
!      if (hdferr /= 0) stop 'Failed to close atom hamiltonian overlap status.'
!   endif
!
!   ! Check on the dipole moment.
!   dimoDone = 0
!   call h5aread_f(atomDMOverlap_aid,H5T_NATIVE_INTEGER,dimoDone,&
!         & attribIntDims,hdferr)
!   if (hdferr /= 0) stop 'Failed to read atom dipole moment overlap status.'
!   if (dimoDone == 1) then
!      write(20,*) "Two-center-ish dipole moment already exists. Skipping."
!      call h5aclose_f(atomDMOverlap_aid,hdferr)
!      if (hdferr /= 0) stop 'Failed to close dipole moment status.'
!   endif
!
!   ! Check on the momentum matrix elements.
!   MMDone = 0
!   call h5aread_f(atomMMOverlap_aid,H5T_NATIVE_INTEGER,MMDone,&
!         & attribIntDims,hdferr)
!   if (hdferr /= 0) stop 'Failed to read momentum matrix status.'
!   if (MMDone == 1) then
!      write(20,*) "Two-center-ish momentum matrix already exists. Skipping."
!      call h5aclose_f(atomMMOverlap_aid,hdferr)
!      if (hdferr /= 0) stop 'Failed to close momentum matrix status.'
!   endif
!
!   ! If all integrals are done, then just return.
!   if ((ovlpDone == 1) .and. (hamDone == 1) .and. &
!         & ((dimoDone == 1) .or. (doDIMO_PSCF == 0)) .and. &
!         & ((MMDone == 1) .or. (doOPTC_PSCF == 0)) then
!      write (20,*) "All requested integrals are already done. Skipping."
!      return
!   endif
!
!   ! Allocate space for locally defined allocatable arrays
!   allocate (currentBasisFns       (maxNumAtomAlphas,maxNumStates,2))
!   allocate (currentAlphas         (maxNumAtomAlphas,2))
!   allocate (currentlmAlphaIndex   (maxNumAtomAlphas,2))
!   allocate (currentlmIndex        (maxNumStates,2))
!   allocate (alphaDist             (maxNumAtomAlphas,maxNumAtomAlphas))
!   allocate (alphaCenter           (maxNumAtomAlphas,maxNumAtomAlphas))
!
!
!   if (doINTG == 1) then
!      allocate (pairXWaveFn2OL     (16,maxNumAtomAlphas,maxNumStates))
!      allocate (pairXWaveFn2Ham    (16,maxNumAtomAlphas,maxNumStates,spin))
!      allocate (pairXWaveFn12OL    (maxNumStates,maxNumStates))
!      allocate (pairXWaveFn12Ham   (maxNumStates,maxNumStates,spin))
!
!      ! Initialize overlap and hamiltonian matrices.
!      pairXWaveFn12OL(:,:)    = 0.0_double
!      pairXWaveFn12Ham(:,:,:) = 0.0_double
!   endif
!
!   do i = 1, numAtomSites
!
!      ! Initialize the HDF5 dataset and dataspace for this atom.
!      call initAtomHDF (i,doINTG,doMOME)
!
!      ! Obtain local copies of key data from larger global data structures for
!      !   the first looped atom.
!      call initializeAtomSite(i,1,currentAtomType,currentElements,&
!            & currentNumTotalStates,currentNumCoreStates,currentNumValeStates,&
!            & currentNumAlphas,currentlmIndex,currentlmAlphaIndex,&
!            & currentPosition,currentAlphas,currentBasisFns)
!
!      ! Begin a loop over the other atoms in the system
!      do j = i, numAtomSites
!
!         ! Obtain local copies of key data from larger global data structures
!         !   for the second looped atom.
!         call initializeAtomSite(j,2,currentAtomType,currentElements,&
!            & currentNumTotalStates,currentNumCoreStates,currentNumValeStates,&
!            & currentNumAlphas,currentlmIndex,currentlmAlphaIndex,&
!            & currentPosition,currentAlphas,currentBasisFns)
!
!         ! At this point all the data that defines the nature of the two atoms
!         !   in this pair have been copied and defined.
!
!         ! Find the lattice point closest to the difference between the two
!         !   atom sites.
!         call findLatticeVector((currentPosition(:,1)-currentPosition(:,2)),&
!               & latticeVector)
!
!         ! Determine the square of the minimum seperation distance between the
!         !   two atoms.
!         atomSiteSepSqrd = sum((currentPosition(:,1) - currentPosition(:,2) - &
!               & latticeVector(:))**2)
!
!         ! Calculate the maximum distance from either atom where the overlap
!         !   is considered to be non-negligable.
!         currentNegligLimit = logBasisFnThresh * (currentAlphas(1,1) + &
!               & currentAlphas(1,2)) / (currentAlphas(1,1) * &
!               & currentAlphas(1,2))
!
!         ! Determine if there are no alpha terms for this atom pair that fall
!         !   within the current negligability limit.  Cycle if there are none.
!         if (atomSiteSepSqrd > currentNegligLimit) cycle
!
!         ! Determine the maximum radius beyond which no lattice point will be
!         !   considered to contribute to the overlap integral for this atom
!         !   pair.
!         maxLatticeRadius = atomSiteSepSqrd + currentNegligLimit + 2.0_double*&
!                          & sqrt(atomSiteSepSqrd * currentNegligLimit)
!
!         ! Check if we have to check more lattice points than were determined
!         !   to be needed by comparing the maxlatticeRadius to the maximum
!         !   radius of the lattice points defined earlier.
!         if (maxLatticeRadius > cellSizesReal(numCellsReal)) then
!            write (20,*) 'More lattice points needed for this atom overlap pair'
!            write (20,*) 'maxLatticeRadius requested=',maxLatticeRadius
!            write (20,*) 'max available=',cellSizesReal(numCellsReal)
!            stop
!         endif
!
!         ! Form an alpha overlap matrix that is later used to determine whether
!         !   a particular alpha pair needs to be considered in the overlap
!         !   matrix or if the alpha pair do not have sufficient overlap.
!         !   NOTE that it may be possible to improve this even more by asking
!         !   some simple question that will determine quickly if the alphas for
!         !   these two atoms are the same.  In that case the alphasOverlap
!         !   matrix is symmetric and can be calculated in about half the time.
!         !   That is not checked for here, but in the future it could be.
!         ! At first I thought that this would be a good idea, but it turns out
!         !   to not be efficient for large cell systems where each atom will
!         !   have an overlap with another atom only once and even then, not
!         !   every alpha pair will be used so this might actually perform worse.
!         !   It may be possible to salvage it in one of two ways.  Save this
!         !   result for the cases of later overlap integral calculations.  This
!         !   would require a lot of extra memory usage though.  OR, have a
!         !   runtime switch that will utilize this method for the case of small
!         !   systems, and the other method for the case of large systems.  That
!         !   is probably not worth it though since small systems will be rather
!         !   fast anyway.  Well, we will see.
!         ! Another possibility (that can be applied to other similar algorithms
!         !   in other subroutines) is to calculate these values based on the
!         !   element as opposed to the current method of calculating them based
!         !   on every atom pair iteration.  In this way the values would only
!         !   have to be calculated once for each element pair.
!         do k = 1,currentNumAlphas(2)
!            alphaCenter(:currentNumAlphas(1),k) = &
!                  &  currentAlphas(:currentNumAlphas(1),1) / &
!                  & (currentAlphas(:currentNumAlphas(1),1) + &
!                  &  currentAlphas(k,2))
!            alphaDist(:currentNumAlphas(1),k) = logBasisFnThresh / &
!                  &  currentAlphas(k,2) / alphaCenter(:currentNumAlphas(1),k)
!         enddo
!
!         ! Begin a loop over all the lattice points to shift the position of
!         !   atom number 2 to all the replicated cells.
!         do k = 1, numCellsReal
!
!            ! Exit the loop when we have exceeded the necessary number of
!            !   lattice points based on distance.
!            if (cellSizesReal(k) > maxLatticeRadius) exit
!
!            ! Obtain the position of atom #2 shifted by the current lattice.
!            shiftedAtomPos(:) = currentPosition(:,2) + latticeVector(:) + &
!                  & cellDimsReal(:,k)
!
!            ! Obtain the seperation vector between atom 1 and the shifted
!            !   position of atom 2.
!            shiftedAtomSiteSep = sum ((currentPosition(:,1) - &
!                  & shiftedAtomPos(:))**2)
!
!            ! Determine if this shifted atom position puts it outside of the
!            !   above determined negligability limit for this atom pair.
!            if (shiftedAtomSiteSep > currentNegligLimit) cycle
!
!
!            ! Initialize the matrices to hold the product of the integrals
!            !    times the atom2 wave functions.
!            if (doINTG == 1) then
!               pairXWaveFn2OL(:,:currentNumAlphas(1),&
!                     &:currentNumTotalStates(2)) = 0.0_double
!               pairXWaveFn2Ham(:,:currentNumAlphas(1),&
!                     &:currentNumTotalStates(2),:spin) = 0.0_double
!            endif
!            if (doMOME == 1) then
!               pairXWaveFn2MomX(:,:currentNumAlphas(1),&
!                     &:currentNumTotalStates(2)) = 0.0_double
!               pairXWaveFn2MomY(:,:currentNumAlphas(1),&
!                     &:currentNumTotalStates(2)) = 0.0_double
!               pairXWaveFn2MomZ(:,:currentNumAlphas(1),&
!                     &:currentNumTotalStates(2)) = 0.0_double
!            endif
!
!            ! Initialize a variable to track the largest atomic alpha used
!            !   from atom 1.
!            maxAlpha1Used = 0
!
!            ! Now we start loops over the alphas of each atomic type.  For each
!            !   alpha pair the overlap is tested to determine if it is
!            !   negligable or not.  The alphas for atom 1 are tested in turn
!            !   with the first alpha of atom 2.  After all the atom 1 alphas
!            !   have been tested or the test fails we test all the alphas
!            !   of atom 2 with the first alpha of atom 1 until all the alphas
!            !   of atom 2 have been used or the test fails.  Then we increment
!            !   the loop counter by 1.  The second iteration of the loop
!            !   obviously works just like the first except that the only
!            !   difference is that instead of starting with the first alphas
!            !   of atom 1 and atom 2, we always start with the second alphas
!            !   of both atoms.
!
!            ! The loop must be up to the smaller number of alphas.
!            do l = 1, min(currentNumAlphas(1),currentNumAlphas(2))
!
!               ! Set the mode to zero to say that we are incrementing through
!               !   neither the atom 1 alpha array nor the atom 2 alpha array.
!               currentMode = 0
!
!               ! Initialize the matrix index values for the two alphas.
!               alphaIndex(:) = l
!
!               ! Check if this alpha pair has any non-negligable contribution.
!               if (alphaDist(l,l) < shiftedAtomSiteSep) exit
!
!               ! Start looping over atomic alphas looking for a negligable
!               !   contribution for each alpha pair.
!               do while (.true.)
!
!                  ! Check if we are going through atom 1 alphas.
!                  if (currentMode == 1) then
!
!                     ! Go to the next atom 1 alpha
!                     alphaIndex(1) = alphaIndex(1) + 1
!
!                     ! Check if there are no alphas left to do for atom 1
!                     if (alphaIndex(1) > currentNumAlphas(1)) then
!
!                        ! Switch mode to go through atom 2 alphas on the
!                        !   current diagonal alpha of atom 1.
!                        currentMode = 2
!                        alphaIndex(1) = l
!                        cycle
!                     endif
!
!                  ! Check if we are going through atom 2 alphas.
!                  elseif (currentMode == 2) then
!
!                     ! Go to the next atom 2 alpha
!                     alphaIndex(2) = alphaIndex(2) + 1
!
!                     ! Check if there are no alphas left to do for atom 2
!                     if (alphaIndex(2) > currentNumAlphas(2)) exit
!                  endif
!
!                  ! Check if this atom alpha pair has any non negligable
!                  !   overlap.
!                  if (alphaDist(alphaIndex(1),alphaIndex(2)) < &
!                        & shiftedAtomSiteSep) then
!
!                     ! Switch mode to go through atom 2 alphas on the current
!                     !   diagonal alpha of atom 1.
!                     if (currentMode == 1) then
!                        currentMode = 2
!                        alphaIndex(1) = l
!                        cycle
!                     else
!                        exit
!                     endif
!                  endif
!
!                  if (doINTG == 1) then
!                     ! At this point a sufficient overlap has been found for
!                     !   the current alpha pair so we start looking through
!                     !   nuclear potentials.
!
!                     ! First we find the center of the overlap between the two
!                     !   alphas.
!                     overlapCenter(:) = &
!                           & alphaCenter(alphaIndex(1),alphaIndex(2)) * &
!                           & (currentPosition(:,1) - shiftedAtomPos(:)) + &
!                           & shiftedAtomPos(:)
!
!                     ! Initialize the result matrix to zero before starting the
!                     !   nuclear potential loop.
!                     potAtomOverlap(:currentlmAlphaIndex(alphaIndex(1),1), &
!                           & :currentlmAlphaIndex(alphaIndex(2),2),:spin) = &
!                           & 0.0_double
!
!
!                     ! Initialize a counter to track which potential
!                     !    coefficient we are currently calculating on.
!                     potCoeffIndex = 0
!
!                     ! Initiate a loop over each potential site for the nuclear
!                     !   potentials.
!                     do m = 1, numPotSites
!
!                        ! Skip equivalent types now.  THey will be accounted
!                        !   for later but we want to do some setup stuff only
!                        !   once for all atoms of the same type.
!                        if (potSites(m)%firstPotType == 0) cycle
!
!                        ! Initialize the parameters for this potential site
!                        !   related to the type of this site.
!                        currentPotType = potSites(m)%potTypeAssn
!                        currentNumPotAlphas = potTypes(currentPotType)%numAlphas
!
!                        ! Initiate a loop over all the potential alphas present
!                        !   for this site including the nuclear contribution.
!                        do n = 1, currentNumPotAlphas + 1
!
!                           ! Assign the potential alpha based on the value of n.
!                           if (n <= currentNumPotAlphas) then
!
!                              ! Apply the case for the atomic potential.
!
!                              ! Increment the index value that indicates which
!                              !   potential coefficient to use of all
!                              !   coefficients in the system.
!                              potCoeffIndex = potCoeffIndex + 1
!
!                              ! Store the current potential coefficient.
!                              currentPotCoeff(:spin) = &
!                                    & potCoeffs(potCoeffIndex,:spin)
!
!                              ! Store the current potential alpha. 
!                              currentPotAlpha = &
!                                    & potTypes(currentPotType)%alphas(n)
!                           else
!
!                              ! Apply the case for the nuclear potential.
!
!                              ! Get nuclear charge associated with this type.
!                              zFactor = potTypes(currentPotType)%nucCharge
!
!                              ! If the zFactor is sufficiently small we
!                              !   consider the effect of the overlap to be
!                              !   negligable and we cycle to the next one.
!                              if (zFactor < smallThresh) cycle
!
!                              ! Get the exponential alpha factor for the
!                              !   nuclear potential.
!                              currentPotAlpha = &
!                                    & potTypes(currentPotType)%nucAlpha
!
!                           endif
!
!                           ! Determine the maximum distance beyond which the
!                           !   overlap of the three gaussians is considered
!                           !   negligable.
!                           threeAlphaDist = logBasisFnThresh * (1 - &
!                                 & shiftedAtomSiteSep / &
!                                 & alphaDist(alphaIndex(1),alphaIndex(2))) * &
!                                 & (currentAlphas(alphaIndex(1),1) + &
!                                 &  currentAlphas(alphaIndex(2),2) + &
!                                 &  currentPotAlpha) / &
!                                 & (currentAlphas(alphaIndex(1),1) + &
!                                 &  currentAlphas(alphaIndex(2),2)) / &
!                                 &  currentPotAlpha
!
!                           ! Loop over the remaining potential sites that are
!                           !   equivalent.
!                           do o = 0, potTypes(currentPotType)%multiplicity-1
!
!                              ! Initialize the parameters for this potential
!                              !   site related to its position.
!                              potPosition(:) = potSites(m+o)%cartPos(:)
!
!                              ! Locate the origin for the potential lattice sum.
!                              call findLatticeVector((overlapCenter(:) - &
!                                    & potPosition(:)), latticeVector2)
!
!                              ! Find the seperation vector and distance between
!                              !   the minimum overlap center, and the origin.
!                              centerOriginVect(:) = overlapCenter(:) - &
!                                    & potPosition(:) - latticeVector2(:)
!                              centerOriginSep = sum(centerOriginVect(:)**2)
!
!                              ! Check if largest potential is negligable or not.
!                              if (centerOriginSep > threeAlphaDist) cycle
!
!                              ! At this point it must be the case that at least
!                              !   one contribution will be calculated.
!
!                              ! First, find the cut-off radius for the potential
!                              !   lattice summation by the triangle inequality.
!                              maxLatticeRadius2 = centerOriginSep + &
!                                 & threeAlphaDist + 2.0_double * &
!                                 & sqrt (centerOriginSep * threeAlphaDist)
!
!
!                              ! Begin loop over lattice points for the site pot.
!                              do p = 1, numCellsReal
!
!                                 ! Check if this lattice point extends beyond
!                                 !   the range of negligability
!                                 if (cellSizesReal(p) > maxLatticeRadius2) exit
!
!                                 ! Get overlap center shifted by lattice point.
!                                 shiftedCenterOriginSep = &
!                                       & sum((centerOriginVect(:) - &
!                                       & cellDimsReal(:,p))**2)
!
!                                 ! Check if shifted seperation between the
!                                 !   center and the origin extends past the
!                                 !   negligability limit.  If so, cycle to the
!                                 !   next cell.
!                                 if (shiftedCenterOriginSep > &
!                                       & threeAlphaDist) cycle
!
!                                 ! Get seperation shifted by the lattice point.
!                                 shiftedPotPos(:) = potPosition(:) + &
!                                       & latticeVector2(:) + cellDimsReal(:,p)
!
!                                 if (n <= currentNumPotAlphas) then
!                                    ! Calculate the opcode to do the correct
!                                    !   set of integrals for the current alpha
!                                    !   pair.
!                                    l1l2Switch = ishft(1,(powerOfTwo(&
!                                         & currentlmAlphaIndex(&
!                                         & alphaIndex(1),1))))+ ishft(16,&
!                                         & (powerOfTwo(currentlmAlphaIndex(&
!                                         & alphaIndex(2),2))))
!
!                                    call electron3CIntg (&
!                                       & currentAlphas(alphaIndex(1),1),&
!                                       & currentAlphas(alphaIndex(2),2),&
!                                       & currentPotAlpha, currentPosition(:,1),&
!                                       & shiftedAtomPos(:), shiftedPotPos(:),&
!                                       & l1l2Switch, oneAlphaSet)
!
!                                    ! Accumulate results returned for alpha set.
!                                    do q = 1, spin
!                                    potAtomOverlap(:currentlmAlphaIndex &
!                                       & (alphaIndex(1),1),:currentlmAlphaIndex&
!                                       & (alphaIndex(2),2),q) = &
!                                       & potAtomOverlap(:currentlmAlphaIndex &
!                                       & (alphaIndex(1),1),:currentlmAlphaIndex&
!                                       & (alphaIndex(2),2),q) + &
!                                       & oneAlphaSet(:currentlmAlphaIndex &
!                                       & (alphaIndex(1),1),:currentlmAlphaIndex&
!                                       & (alphaIndex(2),2)) * &
!                                       & currentPotCoeff(q)
!                                    enddo
!                                 else
!                                   ! Calculate the opcode to do the correct set
!                                   ! of integrals for the current alpha pair
!                                   l1l2Switch = ishft(1,&
!                                     &(powerOfTwo(currentlmAlphaIndex(&
!                                     &   alphaIndex(1),1)))) &
!                                     &+ ishft(16,&
!                                     &(powerOfTwo(currentlmAlphaIndex(&
!                                     &   alphaIndex(2),2))))
!                                   
!                                   call nuclear3CIntg (&
!                                      & currentAlphas(alphaIndex(1),1),&
!                                      & currentAlphas(alphaIndex(2),2),&
!                                      & currentPotAlpha,currentPosition(:,1),&
!                                      & shiftedAtomPos(:),shiftedPotPos(:),&
!                                      & l1l2Switch,oneAlphaSet)
!
!                                    ! Accumulate results returned for alpha set.
!                                    do q = 1, spin
!                                       potAtomOverlap(:currentlmAlphaIndex &
!                                          & (alphaIndex(1),1),&
!                                          & :currentlmAlphaIndex &
!                                          & (alphaIndex(2),2),q) = &
!                                          & potAtomOverlap(:currentlmAlphaIndex&
!                                          & (alphaIndex(1),1),&
!                                          & :currentlmAlphaIndex &
!                                          & (alphaIndex(2),2),q) - &
!                                          & oneAlphaSet(:currentlmAlphaIndex &
!                                          & (alphaIndex(1),1),&
!                                          & :currentlmAlphaIndex &
!                                          & (alphaIndex(2),2)) * zFactor
!                                    enddo
!                                 endif
!                              enddo ! (p numCells)
!                           enddo ! (o multiplicity)
!                        enddo ! (n numCurrentPotAlphas+1)
!                     enddo ! (m numPots (inequivalent))
!
!
!                     ! Calculate the opcode to do the correct set of integrals
!                     ! for the current alpha pair
!                     l1l2Switch = ishft(1,&
!                       &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
!                       &+ ishft(16,&
!                       &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))
!
!                     ! Determine the kinetic energy contribution.
!                     call kinetic2CIntg (&
!                           & currentAlphas(alphaIndex(1),1),&
!                           & currentAlphas(alphaIndex(2),2),&
!                           & currentPosition(:,1), shiftedAtomPos(:),&
!                           & l1l2Switch, oneAlphaSet)
!
!                     ! Accumulate the contribution from this alpha pair
!                     do m = 1, spin
!                        potAtomOverlap(:currentlmAlphaIndex &
!                              & (alphaIndex(1),1),:currentlmAlphaIndex &
!                              & (alphaIndex(2),2),m) = &
!                              & potAtomOverlap(:currentlmAlphaIndex &
!                              & (alphaIndex(1),1),:currentlmAlphaIndex &
!                              & (alphaIndex(2),2),m) + &
!                              & oneAlphaSet(:currentlmAlphaIndex &
!                              & (alphaIndex(1),1),:currentlmAlphaIndex &
!                              & (alphaIndex(2),2))
!                     enddo
!
!                     ! Compute the mass velocity integral if needed for the
!                     !   scalar relativistic calculation.
!                     if (rel == 1) then
!                        ! Calculate the opcode to do the correct set of
!                        !   integrals for the current alpha pair.
!                        l1l2Switch = ishft(1,&
!                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
!                        &+ ishft(16,&
!                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))
!
!                        ! Compute the integral.
!                        call massVel2CIntg (currentAlphas(alphaIndex(1),1),&
!                              & currentAlphas(alphaIndex(2),2),&
!                              & currentPosition(:,1), shiftedAtomPos(:),&
!                              & l1l2Switch, oneAlphaSet)
!                           
!                        ! Accumulate the contribution from this alpha pair.
!                        !   Note: a minus sign is used in the accumulation
!                        !   for the mass velocity integral because it has an
!                        !   overall negative contribution to the Hamiltonian.
!                        do m = 1, spin
!                           potAtomOverlap(:currentlmAlphaIndex &
!                                 & (alphaIndex(1),1),:currentlmAlphaIndex &
!                                 & (alphaIndex(2),2),m) = &
!                                 & potAtomOverlap(:currentlmAlphaIndex &
!                                 & (alphaIndex(1),1),:currentlmAlphaIndex &
!                                 & (alphaIndex(2),2),m) - &
!                                 & oneAlphaSet(:currentlmAlphaIndex &
!                                 & (alphaIndex(1),1),:currentlmAlphaIndex &
!                                 & (alphaIndex(2),2))
!                        enddo
!                     endif
!
!                     ! FIX: It seems unnecessary to recompute the l1l2Switch
!                     !   for OL and MV matrices after it is computed
!                     !   once for the KE.
!                     ! Calculate the opcode to do the correct set of integrals
!                     ! for the current alpha pair.
!                     l1l2Switch = ishft(1,&
!                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
!                        &+ ishft(16,&
!                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))
!                
!                     !print*,l1l2Switch
!                     ! We can proceed with the next step of the calculation.
!                     ! This is the actual integral.
!                     call overlap2CIntg (&
!                           & currentAlphas(alphaIndex(1),1),&
!                           & currentAlphas(alphaIndex(2),2), &
!                           & currentPosition(:,1), shiftedAtomPos(:),&
!                           & l1l2Switch, oneAlphaSet)
!                  endif
!
!
!                  ! Compute the momentum matrix values if requested.
!                  if (doMOME == 1) then
!                     ! Calculate the opcode to do the correct set
!                     ! of integrals for the current alpha pair
!                     l1l2Switch = ishft(1,&
!                        & (powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
!                        & + ishft(16,&
!                        & (powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))
!
!                     call momentum2CIntg (&
!                           & currentAlphas(alphaIndex(1),1),&
!                           & currentAlphas(alphaIndex(2),2),&
!                           & currentPosition(:,1), shiftedAtomPos(:),&
!                           & l1l2Switch, oneAlphaSetMom)
!                  endif
!
!                  ! Collect the results of the overlap of the current alpha
!                  !   times the wave functions for atom 2.
!
!                  if (doINTG == 1) then
!                     ! Potential overlaps first (if any were found).
!                     do q = 1, spin
!                        do m = 1, currentNumTotalStates(2)
!                           pairXWaveFn2Ham(:currentlmAlphaIndex( &
!                                 & alphaIndex(1),1),alphaIndex(1),m,q) = &
!                                 & pairXWaveFn2Ham(:currentlmAlphaIndex &
!                                 & (alphaIndex(1),1),alphaIndex(1),m,q) + &
!                                 & potAtomOverlap(:currentlmAlphaIndex &
!                                 & (alphaIndex(1),1),currentlmIndex(m,2),q) * &
!                                 & currentBasisFns(alphaIndex(2),m,2)
!                        enddo
!                     enddo
!                     ! Atom pair overlaps second.
!                     do m = 1, currentNumTotalStates(2)
!                        pairXWaveFn2OL(:currentlmAlphaIndex(alphaIndex(1),1),&
!                              & alphaIndex(1),m) = &
!                              & pairXWaveFn2OL(:currentlmAlphaIndex &
!                              & (alphaIndex(1),1),alphaIndex(1),m) + &
!                              & oneAlphaSet(:currentlmAlphaIndex &
!                              & (alphaIndex(1),1),currentlmIndex(m,2)) * &
!                              & currentBasisFns(alphaIndex(2),m,2)
!                     enddo
!                  endif
!
!                  ! Momentum matrix last.  (If requested)
!                  if (doMOME == 1) then
!                     do m = 1, currentNumTotalStates(2)
!                        pairXWaveFn2MomX(:currentlmAlphaIndex(alphaIndex(1),&
!                              & 1), alphaIndex(1),m) = &
!                              & pairXWaveFn2MomX(:currentlmAlphaIndex &
!                              & (alphaIndex(1),1),alphaIndex(1),m) + &
!                              & oneAlphaSetMom(:currentlmAlphaIndex &
!                              & (alphaIndex(1),1),currentlmIndex(m,2),1) * &
!                              & currentBasisFns(alphaIndex(2),m,2)
!                        pairXWaveFn2MomY(:currentlmAlphaIndex(alphaIndex(1),&
!                              & 1), alphaIndex(1),m) = &
!                              & pairXWaveFn2MomY(:currentlmAlphaIndex &
!                              & (alphaIndex(1),1),alphaIndex(1),m) + &
!                              & oneAlphaSetMom(:currentlmAlphaIndex &
!                              & (alphaIndex(1),1),currentlmIndex(m,2),2) * &
!                              & currentBasisFns(alphaIndex(2),m,2)
!                        pairXWaveFn2MomZ(:currentlmAlphaIndex(alphaIndex(1),&
!                              & 1), alphaIndex(1),m) = &
!                              & pairXWaveFn2MomZ(:currentlmAlphaIndex &
!                              & (alphaIndex(1),1),alphaIndex(1),m) + &
!                              & oneAlphaSetMom(:currentlmAlphaIndex &
!                              & (alphaIndex(1),1),currentlmIndex(m,2),3) * &
!                              & currentBasisFns(alphaIndex(2),m,2)
!                     enddo
!                  endif
!
!                  ! Update the maximum alpha used from atom 1.
!                  maxAlpha1Used = max(alphaIndex(1),maxAlpha1Used)
!
!                  ! Switch mode from the initial state of no searching to the
!                  !   state of searching along alphas from atom 1.
!                  if (currentMode == 0) then
!                     currentMode = 1
!                  endif
!               enddo
!            enddo  ! min number of alphas between the two atoms l
!
!            if (maxAlpha1Used > 0) then
!
!               ! At this point all the alpha loops are complete and we can form
!               !   a product with the atom 1 basis functions to give the
!               !   overlap and (spin) hamiltonian integral matrix elements.
!               if (doINTG == 1) then
!                  call multWithBasisFn1 (currentBasisFns,pairXBasisFn2OL,&
!                        & pairXBasisFn12,currentlmIndex,currentNumTotalStates,&
!                        & maxAlpha1Used)
!                  do q = 1, spin
!                     call multWithBasisFn1 (currentBasisFns,&
!                           & pairXBasisFn2Ham(:,:,:,q),&
!                           & pairXBasisFn12Ham(:,:,q),currentlmIndex,&
!                           & currentNumTotalStates,maxAlpha1Used)
!                  enddo
!
!                  ! Collect this atom 1, atom 2 basis function overlap matrix
!                  !   for all bloch vectors (kpoints) with phase factors
!                  !   appropriate to the current atom 2 lattice vector. NOTE
!                  !   that the k index is over the number of cells in the
!                  !   superlattice.
!#ifndef GAMMA
!                  call applyPhaseFactors (currentPair,pairXBasisFn12(1:&
!                        & currentNumTotalStates(1),1:&
!                        & currentNumTotalStates(2)),currentNumTotalStates(1),&
!                        & currentNumTotalStates(2),k,0,0)
!                  do q = 1, spin
!                     call applyPhaseFactors (currentPair,pairXBasisFn12Ham(1:&
!                           & currentNumTotalStates(1),1:&
!                           & currentNumTotalStates(2),q),&
!                           & currentNumTotalStates(1),&
!                           & currentNumTotalStates(2),k,0,0)
!                  enddo
!#else
!                  call applyPhaseFactorsGamma (currentPairGamma,&
!                        & pairXBasisFn12(1:currentNumTotalStates(1),&
!                        & 1:currentNumTotalStates(2)),&
!                        & currentNumTotalStates(1),currentNumTotalStates(2),0)
!                  do q = 1, spin
!                     call applyPhaseFactorsGamma (currentPairGamma,&
!                           & pairXBasisFn12Ham(1:&
!                           & currentNumTotalStates(1),1:&
!                           & currentNumTotalStates(2),q),&
!                           & currentNumTotalStates(1),&
!                           & currentNumTotalStates(2),0)
!                  enddo
!#endif
!               endif ! doINTG
!
!               ! The momentum is summed against the wave function 1 if needed.
!               if (doOPTC_PSCF >= 1) then
!                  do q = 1, 3
!                     call multWithBasisFn1 (currentBasisFns,&
!                           & pairXBasisFn2Mom(:,:,q),&
!                           & pairXBasisFn12Mom(:,:,q),&
!                           & currentlmIndex,currentNumTotalStates,&
!                           & maxAlpha1Used)
!                  enddo
!#ifndef GAMMA
!                  do q = 1,3
!                     call applyPhaseFactors (currentPair,pairXBasisFn12Mom(1:&
!                           & currentNumTotalStates(1),1:&
!                           & currentNumTotalStates(2),q),&
!                           & currentNumTotalStates(1),&
!                           & currentNumTotalStates(2),0)
!                  enddo
!#else
!                  do q = 1,3
!                     call applyPhaseFactorsGamma (currentPairGamma,&
!                           & pairXBasisFn12Mom(1:&
!                           & currentNumTotalStates(1),1:&
!                           & currentNumTotalStates(2),q),&
!                           & currentNumTotalStates(1),&
!                           & currentNumTotalStates(2),0)
!                  enddo
!#endif
!               endif ! doOPTC_PSCF
!            endif ! maxAlpha1Used > 0
!
!         enddo !(k superlattice)
!
!      enddo ! (Atom loop #2)
!
!
!      ! Mark the completion of this atom.
!      if (mod(i,10) .eq. 0) then
!         write (20,ADVANCE="NO",FMT="(a1)") "|"
!      else
!         write (20,ADVANCE="NO",FMT="(a1)") "."
!      endif
!      if (mod(i,50) .eq. 0) then
!         write (20,*) " ",i
!      endif
!      call flush (20)
!
!   enddo    ! (Atom loop #1)
!
!   ! Deallocate all arrays and matrices before exiting.
!   deallocate (currentBasisFns)
!   deallocate (currentAlphas)
!   deallocate (currentlmAlphaIndex)
!   deallocate (currentlmIndex)
!   deallocate (alphaDist)   ! Can be saved from before
!   deallocate (alphaCenter) ! Can be saved for later.
!
!   if (doINTG == 1) then
!      deallocate (pairXWaveFn2OL)
!      deallocate (pairXWaveFn2Ham)
!      deallocate (pairXWaveFn12OL)
!      deallocate (pairXWaveFn12Ham)
!   endif
!
!   if (doMOME == 1) then
!      deallocate (pairXWaveFn2MomX)
!      deallocate (pairXWaveFn2MomY)
!      deallocate (pairXWaveFn2MomZ)
!      deallocate (pairXWaveFn12MomX)
!      deallocate (pairXWaveFn12MomY)
!      deallocate (pairXWaveFn12MomZ)
!   endif
!
!   ! Log the date and time we start.
!   call timeStampEnd(20)
!
!
!end subroutine allIntgCombo


subroutine orthoOL

   ! Use necessary modules.
   use HDF5
   use O_Kinds
   use O_AtomicSites, only: coreDim, valeDim
   use O_KPoints, only: numKPoints
   use O_SCFIntegralsHDF5, only: atomOverlap_did, atomDims, atomOverlap_aid
   use O_Orthogonalization

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define local variables.
   integer :: i,j,k
   integer :: hdferr
   integer :: currIndex
   real (kind=double), allocatable, dimension (:,:) :: packedValeVale
   integer(hsize_t), dimension (1) :: attribIntDims ! Attribute dataspace dim

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
   enddo ! numKPoints

   ! Record that the overlap calculation is complete and close the attribute.
   attribIntDims(1) = 1
   call h5awrite_f(atomOverlap_aid,H5T_NATIVE_INTEGER,1,attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to record atom overlap success.'
   call h5aclose_f(atomOverlap_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the atom overlap attribute.'

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
   use HDF5
   use O_Kinds
   use O_PotTypes, only: potTypes
   use O_KPoints, only: numKPoints
   use O_AtomicSites, only: coreDim, valeDim
   use O_SCFIntegralsHDF5, only: atomKEOverlap_did, atomMVOverlap_did, &
         & atomNucOverlap_did, atomPotOverlap_did, atomKEOverlap_aid, &
         & atomNucOverlap_aid, atomMVOverlap_aid, atomPotTermOL_aid,&
         & atomDims
   use O_Orthogonalization

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed dummy arguments.
   integer :: opCode

   ! Define local variables.
   integer :: i,j,k
   integer :: hdferr
   integer :: currIndex
   integer :: potTermIdx
   real (kind=double), allocatable, dimension (:,:) :: packedValeVale
   integer(hsize_t), dimension (1) :: attribIntDims ! Attribute dataspace dim

   ! Identify the index number of the current potential term in the list of
   !   all potential terms (of all types). This only needs to be done for
   !   electron potential integrals.
   if (opCode == 4) then
      potTermIdx = potTypes(currPotTypeNumber)%cumulAlphaSum + currAlphaNumber
   endif

   ! Orthogonalizing against the overlap matrix is unnecessary when the core
   !   dimension is zero.  However, we must still allocate some matrices for
   !   use later.
   if (coreDim /= 0) then
#ifndef GAMMA
      do i = 1, numKPoints
         ! Form product of (valeCoreOL)(coreVale) and (valeCore)(coreValeOL).
         !   Subtract both from the target matrix elements (valeVale).
         call valeCoreCoreVale (valeDim,coreDim,valeVale(:,:,i,1),&
               & coreVale(:,:,i),coreValeOL(:,:,i))
      enddo
#else
      do i = 1, numKPoints
         ! Form a product of (valeCoreOL)(coreVale) and (valeCore)
         !   (coreValeOL).  Subtract both from the target matrix elements
         !   (valeValeGamma).
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
         call h5dwrite_f(atomPotOverlap_did(i,potTermIdx),&
               & H5T_NATIVE_DOUBLE,packedValeVale(:,:),atomDims,hdferr)
         if (hdferr /= 0) stop 'Failed to write electronic potential vale vale'
      case (5) ! Implies that we are doing a scalar rel. calculation.
         call h5dwrite_f(atomMVOverlap_did(i),H5T_NATIVE_DOUBLE,&
               & packedValeVale(:,:),atomDims,hdferr)
         if (hdferr /= 0) stop 'Failed to write mass velocity vale vale'
      end select
   enddo ! numKPoints i

   ! Record that the calculation is complete and close the attribute.
   attribIntDims(1) = 1
   select case (opCode)
   case (2)
      call h5awrite_f(atomKEOverlap_aid,H5T_NATIVE_INTEGER,1,&
            & attribIntDims,hdferr)
      if (hdferr /= 0) stop 'Failed to record atom KE overlap success.'
      call h5aclose_f(atomKEOverlap_aid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atom KE overlap attribute.'
   case (3)
      call h5awrite_f(atomNucOverlap_aid,H5T_NATIVE_INTEGER,1,&
            & attribIntDims,hdferr)
      if (hdferr /= 0) stop 'Failed to record atom nuclear overlap success.'
      call h5aclose_f(atomNucOverlap_aid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atom nuclear overlap attribute.'
   case (4)
      call h5awrite_f(atomPotTermOL_aid(potTermIdx),H5T_NATIVE_INTEGER,1,&
            & attribIntDims,hdferr)
      if (hdferr /= 0) stop 'Failed to record atom pot term overlap success.'
      call h5aclose_f(atomPotTermOL_aid(potTermIdx),hdferr)
      if (hdferr /= 0) stop 'Failed to close atom pot term overlap attribute.'
   case (5)
      call h5awrite_f(atomMVOverlap_aid,H5T_NATIVE_INTEGER,1,&
            & attribIntDims,hdferr)
      if (hdferr /= 0) stop 'Failed to record atom MV overlap success.'
      call h5aclose_f(atomMVOverlap_aid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atom MV overlap attribute.'
   end select

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
   !deallocate (coreValeOL) ! Needed in 3Terms
#else
   deallocate (coreCoreGamma)
   deallocate (valeValeGamma)
   !deallocate (coreValeOLGamma) ! Needed in 3Terms
#endif

end subroutine cleanUpIntegralsSCF


subroutine secondCleanUpIntegralsSCF

   implicit none

   ! If all integrals were skipped, then coreValeOL was never allocated.
   if (someIntegralDone == 0) then
      return
   endif

#ifndef GAMMA
   deallocate (coreValeOL) ! Needed in 3Terms
#else
   deallocate (coreValeOLGamma) ! Needed in 3Terms
#endif

end subroutine secondCleanUpIntegralsSCF


end module O_Integrals
