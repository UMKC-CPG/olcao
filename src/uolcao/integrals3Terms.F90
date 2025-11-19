module O_Integrals3Terms

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
   !   valence sections and their coreVale and valeCore interactions. The
   !   overlap valeVale and coreVale will be pulled from the previous
   !   integrals calculation for the purpose of orthogonalization.

#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:,:) :: coreCore
   complex (kind=double), allocatable, dimension (:,:,:)   :: valeCore
   complex (kind=double), allocatable, dimension (:,:,:,:) :: coreVale
   complex (kind=double), allocatable, dimension (:,:,:,:) :: valeVale
#else
   real    (kind=double), allocatable, dimension (:,:,:)   :: coreCoreGamma
   real    (kind=double), allocatable, dimension (:,:,:)   :: valeCoreGamma
   real    (kind=double), allocatable, dimension (:,:,:)   :: coreValeGamma
   real    (kind=double), allocatable, dimension (:,:,:)   :: valeValeGamma
#endif


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


subroutine allocateIntegrals3Terms(coreDim,valeDim,numKPoints)

   implicit none

   ! Define passed dummy parameters.
   integer, intent(in) :: coreDim
   integer, intent(in) :: valeDim
   integer, intent(in) :: numKPoints

#ifndef GAMMA
   allocate (coreCore (coreDim,coreDim,numKPoints,3))
   allocate (valeVale (valeDim,valeDim,numKPoints,3))
#else
   allocate (coreCoreGamma (coreDim,coreDim,3))
   allocate (valeValeGamma (valeDim,valeDim,3))
#endif

end subroutine allocateIntegrals3Terms



! Three term dipole moment integral.
subroutine gaussOverlapDM(packedVVDims,did,aid)

   ! Import necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Constants, only: dim3
   use O_Input, only: dipoleCenter
   use O_KPoints, only: numKPoints
   use O_GaussianRelations, only: alphaDist
   use O_AtomicSites, only: valeDim, coreDim, numAtomSites, atomSites
   use O_AtomicTypes, only: maxNumAtomAlphas, maxNumStates, atomTypes
   use O_Lattice, only: numCellsReal, cellSizesReal, cellDimsReal, &
         & findLatticeVector
   use O_GaussianIntegrals, only: dipole3CIntg
   use O_Basis, only: initializeAtomSite
   use O_IntgSaving

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define passed parameters.
   integer(hsize_t), dimension(2), intent(in) :: packedVVDims
   integer(hid_t), dimension(numKPoints,3), intent(in) :: did
   integer(hid_t), intent(in) :: aid

   ! Define local variables for logging and loop control
   integer :: i,j,k,l,m,n ! Loop index variables
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
   real (kind=double), dimension (16,16,3) :: oneAlphaSetDM ! One pair in xyz.
   real (kind=double), allocatable, dimension (:,:,:,:) :: pairXBasisFn2 ! The
         ! above overlap times the basis function from atom 2.  This is
         ! accumulated with each iteration of the alpha loop.

   logical :: contrib  ! At least one alpha pair contributes.  For the nuc and
         ! elec calculations this also requires that the pot term contribute.


   ! Variables and data structures that change or are accumulated with each
   !   iteration of the lattice loop.
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:,:) :: currentPair
#else
   real (kind=double), allocatable, dimension (:,:,:) :: currentPairGamma
#endif
   real (kind=double), allocatable, dimension (:,:,:) :: pairXBasisFn12


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

   ! Record the beginning of this phase of the calculation.
   call timeStampStart (29)

   ! Determine if this calculation has already been completed by a previous
   !   OLCAO execution.
   hdf5Status = 0
   attribIntDims(1) = 1
   call h5aread_f(aid,H5T_NATIVE_INTEGER,hdf5Status,attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to read atom DM overlap status.'
   if (hdf5Status == 1) then
      write(20,*) "Three-center DM overlap already exists. Skipping."
      call timeStampEnd(29)
      call h5aclose_f(aid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atom DM overlap status.'
      return
   endif

   ! Allocate space for locally defined allocatable arrays
   allocate (currentBasisFns     (maxNumAtomAlphas,maxNumStates,2))
   allocate (currentAlphas       (maxNumAtomAlphas,2))
   allocate (currentlmAlphaIndex (maxNumAtomAlphas,2))
   allocate (currentlmIndex      (maxNumStates,2))
   allocate (pairXBasisFn2       (16,maxNumAtomAlphas,maxNumStates,3))
   allocate (pairXBasisFn12      (maxNumStates,maxNumStates,3))

#ifndef GAMMA
   allocate (coreVale            (coreDim,valeDim,numKPoints,3))
   allocate (currentPair         (maxNumStates,maxNumStates,numKPoints,3))
#else
   allocate (coreValeGamma       (coreDim,valeDim,3))
   allocate (currentPairGamma    (maxNumStates,maxNumStates,3))
#endif


   ! Initialize key matrices
#ifndef GAMMA
   coreVale      (:,:,:,:) = 0.0_double
   valeVale      (:,:,:,:) = 0.0_double
   coreCore      (:,:,:,:) = 0.0_double
#else
   coreValeGamma (:,:,:)   = 0.0_double
   valeValeGamma (:,:,:)   = 0.0_double
   coreCoreGamma (:,:,:)   = 0.0_double
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
         currentPair(:,:,:,:) = 0.0_double
#else
         currentPairGamma(:,:,:) = 0.0_double
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
            ! times the atom2 basis functions in the XYZ direction.
            pairXBasisFn2(:,:currentNumAlphas(1),&
                  &:currentNumTotalStates(2),:) = 0.0_double

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

                  ! At this point a sufficient overlap between atomic
                  !   Gaussians has been shown to exist. FIX: This appears
                  !   to be a useless term. Similarly, the .eqv. test below
                  !   becomes useless.
                  contrib = .true.

                  ! Calculate the opcode to do the correct set of integrals
                  !   for the current alpha pair
                  l1l2Switch = ishft(1,&
                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
                        &+ ishft(16,&
                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))

                  ! We can proceed with the next step of the calculation.
                  !   This is the actual integral from the Obara-Saika papers.
                  call dipole3CIntg (currentAlphas(alphaIndex(1),1), &
                        & currentAlphas(alphaIndex(2),2), &
                        & currentPosition(:,1), shiftedAtomPos(:), &
                        & dipoleCenter(:), l1l2Switch, oneAlphaSetDM)


                  ! Collect the results of the overlap of the current alpha
                  !   times the basis functions of atom 2.
                  if (contrib .eqv. .true.) then
                     do n = 1, 3
                        do m = 1, currentNumTotalStates(2)
                           pairXBasisFn2(:currentlmAlphaIndex(&
                              & alphaIndex(1),1),alphaIndex(1),m,n) = &
                              & pairXBasisFn2(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),alphaIndex(1),m,n) + &
                              & oneAlphaSetDM(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),currentlmIndex(m,2),n) * &
                              & currentBasisFns(alphaIndex(2),m,2)
                        enddo
                     enddo

                     ! Update the maximum alpha used from atom 1.
                     maxAlpha1Used = max(alphaIndex(1),maxAlpha1Used)
                  endif
               enddo
            enddo  ! min number of alphas between the two atoms l


            ! At this point all the alpha loops are complete and we can form
            !   a product with the atom 1 basis function to give the overlap
            !   integral in a complete basis representation.
            do l = 1, 3
               call multWithBasisFn1 (currentBasisFns,pairXBasisFn2(:,:,:,l),&
                     & pairXBasisFn12(:,:,l),currentlmIndex,&
                     & currentNumTotalStates,maxAlpha1Used)
            enddo


            ! Collect this atom 1, atom 2 basis function overlap matrix for
            !   all bloch vectors (kpoints) with phase factors appropriate
            !   to the current atom 2 lattice vector.  NOTE that the k
            !   index is over the number of cells in the superlattice.
            do l = 1, 3
#ifndef GAMMA
               call applyPhaseFactors (currentPair(:,:,:,l),&
                     & pairXBasisFn12(1:currentNumTotalStates(1),&
                     & 1:currentNumTotalStates(2),l),&
                     & currentNumTotalStates(1),currentNumTotalStates(2),&
                     & k,0)
#else
               call applyPhaseFactorsGamma (currentPairGamma(:,:,l),&
                     & pairXBasisFn12(1:currentNumTotalStates(1),&
                     & 1:currentNumTotalStates(2),l),&
                     & currentNumTotalStates(1),currentNumTotalStates(2),0)
#endif
            enddo
         enddo !(k superlattice)

         ! At this point all the lattice sums for the current atom pair are
         !   complete.

         ! So now we can arrange the data from this atom into a set of three
         !   large matrices.  A valence-valence matrix, a core-valence matrix,
         !   and a core-core matrix.

#ifndef GAMMA
         ! First we must make a correction for the atom 2 lattice origin shift.
         do k = 1, 3
            call kPointLatticeOriginShift (currentNumTotalStates,&
                  & currentPair(:,:,:,k),latticeVector)
            call saveCurrentPair(i,j,numKPoints,currentPair(:,:,:,k),&
                  & valeVale(:,:,:,k),coreVale(:,:,:,k),coreCore(:,:,:,k),0)
         enddo
#else
         do k = 1, 3
            call saveCurrentPairGamma(i,j,currentPairGamma(:,:,k),&
                  & valeValeGamma(:,:,k),coreValeGamma(:,:,k),&
                  & coreCoreGamma(:,:,k))
         enddo
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
   call ortho(6,packedVVDims,did,aid)

   ! Make a finishing time stamp.
   call timeStampEnd (29)

end subroutine gaussOverlapDM


! Three term (xyz) momentum matrix integral.
subroutine gaussOverlapMM(packedVVDims,did,aid)

   ! Import necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Constants, only: dim3
   use O_Input, only: dipoleCenter
   use O_KPoints, only: numKPoints
   use O_GaussianRelations, only: alphaDist
   use O_AtomicSites, only: valeDim, coreDim, numAtomSites, atomSites
   use O_AtomicTypes, only: maxNumAtomAlphas, maxNumStates, atomTypes
   use O_Lattice, only: numCellsReal, cellSizesReal, cellDimsReal, &
         & findLatticeVector, logBasisFnThresh
   use O_GaussianIntegrals, only: momentum2CIntg
   use O_Basis, only: initializeAtomSite
   use O_IntgSaving

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define passed parameters.
   integer(hsize_t), dimension(2), intent(in) :: packedVVDims
   integer(hid_t), dimension(numKPoints,3), intent(in) :: did
   integer(hid_t), intent(in) :: aid

   ! Define local variables for logging and loop control
   integer :: i,j,k,l,m,n ! Loop index variables
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
   real (kind=double), dimension (16,16,3) :: oneAlphaSetMM ! One pair in xyz.
   real (kind=double), allocatable, dimension (:,:,:,:) :: pairXBasisFn2 ! The
         ! above overlap times the basis function from atom 2.  This is
         ! accumulated with each iteration of the alpha loop.

   logical :: contrib  ! At least one alpha pair contributes.  For the nuc and
         ! elec calculations this also requires that the pot term contribute.


   ! Variables and data structures that change or are accumulated with each
   !   iteration of the lattice loop.
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:,:) :: currentPair
#else
   real (kind=double), allocatable, dimension (:,:,:) :: currentPairGamma
#endif
   real (kind=double), allocatable, dimension (:,:,:) :: pairXBasisFn12


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

   ! Record the beginning of this phase of the calculation.
   call timeStampStart (12)

   ! Determine if this calculation has already been completed by a previous
   !   OLCAO execution.
   hdf5Status = 0
   attribIntDims(1) = 1
   call h5aread_f(aid,H5T_NATIVE_INTEGER,hdf5Status,attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to read atom MM overlap status.'
   if (hdf5Status == 1) then
      write(20,*) "Two-center MM overlap already exists. Skipping."
      call timeStampEnd(12)
      call h5aclose_f(aid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atom MM overlap status.'
      return
   endif

   ! Allocate space for locally defined allocatable arrays
   allocate (currentBasisFns     (maxNumAtomAlphas,maxNumStates,2))
   allocate (currentAlphas       (maxNumAtomAlphas,2))
   allocate (currentlmAlphaIndex (maxNumAtomAlphas,2))
   allocate (currentlmIndex      (maxNumStates,2))
   allocate (pairXBasisFn2       (16,maxNumAtomAlphas,maxNumStates,3))
   allocate (pairXBasisFn12      (maxNumStates,maxNumStates,3))

#ifndef GAMMA
   allocate (coreVale            (coreDim,valeDim,numKPoints,3))
   allocate (currentPair         (maxNumStates,maxNumStates,numKPoints,3))
#else
   allocate (coreValeGamma       (coreDim,valeDim,3))
   allocate (currentPairGamma    (maxNumStates,maxNumStates,3))
#endif


   ! Initialize key matrices
#ifndef GAMMA
   coreVale      (:,:,:,:) = 0.0_double
   valeVale      (:,:,:,:) = 0.0_double
   coreCore      (:,:,:,:) = 0.0_double
#else
   coreValeGamma (:,:,:)   = 0.0_double
   valeValeGamma (:,:,:)   = 0.0_double
   coreCoreGamma (:,:,:)   = 0.0_double
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
!         currentNegligLimit = alphaDist(1,1,currentElements(1),&
!               & currentElements(2))
         currentNegligLimit = logBasisFnThresh * (currentAlphas(1,1) + &
               & currentAlphas(1,2)) / (currentAlphas(1,1) * &
               & currentAlphas(1,2))

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
         currentPair(:,:,:,:) = 0.0_double
#else
         currentPairGamma(:,:,:) = 0.0_double
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
            ! times the atom2 basis functions in the XYZ direction.
            pairXBasisFn2(:,:currentNumAlphas(1),&
                  & :currentNumTotalStates(2),:) = 0.0_double

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

                  ! At this point a sufficient overlap between atomic
                  !   Gaussians has been shown to exist. FIX: This appears
                  !   to be a useless term. Similarly, the .eqv. test below
                  !   becomes useless.
                  contrib = .true.

                  ! Calculate the opcode to do the correct set of integrals
                  !   for the current alpha pair
                  l1l2Switch = ishft(1,&
                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
                        &+ ishft(16,&
                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))

                  ! We can proceed with the next step of the calculation.
                  !   This is the actual integral from the Obara-Saika papers.

                  call momentum2CIntg (currentAlphas(alphaIndex(1),1), &
                        & currentAlphas(alphaIndex(2),2), &
                        & currentPosition(:,1), shiftedAtomPos(:), &
                        & l1l2Switch, oneAlphaSetMM)

                  ! Collect the results of the overlap of the current alpha
                  !   times the basis functions of atom 2.
                  if (contrib .eqv. .true.) then
                     do n = 1, 3
                        do m = 1, currentNumTotalStates(2)
                           pairXBasisFn2(:currentlmAlphaIndex(&
                              & alphaIndex(1),1),alphaIndex(1),m,n) = &
                              & pairXBasisFn2(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),alphaIndex(1),m,n) + &
                              & oneAlphaSetMM(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),currentlmIndex(m,2),n) * &
                              & currentBasisFns(alphaIndex(2),m,2)
                        enddo
                     enddo

                     ! Update the maximum alpha used from atom 1.
                     maxAlpha1Used = max(alphaIndex(1),maxAlpha1Used)
                  endif
               enddo
            enddo  ! min number of alphas between the two atoms l


            ! At this point all the alpha loops are complete and we can form
            !   a product with the atom 1 basis function to give the overlap
            !   integral in a complete basis representation.
            do l = 1, 3
               call multWithBasisFn1 (currentBasisFns,pairXBasisFn2(:,:,:,l),&
                     & pairXBasisFn12(:,:,l),currentlmIndex,&
                     & currentNumTotalStates,maxAlpha1Used)
            enddo

            ! Collect this atom 1, atom 2 basis function overlap matrix for
            !   all bloch vectors (kpoints) with phase factors appropriate
            !   to the current atom 2 lattice vector.  NOTE that the k
            !   index is over the number of cells in the superlattice.
            do l = 1, 3
#ifndef GAMMA
               call applyPhaseFactors (currentPair(:,:,:,l),&
                     & pairXBasisFn12(1:currentNumTotalStates(1),&
                     & 1:currentNumTotalStates(2),l),&
                     & currentNumTotalStates(1),currentNumTotalStates(2),&
                     & k,1)
#else
               call applyPhaseFactorsGamma (currentPairGamma(:,:,l),&
                     & pairXBasisFn12(1:currentNumTotalStates(1),&
                     & 1:currentNumTotalStates(2),l),&
                     & currentNumTotalStates(1),currentNumTotalStates(2),1)
#endif
            enddo
         enddo !(k superlattice)

         ! At this point all the lattice sums for the current atom pair are
         !   complete.

         ! So now we can arrange the data from this atom into a set of three
         !   large matrices.  A valence-valence matrix, a core-valence matrix,
         !   and a core-core matrix.

         do k = 1, 3
#ifndef GAMMA
            ! Make a correction for the atom 2 lattice origin shift.
            call kPointLatticeOriginShift (currentNumTotalStates,&
                  & currentPair(:,:,:,k),latticeVector)
            call saveCurrentPair(i,j,numKPoints,currentPair(:,:,:,k),&
                  & valeVale(:,:,:,k),coreVale(:,:,:,k),coreCore(:,:,:,k),0)
#else
            call saveCurrentPairGamma(i,j,currentPairGamma(:,:,k),&
                  & valeValeGamma(:,:,k),coreValeGamma(:,:,k),&
                  & coreCoreGamma(:,:,k))
#endif
         enddo
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
   call ortho(7,packedVVDims,did,aid)

   ! Make a finishing time stamp.
   call timeStampEnd (12)

end subroutine gaussOverlapMM


! Two center K-Overlap integral.
#ifndef GAMMA
subroutine gaussKOverlap(packedVVDims,did,aid)

   ! Import necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Constants, only: dim3
   use O_Input, only: dipoleCenter
   use O_KPoints, only: numKPoints, numAxialKPoints
   use O_GaussianRelations, only: alphaDist
   use O_AtomicSites, only: valeDim, coreDim, numAtomSites, atomSites
   use O_AtomicTypes, only: maxNumAtomAlphas, maxNumStates, atomTypes
   use O_Lattice, only: numCellsReal, cellSizesReal, cellDimsReal, &
         & findLatticeVector, logBasisFnThresh, recipVectors
   use O_GaussianIntegrals, only: KOverlap2CIntg
   use O_Basis, only: initializeAtomSite
   use O_IntgSaving

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define passed parameters.
   integer(hsize_t), dimension(2), intent(in) :: packedVVDims
   integer(hid_t), dimension(numKPoints,3), intent(in) :: did
   integer(hid_t), intent(in) :: aid

   ! Define local variables for logging and loop control
   integer :: axis
   integer :: i,j,k,l,m,n ! Loop index variables
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
   complex (kind=double), dimension (16,16,3) :: oneAlphaSetK ! One pair in xyz.
   complex (kind=double), allocatable, dimension (:,:,:,:) :: pairXBasisFn2 ! The
         ! above overlap times the basis function from atom 2.  This is
         ! accumulated with each iteration of the alpha loop.

   logical :: contrib  ! At least one alpha pair contributes.  For the nuc and
         ! elec calculations this also requires that the pot term contribute.


   ! Variables and data structures that change or are accumulated with each
   !   iteration of the lattice loop.
   complex (kind=double), allocatable, dimension (:,:,:,:) :: currentPair
   complex (kind=double), allocatable, dimension (:,:,:) :: pairXBasisFn12


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

   ! Record the beginning of this phase of the calculation.
   call timeStampStart (32)

   ! Determine if this calculation has already been completed by a previous
   !   OLCAO execution.
   hdf5Status = 0
   attribIntDims(1) = 1
   call h5aread_f(aid,H5T_NATIVE_INTEGER,hdf5Status,attribIntDims,hdferr)
   if (hdferr /= 0) stop 'Failed to read atom KOverlap status.'
   if (hdf5Status == 1) then
      write(20,*) "Two-center KOverlap already exists. Skipping."
      call timeStampEnd(32)
      call h5aclose_f(aid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atom KOverlap status.'
      return
   endif

   ! Allocate space for locally defined allocatable arrays
   allocate (currentBasisFns     (maxNumAtomAlphas,maxNumStates,2))
   allocate (currentAlphas       (maxNumAtomAlphas,2))
   allocate (currentlmAlphaIndex (maxNumAtomAlphas,2))
   allocate (currentlmIndex      (maxNumStates,2))
   allocate (pairXBasisFn2       (16,maxNumAtomAlphas,maxNumStates,3))
   allocate (pairXBasisFn12      (maxNumStates,maxNumStates,3))

   allocate (coreVale            (coreDim,valeDim,numKPoints,3))
   allocate (currentPair         (maxNumStates,maxNumStates,numKPoints,3))


   ! Initialize key matrices
   coreVale      (:,:,:,:) = 0.0_double
   valeVale      (:,:,:,:) = 0.0_double
   coreCore      (:,:,:,:) = 0.0_double

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
!         currentNegligLimit = alphaDist(1,1,currentElements(1),&
!               & currentElements(2))
         currentNegligLimit = logBasisFnThresh * (currentAlphas(1,1) + &
               & currentAlphas(1,2)) / (currentAlphas(1,1) * &
               & currentAlphas(1,2))

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
         currentPair(:,:,:,:) = 0.0_double

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
            ! times the atom2 basis functions in the XYZ direction.
            pairXBasisFn2(:,:currentNumAlphas(1),&
                  & :currentNumTotalStates(2),:) = 0.0_double

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

                  ! At this point a sufficient overlap between atomic
                  !   Gaussians has been shown to exist. FIX: This appears
                  !   to be a useless term. Similarly, the .eqv. test below
                  !   becomes useless.
                  contrib = .true.

                  ! Calculate the opcode to do the correct set of integrals
                  !   for the current alpha pair
                  l1l2Switch = ishft(1,&
                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
                        &+ ishft(16,&
                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))

                  ! Loop over each polarization direction.
                  do axis = 1,3
                     call KOverlap2CIntg (currentAlphas(alphaIndex(1),1), &
                           & currentAlphas(alphaIndex(2),2), &
                           & currentPosition(:,1), shiftedAtomPos(:), &
                           & recipVectors(:,axis) / numAxialKPoints(axis), &
                           & l1l2Switch, oneAlphaSetK(:,:,axis))
                  enddo

                  ! Collect the results of the overlap of the current alpha
                  !   times the basis functions of atom 2.
                  if (contrib .eqv. .true.) then
                     do n = 1, 3
                        do m = 1, currentNumTotalStates(2)
                           pairXBasisFn2(:currentlmAlphaIndex(&
                              & alphaIndex(1),1),alphaIndex(1),m,n) = &
                              & pairXBasisFn2(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),alphaIndex(1),m,n) + &
                              & oneAlphaSetK(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),currentlmIndex(m,2),n) * &
                              & currentBasisFns(alphaIndex(2),m,2)
                        enddo
                     enddo

                     ! Update the maximum alpha used from atom 1.
                     maxAlpha1Used = max(alphaIndex(1),maxAlpha1Used)
                  endif
               enddo
            enddo  ! min number of alphas between the two atoms l


            ! At this point all the alpha loops are complete and we can form
            !   a product with the atom 1 basis function to give the overlap
            !   integral in a complete basis representation.
            do l = 1, 3
               call cmplxMultWithBasisFn1 (currentBasisFns,pairXBasisFn2(:,:,:,l),&
                     & pairXBasisFn12(:,:,l),currentlmIndex,&
                     & currentNumTotalStates,maxAlpha1Used)
            enddo

            ! Collect this atom 1, atom 2 basis function overlap matrix for
            !   all bloch vectors (kpoints) with phase factors appropriate
            !   to the current atom 2 lattice vector.  NOTE that the k
            !   index is over the number of cells in the superlattice.
            do l = 1, 3
               call cmplxApplyPhaseFactors (currentPair(:,:,:,l),&
                     & pairXBasisFn12(1:currentNumTotalStates(1),&
                     & 1:currentNumTotalStates(2),l),&
                     & currentNumTotalStates(1),currentNumTotalStates(2),&
                     & k,0)
            enddo
         enddo !(k superlattice)

         ! At this point all the lattice sums for the current atom pair are
         !   complete.

         ! So now we can arrange the data from this atom into a set of three
         !   large matrices.  A valence-valence matrix, a core-valence matrix,
         !   and a core-core matrix.

         do k = 1, 3
            ! Make a correction for the atom 2 lattice origin shift.
            call kPointLatticeOriginShift (currentNumTotalStates,&
                  & currentPair(:,:,:,k),latticeVector)
            call saveCurrentPair(i,j,numKPoints,currentPair(:,:,:,k),&
                  & valeVale(:,:,:,k),coreVale(:,:,:,k),coreCore(:,:,:,k),0)
         enddo
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


   deallocate (currentPair)

   ! Perform orthogonalization and save the results to disk.
   call ortho(8,packedVVDims,did,aid)

   ! Make a finishing time stamp.
   call timeStampEnd (32)

end subroutine gaussKOverlap
#endif


subroutine ortho (opCode,packedVVDims,did,aid)

   ! Use necessary modules.
   use O_Kinds
   use O_PotTypes, only: potTypes
   use O_KPoints, only: numKPoints
   use O_AtomicSites, only: coreDim, valeDim
   use O_Orthogonalization
   use O_Potential, only: rel
#ifndef GAMMA
   use O_Integrals, only: coreValeOL
#else
   use O_Integrals, only: coreValeOLGamma
#endif

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed dummy parameters.
   integer, intent(in) :: opCode
   integer(hsize_t), dimension(2), intent(in) :: packedVVDims
   integer(hid_t), dimension(numKPoints,3), intent(in) :: did
   integer(hid_t), intent(in) :: aid

   ! Define local variables.
   integer :: i,j,k,l
   integer :: hdferr
   integer :: currIndex
   real (kind=double), allocatable, dimension (:,:,:) :: packedValeVale
   integer(hsize_t), dimension (1) :: attribIntDims ! Attribute dataspace dim

   ! Orthogonalizing against the overlap matrix is unnecessary when the core
   !   dimension is zero.  However, we must still allocate some matrices for
   !   use later.
   if (coreDim /= 0) then
#ifndef GAMMA
      do i = 1, numKPoints
         ! Form product of (valeCoreOL)(coreVale) and (valeCore)(coreValeOL).
         !   Subtract both from the target matrix elements (valeVale).
         do j = 1, 3
            call valeCoreCoreVale (valeDim,coreDim,valeVale(:,:,i,j),&
                  & coreVale(:,:,i,j),coreValeOL(:,:,i))
         enddo
      enddo
#else
      ! Form a product of (valeCoreOL)(coreVale) and (valeCore)
      !   (coreValeOL).  Subtract both from the target matrix elements
      !   (valeValeGamma).
      do j = 1, 3
         call valeCoreCoreValeGamma (valeDim,coreDim,&
               & valeValeGamma(:,:,j),coreValeGamma(:,:,j),&
               & coreValeOLGamma(:,:))
      enddo
#endif
   endif

   ! Allocate space to finish orthogonalization and pack the valeVale matrix.
#ifndef GAMMA
   deallocate (coreVale)
   allocate (valeCore(coreDim,valeDim,3)) ! Pre-transposed format.
   allocate (packedValeVale(2,valeDim*(valeDim+1)/2,3))
#else
   deallocate (coreValeGamma)
   allocate (valeCoreGamma(coreDim,valeDim,3)) ! Pre-transposed format.
   allocate (packedValeVale(1,valeDim*(valeDim+1)/2,3))
#endif

   do i = 1, numKPoints

      ! Orthogonalization is only necessary if the core dimension is non-zero.
      if (coreDim /= 0) then
         do j = 1, 3
#ifndef GAMMA
            ! Form product of (coreValeOL)(coreCore) in temp matrix (valeCore).
            call coreValeCoreCore (valeDim,coreDim,valeCore(:,:,j),&
                  & coreValeOL(:,:,i),coreCore(:,:,i,j))

            ! Finally compute the product of the above
            !   (valeCore)(coreCore) with coreValeOL.
            call makeValeVale (valeDim,coreDim,valeDim,valeCore(:,:,j),&
                  & coreValeOL(:,:,i),valeVale(:,:,i,j),&
                  & packedValeVale(:,:,j),1,0)
#else
            ! Form product of (coreValeOL)(coreCore) in temp matrix (valeCore).
            call coreValeCoreCoreGamma (valeDim,coreDim,valeCoreGamma(:,:,j),&
                  & coreValeOLGamma,coreCoreGamma(:,:,j))

            ! Finally compute the product of the above
            !   (valeCore)(coreCore) with coreValeOL.
            call makeValeValeGamma (valeDim,coreDim,valeDim,&
                  & valeCoreGamma(:,:,j), coreValeOLGamma,&
                  & valeValeGamma(:,:,j), packedValeVale(:,:,j),1,0)
#endif
         enddo
      else

         ! Initialize the index for packing the valeVale matrix.
#ifndef GAMMA
         do j = 1, 3
            currIndex = 0
            do k = 1, valeDim
               do l = 1, k
                  currIndex = currIndex + 1
                  packedValeVale(1,currIndex,j) = &
                        & real(valeVale(l,k,i,j),double)
                  packedValeVale(2,currIndex,j) = aimag(valeVale(l,k,i,j))
               enddo
            enddo
         enddo
#else
         do j = 1, 3
            currIndex = 0
            do k = 1, valeDim
               do l = 1, k
                  currIndex = currIndex + 1
                  packedValeVale(1,currIndex,j) = valeValeGamma(l,k,j)
               enddo
            enddo
         enddo
#endif
      endif

      ! Write the valeVale term onto disk in HDF5 format.
      do j = 1, 3
         call h5dwrite_f(did(i,j),H5T_NATIVE_DOUBLE,packedValeVale(:,:,j),&
               & packedVVDims,hdferr)
         select case (opCode)
         case (6)
            if (hdferr /= 0) stop 'failed to write dipole moment vale vale'
         case (7)
            if (hdferr /= 0) stop 'failed to write momentum matrix vale vale'
         case (8)
            if (hdferr /= 0) stop 'failed to write Koverlap vale vale'
         end select
      enddo
   enddo ! KPoints


   ! Record that the calculation is complete.
   attribIntDims(1) = 1
   call h5awrite_f(aid,H5T_NATIVE_INTEGER,1,attribIntDims,hdferr)
   select case (opCode)
   case (6)
      if (hdferr /= 0) stop 'Failed to record atom DM overlap success.'
   case (7)
      if (hdferr /= 0) stop 'Failed to record atom momentum matrix success.'
   case (8)
      if (hdferr /= 0) stop 'Failed to record atom Koverlap success.'
   end select

   ! Close the attribute.
   call h5aclose_f(aid,hdferr)
   select case (opCode)
   case (6)
      if (hdferr /= 0) stop 'Failed to close atom DM overlap attribute.'
   case (7)
      if (hdferr /= 0) stop 'Failed to close atom momentum matrix attribute.'
   case (8)
      if (hdferr /= 0) stop 'Failed to close atom Koverlap attribute.'
   end select


   ! Deallocate matrices that are no longer necessary.
#ifndef GAMMA
   deallocate (valeCore)
#else
   deallocate (valeCoreGamma)
#endif
   deallocate (packedValeVale)

end subroutine ortho


subroutine cleanUpIntegrals3Terms

   implicit none

#ifndef GAMMA
   deallocate (coreCore)
   deallocate (valeVale)
#else
   deallocate (coreCoreGamma)
   deallocate (valeValeGamma)
#endif

end subroutine cleanUpIntegrals3Terms


end module O_Integrals3Terms
