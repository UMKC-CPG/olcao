module O_IntegralsPSCF

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

subroutine intgAndOrMom(doINTG,doMOME)

   ! Import necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Constants,    only: dim3
   use O_PotTypes,     only: potTypes
   use O_AtomicSites,  only: numAtomSites
   use O_Potential,    only: spin, potCoeffs
   use O_PSCFIntgHDF5, only: targetChunkSize
   use O_Basis,        only: initializeAtomSite
   use O_PotSites,     only: numPotSites, potSites
   use O_AtomicTypes,  only: maxNumAtomAlphas, maxNumStates
   use O_GaussianIntegrals, only: overlapInteg, KEInteg, nucPotInteg, &
         & threeCentInteg, MOMF
   use O_Lattice, only: logBasisFnThresh, numCellsReal, cellSizesReal, &
         & cellDimsReal, findLatticeVector
   use O_ParallelSubs
   use HDF5
   use MPI

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define passed parameters.
   integer :: doINTG  ! (1) = do overlap & hamiltonian integrals; (0) = do not.
   integer, intent(in) :: doMOME

   ! Define variables for MPI
   integer :: maxAtom, minAtom
   integer :: toBalance
   integer :: atomLoop
   integer :: mpiSize, mpiRank, mpiErr

   ! Define local variables for logging, loop control, etc.
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
   real (kind=double), allocatable, dimension (:,:) :: alphaDist ! This
         ! matrix is calculated before the lattice loop is started, but it is
         ! only used within the alpha loop.  It is pre-calculated since it
         ! used multiple times by the same atom pair.
   real (kind=double), allocatable, dimension (:,:) :: alphaCenter ! This
         ! matrix is like the above one except that it provides a factor that
         ! is used to calculate the center point of two overlapping alphas.  It
         ! is also an intermediate step in the calculation of the alphaDist
         ! matrix.
   real (kind=double), dimension (16,16)   :: oneAlphaSet ! The overlap from
         ! one alpha pair. (May be 2-center overlap, 2-center KE, 3-center
         ! nuclear, or 3-center overlap.)
   real (kind=double), dimension (16,16,3) :: oneAlphaSetMom ! The overlap
         ! contribution to the momentum matrix from one alpha pair.
   real (kind=double), dimension (16,16,2) :: potAtomOverlap ! The accumulated
         ! values from all the oneAlphaSet calculations for spin up and down.
   real (kind=double), allocatable, dimension (:,:,:) :: pairXBasesFn2OL ! The
         ! overlap times the basis function from atom 2.  This is accumulated
         ! with each iteration of the alpha loop.
   real (kind=double), allocatable, dimension (:,:,:,:) :: pairXBasisFn2Ham !
         ! The potential overlap times the basis function from atom 2.  This is
         ! accumulated with each iteration of the alpha loop.
   real (kind=double), allocatable, dimension (:,:,:) :: pairXBasisFn2MomX !
         ! The overlap times the basis function from atom 2.  This is
         ! accumulated with each iteration of the alpha loop.
   real (kind=double), allocatable, dimension (:,:,:) :: pairXBasisFn2MomY !
         ! The overlap times the basis function from atom 2.  This is
         ! accumulated with each iteration of the alpha loop.
   real (kind=double), allocatable, dimension (:,:,:) :: pairXBasisFn2MomZ !
         ! The overlap times the basis function from atom 2.  This is
         ! accumulated with each iteration of the alpha loop.


   ! Potential specific variables that change with each potential site 
   !   or potential lattice point iteration.
   integer :: currentPotType
   integer :: currentNumPotAlphas
   integer :: potCoeffIndex
   real (kind=double), dimension (2) :: currentPotCoeff
   real (kind=double) :: currentPotAlpha

   ! Nuclear potential specific variables that will change with each potential
   !   site or potential lattice point iteration.
   real (kind=double) :: zFactor  ! Charge on the current nucleus.
   real (kind=double) :: threeAlphaDist ! A value like the alphaDist above that
         ! can be pre-calculated since part of it will not change with the
         ! actual distance, but only with the choice of alphas (which are known
         ! at input time.)
   real (kind=double), dimension (dim3) :: potPosition ! Original cell
         ! position of the current potential site.
   


   ! Variables and data structures that change or are accumulated with each
   !   iteration of the lattice loop (k).
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:,:) :: currentPair
   complex (kind=double), allocatable, dimension (:,:,:,:) :: currentPairMom
#else
   real (kind=double), allocatable, dimension (:,:,:) :: currentPairGamma
   real (kind=double), allocatable, dimension (:,:,:) :: currentPairGammaMom
#endif
   real (kind=double), allocatable, dimension (:,:,:) :: pairXBasisFn12Ham
   real (kind=double), allocatable, dimension (:,:)   :: pairXBasisFn12OL
   real (kind=double), allocatable, dimension (:,:)   :: pairXBasisFn12MomX
   real (kind=double), allocatable, dimension (:,:)   :: pairXBasisFn12MomY
   real (kind=double), allocatable, dimension (:,:)   :: pairXBasisFn12MomZ



   ! Local position and direction vectors and radii
   real (kind=double), dimension (dim3) :: latticeVector ! Vector to lattice
         ! point closest to the difference between the unit cell positions for
         ! atom 1 and atom 2.
   real (kind=double), dimension (dim3) :: latticeVector2 ! Vector to lattice
         ! point closest to the difference between the unit cell positions for
         ! the center of the alpha1 alpha2 overlap and alpha3 origin.
   real (kind=double), dimension (dim3) :: shiftedAtomPos ! The position of
         ! atom 2 shifted to each relevant lattice point.
   real (kind=double), dimension (dim3) :: shiftedPotPos ! The position of
         ! the nuclear pot shifted to each relevant lattice point.
   real (kind=double), dimension (dim3) :: overlapCenter ! The position of the
         ! center of the overlap between two alphas.  This is used to help get
         ! the distance to the third center (nuclear potential).
   real (kind=double), dimension (dim3) :: centerOriginVect ! The seperation
         ! vector between the center of the two alpha overlaps, and the origin
         ! to be used for the nuclear potential site.
   real (kind=double) :: centerOriginSep ! The sum of squares of the above
         ! vector.
   real (kind=double) :: shiftedCenterOriginSep ! The above value shifted
         ! according to the current lattice point loop.
   real (kind=double) :: atomSiteSepSqrd ! The square of the minimum distance
         ! seperating atom 1 and atom 2 according to their unit cell positions
         ! shifted by the lattice point closest to their difference.
   real (kind=double) :: shiftedAtomSiteSep ! The seperation distance between
         ! atom 1 and the shifted position of atom 2.
   real (kind=double) :: currentNegligLimit ! The distance beyond which all
         ! alpha pairs are considered to have no overlap.
   real (kind=double) :: maxLatticeRadius ! Maximum radius beyond which no
         ! lattice points will be considered for integration.  This is used
         ! for the atomic alpha pairs.
   real (kind=double) :: maxLatticeRadius2 ! Maximum radius beyond which no
         ! lattice points will be considered for integration.  This is used
         ! for the nuclear alpha triple.

   ! Define variables for gauss integrals
   integer :: l1l2Switch
   integer, dimension(16) :: powerOfTwo = (/0,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3/)


   ! Log the date and time we start.
   call timeStampStart(20)


   ! Allocate space for locally defined allocatable arrays
   allocate (currentBasisFns       (maxNumAtomAlphas,maxNumStates,2))
   allocate (currentAlphas         (maxNumAtomAlphas,2))
   allocate (currentlmAlphaIndex   (maxNumAtomAlphas,2))
   allocate (currentlmIndex        (maxNumStates,2))
   allocate (alphaDist             (maxNumAtomAlphas,maxNumAtomAlphas))
   allocate (alphaCenter           (maxNumAtomAlphas,maxNumAtomAlphas))
   if (doINTG == 1) then
#ifndef GAMMA
      allocate (currentPair           (maxNumStates,maxNumStates,numKPoints,2))
#else
      allocate (currentPairGamma      (maxNumStates,maxNumStates,2))
#endif
   endif
   if (doMOME == 1) then
#ifndef GAMMA
      allocate (currentPairMom        (maxNumStates,maxNumStates,numKPoints,3))
#else
      allocate (currentPairGammaMom   (maxNumStates,maxNumStates,3))
#endif
   endif

   if (doINTG == 1) then
      allocate (pairXBasisFn2OL    (16,maxNumAtomAlphas,maxNumStates))
      allocate (pairXBasisFn2Ham   (16,maxNumAtomAlphas,maxNumStates,spin))
      allocate (pairXBasisFn12OL   (maxNumStates,maxNumStates))
      allocate (pairXBasisFn12Ham  (maxNumStates,maxNumStates,spin))

      ! Initialize overlap and hamiltonian matrices.
      pairXBasisFn12OL(:,:)    = 0.0_double
      pairXBasisFn12Ham(:,:,:) = 0.0_double
   endif

   if (doMOME == 1) then
      allocate (pairXBasisFn2MomX   (16,maxNumAtomAlphas,maxNumStates))
      allocate (pairXBasisFn2MomY   (16,maxNumAtomAlphas,maxNumStates))
      allocate (pairXBasisFn2MomZ   (16,maxNumAtomAlphas,maxNumStates))
      allocate (pairXBasisFn12MomX  (maxNumStates,maxNumStates))
      allocate (pairXBasisFn12MomY  (maxNumStates,maxNumStates))
      allocate (pairXBasisFn12MomZ  (maxNumStates,maxNumStates))

      ! Initialize momentum matrix element matrices.
      pairXBasisFn12MomX(:,:) = 0.0_double
      pairXBasisFn12MomY(:,:) = 0.0_double
      pairXBasisFn12MomZ(:,:) = 0.0_double
   endif


   ! Compute the parameters for load balancing the following loop. Each
   !   process needs to have about the same number of atom pairs, but
   !   each atom pair (a,b) only needs to be done once. (I.e. no (b,a)).
   !   Note that the numbers stored in minAtom and maxAtom are from within
   !   the range 1 to numAtomSites*(numAtomSites+1)/2, *not* 1 to numAtomSites.
   toBalance = (numAtomSites * (numAtomSites + 1))/2
   call MPI_Comm_Size(MPI_COMM_WORLD,mpiSize,mpiErr)
   call MPI_Comm_Rank(MPI_COMM_WORLD,mpiRank,mpiErr)
   call loadBalMPI(toBalance,minAtom,maxAtom,mpiRank,mpiSize)

   ! Begin atom-atom loops
   do atomLoop = minAtom, maxAtom

      ! Get the actual atom number of the atoms in this atomLoop pair.
      call finUnpackedIndices(atomLoop,i,j)

      ! Obtain local copies of key data from larger global data structures for
      !   the first looped atom.
      call initializeAtomSite(i,1,currentAtomType,currentElements,&
            & currentNumTotalStates,currentNumCoreStates,currentNumValeStates,&
            & currentNumAlphas,currentlmIndex,currentlmAlphaIndex,&
            & currentPosition,currentAlphas,currentBasisFns)

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

      ! Calculate the maximum distance from either atom where the overlap
      !   is considered to be non-negligable.
      currentNegligLimit = logBasisFnThresh * (currentAlphas(1,1) + &
            & currentAlphas(1,2)) / (currentAlphas(1,1) * &
            & currentAlphas(1,2))

      ! Determine if there are no alpha terms for this atom pair that fall
      !   within the current negligability limit.  Cycle if there are none.
      if (atomSiteSepSqrd > currentNegligLimit) cycle

      ! Determine the maximum radius beyond which no lattice point will be
      !   considered to contribute to the overlap integral for this atom
      !   pair.
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

      ! Form an alpha overlap matrix that is later used to determine whether
      !   a particular alpha pair needs to be considered in the overlap
      !   matrix or if the alpha pair do not have sufficient overlap.
      ! Another possibility (that can be applied to other similar algorithms
      !   in other subroutines) is to calculate these values based on the
      !   element as opposed to the current method of calculating them based
      !   on every atom pair iteration. In this way the values would only
      !   have to be calculated once for each element pair.
      do k = 1,currentNumAlphas(2)
         alphaCenter(:currentNumAlphas(1),k) = &
               &  currentAlphas(:currentNumAlphas(1),1) / &
               & (currentAlphas(:currentNumAlphas(1),1) + &
               &  currentAlphas(k,2))
         alphaDist(:currentNumAlphas(1),k) = logBasisFnThresh / &
               &  currentAlphas(k,2) / alphaCenter(:currentNumAlphas(1),k)
      enddo

      ! Initialize the matrix for this atom pair that will record the
      !   overlap integrals with bloch vector (kpoint) phase factors.
      if (doINTG == 1) then
#ifndef GAMMA
         currentPair(:,:,:,:) = 0.0_double
#else
         currentPairGamma(:,:,:) = 0.0_double
#endif
      endif
      if (doMOME == 1) then
#ifndef GAMMA
         currentPairMom(:,:,:,:) = 0.0_double
#else
         currentPairGammaMom(:,:,:) = 0.0_double
#endif
      endif


      ! Begin a loop over all the lattice points to shift the position of
      !   atom number 2 to all the replicated cells.
      do k = 1, numCellsReal

         ! Exit the loop when we have exceeded the necessary number of
         !   lattice points based on distance.
         if (cellSizesReal(k) > maxLatticeRadius) exit

         ! Obtain the position of atom #2 shifted by the current lattice.
         shiftedAtomPos(:) = currentPosition(:,2) + latticeVector(:) + &
               & cellDimsReal(:,k)

         ! Obtain the seperation vector between atom 1 and the shifted
         !   position of atom 2.
         shiftedAtomSiteSep = sum ((currentPosition(:,1) - &
               & shiftedAtomPos(:))**2)

         ! Determine if this shifted atom position puts it outside of the
         !   above determined negligability limit for this atom pair.
         if (shiftedAtomSiteSep > currentNegligLimit) cycle


         ! Initialize the matrices to hold the product of the integrals
         !    times the atom2 wave functions.
         if (doINTG == 1) then
            pairXBasisFn2OL(:,:currentNumAlphas(1),&
                  &:currentNumTotalStates(2)) = 0.0_double
            pairXBasisFn2Ham(:,:currentNumAlphas(1),&
                  &:currentNumTotalStates(2),:spin) = 0.0_double
         endif
         if (doMOME == 1) then
            pairXBasisFn2MomX(:,:currentNumAlphas(1),&
                  &:currentNumTotalStates(2)) = 0.0_double
            pairXBasisFn2MomY(:,:currentNumAlphas(1),&
                  &:currentNumTotalStates(2)) = 0.0_double
            pairXBasisFn2MomZ(:,:currentNumAlphas(1),&
                  &:currentNumTotalStates(2)) = 0.0_double
         endif

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
            if (alphaDist(l,l) < shiftedAtomSiteSep) exit

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
               if (alphaDist(alphaIndex(1),alphaIndex(2)) < &
                     & shiftedAtomSiteSep) then

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

               if (doINTG == 1) then
                  ! At this point a sufficient overlap has been found for
                  !   the current alpha pair so we start looking through
                  !   nuclear potentials.

                  ! First we find the center of the overlap between the two
                  !   alphas.
                  overlapCenter(:) = &
                        & alphaCenter(alphaIndex(1),alphaIndex(2)) * &
                        & (currentPosition(:,1) - shiftedAtomPos(:)) + &
                        & shiftedAtomPos(:)

                  ! Initialize the result matrix to zero before starting the
                  !   nuclear potential loop.
                  potAtomOverlap(:currentlmAlphaIndex(alphaIndex(1),1), &
                        & :currentlmAlphaIndex(alphaIndex(2),2),:spin) = &
                        & 0.0_double


                  ! Initialize a counter to track which potential
                  !    coefficient we are currently calculating on.
                  potCoeffIndex = 0

                  ! Initiate a loop over each potential site for the nuclear
                  !   potentials.
                  do m = 1, numPotSites

                     ! Skip equivalent types now.  THey will be accounted
                     !   for later but we want to do some setup stuff only
                     !   once for all atoms of the same type.
                     if (potSites(m)%firstPotType == 0) cycle

                     ! Initialize the parameters for this potential site
                     !   related to the type of this site.
                     currentPotType = potSites(m)%potTypeAssn
                     currentNumPotAlphas = potTypes(currentPotType)%numAlphas

                     ! Initiate a loop over all the potential alphas present
                     !   for this site including the nuclear contribution.
                     do n = 1, currentNumPotAlphas + 1

                        ! Assign the potential alpha based on the value of n.
                        if (n <= currentNumPotAlphas) then

                           ! Apply the case for the atomic potential.

                           ! Increment the index value that indicates which
                           !   potential coefficient to use of all
                           !   coefficients in the system.
                           potCoeffIndex = potCoeffIndex + 1

                           ! Store the current potential coefficient.
                           currentPotCoeff(:spin) = &
                                 & potCoeffs(potCoeffIndex,:spin)

                           ! Store the current potential alpha. 
                           currentPotAlpha = &
                                 & potTypes(currentPotType)%alphas(n)
                        else

                           ! Apply the case for the nuclear potential.

                           ! Get nuclear charge associated with this type.
                           zFactor = potTypes(currentPotType)%nucCharge

                           ! If the zFactor is sufficiently small we
                           !   consider the effect of the overlap to be
                           !   negligable and we cycle to the next one.
                           if (zFactor < smallThresh) cycle

                           ! Get the exponential alpha factor for the
                           !   nuclear potential.
                           currentPotAlpha = &
                                 & potTypes(currentPotType)%nucAlpha

                        endif

                        ! Determine the maximum distance beyond which the
                        !   overlap of the three gaussians is considered
                        !   negligable.
                        threeAlphaDist = logBasisFnThresh * (1 - &
                              & shiftedAtomSiteSep / &
                              & alphaDist(alphaIndex(1),alphaIndex(2))) * &
                              & (currentAlphas(alphaIndex(1),1) + &
                              &  currentAlphas(alphaIndex(2),2) + &
                              &  currentPotAlpha) / &
                              & (currentAlphas(alphaIndex(1),1) + &
                              &  currentAlphas(alphaIndex(2),2)) / &
                              &  currentPotAlpha

                        ! Loop over the remaining potential sites that are
                        !   equivalent.
                        do o = 0, potTypes(currentPotType)%multiplicity-1

                           ! Initialize the parameters for this potential
                           !   site related to its position.
                           potPosition(:) = potSites(m+o)%cartPos(:)

                           ! Locate the origin for the potential lattice sum.
                           call findLatticeVector((overlapCenter(:) - &
                                 & potPosition(:)), latticeVector2)

                           ! Find the seperation vector and distance between
                           !   the minimum overlap center, and the origin.
                           centerOriginVect(:) = overlapCenter(:) - &
                                 & potPosition(:) - latticeVector2(:)
                           centerOriginSep = sum(centerOriginVect(:)**2)

                           ! Check if largest potential is negligable or not.
                           if (centerOriginSep > threeAlphaDist) cycle

                           ! At this point it must be the case that at least
                           !   one contribution will be calculated.

                           ! First, find the cut-off radius for the potential
                           !   lattice summation by the triangle inequality.
                           maxLatticeRadius2 = centerOriginSep + &
                              & threeAlphaDist + 2.0_double * &
                              & sqrt (centerOriginSep * threeAlphaDist)


                           ! Begin loop over lattice points for the site pot.
                           do p = 1, numCellsReal

                              ! Check if this lattice point extends beyond
                              !   the range of negligability
                              if (cellSizesReal(p) > maxLatticeRadius2) exit

                              ! Get overlap center shifted by lattice point.
                              shiftedCenterOriginSep = 
                                    & sum((centerOriginVect(:) - &
                                    & cellDimsReal(:,p))**2)

                              ! Check if shifted seperation between the
                              !   center and the origin extends past the
                              !   negligability limit.  If so, cycle to the
                              !   next cell.
                              if (shiftedCenterOriginSep > threeAlphaDist) cycle

                              ! Get seperation shifted by the lattice point.
                              shiftedPotPos(:) = potPosition(:) + &
                                    & latticeVector2(:) + cellDimsReal(:,p)

                              if (n <= currentNumPotAlphas) then
                                 ! Calculate the opcode to do the correct
                                 !   set of integrals for the current alpha
                                 !   pair.
                                 l1l2Switch = ishft(1,(powerOfTwo(&
                                      & currentlmAlphaIndex(&
                                      & alphaIndex(1),1))))+ ishft(16,&
                                      & (powerOfTwo(currentlmAlphaIndex(&
                                      & alphaIndex(2),2))))

                                 call threeCentInteg (&
                                    & currentAlphas(alphaIndex(1),1),&
                                    & currentAlphas(alphaIndex(2),2),&
                                    & currentPotAlpha, currentPosition(:,1),&
                                    & shiftedAtomPos(:), shiftedPotPos(:),&
                                    & l1l2Switch, oneAlphaSet)

                                 ! Accumulate results returned for alpha set.
                                 do q = 1, spin
                                 potAtomOverlap(:currentlmAlphaIndex &
                                    & (alphaIndex(1),1),:currentlmAlphaIndex&
                                    & (alphaIndex(2),2),q) = &
                                    & potAtomOverlap(:currentlmAlphaIndex &
                                    & (alphaIndex(1),1),:currentlmAlphaIndex&
                                    & (alphaIndex(2),2),q) + &
                                    & oneAlphaSet(:currentlmAlphaIndex &
                                    & (alphaIndex(1),1),:currentlmAlphaIndex&
                                    & (alphaIndex(2),2)) * &
                                    & currentPotCoeff(q)
                                 enddo
                              else
                                ! Calculate the opcode to do the correct set
                                ! of integrals for the current alpha pair
                                l1l2Switch = ishft(1,(powerOfTwo(&
                                  &currentlmAlphaIndex(alphaIndex(1),1))))&
                                  &+ ishft(16,&
                                  &(powerOfTwo(&
                                  &currentlmAlphaIndex(alphaIndex(2),2))))
                                
                                call nucPotInteg (&
                                    & currentAlphas(alphaIndex(1),1),&
                                    & currentAlphas(alphaIndex(2),2),&
                                    & currentPotAlpha,currentPosition(:,1),&
                                    & shiftedAtomPos(:),shiftedPotPos(:),&
                                    & l1l2Switch,oneAlphaSet)

                                 ! Accumulate results returned for alpha set.
                                 do q = 1, spin
                                    potAtomOverlap(:currentlmAlphaIndex &
                                       & (alphaIndex(1),1),&
                                       & :currentlmAlphaIndex &
                                       & (alphaIndex(2),2),q) = &
                                       & potAtomOverlap(:currentlmAlphaIndex&
                                       & (alphaIndex(1),1),&
                                       & :currentlmAlphaIndex &
                                       & (alphaIndex(2),2),q) - &
                                       & oneAlphaSet(:currentlmAlphaIndex &
                                       & (alphaIndex(1),1),&
                                       & :currentlmAlphaIndex &
                                       & (alphaIndex(2),2)) * zFactor
                                 enddo
                              endif
                           enddo ! (p numCells)
                        enddo ! (o multiplicity)
                     enddo ! (n numCurrentPotAlphas+1)
                  enddo ! (m numPots (inequivalent))


                  ! Calculate the opcode to do the correct set of integrals
                  ! for the current alpha pair
                  l1l2Switch = ishft(1,&
                       &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
                       &+ ishft(16,&
                       &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))

                  ! Determine the kinetic energy contribution
                  call KEInteg (currentAlphas(alphaIndex(1),1),&
                        & currentAlphas(alphaIndex(2),2),&
                        & currentPosition(:,1), shiftedAtomPos(:),&
                        & l1l2Switch, oneAlphaSet)

                  ! Accumulate the contribution from this alpha pair
                  do m = 1, spin
                     potAtomOverlap(:currentlmAlphaIndex &
                           & (alphaIndex(1),1),:currentlmAlphaIndex &
                           & (alphaIndex(2),2),m) = &
                           & potAtomOverlap(:currentlmAlphaIndex &
                           & (alphaIndex(1),1),:currentlmAlphaIndex &
                           & (alphaIndex(2),2),m) + &
                           & oneAlphaSet(:currentlmAlphaIndex &
                           & (alphaIndex(1),1),:currentlmAlphaIndex &
                           & (alphaIndex(2),2))
                  enddo

                  ! Calculate the opcode to do the correct set of integrals
                  ! for the current alpha pair
                  l1l2Switch = ishft(1,&
                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
                        &+ ishft(16,&
                        &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))
             
                  ! We can proceed with the next step of the calculation.
                  ! This is the actual 2-center overlap integral.
                  call overlapInteg (currentAlphas(alphaIndex(1),1),&
                        & currentAlphas(alphaIndex(2),2), &
                        & currentPosition(:,1), shiftedAtomPos(:),&
                        & l1l2Switch, oneAlphaSet)

               endif ! doINTG


               ! Compute the momentum matrix values if requested.
               if (doMOME == 1) then
                  ! Calculate the opcode to do the correct set
                  ! of integrals for the current alpha pair.
                  l1l2Switch = ishft(1,&
                     & (powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
                     & + ishft(16,&
                     & (powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))

                  call MOMF (currentAlphas(alphaIndex(1),1),&
                        & currentAlphas(alphaIndex(2),2),&
                        & currentPosition(:,1), shiftedAtomPos(:),&
                        & l1l2Switch,oneAlphaSetMom)
               endif

               ! Collect the results of the overlap of the current alpha
               !   times the basis function coefficients for atom 2.

               if (doINTG == 1) then
                  ! Potential overlaps first (if any were found).
                  do q = 1, spin
                     do m = 1, currentNumTotalStates(2)
                        pairXBasisFn2Ham(:currentlmAlphaIndex( &
                              & alphaIndex(1),1),alphaIndex(1),m,q) = &
                              & pairXBasisFn2Ham(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),alphaIndex(1),m,q) + &
                              & potAtomOverlap(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),currentlmIndex(m,2),q) * &
                              & currentBasisFns(alphaIndex(2),m,2)
                     enddo
                  enddo
                  ! Atom pair overlaps second.
                  do m = 1, currentNumTotalStates(2)
                     pairXBasisFn2OL(:currentlmAlphaIndex(alphaIndex(1),1),&
                           & alphaIndex(1),m) = &
                           & pairXBasisFn2OL(:currentlmAlphaIndex &
                           & (alphaIndex(1),1),alphaIndex(1),m) + &
                           & oneAlphaSet(:currentlmAlphaIndex &
                           & (alphaIndex(1),1),currentlmIndex(m,2)) * &
                           & currentBasisFns(alphaIndex(2),m,2)
                  enddo
               endif

               ! Momentum matrix last.  (If requested)
               if (doMOME == 1) then
                  do m = 1, currentNumTotalStates(2)
                     pairXBasisFn2MomX(:currentlmAlphaIndex(alphaIndex(1),&
                           & 1), alphaIndex(1),m) = &
                           & pairXBasisFn2MomX(:currentlmAlphaIndex &
                           & (alphaIndex(1),1),alphaIndex(1),m) + &
                           & oneAlphaSetMom(:currentlmAlphaIndex &
                           & (alphaIndex(1),1),currentlmIndex(m,2),1) * &
                           & currentBasisFns(alphaIndex(2),m,2)
                     pairXBasisFn2MomY(:currentlmAlphaIndex(alphaIndex(1),&
                           & 1), alphaIndex(1),m) = &
                           & pairXBasisFn2MomY(:currentlmAlphaIndex &
                           & (alphaIndex(1),1),alphaIndex(1),m) + &
                           & oneAlphaSetMom(:currentlmAlphaIndex &
                           & (alphaIndex(1),1),currentlmIndex(m,2),2) * &
                           & currentBasisFns(alphaIndex(2),m,2)
                     pairXBasisFn2MomZ(:currentlmAlphaIndex(alphaIndex(1),&
                           & 1), alphaIndex(1),m) = &
                           & pairXBasisFn2MomZ(:currentlmAlphaIndex &
                           & (alphaIndex(1),1),alphaIndex(1),m) + &
                           & oneAlphaSetMom(:currentlmAlphaIndex &
                           & (alphaIndex(1),1),currentlmIndex(m,2),3) * &
                           & currentBasisFns(alphaIndex(2),m,2)
                  enddo
               endif

               ! Update the maximum alpha used from atom 1.
               maxAlpha1Used = max(alphaIndex(1),maxAlpha1Used)

               ! Switch mode from the initial state of no searching to the
               !   state of searching along alphas from atom 1.
               if (currentMode == 0) then
                  currentMode = 1
               endif
            enddo
         enddo  ! min number of alphas between the two atoms l

         ! At this point all the alpha loops are complete and we can form a
         !   product with the atom 1 basis function coeffs to give the overlap
         !   integral in a complete representation.
         if (doINTG == 1) then
            do l = 1, currentNumTotalStates(2)
               do m = 1, currentNumTotalStates(1)
                  pairXBasisFn12OL(m,l) = &
                        & sum(currentBasisFns(:maxAlpha1Used,m,1) * &
                        & pairXBasisFn2OL(currentlmIndex(m,1),:maxAlpha1Used,l))
               enddo
            enddo
            do q = 1, spin
               do l = 1, currentNumTotalStates(2)
                  do m = 1, currentNumTotalStates(1)
                     pairXBasisFn12Ham(m,l,q) = &
                           & sum(currentBasisFns(:maxAlpha1Used,m,1) * &
                           & pairXBasisFn2Ham(currentlmIndex(m,1),&
                           & :maxAlpha1Used,l,q))
                  enddo
               enddo
            enddo
         endif

         ! The momentum is summed against the wave function 1 if needed.
         if (doMOME == 1) then
            do l = 1, currentNumTotalStates(2)
               do m = 1, currentNumTotalStates(1)
               pairXBasisFn12MomX(m,l) = &
                     & sum(currentBasisFns(:maxAlpha1Used,m,1) * &
                     & pairXBasisFn2MomX(currentlmIndex(m,1),:maxAlpha1Used,l))
               pairXBasisFn12MomY(m,l) = &
                     & sum(currentBasisFns(:maxAlpha1Used,m,1) * &
                     & pairXBasisFn2MomY(currentlmIndex(m,1),:maxAlpha1Used,l))
               pairXBasisFn12MomZ(m,l) = &
                     & sum(currentBasisFns(:maxAlpha1Used,m,1) * &
                     & pairXBasisFn2MomZ(currentlmIndex(m,1),:maxAlpha1Used,l))
               enddo
            enddo
         endif

         ! Collect this atom 1, atom 2 basis function overlap matrix for all
         !   bloch vectors (kpoints) with phase factors appropriate to the
         !   current atom 2 lattice vector.  NOTE that the k index is over
         !   the number of cells in the superlattice.
#ifndef GAMMA
         if (doINTG == 1) then
            call applyPhaseFactors (currentPair(:,:,:,1),pairXBasisFn12OL(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),k,0,0)
            call applyPhaseFactors (currentPair(:,:,:,2),pairXBasisFn12Ham(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),k,0,0)
         endif
         if (doMOME == 1) then
            call applyPhaseFactors (currentPairMom(:,:,:,1),&
                  & pairXBasisFn12MomX(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),k,0,0)
            call applyPhaseFactors (currentPairMom(:,:,:,2),&
                  & pairXBasisFn12MomY(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),k,0,0)
            call applyPhaseFactors (currentPairMom(:,:,:,3),&
                  & pairXBasisFn12MomZ(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),k,0,0)
         endif
#else
         if (doINTG == 1) then
            call applyPhaseFactorsGamma (currentPairGamma(:,:,1),&
                  & pairXBasisFn12OL(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),0)
            call applyPhaseFactorsGamma (currentPairGamma(:,:,2),&
                  & pairXBasisFn12Ham(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),0)
         endif
         if (doMOME == 1) then
            call applyPhaseFactorsGamma (currentPairGammaMom(:,:,1),&
                  & pairXBasisFn12MomX(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),0)
            call applyPhaseFactorsGamma (currentPairGammaMom(:,:,2),&
                  & pairXBasisFn12MomY(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),0)
            call applyPhaseFactorsGamma (currentPairGammaMom(:,:,3),&
                  & pairXBasisFn12MomZ(1:&
                  & currentNumTotalStates(1),1:currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),0)
         endif
#endif
      enddo !(k superlattice)

      ! At this point all the lattice sums for the current atom pair are
      !   complete.

      ! So now we can arrange the data from this atom into a set of three
      !   large matrices.  A valence-valence matrix, a core-valence matrix,
      !   and a core-core matrix.

#ifndef GAMMA
      ! First we must make a correction for the atom 2 lattice origin shift.
      if (doINTG == 0) then
         call kPointLatticeOriginShift (currentNumTotalStates,&
               & currentPair(:,:,:,1),latticeVector,numKPoints,0)
         call kPointLatticeOriginShift (currentNumTotalStates,&
               & currentPair(:,:,:,2),latticeVector,numKPoints,0)
         call saveCurrentPair(i,j,numKPoints,currentPair(:,:,:,1),descriptVV,&
               & descriptCC,descriptCV_OL,descriptVC_OL,localVV
               & ga_coreValeOL,ga_coreCore,currentNumTotalStates)
         call saveCurrentPair(i,j,numKPoints,currentPair(:,:,:,2),ga_valeVale,&
               & ga_coreValeOL,ga_coreCore,currentNumTotalStates)
      endif
      if (doMOME == 0) then
         call kPointLatticeOriginShift (currentNumTotalStates,&
               & currentPairMom(:,:,:,1),latticeVector,numKPoints,0)
         call kPointLatticeOriginShift (currentNumTotalStates,&
               & currentPairMom(:,:,:,2),latticeVector,numKPoints,0)
         call kPointLatticeOriginShift (currentNumTotalStates,&
               & currentPairMom(:,:,:,3),latticeVector,numKPoints,0)
         call saveCurrentPair(i,j,numKPoints,currentPairMom(:,:,:,1),&
               & ga_valeVale,ga_coreValeOL,ga_coreCore,currentNumTotalStates)
         call saveCurrentPair(i,j,numKPoints,currentPairMom(:,:,:,2),&
               & ga_valeVale,ga_coreValeOL,ga_coreCore,currentNumTotalStates)
         call saveCurrentPair(i,j,numKPoints,currentPairMom(:,:,:,3),&
               & ga_valeVale,ga_coreValeOL,ga_coreCore,currentNumTotalStates)
      endif
#else
      if (doINTG == 1) then
         call saveCurrentPairGamma(i,j,currentPairGamma(:,:,1),ga_valeVale,&
               & ga_coreValeOL,ga_coreCore,currentNumTotalStates)
         call saveCurrentPairGamma(i,j,currentPairGamma(:,:,2),ga_valeVale,&
               & ga_coreValeOL,ga_coreCore,currentNumTotalStates)
      endif
      if (doMOME == 1) then
         call saveCurrentPairGamma(i,j,currentPairGammaMom(:,:,1),ga_valeVale,&
               & ga_coreValeOL,ga_coreCore,currentNumTotalStates)
         call saveCurrentPairGamma(i,j,currentPairGammaMom(:,:,2),ga_valeVale,&
               & ga_coreValeOL,ga_coreCore,currentNumTotalStates)
         call saveCurrentPairGamma(i,j,currentPairGammaMom(:,:,3),ga_valeVale,&
               & ga_coreValeOL,ga_coreCore,currentNumTotalStates)
      endif
#endif

      ! Mark the completion of this atom (and thus approximately the completion
      !   of one atom by each of the other processes).
      if (mpiRank == 0) then
         do j = 1, mpiSize
            if (mod(i+j-1,10) .eq. 0) then
               write (20,ADVANCE="NO",FMT="(a1)") "|"
            else
               write (20,ADVANCE="NO",FMT="(a1)") "."
            endif
            if (mod(i+j-1,50) .eq. 0) then
               write (20,*) " ",i*50
            endif
            call flush (20)
         enddo
      endif

   enddo    ! (Atom pair loop)


   ! Deallocate all arrays and matrices before exiting.
   deallocate (currentBasisFns)
   deallocate (currentAlphas)
   deallocate (currentlmAlphaIndex)
   deallocate (currentlmIndex)
   deallocate (alphaDist)   ! Can be saved from before
   deallocate (alphaCenter) ! Can be saved for later.
   if (doINTG == 1) then
#ifndef GAMMA
      deallocate (currentPair)
#else
      deallocate (currentPairGamma)
#endif
   endif
   if (doMOME == 1) then
#ifndef GAMMA
      deallocate (currentPairMom)
#else
      deallocate (currentPairGammaMom)
#endif
   endif
   if (doINTG == 1) then
      deallocate (pairXBasisFn2OL)
      deallocate (pairXBasisFn2Ham)
      deallocate (pairXBasisFn12OL)
      deallocate (pairXBasisFn12Ham)
   endif

   if (doMOME == 1) then
      deallocate (pairXBasisFn2MomX)
      deallocate (pairXBasisFn2MomY)
      deallocate (pairXBasisFn2MomZ)
      deallocate (pairXBasisFn12MomX)
      deallocate (pairXBasisFn12MomY)
      deallocate (pairXBasisFn12MomZ)
   endif

   ! Perform orthogonalization and save certain results to disk.
   if (doINTG == 1) then

      ! Orthogonalize the overlap and then the Hamiltonian.
      call orthoOL(ga_valeVale,ga_coreValeOL,ga_coreCore,ga_valeCore,1)
      call ortho(ga_valeVale,ga_coreValeOL,ga_coreCore,ga_valeCore,1)

      ! Save the overlapVV and the overlapCV to disk.
      call saveMatrix ()
      call saveMatrix ()
   endif
   if (doMOME == 1) then
      ! Check if we need to read the coreValeOL orthogonalization coeffs first.
      if (doINTG == 0) then
         call readCoreValeOL ()
      endif
      ! Orthogonalize the momentumXYZ
      call ortho(ga_valeVale,ga_coreValeOL,ga_coreCore,ga_valeCore,1)
      call ortho(ga_valeVale,ga_coreValeOL,ga_coreCore,ga_valeCore,1)
      call ortho(ga_valeVale,ga_coreValeOL,ga_coreCore,ga_valeCore,1)

      ! Save the momentumXYZ to disk
      call saveMatrix ()
      call saveMatrix ()
      call saveMatrix ()
   endif

   ! Log the date and time we finished.
   call timeStampEnd(20)

end subroutine intgAndOrMom


end module O_IntegralsPSCF
