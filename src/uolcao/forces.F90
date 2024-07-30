module O_Force

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

   ! Define the matrix that will hold all integrals of the Hamiltonian.
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:,:,:) :: valeValeF
   complex (kind=double), allocatable, dimension (:,:,:,:,:) :: coreCoreF
   complex (kind=double), allocatable, dimension (:,:,:,:,:) :: coreValeF
#else
   real (kind=double), allocatable, dimension (:,:,:,:) :: valeValeFGamma
   real (kind=double), allocatable, dimension (:,:,:,:) :: coreCoreFGamma
   real (kind=double), allocatable, dimension (:,:,:,:) :: coreValeFGamma
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains


subroutine allocateIntegralsForce(coreDim,valeDim,numKPoints,spin)

   implicit none

   ! Define passed dummy parameters.
   integer, intent(in) :: coreDim
   integer, intent(in) :: valeDim
   integer, intent(in) :: numKPoints
   integer, intent(in) :: spin

#ifndef GAMMA
   allocate (valeValeF (valeDim,valeDim,numKPoints,spin,3))
   allocate (coreCoreF (coreDim,coreDim,numKPoints,spin,3))
   allocate (coreValeF (coreDim,valeDim,numKPoints,spin,3))
#else
   allocate (valeValeFGamma (valeDim,valeDim,spin,3))
   allocate (coreCoreFGamma (coreDim,coreDim,spin,3))
   allocate (coreValeFGamma (coreDim,valeDim,spin,3))
#endif

end subroutine allocateIntegralsForce



subroutine computeForceIntg(totalEnergy)

   ! Import necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Constants,    only: dim3
   use O_KPoints, only: numKPoints
   use O_PotTypes,     only: potTypes
   use O_AtomicSites,  only: numAtomSites, valeDim, coreDim
   use O_Potential,    only: rel, spin, potCoeffs
   use O_Basis,        only: initializeAtomSite
   use O_PotSites,     only: numPotSites, potSites
   use O_AtomicTypes,  only: maxNumAtomAlphas, maxNumStates
   use O_GaussianIntegrals, only: dkinetic2CIntg, &
         & dnuclear3CIntgCB, delectron3CIntgCB, momentum2CIntg
   use O_Lattice, only: logBasisFnThresh, numCellsReal, cellSizesReal, &
         & cellDimsReal, findLatticeVector
   use O_IntgSaving
   use HDF5

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define passed variables.
   real (kind=double) :: totalEnergy

   ! Define local variables for logging, loop control, and file writing.
   integer :: i,j,k,l,m,n,o,p,q,r ! Loop index variables

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
   real (kind=double), dimension (16,16,3) :: oneAlphaSet ! The potential
         ! overlap from one alpha pair.
   real (kind=double), dimension (16,16,2,3) :: potAtomOverlap ! The
         ! accumulated values from all the oneAlphaSet calculations.
   real (kind=double), allocatable, dimension (:,:,:,:,:) :: pairXBasisFn2dHam
         ! The potential overlap times the wave function from atom 2.  This is
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
   

   ! Multiplication of the pairXBasisFn2dHam by the basis coefficients from
   !   basis function 1.
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:,:,:) :: currentPair
#else
   real (kind=double), allocatable, dimension (:,:,:,:) :: currentPairGamma
#endif
   real (kind=double), allocatable, dimension (:,:,:,:) :: pairXBasisFn12dHam


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


#ifndef GAMMA
   allocate (valeValeF (valeDim,valeDim,numKPoints,spin,3))
   allocate (coreCoreF (coreDim,coreDim,numKPoints,spin,3))
   allocate (coreValeF (coreDim,valeDim,numKPoints,spin,3))
#else
   allocate (valeValeFGamma (valeDim,valeDim,spin,3))
   allocate (coreCoreFGamma (coreDim,coreDim,spin,3))
   allocate (coreValeFGamma (coreDim,valeDim,spin,3))
#endif


   ! Allocate space for locally defined allocatable arrays
   allocate (currentBasisFns       (maxNumAtomAlphas,maxNumStates,2))
   allocate (currentAlphas         (maxNumAtomAlphas,2))
   allocate (currentlmAlphaIndex   (maxNumAtomAlphas,2))
   allocate (currentlmIndex        (maxNumStates,2))
   allocate (alphaDist             (maxNumAtomAlphas,maxNumAtomAlphas))
   allocate (alphaCenter           (maxNumAtomAlphas,maxNumAtomAlphas))

   ! Allocate space for the intermediate steps of the integration.
   allocate (pairXBasisFn2dHam (16,maxNumAtomAlphas,maxNumStates,spin,3))
   allocate (pairXBasisFn12dHam (maxNumStates,maxNumStates,spin,3))
#ifndef GAMMA
   allocate (currentPair      (maxNumStates,maxNumStates,numKPoints,spin,3))
#else
   allocate (currentPairGamma (maxNumStates,maxNumStates,spin,3))
#endif

   ! Initialize hamiltonian matrices.
   pairXBasisFn12dHam(:,:,:,:) = 0.0_double

   ! Initialize the main force matrices.
#ifndef GAMMA
   valeValeF(:,:,:,:,:) = 0.0_double
   coreCoreF(:,:,:,:,:) = 0.0_double
   coreValeF(:,:,:,:,:) = 0.0_double
#else
   valeValeFGamma(:,:,:,:) = 0.0_double
   coreCoreFGamma(:,:,:,:) = 0.0_double
   coreValeFGamma(:,:,:,:) = 0.0_double
#endif


   ! Start the outer atom loop.
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

         ! These should probably be computed ahead of time and then used
         !   above and below.
         do k = 1,currentNumAlphas(2)
            alphaCenter(:currentNumAlphas(1),k) = &
                  &  currentAlphas(:currentNumAlphas(1),1) / &
                  & (currentAlphas(:currentNumAlphas(1),1) + &
                  &  currentAlphas(k,2))
            alphaDist(:currentNumAlphas(1),k) = logBasisFnThresh / &
                  &  currentAlphas(k,2) / alphaCenter(:currentNumAlphas(1),k)
         enddo

         ! Initialize the matrix for this atom pair that will record the
         !   integrals with bloch vector (kpoint) phase factors.
#ifndef GAMMA
         currentPair(:,:,:,:,:) = 0.0_double
#else
         currentPairGamma(:,:,:,:) = 0.0_double
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

            ! Obtain the seperation vector between atom 1 and the shifted
            !   position of atom 2.
            shiftedAtomSiteSep = sum ((currentPosition(:,1) - &
                  & shiftedAtomPos(:))**2)

            ! Determine if this shifted atom position puts it outside of the
            !   above determined negligability limit for this atom pair.
            if (shiftedAtomSiteSep > currentNegligLimit) cycle


            ! Initialize the matrices to hold the product of the integrals
            !    times the atom2 wave functions.
            pairXBasisFn2dHam(:,:currentNumAlphas(1),&
                  &:currentNumTotalStates(2),:spin,:) = 0.0_double

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
                        & :currentlmAlphaIndex(alphaIndex(2),2),:spin,:) = &
                        & 0.0_double


                  ! Initialize a counter to track which potential
                  !    coefficient we are currently calculating on.
                  potCoeffIndex = 0

                  ! Initiate a loop over each potential site for the nuclear
                  !   potentials.
                  do m = 1, numPotSites

                     ! Skip equivalent types now.  They will be accounted
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
                              shiftedCenterOriginSep = &
                                    & sum((centerOriginVect(:) - &
                                    & cellDimsReal(:,p))**2)

                              ! Check if shifted seperation between the
                              !   center and the origin extends past the
                              !   negligability limit.  If so, cycle to the
                              !   next cell.
                              if (shiftedCenterOriginSep > &
                                    & threeAlphaDist) cycle

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

                                 call delectron3CIntgCB (&
                                    & currentAlphas(alphaIndex(1),1),&
                                    & currentAlphas(alphaIndex(2),2),&
                                    & currentPotAlpha, currentPosition(:,1),&
                                    & shiftedAtomPos(:), shiftedPotPos(:),&
                                    & l1l2Switch, oneAlphaSet)

                                 ! Accumulate results returned for alpha set.
                                 do q = 1, spin
                                 potAtomOverlap(:currentlmAlphaIndex &
                                    & (alphaIndex(1),1),:currentlmAlphaIndex&
                                    & (alphaIndex(2),2),q,:) = &
                                    & potAtomOverlap(:currentlmAlphaIndex &
                                    & (alphaIndex(1),1),:currentlmAlphaIndex&
                                    & (alphaIndex(2),2),q,:) - &
                                    & oneAlphaSet(:currentlmAlphaIndex &
                                    & (alphaIndex(1),1),:currentlmAlphaIndex&
                                    & (alphaIndex(2),2),:) * &
                                    & currentPotCoeff(q)
                                 enddo
                              else
                                ! Calculate the opcode to do the correct set
                                ! of integrals for the current alpha pair
                                l1l2Switch = ishft(1,&
                                  &(powerOfTwo(currentlmAlphaIndex(&
                                  &   alphaIndex(1),1)))) &
                                  &+ ishft(16,&
                                  &(powerOfTwo(currentlmAlphaIndex(&
                                  &   alphaIndex(2),2))))
                                
                                call dnuclear3CIntgCB (&
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
                                       & (alphaIndex(2),2),q,:) = &
                                       & potAtomOverlap(:currentlmAlphaIndex&
                                       & (alphaIndex(1),1),&
                                       & :currentlmAlphaIndex &
                                       & (alphaIndex(2),2),q,:) + &
                                       & oneAlphaSet(:currentlmAlphaIndex &
                                       & (alphaIndex(1),1),&
                                       & :currentlmAlphaIndex &
                                       & (alphaIndex(2),2),:) * zFactor
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

                  ! Determine the kinetic energy contribution.
                  call dkinetic2CIntg (&
                        & currentAlphas(alphaIndex(1),1),&
                        & currentAlphas(alphaIndex(2),2),&
                        & currentPosition(:,1), shiftedAtomPos(:),&
                        & l1l2Switch, oneAlphaSet)

                  ! Accumulate the contribution from this alpha pair
                  do m = 1, spin
                     potAtomOverlap(:currentlmAlphaIndex &
                           & (alphaIndex(1),1),:currentlmAlphaIndex &
                           & (alphaIndex(2),2),m,:) = &
                           & potAtomOverlap(:currentlmAlphaIndex &
                           & (alphaIndex(1),1),:currentlmAlphaIndex &
                           & (alphaIndex(2),2),m,:) - &
                           & oneAlphaSet(:currentlmAlphaIndex &
                           & (alphaIndex(1),1),:currentlmAlphaIndex &
                           & (alphaIndex(2),2),:)
                  enddo

                  ! Compute the mass velocity integral if needed for the
                  !   scalar relativistic calculation.
                  ! FIX: Do derivative of the relativistic mass velocity
                  !   term to include this effect on forces. (Probably
                  !   a small effect so we ignore now.)
                  !if (rel == 1) then
                  !   ! Calculate the opcode to do the correct set of
                  !   !   integrals for the current alpha pair.
                  !   l1l2Switch = ishft(1,&
                  !   &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
                  !   &+ ishft(16,&
                  !   &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))

                  !   ! Compute the integral.
                  !   call massVel2CIntg (currentAlphas(alphaIndex(1),1),&
                  !         & currentAlphas(alphaIndex(2),2),&
                  !         & currentPosition(:,1), shiftedAtomPos(:),&
                  !         & l1l2Switch, oneAlphaSet)
                  !      
                  !   ! Accumulate the contribution from this alpha pair.
                  !   !   Note: a minus sign is used in the accumulation
                  !   !   for the mass velocity integral because it has an
                  !   !   overall negative contribution to the Hamiltonian.
                  !   do m = 1, spin
                  !      potAtomOverlap(:currentlmAlphaIndex &
                  !            & (alphaIndex(1),1),:currentlmAlphaIndex &
                  !            & (alphaIndex(2),2),m) = &
                  !            & potAtomOverlap(:currentlmAlphaIndex &
                  !            & (alphaIndex(1),1),:currentlmAlphaIndex &
                  !            & (alphaIndex(2),2),m) - &
                  !            & oneAlphaSet(:currentlmAlphaIndex &
                  !            & (alphaIndex(1),1),:currentlmAlphaIndex &
                  !            & (alphaIndex(2),2))
                  !   enddo
                  !endif

                  ! FIX: It seems unnecessary to recompute the l1l2Switch
                  !   for OL and MV matrices after it is computed
                  !   once for the KE.
                  ! Calculate the opcode to do the correct set of integrals
                  ! for the current alpha pair.
                  l1l2Switch = ishft(1,&
                     &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
                     &+ ishft(16,&
                     &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))
                
                  ! Include the pulay correction (by doing the momentum
                  !   integral calculation).
                  call momentum2CIntg (&
                        & currentAlphas(alphaIndex(1),1),&
                        & currentAlphas(alphaIndex(2),2), &
                        & currentPosition(:,1), shiftedAtomPos(:),&
                        & l1l2Switch, oneAlphaSet)

                  ! Accumulate the contribution from this alpha set
                  do m = 1, spin
                     potAtomOverlap(:currentlmAlphaIndex &
                           & (alphaIndex(1),1),:currentlmAlphaIndex &
                           & (alphaIndex(2),2),m,:) = &
                           & potAtomOverlap(:currentlmAlphaIndex &
                           & (alphaIndex(1),1),:currentlmAlphaIndex &
                           & (alphaIndex(2),2),m,:) - 2 * totalEnergy * &
                           & oneAlphaSet(:currentlmAlphaIndex &
                           & (alphaIndex(1),1),:currentlmAlphaIndex &
                           & (alphaIndex(2),2),:)
                  enddo

                  ! Collect the results of the overlap of the current alpha
                  !   times the wave functions for atom 2.

                  ! Potential overlaps first (if any were found).
                  do q = 1, spin
                     do m = 1, currentNumTotalStates(2)
                        pairXBasisFn2dHam(:currentlmAlphaIndex( &
                              & alphaIndex(1),1),alphaIndex(1),m,q,:) = &
                              & pairXBasisFn2dHam(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),alphaIndex(1),m,q,:) + &
                              & potAtomOverlap(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),currentlmIndex(m,2),q,:) * &
                              & currentBasisFns(alphaIndex(2),m,2)
                     enddo
                  enddo

                  ! Update the maximum alpha used from atom 1.
                  maxAlpha1Used = max(alphaIndex(1),maxAlpha1Used)

                  ! Switch mode from the initial state of no searching to the
                  !   state of searching along alphas from atom 1.
                  if (currentMode == 0) then
                     currentMode = 1
                  endif
               enddo
            enddo  ! min number of alphas between the two atoms l

            if (maxAlpha1Used > 0) then

               ! At this point all the alpha loops are complete and we can form
               !   a product with the atom 1 wave function to give the overlap
               !   integral in a complete wavefunction representation.
               do q = 1, spin
                  do r = 1, 3
                     do l = 1, currentNumTotalStates(2)
                        do m = 1, currentNumTotalStates(1)
                           pairXBasisFn12dHam(m,l,q,r) = &
                                 & sum(currentBasisFns(:maxAlpha1Used,m,1) * &
                                 & pairXBasisFn2dHam(currentlmIndex(m,1),&
                                 & :maxAlpha1Used,l,q,r))
                        enddo
                     enddo
                  enddo
               enddo

#ifndef GAMMA
               do q = 1, spin
                  do r = 1, 3
                     call applyPhaseFactors (currentPair(:,:,:,q,r),&
                           & pairXBasisFn12dHam(1:currentNumTotalStates(1),&
                           & 1:currentNumTotalStates(2),q,r),&
                           & currentNumTotalStates(1),&
                           & currentNumTotalStates(2),k,0,0)
                  enddo
               enddo
#else
               do q = 1, spin
                  do r = 1, 3
                     call applyPhaseFactorsGamma (currentPairGamma(:,:,q,r),&
                           & pairXBasisFn12dHam(1:currentNumTotalStates(1),&
                           & 1:currentNumTotalStates(2),q,r),&
                           & currentNumTotalStates(1),&
                           & currentNumTotalStates(2),0)
                  enddo
               enddo
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
         do q = 1, spin
            do r = 1, 3
               call kPointLatticeOriginShift (currentNumTotalStates,&
                     & currentPair(:,:,:,q,r),latticeVector,numKPoints,0)
               call saveCurrentPair(i,j,numKPoints,currentPair(:,:,:,q,r),&
                     & valeValeF(:,:,:,q,r),coreValeF(:,:,:,q,r),&
                     & coreCoreF(:,:,:,q,r),1)
            enddo
         enddo
#else
         do q = 1, spin
            do r = 1, 3
               call saveCurrentPairGamma(i,j,currentPairGamma(:,:,q,r),&
                     & valeValeFGamma(:,:,q,r),coreValeFGamma(:,:,q,r),&
                     & coreCoreFGamma(:,:,q,r))
            enddo
         enddo
#endif

      enddo ! (Atom loop #2)
   enddo    ! (Atom loop #1)

   ! Print the results
#ifndef GAMMA
   do i = 1,3
      do j = 1,spin
         write (97+j,*) "i(x,y,z)=",i
         write (97+j,*) "j(spin)=",j
         do k = 1, numKPoints
            write (97+j,*) "k(kp)=",k
            do l = 1, valeDim
               do m = 1, valeDim
                  write (97+j,advance="no",fmt="(e13.5,a1)") real(valeValeF(m,l,k,j,i)), " "
               enddo
               write (97+j,*)
            !   write (97+j,*) real(sum(valeValeF(:,l,k,j,i)))
            enddo
            write (97+j,*)
            do l = 1, valeDim
               do m = 1, valeDim
                  write (97+j,advance="no",fmt="(e13.5,a1)") aimag(valeValeF(m,l,k,j,i)), " "
               enddo
               write (97+j,*)
            !   write (97+j,*) real(sum(valeValeF(:,l,k,j,i)))
            enddo
            write (97+j,*)
         enddo
      enddo
   enddo
#else
   do i = 1,3
      do j = 1,spin
         write (97+j,*) "i(x,y,z)=",i
         write (97+j,*) "j(spin)=",j
         do l = 1, valeDim
            write (97+j,*) valeValeFGamma(:,l,q,r)
         enddo
      enddo
   enddo
#endif


   ! Deallocate all arrays and matrices before exiting.
   deallocate (currentBasisFns)
   deallocate (currentAlphas)
   deallocate (currentlmAlphaIndex)
   deallocate (currentlmIndex)
   deallocate (alphaDist)   ! Can be saved from before
   deallocate (alphaCenter) ! Can be saved for later.

   deallocate (pairXBasisFn2dHam)
   deallocate (pairXBasisFn12dHam)

!#ifndef GAMMA
!   deallocate (valeValeF)
!   deallocate (coreCoreF)
!   deallocate (coreValeF)
!#else
!   deallocate (valeValeFGamma)
!   deallocate (coreCoreFGamma)
!   deallocate (coreValeFGamma)
!#endif

   ! Log the date and time we start.
   call timeStampEnd(20)

end subroutine computeForceIntg


#ifndef GAMMA
subroutine computeForce(valeValeRho, kPoint, currSpin, currAxis)

   use O_Kinds
   use O_AtomicSites, only : numAtomSites, atomSites, valeDim
   use O_AtomicTypes, only : atomTypes
   use O_PotTypes, only: potTypes
   use O_Potential, only: spin

   implicit none

   ! Define passed parameters.
   complex (kind=double), dimension (valeDim,valeDim,spin), &
         & intent(inout) :: valeValeRho
   integer, intent(in) :: kPoint, currSpin, currAxis

   integer :: i,j,k,l,m
   integer :: counter1, counter2
   real (kind=double), allocatable, dimension(:,:) :: atomSum
   real (kind=double) :: theta
   real (kind=double) :: phi
   real (kind=double) :: sepDist
   real (kind=double) :: nucForceMag
   real (kind=double), dimension(3) :: dxyz
   real (kind=double), dimension(3) :: nucForce
   real (kind=double), dimension(3) :: totalNucForce

   allocate(atomSum(numAtomSites,numAtomSites))

   valeValeRho(:,:,currSpin) = valeValeRho(:,:,currSpin) + &
         & transpose(valeValeRho(:,:,currSpin))
   do l = 1, valeDim
      valeValeRho(l,l,currSpin) = valeValeRho(l,l,currSpin) / 2.0_double
   enddo
   valeValeF(:,:,kPoint,currSpin,currAxis) = &
         & valeValeF(:,:,kpoint,currSpin,currAxis) * &
         & valeValeRho(:,:,currSpin)

   do i = 1, valeDim
      do j = 1, valeDim
         write (6,advance="no",fmt="(e13.5,a1)") &
               & real(valeValeF(j,i,kPoint,currSpin,currAxis)), " "
      enddo
      write (6,*)
   enddo
   write (6,*)
   do i = 1, valeDim
      do j = 1, valeDim
         write (6,advance="no",fmt="(e13.5,a1)") &
               & aimag(valeValeF(j,i,kPoint,currSpin,currAxis)), " "
      enddo
      write (6,*)
   enddo
   write (6,*)

   atomSum(:,:) = 0.0_double
   counter1 = 0
   do i = 1, numAtomSites
      do j = 1, atomTypes(atomSites(i)%atomTypeAssn)%numValeStates
         counter1 = counter1 + 1
         counter2 = 0
         do k = 1, numAtomSites
            do l = 1, atomTypes(atomSites(k)%atomTypeAssn)%numValeStates
               counter2 = counter2 + 1
               atomSum(i,k) = atomSum(i,k) + &
                     & valeValeF(counter1,counter2,kPoint,currSpin,currAxis)
            enddo
         enddo
      enddo
   enddo

   do i = 1, numAtomSites
      do j = i+1, numAtomSites
         atomSum(j,i) = atomSum(i,j) * (-1.0_double)
      enddo
   enddo

   do i = 1, numAtomSites
      totalNucForce(:) = 0.0_double ! Force on atom i
      do j = 1, numAtomSites
         ! Magnitude of force on i from j.
         nucForce(:) = 0.0_double
         if (i /= j) then
            dxyz(:) = atomSites(i)%cartPos(:) - atomSites(j)%cartPos(:)
!write (6,*) "dxyz=", dxyz(:)
            sepDist = sqrt(sum((dxyz)**2))
            theta = acos(dxyz(3) / sepDist)
            phi = sign(1.0_double,dxyz(2)) * acos(dxyz(1) / &
                  & sqrt(dxyz(1)**2 + dxyz(2)**2))
!write (6,*) "theta=", theta
!write (6,*) "phi=", phi
            nucForceMag = potTypes(atomSites(i)%atomTypeAssn)%nucCharge * &
                  potTypes(atomSites(j)%atomTypeAssn)%nucCharge / sepDist**2
            nucForce(1) = nucForceMag * sin(theta)*cos(phi)
            nucForce(2) = nucForceMag * sin(theta)*sin(phi)
            nucForce(3) = nucForceMag * cos(theta)
            totalNucForce(:) = totalNucForce(:) + nucForce(:)
         endif
         write(6,fmt="(2i4,4e13.5)") i,j,atomSum(i,j),nucForce(:)
      enddo
      write (6,fmt="(i4,4e13.5)") i,sum(atomSum(i,:)),totalNucForce(:)
   enddo

   deallocate(atomSum)

end subroutine computeForce

#else

subroutine computeForceGamma(valeValeRhoGamma, currSpin, currAxis)

   use O_Kinds
   use O_AtomicSites, only : numAtomSites, atomSites, valeDim
   use O_AtomicTypes, only : atomTypes
   use O_PotTypes, only: potTypes
   use O_Potential, only: spin

   implicit none

   ! Define passed parameters.
   real (kind=double), dimension (valeDim,valeDim,spin), &
         & intent(inout) :: valeValeRhoGamma
   integer, intent(in) :: currSpin, currAxis

   ! Define local variables.
   integer :: i,j,k,l,m
   integer :: counter1, counter2
   real (kind=double), allocatable, dimension(:,:) :: atomSum
   real (kind=double) :: theta
   real (kind=double) :: phi
   real (kind=double) :: sepDist
   real (kind=double) :: nucForceMag
   real (kind=double), dimension(3) :: dxyz
   real (kind=double), dimension(3) :: nucForce
   real (kind=double), dimension(3) :: totalNucForce

   allocate(atomSum(numAtomSites,numAtomSites))

   valeValeRhoGamma(:,:,currSpin) = valeValeRhoGamma(:,:,currSpin) &
         & + transpose(valeValeRhoGamma(:,:,currSpin))
   do l = 1, valeDim
      valeValeRhoGamma(l,l,currSpin) = valeValeRhoGamma(l,l,currSpin) &
            & / 2.0_double
   enddo
   valeValeFGamma(:,:,currSpin,currAxis) = &
         & valeValeFGamma(:,:,currSpin,currAxis) * &
         & valeValeRhoGamma(:,:,currSpin)

   atomSum(:,:) = 0.0_double
   counter1 = 0
   do i = 1, numAtomSites
      do j = 1, atomTypes(atomSites(i)%atomTypeAssn)%numValeStates
         counter1 = counter1 + 1
         counter2 = 0
         do k = 1, numAtomSites
            do l = 1, atomTypes(atomSites(k)%atomTypeAssn)%numValeStates
               counter2 = counter2 + 1
               atomSum(i,k) = atomSum(i,k) + &
                     & valeValeFGamma(counter1,counter2,currSpin,currAxis)
            enddo
         enddo
      enddo
   enddo

   do i = 1, numAtomSites
      do j = i+1, numAtomSites
         atomSum(j,i) = atomSum(i,j) * -1.0_double
      enddo
   enddo

   do i = 1, numAtomSites
      totalNucForce(:) = 0.0_double ! Force on atom i
      do j = 1, numAtomSites
         ! Magnitude of force on i from j.
         nucForce(:) = 0.0_double
         if (i /= j) then
            dxyz(:) = atomSites(i)%cartPos(:) - atomSites(j)%cartPos(:)
!write (6,*) "dxyz=", dxyz(:)
            sepDist = sqrt(sum((dxyz)**2))
            theta = acos(dxyz(3) / sepDist)
            phi = sign(1.0_double,dxyz(2)) * acos(dxyz(1) / sqrt(dxyz(1)**2 + dxyz(2)**2))
!write (6,*) "theta=", theta
!write (6,*) "phi=", phi
            nucForceMag = potTypes(atomSites(i)%atomTypeAssn)%nucCharge * &
                  potTypes(atomSites(j)%atomTypeAssn)%nucCharge / sepDist**2
            nucForce(1) = nucForceMag * sin(theta)*cos(phi)
            nucForce(2) = nucForceMag * sin(theta)*sin(phi)
            nucForce(3) = nucForceMag * cos(theta)
            totalNucForce(:) = totalNucForce(:) + nucForce(:)
         endif
         write(6,fmt="(2i4,4e13.5)") i,j,atomSum(i,j),nucForce(:)
      enddo
      write (6,fmt="(i4,4e13.5)") i,sum(atomSum(i,:)),totalNucForce(:)
   enddo

   deallocate(atomSum)

end subroutine computeForceGamma
#endif

end module O_Force
