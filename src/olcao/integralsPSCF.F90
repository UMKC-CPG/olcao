module O_IntegralsPSCF

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

   ! Structures that will hold the results from many iterations of the atom-
   !   lattice loops accumulations.
   real (kind=double), allocatable, dimension (:,:)   :: accumulatedOL
   real (kind=double), allocatable, dimension (:,:)   :: accumulatedMomX
   real (kind=double), allocatable, dimension (:,:)   :: accumulatedMomY
   real (kind=double), allocatable, dimension (:,:)   :: accumulatedMomZ
   real (kind=double), allocatable, dimension (:,:,:) :: accumulatedHam

   ! Structure that will track the indices of the loops for the above data.
   integer, allocatable, dimension (:,:) :: loopIndices
   integer, allocatable, dimension (:,:) :: currentLoopIndices

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
   use HDF5

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define passed parameters.
   integer :: doINTG  ! (1) = do overlap & hamiltonian integrals; (0) = do not.
   integer, intent(in) :: doMOME

   ! Define local variables for logging, loop control, and file writing.
   integer :: i,j,k,l,m,n,o,p,q ! Loop index variables
   integer :: totalDataCounter ! Track the number of iterations that have been
                               !   done before a chunk must be saved to disk.
   integer(hsize_t) :: cumulDataSize ! Cumulatively track the size of the data
                               !   computed through the loops.
   integer :: exactFit ! Flag to mark that the last recorded chunk was an
                       !   exact fit and that no more small leftover stubs of
                       !   data will need to be saved.

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
   real (kind=double), dimension (16,16)   :: oneAlphaSet ! The potential
         ! overlap from one alpha pair.
   real (kind=double), dimension (16,16,3) :: oneAlphaSetMom ! The overlap
         ! contribution to the momentum matrix from one alpha pair.
   real (kind=double), dimension (16,16,2) :: potAtomOverlap ! The accumulated
         ! values from all the oneAlphaSet calculations.
   real (kind=double), allocatable, dimension (:,:,:) :: pairXWaveFn2OL ! The
         ! overlap times the wave function from atom 2.  This is accumulated
         ! with each iteration of the alpha loop.
   real (kind=double), allocatable, dimension (:,:,:,:) :: pairXWaveFn2Ham ! The
         ! potential overlap times the wave function from atom 2.  This is
         ! accumulated with each iteration of the alpha loop.
   real (kind=double), allocatable, dimension (:,:,:) :: pairXWaveFn2MomX ! The
         ! overlap times the wave function from atom 2.  This is accumulated
         ! with each iteration of the alpha loop.
   real (kind=double), allocatable, dimension (:,:,:) :: pairXWaveFn2MomY ! The
         ! overlap times the wave function from atom 2.  This is accumulated
         ! with each iteration of the alpha loop.
   real (kind=double), allocatable, dimension (:,:,:) :: pairXWaveFn2MomZ ! The
         ! overlap times the wave function from atom 2.  This is accumulated
         ! with each iteration of the alpha loop.


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
   real (kind=double), allocatable, dimension (:,:,:) :: pairXWaveFn12Ham
   real (kind=double), allocatable, dimension (:,:)   :: pairXWaveFn12OL
   real (kind=double), allocatable, dimension (:,:)   :: pairXWaveFn12MomX
   real (kind=double), allocatable, dimension (:,:)   :: pairXWaveFn12MomY
   real (kind=double), allocatable, dimension (:,:)   :: pairXWaveFn12MomZ



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
      allocate (pairXWaveFn2OL     (16,maxNumAtomAlphas,maxNumStates))
      allocate (pairXWaveFn2Ham    (16,maxNumAtomAlphas,maxNumStates,spin))
      allocate (pairXWaveFn12OL    (maxNumStates,maxNumStates))
      allocate (pairXWaveFn12Ham   (maxNumStates,maxNumStates,spin))

      ! Initialize overlap and hamiltonian matrices.
      pairXWaveFn12OL(:,:)    = 0.0_double
      pairXWaveFn12Ham(:,:,:) = 0.0_double
   endif

   if (doMOME == 1) then
      allocate (pairXWaveFn2MomX   (16,maxNumAtomAlphas,maxNumStates))
      allocate (pairXWaveFn2MomY   (16,maxNumAtomAlphas,maxNumStates))
      allocate (pairXWaveFn2MomZ   (16,maxNumAtomAlphas,maxNumStates))
      allocate (pairXWaveFn12MomX  (maxNumStates,maxNumStates))
      allocate (pairXWaveFn12MomY  (maxNumStates,maxNumStates))
      allocate (pairXWaveFn12MomZ  (maxNumStates,maxNumStates))

      ! Initialize momentum matrix element matrices.
      pairXWaveFn12MomX(:,:) = 0.0_double
      pairXWaveFn12MomY(:,:) = 0.0_double
      pairXWaveFn12MomZ(:,:) = 0.0_double
   endif


   ! Initialize the counter for the number of data sub-chunk contributions have
   !   occured.
   totalDataCounter = 0

   ! Initialize the cumulative data size that holds the amount of data in the
   !   current chunk.
   cumulDataSize = 0

   ! Initialize other variables.
   exactFit = 0

   do i = 1, numAtomSites

      ! Initialize the HDF5 dataset and dataspace for this atom.
      call initAtomHDF (i,doINTG,doMOME)

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

         ! Form an alpha overlap matrix that is later used to determine whether
         !   a particular alpha pair needs to be considered in the overlap
         !   matrix or if the alpha pair do not have sufficient overlap.
         !   NOTE that it may be possible to improve this even more by asking
         !   some simple question that will determine quickly if the alphas for
         !   these two atoms are the same.  In that case the alphasOverlap
         !   matrix is symmetric and can be calculated in about half the time.
         !   That is not checked for here, but in the future it could be.
         ! At first I thought that this would be a good idea, but it turns out
         !   to not be efficient for large cell systems where each atom will
         !   have an overlap with another atom only once and even then, not
         !   every alpha pair will be used so this might actually perform worse.
         !   It may be possible to salvage it in one of two ways.  Save this
         !   result for the cases of later overlap integral calculations.  This
         !   would require a lot of extra memory usage though.  OR, have a
         !   runtime switch that will utilize this method for the case of small
         !   systems, and the other method for the case of large systems.  That
         !   is probably not worth it though since small systems will be rather
         !   fast anyway.  Well, we will see.
         ! Another possibility (that can be applied to other similar algorithms
         !   in other subroutines) is to calculate these values based on the
         !   element as opposed to the current method of calculating them based
         !   on every atom pair iteration.  In this way the values would only
         !   have to be calculated once for each element pair.
         do k = 1,currentNumAlphas(2)
            alphaCenter(:currentNumAlphas(1),k) = &
                  &  currentAlphas(:currentNumAlphas(1),1) / &
                  & (currentAlphas(:currentNumAlphas(1),1) + &
                  &  currentAlphas(k,2))
            alphaDist(:currentNumAlphas(1),k) = logBasisFnThresh / &
                  &  currentAlphas(k,2) / alphaCenter(:currentNumAlphas(1),k)
         enddo


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
               pairXWaveFn2OL(:,:currentNumAlphas(1),&
                     &:currentNumTotalStates(2)) = 0.0_double
               pairXWaveFn2Ham(:,:currentNumAlphas(1),&
                     &:currentNumTotalStates(2),:spin) = 0.0_double
            endif
            if (doMOME == 1) then
               pairXWaveFn2MomX(:,:currentNumAlphas(1),&
                     &:currentNumTotalStates(2)) = 0.0_double
               pairXWaveFn2MomY(:,:currentNumAlphas(1),&
                     &:currentNumTotalStates(2)) = 0.0_double
               pairXWaveFn2MomZ(:,:currentNumAlphas(1),&
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

                                    call threeCentInteg (currentAlphas(alphaIndex(1),1),&
                                       & currentAlphas(alphaIndex(2),2),&
                                       & currentPotAlpha, currentPosition(:,1),&
                                       & shiftedAtomPos(:), shiftedPotPos(:),&
                                       & l1l2Switch, oneAlphaSet)


                                    !call threeCentInteg (oneAlphaSet,&
                                    !   & currentlmAlphaIndex (alphaIndex(1),1),&
                                    !   & currentlmAlphaIndex (alphaIndex(2),2),&
                                    !   & currentAlphas(alphaIndex(1),1),&
                                    !   & currentAlphas(alphaIndex(2),2),&
                                    !   & currentPotAlpha,currentPosition(:,1),&
                                    !   & shiftedAtomPos(:),shiftedPotPos(:))

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
                                   l1l2Switch = ishft(1,&
                                     &(powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
                                     &+ ishft(16,&
                                     &(powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))
                                   
                                   call nucPotInteg (currentAlphas(alphaIndex(1),1),&
                                       & currentAlphas(alphaIndex(2),2),&
                                       & currentPotAlpha,currentPosition(:,1),&
                                       & shiftedAtomPos(:),shiftedPotPos(:),&
                                       & l1l2Switch,oneAlphaSet)



                                    !call nucPotInteg (oneAlphaSet,&
                                    !   & currentlmAlphaIndex (alphaIndex(1),1),&
                                    !   & currentlmAlphaIndex (alphaIndex(2),2),&
                                    !   & currentAlphas(alphaIndex(1),1),&
                                    !   & currentAlphas(alphaIndex(2),2),&
                                    !   & currentPotAlpha,currentPosition(:,1),&
                                    !   & shiftedAtomPos(:),shiftedPotPos(:))

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

                     ! Determine the kinetic energy contribution
                     !call KEInteg (oneAlphaSet,&
                     !      & currentlmAlphaIndex (alphaIndex(1),1),&
                     !      & currentlmAlphaIndex (alphaIndex(2),2),&
                     !      & currentAlphas(alphaIndex(1),1),&
                     !      & currentAlphas(alphaIndex(2),2),&
                     !      & currentPosition(:,1),shiftedAtomPos(:))

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
                
                !print*,l1l2Switch
                ! We can proceed with the next step of the calculation.
                ! This is the actual integral.
                call overlapInteg (currentAlphas(alphaIndex(1),1),&
                  & currentAlphas(alphaIndex(2),2), &
                  & currentPosition(:,1), shiftedAtomPos(:),&
                  & l1l2Switch, oneAlphaSet)


                     ! Compute the atomic overlap for this atom pair.
                     !call overlapInteg (oneAlphaSet,&
                     !      & currentlmAlphaIndex (alphaIndex(1),1),&
                     !      & currentlmAlphaIndex (alphaIndex(2),2),&
                     !      & currentAlphas(alphaIndex(1),1),&
                     !      & currentAlphas(alphaIndex(2),2),&
                     !      & currentPosition(:,1),shiftedAtomPos(:))
                  endif


                  ! Compute the momentum matrix values if requested.
                  if (doMOME == 1) then
                     ! Calculate the opcode to do the correct set
                     ! of integrals for the current alpha pair
                     l1l2Switch = ishft(1,&
                        & (powerOfTwo(currentlmAlphaIndex(alphaIndex(1),1))))&
                        & + ishft(16,&
                        & (powerOfTwo(currentlmAlphaIndex(alphaIndex(2),2))))

                     call MOMF (currentAlphas(alphaIndex(1),1),&
                           & currentAlphas(alphaIndex(2),2),&
                           & currentPosition(:,1), shiftedAtomPos(:),&
                           & l1l2Switch, oneAlphaSetMom)
                     !call MOMF (oneAlphaSetMom,&
                     !      & currentlmAlphaIndex (alphaIndex(1),1),&
                     !      & currentlmAlphaIndex (alphaIndex(2),2),&
                     !      & currentAlphas(alphaIndex(1),1),&
                     !      & currentAlphas(alphaIndex(2),2),&
                     !      & currentPosition(:,1),shiftedAtomPos(:))
                  endif

                  ! Collect the results of the overlap of the current alpha
                  !   times the wave functions for atom 2.

                  if (doINTG == 1) then
                     ! Potential overlaps first (if any were found).
                     do q = 1, spin
                        do m = 1, currentNumTotalStates(2)
                           pairXWaveFn2Ham(:currentlmAlphaIndex( &
                                 & alphaIndex(1),1),alphaIndex(1),m,q) = &
                                 & pairXWaveFn2Ham(:currentlmAlphaIndex &
                                 & (alphaIndex(1),1),alphaIndex(1),m,q) + &
                                 & potAtomOverlap(:currentlmAlphaIndex &
                                 & (alphaIndex(1),1),currentlmIndex(m,2),q) * &
                                 & currentBasisFns(alphaIndex(2),m,2)
                        enddo
                     enddo
                     ! Atom pair overlaps second.
                     do m = 1, currentNumTotalStates(2)
                        pairXWaveFn2OL(:currentlmAlphaIndex(alphaIndex(1),1),&
                              & alphaIndex(1),m) = &
                              & pairXWaveFn2OL(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),alphaIndex(1),m) + &
                              & oneAlphaSet(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),currentlmIndex(m,2)) * &
                              & currentBasisFns(alphaIndex(2),m,2)
                     enddo
                  endif

                  ! Momentum matrix last.  (If requested)
                  if (doMOME == 1) then
                     do m = 1, currentNumTotalStates(2)
                        pairXWaveFn2MomX(:currentlmAlphaIndex(alphaIndex(1),&
                              & 1), alphaIndex(1),m) = &
                              & pairXWaveFn2MomX(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),alphaIndex(1),m) + &
                              & oneAlphaSetMom(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),currentlmIndex(m,2),1) * &
                              & currentBasisFns(alphaIndex(2),m,2)
                        pairXWaveFn2MomY(:currentlmAlphaIndex(alphaIndex(1),&
                              & 1), alphaIndex(1),m) = &
                              & pairXWaveFn2MomY(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),alphaIndex(1),m) + &
                              & oneAlphaSetMom(:currentlmAlphaIndex &
                              & (alphaIndex(1),1),currentlmIndex(m,2),2) * &
                              & currentBasisFns(alphaIndex(2),m,2)
                        pairXWaveFn2MomZ(:currentlmAlphaIndex(alphaIndex(1),&
                              & 1), alphaIndex(1),m) = &
                              & pairXWaveFn2MomZ(:currentlmAlphaIndex &
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
            !   product with the atom 1 wave function to give the overlap
            !   integral in a complete wavefunction representation.
            if (doINTG == 1) then
               do l = 1, currentNumTotalStates(2)
                  do m = 1, currentNumTotalStates(1)
                     pairXWaveFn12OL(m,l) = &
                           & sum(currentBasisFns(:maxAlpha1Used,m,1) * &
                           & pairXWaveFn2OL(currentlmIndex(m,1),&
                           & :maxAlpha1Used,l))
                  enddo
               enddo
               do q = 1, spin
                  do l = 1, currentNumTotalStates(2)
                     do m = 1, currentNumTotalStates(1)
                        pairXWaveFn12Ham(m,l,q) = &
                              & sum(currentBasisFns(:maxAlpha1Used,m,1) * &
                              & pairXWaveFn2Ham(currentlmIndex(m,1),&
                              & :maxAlpha1Used,l,q))
                     enddo
                  enddo
               enddo
            endif

            ! The momentum is summed against the wave function 1 if needed.
            if (doMOME == 1) then
               do l = 1, currentNumTotalStates(2)
                  do m = 1, currentNumTotalStates(1)
                  pairXWaveFn12MomX(m,l) = &
                        & sum(currentBasisFns(:maxAlpha1Used,m,1) * &
                        & pairXWaveFn2MomX(currentlmIndex(m,1), &
                        & :maxAlpha1Used,l))
                  pairXWaveFn12MomY(m,l) = &
                        & sum(currentBasisFns(:maxAlpha1Used,m,1) * &
                        & pairXWaveFn2MomY(currentlmIndex(m,1), &
                        & :maxAlpha1Used,l))
                  pairXWaveFn12MomZ(m,l) = &
                        & sum(currentBasisFns(:maxAlpha1Used,m,1) * &
                        & pairXWaveFn2MomZ(currentlmIndex(m,1), &
                        & :maxAlpha1Used,l))
                  enddo
               enddo
            endif

            ! At this point the current lattice sum for the current atom pair
            !   is complete and the results can be either recorded to disk or
            !   recorded to the memory buffer.

            ! Increment the totalDataCounter to record index information.
            totalDataCounter = totalDataCounter + 1

            ! Copy the current pairXWaveFn12 matrices to the appropriate
            !   accumulation buffers.
            if (doINTG == 1) then
               accumulatedOL(:currentNumTotalStates(1),cumulDataSize + 1:&
                     & cumulDataSize + currentNumTotalStates(2)) = &
                     & pairXWaveFn12OL(:currentNumTotalStates(1),:&
                     & currentNumTotalStates(2))
               accumulatedHam(:currentNumTotalStates(1),cumulDataSize + 1:&
                     & cumulDataSize + currentNumTotalStates(2),:spin) = &
                     & pairXWaveFn12Ham(:currentNumTotalStates(1),:&
                     & currentNumTotalStates(2),:spin)
            endif
            if (doMOME == 1) then
               accumulatedMomX(:currentNumTotalStates(1),cumulDataSize + 1:&
                     & cumulDataSize + currentNumTotalStates(2)) = &
                     & pairXWaveFn12MomX(:currentNumTotalStates(1),:&
                     & currentNumTotalStates(2))
               accumulatedMomY(:currentNumTotalStates(1),cumulDataSize + 1:&
                     & cumulDataSize + currentNumTotalStates(2)) = &
                     & pairXWaveFn12MomY(:currentNumTotalStates(1),:&
                     & currentNumTotalStates(2))
               accumulatedMomZ(:currentNumTotalStates(1),cumulDataSize + 1:&
                     & cumulDataSize + currentNumTotalStates(2)) = &
                     & pairXWaveFn12MomZ(:currentNumTotalStates(1),:&
                     & currentNumTotalStates(2))
            endif

            ! This atom/cell pair will contribute to the overlap and hamiltonian
            !   of the i loop atom.  Record the necessary parameters to track 
            !   that information.
            currentLoopIndices(1,totalDataCounter) = j
            currentLoopIndices(2,totalDataCounter) = k

            ! Increment the cumulative record for the amount of data stored in
            !   memory so far.  This will be compared to the target chunk size
            !   to determine when to store the data.
            cumulDataSize = cumulDataSize + currentNumTotalStates(2)

            ! If the amount of data saved in memory so far exceeds the target
            !   chunk size then we will save the data, and reset the
            !   necessary counters for the next chunk.
            if (cumulDataSize .ge. targetChunkSize) then

               call saveCurrentAccumulation(doINTG,doMOME)

               ! Set a flag if the data was recorded as an exact fit.  It is
               !   possible that the last bit of data will not be large enough
               !   to make the cumulDataSize grow larger than the
               !   targetChunkSize.  In such a case it is necessary to do an
               !   extra saveCurrentAccumulation call after we exit the j loop.
               !   This is the assumed default case.  However, it is possible
               !   that the cumulDataSize will match the targetChunkSize and
               !   there will be no more data accumulated that will need the
               !   extra explicit saveCurrentAccumulation call.  In such a case
               !   we need this flag to be set so that the extra call will not
               !   be performed.
               if (cumulDataSize .eq. targetChunkSize) then
                  exactFit = 1
               else
                  exactFit = 0
               endif

               ! Reset the totalDataCounter and cumulative data size stored.
               totalDataCounter = 0
               cumulDataSize   = 0
            endif
         enddo !(k superlattice)
      enddo ! (Atom loop #2)

      ! Make an explicit call to save the accumulated data if there is any need.
      if (exactFit .eq. 0) then
         call saveCurrentAccumulation(doINTG,doMOME)
      endif

      ! Reset the totalDataCounter and cumulative data size stored.
      totalDataCounter = 0
      cumulDataSize   = 0

      ! Close up the HDF data for this atom.
      call closeAtomHDF (doINTG,doMOME)

      ! Mark the completion of this atom.
      if (mod(i,10) .eq. 0) then
         write (20,ADVANCE="NO",FMT="(a1)") "|"
      else
         write (20,ADVANCE="NO",FMT="(a1)") "."
      endif
      if (mod(i,50) .eq. 0) then
         write (20,*) " ",i
      endif
      call flush (20)

   enddo    ! (Atom loop #1)


   ! Deallocate all arrays and matrices before exiting.
   deallocate (currentBasisFns)
   deallocate (currentAlphas)
   deallocate (currentlmAlphaIndex)
   deallocate (currentlmIndex)
   deallocate (alphaDist)   ! Can be saved from before
   deallocate (alphaCenter) ! Can be saved for later.

   if (doINTG == 1) then
      deallocate (pairXWaveFn2OL)
      deallocate (pairXWaveFn2Ham)
      deallocate (pairXWaveFn12OL)
      deallocate (pairXWaveFn12Ham)
   endif

   if (doMOME == 1) then
      deallocate (pairXWaveFn2MomX)
      deallocate (pairXWaveFn2MomY)
      deallocate (pairXWaveFn2MomZ)
      deallocate (pairXWaveFn12MomX)
      deallocate (pairXWaveFn12MomY)
      deallocate (pairXWaveFn12MomZ)
   endif

   ! Log the date and time we start.
   call timeStampEnd(20)

end subroutine intgAndOrMom



subroutine initAtomHDF (i,doINTG,doMOME)

   ! Import the necessary modules
   use HDF5
   use O_Kinds
   use O_Constants,   only: dim3
   use O_Potential,   only: spin
   use O_Basis,       only: initializeAtomSite
   use O_AtomicSites, only: numAtomSites
   use O_AtomicTypes, only: maxNumAtomAlphas, maxNumStates
   use O_Lattice,     only: logBasisFnThresh, numCellsReal, cellSizesReal, &
         & cellDimsReal, findLatticeVector
   use O_PSCFIntgHDF5 ! Basically it needs almost everything from this so...

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the data structure variables passed to this subroutine
   integer :: i ! Current atom 1 index number
   integer :: doINTG
   integer, intent(in) :: doMOME

   ! Define local variables needed for HDF file computations
   character*30 :: currentName
   integer :: j,k ! Loop control variables
   integer :: hdferr
   integer :: totalDataCounter
   integer :: maxNumDataSizeElements
   integer :: numChunksNeeded
   integer :: traversedDataSize
   integer :: maximalChunkSize
   integer, allocatable, dimension (:) :: cumulDataSize

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

   ! Local position and direction vectors and radii
   real (kind=double), dimension (dim3) :: latticeVector ! Vector to lattice
         ! point closest to the difference between the unit cell positions for
         ! atom 1 and atom 2.
   real (kind=double), dimension (dim3) :: shiftedAtomPos ! The position of
         ! atom 2 shifted to each relevant lattice point.
   real (kind=double) :: atomSiteSepSqrd ! The square of the minimum distance
         ! seperating atom 1 and atom 2 according to their unit cell positions
         ! shifted by the lattice point closest to their difference.
   real (kind=double) :: currentNegligLimit ! The distance beyond which all
         ! alpha pairs are considered to have no overlap.
   real (kind=double) :: shiftedAtomSiteSep ! The seperation distance between
         ! atom 1 and the shifted position of atom 2.
   real (kind=double) :: maxLatticeRadius ! Maximum radius beyond which no
         ! lattice points will be considered for integration.  This is used
         ! for the atomic alpha pairs.

   ! Allocate space for locally defined allocatable arrays
   allocate (currentBasisFns        (maxNumAtomAlphas,maxNumStates,2))
   allocate (currentAlphas         (maxNumAtomAlphas,2))
   allocate (currentlmAlphaIndex   (maxNumAtomAlphas,2))
   allocate (currentlmIndex        (maxNumStates,2))


   ! Allocate sufficient space to hold the loop index information.
   allocate (loopIndices(2,numAtomSites * numCellsReal))
   allocate (cumulDataSize(numAtomSites * numCellsReal))


   ! Obtain local copies of key data from larger global data structures for
   !   the first looped atom.
   call initializeAtomSite(i,1,currentAtomType,currentElements,&
         & currentNumTotalStates,currentNumCoreStates,currentNumValeStates,&
         & currentNumAlphas,currentlmIndex,currentlmAlphaIndex,&
         & currentPosition,currentAlphas,currentBasisFns)

   ! Initialize the loop indices and the counters for the number of data
   !   elements, and the total size that those elements contribute.
   loopIndices(:,:)      = 0
   totalDataCounter      = 0
   loopIndexChunkSize(1) = 2
   loopIndexChunkSize(2) = 0
   totalLoopIndexDims(1) = 2
   totalLoopIndexDims(2) = 0
   dataChunkSize(1)      = currentNumTotalStates(1)
   dataChunkSize(2)      = 0
   totalDataDims(1)      = currentNumTotalStates(1)
   totalDataDims(2)      = 0
   cumulDataSize(:)      = 0
   maximalChunkSize      = 0
   numChunksNeeded       = 0

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

      ! Determine the square of the minimum seperation distance between the two
      !   atoms.
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
      maxLatticeRadius = atomSiteSepSqrd + currentNegligLimit + 2.0_double * &
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

         ! Increment the counter for the number of data components for the i
         !   loop atom.
         totalDataCounter = totalDataCounter + 1

         ! Increase the totalDataDims according to the number of states for the
         !   j loop atom.
         totalDataDims(2) = totalDataDims(2) + &
               & currentNumTotalStates(2)

         ! Keep a cumulative record of the size of the data for each j and k
         !   loop iteration.  Later, if the chunk size needs to be smaller
         !   than the total data for this i atom then this will be used to
         !   help split.
         cumulDataSize(totalDataCounter) = totalDataDims(2)

         ! This atom/cell pair will contribute to the overlap and hamiltonian
         !   of the i loop atom.  Record the necessary parameters to track that
         !   information.
         loopIndices(1,totalDataCounter) = j
         loopIndices(2,totalDataCounter) = k

      enddo
   enddo


   ! Determine the size of the chunk for this i loop atom's data and loop
   !   index data.

   ! What is the upper limit on the number of array elements contributed from
   !   the j and k loop atoms?  maxMemory is the maximum number of elements in
   !   total.  currentNumTotalStates(1) is the number of elements contributed
   !   from the i loop atom.  The value of maxMemory was chosen to try and
   !   force the chunk sizes to be large to reduce the HDF overhead of having
   !   many many chunks.  As an example, consider that the atom in the i loop
   !   has 10 states so currentNumTotalStates(1)=10.  Then the
   !   maxNumDataSizeElements = 1000000 if maxMemory is 10000000.
   maxNumDataSizeElements = (maxMemory / currentNumTotalStates(1))

   ! Correct this limit upwards such that it is an exact multiple of the
   !   currentNumTotalStates(1).  Continuing the example no correction is
   !   needed.
   if (mod(maxMemory,int(currentNumTotalStates(1),hsize_t)) .ne. 0) then
      maxNumDataSizeElements = maxNumDataSizeElements + &
            & currentNumTotalStates(1) - mod(maxMemory,&
            & int(currentNumTotalStates(1),hsize_t))
   endif

   ! Determine the number of chunks needed to hold this data without making
   !   the memory footprint too big.  Increasing the chunk size will cause the
   !   memory needed to hold that chunk before it is written to be larger.
   !   Now, in the example, the number of chunks needed will be the
   !   totalDataDims divided by 1000000.  Essentially, the chunk size is about
   !   one million of the atom i interactions.
   numChunksNeeded = int(totalDataDims(2) / maxNumDataSizeElements)

   ! Correct this number upwards if the division was not an exact multiple.
   if (mod(totalDataDims(2),int(maxNumDataSizeElements,hsize_t)) .ne. 0) then
      numChunksNeeded = numChunksNeeded + 1
   endif

   ! Commenting on the example we first need to understand the origin of the
   !   totalDataDims(2) value.  It is the wave function dimension of all the
   !   atoms from every periodic cell that atom i interacts with.  So, for a
   !   system with a small unit cell we note that there will be a small number
   !   of atoms in the unit cell, but a lot of periodic cell.  For a system
   !   with a large unit cell there will be a large number of atoms in the
   !   unit cell (probably) and a small number of (but >= 27) periodic cells.
   !   These effects somewhat cancel each other and the two cases are similar.

   ! So the number of chunks that are needed will probably be rather small,
   !   maybe around 10.

   ! Now that we know the number of chunks we need, we can determine the exact
   !   size of the chunk.  This is actually rather tricky because the data to
   !   be stored in the chunk has internal structure itself.  We can't cut
   !   one of the internal subChunks in half across two HDF chunks.  The
   !   solution that we will persue is not an optimal solution, but rather an
   !   approximation.  A sub chunk is a set of data corresponding to the
   !   interaction of the loop i atom with a loop j atom in a loop k cell.
   ! 1)Make an approximate guess for the chunk size regardless of subChunks.
   !   Initialize this as the maximalChunkSize.
   ! 2)Traverse the cumulDataSize array to find the subChunk border that is
   !   just larger than the guessed chunk size.  Record this as the maximal
   !   chunk size.  Retain the original guessed chunk size as it was in
   !   dataChunkSize.
   ! 3)Continue the search along the cumulDataSize array using the guessed
   !   chunk size to find the next subChunk border that is just larger than the
   !   guessed chunk size.  If this new chunk is larger than the previous
   !   maximal chunk, then it becomes the maximal chunk size.
   ! 4)Repeat 3 until the cumulDataSize array has been totally traversed.  In
   !   the case of only 1 chunk, the array will be traversed in #2.

   ! Now, in the above process we will also record the exact j and k loop
   !   indices where the chunk will end.  Then, in the computation when those
   !   specific j and k loop indices are reached the chunk will be recorded.
   !   This will leave a small amount of space at the end of each chunk, but
   !   that is an acceptable loss since it should only be a few bytes per
   !   chunk, and the total number of chunks will be small even for large
   !   systems with many atoms.

   ! Make chunk size guess
   dataChunkSize(2) = int((totalDataDims(2) / numChunksNeeded))

   ! Correct this number upwards if the division was not an exact multiple.
   if (mod(int(totalDataDims(2),hsize_t),numChunksNeeded) .ne. 0) then
      dataChunkSize(2) = dataChunkSize(2) + 1
   endif

   ! Initialize the maximal chunk size
   maximalChunkSize = dataChunkSize(2)

   ! Initialize a counter for the amount of data already traversed.
   traversedDataSize = 0

   ! Traverse the cumulDataSize array.
   do j = 1, totalDataCounter
      if (cumulDataSize(j) - traversedDataSize .ge. &
            & dataChunkSize(2)) then

         if (cumulDataSize(j) - traversedDataSize .ge. maximalChunkSize) then

            ! Obtain the maximal chunk size
            maximalChunkSize = cumulDataSize(j) - traversedDataSize

            ! Record the amount of data seen so far in total.
            traversedDataSize = cumulDataSize(j)

            ! Is the number of indices for this chunk the largest number of
            !   indices needed so far?
            if (j - loopIndexChunkSize(2) .gt. loopIndexChunkSize(2)) then

               ! Record the number of indices seen for this chunk since it is
               !   the maximal number.
               loopIndexChunkSize(2) = j - loopIndexChunkSize(2)
            endif
         endif
      endif
   enddo

   ! At this point the maximal chunk size has been determined.  This will be
   !   used as the chunk size for the datasets below.  The basic scheme will
   !   be to have both a target chunk size and a maximal chunk size.  The
   !   maximal chunk size is the upper limit, and the target chunk size is the
   !   lower limit.  As the calculation proceeds, whenever the data equals or
   !   passes the target chunk size the chunk will be stored on disk and it
   !   will occupy a chunk of size maximal chunk size.  (i.e. Each chunk that
   !   is stored will have a small tail of unused data attached to it.)  The
   !   target chunk size is just the initial chunk size guess from above.

   ! Record the target chunk size.
   targetChunkSize = dataChunkSize(2)

   ! Record the actual chunk size to use for the file.
   dataChunkSize(2) = maximalChunkSize

   ! Adjust the totalDataDims to reflect the number of chunks needed and the
   !   maximal chunk size.
   totalDataDims(2)      = numChunksNeeded * maximalChunkSize
   totalLoopIndexDims(2) = numChunksNeeded * loopIndexChunkSize(2)

   ! At this point the number of j loop atoms and corresponding real space
   !   cells necessary for the current i loop atom is known.  We also know the
   !   chunk size that is needed to keep the memory footprint sufficiently
   !   small.  So, we can make the HDF objects and descriptors for the file.


   ! Make the dataspace for the file.
   if (doINTG == 1) then
      do j = 1, spin
         call h5screate_simple_f(2,totalDataDims,fileHamiltonian_dsid(j),hdferr)
         if (hdferr /= 0) stop 'Failed to create file ham dsid'
      enddo
      call h5screate_simple_f (2,totalDataDims,fileOverlap_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to create file overlap dsid'
   endif
   if (doMOME == 1) then
      call h5screate_simple_f (2,totalDataDims,fileMomentumX_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to create file momentum x dsid'
      call h5screate_simple_f (2,totalDataDims,fileMomentumY_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to create file momentum y dsid'
      call h5screate_simple_f (2,totalDataDims,fileMomentumZ_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to create file momentum z dsid'
   endif


   ! Make the property list for the current i atom's integral datasets.
   call h5pcreate_f       (H5P_DATASET_CREATE_F,fileData_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create file data plid'
   call h5pset_layout_f   (fileData_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set file data layout type'
   call h5pset_chunk_f    (fileData_plid,2,dataChunkSize,hdferr)
   if (hdferr /= 0) stop 'Failed to set file data chunk size'
!   call h5pset_shuffle_f (fileData_plid,hdferr)
   call h5pset_deflate_f  (fileData_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set file data for deflation'

   ! Make the datasets for the file.
   if (doINTG == 1) then
      do j = 1, spin
         write (currentName,*) i,j
         currentName = trim (currentName)
         call h5dcreate_f (hamiltonian_gid,currentName,H5T_NATIVE_DOUBLE,&
               & fileHamiltonian_dsid(j),hamiltonian_did(j),&
               & hdferr,fileData_plid)
         if (hdferr /= 0) stop 'Failed to create hamiltonian dataset'
      enddo
      write (currentName,*) i
      currentName = trim (currentName)
      call h5dcreate_f (overlap_gid,currentName,H5T_NATIVE_DOUBLE,&
            & fileOverlap_dsid,overlap_did,hdferr,fileData_plid)
      if (hdferr /= 0) stop 'Failed to create overlap dataset'
   endif
   if (doMOME == 1) then
      write (currentName,*) i
      currentName = trim (currentName)
      call h5dcreate_f (momentumX_gid,currentName,H5T_NATIVE_DOUBLE,&
            & fileMomentumX_dsid,momentumX_did,hdferr,fileData_plid)
      if (hdferr /= 0) stop 'Failed to create momentum x dataset'
      call h5dcreate_f (momentumY_gid,currentName,H5T_NATIVE_DOUBLE,&
            & fileMomentumY_dsid,momentumY_did,hdferr,fileData_plid)
      if (hdferr /= 0) stop 'Failed to create momentum y dataset'
      call h5dcreate_f (momentumZ_gid,currentName,H5T_NATIVE_DOUBLE,&
            & fileMomentumZ_dsid,momentumZ_did,hdferr,fileData_plid)
      if (hdferr /= 0) stop 'Failed to create momentum z dataset'
   endif

   ! Make the dataspace for the memory reference.  This dataspace will be the
   !   same for all the different datasets of the current i atom.
   call h5screate_simple_f (2,dataChunkSize,dataChunk_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create data chunk dsid'

   ! It is only necessary to record the index tracking information when doing
   !   the overlap and hamiltonian integrals.  When doing *only* the momentum
   !   integrals it is *assumed* (certain) that the index information was
   !   already computed for the overlap and hamiltonian (because they should
   !   have been done first).
   if (doINTG == 1) then
      ! Make the dataspace for the file to track the index number of the atom2
      !   loop and the lattice loop for reading in the data later.
      call h5screate_simple_f (2,totalLoopIndexDims,fileLoopIndices_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to create file loop indices dsid'

      ! Make the property list for the index tracking dataset for this i atom.
      call h5pcreate_f       (H5P_DATASET_CREATE_F,fileLoopIndices_plid,hdferr)
      if (hdferr /= 0) stop 'Failed to create file loop indices plid'
      call h5pset_layout_f   (fileLoopIndices_plid,H5D_CHUNKED_F,hdferr)
      if (hdferr /= 0) stop 'Failed to set file loop indices layout'
      call h5pset_chunk_f    (fileLoopIndices_plid,2,loopIndexChunkSize,hdferr)
      if (hdferr /= 0) stop 'Failed to set file loop indices chunk size'
!      call h5pset_shuffle_f (fileLoopIndices_plid,hdferr)
      call h5pset_deflate_f  (fileLoopIndices_plid,1,hdferr)
      if (hdferr /= 0) stop 'Failed to set file loop indices for deflation'


      ! Make a dataset to track the index numbers of the atom2 loop and the
      !   lattice loop for reading in the data later.
      call h5dcreate_f (loopIndices_gid,currentName,H5T_NATIVE_INTEGER,&
            & fileLoopIndices_dsid,loopIndices_did,hdferr,fileLoopIndices_plid)
      if (hdferr /= 0) stop 'Failed to create loop indices data set'

      ! Make the dataspace for the memory reference of the loop index matrix.
      call h5screate_simple_f (2,loopIndexChunkSize,loopIndexChunk_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to create loop index chunk size dsid'
   endif

   ! Initialize the counters for the currentStart of the loopIndex and data
   !   dataset/dataspace for referencing and writing.
   currentStart(1) = 0
   currentStart(2) = 0
   currentIndexStart(1) = 0
   currentIndexStart(2) = 0

   ! Deallocate unnecessary matrices now.
   deallocate (currentBasisFns)
   deallocate (currentAlphas)
   deallocate (currentlmAlphaIndex)
   deallocate (currentlmIndex)
   deallocate (loopIndices)
   deallocate (cumulDataSize)

   ! Allocate current loop index space and initialize it.
   allocate (currentLoopIndices (int(loopIndexChunkSize(1)), &
         & int(loopIndexChunkSize(2))))
   currentLoopIndices(:,:) = 0

   ! Allocate accumulation space for the overlap and hamiltonian if required.
   if (doINTG == 1) then
      allocate (accumulatedOL  (int(dataChunkSize(1)),int(dataChunkSize(2))))
      allocate (accumulatedHam (int(dataChunkSize(1)),&
            & int(dataChunkSize(2)),spin))

      ! Initialize matrices.
      accumulatedOL(:,:)      = 0.0_double
      accumulatedHam(:,:,:)   = 0.0_double
   endif

   ! Allocate space for the momentum matrix components if requested.
   if (doMOME == 1) then

      ! Allocate accumulation space
      allocate (accumulatedMomX (int(dataChunkSize(1)),int(dataChunkSize(2))))
      allocate (accumulatedMomY (int(dataChunkSize(1)),int(dataChunkSize(2))))
      allocate (accumulatedMomZ (int(dataChunkSize(1)),int(dataChunkSize(2))))

      ! Initialize matrices
      accumulatedMomX(:,:)    = 0.0_double
      accumulatedMomY(:,:)    = 0.0_double
      accumulatedMomZ(:,:)    = 0.0_double
   endif

end subroutine initAtomHDF


subroutine closeAtomHDF (doINTG,doMOME)

   ! Import the necessary modules
   use HDF5
   use O_Kinds
   use O_Potential, only: spin
   use O_PSCFIntgHDF5 ! Basically it needs almost everything from here so...

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine
   integer :: doINTG
   integer, intent(in) :: doMOME

   ! Define local variables
   integer :: hdferr
   integer :: i

   ! Now we close the dataspaces used for the file output.
   if (doINTG == 1) then
      do i = 1, spin
         call h5sclose_f (fileHamiltonian_dsid(i),hdferr)
         if (hdferr /= 0) stop 'Failed to close file hamiltonian dsid'
      enddo
      call h5sclose_f (fileOverlap_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to close file overlap dsid'
   endif
   if (doMOME == 1) then
      call h5sclose_f (fileMomentumX_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to close file momentum x dsid'
      call h5sclose_f (fileMomentumY_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to close file momentum y dsid'
      call h5sclose_f (fileMomentumZ_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to close file momentum z dsid'
   endif

   ! Close the data space for the loop indices.
   if (doINTG == 1) then
      call h5sclose_f (fileLoopIndices_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to close file loop indices dsid'
   endif

   ! Close the dataspaces used for the memory representation of the chunk.
   if (doINTG == 1) then
      call h5sclose_f (loopIndexChunk_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to close loop index chunk dsid'
   endif
   call h5sclose_f (dataChunk_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close data chunk dsid'

   ! Close the property list used for the file output.
   if (doINTG == 1) then
      call h5pclose_f (fileLoopIndices_plid,hdferr)
      if (hdferr /= 0) stop 'Failed to close file loop indices plid'
   endif
   call h5pclose_f (fileData_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close file data plid'

   ! Close the datasets used for the file output.
   if (doINTG == 1) then
      do i = 1, spin
         call h5dclose_f (hamiltonian_did(i),hdferr)
         if (hdferr /= 0) stop 'Failed to close hamiltonian did'
      enddo
      call h5dclose_f (overlap_did,hdferr)
      if (hdferr /= 0) stop 'Failed to close overlap did'
   endif
   if (doMOME == 1) then
      call h5dclose_f (momentumX_did,hdferr)
      if (hdferr /= 0) stop 'Failed to close momentum x did'
      call h5dclose_f (momentumY_did,hdferr)
      if (hdferr /= 0) stop 'Failed to close momentum y did'
      call h5dclose_f (momentumZ_did,hdferr)
      if (hdferr /= 0) stop 'Failed to close momentum z did'
   endif

   ! Close the dataset for the loop indices.
   if (doINTG == 1) then
      call h5dclose_f (loopIndices_did,hdferr)
      if (hdferr /= 0) stop 'Failed to close loop indices did'
   endif

   ! Deallocate loop indices.
   deallocate (currentLoopIndices)

   ! Deallocate accumulation matrices.
   if (doINTG == 1) then
      deallocate (accumulatedOL)
      deallocate (accumulatedHam)
   endif
   if (doMOME == 1) then
      deallocate (accumulatedMomX)
      deallocate (accumulatedMomY)
      deallocate (accumulatedMomZ)
   endif

end subroutine closeAtomHDF




subroutine saveCurrentAccumulation(doINTG,doMOME)

   ! Import the necessary modules
   use HDF5
   use O_Kinds
   use O_Potential, only: spin
   use O_PSCFIntgHDF5 ! It needs almost everything from here so...

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the data structure variables passed to this subroutine
   integer :: doINTG
   integer, intent(in) :: doMOME

   ! Define local variables
   integer :: hdferr
   integer :: i

   ! Record loop indices in the doINTG case only.  For the doMOME *only* case,
   !   this *should* already be done and saved in the HDF5 file.
   if (doINTG == 1) then
      ! Choose the hyperslab to write to in the loop indices dataset.
      call h5sselect_hyperslab_f (fileLoopIndices_dsid,H5S_SELECT_SET_F,&
            & currentIndexStart,loopIndexChunkSize,hdferr)
      if (hdferr /= 0) stop 'Failed to select hyperslab for file loop indices'

      ! Now, record the data.
      call h5dwrite_f (loopIndices_did,H5T_NATIVE_INTEGER,&
            & currentLoopIndices(:,:),loopIndexChunkSize,hdferr,&
            & loopIndexChunk_dsid,fileLoopIndices_dsid)
      if (hdferr /= 0) stop 'Failed to write loop indices did'
   endif
   

   ! Perform the extention and recording process for the overlap and
   !   hamiltonian if necessary.
   if (doINTG == 1) then

      ! Choose the hyperslab to write to in overlap and hamiltonian datasets.
      do i = 1, spin
         call h5sselect_hyperslab_f (fileHamiltonian_dsid(i),H5S_SELECT_SET_F,&
               & currentStart,dataChunkSize,hdferr)
         if (hdferr /= 0) stop 'Failed to select hyperslab for file hamiltonian'
      enddo
      call h5sselect_hyperslab_f (fileOverlap_dsid,H5S_SELECT_SET_F,&
            & currentStart,dataChunkSize,hdferr)
      if (hdferr /= 0) stop 'Failed to select hyperslab for file overlap'

      ! Now record the data.
      do i = 1, spin
         call h5dwrite_f (hamiltonian_did(i),H5T_NATIVE_DOUBLE,&
               & accumulatedHam(:,:,i),dataChunkSize,hdferr,dataChunk_dsid,&
               & fileHamiltonian_dsid(i))
         if (hdferr /= 0) stop 'Failed to write hamiltonian did'
      enddo
      call h5dwrite_f (overlap_did,H5T_NATIVE_DOUBLE,accumulatedOL(:,:),&
            & dataChunkSize,hdferr,dataChunk_dsid,fileOverlap_dsid)
      if (hdferr /= 0) stop 'Failed to write overlap did'
   endif


   ! If we need to compute the momentum, then repeat the extention
   !   and recording process for that data too.
   if (doMOME == 1) then

      ! Choose the hyperslab to write to in the file dataset.
      call h5sselect_hyperslab_f (fileMomentumX_dsid,&
            & H5S_SELECT_SET_F,currentStart,dataChunkSize,hdferr)
      if (hdferr /= 0) stop 'Failed to select hyperslab for file momentum x'
      call h5sselect_hyperslab_f (fileMomentumY_dsid,&
            & H5S_SELECT_SET_F,currentStart,dataChunkSize,hdferr)
      if (hdferr /= 0) stop 'Failed to select hyperslab for file momentum y'
      call h5sselect_hyperslab_f (fileMomentumZ_dsid,&
            & H5S_SELECT_SET_F,currentStart,dataChunkSize,hdferr)
      if (hdferr /= 0) stop 'Failed to select hyperslab for file momentum z'

      ! Now record the data.
      call h5dwrite_f (momentumX_did,H5T_NATIVE_DOUBLE,accumulatedMomX(:,:),&
            & dataChunkSize,hdferr,dataChunk_dsid,fileMomentumX_dsid)
      if (hdferr /= 0) stop 'Failed to write momentum x data'
      call h5dwrite_f (momentumY_did,H5T_NATIVE_DOUBLE,accumulatedMomY(:,:),&
            & dataChunkSize,hdferr,dataChunk_dsid,fileMomentumY_dsid)
      if (hdferr /= 0) stop 'Failed to write momentum y data'
      call h5dwrite_f (momentumZ_did,H5T_NATIVE_DOUBLE,accumulatedMomZ(:,:),&
            & dataChunkSize,hdferr,dataChunk_dsid,fileMomentumZ_dsid)
      if (hdferr /= 0) stop 'Failed to write momentum z data'
   endif

   ! Now we need to prepare for the next large chunk.

   ! Increment the start parameters so that we know where to start with the
   !   hyperslab selection next time.
   currentStart(2)      = currentStart(2)      + dataChunkSize(2)
   currentIndexStart(2) = currentIndexStart(2) + loopIndexChunkSize(2)

   ! Clear the current loop index data
   currentLoopIndices(:,:) = 0

   ! Clear the accumulation data.
   if (doINTG == 1) then
      accumulatedOL(:,:)    = 0.0_double
      accumulatedHam(:,:,:) = 0.0_double
   endif

   ! Take care of the momentum data if necessary
   if (doMOME == 1) then
      ! Clear the accumulation data.
      accumulatedMomX(:,:) = 0.0_double
      accumulatedMomY(:,:) = 0.0_double
      accumulatedMomZ(:,:) = 0.0_double
   endif
end subroutine saveCurrentAccumulation


#ifndef GAMMA
subroutine getIntgResults (valeVale,coreValeOL,&
      & currentKPoint,runCode,valeValeBand_did,valeValeBand_dims,&
      & noSaveValeVale,spinDirection)
#else
subroutine getIntgResults (valeValeGamma,coreValeOLGamma,&
      & runCode,valeValeBand_did,valeValeBand_dims,&
      & noSaveValeVale,spinDirection)
#endif

   ! Import the necessary modules
   use HDF5
   use O_Kinds
   use O_Constants, only: dim3
   use O_PSCFIntgHDF5 ! It needs almost everything from here so...
   use O_AtomicSites, only: coreDim, valeDim, numAtomSites, atomSites
   use O_AtomicTypes, only: maxNumStates, atomTypes
   use O_Lattice, only: findLatticeVector
   use O_KPoints
   use O_Orthogonalization
#ifndef GAMMA
   use O_IntgSaving, only: kPointLatticeOriginShift, saveCurrentPair, &
         & applyPhaseFactors
#else
   use O_IntgSaving, only: saveCurrentPairGamma, applyPhaseFactorsGamma
#endif

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the variables passed to this subroutine.
#ifndef GAMMA
   complex (kind=double), dimension (:,:,:) :: valeVale
   complex (kind=double), dimension (:,:,:) :: coreValeOL
   integer :: currentKPoint
#else
   real (kind=double), dimension (:,:) :: valeValeGamma
   real (kind=double), dimension (:,:) :: coreValeOLGamma
#endif
   integer :: runCode  ! 1=overlap; 2=hamiltonian; 3,4,5=MOME x,y,z
   integer (hid_t) :: valeValeBand_did
   integer (hsize_t), dimension (2) :: valeValeBand_dims
   integer :: noSaveValeVale
   integer :: spinDirection

   ! Define local variables for logging and loop control
   integer :: i,j,k ! Loop index variables
   integer :: hdferr

   ! Define local variables for tracking chunks, indices, etc.
   character*30 :: currentName
   integer :: numChunks
   integer(hsize_t) :: numPoints
   integer :: atom2Index
   integer :: accumInitPos

   ! Atom specific variables that change with each atom pair loop iteration.
   integer, dimension (2) :: currentNumTotalStates

   ! Variables and data structures that change or are accumulated with each
   !   iteration of the lattice loop.
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: currentPair
#else
   real (kind=double), allocatable, dimension (:,:)      :: currentPairGamma
#endif
   real (kind=double), allocatable, dimension (:,:)      :: accumulatedIntg
   integer,            allocatable, dimension (:,:)      :: loopIndices

   ! Variables and data structures that are used to calculate the overlap
   !   valeVale matrix.
   integer :: currIndex
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: coreCore
   complex (kind=double), allocatable, dimension (:,:,:) :: coreVale
   complex (kind=double), allocatable, dimension (:,:)   :: valeCore
#else
   real (kind=double), allocatable, dimension (:,:)   :: coreCoreGamma
   real (kind=double), allocatable, dimension (:,:)   :: coreValeGamma
   real (kind=double), allocatable, dimension (:,:)   :: valeCoreGamma
#endif
   real (kind=double), allocatable, dimension (:,:)   :: packedValeVale

   ! Local position and direction vectors and radii
   real (kind=double), dimension (dim3) :: latticeVector ! Vector to lattice
         ! point closest to the difference between the unit cell positions for
         ! atom 1 and atom 2.

   ! Initialize the HDF interface for the components needed here.
   call initReadIntgHDF5 (runCode)

   ! Allocate space for locally defined allocatable arrays.
#ifndef GAMMA
   allocate (currentPair(maxNumStates,maxNumStates,1)) ! Need 1KP at a time.
   allocate (coreCore(coreDim,coreDim,1)) ! Need 1 KPoint at a time only.
   if (runCode /= 1) then
      allocate (coreVale(coreDim,valeDim,1)) ! Need 1 KPoint at a time only.
   endif
# else
   allocate (currentPairGamma(maxNumStates,maxNumStates))
   allocate (coreCoreGamma(coreDim,coreDim))
   if (runCode /= 1) then
      allocate (coreValeGamma(coreDim,valeDim))
   endif
#endif


   ! Initialize key matrices
#ifndef GAMMA
   coreCore(:,:,1) = 0.0_double
   valeVale(:,:,1) = 0.0_double
   if (runCode == 1) then
      coreValeOL(:,:,1) = 0.0_double
   else
      coreVale(:,:,1) = 0.0_double
   endif
#else
   coreCoreGamma(:,:) = 0.0_double
   valeValeGamma(:,:) = 0.0_double
   if (runCode == 1) then
      coreValeOLGamma(:,:) = 0.0_double
   else
      coreValeGamma(:,:) = 0.0_double
   endif
#endif


   ! For each atom read in the requested results of the atom2-lattice loop
   !   matrices from the intg calculation.
   do i = 1, numAtomSites

      ! Initialize the current atom2 index number and the number of total
      !   states for this atom pair.
      atom2Index = i
      currentNumTotalStates(1) = atomTypes(atomSites(i)%atomTypeAssn)%&
            & numCoreStates + atomTypes(atomSites(i)%atomTypeAssn)%numValeStates
      currentNumTotalStates(2) = currentNumTotalStates(1)

      ! Reset the currentPair
#ifndef GAMMA
      currentPair    (:,:,:) = 0.0_double
#else
      currentPairGamma (:,:) = 0.0_double
#endif


      ! Find the lattice point closest to the difference between the two
      !   atom sites.
      call findLatticeVector((atomSites(i)%cartPos(:)-&
            & atomSites(atom2Index)%cartPos(:)),latticeVector)

      ! Obtain the name for the datasets for this atom.
      if (runCode == 2) then
         write (currentName,*) i,spinDirection
         currentName = trim (currentName)

         ! Open the integral dataset for this atom.
         call h5dopen_f (intg_gid,currentName,intg_did,hdferr)
         if (hdferr /= 0) stop 'Can not open integral dataset'

         ! The loopIndex information does not depend on spin direction.
         write (currentName,*) i
         currentName = trim (currentName)

         ! Open the loopIndex dataset for this atom.
         call h5dopen_f (loopIndices_gid,currentName,loopIndices_did,hdferr)
         if (hdferr /= 0) stop 'Can not open loop index dataset'
      else
         write (currentName,*) i
         currentName = trim (currentName)

         ! Open the integral dataset for this atom.
         call h5dopen_f (intg_gid,currentName,intg_did,hdferr)
         if (hdferr /= 0) stop 'Can not open integral dataset'

         ! Open the loopIndex dataset for this atom.
         call h5dopen_f (loopIndices_gid,currentName,loopIndices_did,hdferr)
         if (hdferr /= 0) stop 'Can not open loop index dataset'
      endif

      ! Determine the number of chunks that make up this dataset.

      ! Create the dataspace that will describe the overlap file to be read.
      !   The dataspace covers all the data on the file.  The exact data to
      !   read will be selected in each j-loop iteration via hyperslab.
      call h5dget_space_f (intg_did,fileIntg_dsid,hdferr)
      if (hdferr /= 0) stop 'Can not get space'
      call h5dget_create_plist_f (intg_did,fileIntg_plid,hdferr)
      if (hdferr /= 0) stop 'Can not get create plist'
      call h5pget_chunk_f (fileIntg_plid,2,chunkDims,hdferr)
      if (hdferr == -1) stop 'Failed to get chunk size'
      call h5sget_simple_extent_npoints_f (fileIntg_dsid,numPoints,hdferr)
      if (hdferr == -1) stop 'Failed to get number of points'
      numPoints = numPoints / currentNumTotalStates(1)
      numChunks = numPoints / chunkDims(2)
      if (mod(numPoints,chunkDims(2)) .ne. 0) then
         numChunks = numChunks + 1
      endif

      ! Create the dataspace that will describe the loopIndices file to be
      !   read.  This dataspace covers all the loop index data on the file.
      !   The exact loop indices to read will be chosen in the j-loop.
      call h5dget_space_f (loopIndices_did,fileLoopIndices_dsid,hdferr)
      if (hdferr /= 0) stop 'Can not get loop index space'
      call h5dget_create_plist_f (loopIndices_did,fileLoopIndices_plid,hdferr)
      if (hdferr /= 0) stop 'Can not get loop index create plist'
      call h5pget_chunk_f (fileLoopIndices_plid,2,loopIndexDims,hdferr)
      if (hdferr == -1) stop 'Failed to get loop index chunk size'
      call h5sget_simple_extent_npoints_f (fileLoopIndices_dsid,numPoints,&
            & hdferr)
      if (hdferr == -1) stop 'Failed to get loop index number of points'
      numPoints = numPoints / 2

      ! Create the dataspace that will describe the intg data in memory.
      call h5screate_simple_f (2,chunkDims,dataChunk_dsid,hdferr)
      if (hdferr == -1) stop 'Failed to create dataspace'

      ! Create the dataspace that will describe the index data in memory.
      call h5screate_simple_f (2,loopIndexDims,loopIndices_dsid,hdferr)
      if (hdferr == -1) stop 'Failed to create loop index dataspace'


      ! Allocate space in memory to read the data and indices.
      allocate (accumulatedIntg (int(chunkDims(1)),int(chunkDims(2))))
      allocate (loopIndices (int(loopIndexDims(1)),int(loopIndexDims(2))))

      ! Initialize the starting locations to read from the file.
      currentStart(1)      = 0
      currentStart(2)      = 0
      currentIndexStart(1) = 0
      currentIndexStart(2) = 0


      ! Make a loop to repeatedly read a chunk of data from the datasets for
      !   this atom.  Read both the data and the loop indices.
      do j = 1, numChunks

         ! Choose the hyperslabs to read for this iteration.
         call h5sselect_hyperslab_f (fileIntg_dsid,H5S_SELECT_SET_F,&
               & currentStart,chunkDims,hdferr)
         if (hdferr == -1) stop 'Failed to select hyperslab'
         call h5sselect_hyperslab_f (fileLoopIndices_dsid,H5S_SELECT_SET_F,&
               & currentIndexStart,loopIndexDims,hdferr)
         if (hdferr == -1) stop 'Failed to select loop index hyperslab'

         ! Now read the data

         call h5dread_f (loopIndices_did,H5T_NATIVE_INTEGER,loopIndices,&
               & loopIndexDims,hdferr,loopIndices_dsid,fileLoopIndices_dsid)
         if (hdferr == -1) stop 'Failed to read loop indices'
         call h5dread_f (intg_did,H5T_NATIVE_DOUBLE,accumulatedIntg,&
               & chunkDims,hdferr,dataChunk_dsid,fileIntg_dsid)
         if (hdferr == -1) stop 'Failed to read integral data'

         ! If we are not at the last iteration of the j-loop then we increment
         !   the starting point to the next chunk.
         if (j .ne. numChunks) then
            currentStart(2) = currentStart(2) + chunkDims(2)
            currentIndexStart(2) = currentIndexStart(2) + loopIndexDims(2)
         endif

         ! Initialize the position within the chunk that will be the initial
         !   point to begin operating on the accumulated data.
         accumInitPos = 0

         ! Now that we have read in the data we start to compute on it.
         do k = 1, numPoints

            ! Check if we reached the end of the data within a given chunk.
            if (loopIndices(1,k) == 0) then

               ! Exit the loop over sub-chunks of the current chunk.
               exit
            endif

            ! Check if we are at the next atom from the atom2 loop.
            if (loopIndices(1,k) .ne. atom2Index) then

               ! Collect the currentPair results into the valeVale matrix.
#ifndef GAMMA
               if (runCode .eq. 1) then
                  call kPointLatticeOriginShift (currentNumTotalStates,&
                        & currentPair,latticeVector,1,currentKPoint-1)
                  call saveCurrentPair (i,atom2Index,1,currentPair,&
                        & valeVale,coreValeOL,coreCore)
               else
                  call kPointLatticeOriginShift (currentNumTotalStates,&
                        & currentPair,latticeVector,1,currentKPoint-1)
                  call saveCurrentPair (i,atom2Index,1,currentPair,&
                        & valeVale,coreVale,coreCore)
               endif
#else
               if (runCode .eq. 1) then
                  call saveCurrentPairGamma (i,atom2Index,currentPairGamma,&
                        & valeValeGamma,coreValeOLGamma,coreCoreGamma)
               else
                  call saveCurrentPairGamma (i,atom2Index,currentPairGamma,&
                        & valeValeGamma,coreValeGamma,coreCoreGamma)
               endif
#endif

               ! Update the atom2Index
               atom2Index = loopIndices(1,k)

               ! Reset the currentPair
#ifndef GAMMA
               currentPair    (:,:,:) = 0.0_double
#else
               currentPairGamma (:,:) = 0.0_double
#endif

               ! Update the currentNumTotalStates for the new atom2 index.
               currentNumTotalStates(2) = atomTypes(atomSites(atom2Index)%&
                     & atomTypeAssn)%numCoreStates + atomTypes(atomSites&
                     & (atom2Index)%atomTypeAssn)%numValeStates

               ! Update the latticeVector
               call findLatticeVector((atomSites(i)%cartPos(:)-&
                     & atomSites(atom2Index)%cartPos(:)),latticeVector)

               ! Apply the phase factor for this cell.
#ifndef GAMMA
               call applyPhaseFactors (currentPair,&
                  & accumulatedIntg(1:currentNumTotalStates(1),&
                  & accumInitPos+1:accumInitPos+currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),&
                  & loopIndices(2,k),runCode,currentKPoint)
#else
               call applyPhaseFactorsGamma (currentPairGamma,&
                  & accumulatedIntg(1:currentNumTotalStates(1),&
                  & accumInitPos+1:accumInitPos+currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),runCode)
#endif

            else
               ! Apply the phase factor for this cell.
#ifndef GAMMA
               call applyPhaseFactors (currentPair,&
                  & accumulatedIntg(1:currentNumTotalStates(1),&
                  & accumInitPos+1:accumInitPos+currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),&
                  & loopIndices(2,k),runCode,currentKPoint)
#else
               call applyPhaseFactorsGamma (currentPairGamma,&
                  & accumulatedIntg(1:currentNumTotalStates(1),&
                  & accumInitPos+1:accumInitPos+currentNumTotalStates(2)),&
                  & currentNumTotalStates(1),currentNumTotalStates(2),runCode)
#endif
            endif

            ! Increment the position within the accumulated data where we begin
            !   operating on the data.
            accumInitPos = accumInitPos + currentNumTotalStates(2)
         enddo


         ! Collect the currentPair results into the overlap matrix.
#ifndef GAMMA
         if (runCode .eq. 1) then
            call kPointLatticeOriginShift (currentNumTotalStates,&
                  & currentPair,latticeVector,1,currentKPoint-1)
            call saveCurrentPair (i,atom2Index,1,currentPair,&
                  & valeVale,coreValeOL,coreCore)
         else
            call kPointLatticeOriginShift (currentNumTotalStates,&
                  & currentPair,latticeVector,1,currentKPoint-1)
            call saveCurrentPair (i,atom2Index,1,currentPair,&
                  & valeVale,coreVale,coreCore)
         endif
#else
         if (runCode .eq. 1) then
            call saveCurrentPairGamma (i,atom2Index,currentPairGamma,&
                  & valeValeGamma,coreValeOLGamma,coreCoreGamma)
         else
            call saveCurrentPairGamma (i,atom2Index,currentPairGamma,&
                  & valeValeGamma,coreValeGamma,coreCoreGamma)
         endif
#endif
      enddo

      ! Close the HDF property lists for this atom.
      call h5pclose_f (fileIntg_plid,hdferr)
      if (hdferr == -1) stop 'Failed to close file intg plid'
      call h5pclose_f (fileLoopIndices_plid,hdferr)
      if (hdferr == -1) stop 'Failed to close file loop indices plid'

      ! Close the HDF datasets for this atom.
      call h5dclose_f (intg_did,hdferr)
      if (hdferr == -1) stop 'Failed to close intg did'
      call h5dclose_f (loopIndices_did,hdferr)
      if (hdferr == -1) stop 'Failed to close loop indices did'

      ! Close the HDF dataspaces for this atom.
      call h5sclose_f (fileIntg_dsid,hdferr)
      if (hdferr == -1) stop 'Failed to close intg dsid'
      call h5sclose_f (fileLoopIndices_dsid,hdferr)
      if (hdferr == -1) stop 'Failed to close loop indices dsid'

      ! Deallocate the chunks for this i loop atom.
      deallocate (accumulatedIntg)
      deallocate (loopIndices)
   enddo

   ! Deallocate matrices that are not necessary any more.
#ifndef GAMMA
   deallocate (currentPair)
#else
   deallocate (currentPairGamma)
#endif

   ! Close the dataspace that was used for the data memory dataset.
   call h5sclose_f (dataChunk_dsid,hdferr)
   if (hdferr == -1) stop 'Failed to close data chunk dsid'

   ! Close the dataspace that was used for the loop index memory dataset.
   call h5sclose_f (loopIndices_dsid,hdferr)
   if (hdferr == -1) stop 'Failed to close loop indices dsid'

   ! Close the integral file access property list.
   call h5pclose_f (intg_plid,hdferr)
   if (hdferr == -1) stop 'Failed to close intg plid'

   ! Close the HDF integral file.
   call h5fclose_f (intg_fid,hdferr)
   if (hdferr == -1) stop 'Failed to close intg file'


   ! Form the valeVale matrix
   ! Determine which matrices are to be made and used.  The overlap matrix is
   !   treated slightly differently than the others.
#ifndef GAMMA
   if (runCode .eq. 1) then
      if (coreDim /= 0) then
         ! Obtain the orthogonalization coefficients.  (coreValeOL)
         call valeCoreCoreValeOL (valeDim,coreDim,valeVale(:,:,1),&
               & coreValeOL(:,:,1))

         ! Allocate the valeCore matrix.
         allocate (valeCore (coreDim,valeDim))  ! Pre-transposed format

         ! Form a product of (coreValeOL)(coreCore) in a temp matrix
         !   (valeCore).
         call coreValeCoreCore (valeDim,coreDim,valeCore,coreValeOL(:,:,1),&
               & coreCore(:,:,1))

         ! Deallocate the coreCore and allocate the packedValeVale.
         deallocate (coreCore)
         allocate (packedValeVale(2,valeDim*(valeDim+1)/2))

         ! Finally compute the product of the above (valeCore)(coreCore) with
         !   coreValeOL.  This is added to valeVale to complete
         !   orthogonaliziation.
         call makeValeVale (valeDim,coreDim,valeDim,valeCore, &
               & coreValeOL(:,:,1),valeVale(:,:,1),packedValeVale,1,0)

         ! Save the overlap valeVale only if it is requested.  (This will be
         !   most of the time.  It is only not saved when we are doing a SYBD
         !   type of calculation when we have some 300 kpoints.)  This is
         !   also why this question is not asked for the Gamma kpoint
         !   situation below.
         if (noSaveValeVale == 0) then
            call h5dwrite_f (valeValeBand_did,H5T_NATIVE_DOUBLE,&
                  & packedValeVale(:,:),valeValeBand_dims,hdferr)
            if (hdferr == -1) stop 'Failed to write packed vale vale band'
         endif

         ! Deallocate packedValeVale and valeCore to be ready for next kp.
         deallocate (packedValeVale)
         deallocate (valeCore)
      else  ! No core dim in whole system so orthogonalization is not needed.
         deallocate (coreCore)
         if (noSaveValeVale == 0) then
            allocate (packedValeVale(2,valeDim*(valeDim+1)/2))

            ! Initialize the index counter for packing.
            currIndex = 0

            ! Pack the valeVale matrix.
            do j = 1, valeDim
               do k = 1, j
                  currIndex = currIndex + 1
                  packedValeVale(1,currIndex) = &
                        & real(valeVale(k,j,1),double)
                  packedValeVale(2,currIndex) = aimag(valeVale(k,j,1))
               enddo
            enddo

            call h5dwrite_f (valeValeBand_did,H5T_NATIVE_DOUBLE,&
                  & packedValeVale(:,:),valeValeBand_dims,hdferr)
            if (hdferr == -1) stop 'Failed to write packed vale vale band'
            deallocate (packedValeVale)
         endif
      endif
   else
      if (coreDim /= 0) then
         ! Modify the overlap with the valeCoreCoreValeOL matrix multiplication.
         call valeCoreCoreVale (valeDim,coreDim,valeVale(:,:,1),&
               & coreVale(:,:,1),coreValeOL(:,:,1))

         ! Deallocate the coreVale matrix and allocate the valeCore matrix.
         deallocate (coreVale)
         allocate (valeCore (coreDim,valeDim)) ! Pre-transposed format.

         ! Form a product of (coreValeOL)(coreCore) in a temp matrix (valeCore).
         call coreValeCoreCore (valeDim,coreDim,valeCore,coreValeOL(:,:,1),&
               & coreCore(:,:,1))

         ! Finally compute the product of the above (valeCore)(coreCore) with
         !   coreValeOL.  This is added to valeVale to complete
         !   orthogonaliziation.  Note that the final 1 is used to indicate
         !   that the full matrix should be created for faster cache access
         !   later on.  This is only necessary for the momentum matrix elements.
         allocate (packedValeVale(2,1)) ! Unused but needed for completeness.
         if (runCode == 2) then
            call makeValeVale (valeDim,coreDim,1,valeCore,coreValeOL(:,:,1),&
                  & valeVale(:,:,1),packedValeVale,0,0)
         else
            call makeValeVale (valeDim,coreDim,1,valeCore,coreValeOL(:,:,1),&
                  & valeVale(:,:,1),packedValeVale,0,1)
         endif
         deallocate (packedValeVale)

         ! Deallocate temporary matrices.
         deallocate (valeCore)
         deallocate (coreCore)
      else  ! The valeVale is already made, no need to orthogonalize.
         ! However, the momentum matrix needs to be made into a full matrix,
         !   not just the upper triangle.
         if (runCode /= 2) then
            do i = 1, valeDim
               do j = 1, i
                  valeVale(i,j,1) = conjg(valeVale(j,i,1))
               enddo
            enddo
         endif

         deallocate (coreVale)
         deallocate (coreCore)
      endif
   endif
#else
   if (runCode .eq. 1) then
      if (coreDim /= 0) then
         ! Obtain the orthogonalization coefficients.  (coreValeOL)
         call valeCoreCoreValeOLGamma (valeDim,coreDim,valeValeGamma,&
               & coreValeOLGamma)

         ! Allocate the valeCore matrix in a transposed format for much more
         !   efficient cache usage during the matrix multiplication.
         allocate (valeCoreGamma (coreDim,valeDim))
   
         ! Form a product of (coreValeOL)(coreCore) in a temp matrix
         !   (valeCore).
         call coreValeCoreCoreGamma (valeDim,coreDim,valeCoreGamma,&
               & coreValeOLGamma,coreCoreGamma)

         ! Deallocate the coreCore and allocate the packedValeVale.
         deallocate (coreCoreGamma)
         allocate (packedValeVale(1,valeDim*(valeDim+1)/2)) ! Deallocated
               ! after the getIntgResults subroutine is called from band.f90.
               ! This portion is only run from calls from that program.
   
         ! Finally compute the product of the above (valeCore)(coreCore) with
         !   coreValeOL.  This is added to valeVale to complete
         !   orthogonaliziation.
         call makeValeValeGamma (valeDim,coreDim,valeDim,valeCoreGamma,&
               & coreValeOLGamma,valeValeGamma,packedValeVale,1,0)

         ! Save the overlap valeVale
         call h5dwrite_f (valeValeBand_did,H5T_NATIVE_DOUBLE,&
               & packedValeVale(:,:),valeValeBand_dims,hdferr)
         if (hdferr == -1) stop 'Failed to write vale vale band overlap'

         ! Deallocate temporary matrices to be ready for the next kpoint.
         deallocate (packedValeVale)
         deallocate (valeCoreGamma)
      else  ! No core dimension in the whole system.
         deallocate (coreCoreGamma)
         if (noSaveValeVale == 0) then
            allocate (packedValeVale(1,valeDim*(valeDim+1)/2))

            ! Initialize the index counter for packing.
            currIndex = 0

            ! Pack the valeVale matrix.
            do j = 1, valeDim
               do k = 1, j
                  currIndex = currIndex + 1
                  packedValeVale(1,currIndex) = valeValeGamma(k,j)
               enddo
            enddo

            call h5dwrite_f (valeValeBand_did,H5T_NATIVE_DOUBLE,&
                  & packedValeVale(:,:),valeValeBand_dims,hdferr)
            if (hdferr == -1) stop 'Failed to write packed vale vale band'
            deallocate (packedValeVale)
         endif
      endif
   else
      if (coreDim /= 0) then
         call valeCoreCoreValeGamma (valeDim,coreDim,valeValeGamma,&
               & coreValeGamma,coreValeOLGamma)

         ! Deallocate the coreVale matrix and allocate the valeCore matrix
         !   in a transposed format for much more efficient cache usage
         !   during the matrix multiplication.
         deallocate (coreValeGamma)
         allocate (valeCoreGamma (coreDim,valeDim))

         ! Form a product of (coreValeOL)(coreCore) in a temp matrix
         !   (valeCore).
         call coreValeCoreCoreGamma (valeDim,coreDim,valeCoreGamma,&
               & coreValeOLGamma,coreCoreGamma)

         ! Finally compute the product of the above (valeCore)(coreCore) with
         !   coreValeOL.  This is added to valeVale to complete
         !   orthogonaliziation.  Note that the final 1 is used to indicate
         !   that the full momentum matrix should be stored instead of just
         !   the upper triangle for the purpose of faster cache access later
         !   on.
         allocate (packedValeVale(2,1)) ! Unused, but necessary parameter.
         if (runCode == 2) then
            call makeValeValeGamma (valeDim,coreDim,1,valeCoreGamma,&
                  & coreValeOLGamma,valeValeGamma,packedValeVale,0,0)
         else
            call makeValeValeGamma (valeDim,coreDim,1,valeCoreGamma,&
                  & coreValeOLGamma,valeValeGamma,packedValeVale,0,1)
         endif
         deallocate (packedValeVale)

         ! Deallocate temporary matrices.
         deallocate (valeCoreGamma)
         deallocate (coreCoreGamma)
      else  ! The valeVale is already made, no need to orthogonalize.
         ! However, the momentum matrix needs to be made into a full matrix,
         !   not just the upper triangle.
         if (runCode /= 2) then
            do i = 1, valeDim
               do j = 1, i
                  valeValeGamma(i,j) = valeValeGamma(j,i)
               enddo
            enddo
         endif

         deallocate (coreValeGamma)
         deallocate (coreCoreGamma)
      endif
   endif
#endif

end subroutine getIntgResults


end module O_IntegralsPSCF
