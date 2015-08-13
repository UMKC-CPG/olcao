module O_Bond

   ! Import the necessary modules.
   use O_Kinds

   ! Make sure that no variables are declared accidentally.
   implicit none

   contains

subroutine computeBond

   ! Import necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Potential,   only: spin
   use O_CommandLine, only: doBond
   use O_Lattice,     only: realVectors
   use O_AtomicTypes, only: numAtomTypes, atomTypes
   use O_KPoints,     only: numKPoints, kPointWeight
   use O_AtomicSites, only: valeDim, numAtomSites, atomSites
   use O_Constants,   only: pi, hartree, bigThresh, smallThresh
   use O_Input,       only: numStates, sigmaBOND, eminBOND, emaxBOND, &
         & deltaBOND, maxBL, detailCodeBond, excitedAtomPACS, maxNumNeighbors
#ifndef GAMMA
   use O_SecularEquation, only: valeVale, valeValeOL, readDataSCF, &
         & readDataPSCF, energyEigenValues
#else
   use O_SecularEquation, only: valeValeGamma, valeValeOLGamma, readDataSCF, &
         & readDataPSCF, energyEigenValues
#endif

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables.
   integer :: h,i,j,k,l,m,n,o ! Loop index variables
   integer, allocatable, dimension (:) :: numTypeBonds
   integer, allocatable, dimension (:) :: numChargedAtoms
   integer, allocatable, dimension (:) :: numOrbIndex
   integer, allocatable, dimension (:) :: bondIndex
   integer, allocatable, dimension (:) :: numAtomStates
   integer, allocatable, dimension (:) :: numAtomsBonded
   integer, allocatable, dimension (:) :: bondedAtomID
   integer, allocatable, dimension (:) :: bondedTypeID
   integer, allocatable, dimension (:) :: currentBonds
   integer, allocatable, dimension (:) :: numBondAngles
   integer, allocatable, dimension (:,:,:) :: bondAngleAtoms
   integer, dimension (2) :: atom1Index
   integer, dimension (2) :: atom2Index
   integer :: numSOrbitals
   integer :: numPOrbitals
   integer :: numDOrbitals
   integer :: numFOrbitals
   integer :: maxNumBondAngles
   integer :: maxNumBonds
   integer :: currentNumBonds
   integer :: currentType
   integer :: valeDimIndex
   integer :: numSystemBonds
   integer :: minDistCount
   integer :: bondedAtomCounter
   integer :: numEnergyPoints
   real (kind=double) :: oneValeRealAccum
   real (kind=double) :: systemCharge
   real (kind=double) :: systemBondOrder
   real (kind=double) :: minDistMag
   real (kind=double) :: currentDistMag
   real (kind=double) :: expAlpha
   real (kind=double) :: expFactor
   real (kind=double) :: sigmaSqrtPi
   real (kind=double), dimension (3) :: currentDistance
   real (kind=double), allocatable, dimension (:)      :: atomCharge
   real (kind=double), allocatable, dimension (:)      :: atomChargeTotal
   real (kind=double), allocatable, dimension (:)      :: bondOrderTotal
   real (kind=double), allocatable, dimension (:)      :: bondLengthTotal
   real (kind=double), allocatable, dimension (:)      :: energyScale
   real (kind=double), allocatable, dimension (:)      :: bondOrderEnergyDep !
         ! Energy dep. bond order for the *target atom* in a PACS calculations.
         ! Index range is maxNumNeighbors.
   real (kind=double), allocatable, dimension (:)      :: bondOrderAllEnergy
   real (kind=double), allocatable, dimension (:)      :: totalCalcBondOrder
   real (kind=double), allocatable, dimension (:,:)    :: bondLength
   real (kind=double), allocatable, dimension (:,:)    :: bondOrder
   real (kind=double), allocatable, dimension (:,:)    :: bondAngle
   real (kind=double), allocatable, dimension (:,:)    :: bondedDist
   real (kind=double), allocatable, dimension (:,:)    :: bondCompleteAtom
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:) :: waveFnSqrd
#else
   real (kind=double), allocatable, dimension (:)    :: waveFnSqrdGamma
#endif

   ! Log the date and time we start.
   call timeStampStart (22)


   ! Compute local factors for broadening the energy dependent bond order.
   sigmaSqrtPi = sqrt(pi) * sigmaBOND

   ! Allocate arrays and matrices for this computation.
   allocate (numOrbIndex     (numAtomTypes + 1))
   allocate (bondIndex       (valeDim))
   allocate (numAtomStates   (numAtomSites))
#ifndef GAMMA
   allocate (waveFnSqrd (valeDim))
   if (doBond .eq. 0) then
      allocate (valeVale   (valeDim,numStates,1,1))
      allocate (valeValeOL (valeDim,valeDim,1,1))
   endif
#else
   allocate (waveFnSqrdGamma (valeDim))
   if (doBond .eq. 0) then
      allocate (valeValeGamma   (valeDim,numStates,1))
      allocate (valeValeOLGamma (valeDim,valeDim,1))
   endif
#endif


   ! Initialize counter to index the cumulative sum of orbitals of all types.
   numOrbIndex(1) = 0

   ! Loop to record the index number for each type's orbitals.
   do i = 1, numAtomTypes
      numOrbIndex (i+1) = numOrbIndex(i) + &
            & sum(atomTypes(i)%numQN_lValeRadialFns(:))
   enddo


   ! Initialize the index number in valeDim that each state is at.
   valeDimIndex = 0

   ! Loop over every atom in the system to index where the values for each
   !   atom should be stored.
   do i = 1, numAtomSites

      ! Obtain the type of the current atom.
      currentType = atomSites(i)%atomTypeAssn

      ! Identify and store the number of valence states for this atom.
      numAtomStates(i) = atomTypes(currentType)%numValeStates

      numSOrbitals = atomTypes(currentType)%numQN_lValeRadialFns(1)
      numPOrbitals = atomTypes(currentType)%numQN_lValeRadialFns(2)
      numDOrbitals = atomTypes(currentType)%numQN_lValeRadialFns(3)
      numFOrbitals = atomTypes(currentType)%numQN_lValeRadialFns(4)

      do j = 1, numSOrbitals
         valeDimIndex = valeDimIndex + 1
         bondIndex(valeDimIndex) = numOrbIndex(currentType) + j
      enddo
      do j = 1, numPOrbitals
         do k = 1,3
            valeDimIndex = valeDimIndex + 1
            bondIndex(valeDimIndex) = numOrbIndex(currentType) + numSOrbitals +j
         enddo
      enddo
      do j = 1, numDOrbitals
         do k = 1,5
            valeDimIndex = valeDimIndex + 1
            bondIndex(valeDimIndex) = numOrbIndex(currentType) + numSOrbitals +&
                  & numPOrbitals + j
         enddo
      enddo
      do j = 1, numFOrbitals
         do k = 1,5
            valeDimIndex = valeDimIndex + 1
            bondIndex(valeDimIndex) = numOrbIndex(currentType) + numSOrbitals +&
                  & numPOrbitals + numDOrbitals + j
         enddo
      enddo
   enddo

   ! Determine the number of energy points to be computed for in the case of
   !   the energy dependent bond order calculation.
   numEnergyPoints = (emaxBOND - eminBOND ) / deltaBOND

   ! Allocate space to hold the bondorder and charge results
   allocate (bondOrderEnergyDep (maxNumNeighbors))
   allocate (bondOrder          (numAtomSites,numAtomSites))
   allocate (bondLength         (numAtomSites,numAtomSites))
   allocate (atomCharge         (numAtomSites))
   if (excitedAtomPACS .ne. 0) then
      allocate (energyScale         (numEnergyPoints))
      allocate (bondCompleteAtom    (numEnergyPoints,maxNumNeighbors))
      allocate (bondOrderAllEnergy  (maxNumNeighbors))
      allocate (bondedAtomID        (maxNumNeighbors))
      allocate (bondedTypeID        (maxNumNeighbors))
      allocate (totalCalcBondOrder  (maxNumNeighbors))
   endif


   ! Assign values to the energy scale if necessary.
   if (excitedAtomPACS .ne. 0) then
      do i = 1, numEnergyPoints
         energyScale(i) = eminBOND + (i-1) * deltaBOND
      enddo
   endif

   do h = 1, spin

      ! Initialize various arrays and matrices.
      bondOrder          (:,:) = 0.0_double
      bondLength         (:,:) = 0.0_double
      atomCharge         (:)   = 0.0_double
      if (excitedAtomPACS .ne. 0) then
         bondCompleteAtom    (:,:) = 0.0_double
         bondOrderAllEnergy  (:)   = 0.0_double
         bondedAtomCounter       = 0
         bondedAtomID        (:) = 0
         bondedTypeID        (:) = 0
      endif

      ! Begin accumulating the BOND values over each kpoint
      do i = 1, numKPoints


         ! Determine if we are doing bond+Q* in a post-SCF calculation, or
         !   within an SCF calculation.  (The doBond flag is only set to
         !   true within an SCF calculation because that is the only place that
         !   this option makes sense to use.  In a post-SCF case, the bond job
         !   was explicitly asked for by the user and so this request flag is
         !   not used.)
         if (doBond == 1) then
            ! Read necessary data from SCF (setup,main) data structures.
            call readDataSCF(h,i,numStates)
         else
            ! Read necessary data from post SCF (intg,band) data structures.
            call readDataPSCF(h,i,numStates)
         endif


         do j = 1, numStates

            ! Initialize the bond order for energy dependency.
            bondOrderEnergyDep (:) = 0.0_double

            ! Begin with the energy dependent computation including all energy
            !   states.  Note that only the chosen excited atom is considered.
            !   If the excited atom in the bond order input is set to zero, then
            !   this part of the computation does not proceed.

            if (excitedAtomPACS .ne. 0) then

               ! Initialize the indices of the first atom. (1) = init, (2) = fin
               atom1Index(1) = 0
               atom1Index(2) = 0

               ! Loop over all atoms up to the excited atom to determine the
               !   atom1Index of the excited atom.
               do k = 1, excitedAtomPACS

                  ! Identify the indices of the current atom from loop (k).
                  atom1Index(1) = atom1Index(2) + 1
                  atom1Index(2) = atom1Index(1) + numAtomStates(k) - 1
               enddo

               ! Initialize the indices of the second atom. 1 = init, 2 = fin
               atom2Index(1) = 0
               atom2Index(2) = 0

               ! Initialize the index for the number of bonded atoms.
               bondedAtomCounter = 0

               ! Begin a loop over all the other atoms.
               do l = 1, numAtomSites

                  ! Identify the indices of the current atom from loop (l).
                  atom2Index(1) = atom2Index(2) + 1
                  atom2Index(2) = atom2Index(1) + numAtomStates(l) - 1


                  ! Determine the smallest distance between the atoms of the
                  !   two loops (k) and (l).

                  ! Initialize the min distance to an impossibly large number.
                  minDistMag = bigThresh

                  ! Initialize a counter for the number of atoms that are
                  !   within the max bond length parameter.
                  minDistCount = 0

                  ! Search the original cell and neighboring 124 cells to find
                  !   the minimum distance between atom (k) and atom (l).
                  do m = -2,2
                  do n = -2,2
                  do o = -2,2

                     ! Compute seperation vector for the atoms in this combo.
                     currentDistance(:) = &
                           & atomSites(excitedAtomPACS)%cartPos(:) - &
                           & atomSites(l)%cartPos(:) - &
                           & m * realVectors(:,1) - &
                           & n * realVectors(:,2) - &
                           & o * realVectors(:,3)

                     ! Compute the distance for the atoms of this cell combo.
                     currentDistMag = sqrt(sum(currentDistance(:)**2))

                     ! Compare the current magnitude to the current minimum.
                     if (currentDistMag < minDistMag .and. &
                           & abs(currentDistMag) > smallThresh) then
                        minDistMag = currentDistMag
                     endif

                     ! Check if this currentDistMag is less than the cut-off
                     !   radius.  If so, then increment the counter for this
                     !   atom pair.  This is done in case a replicated atom has
                     !   more than one configuration where it is sufficiently
                     !   close to the current target atom that it regesters an
                     !   overlap.  This effect is not really implemented yet
                     !   because the minDistMag will have to track all the
                     !   nearest atoms, not just the closest one.
                     if (currentDistMag < maxBL .and. &
                           & abs(currentDistMag) > smallThresh) then
                        minDistCount = minDistCount + 1
                     endif
                  enddo
                  enddo
                  enddo

                  ! Determine if any image atoms are within the cut-off radius.
                  if (minDistCount > 0) then

                     ! Track indices and IDs for the atoms that bond to the
                     !    excited atom.

                     ! Increment the number of atoms that are bonded with the
                     !   excited atom.
                     bondedAtomCounter = bondedAtomCounter + 1

                     ! Record the atom ID number for this atom that bonds with
                     !   the excited atom.
                     bondedAtomID(bondedAtomCounter) = l

                     ! Store the determined minimum bond length for this pair.
                     bondLength(excitedAtomPACS,l) = minDistMag

                     ! Loop over the atom 1 states against the atom 2 states
                     !   for this band.
                     do m = atom1Index(1),atom1Index(2)

#ifndef GAMMA
                        ! Compute ^2 of the wave function for each element.
                        waveFnSqrd(atom2Index(1):atom2Index(2)) = &
                              & conjg(valeVale(m,j,1,1)) * &
                              & valeVale(atom2Index(1):atom2Index(2),j,1,1)

                        ! Compute the effects of overlap for the real part
                        !   only (real*real) + (imag*imag).
                        oneValeRealAccum = sum(&
                              & real (waveFnSqrd(atom2Index(1):&
                              & atom2Index(2)),double) * &
                              & real (valeValeOL(atom2Index(1):&
                              & atom2Index(2),m,1,1),double) + &
                              & aimag(waveFnSqrd(atom2Index(1):&
                              & atom2Index(2)))*&
                              & aimag(valeValeOL(atom2Index(1):&
                              & atom2Index(2),m,1,1)))
#else
                        ! Compute ^2 of the wave function for each element.
                        waveFnSqrdGamma(atom2Index(1):atom2Index(2)) = &
                              & valeValeGamma(m,j,1) * &
                              & valeValeGamma(atom2Index(1):atom2Index(2),&
                              & j,1)

                        ! Compute the effects of overlap for the real part
                        !   only (real*real) + (imag*imag).
                        oneValeRealAccum = sum(&
                              & waveFnSqrdGamma(atom2Index(1):&
                              & atom2Index(2)) * &
                              & valeValeOLGamma(atom2Index(1):&
                              & atom2Index(2),m,1))
#endif

                        ! Store the atom pair bond order contribution for this
                        !   energy level.
                        bondOrderEnergyDep(bondedAtomCounter) = &
                              & bondOrderEnergyDep(bondedAtomCounter) + &
                              & oneValeRealAccum * kPointWeight(i) / &
                              & real(spin,double)

                        ! Store the atom pair bond order contribution.
                        bondOrderAllEnergy(bondedAtomCounter) = &
                              & bondOrderAllEnergy(bondedAtomCounter) + &
                              & oneValeRealAccum * kPointWeight(i) / &
                              & real(spin,double)

                     enddo
                  endif
               enddo ! (l atom2)
               do l = 1, numEnergyPoints

                  ! Compute the exponential alpha for broadening this point.
                  expAlpha = ((energyEigenValues(j,i,h) - &
                        & energyScale(l)) / sigmaBOND)**2

                  ! If the exponential alpha is less than 50 apply broadening.
                  if (expAlpha < 50) then

                     ! Compute the exponential factor.
                     expFactor = exp(-expAlpha) / sigmaSqrtPi

                     ! Broaden and store the complete energy dep. bond order.
                     bondCompleteAtom(l,:bondedAtomCounter) = &
                           & bondCompleteAtom(l,:bondedAtomCounter) + &
                           & bondOrderEnergyDep(:bondedAtomCounter) * expFactor
                  endif
               enddo
            endif ! (excited atom .ne. 0)

            ! Skip energy states that are not occupied.
            if (energyEigenValues(j,i,h) > 0.0_double) exit

            ! Initialize the indices of the first atom. (1) = init, (2) = fin
            atom1Index(1) = 0
            atom1Index(2) = 0

            ! Begin the first loop over all atoms
            do k = 1, numAtomSites

               ! Identify the indices of the current atom from loop (k).
               atom1Index(1) = atom1Index(2) + 1
               atom1Index(2) = atom1Index(1) + numAtomStates(k) - 1

               ! Loop over the atom 1 indexed states against all other states
               !   in this band (j).
               do l = atom1Index(1), atom1Index(2)

#ifndef GAMMA
                  ! Compute the square of the wave function for each element.
                  waveFnSqrd(:valeDim) = conjg(valeVale(l,j,1,1)) * &
                        & valeVale(:valeDim,j,1,1)

                  ! Compute the effects of overlap for the real part only
                  !   (real*real) + (imag*imag).
                  oneValeRealAccum = sum( &
                        & real(waveFnSqrd(:valeDim),double) * &
                        & real (valeValeOL(:valeDim,l,1,1),double) + &
                        & aimag(waveFnSqrd(:valeDim)) * &
                        & aimag(valeValeOL(:valeDim,l,1,1)))
#else
                  ! Compute the square of the wave function for each element.
                  waveFnSqrdGamma(:valeDim) = valeValeGamma(l,j,1) * &
                        & valeValeGamma(:valeDim,j,1)

                  ! Compute the effects of overlap for the real part only
                  !   (real*real) + (imag*imag).
                  oneValeRealAccum = sum(waveFnSqrdGamma(:valeDim) * &
                        & valeValeOLGamma(:valeDim,l,1))
#endif

                  ! Store the atom charge contribution from this band (j) and
                  !   kpoint (i).
                  atomCharge(k) = atomCharge(k) + oneValeRealAccum * &
                        & kPointWeight(i) / real(spin,double)
               enddo

               ! Initialize the indices of the second atom. (1)=init, (2)=fin
               atom2Index(1) = 0
               atom2Index(2) = 0

               ! Begin the second loop over all atoms
               do l = 1, numAtomSites

                  ! Identify the indices of the current atom from loop (l).
                  atom2Index(1) = atom2Index(2) + 1
                  atom2Index(2) = atom2Index(1) + numAtomStates(l) - 1


                  ! Determine the smallest distance between the atoms of the
                  !   two loops (k) and (l).

                  ! Initialize the min distance to an impossibly large number.
                  minDistMag = bigThresh

                  ! Initialize a counter for the number of atoms that are
                  !   within the max bond length parameter.
                  minDistCount = 0

                  ! Search the original cell and the neighboring 124 cells to
                  !   find the minimum distance between atom (k) and atom (l).
                  do m = -2,2
                  do n = -2,2
                  do o = -2,2

                     ! Compute the seperation vector for the atoms in this
                     !   combo.
                     currentDistance(:) = atomSites(k)%cartPos(:) - &
                           & atomSites(l)%cartPos(:) - &
                           & m * realVectors(:,1) - &
                           & n * realVectors(:,2) - &
                           & o * realVectors(:,3)

                     ! Compute the distance for the atoms of this cell combo.
                     currentDistMag = sqrt(sum(currentDistance(:)**2))

                     ! Compare the current magnitude to the current minimum.
                     if (currentDistMag < minDistMag .and. &
                           & abs(currentDistMag) > smallThresh) then
                        minDistMag = currentDistMag
                     endif

                     ! Check if this minDistMag is less than the cut-off radius.
                     !   radius.  If so, then increment the counter for this
                     !   atom pair.  This is done in case a replicant atom has
                     !   more than one configuration where it is sufficiently
                     !   close to the current target atom that it regesters an
                     !   overlap.  This effect is not really implemented yet
                     !   since the minDistMag will have to track all the
                     !   nearest atoms, not just the closest one.
                     if (currentDistMag < maxBL .and. &
                           & abs(currentDistMag) > smallThresh) then
                        minDistCount = minDistCount + 1
                     endif
                  enddo
                  enddo
                  enddo

                  ! Determine if any image atoms are within the cut-off radius.
                  if (minDistCount > 0) then

                     ! Store the determined minimum bond length for this pair.
                     bondLength(k,l) = minDistMag

                     ! Loop over the atom 1 states against the atom 2 states
                     !   for this band.
                     do m = atom1Index(1),atom1Index(2)

#ifndef GAMMA
                        ! Compute ^2 of the wave function for each element.
                        waveFnSqrd(atom2Index(1):atom2Index(2)) = &
                              & conjg(valeVale(m,j,1,1)) * &
                              & valeVale(atom2Index(1):atom2Index(2),j,1,1)

                        ! Compute the effects of overlap for the real part
                        !   only (real*real) + (imag*imag).
                        oneValeRealAccum = sum(&
                              & real (waveFnSqrd(atom2Index(1):&
                              &   atom2Index(2)),double)*&
                              & real (valeValeOL(atom2Index(1):&
                              &   atom2Index(2),m,1,1),double) +&
                              & aimag(waveFnSqrd(atom2Index(1):&
                              &   atom2Index(2)))*&
                              & aimag(valeValeOL(atom2Index(1):&
                              &   atom2Index(2),m,1,1)))
#else
                        ! Compute ^2 of the wave function for each element.
                        waveFnSqrdGamma(atom2Index(1):atom2Index(2)) = &
                              & valeValeGamma(m,j,1) * &
                              & valeValeGamma(atom2Index(1):atom2Index(2),j,1)

                        ! Compute the effects of overlap.
                        oneValeRealAccum = sum(&
                              & waveFnSqrdGamma(atom2Index(1):&
                              & atom2Index(2)) * &
                              & valeValeOLGamma(atom2Index(1):&
                              & atom2Index(2),m,1))
#endif

                        ! Store the atom pair bond order contribution.
                        bondOrder(k,l) = bondOrder(k,l) + oneValeRealAccum * &
                              & kPointWeight(i) / real(spin,double)

                     enddo
                  endif
               enddo ! (l atom2)
            enddo ! (k atom1)
         enddo ! (j states)

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

      enddo ! (i kpoints)



      ! Begin accumulating statistics and recording output to disk.



      ! Compute the maximum number of bonds present for an atom in the system.
      !   Note that we consider only those bonds with a positive bond order.
      !   (i.e. no anti-bonding considered).
      maxNumBonds = 0
      do i = 1, numAtomSites

         ! Initialize a counter for the number of bonds for this (i) atom.
         currentNumBonds = 0

         ! Loop over all other atoms to accumulate how many are bonded to (i).
         do j = 1, numAtomSites
            if (bondOrder(j,i) > 0.0_double) then
               currentNumBonds = currentNumBonds + 1
            endif
         enddo

         ! Compare the number of bonds for atom (i) to the current maximum.
         if (currentNumBonds > maxNumBonds) then
            maxNumBonds = currentNumBonds
         endif
      enddo

      ! From the max number of bonds (that some atom must have), compute the
      !   max number of bond angles (which must apply to that same atom).  The
      !   first bond can be angled with the other N-1.  The second can be
      !   angled with the other N-2, the third with N-3, and so on.  Hence,
      !   the number of bond angles is (N-1) + (N-2) + (N-3) + ... (N-N) =
      !   Sum(0 to N-1) of i = N * (N-1) / 2.
      maxNumBondAngles = maxNumBonds * (maxNumBonds-1) / 2

      ! Allocate space for statistics holding arrays and matrices.
      allocate (bondAngleAtoms     (2,maxNumBondAngles,numAtomSites))
      allocate (bondAngle          (maxNumBondAngles,numAtomSites))
      allocate (currentBonds       (maxNumBonds))
      allocate (numBondAngles      (numAtomSites))
      allocate (numChargedAtoms    (numAtomTypes))
      allocate (numTypeBonds       (numAtomTypes))
      allocate (atomChargeTotal    (numAtomTypes))
      allocate (bondOrderTotal     (numAtomTypes))
      allocate (bondLengthTotal    (numAtomTypes))
      allocate (numAtomsBonded     (numAtomSites))

      ! Compute the bond triplets for the bond angles.  i is always the middle
      !   component of the triplet.
      do i = 1, numAtomSites

         ! Obtain a convenient list of all the atom numbers that this atom
         !   bonds to.
         currentBonds(:) = 0
         currentNumBonds = 0
         do j = 1, numAtomSites
            if (bondOrder(j,i) > 0.0_double) then
               currentNumBonds = currentNumBonds + 1
               currentBonds(currentNumBonds) = j
            endif
         enddo

         ! Find the positions of the atoms that bond to this atom (i) relative
         !   to the position of (i).

         ! Allocate space to hold the positions of each atom.
         allocate (bondedDist (3,currentNumBonds))

         do j = 1, currentNumBonds

            ! Initialize the min distance to an impossibly large number.
            minDistMag = bigThresh

            ! Loop through 125 cells to find the minimal distance
            do k = -2, 2
            do l = -2, 2
            do m = -2, 2

               ! Compute seperation vector for the atoms in this combo.
               currentDistance(:) = atomSites(i)%cartPos(:) - &
                     & atomSites(currentBonds(j))%cartPos(:) - &
                     & k * realVectors(:,1) - &
                     & l * realVectors(:,2) - &
                     & m * realVectors(:,3)

               ! Compute the distance for the atoms of this cell combo.
               currentDistMag = sqrt(sum(currentDistance(:)**2))

               ! Compare the current magnitude to the current minimum.
               if (currentDistMag < minDistMag .and. &
                     & abs(currentDistMag) > smallThresh) then
                  minDistMag = currentDistMag
               endif
            enddo
            enddo
            enddo

            ! Record the minimal distance.
            bondedDist(:,j) = currentDistance(:)
         enddo

         ! Traverse the list of bonding atoms to record the bond triplets and
         !   compute the angles.
         numBondAngles(i) = 0
         do j = 1, currentNumBonds
            do k = j+1, currentNumBonds

               ! Increment the current number of bonds.
               numBondAngles(i) = numBondAngles(i) + 1

               ! Record which atoms are present for this bond angle.
               bondAngleAtoms(1,numBondAngles(i),i) = currentBonds(j)
               bondAngleAtoms(2,numBondAngles(i),i) = currentBonds(k)

               ! Compute the bond angle as arccos(v1.v2 / (|v1| * |v2|)).
               bondAngle(numBondAngles(i),i) = acos(dot_product(&
                     & bondedDist(:,j),bondedDist(:,k)) / &
                     & sqrt(dot_product(bondedDist(:,j),bondedDist(:,j))) / &
                     & sqrt(dot_product(bondedDist(:,k),bondedDist(:,k))))
            enddo
         enddo

         deallocate (bondedDist)
      enddo


      ! Compute the total bond order for the excited atom as an energy
      !   dependent function for each species it is bonded to.
      if (excitedAtomPACS .ne. 0) then

         ! First the total bond order is computed for each bonded atom as an
         !   integration of the bond energy dependent bond order.  Simply
         !   accumulate the area under the curve by the areas of rectangles.
         do i = 1, bondedAtomCounter
            totalCalcBondOrder(i) = sum(&
                  & bondCompleteAtom(1:numEnergyPoints-1,i) + &
                  & bondCompleteAtom(2:numEnergyPoints,i)) * deltaBOND * &
                  & 0.5_double
         enddo

         ! Then the complete bond order for each bonded atom is normalized.
         do i = 1, bondedAtomCounter
            bondCompleteAtom(:,i) = bondCompleteAtom(:,i) * &
                  & bondOrderAllEnergy(i) / totalCalcBondOrder(i)
         enddo
      endif


      ! Initialize counters for totals.
      atomChargeTotal (:) = 0.0_double
      bondOrderTotal  (:) = 0.0_double
      bondLengthTotal (:) = 0.0_double
      numChargedAtoms (:) = 0
      numTypeBonds    (:) = 0
      numAtomsBonded  (:) = 0
      systemBondOrder     = 0.0_double
      numSystemBonds      = 0
      systemCharge        = 0

      ! Begin a loop over all atoms.
      do i = 1, numAtomSites

         ! Determine the type of the current atom.
         currentType = atomSites(i)%atomTypeAssn

         ! Accumulate the charge for this atom to its atom type total.
         atomChargeTotal(currentType) = atomChargeTotal(currentType) + &
               & atomCharge(i)
         numChargedAtoms(currentType) = numChargedAtoms(currentType) + 1

         ! Accumulate the charge for this atom to the total system charge.
         systemCharge = systemCharge + atomCharge(i)

         ! Begin a second loop over the other atoms.
         do j = 1, numAtomSites

            ! Only include bonds that exist in the statistics.
            if (bondLength(i,j) /= 0.0_double) then

               ! Accumulate the bond length and bond order for this atom pair
               !   to the type of the (i) atom.
               bondLengthTotal(currentType) = bondLengthTotal(currentType) + &
                     & bondLength(i,j)
               bondOrderTotal(currentType) = bondOrderTotal(currentType) + &
                     & bondOrder(i,j)
               numTypeBonds(currentType) = numTypeBonds(currentType) + 1

               ! Only consider bonds where i < j to prevent double counting.
               if (i < j) then
                  ! Accumulate BO for this atom pair to the total system bond
                  !   order.
                  systemBondOrder = systemBondOrder + bondOrder(i,j)

                  ! Accumulate the number of bonds in the system.
                  numSystemBonds = numSystemBonds + 1

                  ! Record the number of atoms in the j loop that pair with
                  !   atoms in the i loop.
                  numAtomsBonded(i) = numAtomsBonded(i) + 1
               endif
            endif
         enddo
      enddo

      ! Begin recording the energy dependent bond order results to disk, making
      !   sure to convert the energy scale to eV.
      if (excitedAtomPACS .ne. 0) then
         write (16,ADVANCE="NO",FMT="(a12,a12)") "Energy      ","Total       "
         do i = 1, bondedAtomCounter
            write (16,ADVANCE="NO",FMT="(a2,a10)") "  ",&
                  & atomTypes(i)%typeLabel
         enddo
         write (16,*)
         do i = 1, numEnergyPoints
            write (16,fmt="(30f12.8)") energyScale(i)*hartree,&
                  & sum(bondCompleteAtom(i,:bondedAtomCounter)),&
                  & bondCompleteAtom(i,:bondedAtomCounter)
         enddo
      endif

      ! Begin recording general bond order results to disk.
      if (detailCodeBond == 0) then
         do i = 1, numAtomTypes

            write (9+h,fmt="(2x,a10,i5,a3,a18)") 'Atom type:',i,' - ',&
                  & atomTypes(i)%typeLabel
            write (9+h,fmt="(4x,a28)") 'Individual Effective Charge:'

            ! Loop over all atoms looking for those of the current type.
            do j = 1, numAtomSites

               ! Only record atoms of the current type.
               if (atomSites(j)%atomTypeAssn == i) then
                  write (9+h,fmt="(4x,i5,10x,f10.6)") j,atomCharge(j)
               endif
            enddo

            ! Record text notes for total bond order output in fort.10.
            write (9+h,fmt="(4x,a17,i5)") 'Number of Atoms: ',numChargedAtoms(i)
            write (9+h,fmt="(4x,a16,f12.4)") 'Average Charge: ',&
                  & atomChargeTotal(i)/numChargedAtoms(i)

            write (9+h,*)
            write (9+h,fmt="(4x,a36)") 'Individual Bondlength  &  Bondorder:'

            ! Begin nested loops over all atoms for the bond order.
            do j = 1, numAtomSites

               ! Only record atoms of the current type.
               if (atomSites(j)%atomTypeAssn == i) then
                  do k = 1, numAtomSites

                     ! Only include bonds that exist in the statistics.
                     if (bondLength(j,k) /= 0.0_double) then

                        ! Record the first atom and second atom, the bond
                        !   length, and the bondOrder.
                        write (9+h,fmt="(4x,2i5,5x,f8.4,2x,f8.4)") j,k,&
                              & bondLength(j,k),bondOrder(j,k)
                     endif
                  enddo
               endif
            enddo

            write (9+h,fmt="(4x,a17,i5)") 'Number of Bonds: ',numTypeBonds(i)
            if (numTypeBonds(i) > 0) then
               write (9+h,fmt="(4x,a20,f12.4)") 'Average Bondlength: ',&
                     & bondLengthTotal(i) / numTypeBonds(i)
               write (9+h,fmt="(4x,a19,f12.4)") 'Average BondOrder: ',&
                     & bondOrderTotal(i) / numTypeBonds(i)
            endif
         enddo

         write (9+h,fmt="(4x,a23,i6)") 'Total Number of Bonds: ',numSystemBonds
         write (9+h,fmt="(4x,a17,f12.4)") 'Total Bondorder: ',systemBondOrder
         write (9+h,fmt="(4x,a14,f12.4)") 'Total Charge: ',systemCharge
      else
         ! Print out useful header information.
         write (9+h,fmt=200) "NUM_ATOMS        ",numAtomSites

         do i = 1, numAtomSites
            currentType = atomSites(i)%atomTypeAssn
            write (9+h,fmt=200) "ATOM_NUM         ", i
            write (9+h,fmt=200) "SYSTEM_NUM       ", 1
            write (9+h,fmt=400) "ELEMENT_NAME     ", &
                  & atomTypes(currentType)%elementName
            write (9+h,fmt=200) "ELEMENT_ID       ", &
                  & atomTypes(currentType)%elementID
            write (9+h,fmt=200) "SPECIES_ID       ", &
                  & atomTypes(currentType)%speciesID
            write (9+h,fmt=200) "TYPE_ID          ", &
                  & atomTypes(currentType)%typeID
            write (9+h,fmt=300) "ATOM_CHARGE      ", atomCharge(i)
            write (9+h,fmt=200) "NUM_BONDED_ATOMS ", numAtomsBonded(i)
            do j = i+1, numAtomSites
               ! Only include bonds that exist in the statistics.
               if (bondLength(i,j) /= 0.0_double) then

                  ! Record the bonded atom, the bond length, and the bond order.
                  write (9+h,fmt="(i5,5x,f8.4,2x,f8.4)") j,bondLength(i,j),&
                        & bondOrder(i,j)
               endif
            enddo
            write (9+h,fmt=200) "NUM_BOND_ANGLES  ", numBondAngles(i)
            do j = 1, numBondAngles(i)
               ! Record the triplet of atoms for this bond angle and the bond
               !   angle converted into degrees.
               write (9+h,fmt="(i5,5x,i5,5x,i5,f8.4)"), bondAngleAtoms(1,j,i),&
                     & i,bondAngleAtoms(2,j,i), bondAngle(j,i)*180.0_double/pi
            enddo
                  
         enddo
      endif

      ! Deallocate unnecessary matrices for this spin iteration
      deallocate (numChargedAtoms)
      deallocate (numTypeBonds)
      deallocate (atomChargeTotal)
      deallocate (bondOrderTotal)
      deallocate (bondLengthTotal)
      deallocate (numAtomsBonded)
      deallocate (numBondAngles)
      deallocate (bondAngleAtoms)
      deallocate (bondAngle)
      deallocate (currentBonds)
   enddo ! (h spin)

   ! Deallocate matrices that are no longer needed.
   deallocate (numOrbIndex)
   deallocate (bondIndex)
   deallocate (numAtomStates)
   deallocate (bondOrderEnergyDep)
   deallocate (bondOrder)
   deallocate (bondLength)
   deallocate (atomCharge)
#ifndef GAMMA
   deallocate (waveFnSqrd)
   if (doBond == 0) then
      deallocate (valeVale)
      deallocate (valeValeOL)
   endif
#else
   deallocate (waveFnSqrdGamma)
   if (doBond == 0) then
      deallocate (valeValeGamma)
      deallocate (valeValeOLGamma)
   endif
#endif
   if (excitedAtomPACS /= 0) then
      deallocate (energyScale)
      deallocate (bondCompleteAtom)
      deallocate (bondOrderAllEnergy)
      deallocate (bondedAtomID)
      deallocate (bondedTypeID)
      deallocate (totalCalcBondOrder)
   endif


   ! Log the date and time we end.
   call timeStampEnd (22)


   ! Define all the formatted output styles.
   200 format (a17,i5)
   300 format (a17,f12.8)
   400 format (a17,a3)

end subroutine computeBond

end module O_Bond
