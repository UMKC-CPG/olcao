module O_Bond3C

   ! Import the necessary modules.
   use O_Kinds

   ! Make sure that no variables are declared accidentally.
   implicit none

   ! Define access
   public

   ! Begin list of module derived types.
   type atom3Node
      type(atom3Node), pointer :: next
      integer :: atomSiteNumber
      integer :: valeDimRangeMag
      integer, dimension (2) :: atom3Index ! valeDim beginning=1 ending=2
      integer, dimension (3) :: latticeIndices ! Index values that indicate
            ! which replicated lattice this atom belongs to.
      real (kind=double) :: distToAtom1
      real (kind=double) :: distToAtom2
      real (kind=double), dimension (4) :: bondOrder3C ! Index 1 = BO between
            ! atoms 1,2; Index 2 = BO between atoms 1,3; Index 3 = BO between
            ! atoms 2,3; INdex 4 = Total three center bond order.
      real (kind=double), dimension (3) :: centroid ! Geometric centroid
      real (kind=double), dimension (3) :: bondCentroid ! Centroid weighted by
            ! the magnitudes of the BO between atom pairs.  (This should
            ! identify the centroid of the three center bond.)
      real (kind=double), dimension (3) :: bondCentroidInCell ! This is the
            ! same as the bondCentroid except that (if necessary) it has been
            ! shifted to be inside the unit cell.
   end type atom3Node

   type atom2Node
      type(atom2Node), pointer :: next
      type(atom3Node), pointer :: head
      integer :: atomSiteNumber
      integer :: valeDimRangeMag
      integer, dimension (2) :: atom2Index ! valeDim beginning=1 ending=2
      integer, dimension (3) :: latticeIndices ! Index values that indicate
            ! which replicated lattice this atom belongs to.
      real (kind=double) :: distToAtom1
   end type atom2Node

   type atom1Node
      type(atom2Node), pointer :: head
      integer :: valeDimRangeMag
      integer, dimension (2) :: atom1Index ! valeDim beginning=1 ending=2
   end type atom1Node

   ! Begin list of module data.
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:) :: waveFnSqrd
#else
   real (kind=double), allocatable, dimension (:)    :: waveFnSqrdGamma
#endif

   contains

subroutine computeBond3C

   ! Import necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Potential,   only: spin
   use O_CommandLine, only: doBond
   use O_Input,       only: numStates
   use O_Constants,   only: smallThresh
   use O_Lattice,     only: realVectors
   use O_Populate,    only: electronPopulation
   use O_KPoints,     only: numKPoints, kPointWeight
   use O_AtomicSites, only: valeDim, numAtomSites, atomSites
   use O_AtomicTypes, only: numAtomTypes, atomTypes, maxNumValeStates
#ifndef GAMMA
   use O_SecularEquation, only: valeVale, valeValeOL, readDataSCF, readDataPSCF
#else
   use O_SecularEquation, only: valeValeGamma, valeValeOLGamma, readDataSCF, &
         & readDataPSCF
#endif

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables.
   integer :: h,i,j,k,l,m,n,o ! Loop index variables
   integer, allocatable, dimension (:) :: numTypeBonds
   integer, allocatable, dimension (:) :: numOrbIndex
   integer, allocatable, dimension (:) :: bondIndex
   integer, allocatable, dimension (:) :: numAtomStates
   integer, allocatable, dimension (:) :: numAtomsBonded
   integer, allocatable, dimension (:) :: currentBonds
   integer, allocatable, dimension (:) :: numBondAngles
   integer, allocatable, dimension (:,:,:) :: bondAngleAtoms
   integer, dimension (2) :: atom1Index
   integer, dimension (2) :: atom2Index
   integer, dimension (2) :: atom3Index
   integer, dimension (3,3) :: latticeIndices
   integer :: num3CB ! A count of the number of 3C bonds in the system.
   integer :: orderedIndex ! An index map between a linear list of states and
         ! the spin-kpoint-state loop variables.
   integer :: numSOrbitals
   integer :: numPOrbitals
   integer :: numDOrbitals
   integer :: numFOrbitals
   integer :: currentType
   integer :: valeDimIndex
   integer :: minDistCount
   real (kind=double) :: systemBondOrder
   real (kind=double) :: currentDistMag
   real (kind=double) :: expAlpha
   real (kind=double) :: expFactor
   real (kind=double) :: chargeScaleFactor ! A # from 0-1 that scales
         ! contributions to the bondOrder according to the electron population
         ! in the current state.
   real (kind=double), dimension (2) :: currentMinDist
   real (kind=double), allocatable, dimension (:)      :: bondOrderTotal
   real (kind=double), allocatable, dimension (:)      :: bondLengthTotal
   real (kind=double), allocatable, dimension (:,:)    :: bondAngle
   real (kind=double), allocatable, dimension (:,:)    :: bondedDist
   type (atom1Node), allocatable, dimension (:) :: BO3C
   type (atom2Node), pointer :: currentNode2
   type (atom3Node), pointer :: currentNode3


   ! Log the date and time we start.
   call timeStampStart (26)


   ! Allocate arrays and matrices for this computation.
   allocate (BO3C            (numAtomSites))
   allocate (numOrbIndex     (numAtomTypes + 1))
   allocate (bondIndex       (valeDim))
   allocate (numAtomStates   (numAtomSites))
#ifndef GAMMA
   allocate (waveFnSqrd (maxNumValeStates))
   if (doBond .eq. 0) then
      allocate (valeVale   (valeDim,numStates,1,1))
      allocate (valeValeOL (valeDim,valeDim,1,1))
   endif
#else
   allocate (waveFnSqrdGamma (maxNumValeStates))
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

   ! Create the bond data structure and populate with all information that will
   !   not change for any given spin-kpoint-state combination (e.g. bond
   !   distances and valence dimension indices for the bonded atoms).
   ! Note that there is a bit of tricky business regarding the record of the
   !   valence dimension indices of the bonded atoms. In executing this
   !   algorithm we must know the valeDim range of whatever atom we are on in
   !   the k,l,m loops, but we only want to record the range for atoms that
   !   are actually bonded in the linked list. Therefore, we will keep track
   !   of the ranges with atom1Index(:) types of variables and we will record
   !   into ~~~%atom1Index(:) variables inside the linked list.  That is, there
   !   will be two sets of the same kind of varaible appearing in the algorithm.
   !   I will try not to be confusing about it by using "Identify" and
   !   "Initialize" in the comments relating to the first mission and "Store"
   !   in the comments related to the second.
   ! Note also that this three center bond order calculation will be operating
   !   under the assumption that a given atom will have only one bond with any
   !   other atom *including* any replicated copies of that other atom.  So,
   !   consider that atom 1 is bonded to atom 2 and that both are inside the
   !   unit cell. (This is not a requirement, but it is just for example. We
   !   could also consider atom 1 as being inside the unit cell and atom 2 as
   !   being in a replicated periodic cell. The same logic will apply in that
   !   case.) We assume that there are no replicated copies of atom 2 that
   !   are in neighboring cells that atom 1 will be close enough for atom 1 to
   !   bond with. FURTHER, we assume that the final atom of the three center
   !   bond (atom 3) will follow the same rule regarding its relation to atom 1.
   !   FINALLY, we assume that the relationship between atom 2 and atom 3 will
   !   follow a similar rule such that *the* atom 2 that atom 3 is bonded to
   !   will be the *same* atom 2 that atom 1 is bonded to (meaning it will not
   !   be a different replicated copy).  We will compute the lattice indices
   !   for atom2 and atom3 from the prospective that atom1 will always be in
   !   the central unit cell.
   ! It is a bit tricky to check the appropriate indices for which replicated
   !   cell an atom is in to be sure that the above assumption is not violated.
   !   This will have to be done in the future, but it isn't now. The problem
   !   is that the lattice indices for atom2 and atom3 are obtained assuming
   !   that atom1 is in the central unit cell. Therefore if we want to check
   !   the atom2 to atom3 relationship we might have the case that the atom2
   !   that bonds to atom1 is in a given cell (say 1,2,2) but then atom3 is
   !   also in another cell (say 2,2,2) so that if we looked for the cell
   !   indices of atom2 with atom3 being in the central unit cell we would
   !   find cell (1,-1,-1) or something. (The values of the numbers are
   !   meaningless and just given as mathematically *incorrect* examples.)
   ! Essentially, this whole three center bond order calculation will not work
   !   if the cell is small such that atoms bond with multiple copies of other
   !   atoms. You should expect broken and uniterpretable results if you
   !   attempt such a small cell calculation.

   ! Initialize the valeDim indices of the first atom. (1)=init, (2)=fin
   atom1Index(:) = 0

   ! Initialize the head pointers in the first atom array.
   do k = 1, numAtomSites
      BO3C(k)%head => null()
   enddo

   ! Begin the main algorithm for setting up the list of bonded atoms. There
   !   are three loops over numAtomSites (k,l,m) that will probe every possible
   !   grouping of three atoms. The sets that are all mutually within range of
   !   each other will be recorded and will later have their three center bond
   !   order computed.
   do k = 1, numAtomSites

      ! Identify the valeDim indices of the current atom from loop (k).
      atom1Index(1) = atom1Index(2) + 1
      atom1Index(2) = atom1Index(1) + numAtomStates(k) - 1

      ! Store the valeDim indices of the current atom from loop (k).
      BO3C(k)%atom1Index(1) = atom1Index(1)
      BO3C(k)%atom1Index(2) = atom1Index(2)

      ! Store the magnitude of the valeDim range.
      BO3C(k)%valeDimRangeMag = atom1Index(2) - atom1Index(1) + 1

      ! Initialize the valeDim indices of the *second* atom. Because the loop
      !   indices will not double count, we have to count up to the current
      !   index for atom2.
      atom2Index(:) = 0
      do l = 1, k
         atom2Index(1) = atom2Index(2) + 1
         atom2Index(2) = atom2Index(1) + numAtomStates(l) - 1
      enddo

      do l = k+1, numAtomSites

         ! Identify the indices of the current atom from loop (l).
         atom2Index(1) = atom2Index(2) + 1
         atom2Index(2) = atom2Index(1) + numAtomStates(l) - 1

         ! Compute distance from atom1 to atom2 and cycle if out of range.
         call getMinDist(k,l,minDistCount,latticeIndices(:,1),&
               & currentMinDist(1))
         if (minDistCount == 0) cycle

         ! We found an atom in range of atom1. Allocate space to store
         !   information about it and then point to it.
         if (.not.associated(BO3C(k)%head)) then
            allocate(BO3C(k)%head)
            currentNode2 => BO3C(k)%head
         else
            allocate(currentNode2%next)
            currentNode2 => currentNode2%next
         endif

         ! Store the valeDim indices of the current atom from the loop (l).
         currentNode2%atom2Index(1) = atom2Index(1)
         currentNode2%atom2Index(2) = atom2Index(2)

         ! Store the magnitude of the valeDim range.
         currentNode2%valeDimRangeMag = atom2Index(2) - atom2Index(1) + 1

         ! Store the replicated lattice cell indices for atom2.
         currentNode2%latticeIndices(:) = latticeIndices(:,1)

         ! Store the minimum distance from this atom2 (l) to atom1 (k).
         currentNode2%distToAtom1 = currentMinDist(1)

         ! Store the site number for this atom.
         currentNode2%atomSiteNumber = l

         ! Prepare this node for storing a list of atom3s
         currentNode2%head => null()

         ! Establish that this is the current end of the list for atom2.
         currentNode2%next => null()

         ! Initialize the valeDim indices of the *third* atom.
         atom3Index(:) = 0
         do m = 1, l
            atom3Index(1) = atom3Index(2) + 1
            atom3Index(2) = atom3Index(1) + numAtomStates(m) - 1
         enddo

         do m = l+1, numAtomSites

            ! Identify the valeDim indices of the current atom from loop (m).
            atom3Index(1) = atom3Index(2) + 1
            atom3Index(2) = atom3Index(1) + numAtomStates(m) - 1

            ! Compute distance from atom1 to atom3 and cycle if out of range.
            call getMinDist(k,m,minDistCount,latticeIndices(:,2),&
                  & currentMinDist(1))
            if (minDistCount == 0) cycle

            ! Compute distance from atom2 to atom3 and cycle if out of range.
            call getMinDist(l,m,minDistCount,latticeIndices(:,3),&
                  & currentMinDist(2))
            if (minDistCount == 0) cycle

            ! We found an atom in range of atom1 and atom2. Allocate space to
            !   store information about it and then point to it.
            if (.not.associated(currentNode2%head)) then
               allocate(currentNode2%head)
               currentNode3 => currentNode2%head
            else
               allocate(currentNode3%next)
               currentNode3 => currentNode3%next
            endif

            ! Store the valeDim indices of the current atom from the loop (m).
            currentNode3%atom3Index(1) = atom3Index(1)
            currentNode3%atom3Index(2) = atom3Index(2)

            ! Store the magnitude of the valeDim range.
            currentNode3%valeDimRangeMag = atom3Index(2) - atom3Index(1) + 1

            ! Store the minimum distance from this atom (m) to atom 1 (k) and
            !   atom 2 (l).
            currentNode3%distToAtom1 = currentMinDist(1)
            currentNode3%distToAtom2 = currentMinDist(2)

            ! Store the replicated lattice cell indices for atom3.
            currentNode3%latticeIndices(:) = latticeIndices(:,2)

            ! Store the site number for this atom.
            currentNode3%atomSiteNumber = m

            ! Compute and store the centroid for this set of three atoms.
            currentNode3%centroid(:) = (atomSites(k)%cartPos(:) + &
                  & atomSites(l)%cartPos(:) + &
                  & latticeIndices(1,1) * realVectors(:,1) + &
                  & latticeIndices(2,1) * realVectors(:,2) + &
                  & latticeIndices(3,1) * realVectors(:,3) + &
                  & atomSites(m)%cartPos(:) + &
                  & latticeIndices(1,2) * realVectors(:,1) + &
                  & latticeIndices(2,2) * realVectors(:,2) + &
                  & latticeIndices(3,2) * realVectors(:,3)) / 3.0_double

            ! Initialize the value of the bond order.
            currentNode3%bondOrder3C(:) = 0.0_double

            ! Establish that this is the current end of the list for atom3.
            currentNode3%next => null()

         enddo ! m = l, numAtomSites
      enddo ! l = k, numAtomSites
   enddo ! k = 1, numAtomSites

   ! Begin the calculation of actual bond values for the set of bonded atoms
   !   stored in the linked list data structure.
   do h = 1, spin

      ! Begin accumulating the values over each kpoint
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

         ! Begin collecting the three center bond order.
         do j = 1, numStates

            ! Skip energy states that are not occupied. We do this by first
            !   determining the array index value of the current spin-kpoint-
            !   state as defined by the tempEnergyEigenValues loop near the
            !   beginning of the population subroutine. Then we check if the
            !   population of this state is less than smallThresh. If so, then
            !   it is deemed a non-populated state and we can exit the j loop.
            orderedIndex = j+numStates*(h-1)+numStates*spin*(i-1)
            if (electronPopulation(orderedIndex) < smallThresh) exit

            ! We have found a state with at least some electron population in
            !   it. The bond order calculation will be scale according to the
            !   amount of charge.
            chargeScaleFactor = electronPopulation(orderedIndex) / &
                  & (kPointWeight(i)/real(spin,double))

            ! Loop across all atoms in the BO3C array and traverse the attached
            !   linked list data structure.
            do k = 1, numAtomSites

               ! Find an atom that has a bond to some other atom (or itself).
               if (.not.associated(BO3C(k)%head)) cycle

               ! Point at the head of the list.
               currentNode2 => BO3C(k)%head

               do while (associated(currentNode2))

                  ! Traverse the atom3 list.
                  currentNode3 => currentNode2%head

                  do while (associated(currentNode3))

                     ! Accumulate the three center bond order for the current
                     !   set of three atoms.  Also compute the two center bond
                     !   order between each pair of atoms: 1=atoms 1,2;
                     !   2=atoms 1,3; 3=atoms 2,3. This will be used to compute
                     !   the weighted centroid.
                     call compute3CBO(kPointWeight(i),real(spin,double),&
                           & chargeScaleFactor,j,BO3C(k)%valeDimRangeMag,&
                           & BO3C(k)%atom1Index(:),currentNode2%atom2Index(:),&
                           & currentNode3%bondOrder3C,1)
                     call compute3CBO(kPointWeight(i),real(spin,double),&
                           & chargeScaleFactor,j,BO3C(k)%valeDimRangeMag,&
                           & BO3C(k)%atom1Index(:),currentNode3%atom3Index(:),&
                           & currentNode3%bondOrder3C,2)
                     call compute3CBO(kPointWeight(i),real(spin,double),&
                           & chargeScaleFactor,j,currentNode2%valeDimRangeMag,&
                           & currentNode2%atom2Index(:),&
                           & currentNode3%atom3Index(:),&
                           & currentNode3%bondOrder3C,3)

                     ! Point at the next node on the list.
                     currentNode3 => currentNode3%next

                  end do ! CurrentNode3

                  ! Point at the next node on the list.
                  currentNode2 => currentNode2%next

               end do ! CurrentNode2
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

      ! Create a blank line.
      write (20,*)

      ! Now that all the bond information for this spin orientation has been
      !   computed (including all states and kpoints) we can determine the
      !   bondCentroid and print the results.

      ! Initialize a counter of the number of three center bonds.
      num3CB = 0

      ! Print information for the system as a whole.
      write (11+h,fmt="(a23)") "LATTICE_VECTORS"
      write (11+h,300) "ax ay az ",realVectors(:,1)
      write (11+h,300) "bx by bz ",realVectors(:,2)
      write (11+h,300) "cx cy cz ",realVectors(:,3)
      do i = 1, numAtomSites

         ! Find an atom that has a bond to some other atom (or itself).
         if (.not.associated(BO3C(i)%head)) cycle

         ! Point at the head of the list.
         currentNode2 => BO3C(i)%head

         do while (associated(currentNode2))

            ! Traverse the atom3 list.
            currentNode3 => currentNode2%head

            do while (associated(currentNode3))

! Establish rejection criteria.
if ((currentNode3%bondOrder3C(1) < 0.0_double) .or. &
  & (currentNode3%bondOrder3C(2) < 0.0_double) .or. &
  & (currentNode3%bondOrder3C(3) < 0.0_double)) then
currentNode3 => currentNode3%next
cycle
endif

               ! Increment the counter of the number of three center bonds
               !   that meet acceptance criteria (e.g. all three contributing
               !   two center bonds are positive).
               num3CB = num3CB + 1

               ! Compute the centroid of the bond between this set of
               !   atoms. This like the normal centroid except it is
               !   weighted according to the strength of the interaction
               !   between atoms. The deviation from the geometric
               !   centroid will speak to the degree of evenness in the
               !   three centered bond.
               call computeBondCentroid(i,currentNode2%atomSiteNumber,&
                     & currentNode3%atomSiteNumber,currentNode2%latticeIndices,&
                     & currentNode3%latticeIndices,currentNode3%bondOrder3C,&
                     & currentNode3%centroid,currentNode3%bondCentroid,&
                     & currentNode3%bondCentroidInCell)

               ! Print descriptors for this three center bond.
               write (11+h,200) "ATOMS_INVOLVED ",i,&
                     & currentNode2%atomSiteNumber,currentNode3%atomSiteNumber
               write (11+h,200) "LATTICE_INDICES_1 ",0,0,0
               write (11+h,200) "LATTICE_INDICES_2 ",currentNode2%latticeIndices
               write (11+h,200) "LATTICE_INDICES_3 ",currentNode3%latticeIndices
               write (11+h,300) "COORDS_1 ",atomSites(i)%cartPos(:)
               write (11+h,300) "COORDS_2 ",&
                     & atomSites(currentNode2%atomSiteNumber)%cartPos(:)
               write (11+h,300) "COORDS_3 ",&
                     & atomSites(currentNode3%atomSiteNumber)%cartPos(:)
               write (11+h,300) "CENTROID ",currentNode3%centroid(:)
               write (11+h,300) "BOND_CENTROID ",currentNode3%bondCentroid(:)
               write (11+h,300) "BOND_CENTROID_IN_CELL ",&
                     & currentNode3%bondCentroidInCell(:)
               write (11+h,300) "2CBO_PAIRS ",currentNode3%bondOrder3C(1:3)
               write (11+h,400) "3CBO ",currentNode3%bondOrder3C(4)

               ! Move to the next member of the node 3 list.
               currentNode3 => currentNode3%next

            end do ! List of atom3 (while loop)

            ! Move to the next member of the node 2 list.
            currentNode2 => currentNode2%next

         end do ! List of atom2 (while loop)
      enddo ! i=1, numAtomSites

      ! Record the number of three center bonds.
      write (11+h,fmt="(a17,i5)") "NUM_3C_BONDS ",num3CB

   enddo ! (h spin)

   ! Deallocate matrices that are no longer needed.
   deallocate (numOrbIndex)
   deallocate (bondIndex)
   deallocate (numAtomStates)
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


   ! Log the date and time we end.
   call timeStampEnd (26)


   ! Define all the formatted output styles.
   200 format (a23,3i5)
   300 format (a23,3e16.8)
   400 format (a23,e16.8)

end subroutine computeBond3C


subroutine getMinDist (i,j,minDistCount,latticeIndices,minDist)

   ! Import necessary modules.
   use O_Kinds
   use O_Input,       only: maxBL
   use O_AtomicSites, only: atomSites
   use O_Lattice,     only: realVectors
   use O_Constants,   only: bigThresh, smallThresh

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define passed parameters
   integer, intent (in)  :: i
   integer, intent (in)  :: j
   integer, intent (out) :: minDistCount
   integer, dimension (:), intent (out) :: latticeIndices
   real (kind=double), intent (out) :: minDist

   ! Define local variables
   integer :: m,n,o
   real (kind=double) :: currentDist
   real (kind=double), dimension (3) :: displacement

   ! Determine the smallest distance between the atoms (i) and (j).

   ! Initialize the min distance to an impossibly large number.
   minDist = bigThresh

   ! Initialize a counter for the number of atoms that are within the max bond
   !   length parameter.
   minDistCount = 0

   ! Search the original cell and the neighboring 124 cells to find the minimum
   !   distance between atom (i) and atom (j).
   do m = -2,2
   do n = -2,2
   do o = -2,2

      ! Compute the seperation vector for the atoms in this combo.
      displacement(:) = atomSites(i)%cartPos(:) - atomSites(j)%cartPos(:) - &
            & m * realVectors(:,1) - &
            & n * realVectors(:,2) - &
            & o * realVectors(:,3)

      ! Compute the distance for the atoms of this cell combo.
      currentDist = sqrt(sum(displacement(:)**2))

      ! Compare the current magnitude to the current minimum.
      if ((currentDist < minDist) .and. (abs(currentDist) > smallThresh)) then
         minDist = currentDist
         latticeIndices(1) = m
         latticeIndices(2) = n
         latticeIndices(3) = o
      endif

      ! Check if this minDist is less than the cut-off radius. If so, then
      !   increment the counter for this atom pair.  This is done in case a
      !   replicant atom has more than one configuration where it is
      !   sufficiently close to the current target atom that it regesters an
      !   overlap.  This effect is not really implemented yet since the
      !   minDist will have to track all the nearest atoms, not just the
      !   closest one.
      if ((currentDist < maxBL) .and. &
            & (abs(currentDist) > smallThresh)) then
         minDistCount = minDistCount + 1
      endif
   enddo
   enddo
   enddo

end subroutine getMinDist


subroutine compute3CBO (kPointWeight,spin,scaleFactor,j,index2Mag,index1,&
      & index2,bondOrder3C,atomPairCode)

   ! Import necessary modules.
   use O_Kinds
#ifndef GAMMA
   use O_SecularEquation, only: valeVale, valeValeOL
#else
   use O_SecularEquation, only: valeValeGamma, valeValeOLGamma
#endif

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define passed parameters.
   real (kind=double), intent (in) :: kPointWeight
   real (kind=double), intent (in) :: spin
   real (kind=double), intent (in) :: scaleFactor
   integer, intent (in) :: j ! State index number
   integer, intent (in) :: index2Mag ! Magnitude of the valeDim index2 range
   integer, dimension (:), intent (in) :: index1 ! atom1 valeDim index range
   integer, dimension (:), intent (in) :: index2 ! atom2 valeDim index range
   real (kind=double), dimension (4), intent (inout) :: bondOrder3C ! Bond order
   integer, intent (in) :: atomPairCode

   ! Define local variables.
   integer :: n
   real (kind=double) :: oneValeRealAccum

   ! Compute the three center bond order. Loop over the states of the atom at
   !   site 1 against the states of the atom at site 2 for this spin-kpoint-
   !   state (band index).
   do n = index1(1),index1(2)

#ifndef GAMMA
      ! Compute ^2 of the wave function coefficients for each element.
      waveFnSqrd(1:index2Mag) = conjg(valeVale(n,j,1,1)) * &
            & valeVale(index2(1):index2(2),j,1,1)

      ! Multiply by the overlap to get the electron number associated with
      !   these two sites in this state.  (Note that we only need to look at
      !   the real part: (real*real) + (imag*imag).)
      oneValeRealAccum = sum(&
            & real (waveFnSqrd(1:index2Mag),double) * &
            & real (valeValeOL(index2(1):index2(2),n,1,1),double) + &
            & aimag(waveFnSqrd(1:index2Mag)) * &
            & aimag(valeValeOL(index2(1):index2(2),n,1,1)))
#else
      ! Compute ^2 of the wave function for each element.
      waveFnSqrdGamma(1:index2Mag) = valeValeGamma(n,j,1) * &
            & valeValeGamma(index2(1):index2(2),j,1)

      ! Compute the effects of overlap.
      oneValeRealAccum = sum(&
            & waveFnSqrdGamma(1:index2Mag) * &
            & valeValeOLGamma(index2(1):index2(2),n,1))
#endif

      ! Store the contribution to the atom pair bond order and the atom
      !   triplet bond order.
      bondOrder3C(atomPairCode) = bondOrder3C(atomPairCode) + &
            & oneValeRealAccum * kPointWeight / spin * scaleFactor
      bondOrder3C(4) = bondOrder3C(4) + &
            & oneValeRealAccum * kPointWeight / spin * scaleFactor
   enddo

end subroutine compute3CBO

subroutine computeBondCentroid(k,l,m,latticeIndices2,latticeIndices3,&
      & currentBO3C,centroid,bondCentroid,bondCentroidInCell)

   ! Import necessary modules.
   use O_Kinds
   use O_AtomicSites, only: atomSites
   use O_Lattice, only: realVectors, invRealVectors

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define passed parameters.
   integer, intent (in) :: k ! atom site number of the first atom
   integer, intent (in) :: l ! atom site number of the second atom
   integer, intent (in) :: m ! atom site number of the third atom
   integer, dimension (3), intent (in) :: latticeIndices2
   integer, dimension (3), intent (in) :: latticeIndices3
   real (kind=double), dimension (4), intent (in) :: currentBO3C
   real (kind=double), dimension (3), intent (in) :: centroid
   real (kind=double), dimension (3), intent (inout) :: bondCentroid
   real (kind=double), dimension (3), intent (out) :: bondCentroidInCell

   ! Define local variables.
   integer :: abc
   integer :: xyz
   real (kind=double) :: magnitude12 ! Mag of displacement from mid to centroid
   real (kind=double) :: magnitude13 ! Mag of displacement from mid to centroid
   real (kind=double) :: magnitude23 ! Mag of displacement from mid to centroid
   real (kind=double) :: shiftMag12 ! Mag of the shift to be given to midpoint.
   real (kind=double) :: shiftMag13 ! Mag of the shift to be given to midpoint.
   real (kind=double) :: shiftMag23 ! Mag of the shift to be given to midpoint.
   real (kind=double), dimension (3) :: bondCentroidABC ! Fractional abc coords
   real (kind=double), dimension (3) :: atom1Pos
   real (kind=double), dimension (3) :: atom2Pos
   real (kind=double), dimension (3) :: atom3Pos
   real (kind=double), dimension (3) :: midPoint12     ! 2CBond midpoint
   real (kind=double), dimension (3) :: midPoint13     ! 2CBond midpoint
   real (kind=double), dimension (3) :: midPoint23     ! 2CBond midpoint
   real (kind=double), dimension (3) :: displacement12 ! From mid to centroid
   real (kind=double), dimension (3) :: displacement13 ! From mid to centroid
   real (kind=double), dimension (3) :: displacement23 ! From mid to centroid
   real (kind=double), dimension (3) :: unitShift12    ! Unit vector
   real (kind=double), dimension (3) :: unitShift13    ! Unit vector
   real (kind=double), dimension (3) :: unitShift23    ! Unit vector
   real (kind=double), dimension (3) :: shiftedPos1    ! x,y,z for atom 1
   real (kind=double), dimension (3) :: shiftedPos2    ! x,y,z for atom 2
   real (kind=double), dimension (3) :: shiftedPos3    ! x,y,z for atom 3

   ! Determine the positions of the atoms (including periodic lattice effects).
   atom1Pos(:) = atomSites(k)%cartPos(:) ! Always in the central cell.
   atom2Pos(:) = atomSites(l)%cartPos(:) + &
         & latticeIndices2(1) * realVectors(:,1) + &
         & latticeIndices2(2) * realVectors(:,2) + &
         & latticeIndices2(3) * realVectors(:,3)
   atom3Pos(:) = atomSites(m)%cartPos(:) + &
         & latticeIndices3(1) * realVectors(:,1) + &
         & latticeIndices3(2) * realVectors(:,2) + &
         & latticeIndices3(3) * realVectors(:,3)

   ! Compute the midpoint of each bond.
   midPoint12(:) = (atom1Pos(:) + atom2Pos(:))/2.0_double
   midPoint13(:) = (atom1Pos(:) + atom3Pos(:))/2.0_double
   midPoint23(:) = (atom2Pos(:) + atom3Pos(:))/2.0_double

   ! Move each midpoint a distance perpendicular to the bond that is equal to
   !   the inverse of the ratio of the 2-center bond order for this bond to the
   !   3-center bond order for this group of atoms.

   ! To do this we first compute the direction that needs to be moved for each
   !   midpoint. This is obtained from the displacement vector between a given
   !   midPoint and the centroid scaled by the magnitude of the displacement.
   !   (Essentially a unit vector which will be called unitShift.)
   displacement12(:) = midPoint12(:) - centroid(:)
   displacement13(:) = midPoint13(:) - centroid(:)
   displacement23(:) = midPoint23(:) - centroid(:)
   magnitude12 = sqrt(sum(displacement12(:)*displacement12(:)))
   magnitude13 = sqrt(sum(displacement13(:)*displacement13(:)))
   magnitude23 = sqrt(sum(displacement23(:)*displacement23(:)))
   unitShift12(:) = displacement12(:) / magnitude12
   unitShift13(:) = displacement13(:) / magnitude13
   unitShift23(:) = displacement23(:) / magnitude23

   ! Now we compute the magnitude of the shift to apply
   shiftMag12 = currentBO3C(4) / currentBO3C(1)
   shiftMag13 = currentBO3C(4) / currentBO3C(2)
   shiftMag23 = currentBO3C(4) / currentBO3C(3)

   ! Apply the shift in the appropriate direction.
   shiftedPos1(:) = atom1Pos(:) + shiftMag12 * unitShift12(:)
   shiftedPos2(:) = atom2Pos(:) + shiftMag13 * unitShift13(:)
   shiftedPos3(:) = atom3Pos(:) + shiftMag23 * unitShift23(:)

   ! Compute the centroid of this new triangle (where the vertices represent
   !   bond midpoints that have been shifted away from the geometric centroid
   !   by an amount proportional to the 2C bond order of the bond that the
   !   midpoint represents). Note that this position may be outside of the
   !   central cell and must therefore be shifted back inside of it before it
   !   can be used as a "atomic" position.
   bondCentroid(:) = (shiftedPos1(:) + shiftedPos2(:) + shiftedPos3(:)) / &
         & 3.0_double

   ! Compute the fractional abc coordinates of the bondCentroid.
   do abc = 1,3 ! a,b,c axes
      bondCentroidABC(abc) = 0.0_double
      do xyz = 1,3 ! x,y,z axes
         bondCentroidABC(abc) = bondCentroidABC(abc) + bondCentroid(xyz) * &
               invRealVectors(xyz,abc)
      enddo
   enddo

   ! If any value of the bondCentroid vector is > 1 then we reduce it by 1.
   !   Similarly, if any value is < 0 then we increase it by 1. Some special
   !   care needs to be taken for the case that the bondCentroid is on the edge
   !   of the unit cell.
   do abc = 1,3
      if (bondCentroidABC(abc) > 1.0_double) then
         bondCentroidABC(abc) = bondCentroidABC(abc) - 1
      elseif (bondCentroidABC(abc) < 0.0_double) then
         bondCentroidABC(abc) = bondCentroidABC(abc) + 1
      endif
   enddo

   ! Now we convert back to direct space x,y,z coordinates.
   do xyz = 1,3
      bondCentroidInCell(xyz) = 0.0_double
      do abc = 1,3
         bondCentroidInCell(xyz) = bondCentroidInCell(xyz) + &
               & bondCentroidABC(abc) * realVectors(xyz,abc)
      enddo
   enddo

end subroutine computeBondCentroid


end module O_Bond3C
