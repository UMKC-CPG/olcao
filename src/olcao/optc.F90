module O_OptcTransitions

   ! Import necessary modules
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define module data.
   real (kind=double) :: orbitalDiff ! Energy difference between highest
         !   occupied state and lowest unoccupied state. (Normal energy onset.)
   real (kind=double) :: energyCutoff ! The energy cut off value determined
         !   from the input data and the type of calculation.  This represents
         !   the highest energy level above 0 included in the calc.
   real (kind=double) :: energyMin ! The minimum energy in the energy window
         !   to be computed for.
   real (kind=double) :: maxTransEnergy ! The maximum transition energy
         !   determined from the input data and the type of calculation.  This
         !   represents the largest transition energy included in the calc.
   real (kind=double), allocatable, dimension (:) :: energyScale ! The fine
         !   scale used for the energy axis.

   integer :: maxPairs ! The largest value in the transCounter array defined
         !   below.
   integer, allocatable, dimension (:,:) :: transCounter ! A count of the
         !   number of transitions present for each kpoint of both spin dirs.

   integer, allocatable, dimension (:,:) :: firstOccupiedState ! Array with the
         !   index number of the first occupied state for each kpoint and spin.
   integer, allocatable, dimension (:,:) :: lastOccupiedState ! Array with the
         !   index number of the last occupied state for each kpoint and spin.
   integer, allocatable, dimension (:,:) :: firstUnoccupiedState ! Array with
         !   the index number of the first unoccupied state for each
         !   kpoint and spin.
   integer, allocatable, dimension (:,:) :: lastUnoccupiedState ! Array with the
         !   index number of the last unoccupied state for each kpoint and spin.

   real (kind=double), allocatable, dimension (:) :: indirectGap ! This is the
         !   smallest energy difference between the highest occupied state and
         !   the lowest unoccupied state for all kpoints.
   real (kind=double), allocatable, dimension (:,:) :: directGap ! This is the
         !   smallest energy difference between the highest occupied state and
         !   the lowest unoccupied state for any single kpoint.
   real (kind=double), allocatable, dimension (:,:) :: indirectGapEnergies !
         !   This is the list of upper and lower energy values for the
         !   directGap of each kpoint and spin.

   real (kind=double), allocatable, dimension (:,:,:,:) :: transitionProb
   real (kind=double), allocatable, dimension (:,:,:)   :: energyDiff

#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:,:) :: valeValeMom
   complex (kind=double), allocatable, dimension (:,:,:)   :: coreValeOL
#else
   real (kind=double), allocatable, dimension (:,:,:) :: valeValeMomGamma
   real (kind=double), allocatable, dimension (:,:)   :: coreValeOLGamma
#endif

   real (kind=double), allocatable, dimension (:,:) :: sigmaEAccumulator

contains

subroutine getEnergyStatistics

   ! Import necessary modules.
   use O_Kinds
   use O_Potential,       only: spin
   use O_CommandLine,     only: stateSet
   use O_SecularEquation, only: energyEigenValues
   use O_Populate,        only: electronPopulation
   use O_KPoints,         only: numKPoints, kPointWeight
   use O_Constants,       only: dim3, smallThresh, bigThresh, hartree
   use O_Input, only: numStates, cutoffEnOPTC, maxTransEnOPTC, &
         & totalEnergyDiffPACS, energyWindowPACS, firstInitStatePACS, &
         & lastInitStatePACS, onsetEnergySlackPACS, cutoffEnSIGE, &
         & maxTransEnSIGE

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables.
   integer :: h,i,j,k ! Loop index variables.
   real (kind=double) :: maxOccupiedEnergy
   real (kind=double) :: minUnoccupiedEnergy
   real (kind=double) :: currentEnergyDiff
   real (kind=double) :: currentGap
   integer :: firstInit
   integer :: lastInit
   integer :: firstFin
   integer :: lastFin
   integer :: orderedIndex

   ! Pull variables out of imported modules.
   if (stateSet == 0) then      ! Standard optical properties calculation.
      ! The energy onset for standard optical properties calculations is the
      !   band gap width in eV.  The energy scale for evaluating the
      !   accumulated (and broadened) transitions will begin as close to 0 eV
      !   as conveniently possible.
      energyCutoff   = cutoffEnOPTC
      maxTransEnergy = maxTransEnOPTC
   elseif (stateSet == 1) then  ! PACS type calculation
      ! PACS calculations have the interesting feature that the energy onset is
      !   at some very high energy that is dependent on the particular target
      !   atom being excited (elemental dependency).  This is because the
      !   transitions we are considering are from a deep core state to the
      !   conduction band.  We want that energy onset to be defined by the
      !   difference in total energy between the ground state and the excited
      !   state instead of the orbital energy difference.  (WHY?)

      !   To make that happen we must first subtract away the original energy
      !   onset which is defined by the orbital energy difference between the
      !   initial state and the lowest energy state that the electron is
      !   excited into (bottom of the conduction band).
      energyCutoff   = bigThresh

      ! Compute the max transition energy as the total energy difference
      !   between ground and excited states.  Since we want the first value to
      !   be a nice round number (multiple of 5) we subtract out so-called
      !   "slack" which is the remainder of an integer division by 5.  This is
      !   the lowest energy in the output data set (no transitions at this
      !   energy though).  Then, we add the energy window we want to compute
      !   for to get the maximum transition energy.
      energyMin = totalEnergyDiffPACS - &
            & mod(totalEnergyDiffPACS,onsetEnergySlackPACS)
      maxTransEnergy = energyMin + energyWindowPACS
!write (20,*) maxTransEnergy,energyMin,deltaPACS
!call flush (20)
   else                         ! Sigma(E) type calculation
      energyCutoff   = cutoffEnSIGE
      maxTransEnergy = maxTransEnSIGE
   endif


   ! Allocate arrays
   allocate (transCounter         (numKPoints,spin))
   allocate (firstOccupiedState   (numKPoints,spin))
   allocate (lastOccupiedState    (numKPoints,spin))
   allocate (firstUnoccupiedState (numKPoints,spin))
   allocate (lastUnoccupiedState  (numKPoints,spin))
   allocate (directGap            (numKPoints,spin))
   allocate (indirectGapEnergies  (spin,2)) ! Last index: Upper=1, Lower=2
   allocate (indirectGap          (spin))

   ! Initialize the arrays.
   transCounter(:,:)         = 0
   firstOccupiedState(:,:)   = 0
   lastOccupiedState(:,:)    = 0
   firstUnoccupiedState(:,:) = 0
   lastUnoccupiedState(:,:)  = 0
   directGap(:,:)            = bigThresh
   indirectGapEnergies(:,1)  = bigThresh
   indirectGapEnergies(:,2)  = -bigThresh
   inDirectGap(:)            = 0.0_double

   ! The purpose of this subroutine is to gather important statistics and
   !   indices for use later on.  The important values that will be determined
   !   are:  1) The first and last occupied state and the first and last 
   !   unoccupied state for each kpoint and each spin that each different type
   !   of calculation cares about (e.g. PACS transitions are from one core
   !   state to numerous CB states; OPTC transitions are from many VB states to
   !   numerous CB states; and SIGE considers a finite range near the fermi
   !   level); 2) The direct band gap for each spin; 3) The minimum indirect
   !   band gap for each spin; 4) The number of transitions for each kpoint and
   !   spin so that we can allocate memory easily; 5) whatever ...
   ! We must take some extra consideration of the situation in which thermal
   !   smearing is present. It should be understood that in this case the
   !   population statistics are a bit different than the 0K case where the
   !   Fermi edge is a flat step function. For thermal smearing we have the
   !   condition that states near the Fermi level that are fully occupied at
   !   0K will be partially occupied at a finite temperature. Further, states
   !   that are totally unoccupied at 0K will also be partially occupied at a
   !   finite temperature. Thus we must consider the case of transitions
   !   between occupied and unoccupied states a bit carefully.
   ! Consider first the PACS case. Here we have the initial state(s) being core
   !   states that will always be fully occupied in the ground state regardless
   !   of temperature (becaues they are so deep). The electron from the core
   !   state may transition into any state with partial or zero occupation.
   !   Thus, when making a list of the states to transition into, we must find
   !   the first one that has a non-negligable portion that is unoccupied (i.e.
   !   a mostly occupied VB state that is a bit far-ish from the traditional
   !   top of the VB.) The core electron may transition into this state, but
   !   the probability of doing so must be affected by (of course) the
   !   momentum matrix element that represents conservation of angular momentum
   !   selection rules (i.e. s->p and p->s,d etc) *and* the fact that the
   !   state at that energy is mostly occupied already (e.g. 99%). Therefore,
   !   the transition calculation proceeds as normal except that the intensity
   !   must be scaled by a factor of 0.01. Without thermal smearing, the
   !   probability of transitioning into this state would be 0%. As the higher
   !   energy states are considered, the so-called occupation scaling factor
   !   will increase to a maximum of 1.0 in the case that the core electron is
   !   transitioning into a state that has zero initial occupation.
   ! Consider second the case of a traditional VB optical properties
   !   calculation. Now the situation is even more complicated because the
   !   initial state may be partially occupied as well as the final state.
   !   However, similar physical principles will apply and the probability of
   !   making a transition will be scaled by two occupation scaling factors,
   !   one for the originating state and one for the final state.

   do h = 1, spin

      ! Initialize the values used to find the indirect gap and the direct gap.
      maxOccupiedEnergy = -bigThresh
      minUnoccupiedEnergy = bigThresh

      ! Pacs calculations need to have the resultant spectra shifted according
      !   to the difference in orbital energies from the ground and excited
      !   states.  To find the amount of the shift we initialize the search
      !   number.
      if (stateSet == 1) then  ! Doing PACS calculation.
         orbitalDiff = bigThresh
      endif

      do i = 1, numKPoints
         do j = 2, numStates

            ! Find the last occupied state for this KPoint and spin. Also note
            !   that this will deal with degenerate highest occupied states.
            ! Note: In the case that thermal smearing is turned on, then the
            !   zero of energy is the Fermi level (which may be in the gap for
            !   insulators). There may be paritially occupied orbitals above
            !   the Fermi level and partially occupid orbitals below the Fermi
            !   level because of the smearing.

            ! Transitions are from the core states into any of the partially
            !   occupied states scaled by the degree of population for PACS
            !   calculations.
            ! Transitions are from the valence states (including any partially
            !   occupied states above or below the fermi level) into any of the
            !   partially occupied states above or below the fermi level (that
            !   are also *above* the initial state) for traditional optical
            !   properties calculations.
            ! Thermal excitations in conducting materials are taken to occur
            !   between any occupied (or partially occupied) state and any
            !   unoccupied (or partially occupied state of higher energy) that
            !   is within a designated range (usually on the order of a few eV
            !   or tenths of an eV). This is for sigma(E) calculations.

            ! Determine the array index value of the current spin-kpoint-state
            !   as defined by the tempEnergyEigenValues loop near the beginning
            !   of the population subroutine.
            orderedIndex = j+numStates*(h-1)+numStates*spin*(i-1)

            ! Note that we should anticipate that the following "if" statement
            !   will never be true when j=1. This means that the electron
            !   population of the first (lowest) state will always be the
            !   maximum possible (the kPointWeight). This condition might be
            !   broken if the thermal smearing was set so insanely high that
            !   it smeared all the way to the lowest state (some -20 eV). This
            !   would correspond to rediculous temperatures. Just to be sure
            !   though, the j-loop starts at 2. (The first state will not be
            !   the first (partially) unoccupied state.) This is important
            !   because some j-1 and j-k with k/=0 type calculations are done
            !   inside this if-statement.
            ! This "if" statement will be true any time that we find a state
            !   with less than full occupation.
            if (abs(electronPopulation(orderedIndex)-kPointWeight(i) / &
                  & real(spin,double)) > smallThresh) then
!            if (energyEigenValues(j,i,h) > 0.0_double) then
!write (20,*) "ePop,ordIndex,kPWeight,j,i,h",electronPopulation(orderedIndex),&
!   & orderedIndex,kPointWeight(i),j,i,h
!call flush (20)

               ! In the condition that we are doing a PACS calculation then
               !   the initial state core state(s) *will* be partially
               !   occupied because we have pulled an electron out of them.
               !   These states will never be final states and so we have to
               !   skip them.
               if (stateSet == 1) then
                  if (j <= lastInitStatePACS) cycle
               endif

               ! We have found a state that is either partially or fully
               !   unoccupied. Therefore, this state may be the lowest energy
               !   state that an electron *could* transfer into (indirect gap)
               !   considering that each spin and kpoint has its own set of
               !   states. If this state is lower in energy than any other
               !   previously found state with non-occupation then we record it
               !   as the final state for the indirect gap for this spin. (Note
               !   that for metals with thermal smearing turned on, this state
               !   will be below the fermi level while there will be partially
               !   occupied states above the fermi level. Thus, there will be
               !   no gap. Metals should have no gap.)
               if (indirectGapEnergies(h,1) > energyEigenValues(j,i,h)) then
                  indirectGapEnergies(h,1) = energyEigenValues(j,i,h)
               endif

               ! Store the index values for the first unoccupied state for this
               !   kpoint and spin.
               if (firstUnoccupiedState(i,h) == 0) then
                  firstUnoccupiedState(i,h) = j
               endif

               ! In the event that the previous state is partially or fully
               !   occupied then we want to check the energy difference
               !   between it and the current state (which is partially
               !   occupied or fully unoccupied) to see if it is a "gap". For
               !   the 0K case with insulators, this is an obvious calculation,
               !   but for the finite temperature case with thermal smearing
               !   the situation is more complex.  Essentailly, for the current
               !   kpoint, we are looking for the largest separation between a
               !   state with some electrons in it and a state that isn't fully
               !   occupied. At the moment we are at a state that has less than
               !   full occupation and we are checking the previous state to
               !   see if it has any electrons.
               if (electronPopulation(orderedIndex-1) > 0.0_double) then

                  ! Obviously, this is also an occupied state and the last time
                  !   that we get inside this "if" statement (inside the upper
                  !   one too) we will have found the last occupied state.
                  !   Thus, we assume that every time is the last time and
                  !   eventually it will be correct.  (This is only useful for
                  !   traditional VB optical properties calculations and
                  !   sigma(E) calculations. It is not useful for PACS
                  !   calculatoins because they will use the user specified
                  !   core states by overriding this determination later on.)
                  lastOccupiedState(i,h) = j-1

                  ! Similarly, we will find the highest occupied state energy.
                  if (indirectGapEnergies(h,2)<energyEigenValues(j-1,i,h)) then
                     indirectGapEnergies(h,2) = energyEigenValues(j-1,i,h)
                  endif

                  ! Compute the size of the energy difference between the
                  !   current partially occupied or fully unoccupied state and
                  !   the previous (j-1) state which is either partially
                  !   occupied or fully occupied. Note that we had to check
                  !   that the previous state was at least partially occupied.
                  currentGap = abs(energyEigenValues(j,i,h) - &
                        & energyEigenValues(j-1,i,h))

                  ! This so-called currentGap will become the direct gap *for
                  !   this kpoint and spin* if it is minimal compared to other
                  !   previous currentGap values for this kpoint (and spin).
                  !   Note that the comparison between kpoints will be done
                  !   later. Also, just remember that we are always only
                  !   comparing the energy difference between adjacent states
                  !   for the current kpoint.
                  if (currentGap < directGap(i,h)) then
                     directGap(i,h) = currentGap
!write (20,*) "currGap,dirGap,i,h",currentGap,directGap(i,h),i,h
!call flush(20)
                     ! The largest currentGap for this kpoint might be a
                     !   directGap when compared with other directGap values
                     !   from other kpoints (to be determined later).
                  endif
               else ! The previous state has no electrons in it.
                  ! No need to search through any higher states in the outer j
                  !   loop.  All relevant information has been obtained.
                  exit
               endif ! The previous state has at least some electrons in it.

            endif ! (Found a partially occupied or fully unoccupied state)
         enddo ! (j = numStates)

         ! If we are doing a sige calculation then we are only concerned
         !   with states that are a few eV from the Fermi energy.  This
         !   means that in addition to having to seek out the first
         !   unoccupied state and the last occupied state we also need to
         !   identify the last unoccupied state (which will be *just*
         !   above the last occupied state according to the
         !   maxTransEnergy given in the input file) and the first
         !   occupied state (which will be *just* below the first
         !   unoccupied state according to the maxTransEnergy given in
         !   the input file (olcao.dat)).
         ! Note that this search only needs to be done once per kpoint so this
         !   is why we are doing it after the numStates loop.
         if (stateSet == 2) then

            ! Loop higher than the last occupied state to find the
            !   lowest unoccupied state that is *greater* than the maximum
            !   transition energy from the Fermi level.
            do k = lastOccupiedState(i,h), numStates
               if (energyEigenValues(k,i,h) > maxTransEnergy) then
                   lastUnoccupiedState(i,h) = k
                   exit
               endif
               if (k == numStates) then
                  lastUnoccupiedState(i,h) = numStates
               endif
            enddo

            ! Loop lower than the first unoccupied state to find the
            !   highest occupied state that is *less* than the maximum
            !   transition energy from the Fermi level.
            do k = 1, firstUnoccupiedState(i,h)-1
               if (abs(energyEigenValues(firstUnoccupiedState(i,h)-k,i,h)) > &
                        & maxTransEnergy) then
                  firstOccupiedState(i,h) = firstUnoccupiedState(i,h)-k
                  exit
               endif
               if (k == firstUnoccupiedState(i,h)-1) then
                  firstOccupiedState(i,h) = 1
               endif
            enddo
         endif ! stateSet 2

         ! Initialize the counter for the number of transitions for this kpoint.
         transCounter(i,h) = 0

         ! For normal optical properties calculations the last unoccupied state
         !   we care about is the last (highest) state in the calculation, and
         !   the first occupied state is always the first (lowest) state in the
         !   calculation.  For PACS, these values depend on which core state
         !   has been excited.  For Sigma(E) these values depend on the range
         !   around the Fermi level to consider (defined by the user) and the
         !   determination above.

         if (stateSet == 0) then ! Normal optical properties calculation.
            firstOccupiedState(i,h) = 1
            ! lastOccupiedState determined above.
            ! firstUnoccupiedState determined above.
            lastUnoccupiedState(i,h) = numStates
         elseif (stateSet == 1) then ! Doing a PACS calculation.
            firstOccupiedState(i,h) = firstInitStatePACS ! From O_Input
            lastOccupiedState(i,h)  = lastInitStatePACS  ! From O_Input
            ! firstUnoccupiedState determined above.
            lastUnoccupiedState(i,h) = numStates
         elseif (stateSet == 2) then ! Doing a Sigma(E) calculation.
            ! firstOccupiedState determined above.
            ! lastOccupiedState determined above.
            ! firstUnoccupiedState determined above.
            ! lastUnoccupiedState determined above.
         else
            ! Error, no other options.
            stop "Check optical properties command line parameter: stateSet"
         endif

         ! Store the state variables for temporary use as loop indices.
         firstInit = firstOccupiedState(i,h)
         lastInit  = lastOccupiedState(i,h)
         firstFin  = firstUnoccupiedState(i,h)
         lastFin   = lastUnoccupiedState(i,h)

         if (stateSet == 1) then ! Doing PACS calculation.
            ! Determine the orbital energy difference for this kpoint and
            !   compare it to the smallest difference yet obtained.
            orbitalDiff = min(orbitalDiff,&
                  & energyEigenValues(firstFin,i,h) - &
                  & energyEigenValues(lastInit,i,h))
         endif

         ! Loop over all the possible transitions for this kpoint to determine
         !   the number of accepted transitions for this kpoint.  We will also
         !   refine the value for the lastUnoccupied state to consider.
         do j = firstInit, lastInit
            do k = firstFin, lastFin

               ! An important note is that with thermal smearing turned on,
               !   there is the potential for a state to be *both* an initial
               !   and a finel state because it is partially occupied. If we
               !   encounter any such states we will not count the case where
               !   the final state is lower in energy that the initial state.
               if (j >= k) cycle

               ! If the energy of the final state is higher than the requested
               !   cut-off then we adjust the record for the last unoccupied
               !   state and go to the next initial state because all the
               !   remaining final states will be greater.
               if (energyEigenValues(k,i,h) > energyCutoff) then
                  lastUnoccupiedState(i,h) = k-1
                  exit
               endif

               ! Compute the energy of transition between the current states.
               currentEnergyDiff = energyEigenValues(k,i,h) - &
                     & energyEigenValues(j,i,h)

               ! Check if the energy difference is less than the maximum
               !   transition energy that the input file requested computation
               !   for.  If it fails, then we adjust the record for the last
               !   unoccupied state and go to the next initial state because
               !   all the remaining final states for this energy will be
               !   greater.
               if (currentEnergyDiff > maxTransEnergy) then
                  lastUnoccupiedState(i,h) = k-1
                  exit
               endif

               ! At this point the transition is valid and one we would want to
               !   compute.  Unfortunately there isn't a good way to save the
               !   energyDiff computation that we did above for later use, and
               !   it will have to be done again.  (BOO)  (Unless we do some
               !   static memory allocation. (BOO))
               transCounter(i,h) = transCounter(i,h) + 1
            enddo ! (k First to Last Fin)
         enddo ! (j First to Last Init)
      enddo ! (i kpoints)

      ! Determine the indirect band gap.
      indirectGap(h) = indirectGapEnergies(h,1) - indirectGapEnergies(h,2)
      if (indirectGap(h) < 0.0_double) then
         indirectGap(h) = 0.0_double
      endif

      ! Write the indirect band gap and determine+write the direct band gap.
      if (spin == 1) then
         write (20,*) "Indirect Band Gap(eV) = ",indirectGap(h)*hartree
         write (20,*) "Direct Band Gap(eV)   = ",minval(directGap(:,h))*hartree
         call flush (20)
      elseif (h == 1) then
         write (20,*) "(Up) Indirect Band Gap(eV) = ",indirectGap(h)*hartree
         write (20,*) "(Up) Direct Band Gap(eV)   = ",minval(directGap(:,h))* &
               & hartree
         call flush (20)
      else ! spin == 2 and h == 2
         write (20,*) "(Dn) Indirect Band Gap(eV) = ",indirectGap(h)*hartree
         write (20,*) "(Dn) Direct Band Gap(eV)   = ",minval(directGap(:,h))* &
               & hartree
         call flush (20)
      endif

      ! Determine the maximum number of transition pairs of all the kpoints and
      !   for both spins.
      maxPairs = max(maxPairs,maxval(transCounter(:,h)))
!write (20,*) "maxPairs=",maxPairs
!call flush (20)
   enddo ! (h spin)

   ! Now that the number of transitions for each kpoint are known we can 
   !   allocate space to hold information based on the number transitions.
   if (stateSet /= 2) then  ! Not doing a Sigma(E) calculation.
      allocate (energyDiff (maxPairs,numKPoints,spin))
      allocate (transitionProb  (dim3,maxPairs,numKPoints,spin))

      ! Initialize these arrays.
      energyDiff(:maxPairs,:numKPoints,:) = 0.0_double
      transitionProb(:dim3,:maxPairs,:numKPoints,:) = 0.0_double
   endif

end subroutine getEnergyStatistics



subroutine computeTransitions


   ! Import the necessary modules.
   use HDF5
   use O_Kinds
   use O_TimeStamps
   use O_Potential,     only: spin
   use O_KPoints,       only: numKPoints
   use O_IntegralsPSCF, only: getIntgResults
   use O_AtomicSites,   only: coreDim, valeDim
   use O_CommandLine,   only: stateSet, serialXYZ
   use O_Input,         only: numStates, totalEnergyDiffPACS
   use O_PSCFBandHDF5,  only: valeValeBand, coreValeBand, valeStatesBand, &
         & eigenVectorsBand_did, eigenVectorsBand2_did, coreValeBand_did, &
         & valeValeBand_did
#ifndef GAMMA
   use O_SecularEquation, only: valeVale
   use O_MatrixSubs,      only: readMatrix, readPartialWaveFns
#else
   use O_SecularEquation, only: valeValeGamma
   use O_MatrixSubs,      only: readMatrixGamma, readPartialWaveFnsGamma
#endif

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define local variables.
   integer :: h,i,j ! Loop index variables
   real (kind=double) :: energyShift
   real (kind=double), allocatable, dimension (:,:) :: tempRealValeVale
#ifndef GAMMA
   real (kind=double), allocatable, dimension (:,:) :: tempImagValeVale
#endif


   ! Log the date and time we start.
   call timeStampStart (23)

   ! Allocate the matrix to hold the wave functions and initialize it.
#ifndef GAMMA
   allocate (valeVale (valeDim,numStates,1,1))
   valeVale(:,:,1,1) = cmplx(0.0_double,0.0_double,double)
#else
   allocate (valeValeGamma(valeDim,numStates,1))
   valeValeGamma(:,:,1) = 0.0_double
#endif


   do h = 1, spin

      ! Record the fact that we are starting the k-point loop so that when
      !   someone looks at the output file as the job is running they will
      !   know that all the little dots represent a count of the number of
      !   k-points. Then they can figure out the progress and progress rate.
      write (20,*) "Beginning k-point loop."
      if (numKPoints > 1) write (20,*) "Expecting ",numKPoints," iterations."
      call flush (20)

      ! Begin a loop over the number of kpoints
      do i = 1, numKPoints


         ! Allocate temporary reading matrices.
#ifndef GAMMA
         allocate (tempRealValeVale(valeDim,numStates))
         allocate (tempImagValeVale(valeDim,numStates))
#else
         allocate (tempRealValeVale(valeDim,numStates))
#endif

         if (stateSet /= 1) then  ! Not doing a PACS calculation.

            ! Read the datasets for this kpoint.
#ifndef GAMMA
            call readMatrix(eigenVectorsBand_did(:,i,h),valeVale(:,:,1,1),&
                  & tempRealValeVale(:,:),tempImagValeVale(:,:),&
                  & valeStatesBand,valeDim,numStates)
#else
            call readMatrixGamma(eigenVectorsBand_did(1,i,h),&
                  & valeValeGamma(:,:,1),valeStatesBand,valeDim,numStates)
#endif
         else

#ifndef GAMMA
            ! Read the data for the ground state for this kpoint.
            call readPartialWaveFns(eigenVectorsBand_did(:,i,h),&
                  & valeVale(:,:,1,1),tempRealValeVale(:,:),&
                  & tempImagValeVale(:,:),valeStatesBand,&
                  & firstOccupiedState(i,h),lastOccupiedState(i,h),&
                  & valeDim,numStates)

            ! Read the data for the excited state for this kpoint.
            call readPartialWaveFns(eigenVectorsBand2_did(:,i,h),&
                  & valeVale(:,:,1,1),tempRealValeVale(:,:),&
                  & tempImagValeVale(:,:),valeStatesBand,&
                  & lastOccupiedState(i,h)+1,lastUnoccupiedState(i,h),&
                  & valeDim,numStates)
#else
            ! Read the data for the ground state for this kpoint.
            call readPartialWaveFnsGamma(eigenVectorsBand_did(1,i,h),&
                  & valeValeGamma(:,:,1),tempRealValeVale(:,:),valeStatesBand,&
                  & firstOccupiedState(i,h),lastOccupiedState(i,h),&
                  & valeDim,numStates)

            ! Read the data for the excited state for this kpoint.
            call readPartialWaveFnsGamma(eigenVectorsBand2_did(1,i,h),&
                  & valeValeGamma(:,:,1),tempRealValeVale(:,:),valeStatesBand,&
                  & lastOccupiedState(i,h)+1,lastUnoccupiedState(i,h),&
                  & valeDim,numStates)
#endif
         endif


         ! Read the orthogonalization coefficients after allocating space to
         !   hold them (and the temp reading matrix).
#ifndef GAMMA
         deallocate (tempRealValeVale)
         deallocate (tempImagValeVale)
         allocate   (tempRealValeVale (coreDim,valeDim))
         allocate   (tempImagValeVale (coreDim,valeDim))
         allocate   (coreValeOL (coreDim,valeDim,1))
         if (coreDim /= 0) then
            call readMatrix(coreValeBand_did(:,i),coreValeOL(:,:,1),&
                  & tempRealValeVale(:,:),tempImagValeVale(:,:),&
                  & coreValeBand,coreDim,valeDim)
         endif
         deallocate (tempRealValeVale)
         deallocate (tempImagValeVale)
#else
         deallocate (tempRealValeVale)
         allocate   (coreValeOLGamma  (coreDim,valeDim))
         if (coreDim /= 0) then
            call readMatrixGamma(coreValeBand_did(1,i),coreValeOLGamma(:,:),&
                  & coreValeBand,coreDim,valeDim)
         endif
#endif


         ! Perform the computations in serial or all together.
         if (serialXYZ == 0) then
#ifndef GAMMA

            allocate   (valeValeMom (valeDim,valeDim,1,3))
            ! Get the integral results for the x, y, z momentum matrices.
            ! Runcode:  3 = XMom, 4 = YMom, 5 = ZMom; (j+2)
            do j = 1, 3
               call getIntgResults (valeValeMom(:,:,:,j),coreValeOL,&
                     & i,j+2,valeValeBand_did(i),valeValeBand,0,1)
            enddo
            ! Deallocate matrices that are no longer needed in this iteration to
            !   make space for those that are needed.
            deallocate (coreValeOL)
#else
            allocate   (valeValeMomGamma (valeDim,valeDim,3))
            ! Get the integral results for the x, y, z momentum matrices.
            ! Runcode:  3 = XMom, 4 = YMom, 5 = ZMom; (j+2)
            do j = 1, 3
               call getIntgResults (valeValeMomGamma(:,:,j),coreValeOLGamma,&
                     & j+2,valeValeBand_did(i),valeValeBand,0,1)
            enddo
            ! Deallocate matrices that are no longer needed in this iteration to
            !   make space for those that are needed.
            deallocate (coreValeOLGamma)
#endif

            if (stateSet /= 2) then  ! Not doing a Sigma(E) calculation.
               call computePairs (i,0,h)
            else
               call computeSigmaE (i,0,h)
            endif

         else
            do j = 1, 3
#ifndef GAMMA
               allocate   (valeValeMom (valeDim,valeDim,1,1))
               ! Runcode:  3 = XMom, 4 = YMom, 5 = ZMom; (j+2)
               call getIntgResults (valeValeMom(:,:,:,1),coreValeOL,&
                     & i,j+2,valeValeBand_did(i),valeValeBand,0,1)
#else
               allocate   (valeValeMomGamma (valeDim,valeDim,1))
               ! Runcode:  3 = XMom, 4 = YMom, 5 = ZMom; (j+2)
               call getIntgResults (valeValeMomGamma(:,:,1),coreValeOLGamma,&
                     & j+2,valeValeBand_did(i),valeValeBand,0,1)
#endif

               if (stateSet /= 2) then  ! Not doing a Sigma(E) calc.
                  call computePairs (i,j,h)
               else
                  call computeSigmaE (i,j,h)
               endif
            enddo

#ifndef GAMMA
            deallocate (coreValeOL)
#else
            deallocate (coreValeOLGamma)
#endif
         endif


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

      enddo

      ! Add a final return if one was not already made.
      if (mod(numKPoints,10) /= 0) then
         write (20,*)
      endif


      ! In the case of XANES calculations we need to shift the calculated
      !   energy values.  The amount of the shift is equal to the difference
      !   between (the total energy difference between the ground and excited
      !   states) and (the calculated LUMO,core difference).
      if (stateSet == 1) then  ! Doing PACS calculation.
         energyShift = totalEnergyDiffPACS - orbitalDiff
         energyDiff(:,:,h) = energyDiff(:,:,h) + energyShift
      endif
   enddo ! (h spin)

   ! Deallocate unnecessary matrices and arrays
#ifndef GAMMA
   deallocate (valeVale)
#else
   deallocate (valeValeGamma)
#endif
   deallocate (firstOccupiedState)
   deallocate (lastOccupiedState)
   deallocate (firstUnoccupiedState)
   deallocate (lastUnoccupiedState)
   deallocate (directGap)
   deallocate (indirectGapEnergies)
   deallocate (indirectGap)

   ! Log the date and time we end.
   call timeStampEnd (23)

end subroutine computeTransitions



subroutine computePairs (currentKPoint,xyzComponents,spinDirection)

   ! Import the necessary modules.
   use O_Kinds
   use O_Constants,   only: dim3
   use O_Potential,   only: spin
   use O_AtomicSites, only: valeDim
   use O_CommandLine, only: stateSet
   use O_SortSubs,    only: mergeSort
   use O_Input,       only: numStates
   use O_KPoints,     only: kPointWeight
   use O_Populate,    only: electronPopulation
#ifndef GAMMA
   use O_SecularEquation, only: energyEigenValues, valeVale
#else
   use O_SecularEquation, only: energyEigenValues, valeValeGamma
#endif

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer :: currentKPoint
   integer :: xyzComponents ! 0=all, 1=x, 2=y, 3=z
   integer :: spinDirection

   ! Define local variables.
   integer :: i,j,k ! Loop index variables
   integer :: initComponent
   integer :: finComponent
   integer :: transPairCount
   integer :: firstInit
   integer :: lastInit
   integer :: firstFin
   integer :: lastFin
   integer :: finalStateIndex
   integer :: orderedIndex
   integer, allocatable, dimension (:) :: sortOrder
   integer, allocatable, dimension (:) :: segmentBorders
   real    (kind=double) :: initStateFactor
   real    (kind=double) :: finStatefactor
   real    (kind=double) :: currentEnergyDiff
   real    (kind=double), allocatable, dimension (:)     :: energyDiffTemp
   real    (kind=double), allocatable, dimension (:,:)   :: transitionProbTemp
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: conjWaveMomSum
   complex (kind=double),              dimension (dim3)  :: valeValeXMom
#else
   real    (kind=double), allocatable, dimension (:,:,:) :: conjWaveMomSumGamma
   real    (kind=double),              dimension (dim3)  :: valeValeXMomGamma
#endif

   ! Initialize a counter for the current number of transition pairs
   transPairCount = 0

   ! Make shorthand for the state indices.
   firstInit = firstOccupiedState(currentKPoint,spinDirection)
   lastInit  = lastOccupiedState(currentKPoint,spinDirection)
   firstFin  = firstUnoccupiedState(currentKPoint,spinDirection)
   lastFin   = lastUnoccupiedState(currentKPoint,spinDirection)

   ! Determine the range of components (xyz) that should be considered.
   if (xyzComponents == 0) then
      initComponent = 1
      finComponent = 3
   else
      initComponent = 1
      finComponent = 1
   endif


#ifndef GAMMA

   ! Allocate space to hold the sum(conjg(valeVale(:,j)) * valeVale_Mom(:,k,1))
   !   for each of the possible final states.  This is done since the values
   !   are independent of the initial states.  The 1 for the valeVale is
   !   for the 1 kpoint.  The finComponent is 3 for all three at once, and
   !   1 for when X, Y, Z are done separately.
   allocate (conjWaveMomSum (valeDim,lastFin-firstFin+1,finComponent))

   ! Compute the sum over the final states.
   do i = initComponent, finComponent
      finalStateIndex = 0
      do j = firstFin, lastFin
         ! Define the final index for conjWaveMomSum
         finalStateIndex = finalStateIndex + 1
         do k = 1, valeDim
            conjWaveMomSum(k,finalStateIndex,i) = &
                  & sum(conjg(valeVale(:,j,1,1)) * valeValeMom(:,k,1,i))
         enddo
      enddo
   enddo

   ! Deallocate the momentum matrix elements.
   deallocate (valeValeMom)

#else
   ! Documentation similar to the above non-gamma case.
   allocate (conjWaveMomSumGamma (valeDim,lastFin-firstFin+1,finComponent))


   ! Compute the sum over the final states.
   do i = initComponent, finComponent

      ! Make the upper triangle correct for Hermiticity.  Recall that for
      !   the Gamma K Point all the matrices are real (except the momentum
      !   matrix which was multiplied by a -i and is hence imaginary).
      !   Since it must be Hermitian we need to apply that now.
      do j = 1, valeDim
         valeValeMomGamma(1:j,j,i) = -valeValeMomGamma(1:j,j,i)
      enddo

      finalStateIndex = 0
      do j = firstFin, lastFin

         ! Increment the finalStateIndex for conjWaveMomSum
         finalStateIndex = finalStateIndex + 1

         do k = 1, valeDim
            conjWaveMomSumGamma(k,finalStateIndex,i) = &
                  & sum(valeValeGamma(:,j,1) * valeValeMomGamma(:,k,i))
         enddo
      enddo
   enddo

   ! Deallocate the momentum matrix elements.
   deallocate (valeValeMomGamma)
#endif


   ! Allocate space for the energy difference.
   allocate (energyDiffTemp (maxPairs))
   allocate (transitionProbTemp  (finComponent,maxPairs))

   ! Initialize the temporary energy transition array.
   energyDiffTemp(:) = 0.0_double

   ! Allocate space to hold the indices for each segment of the energyDiff
   !   array.
   allocate (segmentBorders (lastInit-firstInit+2))

   ! Initialize the first index since it will always be 0.
   segmentBorders(1) = 0

!write (20,*) "firstInit, lastInit=",firstInit,lastInit
!write (20,*) "firstFin, lastFin=",firstFin,lastFin
!write (20,*) "energyCutoff=",energyCutoff
!call flush (20)
   ! Begin the double loop to determine the transition energies.
   do i = firstInit, lastInit
      finalStateIndex = 0
      do j = firstFin, lastFin

         ! Recall that thermal smearing may allow some states to be both
         !   initial and final. We do not consider transitions where the final
         !   state has an energy less than the initial.
         if (i >= j) cycle
!write (20,*) "energyEigenValue(",j,")=",energyEigenValues(j,currentKPoint,spinDirection)
!call flush (20)

         ! If the energy of the final state is higher than the requested
         !   cut-off we go to the next initial state.
         if (energyEigenValues(j,currentKPoint,spinDirection) > &
               & energyCutoff) exit

         ! Compute the energy of the transition from the current states.
         currentEnergyDiff = energyEigenValues(j,currentKPoint,spinDirection)-&
               & energyEigenValues(i,currentKPoint,spinDirection)

!write (20,*) "eDiff,transE",currentEnergyDiff,maxTransEnergy
!call flush (20)
         ! Check if the energy difference is less than the maximum
         !   transition energy that the input file requested computation
         !   for.  If it fails, then we go to the next initial state because
         !   all the remaining final states for this energy will be greater.
         if (currentEnergyDiff > maxTransEnergy) exit

         ! Increment the number of transition pairs counted so far.
         transPairCount = transPairCount + 1
!write (20,*) "transPairCount=",transPairCount,maxPairs
!call flush (20)

         ! Store the transition energy for the current pair.
         energyDiffTemp(transPairCount) = currentEnergyDiff

         ! Increment the final state index for the conjWaveMomSum
         finalStateIndex = finalStateIndex + 1

         ! In the event that thermal smearing is turned on. The state that the
         !   e- comes from and goes into may be fully, partially, or not
         !   occupied. We will scale the probability of a transition linearly
         !   according to the percent occupation of both the initial and final
         !   states.

         ! Determine the array index value of the current initial (index i)
         !   spin-kpoint-state as defined by the tempEnergyEigenValues loop
         !   near the beginning of the population subroutine.
         orderedIndex = i + numStates*(spinDirection-1) + &
               & numStates*spin*(currentKPoint-1)

         ! Use the normal state factor for non-PACS calculations. For PACS
         !   calculations the initStateFactor is always 1 even though the
         !   initial core state(s) will have an electron missing.
         if (stateSet /= 1) then
            initStateFactor = electronPopulation(orderedIndex) / &
                  & (kPointWeight(currentKPoint)/real(spin,double))
         else
            initStateFactor = 1.0_double
         endif

         ! Determine the array index value of the current final (index j)
         !   spin-kpoint-state as defined by the tempEnergyEigenValues loop
         !   near the beginning of the population subroutine.
         orderedIndex = j + numStates*(spinDirection-1) + &
               & numStates*spin*(currentKPoint-1)

         finStateFactor = 1.0_double - electronPopulation(orderedIndex) / &
               & (kPointWeight(currentKPoint)/real(spin,double))
!initStateFactor = 1.0_double
!finStateFactor = 1.0_double
!write (20,*) "i,j,stateFactors",i,j,initStateFactor,finStateFactor
!call flush (20)

#ifndef GAMMA

         ! Loop to obtain the wave function times the momentum integral.
         do k = initComponent,finComponent
             valeValeXMom(k) = sum(valeVale(:,i,1,1) * &
                   & conjWaveMomSum(:,finalStateIndex,k))

            ! Compute the imaginary component of the square of the valeValeXMom
            !   to obtain the transition probability.  Note that the reason it
            !   is imaginary is because of the negative sign included in the
            !   getIntgResults subroutine for the momentum matrix. (I.e. this
            !   is stored as a real number, but it represents the y in (x+iy).)
            transitionProbTemp(k,transPairCount) = &
                  & real(valeValeXMom(k),double)**2+aimag(valeValeXMom(k))**2 *&
                  & initStateFactor*finStateFactor
         enddo
#else

         ! Loop to get the wave function times the momentum matrix element.
         do k = initComponent,finComponent
            valeValeXMomGamma(k) = sum(valeValeGamma(:,i,1) * &
                  & conjWaveMomSumGamma(:,finalStateIndex,k))

            ! Compute the real component of the square of the valeValeXMom to
            !   obtain the transition probability.
            transitionProbTemp(k,transPairCount) = valeValeXMomGamma(k)**2 * &
                  & initStateFactor*finStateFactor
         enddo
#endif
      enddo ! Fin loop j

      ! Save the index for the end border of this segment.
      segmentBorders(i - firstInit + 2) = transPairCount
   enddo ! Init loop i

   ! Deallocate unnecessary matrix
#ifndef GAMMA
   deallocate (conjWaveMomSum)
#else
   deallocate (conjWaveMomSumGamma)
#endif

   ! Determine if there was only one segment.  In this case we don't have to
   !   sort anything.

   ! Sort energyDiffTemp into energyDiff, and obtain the indices for the
   !   correct sorted order of energyDiff so that we can copy the energy
   !   momentum directly.

   allocate (sortOrder (transPairCount))

   call mergeSort (energyDiffTemp,energyDiff(:,currentKPoint,spinDirection),&
         & sortOrder,segmentBorders,transPairCount)

   ! Copy transitionProbTemp to the real transitionProb data structure using
   !   the sorting order determined in the mergeSort subroutine.
   if (xyzComponents == 0) then
      do i = 1, transPairCount
         transitionProb(:,i,currentKPoint,spinDirection) = &
               & transitionProbTemp(:,sortOrder(i))
      enddo
   else
      do i = 1, transPairCount
         transitionProb(xyzComponents,i,currentKPoint,spinDirection) = &
               & transitionProbTemp(1,sortOrder(i))
      enddo
   endif

   ! Deallocate unnecessary arrays and matrices
   deallocate (energyDiffTemp)
   deallocate (transitionProbTemp)
   deallocate (segmentBorders)
   deallocate (sortOrder)

end subroutine computePairs



subroutine computeSigmaE (currentKPoint,xyzComponents,spinDirection)

   ! Import the necessary data modules.
   use O_Kinds
   use O_Constants, only: dim3, pi, auTime, lightFactor, hartree
   use O_KPoints, only: numKPoints, kPointWeight
   use O_SortSubs
   use O_Input, only: maxTransEnSIGE, deltaSIGE, sigmaSIGE
   use O_Lattice, only: realCellVolume
   use O_Potential, only: spin
   use O_AtomicSites, only: valeDim
#ifndef GAMMA
   use O_SecularEquation, only: valeVale, energyEigenValues
#else
   use O_SecularEquation, only: valeValeGamma, energyEigenValues
#endif

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer :: currentKPoint
   integer :: xyzComponents
   integer :: spinDirection

   ! Define local variables.
   integer :: i,j,k ! Loop index variables
   integer :: initComponent
   integer :: finComponent
   integer :: firstInit
   integer :: lastInit
   integer :: firstFin
   integer :: lastFin
   integer :: numInitEnergyStates
   integer :: numFinEnergyStates
   integer :: numTotalEnergyStates
   integer :: finalStateIndex
   integer :: initialStateIndex
   integer :: numEnergyPoints
   real (kind=double) :: alphaFactor
   real (kind=double) :: kPointFactor
   real (kind=double) :: sigmaSqrt2Pi
   real (kind=double) :: conversionFactor
   real (kind=double) :: totalSigma
   real    (kind=double), allocatable, dimension (:,:,:) :: transitionProb
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: conjWaveMomSum
   complex (kind=double)                                 :: valeValeXMom
#else
   real    (kind=double), allocatable, dimension (:,:,:) :: conjWaveMomSumGamma
   real    (kind=double)                                 :: valeValeXMomGamma
#endif

   ! Compute the number of energy points to evaluate.
   numEnergyPoints = int((maxTransEnSIGE * 2.0_double) / deltaSIGE + 1)

   ! Define constants for normalizing the broadening gaussian.
   sigmaSqrt2Pi = sigmaSIGE * sqrt(2.0_double * pi)

   ! Obtain the state bounds for this kpoint.
   firstInit =   firstOccupiedState(currentKPoint,spinDirection)
   lastInit  =    lastOccupiedState(currentKPoint,spinDirection)
   firstFin  = firstUnoccupiedState(currentKPoint,spinDirection)
   lastFin   =  lastUnoccupiedState(currentKPoint,spinDirection)

   ! Get the number of energy states in the system.
   numFinEnergyStates   = lastFin  - firstFin  + 1
   numInitEnergyStates  = lastInit - firstInit + 1
   numTotalEnergyStates = lastFin  - firstInit + 1

   ! Allocate space for the resulting sigmaE and initialize it to zero during
   !   the first kpoint iteration.  Also create the energy scale and initialize
   !   it.  NOTE:  The sigmaEAccumulator
   if (.not. allocated(sigmaEAccumulator)) then
      allocate (sigmaEAccumulator (numEnergyPoints,dim3))
      sigmaEAccumulator (:,:) = 0.0_double

      allocate (energyScale (numEnergyPoints))
      do i = 1, numEnergyPoints
         energyScale(i) = -maxTransEnSIGE + deltaSIGE * (i-1)
      enddo
   endif

   ! Determine the range of components (xyz) that should be considered.
   if (xyzComponents == 0) then
      initComponent = 1
      finComponent = 3
   else
      initComponent = 1
      finComponent = 1
   endif

#ifndef GAMMA
   ! Allocate space to hold the sum(conjg(valeVale(:,j)) * valeVale_Mom(:,k,1))
   !   for each of the possible energy states.
   allocate (conjWaveMomSum (valeDim,numFinEnergyStates,finComponent))

   ! Compute the sum over all the energy states.
   do i = initComponent,finComponent

      finalStateIndex = 0
      do j = firstFin, lastFin

         ! Increment the final state index for conjWaveMomSum. Note that this
         !   will compute the conjWaveMomSumGamma for every possible final
         !   state wihtin the firstFin - lastFin range. Later on, we may find
         !   that for certain initial states we don't actually need every final
         !   state. In those cases we will cycle past the final states.
         finalStateIndex = finalStateIndex + 1

         do k = 1, valeDim
            conjWaveMomSum(k,finalStateIndex,i) = &
                  & sum(conjg(valeVale(:,j,1,1) * valeValeMom(:,k,1,i)))
         enddo
      enddo
   enddo

   ! Deallocate the momentum matrix elements
   deallocate (valeValeMom)

#else

   allocate (conjWaveMomSumGamma (valeDim,numFinEnergyStates,finComponent))

   ! Compute the sum over all the energy states.
   do i = initComponent,finComponent

      finalStateIndex = 0
      do j = firstFin, lastFin

         ! Increment the final state index for conjWaveMomSumGamma. Note that
         !   this will compute the conjWaveMomSumGamma for every possible final
         !   state wihtin the firstFin - lastFin range. Later on, we may find
         !   that for certain initial states we don't actually need every final
         !   state. In those cases we will cycle past the final states.
         finalStateIndex = finalStateIndex + 1

         do k = 1, valeDim
            conjWaveMomSumGamma(k,finalStateIndex,i) = &
                  & sum(valeValeGamma(:,j,1) * valeValeMomGamma(:,k,i))
         enddo
      enddo
   enddo

   ! Deallocate the momentum matrix elements
   deallocate (valeValeMomGamma)
#endif


   ! Allocate space for the fully computed set of transition probabilties
   !   between all of the involved states.  FinComponent is either 1 or 3.
   allocate (transitionProb (finComponent,numTotalEnergyStates,&
         & numTotalEnergyStates))

   ! Begin the double loop to determine the momentum matrix elements between
   !   each possible state.
   ! Initialize the counter for the initial state index number. (The point is
   !   that the firstInit and lastInit values could be something like 345 and
   !   400.  We want to index an array from 1 to 56 so we use counters like
   !   this.)
   initialStateIndex = 0
   do i = firstInit, lastInit

      ! Increment the index number for the initial states.
      initialStateIndex = initialStateIndex + 1

      ! Initialize the counter for the final state index number.
      finalStateIndex = 0
      do j = firstFin, lastFin

         ! Define the indices for the conjWaveMomSum and the momentum matrix
         !   elements. Note that we need to increment the index for every
         !   j-loop iteration because we computed the conjWaveMomSum and
         !   conjWaveMomSumGamma for every final state. We have to increment
         !   the counter for every case, even if we don't use it.  (This note
         !   is here because this was a point of confusion in the past when the
         !   code was being developed.)
         finalStateIndex = finalStateIndex + 1

         ! We don't want double counting or self interactions so we skip the
         !   i>=j cases.
         if (i >= j) cycle

#ifndef GAMMA
         ! Loop to obtain the wave function times the momentum integral.
         do k = initComponent,finComponent
             valeValeXMom = sum(valeVale(:,i,1,1) * &
                   & conjWaveMomSum(:,finalStateIndex,k))

            ! Compute the imaginary component of the square of the
            !   valeValeXMom to obtain the transition probabilty element.
            ! Note that the reason it is imaginary is because of the negative
            !   sign included in the getIntgResults subroutine for the
            !   momentum matrix. (See notes in that code.)
            transitionProb(k,initialStateIndex,finalStateIndex) = &
                  &  real(valeValeXMom,double)**2 + &
                  & aimag(valeValeXMom)**2
         enddo
#else

         ! Loop to obtain the wave function times the momentum integral.
         do k = initComponent, finComponent
            valeValeXMomGamma = sum(valeValeGamma(:,i,1) * &
                  & conjWaveMomSumGamma(:,finalStateIndex,k))

            ! Compute the real component of the square of the valeValeXMom to
            !   obtain the transition probability.
            transitionProb(k,initialStateIndex,finalStateIndex) = &
                  & valeValeXMomGamma**2
         enddo
#endif
      enddo
   enddo

   ! Determine the weighting effect of this kpoint and include the normalizaion
   !   factor for the gaussian. We will actually be multiplying two Gaussians
   !   together so we need this squared.
   kPointFactor = kPointWeight(currentKPoint)/real(spin,double) / &
         & (sigmaSqrt2Pi)**2

   ! Now compute the exponential alpha factor which is -1/(2 * sigma^2).
   alphaFactor = -1.0_double / (2.0_double * sigmaSIGE**2)

   ! Now we fill up the sigmaEAccumulator.  There are two scenarios in which
   !   this may happen, the xyz all at once (xyzComponents==0) case, and the
   !   each axis (x, y, z) one at a time case (xyzComponents/=0). The only
   !   real difference is the k-loop in the first case and the specific index
   !   access in the second.
   ! The basic approach is to loop over all pairs of states and for each pair
   !   create a Gaussian (of unit area) for each, multiply them together to
   !   get a new Gaussian between them, and then numerically evaluate and
   !   accumulate those on a mesh.
   if (xyzComponents == 0) then

      ! Initialize the counter for the index number of the initial states.
      initialStateIndex = 0

      ! Loop over the set of initial states from first to last. Recall that
      !   this list will include fully occupied states a little below the first
      !   unoccupied state up to the last state with any occupation (as long
      !   as it is within range of the fermi level).
      do i = firstInit, lastInit

         ! Increment the index number for the initial states.
         initialStateIndex = initialStateIndex + 1

         ! Initialize the counter for the index number of the final states.
         finalStateIndex = 0
         do j = firstFin, lastFin

            ! Increment the index number for the final states. Note that we
            !   need to increment the index for every j-loop iteration because
            !   we computed the conjWaveMomSum and conjWaveMomSumGamma for
            !   every final state. We have to increment the counter for every
            !   case, even if we don't use it.  (This note is here because this
            !   was a point of confusion in the past when the code was being
            !   developed.)
            finalStateIndex = finalStateIndex + 1

            ! We don't want double counting or self interactions so we skip the
            !   i>=j cases.
            if (i >= j) cycle

            ! Here we compute the product of two gaussians evaluated on a mesh.
            !   The kPointFactor is the square of (the kpoint weight divided by
            !   the spin (2 or 1 depending on spin polarized or not), divided
            !   by sqrt(2Pi)*sigmaSIGE broadening factor). The exponential has
            !   the alphaFactor of -1/(2*sigmaSIGE^2) times the r^2 value where
            !   the r^2 value is essentially the distance between the current
            !   mesh point and the center of the new Gaussian that is formed
            !   by the product of the Gaussians for each transition state. The
            !   r^2 can be expressed as r1^2 +r2^2 where r1=e-e1 and r2=e-e2
            !   with e=the current mesh point energy and e1=the current initial
            !   state energy value and e2=the current final state energy value.
            do k = 1, 3
               sigmaEAccumulator(:,k) = sigmaEAccumulator(:,k) + &
                     & kPointFactor * exp(alphaFactor * ( &
                     & (energyScale(:) - energyEigenValues(i,currentKPoint,&
                     & spinDirection))**2 + &
                     & (energyScale(:) - energyEigenValues(j,currentKPoint,&
                     & spinDirection))**2)) * transitionProb(k,&
                     & initialStateIndex,finalStateIndex)
            enddo
         enddo
      enddo
   else

      ! Now we create a Gaussin of unit area for each pair of states and
      !   multiply them together to create a new Gaussian (somewhere between
      !   them). These product Gaussians will be evaluated and accumulated on a
      !   numerical mesh.
      ! Initialize the counter for the index number of the initial states.
      initialStateIndex = 0
      do i = firstInit, lastInit

         ! Increment the index number for the initial states.
         initialStateIndex = initialStateIndex + 1

         ! Initialize the counter for the index number of the final states.
         finalStateIndex = 0
         do j = firstFin, lastFin

            ! We don't want double counting or self interactions so we skip the
            !   i>=j cases.
            if (i >= j) cycle

            ! Increment the index number for the final states. Note that we
            !   need to increment the index for every j-loop iteration because
            !   we computed the conjWaveMomSum and conjWaveMomSumGamma for
            !   every final state. We have to increment the counter for every
            !   case, even if we don't use it.  (This note is here because this
            !   was a point of confusion in the past when the code was being
            !   developed.)
            finalStateIndex = finalStateIndex + 1

            ! See notes for the above case.
            sigmaEAccumulator(:,xyzComponents) = &
                  & sigmaEAccumulator(:,xyzComponents) + &
                  & kPointFactor * exp(alphaFactor * ( &
                  & (energyScale(:) - energyEigenValues(i,currentKPoint,&
                  & spinDirection))**2 + &
                  & (energyScale(:) - energyEigenValues(j,currentKPoint,&
                  & spinDirection))**2)) * transitionProb(1,&
                  & initialStateIndex,finalStateIndex)
         enddo
      enddo
   endif

   ! Print the results during the last KPoint iteration for the last component
   !   or set of components.
   if ((currentKPoint == numKPoints) .and. ((xyzComponents == 0) .or. &
         & (xyzComponents == 3))) then

      ! This is a bit tricky and so it will help a lot to know something
      !   about cgs esu units before you start.  The Kubo-Greenwood formula
      !   (KGF) given in this reference:  Greenwood, Proc. Phys. Soc. London,
      !   71,585, (1958) is in cgs esu units.  We did the calculation in a.u.
      !   and the result must be presented in SI.  Arrrg.  Note some deviations
      !   from the formula as presented and as we calculated.  (1)  The
      !   formula used the velocity operator while this calculation used the
      !   momentum operator, so our equations have the electron mass squared in
      !   the denominator to compensate.  (2)  The formula should be considered
      !   as a vector outer product that produces a rank-2 tensor.  Our
      !   computed result is an x, y, z vector that is then averaged when it is
      !   printed out.  That is where the 1/3 factor comes from in our equation.

      ! The equation used is as follows:
      !    2 * Pi * hbar * e^2 / (3 * m^2 * Omega) *
      !    Sum ( |Pnm|^2 * delta(E-En) * delta(E-Em))

      ! In a.u. hbar = 1, e = 1, and m = 1 (They still have units though)

      ! In a.u. we have E=energy, T=time, M=mass, L=length.

      ! The first step is to apply the terms in the coefficient that are not
      !   equal to one.  (Note that we do not divide by three here since we
      !   will perform the averaging later.  We also do not change the units of
      !   the cell volume since that will be accounted for next.)
      conversionFactor = 2.0_double * pi / realCellVolume


      ! The units of the KGF in a.u. are seen as:
      ! SigmaE = (ET EL) / (M^2 L^3) * (M^2 L^2)/(T^2) * 1/(E^2)
      ! SigmaE = 1/T
      ! hBar is ET, charge^2 (Q^2) is ML^3/T^2 = EL, we divide by the volume
      !   explicitly unlike in the KGF.  The tricky part is the charge because
      !   it is not a separate fundamental unit in cgs.  In SI it is in a
      !   separate fundamental unit called Coulombs, but in cgs it is not.  In
      !   cgs the units of Q are sqrt(g cm^3 / s^2).

      ! The second step is then to convert the 1/T in a.u. to 1/sec in cgs.
      !   This can be done by the conversion factor of 2.418884326505x10-17 s
      !   equals one atomic unit of time taken from NIST (auTime).  Note that
      !   the factor of 1e-17 is not included in the constant auTime and it
      !   must be accounted for later.

      ! This will convert the result from (a.u. 1/T) to (cgs emu 1/s).  Note
      !   that it still needs to be multiplied by 1e+17.
      conversionFactor = conversionFactor / auTime

      ! The third step is then to convert the 1/sec in cgs to 1/(ohm m) in
      !   SI.  This can be done with the knowledge that ohm = ET/Q^2 = Js/C^2.
      !   Writing that in cgs esu we have:  (g cm^2 / s^2) s s^2/(g cm^3) = 
      !   s / cm.  Then 1/(ohm cm) = cm / (s cm) = 1/s in cgs units.  This
      !   shows that the two units are equivalent.  Now we just need the
      !   conversion factor between them.
      !   We have from Jackson:  1 / (ohm m) = (1e-16 * c^2) * 1e9 1/s
      !   Therefor: 1/s (cgs emu) = 1/(1e-16 * c^2) * 1e-9 1/(ohm m)
      !   The value of 1e-16 * c^2 is DEFINED as:  8.9875517873681764.

      ! This will convert the result from (cgs emu 1/s) to (SI 1/(ohm m)).
      !   Note that it still needs to be multiplied by an additional 1e-9.
      !   This creates a total exponent multiplication factor of 1e+8 so far.
      conversionFactor = conversionFactor / lightFactor

      ! Finally, we want to express the result in (micro ohm cm)^-1 which is
      !   equal to (1e-8 ohm m)^-1 = 1e8(ohm m)^-1.  This creates a total
      !   exponent multiplication factor of 1 so we don't have to make any more
      !   modifications.

      ! Adjust the sigmaE to have the correct units.
      sigmaEAccumulator(:,:) = sigmaEAccumulator(:,:) * conversionFactor

      do i = 1, numEnergyPoints

         ! Write the averaged total, and the x, y, z components, making sure to
         !   convert the energy scale into eV.
         write (49+spinDirection,fmt="(5e20.8e4)") energyScale(i)*hartree,&
            & sum(sigmaEAccumulator(i,:))/3.0_double,sigmaEAccumulator(i,:)
      enddo

      ! Accumulate the total electronic contribution to the thermal conductivity.
      totalSigma = sum(sigmaEAccumulator(:,:)) / 3.0_double
      write (20,*) "The total sigma is: ",totalSigma

      deallocate (sigmaEAccumulator)
      deallocate (energyScale)
   endif

   ! Deallocate unnecessary matrices.
#ifndef GAMMA
   deallocate (conjWaveMomSum)
#else
   deallocate (conjWaveMomSumGamma)
#endif
   deallocate (transitionProb)

end subroutine computeSigmaE


end module O_OptcTransitions
