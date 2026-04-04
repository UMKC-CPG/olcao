module O_Populate

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define module data.
   integer :: coreStructInit ! Flag for whether or not the core state
         ! structures have been initialized in an XANES/ELNES type calculation.
   integer :: occupiedEnergyIndex ! The index value of the highest occupied
         ! state in the sorted list.
   real (kind=double) :: occupiedEnergy ! The energy (in a.u.) that is the
         ! highest occupied state.
   real (kind=double), allocatable, dimension (:) :: sortedEnergyEigenValues
         ! This array holds all the energy eigen values of all kpoints and all
         ! spin states in sorted order from lowest to highest.
   real (kind=double), allocatable, dimension (:) :: electronPopulation
         ! This array holds the number of electrons in each of the kpoint,
         !   spin states of the original unsorted energy eigen values.
   real (kind=double), allocatable, dimension (:,:,:) :: &
         & electronPopulation_LAT
         ! LAT (Linear Analytic Tetrahedron) analog of electronPopulation.
         !   Each entry gives the fractional electron occupation of state
         !   (band, kpoint, spin) as determined by the Bloechl tetrahedron
         !   method (PRB 49, 16223, 1994). This replaces electronPopulation
         !   when the LAT integration method is active (kPointIntgCode==1).
         !   While electronPopulation uses Gaussian broadening and Fermi
         !   filling to assign weights, this array uses analytic corner
         !   integration weights that decompose each tetrahedron's occupied
         !   volume among its four corner k-points. The result is a
         !   broadening-parameter-free occupation weight for each state.
         !   Dimensions: (numStates, numKPoints, spin).
         !   Allocated and filled by computeElectronPopulation_LAT.
   integer, allocatable, dimension (:) :: indexEnergyEigenValues
         ! This array holds a mapping between the original energy eigen values
         ! and a sorted list of the energy eigen values.  The index number in
         ! this array corresponds to the index number of the sorted energy
         ! eigen value array.  The value in this array corresponds to the index
         ! of the original unsorted array.  (Note that while the original
         ! array is not sorted, it is grouped by kpoint, spin, and then state.
         ! Thus, given a kpoint, spin, and state it is possible to deduce the
         ! appropriate index number to use.)
   integer, dimension(6,4) :: QN_nlOrderIndex ! The first index is the QN_n,
         ! and the second index is the QN_l+1 (because Fortran starts counting
         ! array indices at 1).  The value is the order of the state from
         ! lowest to higher energy in the ideal atom case.  The mixed filling
         ! of various states is not considered here (but it should be in the
         ! future especially since some basis set atoms may not fill in the
         ! ideal way.
   integer, dimension(14) :: numOrbitalStates ! For each of the above
         ! QN_nlOrderIndex values this records how many states there are for
         ! that orbital (including spin).
   integer :: excitedCoreStateIndex ! The band index number of the core state
         ! that is to be excited if a ELNES/XANES type calculation is being
         ! done.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


! The population is done in two phases.  The first phase uses the old
!   standard technique and does not apply any thermal smearing to the
!   populated levels.  This phase is done ONLY to obtain the minimum
!   energy level that should be populated (which is trivial), and the
!   maximum energy level that should be populated (which is not trivial in
!   multi kpoint and excited state environments).

! The basic problem of doing spin polarized populations and also multi-
!   kpoint populations is that you don't know what the order of states is.
!   Is the lowest energy from kpoint #1 or kpoint #2?  How about the second
!   or third lowest?  Is the spin up from state #1 of kpoint #1 higher
!   energy than the spin down from state #2 of kpoint #1?  What about mixed
!   kpoints?

! The solution used here is to sort all the energy states and keep track of
!   where they came from.

! 1)  Split the energy levels on the first iteration of the SCF cycle for
!     spin polarized calculations only.
! 2)  Copy all the energy values into a 1 dimensional temp array.  The order
!     of the indices is important.  The values are put in by state first,
!     then spin, and finally kpoints.  So for a 4 state, 2 spin, 2 kpoint
!     system we would have:  111, 211, 311, 411, 121, 221, 321, 421, 112,
!     212, 312, 412, 122, 222, 322, 422.  Where the first digit of the
!     triplet is the state, the second is the spin, and the third is the
!     kpoint.   In the linear array the values are grouped first by their
!     kpoint (e.g. the first 8 are kpoint 1, and the second 8 are kpoint 2.)
!     Then, within the kpoint groups the values are grouped by spin.  (e.g.
!     the first 4 are spin up and the second 4 are spin down.)  Finally,
!     within the spin groups, the values are ordered by their state number.
!     This list is grouped, but not sorted.
! 3)  Sort the energy eigen values into a new array:
!     sortedEnergyEigenValues.  A mapping is created between the sorted
!     values and the original values grouped by kpoints, then spin, and then
!     state.  The mapping is held in indexEnergyEigenValues.  The index
!     number for this array corresponds to the index number for the sorted
!     array.  The value in the array is the index of the original unsorted
!     array.  (e.g. sortedArray = (2, 3, 6, 8) originalArray = (3, 6, 8, 2),
!     indexArray(1) = 4, indexArray(2) = 1, indexArray(3) = 2, and
!     indexArray(4) = 3.
! 4)  Now, the electron population for each state/kpoint/spin set is
!     computed.  The electron population array is one dimensional and the
!     indices for this array match the indices for the original unsorted
!     energy eigen value array.  The electron population for each triplet is
!     is then assigned IN THE ORDER OF THE INDEX ARRAY.  This will make the
!     lowest energy triplit the first to be populated, and the second lowest
!     triplet the next to be populated, etc.  ONE CATCH!  The amount of
!     electron popultion to go into a given triplet depends on which kpoint
!     that triplet includes.  FORTUNATELY, we thought ahead and grouped the
!     values according to their kpoint.  So, we just have to do a simple
!     integer division to get the kpoint number:
!     e.g. 1+(indexArray(i)-1) / (numStates*spin).  Once the cumulative
!     number of electrons is the same as the number of electrons in the
!     system, we abort the loop.
! 5)  There is one issue with the highest occupied state.  If it turns out
!     that this state is degenerate, then we need to distribute the charge
!     in all degenerate states evenly.  This will prevent the case of some
!     states being occupied while others at the same energy level are
!     totally unoccupied.  We find which states have the same energy and
!     then distribute the charge evenly among them by summing the total
!     charge that exists in the degenerate states and summing the total
!     charge that could be put in the degenerate states if they were all
!     completely filled.  Then, each state gets an amount of charge equal to
!     the kPoint weight of that state times the ratio of total existing over
!     total possible charge.  (NOTE that this scheme is NOT applied to the
!     thermal smearing phase.)
! 6)  The thermal smearing scheme is basically the same except that the
!     population is modified by the fermi function.
! 7)  The case of xanes calculation is also basically the same except that
!     we consider the system to have one extra electron, and then we remove
!     an electron from a core orbital after the system was completely
!     populated.


subroutine populateStates

   ! Import necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Constants, only: hartree
   use O_Input, only: thermalSigma
   use O_CommandLine, only: excitedQN_n

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Log the date and time we start.
   call timeStampStart (16)

   ! In the case that this is a XANES/ELNES type of calculation then we must
   !   initialize the core state structures if they have not yet been set up.
   if ((excitedQN_n /= 0) .and. (coreStructInit == 0)) then
      call initCoreStateStructures
   endif

   ! Always populate the states in the standard way.
   call populateStandard

   ! If the thermal smearing parameter is non-zero, then we continue the
   !   population using the thermal smearing scheme.
   if (thermalSigma /= 0.0_double) then
      call populateSmearing
   endif

   ! Write out the calculated occupied energy (Fermi energy for metals) in eV.
   write (20,*) 'Highest Occupied Energy=',occupiedEnergy * hartree,' eV'
   call flush (20)

   ! Log the date and time we end.
   call timeStampEnd (16)

end subroutine populateStates


subroutine populateStandard

   ! Import necessary modules.
   use O_Kinds
   use O_Potential,       only: spin
   use O_Constants,       only: smallThresh
   use O_SortSubs,        only: mergeSort
   use O_CommandLine,     only: excitedQN_n
   use O_SecularEquation, only: energyEigenValues
   use O_Input,           only: numStates, numElectrons
   use O_KPoints,         only: numKPoints, kPointWeight

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the local variables used in this subroutine.
   integer :: i,j,k ! Loop index variables
   integer :: tripletCounter ! Counts the total number of states in the system
         ! and is primarily used for tracking the index of various arrays.
   integer :: numSegments ! Records the number of segments in the energy


   ! Define variables used to populate the electron levels
   integer, allocatable, dimension (:) :: segmentBorders
   real (kind=double) :: populatedElectrons ! Tracks the total number of
         ! electrons (including fractional components) that have been assigned
         ! to energy levels ordered from lowest to highest over all kpoints.
   real (kind=double), allocatable, dimension (:) :: tempEnergyEigenValues

   ! Define variables to evenly fill the highest occupied degenerate states.
   integer :: minStateIndex
   integer :: maxStateIndex
   real (kind=double) :: possibleDegenCharge
   real (kind=double) :: degenerateCharge


   ! Allocate arrays and matrices for populating the electron levels
   ! These arrays will be deallocated in makeValenceRho after the kpoint
   !   loop because they are not needed until the next scf iteration after that
   !   point.
   allocate (electronPopulation      (numStates*numKPoints*spin))
   allocate (sortedEnergyEigenValues (numStates*numKPoints*spin))
   allocate (indexEnergyEigenValues  (numStates*numKPoints*spin))


   ! These arrays are used only for sorting and will be deallocated here.
   allocate (segmentBorders        (numStates*numKPoints*spin+1))
   allocate (tempEnergyEigenValues (numStates*numKPoints*spin))

   ! Sort the energy eigen values in ascending order for populating.  This is
   !   done once at the beginning of this subroutine and the sorted order is
   !   used to populate the levels for both phases.

   ! Copy the energyEigenValues into the tempEnergyEigenValues single dimension
   !   array.  It is important to have the kpoints be the outer loop even if it
   !   is slower because we need the energyEigenValues sorted into groups based
   !   on their kpoint to assign the electron population easily later.
   tripletCounter = 0
   do i = 1, numKPoints
      do j = 1, spin
         do k = 1, numStates
            tripletCounter = tripletCounter + 1
            tempEnergyEigenValues(tripletCounter) = energyEigenValues(k,i,j)
         enddo
      enddo
   enddo


   ! Initialize the segment borders to each individual array slot.
   numSegments = numStates * numKPoints * spin
   segmentBorders(1) = 0
   do i = 1, numSegments
      segmentBorders (i+1) = i
   enddo


   ! Call the sorting subroutine
   call mergeSort (tempEnergyEigenValues,sortedEnergyEigenValues,&
         & indexEnergyEigenValues,segmentBorders,numSegments)

   ! Deallocate the unnecessary arrays from the sorting procedure.
   deallocate (segmentBorders)
   deallocate (tempEnergyEigenValues)


   ! Initialize the electron population of each level to 0.
   electronPopulation(:) = 0.0_double

   ! Initialize local variables
   populatedElectrons = 0.0_double ! Total number of populated electrons.


   ! The basic method is to fill the energy levels in order from the lowest
   !   energy level to the highest.  This is simple since we have just sorted
   !   the levels from lowest to highest.  Now we just fill them in with one
   !   exception.

   ! In the case where the level number is equal to the level number that has
   !   been depopulated by one electron due to xanes excitation.  The simplist
   !   method to deal with xanes calculations here is to populate all the
   !   levels in order as described above except consider the numElectron
   !   number to be larger by 1 than the real numElectron number.  Then we
   !   remove one electron worth of population value from each kpoint at the
   !   core excitation level.  This will put the excited electron in the lowest
   !   conduction band state and leave a core hole.
   

   ! First we determine if it is necessary to artificially increment the
   !   numElectrons in the system by one due to xanes excitation.
   if (excitedQN_n /= 0) then
      numElectrons = numElectrons + 1
   endif


   ! Consider each energy value and determine which KPoint it came from.  Then
   !   assign the appropriately weighted number of electrons to that state.
   do i = 1, numSegments

      ! Weight the electron population for this state by the kpoint weighting
      !   factor and divide by the spin value so that a spin-polarized
      !   calculation has one electron per state while a non-spin-polarized
      !   calculation has two electrons per state.  The value in the (x/y) will
      !   be a number from 0 to numKPoints-1 so we add one to that result to
      !   get the index for the kPointWeight.  The reason this works is because
      !   the energy values were grouped by their kpoints before sorting.
      electronPopulation(indexEnergyEigenValues(i)) = kPointWeight &
            & (1+(indexEnergyEigenValues(i)-1)/(numStates*spin)) / &
            & real(spin,double)

      ! Accumulate the electron population
      populatedElectrons = populatedElectrons + &
            & electronPopulation(indexEnergyEigenValues(i))

      ! Determine if the number of electrons that have been populated
      !   exceeds the number of electrons in the system.  If so, reduce the
      !   last populated level so that the total number of populated electrons
      !   equals the total number of electrons in the system.
      if (populatedElectrons > real(numElectrons,double)) then

         ! First adjust the actual population.
         electronPopulation(indexEnergyEigenValues(i)) = kPointWeight &
            & (1+(indexEnergyEigenValues(i)-1)/(numStates*spin)) / &
            & real(spin,double) + real(numElectrons,double) - &
            & populatedElectrons

         ! Then adjust the record of the number of populated electrons.
         populatedElectrons = real(numElectrons,double)

      endif

      ! Abort the loop when the last electron has been populated, recording
      !   the occupied energy (fermi energy for metals) as we leave.
      if (abs(real(numElectrons,double) - populatedElectrons) < &
            & smallThresh) then
         occupiedEnergyIndex = i
         occupiedEnergy = sortedEnergyEigenValues(occupiedEnergyIndex)
         exit
      endif
   enddo

   ! Check for degeneracy of the highest occupied energy level.  If it is
   !   degenerate we need to find all the states that have the same energy and
   !   then distribute the electron population evenly to all of them.  If we
   !   do not do this in, for example, an isolated Fe atom then some of the d
   !   electron states would be occupied and others would not and the charge
   !   distribution would be very distorted.  The charge would never settle
   !   down and would keep flushing back and forth between the different d
   !   angular orientation sub-orbitals.  This will help with the convergence
   !   of metals too.


   ! Initialize to the case of no degenerate states.
   minStateIndex = occupiedEnergyIndex
   maxStateIndex = occupiedEnergyIndex
   degenerateCharge = electronPopulation(indexEnergyEigenValues( &
         & occupiedEnergyIndex))
   possibleDegenCharge = kPointWeight(1+(indexEnergyEigenValues &
         & (occupiedEnergyIndex)-1)/(numStates*spin)) / real(spin,double)

   ! The way to check for degeneracy is to first assume that the last populated
   !   state is not degenerate and then search higher and lower energy states
   !   to see if they are at a similar energy.
   do while (abs(occupiedEnergy - sortedEnergyEigenValues(maxStateIndex+1)) <= &
         & smallThresh)

      ! Increment the index counter.
      maxStateIndex = maxStateIndex + 1

      ! Accumulate the amount of charge that *could* be held in the degenerate
      !   states.  (The possible charge.)
      possibleDegenCharge = possibleDegenCharge + kPointWeight &
            & (1+(indexEnergyEigenValues(maxStateIndex)-1) / &
            & (numStates*spin)) / real(spin,double)

      ! Abort if we reach the bounds of the calculation.
      if (maxStateIndex == numKPoints * numStates * spin) exit
   enddo

   ! If we are doing a single h atom, then minStateIndex will equal 1.
   if (minStateIndex > 1) then
      do while (abs(occupiedEnergy &
            & - sortedEnergyEigenValues(minStateIndex-1)) <= smallThresh)

         ! Decrement the index counter.
         minStateIndex = minStateIndex - 1
   
         ! Accumulate the amount of charge that *could* be held in the degenerate
         !   states.  (The possible charge.)
         possibleDegenCharge = possibleDegenCharge + kPointWeight &
               & (1+(indexEnergyEigenValues(minStateIndex)-1) / &
               & (numStates*spin)) / real(spin,double)
   
         ! Accumulate the amount of existing charge.
         degenerateCharge = degenerateCharge + kPointWeight &
               & (1+(indexEnergyEigenValues(minStateIndex)-1) / &
               & (numStates*spin)) / real(spin,double)
   
          ! Abort if we reach the bounds of the calculation.
         if (minStateIndex == 1) exit
      enddo
   endif

   
   ! For each of the degenerate states redistribute the charge evenly.
   do i = minStateIndex, maxStateIndex
      electronPopulation(indexEnergyEigenValues(i)) = degenerateCharge / &
            & possibleDegenCharge * kPointWeight &
            & (1+(indexEnergyEigenValues(i)-1)/(numStates*spin)) / &
            & real(spin,double)

   enddo

   ! Adjust the occupied energy and occupied energy index.
   occupiedEnergyIndex = maxStateIndex
   occupiedEnergy = sortedEnergyEigenValues(occupiedEnergyIndex)


   ! Now that all the levels have been populated we remove an electron from
   !   the level that had a xanes excitation applied to it (if applicable).
   if (excitedQN_n /= 0) then
      call correctCorePopulation
   endif

end subroutine populateStandard


subroutine populateSmearing

   ! Import necessary modules.
   use O_Kinds
   use O_Constants,   only: smallThresh
   use O_Potential,   only: spin
   use O_CommandLine, only: excitedQN_n
   use O_MathSubs,    only: stepFunction
   use O_KPoints,     only: numKPoints, kPointWeight
   use O_Input,       only: numStates, numElectrons, fermiSearchLimit, &
         & thermalSigma, maxThermalRange

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the local variables used in this subroutine.
   integer :: i,j
   integer :: numSegments
   real (kind=double) :: fermiFunction
   real (kind=double) :: minEnergy
   real (kind=double) :: maxEnergy
   real (kind=double) :: chargeSum

   ! Compute the number of energy eigen value terms.
   numSegments = numStates * numKPoints * spin

   ! So the basic idea here is to make two posts that bound the possible
   !   values of the fermi energy. The first guess for the fermi energy is in
   !   the exact middle between the two posts. (Which will be in the middle of
   !   the gap for an insulator.) Then we apply the smearing function and if
   !   the integrated charge is too small we bring the minimum post up to the
   !   current guess and if the integrated charge is too large we bring the
   !   maximum post down to the current guess. Then we make a new guess between
   !   the posts and evaluate the smearing function all over again.

   ! The first guess will be pretty close to correct for most systems. In the
   !   case of most insulators with a gap this will put the fermi level right
   !   in the middle of the gap and it will be correct on the first iteration.
   !   In the case of metals or small gap semiconductors it may be close but
   !   wrong on the first iteration. Then one border post will shift a lot (by
   !   moving to the point of the zeroth fermi level guess) and the integrated
   !   charge will be wrong for a few iterations while the two border posts
   !   move closer to each other.

   ! The process may be envisioned in the following sequence of interations in
   !   the loop over variable i. Consider that the first guess integrates to a
   !   charge that is too high. The high post moves in and the low post stays.
   !   The new guess (1) will probably be way too low (the integrated charge
   !   will be < the number of electrons). The low post moves up and the new
   !   guess is made at (2) for the fermi level. In the figures we assume that
   !   each of the subsequent guesses are too low (producing an integrated
   !   charge that is too small). Thus, steps 3 and 4 are moves of the lower
   !   post. Eventually, it will find a situation that is neither too high or
   !   too low. This is the appropriate fermi level.
   !   |                    0                   |
   !   |          1         |
   !              |    2    |
   !                   |  3 |
   !                      |4|

   ! The Fermi search range is given in the input (olcao.dat) file with a
   !   default value of 13.6 eV (which is ~0.5 a.u.). This should realistically
   !   encompass any appropriate search space for a fermi level because it is
   !   13.6 eV lower than the highest occupied state and 13.6 eV higher than
   !   the lowest unoccupied state. (By this time the value is stored in a.u.)

   ! Initialize the min and max boundaries acording to the given thermal
   !   smearing range and the highest occupied and lowest unoccupied states.
   !   That is, the min boundary starts at the highest occupied minus some
   !   defined smearing range while the max boundary starts at the lowest
   !   occupied plus the same smearing range. (Thus the initial energy for the
   !   fermi level will be the fermi energy (occupiedEnergy) determined above
   !   in the non-smearing population subroutine that is always performed
   !   before any smearing is done.
   minEnergy = sortedEnergyEigenValues(occupiedEnergyIndex)   - &
         & fermiSearchLimit
   maxEnergy = sortedEnergyEigenValues(occupiedEnergyIndex+1) + &
         & fermiSearchLimit

   ! Begin a search for the Fermi energy.  500 iterations should be enough. :)
   do i = 1,500

      ! Now we have to repeat something similar to the first phase.  We will
      !   fill the energy levels in order again.

      ! So first we determine if it is necessary to artificially increment the
      !   numElectrons in the system by one due to xanes excitation.
      if (excitedQN_n /= 0) then
         numElectrons = numElectrons + 1
      endif

      ! Make an initial guess for the Fermi energy (occupiedEnergy) as being
      !   between the two border points.
      occupiedEnergy = 0.5_double * (minEnergy + maxEnergy)

      ! Initialize the counter of the number of electrons populated so far.
      chargeSum = 0.0_double

      ! Initialize the electron population.
      electronPopulation (:) = 0.0_double

      ! Start populating every state for every kpoint.
      do j = 1, numSegments

         ! Compute fermi function multiplicative factor for this energy value.
         fermiFunction = stepFunction((sortedEnergyEigenValues(j) - &
               & occupiedEnergy) / thermalSigma,maxThermalRange)

         ! Abort in the event that the fermi function is zero because all the
         !   remaining energy values will also be zero.
         if (fermiFunction == 0.0_double) then
            exit
         endif

         ! Weight the electron population for this state by the kpoint 
         !   weighting factor and divide by the spin value so that a
         !   spin-polarized calculation has one electron per state while a
         !   non-spin-polarized calculation has two electrons per state.  Then
         !   multiply by a smeared step function to smear the population.
         electronPopulation(indexEnergyEigenValues(j)) = (kPointWeight &
               & (1+(indexEnergyEigenValues(j)-1)/(numStates*spin)) / &
               & real(spin,double)) * fermiFunction

         ! Accumulate the electron population.
         chargeSum = chargeSum + electronPopulation(indexEnergyEigenValues(j))

      enddo

      ! Now that all the levels have been populated we remove an electron from
      !   the core level that had an excitation applied to it (if applicable).
      !   We must also reduce the chargeSum by one because it also counted that
      !   one core electron.
      if (excitedQN_n /= 0) then  ! If not a ground state type calculation.
         call correctCorePopulation
         chargeSum = chargeSum - 1
      endif

      ! Record the integrated charge for posterity.
!      write (20,fmt="(a8,e16.8,a10,i4,a11,e16.8)") "IntgChg=",chargeSum,&
!            & " for iter ",i," occEnergy=",occupiedEnergy*hartree
!      call flush (20)

      ! If the difference in charge between the correct numElectrons and
      !   the just assigned chargeSum is sufficiently small then we can exit
      !   the population outer (i) loop since we have found the correct
      !   population values, and Fermi level.
      if (abs(chargeSum - real(numElectrons,double)) < smallThresh) then
         exit
      endif

      ! Move one or the other border value to the position of the last
      !   Fermi energy.  The Fermi energy will be re-guessed in the next
      !   iteration.
      if (chargeSum < real(numElectrons,double)) then
         minEnergy = occupiedEnergy
      else
         maxEnergy = occupiedEnergy
      endif
   enddo

end subroutine populateSmearing


subroutine correctCorePopulation

   ! Import necessary modules.
   use O_Kinds
   use O_Potential,   only: spin
   use O_KPoints,     only: numKPoints
   use O_Input,       only: numElectrons
   use O_CommandLine, only: excitedQN_n, excitedQN_l

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the local variables used in this subroutine.
   integer :: i
   integer :: initialState
   integer :: finalState

   ! Core occupation should proceed in the following order:
!   !   1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, 4d, 5p, 6s, 4f
   !   1s, 2s, 2p, 3s, 3p, 3d, 4s, 4p, 4d, 4f, 5s, 5p, 5d, 6s
   ! This was assumed when assigning the excitedCoreStateIndex and the values
   !   in the numOrbitalStates(:,:) array.

   ! The initial state must begin at the next state after all the states lower
   !   than the excitedCoreStateIndex.  (excitedCoreStateIndex-1) is the band
   !   index of the state just below the excited core state.  The numKPoints
   !   will get us past all copies of this state at each kPoint.  The spin will
   !   get us past the spin degenerate core states (this will not be true if
   !   we include spin-orbit coupling so it will be necessary to devise a
   !   better way to select the states to be excited).  The +1 is then the
   !   next higher state.
   initialState = numKPoints * (excitedCoreStateIndex-1) * spin + 1

   ! The final state must be the last of the states of the same orbital type as
   !   the initial state.  We simply add to the initial state the number of
   !   states in this orbital type and multiply by the kpoint and spin factors.
   !   The -1 is the last state (without it we would mark the next higher
   !   orbital).
   finalState = initialState + numKPoints * spin * &
         & numOrbitalStates(QN_nlOrderIndex(excitedQN_n,excitedQN_l+1)) - 1

   ! Remove a fraction of an electron from every kpoint/spin pair for the
   !   choice of QN_n, QN_l states.
   do i = initialState, finalState
      electronPopulation(indexEnergyEigenValues(i)) = &
            & electronPopulation(indexEnergyEigenValues(i)) - &
            & electronPopulation(indexEnergyEigenValues(i)) * &
            & real(spin,double) / 2.0_double / &
            & numOrbitalStates(QN_nlOrderIndex(excitedQN_n,excitedQN_l+1))
   enddo

   ! We also restore the local numElectrons variable for the next part of
   !   the OLCAO calculation that may rely on it.
   numElectrons = numElectrons - 1

end subroutine correctCorePopulation


subroutine initCoreStateStructures

   ! Import necessary modules.
   use O_Kinds
   use O_Potential, only: spin
   use O_CommandLine, only: excitedQN_n, excitedQN_l

   ! Make sure that nothing funny is accidentally used.
   implicit none

   ! Define local variables.
   integer :: i

   ! Core occupation should proceed in the following order:
!   !   1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, 4d, 5p, 6s, 4f
   !   1s, 2s, 2p, 3s, 3p, 3d, 4s, 4p, 4d, 4f, 5s, 5p, 5d, 6s

   ! Associate each possible QN_n,QN_l combination with an index number that
   !   represents the order that this pair is supposed to appear from lowest
   !   energy to highest.  (Applies only to states up to 4f as this is the
   !   highest "core" state considered.)
!   QN_nlOrderIndex(1,1) = 1
!   QN_nlOrderIndex(2,1) = 2
!   QN_nlOrderIndex(2,2) = 3
!   QN_nlOrderIndex(3,1) = 4
!   QN_nlOrderIndex(3,2) = 5
!   QN_nlOrderIndex(4,1) = 6
!   QN_nlOrderIndex(3,3) = 7
!   QN_nlOrderIndex(4,2) = 8
!   QN_nlOrderIndex(5,1) = 9
!   QN_nlOrderIndex(4,3) = 10
!   QN_nlOrderIndex(5,2) = 11
!   QN_nlOrderIndex(6,1) = 12
!   QN_nlOrderIndex(4,4) = 13
   QN_nlOrderIndex(1,1) = 1
   QN_nlOrderIndex(2,1) = 2
   QN_nlOrderIndex(2,2) = 3
   QN_nlOrderIndex(3,1) = 4
   QN_nlOrderIndex(3,2) = 5
   QN_nlOrderIndex(3,3) = 6
   QN_nlOrderIndex(4,1) = 7
   QN_nlOrderIndex(4,2) = 8
   QN_nlOrderIndex(4,3) = 9
   QN_nlOrderIndex(4,4) = 10
   QN_nlOrderIndex(5,1) = 11
   QN_nlOrderIndex(5,2) = 12
   QN_nlOrderIndex(5,3) = 13
   QN_nlOrderIndex(6,1) = 14

   ! Record the number of states associated with each QN_n,QN_l by index number.
!   numOrbitalStates(1)  = 1
!   numOrbitalStates(2)  = 1
!   numOrbitalStates(3)  = 3
!   numOrbitalStates(4)  = 1
!   numOrbitalStates(5)  = 3
!   numOrbitalStates(6)  = 1
!   numOrbitalStates(7)  = 5
!   numOrbitalStates(8)  = 3
!   numOrbitalStates(9)  = 1
!   numOrbitalStates(10) = 5
!   numOrbitalStates(11) = 3
!   numOrbitalStates(12) = 1
!   numOrbitalStates(13) = 7
   numOrbitalStates(1)  = 1
   numOrbitalStates(2)  = 1
   numOrbitalStates(3)  = 3
   numOrbitalStates(4)  = 1
   numOrbitalStates(5)  = 3
   numOrbitalStates(6)  = 5
   numOrbitalStates(7)  = 1
   numOrbitalStates(8)  = 3
   numOrbitalStates(9)  = 5
   numOrbitalStates(10) = 7
   numOrbitalStates(11) = 1
   numOrbitalStates(12) = 3
   numOrbitalStates(13) = 5
   numOrbitalStates(14) = 1
   numOrbitalStates(:) = numOrbitalStates(:) * spin

   ! Determine the band index number of the requested core excitation
   !   (excluding the ground state of course).
   if (excitedQN_n /= 0) then

      ! Initialize the core state index number.
      excitedCoreStateIndex = 0

      ! Accumulate all the states except the core state to be excited.
      do i = 1, QN_nlOrderIndex(excitedQN_n,excitedQN_l+1) - 1
         excitedCoreStateIndex = excitedCoreStateIndex + numOrbitalStates(i)
      enddo

      ! Increment by 1 to have the index number of the core state.
      excitedCoreStateIndex = excitedCoreStateIndex + 1
   endif

   ! Turn on the flag indicating that the core state structures have been
   !   initialized.
   coreStructInit = 1

end subroutine initCoreStateStructures


! Compute the LAT (Linear Analytic Tetrahedron) analog of the
!   electronPopulation array. While the standard electronPopulation
!   uses Gaussian broadening and Fermi filling to assign a BZ-
!   integration weight to each state (band, kpoint, spin), this
!   subroutine uses the Bloechl tetrahedron method (PRB 49, 16223,
!   1994) to determine these weights analytically.
!
! The physical question both arrays answer is: "how occupied is this
!   state, weighted for Brillouin zone integration?" For the standard
!   method, the answer depends on kPointWeight and the Fermi function.
!   For the LAT method, the answer comes from the Bloechl corner
!   integration weights, which decompose each tetrahedron's occupied
!   volume among its four corners based on the geometry of the linear
!   eigenvalue interpolation within the tetrahedron.
!
! The resulting array electronPopulation_LAT(n, k, spin) can replace
!   electronPopulation in the bond order and effective charge
!   computations, providing a broadening-parameter-free alternative
!   for these integrated properties.
!
! PRECONDITION: This subroutine must be called AFTER the call to
!   shiftEnergyEigenValues(occupiedEnergy), so that the Fermi level
!   sits at E = 0 in the shifted eigenvalue spectrum.
!
! Algorithm (see PSEUDOCODE.md sections 3 and 3a for derivation):
!   For each band n and each tetrahedron T:
!   1. Look up eigenvalues at the four corners via fullToIBZMap
!      (unfolding the IBZ to the full mesh).
!   2. Sort the four eigenvalues ascending; track the permutation
!      so we can map corner weights back to the correct k-points.
!   3. Compute corner integration weights at E_Fermi = 0 using the
!      Bloechl formulas. Four cases arise depending on where E_Fermi
!      falls among the sorted eigenvalues (see Case comments below).
!   4. Accumulate each corner's weight into the output array, using
!      the sort permutation to attribute weight to the correct
!      k-point.
subroutine computeElectronPopulation_LAT

   ! Import the necessary modules.
   use O_Kinds
   use O_Potential,       only: spin
   use O_Input,           only: numStates
   use O_KPoints,         only: numKPoints, &
         & numTetrahedra, tetraVol, tetrahedra, &
         & fullToIBZMap
   use O_SecularEquation, only: energyEigenValues

   ! Make sure that no variables are declared accidentally.
   implicit none

   ! -------------------------------------------------
   ! Local variable declarations.
   ! -------------------------------------------------

   ! Loop indices: h = spin orientation, n = band (state),
   !   t = tetrahedron, i and j = corner indices used for
   !   eigenvalue lookup and sorting.
   integer :: h, n, t, i, j

   ! Index of the minimum eigenvalue during the selection sort of
   !   the four corner eigenvalues.
   integer :: minIdx

   ! The IBZ k-point index of a given corner, obtained by following
   !   the sort permutation back to the original corner and then
   !   through the fullToIBZMap.
   integer :: ibzKP

   ! Temporary variables for swapping values during the selection
   !   sort (real for eigenvalues, integer for the permutation).
   real (kind=double) :: tempVal
   integer :: tempInt

   ! Tolerance for detecting degenerate tetrahedra where two or more
   !   corner eigenvalues coincide. When denominators in the Bloechl
   !   formulas fall below this threshold, the tetrahedron's
   !   contribution is treated as negligible and skipped.
   real (kind=double), parameter :: tol = 1.0d-12

   ! After shiftEnergyEigenValues, the Fermi level sits at exactly
   !   zero in the shifted spectrum. All Bloechl corner weights are
   !   evaluated at this energy.
   real (kind=double), parameter :: eFermi = 0.0_double

   ! The four sorted corner eigenvalues in ascending order. Named
   !   e1..e4 matching the notation in PSEUDOCODE.md sections 2
   !   and 3a.
   real (kind=double) :: e1, e2, e3, e4

   ! Occupied fraction f (volume of the occupied sub-region divided
   !   by the tetrahedron volume), unoccupied fraction f_un
   !   (complement, used in Case 3), and denom (denominator for
   !   degenerate guards).
   real (kind=double) :: f, f_un, denom

   ! Case 1 edge-intersection parameters. Each t_j measures how far
   !   along edge 1->j the iso-energy surface E_Fermi cuts, as a
   !   fraction from 0 to 1. The occupied sub-tetrahedron near
   !   corner 1 has volume ratio f = t2 * t3 * t4.
   real (kind=double) :: t2, t3, t4

   ! Case 3 edge-intersection parameters. Each s_j measures how far
   !   along edge 4->j the iso-energy surface cuts. The unoccupied
   !   sub-tetrahedron near corner 4 has volume f_un = s1 * s2 * s3.
   real (kind=double) :: s1, s2, s3

   ! Case 2 (middle range) intersection parameters. The iso-energy
   !   surface cuts four edges of the tetrahedron:
   !     a on edge 1->3,  b on edge 1->4,
   !     c on edge 2->3,  d on edge 2->4.
   !   These are fractions in [0,1] giving the location of the
   !   intersection points A, B, C, D that bound the pentahedral
   !   occupied region.
   real (kind=double) :: a, b, c, d

   ! Eigenvalue differences used as denominators in the Case 2
   !   intersection parameters. Named to match the PSEUDOCODE
   !   notation: e31 = e3 - e1, e41 = e4 - e1, e32 = e3 - e2,
   !   e42 = e4 - e2.
   real (kind=double) :: e31, e41, e32, e42

   ! Case 2 sub-tetrahedra volume ratios. The occupied pentahedron
   !   decomposes into three sub-tetrahedra (see PSEUDOCODE.md
   !   section 3a for the geometric construction and derivation):
   !     T_I   = (corner1, corner2, A, B)
   !     T_II  = (corner2, A, B, D)
   !     T_III = (corner2, A, C, D)
   real (kind=double) :: v_I, v_II, v_III

   ! The four corner eigenvalues in their original (unsorted) order,
   !   before the selection sort. One per tetrahedron corner.
   real (kind=double), dimension(4) :: cornerEigenVals

   ! The IBZ k-point index for each of the four tetrahedron corners,
   !   obtained by mapping the full-mesh tetrahedron corner index
   !   through fullToIBZMap.
   integer, dimension(4) :: ibzCorner

   ! Sort permutation array. After sorting, sortPerm(i) gives the
   !   original corner index (1-4) of the i-th smallest eigenvalue.
   !   This is essential for mapping the sorted corner weights back
   !   to the correct k-points in the accumulation step.
   integer, dimension(4) :: sortPerm

   ! The four Bloechl corner integration weights computed for one
   !   tetrahedron. cornerWeights(i) is the fraction of the
   !   tetrahedron's occupied volume attributed to the i-th sorted
   !   corner. These four values sum to f(E_Fermi), the total
   !   occupied fraction of the tetrahedron for this band.
   real (kind=double), dimension(4) :: cornerWeights

   ! -------------------------------------------------
   ! Allocate and initialize the output array.
   ! -------------------------------------------------
   ! Dimensions: (band, kpoint, spin). Each entry will accumulate
   !   contributions from every tetrahedron that shares that k-point
   !   as a corner.

   if (allocated(electronPopulation_LAT)) &
         & deallocate(electronPopulation_LAT)
   allocate(electronPopulation_LAT( &
         & numStates, numKPoints, spin))
   electronPopulation_LAT(:,:,:) = 0.0_double

   write (20, *) "Computing LAT electron population."
   call flush(20)

   ! -------------------------------------------------
   ! Main loop: spin x bands x tetrahedra.
   ! -------------------------------------------------
   ! For each (band, tetrahedron) pair, we compute the Bloechl
   !   corner integration weights at the Fermi energy. These weights
   !   determine how much of the tetrahedron's occupied volume should
   !   be attributed to each of the four corner k-points. The weights
   !   are accumulated into the output array, building up the total
   !   occupation for each (band, kpoint, spin) triple.

   do h = 1, spin

      do n = 1, numStates

         do t = 1, numTetrahedra

            ! ----------------------------------------
            ! Step 1: Look up corner eigenvalues.
            ! ----------------------------------------
            ! The tetrahedra array stores full-mesh k-point indices
            !   (1..numFullMeshKP). We map each corner to its IBZ
            !   representative via fullToIBZMap, then look up the
            !   eigenvalue for this band and spin. The sort
            !   permutation is initialized to the identity.
            do i = 1, 4
               ibzCorner(i) = fullToIBZMap( &
                     & tetrahedra(i, t))
               cornerEigenVals(i) = &
                     & energyEigenValues( &
                     & n, ibzCorner(i), h)
               sortPerm(i) = i
            enddo

            ! ----------------------------------------
            ! Step 2: Sort eigenvalues ascending.
            ! ----------------------------------------
            ! We need the four eigenvalues in sorted order for the
            !   Bloechl formulas. The sort permutation tracks which
            !   original corner each sorted position came from, so
            !   that the corner weights can be mapped back to the
            !   correct k-points. Selection sort is used (optimal
            !   for exactly 4 elements).
            do i = 1, 3
               minIdx = i
               do j = i + 1, 4
                  if (cornerEigenVals(j) < &
                        & cornerEigenVals(minIdx)) &
                        & then
                     minIdx = j
                  endif
               enddo
               if (minIdx /= i) then
                  ! Swap eigenvalues.
                  tempVal = cornerEigenVals(i)
                  cornerEigenVals(i) = &
                        & cornerEigenVals(minIdx)
                  cornerEigenVals(minIdx) = tempVal
                  ! Swap permutation indices.
                  tempInt = sortPerm(i)
                  sortPerm(i) = sortPerm(minIdx)
                  sortPerm(minIdx) = tempInt
               endif
            enddo

            e1 = cornerEigenVals(1)
            e2 = cornerEigenVals(2)
            e3 = cornerEigenVals(3)
            e4 = cornerEigenVals(4)

            ! ----------------------------------------
            ! Step 3: Compute Bloechl corner weights
            !   at E_Fermi = 0.
            ! ----------------------------------------
            ! The Bloechl formulas decompose the tetrahedron's
            !   occupied volume into four corner contributions.
            !   Four cases arise depending on where E_Fermi falls
            !   relative to the sorted corner eigenvalues
            !   e1 <= e2 <= e3 <= e4. See PSEUDOCODE.md section 3a
            !   for the full derivation from the vertex-averaging
            !   property of linear functions over tetrahedra.

            ! Case 0a: Fermi level below all corners. No states in
            !   this tetrahedron are occupied for this band -- skip.
            if (eFermi < e1) cycle

            ! Case 0b: Fermi level above all corners. The entire
            !   tetrahedron is occupied for this band. By the vertex
            !   averaging property, the integral of each barycentric
            !   coordinate over the full tetrahedron is V_T/4, so
            !   each corner gets an equal share of 1/4.
            if (eFermi >= e4) then
               cornerWeights(:) = 0.25_double

            ! Case 1: e1 <= E_Fermi < e2. Only corner 1 lies below
            !   the Fermi level. The occupied region is a small sub-
            !   tetrahedron with apex at corner 1. The iso-energy
            !   surface epsilon = E_F cuts the three edges from
            !   corner 1 at fractional distances t2, t3, t4. The
            !   sub-tet volume ratio is f = t2*t3*t4, and the corner
            !   weights follow from vertex averaging over the four
            !   vertices of the sub-tetrahedron.
            elseif (eFermi < e2) then
               denom = (e2-e1) * (e3-e1) * (e4-e1)
               if (abs(denom) < tol) cycle
               t2 = (eFermi - e1) / (e2 - e1)
               t3 = (eFermi - e1) / (e3 - e1)
               t4 = (eFermi - e1) / (e4 - e1)
               f = t2 * t3 * t4
               cornerWeights(2) = &
                     & f * t2 / 4.0_double
               cornerWeights(3) = &
                     & f * t3 / 4.0_double
               cornerWeights(4) = &
                     & f * t4 / 4.0_double
               cornerWeights(1) = f &
                     & - cornerWeights(2) &
                     & - cornerWeights(3) &
                     & - cornerWeights(4)

            ! Case 2: e2 <= E_Fermi < e3 (middle range). Corners 1
            !   and 2 lie below the Fermi level; corners 3 and 4 lie
            !   above. The iso-energy surface cuts four edges,
            !   creating a pentahedral occupied region with vertices
            !   at corners 1 and 2 plus the four intersection points
            !   A, B, C, D. This pentahedron decomposes into three
            !   sub-tetrahedra (T_I, T_II, T_III) with volume ratios
            !   v_I, v_II, v_III. The corner weights come from
            !   applying vertex averaging to each sub-tet and summing.
            elseif (eFermi < e3) then
               e31 = e3 - e1
               e41 = e4 - e1
               e32 = e3 - e2
               e42 = e4 - e2
               if (e31 * e41 < tol .or. &
                     & e32 * e42 < tol) cycle

               ! Intersection parameters: fractional positions where
               !   the iso-energy surface cuts each edge.
               a = (eFermi - e1) / e31
               b = (eFermi - e1) / e41
               c = (eFermi - e2) / e32
               d = (eFermi - e2) / e42

               ! Sub-tetrahedra volume ratios (as fractions of the
               !   full tetrahedron volume).
               v_I   = a * b
               v_II  = a * d * (1.0_double - b)
               v_III = (1.0_double - a) * c * d

               ! Corner weights from vertex averaging over the three
               !   sub-tetrahedra. For each sub-tet, we sum the
               !   barycentric coordinate lambda_j at its four
               !   vertices, multiply by v_k/4, and add contributions
               !   from all three sub-tetrahedra.
               cornerWeights(1) = ( &
                     & v_I * (3.0_double - a - b) &
                     & + v_II*(2.0_double - a - b) &
                     & + v_III * (1.0_double - a) &
                     & ) / 4.0_double
               cornerWeights(2) = ( &
                     & v_I &
                     & + v_II*(2.0_double - d) &
                     & + v_III*(3.0_double-c-d) &
                     & ) / 4.0_double
               cornerWeights(3) = ( &
                     & v_I * a + v_II * a &
                     & + v_III * (a + c) &
                     & ) / 4.0_double
               cornerWeights(4) = ( &
                     & v_I * b &
                     & + v_II * (b + d) &
                     & + v_III * d &
                     & ) / 4.0_double

            ! Case 3: e3 <= E_Fermi < e4. Only corner 4 lies above
            !   the Fermi level. The unoccupied region is a small
            !   sub-tetrahedron near corner 4 (complement of Case 1).
            !   The occupied weights are the whole-tetrahedron weights
            !   (1/4 each) minus the unoccupied sub-tet contributions.
            else
               denom = (e4-e1) * (e4-e2) * (e4-e3)
               if (abs(denom) < tol) then
                  cornerWeights(:) = 0.25_double
               else
                  s1 = (e4-eFermi) / (e4-e1)
                  s2 = (e4-eFermi) / (e4-e2)
                  s3 = (e4-eFermi) / (e4-e3)
                  f_un = s1 * s2 * s3
                  cornerWeights(1) = 0.25_double &
                        & - f_un*s1 / 4.0_double
                  cornerWeights(2) = 0.25_double &
                        & - f_un*s2 / 4.0_double
                  cornerWeights(3) = 0.25_double &
                        & - f_un*s3 / 4.0_double
                  cornerWeights(4) = &
                        & (1.0_double - f_un) &
                        & - cornerWeights(1) &
                        & - cornerWeights(2) &
                        & - cornerWeights(3)
               endif
            endif

            ! ----------------------------------------
            ! Step 4: Accumulate into the output.
            ! ----------------------------------------
            ! The sort permutation maps each sorted position i back
            !   to its original corner index sortPerm(i). The IBZ
            !   k-point for that corner is ibzCorner(sortPerm(i)).
            !   Multiple tetrahedra sharing a k-point accumulate
            !   their contributions additively, building the total
            !   occupation weight for each (band, kpoint, spin).
            do i = 1, 4
               ibzKP = ibzCorner(sortPerm(i))
               electronPopulation_LAT(n,ibzKP,h) &
                     & = electronPopulation_LAT( &
                     & n, ibzKP, h) &
                     & + cornerWeights(i) * tetraVol
            enddo

         enddo ! t (tetrahedra)

         ! Progress indicator for long computations.
         if (mod(n, 50) == 0) then
            write (20, ADVANCE="NO", &
                  & FMT="(a1)") "."
            call flush(20)
         endif

      enddo ! n (states/bands)

   enddo ! h (spin)

   write (20, *) ""

   ! Diagnostic: print total occupation per spin. For a fully
   !   occupied band, the sum over all k-points equals 1.0. The
   !   total over all occupied bands and spins should be consistent
   !   with the electron count from populateStates.
   do h = 1, spin
      write (20, fmt="(a,i1,a,f12.6)") &
            & "  LAT electron pop spin ", h, &
            & " total: ", &
            & sum(electronPopulation_LAT(:,:,h))
   enddo
   call flush(20)

end subroutine computeElectronPopulation_LAT


subroutine cleanUpPopulation

   implicit none

   deallocate (electronPopulation)
   deallocate (sortedEnergyEigenValues)
   deallocate (indexEnergyEigenValues)

   ! The LAT population array is only allocated when the LAT
   !   integration method was active (kPointIntgCode == 1). Guard
   !   the deallocation with an allocated() check.
   if (allocated(electronPopulation_LAT)) then
      deallocate(electronPopulation_LAT)
   endif

end subroutine cleanUpPopulation

end module O_Populate
