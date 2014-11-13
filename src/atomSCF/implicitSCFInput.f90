module implicitSCFInputSubs

private
public :: implicitSCFInput

contains

subroutine allocateMemory

   ! Import necessary definitions.
   use O_Kinds

   ! Import necessary data modules.
   use EnergyData
   use AtomData

   implicit none

   ! Define local variables.

   ! Allocate space for various arrays.
   allocate (eigenValues(numOrbitals))
   allocate (orbitalKE(numOrbitals))
   allocate (orbitalPE(numOrbitals))

end subroutine allocateMemory



subroutine coreOrbInit

   ! Import necessary definitions.
   use O_Kinds

   ! Import necessary data modules.
   use AtomData
   use ExecutionData

   implicit none

   ! Define local variables
   integer :: i,j ! Loop index variables.
   integer :: currentOrbital  ! Counter for the number of orbitals.
   integer, dimension (15) :: corePrincipleQN  ! Order of n Quantum number.
   integer, dimension (15) :: coreAngMomQN     ! Order of l Quantum number.
   real (kind=double) :: currentCoreCharge ! Accumulator for the core charge
         ! in the system.
   real (kind=double) :: orbitalSpinAngMom

   ! Initialize an ordered list of the core orbital's principle quantum
   !   number (n).  The ordering is based on energy of the core orbitals of
   !   elements from the periodic table.  The valence orbitals will have to
   !   be specified specifically because in some cases (e.g. K and Ca) the
   !   higher energy orbitals occur out of this order.
   corePrincipleQN(1)  = 1
   corePrincipleQN(2)  = 2
   corePrincipleQN(3)  = 2
   corePrincipleQN(4)  = 3
   corePrincipleQN(5)  = 3
   corePrincipleQN(6)  = 3
   corePrincipleQN(7)  = 4
   corePrincipleQN(8)  = 4
   corePrincipleQN(9)  = 4
   corePrincipleQN(10) = 5
   corePrincipleQN(11) = 5
   corePrincipleQN(12) = 4
   corePrincipleQN(13) = 5
   corePrincipleQN(14) = 6
   corePrincipleQN(15) = 6

   ! Initialize an ordered list of the core orbital's angular momentum
   !   quantum number (l).  This list corresponds to the list above and
   !   so applies only to the core orbitals.
   coreAngMomQN(1)  = 0
   coreAngMomQN(2)  = 0
   coreAngMomQN(3)  = 1
   coreAngMomQN(4)  = 0
   coreAngMomQN(5)  = 1
   coreAngMomQN(6)  = 2
   coreAngMomQN(7)  = 0
   coreAngMomQN(8)  = 1
   coreAngMomQN(9)  = 2
   coreAngMomQN(10) = 0
   coreAngMomQN(11) = 1
   coreAngMomQN(12) = 3
   coreAngMomQN(13) = 2
   coreAngMomQN(14) = 0
   coreAngMomQN(15) = 1

   ! Determine the quantum numbers for each core orbital in the sys and the
   !   amount of charge that occupies each orbital.  For systems that are
   !   not spin-polarized, the number of core orbitals given in the
   !   input file is the number of core orbitals in the system and each
   !   orbital will be doubly occupied.  The total spin of each orbital
   !   will be zero.  For systems that are spin-polarized, the number of
   !   core orbitals will be double that given in the input file and each
   !   orbital will be singly occupied.  The spin of each orbital will be
   !   either -0.5 or +0.5.  If the calculation is relativistic then it is
   !   forced to be spin-polarized and the occupation is modified so that
   !   one spin orbital (up) is over occupied while the other is under
   !   occupied.

   ! Initialize a counter for the number of core orbitals.
   currentOrbital = 0

   ! Initialize an accumulator for the charge in the system.
   currentCoreCharge = 0.0_double

   ! We can abort this subroutine now if there are no core orbitals (e.g.
   !   we are computing for the H atom).
   if (numCoreOrb == 0) then
      coreCharge = currentCoreCharge
      return
   endif

   ! Initialize the factor that is used to assign electron spin orientation
   !   for a given orbital.  In the non spin-polarized case this will
   !   remain at 0.  In the spin-polarized case it will oscillate between
   !   -0.5 and 0.5 to indicate spin down and spin up.
   if (doSpinPol == 1) then
      orbitalSpinAngMom = -0.5_double
   else
      orbitalSpinAngMom =  0.0_double
   endif

   do i = 1, numCoreOrb
      do j = 0, doSpinPol ! This will only do one iteration for non-spinpol.

         ! Increment the total number of orbitals.
         currentOrbital = currentOrbital + 1

         ! Assign the principle, angular, and spin quantum numbers to this
         !   orbital from the predefined data.
         orbitalQNn(currentOrbital) = corePrincipleQN(i)
         orbitalQNl(currentOrbital) = coreAngMomQN(i)
         orbitalQNs(currentOrbital) = orbitalSpinAngMom

         ! Assign the charge (occupation number) to this orbital.  This
         !   statement assumes that we are doing a spin-polarized
         !   calculation.  (e.g. 1s has an up and a down orbital.)
         orbitalCharge(currentOrbital) = 2 * coreAngMomQN(i) + 1

         ! In the case that the calculation is not spin polarized we need to
         !   adjust the occupation by doubling it.  (e.g. 1s has only one
         !   state for both up and down spin so it must be doubly occupied.)
         if (doSpinPol == 0) then
            orbitalCharge(currentOrbital) = &
                  & 2 * orbitalCharge(currentOrbital)
         endif

         ! In the case that the calculation is relativistic we need to
         !   adjust the occupation by shifting it to over occupy the spin
         !   up states.  (I'm not sure what this is all about yet.)
         if (doRelativistic == 1) then
            orbitalCharge(currentOrbital) = 2 * &
                  & (coreAngMomQN(i) + orbitalSpinAngMom) + 1
         endif

         ! Accumulate the charge assigned to an orbital in this iteration.
         currentCoreCharge = currentCoreCharge + &
               & orbitalCharge(currentOrbital)

         ! In the case where no charge was assigned to an orbital (which
         !   could happen for the orbitalQNs = -0.5 and orbitalQNl = 0
         !   case in the relativistic calculation) we do not consider
         !   that orbital to exist and we must reduce the orbital counter
         !   by one.
         if (abs(orbitalCharge(currentOrbital)) <= 0.0_double) then
            currentOrbital = currentOrbital - 1
         endif

         ! In the case that the calculation is spin-polarized we need to
         !   make the spin angular momentum QN oscillate for the next j
         !   loop iteration.
         if (doSpinPol == 1) then
            orbitalSpinAngMom = -orbitalSpinAngMom
         endif
      enddo
   enddo

   ! Record the new value for the number of core orbitals.
   numCoreOrb = currentOrbital

   ! Record the total accumulated core charge.
   coreCharge = currentCoreCharge

end subroutine coreOrbInit



subroutine valeOrbInit

   ! Import necessary definitions.
     use O_Kinds

   ! Import necessary data modules.
   use AtomData
   use ExecutionData
   use TempOrbitalData

   implicit none

   ! Define local variables
   integer :: i,j ! Loop index variables.
   integer :: currentOrbital  ! Counter for the number of orbitals.
   real (kind=double) :: currentValeCharge ! Accumulator for the valence
         ! charge in the system.
   real (kind=double) :: orbitalSpinAngMom

   ! Determine the quantum numbers for each valence orbital and the
   !   amount of charge that occupies each orbital.  For systems that are
   !   not spin-polarized, the number of valence orbitals given in the
   !   input file is the number of valence orbitals in the system and each
   !   orbital will be occupied with the combined spin-up and spin-down
   !   charge.  The total spin of each orbital will be the sum of the two,
   !   (often zero).  For systems that are spin-polarized, the number of
   !   valence orbitals will be double that given in the input file and each
   !   orbital will be occupied with either the spin-up charge or the
   !   spin-down charge depending on the current iteration.  The spin of
   !   each orbital will be  either -0.5 or +0.5.  If the calculation is
   !   relativistic then it is forced to be spin-polarized and the
   !   occupation is modified so that one spin orbital (up) is over occupied
   !   while the other (down) is under occupied (for some odd reason).

   ! Initialize a counter for the number of orbitals, this is a continuation
   !   of the core orbital assignments that just completed.
   currentOrbital = numCoreOrb

   ! Initialize an accumulator for the valence charge in the system.
   currentValeCharge = 0.0_double

   ! We can abort this subroutine now if there are no valence orbitals (I
   !   have no idea why someone would do this).
   if (numValeOrb == 0) then
      valeCharge = currentValeCharge
      return
   endif

   ! Initialize the factor that is used to assign electron spin orientation
   !   for a given orbital.  In the non spin-polarized case this will
   !   remain at 0.  In the spin-polarized case it will oscillate between
   !   -0.5 and 0.5 to indicate spin down and spin up.
   if (doSpinPol == 1) then
      orbitalSpinAngMom = -0.5_double
   else
      orbitalSpinAngMom =  0.0_double
   endif

   do i = 1, numValeOrb
      do j = 0, doSpinPol ! This will only do one iteration for non-spinpol.

         ! Increment the total number of orbitals.
         currentOrbital = currentOrbital + 1

         ! Assign the principle, angular, and spin quantum numbers to this
         !   orbital from the data read in from the input file and stored
         !   in the temporary arrays.
         orbitalQNn(currentOrbital) = tempValeQNn(i)
         orbitalQNl(currentOrbital) = tempValeQNl(i)
         orbitalQNs(currentOrbital) = orbitalSpinAngMom

         ! Assign the charge (occupation number) to this orbital.  This
         !   statement assumes that we are doing a non spin-polarized
         !   calculation.  (This is the opposite assumption from the
         !   core version of this subroutine.)
         orbitalCharge(currentOrbital) = tempValeOrbChargeDn(i) + &
               & tempValeOrbChargeUp(i)

         ! In the case that the calculation is spin polarized we need to
         !   adjust the occupation by assigning the appropriate up or down
         !   spin.
         if (doSpinPol == 1) then
            if     (orbitalSpinAngMom == -0.5_double) then
               orbitalCharge(currentOrbital) = tempValeOrbChargeDn(i)
            elseif (orbitalSpinAngMom ==  0.5_double) then
               orbitalCharge(currentOrbital) = tempValeOrbChargeUp(i)
            endif

            ! In the case that the calculation is also relativistic, then
            !   we will adjust the occupation by shifting it to over occupy
            !   the spin up states.  (Don't know what this is about yet.)
            if (doRelativistic == 1) then
               orbitalCharge(currentOrbital) = &
                     & orbitalCharge(currentOrbital) * (2.0_double * &
                     & (tempValeQNl(i) + orbitalSpinAngMom) + 1.0_double) / &
                     & (4.0_double * tempValeQNl(i) + 2.0_double)
            endif
         endif

         ! Accumulate the charge assigned to an orbital in this iteration.
         currentValeCharge = currentValeCharge + &
               & orbitalCharge(currentOrbital)

         ! In the case where l+s=j is less than 0 for relativistic calculations
         !   we need to remove this orbital.  (Not sure why.)
         if (((orbitalQNl(currentOrbital)+orbitalQNs(currentOrbital)) <= &
               & 0.0_double) .and. (doRelativistic == 1)) then
            currentOrbital = currentOrbital - 1
         endif

         ! In the case that the calculation is spin-polarized we need to
         !   make the spin angular momentum QN oscillate for the next j
         !   loop iteration.
         if (doSpinPol == 1) then
            orbitalSpinAngMom = -orbitalSpinAngMom
         endif
      enddo
   enddo

   ! Record the new value for the number of valence orbitals.
   numValeOrb = currentOrbital - numCoreOrb

   ! Record the total accumulated vale charge.
   valeCharge = currentValeCharge

   ! Record the total number of orbitals.
   numOrbitals = currentOrbital

   ! At this point we can deallocate the arrays that were used to form the
   !   valence orbitals.
   deallocate (tempValeQNn)
   deallocate (tempValeQNl)
   deallocate (tempValeOrbChargeDn)
   deallocate (tempValeOrbChargeUp)

end subroutine valeOrbInit



subroutine finalizeOrbitals

   ! Import necessary definitions.
   use O_Kinds
   use O_Constants

   ! Import necessary data modules.
   use AtomData

   implicit none

   ! Define local variables
   integer :: i,j ! Loop index variables.

   ! Determine if the core orbitals and the valence orbitals are the same.
   !   If they are, then we abort the computation.  This is simply a check
   !   to make sure that no two states have all the same quantum numbers.
   do i = 1, numOrbitals-1
      do j = i+1, numOrbitals

         ! Compare the two principle quantum numbers.
         if (orbitalQNn(i) /= orbitalQNn(j)) then
            cycle
         endif

         ! Compare the orbital angular momentum quantum numbers.
         if (orbitalQNl(i) /= orbitalQNl(j)) then
            cycle
         endif

         ! Compare the spin angular momentum quantum numbers.
         if (abs(orbitalQNs(i) - orbitalQNs(j)) > smallThresh) then
            cycle
         endif

         ! If we made it this far without cycling to the next iteration,
         !   then all the quantum numbers for two different orbitals are
         !   the same and we should abort.
         write (20,*) "All quantum numbers for two orbitals are the same.",&
               & i, j
         stop

      enddo
   enddo

   ! Compute the total electron charge as the sum of the core and valence
   !   charges.
   electronCharge = valeCharge + coreCharge

   ! Compute the ionicity (ioninc charge) as the difference between the
   !   electron charge and the nuclear charge (Given by the atomic number).
   ionicCharge = atomicNumber - electronCharge

end subroutine finalizeOrbitals



subroutine implicitSCFInput

   ! Import necessary definitions.
   use O_Kinds

   ! Import necessary data modules.
   use AtomData
   use ExecutionData

   ! Import necessary external modules.
   use O_RadialGrid

   implicit none

   ! Demand spin polarization for the case of a relativistic calculation.
   if (doRelativistic == 1) then
      doSpinPol = 1
   endif

   ! Set up the radial grid.
   call setupRadialGrid(atomicNumber)

   ! Compute the occupation numbers and orbital energies for the core.
   call coreOrbInit

   ! Compute the occupation numbers and orbital energies for the valence.
   call valeOrbInit

   ! Make final analysis and assignments based on all the orbital info.
   call finalizeOrbitals

   ! Allocate space for data necessary throughout the program.
   call allocateMemory

   ! Initialize the current iteration, converged flag, initial choice of
   !   eigenValue solver, and the value for the previous maximum deviation of
   !   the new potential from the old.
   iteration      = 0 ! Iteration starts at 0.
   converged      = 0 ! System starts unconverged.
   solverChoice   = 1 ! This solver does not need an initial guess.
   oldMaxDeltaPot = 1.0_double ! Assume the last deviation was huge.

end subroutine implicitSCFInput


end module implicitSCFInputSubs
