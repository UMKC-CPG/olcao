module O_Basis

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Index variables
   integer, allocatable, dimension (:,:) :: lmAlphaIndex ! The first
         !   dimension references the maximum number of atomic alphas in
         !   the system.  The second dimension references the number of
         !   atomic types in the system.  Each array index for a given
         !   atomic type represents the total number of l,m combinations
         !   that this alpha will be used for.  E.g. if the alpha is only
         !   used for s orbitals, then the total number of l,m combos is
         !   1 (l=0,m=0).  If the alpha is used for s and p orbitals then
         !   the total number of l,m combos is l=0,m=0, l=1, m=-1,0,1 for
         !   1+3=4.  s,p,d leads to 9 ...  The idea is for each alpha of
         !   each atomic type to be assigned a number indicating how many
         !   different (spin degenerate) states it is used for.  *NOT* the
         !   total number of states for each alpha, but rather the number
         !   of different kinds of states.
   integer, allocatable, dimension (:,:) :: lmIndex ! The first dimension
         !   references the maximum number of states for an atomic type out
         !   of all available atomic types, the second dimension references
         !   the number of atomic types in the system.  Each array index for
         !   a given atomic type is an integer that uniquely represents a
         !   particular l,m combination.  1 = l=0,m=0; 2 = l=1,m=-1; 3 = 
         !   l=1,m=0; 4 = l=1,m=1; 5 = l=2,m=-2 ... 9 = l=2,m=2; 10 =
         !   l=3,m=-3; 11 = l=3,m=-2; ... 15 = l=3,m=2; 16 = l=3,m=3.  In this
         !   way each state of each atomic type can be identified as to what
         !   its QN_l and QN_m are as part of a linear ordered list.


   ! The system basis functions.
   real (kind=double), allocatable, dimension (:,:,:) :: basis ! The last
         !   dimension is over the number of types in the system.  The 
         !   middle dimension is over the total number of states for the
         !   given type.  The first dimension is over the number of alphas
         !   necessary for the state indexed in the second dimension.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine renormalizeBasis

   ! Include the modules we need
   use O_Kinds
   use O_AtomicTypes, only: numAtomTypes, maxNumStates, maxNumAtomAlphas, &
         & atomTypes

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the local variables for this subroutine
   integer :: i,j,k ! Loop index variables.
   integer :: totalmCount ! Counter of the number of time that the m quantum
         ! number can be changed for a given atomic type.
   integer :: currentNumAlphas ! The current number of atomic alphas for this
         ! l quantum number and this type.
   integer :: currentOrbType ! The l quantum number of the current wave
         ! function for either valence or core (depending the loop) orbitals.
   real (kind=double), dimension (16) :: angNorm


   ! Allocate space to hold the necessary data to make a matrix out of the
   !   wave functions and properly index them.
   allocate (lmIndex      (maxNumStates,numAtomTypes))
   allocate (lmAlphaIndex (maxNumAtomAlphas,numAtomTypes))
   allocate (basis        (maxNumAtomAlphas,maxNumStates,numAtomTypes))

   ! Initialize these matrices to zero.
   lmIndex(:,:)      = 0
   lmAlphaIndex(:,:) = 0
   basis(:,:,:)      = 0.0_double

   ! Initialize the values in angNorm.  These are angular normalizations.  The
   !   coefficient and the real spherical harmonic function in cartesian
   !   coordinates follows.  Computed according to Romanowski, et al., 
   !   ACTA PHYSICA POLONICA B, 39, p1985, (2008).  (And references therein.)
   angNorm(1)  = 0.2820947917738782_double ! sqrt(1/4pi)    * 1
   angNorm(2)  = 0.4886025119029199_double ! sqrt(3/4pi)    * x/r
   angNorm(3)  = 0.4886025119029199_double ! sqrt(3/4pi)    * y/r
   angNorm(4)  = 0.4886025119029199_double ! sqrt(3/4pi)    * z/r
   angNorm(5)  = 1.0925484305920790_double ! sqrt(15/4pi)   * xy/r^2
   angNorm(6)  = 1.0925484305920790_double ! sqrt(15/4pi)   * xz/r^2
   angNorm(7)  = 1.0925484305920790_double ! sqrt(15/4pi)   * yz/r^2
   angNorm(8)  = 0.5462742152960396_double ! sqrt(15/16pi)  * (x^2-y^2)/r^2
   angNorm(9)  = 0.3153915652525200_double ! sqrt(5/16pi)   * (3z^2-r^2)/r^2
   angNorm(10) = 2.8906114426405543_double ! sqrt(105/4pi)  * xyz/r^3
   angNorm(11) = 1.4453057213202771_double ! sqrt(105/16pi) * z(x^2-y^2)/r^3
   angNorm(12) = 0.5900435899266435_double ! sqrt(35/32pi)  * x(x^2-3y^2)/r^3
   angNorm(13) = 0.5900435899266435_double ! sqrt(35/32pi)  * y(y^2-3x^2)/r^3
   angNorm(14) = 0.3731763325901154_double ! sqrt(7/16pi)   * z(5z^2-3r^2)/r^3
   angNorm(15) = 0.4570457994644657_double ! sqrt(21/32pi)  * x(5z^2-r^2)/r^3
   angNorm(16) = 0.4570457994644657_double ! sqrt(21/32pi)  * y(5z^2-r^2)/r^3


   ! This will store the index number for the maximum index that each
   !   atomic alpha will be used for.  So consider the first atomic alpha
   !   for a given type.  All orbital kinds (0, 1, 2 l quantum numbers)
   !   will use this alpha if we assume that the alphas requested for each
   !   orbital kind are 22 for l=0, 20 for l=1, and 16 for l=2.  Here,
   !   the maximum index will be 9 since the l=0 case will use 1 array
   !   slot, the l=1 case will use 3 additional slots, and finally, the
   !   l=2 case will use 5 slots for a total of 9.  Now, consider the
   !   case of the 17th atomic alpha for a given type using the same
   !   (22,20,16) number of alpha requests (l=0,1,2 respectively).  Here
   !   the max index will be 4 since l=0 will be at index 1, and l=2
   !   will use indices 2, 3, and 4 for m = -1,0,1.  So, this stores the
   !   maximum number of possible quantum number states for a given alpha.

   ! Loop over each type.
   do i = 1, numAtomTypes

      ! Loop over the maximum number of alphas for this atomic type.
      do j = 1, atomTypes(i)%numOrbAlphas(1)

         ! The number of alphas for each type of orbital (l=0,1,2,3) are
         !   compared to the current j index.   If the atomic alpha indexed
         !   by j will be used by all four orbital types, then the number of
         !   wave functions that will use that alpha will be sixteen.  If the
         !   atomic alpha indexed by j will be used only by the l=0, and l=1
         !   orbital types, then only four wave functions will use that
         !   alpha.  Etc.  The number comes from the m=0,+-1..+-l quantum
         !   number.  An l=0 and an l=1 pair of orbitals will contribute
         !   1 + 3 to get 4.
         if (j <= atomTypes(i)%numOrbAlphas(4)) then
            lmAlphaIndex(j,i) = 16
         elseif (j <= atomTypes(i)%numOrbAlphas(3)) then
            lmAlphaIndex(j,i) = 9
         elseif (j <= atomTypes(i)%numOrbAlphas(2)) then
            lmAlphaIndex(j,i) = 4
         else
            lmAlphaIndex(j,i) = 1
         endif
      enddo
   enddo

   ! Store all the wave functions (core and valence) into one large three
   !   dimensional matrix.  The first dimension is basically all the
   !   wave functions for a given atomic type put into one long array.
   !   The second dimension indexes the wave functions of the current type
   !   starting with the valence wave functions, and ending with the core
   !   wave functions.  The trick of course is that each wave function is
   !   repeated 2*l+1 times (once for each available m quantum number for
   !   the orbital's l quantum number).  In each case, the wavefunction is
   !   renormalized according to a specific factor.  The last dimension is
   !   the index of the atomic type under consideration.


   ! Loop over each atomic type.
   do i = 1, numAtomTypes

      ! Initialize the counter of all the m quantum numbers possible for this
      !   atomic type over all the orbitals for this atomic type.  This is
      !   essentially the number of states for the atomic type.
      totalmCount = 0

      ! Loop over all the valence radial functions for this atomic type.
      do j = 1, atomTypes(i)%numValeRadialFns

         ! Get the current orbital type (l quantum number)
         currentOrbType = atomTypes(i)%valeQN_lList(j)

         ! Get the current number of alphas needed for the current l
         !   quantum number of the present j indexed valence radial function.
         currentNumAlphas = atomTypes(i)%numOrbAlphas(currentOrbType + 1)

         ! Loop over the number of possible m quantum numbers for the
         !   current l quantum number orbital type.
         do k = 1, 2 * currentOrbType + 1

            ! Increment the counter of the total number of m quantum numbers
            !   used for this atomic type.
            totalmCount = totalmCount + 1

            ! Save the index number that references which renormalization
            !   factor to use for the current l,m combination.  Essentially,
            !   each state of the current atomic type has an associated number
            !   that denotes what the l,m combination is.  For example: an l=0,
            !   and m=0 for the 1s orbital would make lmIndex(1,i) = 0+1.
            !   Then, the l=0 and m=0 for 2s orbital would make lmIndex(2,i) = 
            !   0+1.  Then, for l=1, m=-1 for the 2p orbital, lmIndex(3,i) = 
            !   1+1=2.  Then, l=1,m=0 for 2p again makes lmIndex(4,i) = 1+2=3.
            !   Finally, l=1,m=1 to complete 2p makes lmIndex(5,i) = 1+3=4.
            !   lmIndex(6,i) = 1 for 3s, lmIndex(7,i) = 2, lmIndex(8,i) = 3,
            !   lmIndex(9,i) = 4 for 3p (m=-1,0,1).  For 3d we have
            !   lmIndex(10,i) = 4+1=5, 4+2, 4+3, 4+4, 4+5=9.
            lmIndex(totalmCount,i) = currentOrbType**2 + k

            ! Record the radial function of this atomic type multiplied by the
            !   l,m determined renormalization factor.  This is now a system
            !   basis function.
            basis(:currentNumAlphas,totalmCount,i) = &
               & atomTypes(i)%valeRadialFns(:currentNumAlphas,j) * &
               & angNorm(lmIndex(totalmCount,i))
         enddo
      enddo

      ! Now, loop over all the core radial functions for this atomic type.
      do j = 1, atomTypes(i)%numCoreRadialFns

         ! Get the current orbital type (l quantum number)
         currentOrbType = atomTypes(i)%coreQN_lList(j)

         ! Record the current number of alphas needed for the current l
         !   quantum number of the present j indexed core wave function.
         currentNumAlphas = atomTypes(i)%numOrbAlphas(currentOrbType + 1)

         ! Loop over the number of possible m quantum numbers for the
         !   current l quantum number orbital type.
         do k = 1, 2 * currentOrbType + 1

            ! Increment the counter of the total number of m quantum numbers
            !   used for this atomic type.  This essentially counts the total
            !   number of states.
            totalmCount = totalmCount + 1

            ! Save the index number that references which renormalization
            !   factor to use for the current l,m combination.
            lmIndex(totalmCount,i) = currentOrbType**2 + k

            ! Record the radial function of this atomic type multiplied by the
            !   l,m determined renormalization factor.  This is now a system
            !   basis function.
            basis(:currentNumAlphas,totalmCount,i) = &
               & atomTypes(i)%coreRadialFns(:currentNumAlphas,j) * &
               & angNorm(lmIndex(totalmCount,i))
         enddo
      enddo
   enddo
end subroutine renormalizeBasis


subroutine cleanUpBasis

   implicit none

   deallocate (lmIndex)
   deallocate (lmAlphaIndex)
   deallocate (basis)

end subroutine cleanUpBasis


subroutine initializeAtomSite (atomNumber,atomIndex,atomType,currentElements,&
      & numTotalStates,numCoreStates,numValeStates,numAlphas,lmIndexCurrent,&
      & lmAlphaIndexCurrent,cartPos,alphas,currentBasisFns)

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3
   use O_AtomicSites, only: numAtomSites, atomSites
   use O_AtomicTypes, only: maxNumAtomAlphas, maxNumStates, numAtomTypes, &
         & atomTypes

   ! Input variables that request which atom to extract info from, and where
   !   to put it.
   integer :: atomNumber ! Atom number
   integer :: atomIndex  ! Array index to store.

   ! Output varibales the will be filled in atomIndex array slot with data
   !   from atom number atomNumber.
   integer, dimension (2)    :: atomType
   integer, dimension (2)    :: currentElements
   integer, dimension (2)    :: numTotalStates
   integer, dimension (2)    :: numCoreStates
   integer, dimension (2)    :: numValeStates
   integer, dimension (2)    :: numAlphas
   integer, dimension (:,:)  :: lmIndexCurrent
   integer, dimension (:,:)  :: lmAlphaIndexCurrent
   real (kind=double), dimension (dim3,2) :: cartPos
   real (kind=double), dimension (:,:)    :: alphas
   real (kind=double), dimension (:,:,:)  :: currentBasisFns


   ! Fill the requested information for the atomNumber provided into the
   !   array index requested by atomIndex.

   ! Get the type number of the requested atom.
   atomType(atomIndex) = atomSites(atomNumber)%atomTypeAssn

   ! Get the element ID of the requested atom.
   currentElements(atomIndex) = atomTypes(atomType(atomIndex))%elementID

   ! Get the number of alphas for the requested atom.  The first number of
   !   orbital alphas is always the largest.
   numAlphas(atomIndex) = atomTypes(atomType(atomIndex))%numOrbAlphas(1)

   ! Get the alphas for the requested atom.
   alphas(:numAlphas(atomIndex),atomIndex) = atomTypes &
         & (atomType(atomIndex))%alphas(:numAlphas(atomIndex))

   ! Get the number of core states for this atom.
   numCoreStates(atomIndex) = atomTypes(atomType(atomIndex))%numCoreStates

   ! Get the number of valence states for this atom.
   numValeStates(atomIndex) = atomTypes(atomType(atomIndex))%numValeStates

   ! Determine the total number of states for this atom.
   numTotalStates(atomIndex) = numCoreStates(atomIndex) + &
         & numValeStates(atomIndex)

   ! Get the basis functions for this atom.
   currentBasisFns(:numAlphas(atomIndex), &
         & :numTotalStates(atomIndex),atomIndex) = &
         & basis(:numAlphas(atomIndex),:numTotalStates(atomIndex), &
         & atomType(atomIndex))

   ! Get the number of possible l,m combinations that a given
   !   atomic alpha will be used for.  (If an alpha for an atomic type i
   !   is used for s, p, and d orbital types, then the value at that matrix
   !   index is 9.  (Sum # states l=0,m=0; l=1,m=-1,0,1; l=2,m=-2,-1,0,1,2).)
   lmAlphaIndexCurrent(:numAlphas(atomIndex),atomIndex) = &
         & lmAlphaIndex(:numAlphas(atomIndex),atomType(atomIndex))

   ! Get the integer number uniquely identifying the kind of l,m
   !   combination that each state in the requested atom is in.  (If the
   !   state is l=0,m=0, then the value is 1.  If the state is l=2,m=2, then
   !   the value is 9.  If the state is l=1,m=0 then the value is 3.  If the
   !   state is l=1,m=1 then the value is 4.  etc.)
   lmIndexCurrent(:numTotalStates(atomIndex),atomIndex) = &
         & lmIndex(:numTotalStates(atomIndex),atomType(atomIndex))

   ! Get the x,y,z position of the requested atom.
   cartPos(:dim3,atomIndex) = atomSites(atomNumber)%cartPos(:dim3)

end subroutine initializeAtomSite

subroutine initializePotSite (potNumber,potType,potElement,zFactor,nucAlpha,&
      & cartPos)

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3
   use O_PotSites, only: numPotSites, potSites
   use O_PotTypes, only: numPotTypes, potTypes

   ! Input variables that request which atom to extract info from, and where
   !   to put it.
   integer :: potNumber ! Potential number

   ! Output varibales the will be filled in atomIndex array slot with data
   !   from atom number atomNumber.
   integer :: potType
   integer :: potElement
   real (kind=double) :: zFactor
   real (kind=double) :: nucAlpha
   real (kind=double), dimension (dim3) :: cartPos


   ! Fill the requested information for the potNumber provided.

   ! Get the type number of the requested potential nucleus.
   potType = potSites(potNumber)%potTypeAssn

   ! Get the element ID number of the requested potential nucleus.
   potElement = potTypes(potType)%elementID

   ! Get the nuclear charge associated with this type.
   zFactor = potTypes(potType)%nucCharge

   ! Get the exponential alpha factor for the nuclear potential.
   nucAlpha = potTypes(potType)%nucAlpha

   ! Get the position in the system cell for this nucleus.
   cartPos(:) = potSites(potNumber)%cartPos(:)

end subroutine initializePotSite

end module O_Basis
