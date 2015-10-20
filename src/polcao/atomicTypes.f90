module O_AtomicTypes

   ! Import necessary modules.
   use O_Kinds
   use O_Constants

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public
   private :: readRadialFns, countStatesAndFns

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module derived types.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Data structure to contain all the information relevant to an individual
   !   atomic type.  This structure is in turn indexed in a array of atom types
   !   in a higher level data structure.
   type AtomTypeType
      character*70 :: typeLabel  ! A label for this atomic type.  The label
            !   denotes the key information about this atom type and how it
            !   was defined.  If it is a core atom a 'C' will prefix the
            !   elemental name.  Otherwise it starts with the element name
            !   from the periodic table with appropriate capitalization (e.g.
            !   silicon is 'Si' not 'si').  This is followed immediately by
            !   a number indicating the species number of the atom type.  This
            !   is then followed by an '_' and then another number indicating
            !   the type number FOR THIS SPECIES (e.g. the first silicon of
            !   species 2 will look like Si2_1 and the second silicon of
            !   species 2 will look like Si2_2.  The 1 and 2 may not correspond
            !   to the overall type number of atom types in the system.)

      ! Identity information
      character*3 :: elementName
      integer :: elementID
      integer :: speciesID
      integer :: typeID

      ! Basis wave function information
      integer :: numCoreRadialFns ! The number of core radial wave functions.
      integer :: numValeRadialFns ! The number of valence radial wave functions.
      real (kind=double), pointer, dimension (:,:) :: coreRadialFns ! The radial
            !   wave functions of the core orbitals.  The first dimension holds
            !   the radial function coefficients while the second dimension
            !   indexes the number of core radial functions.
      real (kind=double), pointer, dimension (:,:) :: valeRadialFns ! The radial
            !   wave functions of the valence orbitals.  The first dimension
            !   holds the radial function coefficients while the second
            !   dimension indexes the number of valence radial functions.


      ! Alpha information
      integer, dimension(maxOrbitals) :: numOrbAlphas ! This holds the number
            !   gaussian terms (alphas) that each type of orbital
            !   (s,p,d,f,...) will use.
      real (kind=double), pointer, dimension (:) :: alphas ! The alphas
            !   used to describe the atomic wave function basis functions of
            !   all orbitals of this atomic type.  From exp(-alpha*r^2)


      ! Radial function information
      integer :: maxCoreQN_l  ! The highest angular momentum QN in all
            !   the core radial functions.
      integer :: maxValeQN_l  ! The highest angular momentum QN in all
            !   the valence orbitals.
      integer, pointer, dimension (:) :: coreQN_nList ! A list of the n
            !   principle QN for each core radial function.
      integer, pointer, dimension (:) :: valeQN_nList ! A list of the n
            !   principle QN for each valence radial function.
      integer, pointer, dimension (:) :: coreQN_lList ! A list of the l
            !   angular momentum values for each core radial function.
      integer, pointer, dimension (:) :: valeQN_lList ! A list of the l
            !   angular momentum values for each valence radial function.
      integer, pointer, dimension (:) :: coreQN_2jList ! A list of the 2*j
            !   QN for each core radial function.
      integer, pointer, dimension (:) :: valeQN_2jList ! A list of the 2*j
            !   QN for each valence radial function.
      integer, dimension(maxOrbitals) :: numQN_lCoreRadialFns ! This holds the
            !   number of s, p, d, and f radial functions in the core
            !   part of the wave function description.
      integer, dimension(maxOrbitals) :: numQN_lValeRadialFns ! This holds the
            !   number of s, p, d, and f radial functions in the valence
            !   part of the basis description.


      ! State information
      integer :: numCoreStates ! This is the number of states that are possible
            !   in the core part of the basis function.  It is basically the
            !   l quantum number times 2 plus 1 summed over the number of core
            !   radial basis functions.
      integer :: numValeStates ! The same as the above except for the valence.
      integer :: numTotalStates ! The sum of the number of states in the
            !   valence and the core.

   end type AtomTypeType

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer :: numAtomTypes ! The number of atomic types in the calc.
   type (AtomTypeType), allocatable, dimension (:) :: atomTypes ! Array of
         !   atomic type information for each type in the system.

   ! Statistics relating to atomic types.
   integer :: maxNumAtomAlphas ! This is the largest number of atomic
         !   alphas used by any atomic type in the system.
   integer :: maxNumCoreRadialFns ! This is the largest number of core
         !   radial functions used by any atomic type in the system.
   integer :: maxNumValeRadialFns ! This is the largest number of valence
         !   radial functions used by any atomic type in the system
   integer :: maxNumValeStates ! This is the largest number of valence states
         !   (2*l+1 for each radial function) used by any atomic type in
         !   the system.
   integer :: maxNumCoreStates ! This is the largest number of core states
         !   states (as above) used by any atomic type in the system.
   integer :: maxNumStates ! This is the largest total number of states
         !   (core and valence as above) used by any atomic type in the
         !   system.
   real (kind=double) :: minAtomicAlpha ! Smallest alpha from all the
         !   atomic alphas of all the atomic types.

   ! Statistics relating to the system as a whole.
   integer :: numElements

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine readAtomicTypes(readUnit,writeUnit)

   ! Bring in necessary modules.
   use O_Kinds
   use O_Constants, only: maxOrbitals
   use O_CommandLine, only: basisCode

   ! Import necessary subroutine modules.
   use O_ReadDataSubs

   ! Make sure that no funny variables are defined.
   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Define local variables.
   integer :: i,j
   integer, dimension (4) :: typeIDs
   integer, dimension (3) :: tempIntArray

   ! Read the number of atomic types.
   call readData(readUnit,writeUnit,numAtomTypes,len('NUM_ATOM_TYPES'),&
         & 'NUM_ATOM_TYPES')

   ! At this point we can allocate space for the atom type array and point the
   !   derived type pointer to it.
   allocate (atomTypes(numAtomTypes))

   ! Run a loop to read in all the atomic basis function information and
   !   allocate the associated arrays.
   do i = 1, numAtomTypes

      ! It it important to note that the type number is given based on the
      !   order that the types are given.  (i.e. The first type given is type
      !   number one...)

      ! Read the element/species/type ID numbers for this atom type.
      call readData(readUnit,writeUnit,4,typeIDs,&
            & len('ATOM_TYPE_ID__SEQUENTIAL_NUMBER'),&
            & 'ATOM_TYPE_ID__SEQUENTIAL_NUMBER')
      atomTypes(i)%elementID = typeIDs(1)
      atomTypes(i)%speciesID = typeIDs(2)
      atomTypes(i)%typeID    = typeIDs(3)

      ! The last value (typeIDs(4)) equals loop index i.
      if (typeIDs(4) /= i) then
         write (20,*) 'Types numbered out of order:  '
         write (20,*) 'i=',i,' typeIDs(4)=',typeIDs(4)
         stop
      endif

      ! Read the type label for the atom type and adjust the spacing.
      call readData(readUnit,writeUnit,70,atomTypes(i)%typeLabel,&
            & len('ATOM_TYPE_LABEL'),'ATOM_TYPE_LABEL')
      atomTypes(i)%typeLabel = trim (atomTypes(i)%typeLabel)

      ! Read the number of gaussian terms for the s,p,d,f radial basis
      !   functions.  They are indexed as 'alphas' from exp(-alpha*r^2).
      call readData(readUnit,writeUnit,maxOrbitals,&
            & atomTypes(i)%numOrbAlphas(1:maxOrbitals),&
            & len('NUM_ALPHA_S_P_D_F'),'NUM_ALPHA_S_P_D_F')

      ! Check that the number of alphas is monotonically decreasing.
      do j=1,maxOrbitals-1
         if (atomTypes(i)%numOrbAlphas(j) < &
               & atomTypes(i)%numOrbAlphas(j+1)) then
            write (20,*) 'Num of alphas not monotonically decreasing.'
            stop
         endif
      enddo

      ! Now we know how many alphas are going to be used for this type in the
      !   maximum case (numOrbAlphas(1)) so we can allocate space to hold them.
      allocate (atomTypes(i)%alphas(atomTypes(i)%numOrbAlphas(1)))

      ! Read the alphas for this atomic type.
      call readData(readUnit,writeUnit,atomTypes(i)%numOrbAlphas(1),&
            & atomTypes(i)%alphas(:atomTypes(i)%numOrbAlphas(1)),&
            & len('ALPHAS'),'ALPHAS')


      ! Read in the number of core radial basis functions for this atomic type.
      call readData(readUnit,writeUnit,3,tempIntArray(1:3),&
            & len('NUM_CORE_RADIAL_FNS'),'NUM_CORE_RADIAL_FNS')
      atomTypes(i)%numCoreRadialFns = tempIntArray(basisCode)

      ! Allocate space to hold the QN_n, QN_l, & QN_2j for each radial function.
      allocate (atomTypes(i)%coreQN_nList(atomTypes(i)%numCoreRadialFns))
      allocate (atomTypes(i)%coreQN_lList(atomTypes(i)%numCoreRadialFns))
      allocate (atomTypes(i)%coreQN_2jList(atomTypes(i)%numCoreRadialFns))

      ! Allocate space to hold the list of radial functions.
      allocate (atomTypes(i)%coreRadialFns(atomTypes(i)%numOrbAlphas(1),&
              & atomTypes(i)%numCoreRadialFns))

      ! Read in the core radial basis functions.  We pass this off to a
      !   subroutine because the procedure for reading in the functions
      !   for the core and valence is the same.  They just record the data
      !   to different data structures.
      call readRadialFns(readUnit, writeUnit,tempIntArray(3),&
                       & atomTypes(i)%numOrbAlphas,&
                       & atomTypes(i)%coreQN_nList,&
                       & atomTypes(i)%coreQN_lList,&
                       & atomTypes(i)%coreQN_2jList,&
                       & atomTypes(i)%coreRadialFns)

      ! Read in the num of valence radial basis functions for this atomic type.
      call readData(readUnit,writeUnit,3,tempIntArray(1:3),&
            & len('NUM_VALE_RADIAL_FNS'),'NUM_VALE_RADIAL_FNS')
      atomTypes(i)%numValeRadialFns = tempIntArray(basisCode)

      ! Allocate space to hold the QN_l for each radial function.
      allocate (atomTypes(i)%valeQN_nList(atomTypes(i)%numValeRadialFns))
      allocate (atomTypes(i)%valeQN_lList(atomTypes(i)%numValeRadialFns))
      allocate (atomTypes(i)%valeQN_2jList(atomTypes(i)%numValeRadialFns))

      ! Allocate space to hold the list of radial functions.
      allocate (atomTypes(i)%valeRadialFns(atomTypes(i)%numOrbAlphas(1),&
              & atomTypes(i)%numValeRadialFns))

      ! Read in the valence radial basis functions.
      call readRadialFns(readUnit,writeUnit,tempIntArray(3),&
                       & atomTypes(i)%numOrbAlphas,&
                       & atomTypes(i)%valeQN_nList,&
                       & atomTypes(i)%valeQN_lList,&
                       & atomTypes(i)%valeQN_2jList,&
                       & atomTypes(i)%valeRadialFns)
   enddo

end subroutine readAtomicTypes


subroutine readRadialFns(readUnit,writeUnit,numRadialFns,numOrbAlphas,&
            QN_nList,QN_lList,QN_2jList,radialFns)

   ! Import the necessary modules.
   use O_Kinds
   use O_Constants, only: maxOrbitals
   use O_CommandLine, only: basisCode

   ! Import necessary subroutine modules.
   use O_ReadDataSubs

   ! Make sure that there are not accidental variable declarations.
   implicit none


   ! Define dummy variables passed to this subroutine, in order.
   integer                             :: numRadialFns
   integer, dimension(maxOrbitals)     :: numOrbAlphas
   integer, dimension (:)              :: QN_nList
   integer, dimension (:)              :: QN_lList
   integer, dimension (:)              :: QN_2jList
   real (kind=double), dimension (:,:) :: radialFns
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.


   ! Define the local variables used in this routine
   integer                             :: i,j ! Loop counting integer
   integer, dimension(5)               :: tempIntArray
   real (kind=double), allocatable, dimension (:) :: tempRealArray

   ! Allocate space to hold the unused coefficients.
   allocate (tempRealArray(numOrbAlphas(1)))

   ! Next we read in the radial functions.  However, we must
   !   account for the case where there are no orbitals for this type.
   !   This may occur for an excited atom (no core) or a vacancy (no core and
   !   no valence).

   if (numRadialFns /= 0) then

      call readAndCheckLabel(readUnit,writeUnit,len('NL_RADIAL_FUNCTIONS'),&
            & 'NL_RADIAL_FUNCTIONS')

      i = 0
      do j = 1, numRadialFns
         ! Read the number of components for this radial function and the
         !   identifier for which basis set this radial function is a part of.
         !   At present, there can be only 1 component, and the basis set
         !   numbers are used as follows:  1=mb,fb,eb; 2=fb,eb; 3=eb
         call readData(readUnit,writeUnit,2,tempIntArray(1:2),0,'')

         ! Check if this radial function is in the requested basis.
         if (tempIntArray(2) <= basisCode) then

            ! Increment the counter of the number of radial functions read in.
            i = i + 1

            ! Read in the identifiers for this radial component.
            !  (QN_n, QN_l, QN_j*2, numStates in component, component index.)
            !  At present, only the QN_n and QN_l are used.
            call readData(readUnit,writeUnit,5,tempIntArray(1:5),0,'')
            QN_nList(i)  = tempIntArray(1)
            QN_lList(i)  = tempIntArray(2)
            QN_2jList(i) = tempIntArray(3)

            ! Read in the radial function coeffs for this n,l combination.
            call readData(readUnit,writeUnit,numOrbAlphas(QN_lList(i)+1),&
                  & radialFns(1:numOrbAlphas(QN_lList(i)+1),i),0,'')
         else
            ! Read past this radial function information and do not record it
            !   because it is used for a "higher" basis.
            call readData(readUnit,writeUnit,5,tempIntArray(1:5),0,'')
            call readData(readUnit,writeUnit,numOrbAlphas(tempIntArray(2)+1),&
                  & tempRealArray(1:numOrbAlphas(QN_lList(i)+1)),0,'')
         endif
      enddo
   endif

   deallocate (tempRealArray)
   
end subroutine readRadialFns


subroutine getAtomicTypeImplicitInfo

   ! Include necessary modules.
   use O_StringSubs

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the local variables used in this subroutine.
   integer :: i ! Loop index variables

   ! Initialize variables
   minAtomicAlpha      = 100 ! Some large integer
   maxNumAtomAlphas    = 0
   maxNumCoreRadialFns = 0
   maxNumValeRadialFns = 0
   maxNumCoreStates    = 0
   maxNumValeStates    = 0
   maxNumStates        = 0
   numElements         = 0

   do i = 1, numAtomTypes

      ! Obtain the element name for each atomic type.
      atomTypes(i)%elementName = getLeadChars(atomTypes(i)%typeLabel)

      ! Compare the smallest alpha for this type to the smallest alpha seen yet.
      minAtomicAlpha = min(minAtomicAlpha,atomTypes(i)%alphas(1))


      ! Compare the number of alphas for this type to the largest number of
      !   alphas seen for any atomic type yet.
      if (atomTypes(i)%numOrbAlphas(1) > maxNumAtomAlphas) then
         maxNumAtomAlphas = atomTypes(i)%numOrbAlphas(1)
      endif


      ! Track the number of s,p,d,f radial functions for the core and valence
      !   of this type.
      call countStatesAndFns (atomTypes(i)%numQN_lCoreRadialFns(:),&
            & atomTypes(i)%maxCoreQN_l, atomTypes(i)%numCoreStates,&
            & atomTypes(i)%numCoreRadialFns,atomTypes(i)%coreQN_lList(:))
      call countStatesAndFns (atomTypes(i)%numQN_lValeRadialFns(:),&
            & atomTypes(i)%maxValeQN_l, atomTypes(i)%numValeStates,&
            & atomTypes(i)%numValeRadialFns,atomTypes(i)%valeQN_lList(:))

      ! Sum the results of the counters for the number of valence and core
      !   states into one value.
      atomTypes(i)%numTotalStates = atomTypes(i)%numValeStates + &
            & atomTypes(i)%numCoreStates

      ! Compare results for this atomic type to the maximum values yet seen.

      ! Core first
      maxNumCoreStates = max(maxNumCoreStates,atomTypes(i)%numCoreStates)

      ! Valence next
      maxNumValeStates = max(maxNumValeStates,atomTypes(i)%numValeStates)

      ! Total last
      maxNumStates = max(maxNumStates,atomTypes(i)%numTotalStates)


      ! Compare the number of core radial functions for this type to the
      !   maximum value yet seen.
      if (atomTypes(i)%numCoreRadialFns > maxNumCoreRadialFns) then
         maxNumCoreRadialFns = atomTypes(i)%numCoreRadialFns
      endif

      ! Compare the number of valence radial functions for this type to the
      !   maximum value yet seen.
      if (atomTypes(i)%numValeRadialFns > maxNumValeRadialFns) then
         maxNumValeRadialFns = atomTypes(i)%numValeRadialFns
      endif

      ! Compare the element ID number for this type to the highest yet seen.
      if (atomTypes(i)%elementID > numElements) then
         numElements = atomTypes(i)%elementID
      endif
   enddo

   write (20,fmt="(a,e12.5)") 'Smallest Atomic    Alpha = ', minAtomicAlpha

end subroutine getAtomicTypeImplicitInfo


subroutine countStatesAndFns(numQN_lRadialFns,maxQN_l,numStates,numRadialFns,&
      & QN_lList)

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed dummy variables.
   integer, dimension (:) :: numQN_lRadialFns
   integer                :: maxQN_l
   integer                :: numStates
   integer                :: numRadialFns
   integer, dimension (:) :: QN_lList

   ! Define the local variables used in this subroutine.
   integer :: i  ! Loop index variable

   ! Initialize the s,p,d,f radial function counts to 0
   numQN_lRadialFns(:) = 0

   ! Initialize the maximum angular momentum value to an impossible #.
   maxQN_l = -1

   ! Initialize the number of states for this type to zero.
   numStates = 0

   ! Loop over the radial functions for this given type (core or valence).
   do i = 1, numRadialFns

      ! Increment the valence s,p,d,f function counter.
      numQN_lRadialFns(QN_lList(i)+1) = numQN_lRadialFns(QN_lList(i)+1) + 1

      ! Increment the total state counter.
      numStates = numStates + 2*QN_lList(i)+1

      ! Determine if this is the highest angular momentum value yet seen
      !   for this atomic type (core or valence).
      if (QN_lList(i) > maxQN_l) then
         maxQN_l = QN_lList(i)
      endif
   enddo

end subroutine countStatesAndFns


subroutine cleanUpRadialFns

   implicit none

   ! Define local variables.
   integer :: i

   do i = 1, numAtomTypes
      deallocate (atomTypes(i)%coreRadialFns)
      deallocate (atomTypes(i)%valeRadialFns)
   enddo

end subroutine cleanUpRadialFns


subroutine cleanUpAtomTypes

   implicit none

   ! Define local variables.
   integer :: i

   ! Deallocate all components of each atom type.
   do i = 1, numAtomTypes
      deallocate (atomTypes(i)%alphas)
      deallocate (atomTypes(i)%coreQN_nList)
      deallocate (atomTypes(i)%coreQN_lList)
      deallocate (atomTypes(i)%coreQN_2jList)
      deallocate (atomTypes(i)%valeQN_nList)
      deallocate (atomTypes(i)%valeQN_lList)
      deallocate (atomTypes(i)%valeQN_2jList)
   enddo

   ! Deallocate the list of atom types.
   deallocate (atomTypes)

end subroutine cleanUpAtomTypes


end module O_AtomicTypes
