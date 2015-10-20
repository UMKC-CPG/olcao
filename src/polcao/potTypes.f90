module O_PotTypes

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public
   private :: makePotAlphas

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module derived types.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   type PotTypeType
      character*70 :: typeLabel  ! This label functions just like the atomic
            !   type label.

      ! Identity information
      character*3 :: elementName
      integer :: elementID
      integer :: speciesID
      integer :: typeID


      ! Nuclear coulomb potential information
      real (kind=double) :: nucCharge ! The nuclear charge for this
            !   potential type according to the periodic table of the elements
            !   is stored here.
      real (kind=double) :: nucAlpha ! The nuclear potential is a Gaussian type
            !   function (exp(-alpha*r^2)) with this value as the alpha.


      ! Electronic coulomb potential information
      integer :: numAlphas ! This is just the number of potential terms for
            !   this potential type indexed by alphas in exp(-alpha*r^2).
      real (kind=double), pointer, dimension (:) :: alphas ! The input file
            !   defines a maximum and a minimum value for the alphas here.
            !   Then based on the number of alphas a geometric expression
            !   calculates the range of alphas that are stored here.
      integer :: cumulAlphaSum ! As with the cumulative sum of states that are
            !   seen in the atomic site description we also encounter a similar
            !   value here.  In this case the cumulative sum is stored with the
            !   potential types instead of the potential sites because it is
            !   not affected by the multiplicity (defined below).


      ! Exchange correlation potential information.
      real (kind=double) :: covalentRadius ! This is a value that should define
            !   the range at which covalent bonding is important for this
            !   potential type.  An accurate knowledge of this number for each
            !   potential type will help in determining the mesh for the
            !   exchange correlation potential, but for now it is impossible to
            !   know so this is always 1.0.


      ! Misc. information used commonly.
      integer :: multiplicity ! This stores the number of potential sites that
            !   share this particular potential type.

   end type PotTypeType

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer :: numPotTypes  ! The number of potential types in the calc.
   type (PotTypeType), allocatable, dimension (:) :: potTypes ! Array of
         !   potential type information for each type in the system.

   ! Statistics relating to potential types.
   integer :: maxNumPotAlphas ! This is the largest number of potential
         !   alphas used by any potential type in the system.
   real (kind=double) :: minPotAlpha ! Smallest alpha from all the
         !   potential alphas of all the potential types.
   real (kind=double) :: maxPotAlpha ! Largest alpha from all the potential
         !   alphas of all the potential types.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine readPotTypes(readUnit, writeUnit)

   ! Bring in necessary modules.
   use O_Kinds
   use O_ReadDataSubs

   ! Make sure that no funny variables are defined.
   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Define local variables.
   integer :: i
   integer, dimension (4) :: typeIDs
   real (kind=double), dimension (2) :: tempRealArray

   ! Read the number of potential types.
   call readData(readUnit,writeUnit,numPotTypes,len('NUM_POTENTIAL_TYPES'),&
         & 'NUM_POTENTIAL_TYPES')

   ! At this point we can allocate space for the potential type array.
   allocate (potTypes(numPotTypes))

   ! Run a loop to read in all the potential types.
   do i = 1, numPotTypes

      ! It it important to note that the type number is given based on the
      !   order that the types are given.  (i.e. The first type given is type
      !   number one...)

      ! Read the element/species/type ID numbers for this potential type.
      call readData(readUnit,writeUnit,4,typeIDs,&
            & len('POTENTIAL_TYPE_ID__SEQUENTIAL_NUMBER'),&
            & 'POTENTIAL_TYPE_ID__SEQUENTIAL_NUMBER')
      potTypes(i)%elementID = typeIDs(1)
      potTypes(i)%speciesID = typeIDs(2)
      potTypes(i)%typeID    = typeIDs(3)

      ! The last value (typeIDs(4)) equals loop index i.
      if (typeIDs(4) /= i) then
         write (20,*) 'Types numbered out of order:  '
         write (20,*) 'i=',i,' typeIDs(4)=',typeIDs(4)
         stop
      endif

      ! Read the type label for the potential type and adjust the spacing.
      call readData(readUnit,writeUnit,70,potTypes(i)%typeLabel,&
            & len('POTENTIAL_TYPE_LABEL'),&
            & 'POTENTIAL_TYPE_LABEL')
      potTypes(i)%typeLabel = trim (potTypes(i)%typeLabel)

      ! Read the nuclear charge and alpha for this type.
      call readData(readUnit,writeUnit,2,tempRealArray(1:2),&
            & len('NUCLEAR_CHARGE__ALPHA'),'NUCLEAR_CHARGE__ALPHA')
      potTypes(i)%nucCharge = tempRealArray(1)
      potTypes(i)%nucAlpha  = tempRealArray(2)

      ! Read the covalent radius for this type.
      call readData(readUnit,writeUnit,potTypes(i)%covalentRadius,&
            & len('COVALENT_RADIUS'),'COVALENT_RADIUS')

      ! Read the number of electronic potential alphas for this type.
      call readData(readUnit,writeUnit,potTypes(i)%numAlphas,&
            & len('NUM_ALPHAS'),'NUM_ALPHAS')

      ! Allocate space needed for the alphas for this type.
      allocate (potTypes(i)%alphas(potTypes(i)%numAlphas))

      ! Read the min and max alphas for this type into local temp variables.
      call readData(readUnit, writeUnit,2,tempRealArray(1:2),&
            & len('ALPHAS'),'ALPHAS')

      ! Create the list of alphas using a geometric series.
      call makePotAlphas(potTypes(i)%numAlphas,potTypes(i)%alphas,&
            & tempRealArray(1),tempRealArray(2))
   enddo

end subroutine readPotTypes


subroutine makePotAlphas(numAlphas,alphas,minAlpha,maxAlpha)

   ! Bring in necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define dummy parameters.
   integer                           :: numAlphas
   real (kind=double), dimension (:) :: alphas
   real (kind=double)                :: minAlpha
   real (kind=double)                :: maxAlpha

   ! Define local variables.
   integer :: i

   do i = 1, numAlphas

      ! Make the alpha value: minA*[(maxA/minA)**(1/numAlphas-1)]**(i-1)
      alphas(i) = minAlpha * ((maxAlpha/minAlpha) ** &
         & (1.0_double/real(numAlphas-1,double))) ** (i-1)
   enddo

end subroutine makePotAlphas


subroutine getPotTypeImplicitInfo

   ! Include necessary modules.
   use O_StringSubs
   use O_PotSites ! For getMultiplicities subroutine

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the local variables used in this subroutine.
   integer :: i
   integer :: potDimTemp  ! potDim is recomputed later in the potential module.
   integer, allocatable, dimension(:) :: multiplicities

   ! Initialize variables
   maxNumPotAlphas = 0
   minPotAlpha     = 10000 ! Some large integer
   maxPotAlpha     = 0
   potDimTemp      = 0

   ! Get the multiplicity for all potential types with one pass.  Later, this
   !   value will be copied to individual potential types.
   allocate (multiplicities(numPotTypes))
   call getMultiplicities(multiplicities)

   ! Obtain the following information:
   !   Element Name for each type.
   !   Number of potential sites for each type.
   !   The total number of electronic potential terms (alphas) in the system.
   !   The total number of terms prior to each type in the total list of types.
   !   The max number of potential terms (alphas) of any type.
   !   The smallest potential alpha of any type.
   !   The largest potential alpha of any type.
   do i = 1, numPotTypes

      ! Obtain the element name for each potential type.
      potTypes(i)%elementName = getLeadChars(potTypes(i)%typeLabel)

      ! Get the multiplicity of this potential type.
      potTypes(i)%multiplicity = multiplicities(i)

      ! Record the cumulative value of the number of potential alphas defined.
      potTypes(i)%cumulAlphaSum = potDimTemp

      ! Add to the summation by the amount of alphas used for this type
      potDimTemp = potDimTemp + potTypes(i)%numAlphas

      ! Check if the number of alphas for this type is the most yet seen.
      if (potTypes(i)%numAlphas > maxNumPotAlphas) then
         maxNumPotAlphas = potTypes(i)%numAlphas
      endif

      ! Compare the smallest alpha for this type to the smallest alpha seen yet.
      minPotAlpha = min(minPotAlpha, potTypes(i)%alphas(1))

      ! Compare the largest alpha for this type to the largest alpha seen yet.
      maxPotAlpha = max(maxPotAlpha,potTypes(i)%alphas(potTypes(i)%numAlphas))
   enddo

   deallocate (multiplicities)

   ! Log the important statistics
   write (20,fmt="(a,i5)")    'Potential Dimension      = ', potDimTemp
   write (20,fmt="(a,e12.5)") 'Smallest Potential Alpha = ', minPotAlpha
   write (20,fmt="(a,e12.5)") 'Largest  Potential Alpha = ', maxPotAlpha

end subroutine getPotTypeImplicitInfo


subroutine cleanUpPotTypes

   implicit none

   ! Define local variables.
   integer :: i

   ! Deallocate all potential type components.
   do i = 1, numPotTypes
      deallocate (potTypes(i)%alphas)
   enddo

   ! Deallocate the structure.
   deallocate (potTypes)

end subroutine cleanUpPotTYpes


end module O_PotTypes
