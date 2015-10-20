module O_PotSites

   ! Import necessary modules.
   use O_Kinds
   use O_Constants, only: dim3

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module derived types.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   type PotSiteType
      integer :: potTypeAssn ! Defines which potential type this potential site
            !   is associated with.
      integer :: firstPotType ! Defines if this potential site is the first
            !   site with of its potential type or not.  If it is, then
            !   this value is 1.  If it is not, then this value is 0.
      real (kind=double), dimension (dim3) :: cartPos ! This holds the three
            !   cartesean coordinates of the current potential site.
   end type PotSiteType

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer :: numPotSites  ! The number of potential sites in the calc.
   type (PotSiteType), allocatable, dimension (:) :: potSites ! Array of
         !   potential site information for each site in the system.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine readPotSites(readUnit, writeUnit)

   ! Import necessary modules.
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
   integer :: counter
   character*2 :: atomName

   ! Read the number of potential sites.
   call readData(readUnit,writeUnit,numPotSites,len('NUM_POTENTIAL_SITES'),&
         & 'NUM_POTENTIAL_SITES')

   ! The number of potential sites is known so we can allocate space to hold
   !   the potential type associations, equivalencies, and coordinates.
   allocate (potSites(numPotSites))

   call readAndCheckLabel(readUnit,writeUnit,len('NUM_TYPE_X_Y_Z_ELEM'),&
         & 'NUM_TYPE_X_Y_Z_ELEM')

   do i = 1, numPotSites
      read (4,*)     counter,potSites(i)%potTypeAssn,potSites(i)%cartPos,&
            & atomName
      write (20,900) counter,potSites(i)%potTypeAssn,potSites(i)%cartPos,&
            & atomName

      ! Check that the list is in order
      if (counter /= i) then
         write (20,*) 'Potential site list is out of order at i = ',i
         stop
      endif
   enddo

   900 format (2(i5,1x),3(f18.8,1x),a2)         ! Site coordinates

end subroutine readPotSites


subroutine getPotSiteImplicitInfo

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define local variables.
   integer :: i
   integer :: currentType

   ! Determine which potential sites are the first ones in the ordered list
   !   to be of a particular type.  (Recall that the input is already sorted
   !   by element/species/type.)
   potSites(1)%firstPotType = 1
   currentType = 1

   if (numPotSites > 1) then
      do i = 2, numPotSites

         if (potSites(i)%potTypeAssn > currentType) then
            currentType = currentType + 1
            potSites(i)%firstPotType = 1
         else
            potSites(i)%firstPotType = 0
         endif
      enddo
   endif

end subroutine getPotSiteImplicitInfo


subroutine getMultiplicities (multiplicities)

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed parameters.
   integer, dimension (:) :: multiplicities

   ! Define local variables.
   integer :: i

   ! Initialize the multiplicity counters
   multiplicities(:) = 0

   ! For each site increment the count of the multiplicity for the associated
   !   potential type.
   do i = 1, numPotSites
      multiplicities(potSites(i)%potTypeAssn) = &
            & multiplicities(potSites(i)%potTypeAssn) + 1
   enddo

end subroutine getMultiplicities


subroutine cleanUpPotSites

   implicit none

   deallocate (potSites)

end subroutine cleanUpPotSites


end module O_PotSites
