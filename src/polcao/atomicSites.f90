module O_AtomicSites

   ! Import necessary modules.
   use O_Kinds
   use O_Constants

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module derived types.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   type AtomSiteType

      ! Atom data
      integer :: atomTypeAssn ! This defines which atomic type this atom site
            !   is associated with.
      real (kind=double), dimension (dim3) :: cartPos ! This holds the three
            !   cartesean coordinates of the current atomic site.

      ! Cumulative state counters.
      integer :: cumulCoreStates ! The number of core states is summed over
            !   each atom in the system.  The type of each atom site is
            !   determined and the number of core states for that type is added
            !   to a running total.  The running total for each atom index is
            !   stored here.
      integer :: cumulValeStates ! The same as the above except for the valence.
   end type AtomSiteType

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer :: numAtomSites  ! The number of atomic sites in the calc.
   type (AtomSiteType), allocatable, dimension (:) :: atomSites ! Array of
         !   atomic site information for each site in the system.

   ! Statistics relating to atoms.
   integer :: valeDim, coreDim ! Valence and core dimension determined as
         !   s=1,p=3,d=5,f=7 summed orbital types for all atoms.
   integer :: maxDim ! The maximum of the above two dimensions.


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine readAtomicSites(readUnit,writeUnit)

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
   integer :: i
   integer :: counter
   character*2 :: atomName

   ! Read the number of atomic sites.
   call readData(readUnit,writeUnit,numAtomSites,len('NUM_ATOM_SITES'),&
         & 'NUM_ATOM_SITES')

   ! The number of atomic sites is known so we can allocate space to hold
   !   the site data structure.
   allocate (atomSites(numAtomSites))

   call readAndCheckLabel(readUnit,writeUnit,len('NUM_TYPE_X_Y_Z_ELEM'),&
         & 'NUM_TYPE_X_Y_Z_ELEM')

   do i = 1, numAtomSites
      read (4,*)     counter,atomSites(i)%atomTypeAssn,atomSites(i)%cartPos,&
            & atomName
      write (20,900) counter,atomSites(i)%atomTypeAssn,atomSites(i)%cartPos,&
            & atomName

      ! Check that the list is in order
      if (counter /= i) then
         write (20,*) 'Potential site list is out of order at i = ',i
         stop
      endif
   enddo

   900 format (2(i5,1x),3(f18.8,1x),a2)         ! Site coordinates

end subroutine readAtomicSites


subroutine getAtomicSiteImplicitInfo

   ! Include necessary object modules.
   use O_AtomicTypes, only: atomTypes

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the local variables used in this subroutine.
   integer :: i ! Loop index variables
   integer :: currentType

   ! Initialize the core and valence dimensions defined by the number of core
   !   or valence states for each atomic type multiplied by the number of
   !   atoms for that type.
   valeDim = 0
   coreDim = 0
   maxDim  = 0

   ! Loop over the number of atoms in the system.  At each iteration:
   ! 1) Increment the coreDim and valeDim by the number of states for the
   !    current atom.
   ! 2) Store that value with that atom's data structure.
   do i = 1, numAtomSites

      ! Get the current type of this atom
      currentType = atomSites(i)%atomTypeAssn

      ! Store the current summation value with the current atom's data.  Note
      !   that the cumulative sum is off by one atom so that atom 1 has a sum
      !   of zero.
      atomSites(i)%cumulValeStates = valeDim
      atomSites(i)%cumulCoreStates = coreDim

      ! Cumulatively add the valence and core dimensions of each atom
      valeDim = valeDim + atomTypes(currentType)%numValeStates
      coreDim = coreDim + atomTypes(currentType)%numCoreStates

   enddo

   ! Determine the maximum dimension between the valence and core.
   maxDim = max(valeDim,coreDim)

   ! Log these statistics
   write (20,fmt="(a,i5)") 'Valence Dimension   = ',valeDim
   write (20,fmt="(a,i5)") 'Core    Dimension   = ',coreDim
   write (20,fmt="(a,i5)") 'Maximum Dimension   = ',maxDim

end subroutine getAtomicSiteImplicitInfo


subroutine cleanUpAtomSites

   implicit none

   deallocate (atomSites)

end subroutine cleanUpAtomSites


end module O_AtomicSites
