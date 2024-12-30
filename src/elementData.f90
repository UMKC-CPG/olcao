module O_ElementData

   ! Use any necessary modules.
   use O_Kinds

   ! Make sure no implicit variables are defined.
   implicit none

   ! Define module data.
   character*3, allocatable, dimension (:) :: elementNames ! The
         ! abbreviated names of the elements from the periodic table.
   character*12, allocatable, dimension (:) :: elementFullNames ! The full
         ! unabbreviated names of the elements from the periodic table.
   integer, allocatable, dimension (:) :: numUJElectrons ! Number
         ! of electrons in the highest occupied d or f orbital. This data was
         ! taken from a periodic table of the elements and should be improved
         ! upon in some way for more complicated scenarios.
   real (kind=double), allocatable, dimension (:) :: atomicMass ! Atomic mass
         ! of each unique element.
   real (kind=double), allocatable, dimension (:) :: covalRadii ! Covalent
         ! radius of each unique element.
   real (kind=double), allocatable, dimension (:) :: atomicRadii ! Atomic
         ! radius of each unique element.
   real (kind=double), allocatable, dimension (:) :: neutScatt ! Neutron
         ! scattering factor of each unique element.
   real (kind=double), allocatable, dimension (:,:) :: coreCharge ! Number of
         ! core electrons in each s,p,d,f orbital of each unique element.
   real (kind=double), allocatable, dimension (:,:) :: valeCharge ! Number of
         ! valence electrons in each s,p,d,f orbital of each unique element.
   real (kind=double), allocatable, dimension (:,:) :: colorVTK ! Default
         ! color values to use in VTK, based on CPK chemistry color scheme.
   real (kind=double), allocatable, dimension (:) :: colorDX ! Default color
         ! values to use in OpenDX on a scale of 1-100 for each unique element.
   real (kind=double), allocatable, dimension (:) :: greyDX ! Default grey
         ! scale values to use in OpenDX on a scale of 1-100 for each unique
         ! element.
   integer :: numUniqueElements ! Number of elements with some physical data
         ! from the periodic table of the elements that are included in OLCAO.

   contains

! This subroutine will read a bunch of data from the elements.dat file. Much
!   of the data in elements.dat is already known by the olcao program upon
!   reading the input file or it is not specifically needed by the olcao
!   program. That information will be skipped. It is important that any new
!   information that is added to elements.dat not disrupt the order present in
!   that file. If such a disruption is caused, then this subroutine and other
!   programs (e.g. perl scripts) that read elements.dat will need to be
!   updated.
subroutine initElementData

   ! Use any necessary modules.
   use O_Kinds
   use O_Constants, only: lAngMomCount

   ! Make sure no implicit variables are defined.
   implicit none

   ! Define local variables.
   character*100 :: dataDirectory
   character*100 :: elementDataFile
   integer   :: info
   integer   :: i

   ! Determine if the elemental data has already been read in or not. If this
   !   data has already been read in, then we abort and do not read it in
   !   again.
   if (allocated(elementNames)) then
      return
   endif

   ! Determine the file name that contains the elemental data information.
   call get_environment_variable(NAME="OLCAO_DATA",VALUE=dataDirectory,&
         & STATUS=info)
   if (info /= 0) then
      write (20,*) "info=",info
      stop "Failed to get OLCAO_DATA environment variable."
   endif
   elementDataFile=trim(dataDirectory)//"/elements.dat"

   ! Open the elementDataFile for reading.
   open (unit=313,file=elementDataFile,form='formatted',status='old',&
         & IOSTAT=info)
   if (info /= 0) then
      stop "Failed to open elementDataFile elements.dat"
   endif

   ! Read the number of elements in the data base.
   read (313,*)
   read (313,*) numUniqueElements

   ! Read the maximum QN_L quantum number, but we ignore this information for
   !   now because this number is already known in "O_Constants".
   read (313,*)
   read (313,*)

   ! Allocate arrays for holding data.
   allocate (elementNames(numUniqueElements))
   allocate (elementFullNames(numUniqueElements))
   allocate (atomicMass(numUniqueElements))
   allocate (covalRadii(numUniqueElements))
   allocate (atomicRadii(numUniqueElements))
   allocate (neutScatt(numUniqueElements))
   allocate (numUJElectrons(numUniqueElements))
   allocate (coreCharge(lAngMomCount,numUniqueElements))
   allocate (valeCharge(lAngMomCount,numUniqueElements))
   allocate (colorVTK(4,numUniqueElements)) ! RGB + alpha
   allocate (colorDX(numUniqueElements))
   allocate (greyDX(numUniqueElements))

   ! Assign data values for the element names.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*) elementNames(i)
   enddo

   ! Assign data values for the full element names.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*) elementFullNames(i)
   enddo

   ! Atomic mass of each element.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*) atomicMass(i)
   enddo

   ! Covalent Radii of each element.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*) covalRadii(i)
   enddo

   ! Atomic Radii of each element.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*) atomicRadii(i)
   enddo

   ! Neutron scattering factor of each element.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*) neutScatt(i)
   enddo

   ! Number of UJ electrons (highest d or f orbital).
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*) numUJElectrons(i)
   enddo

   ! Read past the LJ Pair Coefficients for each element.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*)
   enddo

   ! Read past the definition of the number of core orbitals of each spdf type.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*)
   enddo

   ! Read past the definition of the number of minimal basis valence orbitals
   !   of each spdf type beyond the core.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*)
   enddo

   ! Read past the definition of the number of full basis valence orbitals
   !   of each spdf type beyond the minimal basis.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*)
   enddo

   ! Read past the definition of the number of extended basis valence orbitals
   !   of each spdf type beyond the full basis.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*)
   enddo

   ! Define the number of electrons in the occupied core orbitals.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*) coreCharge(:,i)
   enddo

   ! Define the number of electrons in the occupied orbitals above the core.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*) valeCharge(:,i)
   enddo

   ! Read past the definition of the number of radial wave function terms for
   !   each element.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*)
   enddo

   ! Read past the definition of the minimum exponential alpha for the radial
   !   wave function terms for each element.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*)
   enddo

   ! Read past the definition of the maximum exponential alpha for the radial
   !   wave function terms for each element.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*)
   enddo

   ! Read past the definition of the number of potential terms for each element.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*)
   enddo

   ! Read past the definition of the minimum exponential alpha for the
   !   potential function terms for each element.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*)
   enddo

   ! Read past the definition of the maximum exponential alpha for the
   !   potential function terms for each element.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*)
   enddo

   ! Read past the definition of which terms to use in the Gaussian expansion
   !   of each orbital type (spdf) for each element.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*)
      read (313,*)
      read (313,*)
      read (313,*)
   enddo

   ! Read the definition of the VTK color+alpha values for each element.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*) colorVTK(:,i)
   enddo

   ! Read the definition of openDX color values for each element.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*) colorDX(i)
   enddo

   ! Read the definition of openDX grey scale values for each element.
   read (313,*)
   do i = 1, numUniqueElements
      read (313,*) greyDX(i)
   enddo

   ! Close the element data file.
   close (313)

end subroutine initElementData



subroutine deallocateElementData

   ! Make sure that no implicit variables are defined.
   implicit none

   ! Deallocate arrays used to hold data.
   deallocate (elementNames)
   deallocate (atomicMass)
   deallocate (covalRadii)
   deallocate (atomicRadii)
   deallocate (neutScatt)
   deallocate (numUJElectrons)
   deallocate (coreCharge)
   deallocate (valeCharge)
   deallocate (colorVTK)
   deallocate (colorDX)
   deallocate (greyDX)

end subroutine deallocateElementData


end module O_ElementData
