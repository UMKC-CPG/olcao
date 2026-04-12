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

   ! Ionic data.
   real (kind=double), dimension(dim3) :: xyzIonMoment
   real (kind=double), dimension(dim3) :: abcIonMoment

   ! Atom permutation table for IBZ unfolding.
   integer, allocatable, dimension (:,:) :: atomPerm
         !   For each point group operation R and atom
         !   site A, atomPerm(R, A) = B gives the atom
         !   site B that operation R maps A onto. Used by
         !   computeBond to distribute eigenvector-dependent
         !   properties (effective charge, bond order) across
         !   the star of each IBZ k-point via the atom
         !   permutation fix (DESIGN 2.4). Built by
         !   buildAtomPerm after point group operations and
         !   atom positions are both available. Dimensions:
         !   (numPointOps, numAtomSites). Only allocated for
         !   k-point style codes 1 and 2.

   ! Inverse atom permutation for IBZ unfolding.
   integer, allocatable, dimension (:,:) :: invAtomPerm
         !   For each point group operation R and atom
         !   site B, invAtomPerm(R, B) = A gives the atom
         !   site A such that atomPerm(R, A) = B, i.e.,
         !   A = R^{-1}(B). Used during LAT PDOS tetrahedron
         !   corner assembly to map channel indices from
         !   full-mesh k-points back to their IBZ
         !   representatives (DESIGN 1.4, PSEUDOCODE 4a).
         !   Built by buildInvAtomPerm immediately after
         !   buildAtomPerm. Dimensions: (numPointOps,
         !   numAtomSites). Only allocated when atomPerm
         !   is allocated.


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine readAtomicSites(readUnit,writeUnit)

   ! Import necessary subroutine modules.
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


subroutine computeIonicMoment(incChargePerVol)

   use O_Kinds
   use O_Constants, only: eCharge,pi
   use O_PotTypes, only: potTypes
   use O_Lattice, only: getRealmagnitudes,realCellVolume,invRealVectors,realMag
   use O_ElementData, only: coreCharge

   implicit none

   ! Define passed parameters.
   integer, intent(in) :: incChargePerVol

   ! Define local variables.
   integer :: i,j,k
   real (kind=double) :: currCoreCharge

   call getRealMagnitudes

   xyzIonMoment(:) = 0.0_double
   abcIonMoment(:) = 0.0_double

   do i = 1, numAtomSites
      if (coreDim == 0) then
         currCoreCharge = 0.0_double
      else
         currCoreCharge = sum(coreCharge(:,&
               & int(potTypes(atomSites(i)%atomTypeAssn)%nucCharge)))
      endif

!write(20,*) "eCharge/realCellVol = ", eCharge / realCellVolume
!         xyzIonMoment(:) = xyzIonMoment(:) + eCharge / realCellVolume * &
!               & (potTypes(atomSites(i)%atomTypeAssn)%nucCharge - &
!               & currCoreCharge) * atomSites(i)%cartPos(:)
      xyzIonMoment(:) = xyzIonMoment(:) + &
            & (potTypes(atomSites(i)%atomTypeAssn)%nucCharge - &
            & currCoreCharge) * atomSites(i)%cartPos(:)
!write(20,*) "Z = ", (potTypes(atomSites(i)%atomTypeAssn)%nucCharge - &
!            & currCoreCharge)
!write(20,*) "xyzPos = ", atomSites(i)%cartPos(:)
!write(20,*) "xyzIonMom = ", xyzIonMoment(:)
   enddo

   if (incChargePerVol == 1) then
      xyzIonMoment(:) = xyzIonMoment(:) * eCharge / realCellVolume
   endif

   ! Convert xyz moment into abc coordinates.
   do j = 1,3 ! abc axes
      abcIonMoment(j) = 0.0_double
      do k = 1,3 ! xyz axes
         abcIonMoment(j) = abcIonMoment(j) + &
               & xyzIonMoment(k)*invRealVectors(k,j)
      enddo
      abcIonMoment(j) = abcIonMoment(j)*realMag(j)
   enddo
!write(20,*) "abcIonMom = ", abcIonMoment(:)

   ! Convert from atomic units of distance (bohr radii).

end subroutine computeIonicMoment


subroutine buildAtomPerm

   ! Build the atom permutation table atomPerm(R, A) = B, which records that
   ! point group operation R maps atom A onto atom B. This is the single piece
   ! of infrastructure needed for correct IBZ unfolding of all shell-summed
   ! quantities: effective charge (Q*), bond order, and PDOS modes 0-2. See
   ! PSEUDOCODE section 4 and DESIGN section 2.4.
   !
   ! The algorithm works in fractional (abc) coordinates because the point
   ! group operations are stored in that basis. Cartesian atom positions are
   ! converted to fractional using invRealVectors (= recipVectors / 2*pi).
   ! For each operation R and each atom A, the rotation is applied to A's
   ! fractional position, the result is wrapped into [0,1), and the matching
   ! atom B of the same species is identified by minimum-image comparison.

   ! Import the point group rotation matrices and the fractional translation
   !   vectors needed for IBZ symmetry reduction and real-space atom mapping.
   use O_KPoints, only: numPointOps, abcPointOps, abcFracTrans

   ! Import the inverse real-space lattice vectors for the Cartesian-to-
   !   fractional coordinate conversion: abc = invRealVectors * xyz.
   use O_Lattice, only: invRealVectors

   ! Make sure no accidental variables are defined.
   implicit none

   ! Tolerance for fractional coordinate matching. Atoms whose wrapped
   !   fractional positions differ by less than this threshold on every
   !   axis are considered the same crystallographic site.
   real (kind=double), parameter :: posTol = 1.0e-5_double

   ! Define local variables.
   integer :: opIdx      ! Loop index over point group operations.
   integer :: atomA      ! Loop index: source atom being rotated.
   integer :: atomB      ! Loop index: candidate target atom.
   integer :: axis       ! Loop index over fractional coordinate axes.
   logical :: matchFound ! True when atom B matching R(atomA) is found.

   ! Fractional (abc) coordinates of every atom site, converted from the
   !   Cartesian positions stored in atomSites(:)%cartPos. The conversion
   !   uses the same convention as computeIonicMoment: abc(i) = sum_j
   !   invRealVectors(i,j) * xyz(j). Dimensions: (3, numAtomSites).
   real (kind=double), allocatable, dimension (:,:) :: abcAtomPos

   ! The fractional position of atom A after the full space group operation
   !   {R|t}: rotPos(i) = sum_j abcPointOps(i,j,R) * abcAtomPos(j,A) +
   !   abcFracTrans(i,R).
   real (kind=double), dimension (dim3) :: rotPos

   ! Difference between the rotated position and a candidate atom's
   !   fractional position, wrapped into [-0.5, 0.5) via nint to enforce
   !   the minimum-image convention under periodic boundaries.
   real (kind=double), dimension (dim3) :: diff

   ! -----------------------------------------------------------------
   ! Step 1: Convert all atom positions from Cartesian (xyz) to
   !   fractional (abc) coordinates. The transformation is:
   !     abcAtomPos(i, A) = sum_j invRealVectors(i,j) * cartPos(j)
   !   This is the same convention used in computeIonicMoment.
   ! -----------------------------------------------------------------

   allocate (abcAtomPos(dim3, numAtomSites))

   do atomA = 1, numAtomSites
      do axis = 1, dim3
         abcAtomPos(axis, atomA) = &
               & sum(invRealVectors(axis,:) &
               &     * atomSites(atomA)%cartPos(:))
      enddo
   enddo

   ! -----------------------------------------------------------------
   ! Step 2: Allocate the permutation table. Each entry will be filled
   !   in Step 3. The table has one row per point group operation and
   !   one column per atom site.
   ! -----------------------------------------------------------------

   allocate (atomPerm(numPointOps, numAtomSites))

   ! -----------------------------------------------------------------
   ! Step 3: For each operation R and each atom A, apply R to the
   !   fractional position of A, wrap the result into [0,1), and find
   !   the matching atom B of the same species. Point group operations
   !   preserve species (they are pure rotations/reflections of the
   !   crystal), so only atoms sharing the same atomTypeAssn can match.
   ! -----------------------------------------------------------------

   do opIdx = 1, numPointOps
      do atomA = 1, numAtomSites

         ! Apply the full space group operation {R|t} in fractional coordinates:
         !   rotPos = R*abc + t (R = rotation matrix, t = fractional translation).
         do axis = 1, dim3
            rotPos(axis) = &
                  & sum(abcPointOps(axis,:,opIdx) &
                  &     * abcAtomPos(:,atomA)) &
                  & + abcFracTrans(axis, opIdx)
         enddo

         ! Wrap the rotated position into [0, 1). The Fortran modulo
         !   intrinsic correctly handles negative arguments (unlike mod).
         do axis = 1, dim3
            rotPos(axis) = modulo(rotPos(axis), 1.0_double)
         enddo

         ! Search for the atom at the rotated position. Only atoms of
         !   the same type can match (species is conserved by symmetry).
         matchFound = .false.
         do atomB = 1, numAtomSites

            ! Skip atoms of a different species.
            if (atomSites(atomB)%atomTypeAssn /= &
                  & atomSites(atomA)%atomTypeAssn) cycle

            ! Compute the minimum-image difference on each fractional
            !   axis: subtract, then wrap into [-0.5, 0.5) via nint.
            do axis = 1, dim3
               diff(axis) = rotPos(axis) &
                     & - abcAtomPos(axis, atomB)
               diff(axis) = diff(axis) &
                     & - nint(diff(axis))
            enddo

            ! Check whether all components fall within tolerance.
            if (all(abs(diff(:)) < posTol)) then
               atomPerm(opIdx, atomA) = atomB
               matchFound = .true.
               exit
            endif

         enddo ! atomB

         ! Safety check: every atom must map to exactly one partner
         !   under every point group operation. A missing match means
         !   the point group operations and the atom positions are
         !   inconsistent -- this is a fatal error.
         if (.not. matchFound) then
            write (20,*) 'ERROR in buildAtomPerm:'
            write (20,*) '  No matching atom for site ', &
                  & atomA, ' under operation ', opIdx
            write (20,*) '  Rotated fractional pos = ', &
                  & rotPos(:)
            stop 'buildAtomPerm: no atom match found'
         endif

      enddo ! atomA
   enddo ! opIdx

   ! Clean up the temporary fractional coordinate array.
   deallocate (abcAtomPos)

end subroutine buildAtomPerm


! Build the inverse atom permutation table. For each
!   point group operation R and atom site B,
!   invAtomPerm(R, B) = A where atomPerm(R, A) = B.
!   This gives A = R^{-1}(B): the atom at the IBZ
!   k-point whose image under R is atom B at the
!   full-mesh k-point. Required for LAT PDOS channel
!   permutation (DESIGN 1.4, PSEUDOCODE 4a).
!
!   Must be called after buildAtomPerm, which
!   populates the forward table atomPerm. The inverse
!   is constructed by scanning atomPerm and swapping
!   the role of A and B.
subroutine buildInvAtomPerm

   use O_KPoints, only: numPointOps

   implicit none

   ! Local loop indices.
   integer :: opIdx  ! Point group operation index.
   integer :: atomA  ! Source atom (forward map).
   integer :: atomB  ! Target atom (forward map).

   ! Allocate the inverse table with the same shape
   !   as atomPerm: (numPointOps, numAtomSites).
   allocate (invAtomPerm(numPointOps, numAtomSites))

   ! Invert: if atomPerm(R, A) = B then
   !   invAtomPerm(R, B) = A.
   do opIdx = 1, numPointOps
      do atomA = 1, numAtomSites
         atomB = atomPerm(opIdx, atomA)
         invAtomPerm(opIdx, atomB) = atomA
      enddo
   enddo

end subroutine buildInvAtomPerm


subroutine cleanUpAtomSites

   implicit none

   deallocate (atomSites)

   ! Only allocated for k-point style codes 1 and 2.
   if (allocated(atomPerm)) then
      deallocate (atomPerm)
   endif
   if (allocated(invAtomPerm)) then
      deallocate (invAtomPerm)
   endif

end subroutine cleanUpAtomSites


end module O_AtomicSites
