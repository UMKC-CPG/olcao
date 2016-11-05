module O_SpaceGroupOperations

   ! Import necessary modules
   use O_Kinds

   ! Make sure no funny variables are defined.
   implicit none

   ! Begin list of module data.
   character*1 :: spaceLattice
   integer :: spaceGroupNumber
   integer :: numSpaceOps
   integer :: numShifts ! This variable is essentially unused in this program
         ! because all space group operations are treated equally and the
         ! concept of a shift or no shift is not important.
   real (kind=double), allocatable, dimension (:,:,:) :: spaceOps ! Index 1 =
         ! abc contributions to one axis (1:3).  Index 2 = the axes contributed
         ! to by the values in index 1 (1:3).  Index 3 = the number of space
         ! group operations.  (i.e. the new 'a' (index 2, slot 1) coordinate
         ! depends on the values in the old abc according to the relationships
         ! specified in index 1.
   real (kind=double), allocatable, dimension (:,:)   :: spaceShifts ! Index 1=
         ! the shift applied to each abc axis.  Index 2 = the number of space
         ! group operations.

   ! Begin list of module functions.
   contains


   subroutine readSpaceOps (fileUnit)

      ! Make sure no funny variables are defined.
      implicit none

      ! Define dummy variables.
      integer, intent (in) :: fileUnit

      ! Define local variables.
      integer :: i

      ! Read the character identifier of the space group lattice type.
      read (fileUnit,fmt="(a1)") spaceLattice

      ! Read the root space group number for this space group.
      read (fileUnit,*) spaceGroupNumber

      ! Read the number of space group operations and number of shifted
      !   repetitions that exist in that list.
      read (fileUnit,*) numSpaceOps,numShifts

      ! Initialize storage space for the space group operations and the shift
      !   operations associated with each space group operation.
      allocate (spaceOps(3,3,numSpaceOps))
      allocate (spaceShifts(3,numSpaceOps))

      ! Read the cartesian form of the space group operations from stdin.
      do i = 1, numSpaceOps
         read (fileUnit,*) ! Read empty line preceeding operations.
         read (fileUnit,*) spaceOps(:,1,i)  ! a coordinate components
         read (fileUnit,*) spaceOps(:,2,i)  ! b coordinate components
         read (fileUnit,*) spaceOps(:,3,i)  ! c coordinate components
         read (fileUnit,*) spaceShifts(:,i) ! shifts to apply to a,b,c
      enddo
   end subroutine readSpaceOps


   subroutine cleanSpaceGroupOperations

      ! Make sure no funny variables are defined.
      implicit none

      ! Deallocate memory that is no longer needed.
      deallocate (spaceOps)
      deallocate (spaceShifts)
   end subroutine cleanSpaceGroupOperations

end module O_SpaceGroupOperations

module O_CrystalSystem

   ! Import necessary parameter modules.
   use O_Kinds

   ! Make sure no funny variables are defined.
   implicit none

   ! Begin list of module data.
   integer :: makeFull     ! Flag to make non-primitive cell 1=full;0=primitive
   integer :: numAtoms     ! Initial number of atoms before symmetry.
   integer :: numSymmAtoms ! Number of atoms after symmetry operations applied.
   integer, allocatable, dimension (:) :: atomElementID
   integer, allocatable, dimension (:) :: atomSpeciesID
   integer, allocatable, dimension (:) :: symmElementID
   integer, allocatable, dimension (:) :: symmSpeciesID
   real (kind=double), allocatable, dimension (:,:) :: atomFractABC ! This is
         ! the initial list of atoms before any symmetry operations have been
         ! applied.
   real (kind=double), allocatable, dimension (:,:) :: symmFractABC ! This is
         ! the final list of atoms after all symmtery operations have been
         ! applied and duplicate atomic positions ignored.
   real (kind=double), dimension (3,3) :: realLattice    ! [x,y,z][a,b,c]

   ! Begin list of subroutines in this module.
   contains


   ! This subroutine will read in all the atomic positions before any symmetry
   !   operations have been applied.  The coordinates are in fractional abc
   !   form.  It also reads in the atomic element ID and species ID numbers for
   !   each atom so that this information is also included for the symmetry
   !   created atoms.
   subroutine readAtomicData(fileUnit)

      ! Import necessary parameter modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define dummy variables.
      integer, intent (in) :: fileUnit

      ! Define local variables.
      integer :: i

      ! Read the number of atoms.
      read (fileUnit,*) numAtoms

      ! Allocate space to hold the element ID, species ID, and the atomic abc
      !   fractional coordinates.
      allocate (atomFractABC(3,numAtoms))
      allocate (atomElementID(numAtoms))
      allocate (atomSpeciesID(numAtoms))

      ! Read the atomic abc fractional coordinates.
      do i = 1, numAtoms
         read (fileUnit,*) atomElementID(i),atomSpeciesID(i),atomFractABC(:,i)
      enddo
   end subroutine readAtomicData


   subroutine readLattice(fileUnit)

      ! Import necessary parameter modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define dummy variables.
      integer, intent (in) :: fileUnit

      ! Read the flag that requests a full (non-primitive) cell.
      read (fileUnit,*) makeFull

      ! Read the lattice vectors.
      read(fileUnit,*) realLattice(:,:) ! (x,y,z)(a,b,c)

   end subroutine readLattice


   subroutine writeLattice(fileUnit)

      ! Import necessary parameter modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define dummy variables.
      integer, intent (in) :: fileUnit

      ! Write the lattice vectors.
      write(fileUnit,fmt="(3e20.8)") realLattice(:,1) ! (x,y,z)(a)
      write(fileUnit,fmt="(3e20.8)") realLattice(:,2) ! (x,y,z)(b)
      write(fileUnit,fmt="(3e20.8)") realLattice(:,3) ! (x,y,z)(c)

   end subroutine writeLattice


   subroutine writeAtomicData(fileUnit)

      ! Import necessary parameter modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define dummy variables.
      integer, intent (in) :: fileUnit

      ! Define local variables.
      integer :: i

      ! Print the number of atoms in the system.  If the system was reduced,
      !   this also would print the correct number because that reduced number
      !   of atoms was copied to numSymmAtoms.
      write (fileUnit,*) numSymmAtoms

      ! Print the symmetry created atoms.  Note that this will also print
      !   the atoms in the reduced (primitive) cell because when it was
      !   created the atoms were copied back to these data structures after
      !   determining which should remain.  The numSymmAtoms would also be
      !   be changed in the process if the reduction was done.
      do i = 1, numSymmAtoms
         write (fileUnit,fmt="(2i5,3e20.8)") symmElementID(i),symmSpeciesID(i),&
               & symmFractABC(:,i)
      enddo
   end subroutine writeAtomicData

   
   subroutine prepareSymmAtomicCoords

      ! Import necessary parameter modules.
      use O_Kinds

      ! Import necessary object modules.
      use O_SpaceGroupOperations

      ! Make sure no funny variables are defined.
      implicit none

      ! Allocate space to hold the maximum possible number of symmetry created
      !   atomic positions, new element IDs, and new species ID, assuming that
      !   no duplicate are created.  It is most likely that not all of these
      !   available positions will be used though.
      allocate (symmFractABC(3,numAtoms*numSpaceOps))
      allocate (symmElementID(numAtoms*numSpaceOps))
      allocate (symmSpeciesID(numAtoms*numSpaceOps))

      ! Initialize the number of symmetry created atoms.
      numSymmAtoms = 0

   end subroutine prepareSymmAtomicCoords


   subroutine applySymmetry

      ! Import necessary parameter modules.
      use O_Kinds

      ! Import necessary object modules.
      use O_SpaceGroupOperations

      ! Make sure no funny variables are defined.
      implicit none

      ! Define local variables.
      integer :: i,j,k,l
      real (kind=double), dimension (3) :: newPosition

      ! For each atom apply all symmetry operations
      do i = 1, numAtoms
         do j = 1, numSpaceOps

            ! Create the new atomic position base upon the old atomic position
            !   and the current space group operation.
            do k = 1,3 ! abc axes of the new atomic position.

               ! Initialize the new position for this axis.
               newPosition(k) = 0.0_double

               ! Add together the contributions to this k from each abc axis.
               do l = 1,3
                  newPosition(k) = newPosition(k) + atomFractABC(l,i) * &
                        & spaceOps(l,k,j)
               enddo

               ! Apply the axis shift.
               newPosition(k) = newPosition(k) + spaceShifts(k,j)
            enddo

            ! Ensure that this new atomic position is within the cell.
            call shiftToCell(newPosition)

            ! Add this new atom position only if it is not a duplicate to other
            !   previously determined atomic positions.
            if (findMatch(newPosition) == 0) then

               ! Increment the number of atoms in the system.
               numSymmAtoms = numSymmAtoms + 1

               ! Save the atomic location.
               symmFractABC(:,numSymmAtoms) = newPosition(:)

               ! Save the element ID number and Species ID number.
               symmElementID(numSymmAtoms) = atomElementID(i)
               symmSpeciesID(numSymmAtoms) = atomSpeciesID(i)
            endif
         enddo
      enddo
   end subroutine applySymmetry


   subroutine shiftToCell(newPosition)

      ! Import necessary parameter modules.
      use O_Kinds
      use O_Constants

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed parameters.
      real (kind=double), dimension (3) :: newPosition

      ! Define local variables.
      integer :: i

      ! This is actually a bit tricky.  We will demand that all atoms that are
      !   sufficiently close to a border line be made to rest *exactly* on that
      !   border line.  Further, we will make any atom resting on a border line
      !   rest on the 0 value border as opposed to the 1.0 valued border.
      do i = 1,3
         if (abs(newPosition(i)) < smallThresh) then
            newPosition(i) = 0.0_double
         elseif (abs(newPosition(i) - 1.0_double) < smallThresh) then
            newPosition(i) = 0.0_double
         endif
      enddo

      ! Then the remaining atoms that exist beyond any border line will be
      !   shifted inside the cell.
      do i = 1,3
         if (newPosition(i) < 0.0_double) then
            newPosition(i) = newPosition(i) + 1.0_double
         elseif (newPosition(i) > 1.0_double) then
            newPosition(i) = newPosition(i) - 1.0_double
         endif
      enddo

   end subroutine shiftToCell


   function findMatch(newPosition)

      ! Import necessary parameter modules.
      use O_Kinds
      use O_Constants

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed parameters.
      real (kind=double), dimension (3) :: newPosition

      ! Define local variables.
      integer :: i
      integer :: findMatch

      ! Assume that no match will be found.
      findMatch = 0

      ! Compare this new position to all previously created symmetric positions.
      do i = 1, numSymmAtoms
         if (sum((newPosition(:)-symmFractABC(:,i))**2) < smallThresh) then
            findMatch = 1
            exit
         endif
      enddo

   end function findMatch


   subroutine reduceCell (spaceLattice)

      ! Import necessary parameter modules.
      use O_Kinds
      use O_Constants

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed parameters.
      character*1 :: spaceLattice

      ! Some of the members of this set of local variables will be copied back
      !   to the "symm" set of module variables for easy printing.
      integer :: i,j,k,l,m
      integer :: numReduceAtoms
      integer, allocatable, dimension (:) :: reduceElementID
      integer, allocatable, dimension (:) :: reduceSpeciesID
      real (kind=double), dimension(3) :: shiftVector
      real (kind=double), dimension(3) :: shiftedReduceFractABC
      real (kind=double), allocatable, dimension (:,:) :: reduceFractABC
      real (kind=double), dimension (3,3) :: reduceLattice ! [x,y,z][a,b,c]
      real (kind=double), dimension (3,3) :: reduceLatticeInv

      ! In the case that a full (non-primitive) cell was requested or if the
      !   cell is already a primitive cell (and therefor can not be reduced)
      !   we have nothing to do here.  I will have to figure out later how to
      !   reduce the H centered trigonal cells.  Until then, they will have to
      !   be done as non-primitive cells only.
      if ((makeFull == 1) .or. (spaceLattice == 'P') .or. &
            & (spaceLattice == 'H')) then
         return
      endif

      ! Allocate space to hold the reduced set of atoms and atom ID numbers
      !   which could be identical to the symm set but not bigger.
      allocate (reduceFractABC(3,numSymmAtoms))
      allocate (reduceElementID(numSymmAtoms))
      allocate (reduceSpeciesID(numSymmAtoms))

      ! First we will reduce the lattice vectors.
      call reduceLatticeVectors(spaceLattice,reduceLattice,reduceLatticeInv)

      ! Initialize the number of reduced atoms.  Only the atoms that remain
      !   inside the new reduced cell will be kept.
      numReduceAtoms = 0

      ! Now we convert atomic coordinates from their representation in
      !   fractional abc units where abc are the lattice vectors of the
      !   original non-primitive cell to the primitive cell lattice vectors
      !   (in reduceLattice) called abc' (abcPrime).
      do i = 1, numSymmAtoms

         ! I know it is a bit ugly but watch carefully what is done to the
         !   reduceFractABC.  The atom position in the abc' is computed and
         !   saved in index i which is *always* >= numReduceAtoms.  The
         !   'accepted' atoms are saved in the reduceFractABC indices that are
         !   *always* = numReduceAtoms which is <= i and numSymmAtoms.  (i.e.
         !   this array is used for two purposes, computing all new positions,
         !   and saving only the ones that qualify.)  There is a ** by the
         !   tricky spot.  This was done to make the isUniqueAtom a little
         !   easier.

         ! Convert the coordinates into fractional abcPrime units.
         call getFractABCPrime(symmFractABC(:,i),realLattice,&
               & reduceFractABC(:,i),reduceLatticeInv)

         ! Now, it is the case for (at least) the conversion from the trigonal
         !   full-cell (which is hexagonal) to the associated primitive cell
         !   (which is rhomohedral) and from the body-centered-cubic full cell
         !   to the associated primitive cell that the primitive cell is *not*
         !   spacially contained within the full cell and thus simply checking
         !   to see if the full cell atoms are also inside the primitive cell
         !   will not be sufficient to "populate" the primitive cell.  It is
         !   necessary rather to check the atoms in their full cell positions
         !   *and* those same atoms after the full cell has been shifted to
         !   any of the neighboring or diagonally neighboring cells.

         do j = -1, 1
         do k = -1, 1
         do l = -1, 1

            ! Create a vector to shift the atom's position by.
            shiftVector = (/j,k,l/)

            ! In most cases we now simply shift the fractional position, but
            !   in the case that one of the coordinates is equal to either
            !   zero or one we don't apply the shift along that axis because
            !   this would simply create dumplicate atoms.
            do m = 1, 3
               !  If the coordinate is not 1, -1, or 0 then we add the shift
               !   vector to it.  Otherwise we just keep the coordiante as is.
               if (.not.((abs(abs(reduceFractABC(m,i)) - 1.0_double) < &
                     & smallThresh) .or. (abs(reduceFractABC(m,i)) < &
                     & smallThresh))) then
                  shiftedReduceFractABC(m) = reduceFractABC(m,i) + &
                        & shiftVector(m)
               else
                  shiftedReduceFractABC(m) = reduceFractABC(m,i)
               endif
            enddo

            ! Determine if this atom is inside the new primitive cell or if it
            !   will be in the primitive cell after a translation of the full
            !   cell.  If so, (and if it does not duplicate another already
            !   present atom) then we save it.  Otherwise we ignore it.
            if ((inPrimitiveCell(shiftedReduceFractABC(:)) == 1) .and. &
                  & (isUniqueAtom(numReduceAtoms,reduceFractABC(:,:),&
                  &  shiftedReduceFractABC(:)) == 1)) then
               numReduceAtoms = numReduceAtoms + 1
               reduceFractABC(:,numReduceAtoms) = shiftedReduceFractABC(:)!**
               reduceElementID(numReduceAtoms) = symmElementID(i)
               reduceSpeciesID(numReduceAtoms) = symmSpeciesID(i)
            endif
         enddo
         enddo
         enddo
         
      enddo

      ! Once the conversion is complete we will copy the results back over the
      !   original symm set for printing.  The components that are not
      !   overwritten will be ignored because numSymmAtoms will be equal to
      !   numReduceAtoms.
      do i = 1, numReduceAtoms
         symmFractABC(:,i) = reduceFractABC(:,i)
         symmElementID(i)  = reduceElementID(i)
         symmSpeciesID(i)  = reduceSpeciesID(i)
      enddo
      realLattice(:,:) = reduceLattice(:,:)
      numSymmAtoms = numReduceAtoms

      ! Deallocate memory that is no longer needed.
      deallocate (reduceFractABC)
      deallocate (reduceElementID)
      deallocate (reduceSpeciesID)
   end subroutine reduceCell


   subroutine reduceLatticeVectors (spaceLattice,reduceLattice,reduceLatticeInv)

      ! Import necessary parameter modules.
      use O_Kinds
      use O_Constants

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed parameters.
      character*1 :: spaceLattice
      real (kind=double), dimension (3,3) :: reduceLattice
      real (kind=double), dimension (3,3) :: reduceLatticeInv ! Inverse

      ! Define local variables.
      real (kind=double) :: twoThrd
      real (kind=double) :: oneThrd

      ! Define local variables for matrix inversion via lapack dgetrf,dgetri.
      integer :: info
      integer :: lwork
      integer, dimension (3) :: pivotIndices
      real (kind=double), allocatable, dimension (:) :: work

      ! Define the interface to the LAPACK subroutines.
      interface
         integer function ilaenv (ISPEC,NAME,OPTS,N1,N2,N3,N4)
            use O_Kinds
            integer       :: ISPEC
            character*(*) :: NAME
            character*(*) :: OPTS
            integer       :: N1
            integer       :: N2
            integer       :: N3
            integer       :: N4
         end function ilaenv
         subroutine dgetrf (M,N,A,LDA,IPIV,INFO)
            use O_Kinds
            integer :: INFO
            integer :: LDA
            integer :: M
            integer :: N
            integer, dimension(min(M,N)) :: IPIV
            real (kind=double), dimension (LDA,N) :: A
         end subroutine dgetrf
         subroutine dgetri (N,A,LDA,IPIV,WORK,LWORK,INFO)
            use O_Kinds
            integer :: INFO
            integer :: LDA
            integer :: LWORK
            integer :: N
            integer, dimension (N) :: IPIV
            real (kind=double), dimension (LDA,N) :: A
            real (kind=double), dimension (LWORK) :: WORK
         end subroutine dgetri
      end interface

      ! Initialize constants.
      twoThrd = 2.0_double / 3.0_double
      oneThrd = 1.0_double / 3.0_double

      ! Initialize variables for lapack inversion.
      lwork = 3*ilaenv(1,'dgetri','',3,3,-1,-1)
      allocate (work(lwork))

      ! Compute the new lattice parameters based on the spaceLattice and the
      !   old lattice parameters.  This information is from "Symmetry
      !   Principles in Solid State and Molecular Physics" by Melvin Lax,
      !   Dover, New York, 2001.
      if (spaceLattice == 'A') then
         reduceLattice(:,1) = realLattice(:,1)
         reduceLattice(:,2) = realLattice(:,2)
         reduceLattice(:,3) = 0.5_double*(realLattice(:,2)+realLattice(:,3))
      elseif (spaceLattice == 'B') then
         reduceLattice(:,1) = realLattice(:,1)
         reduceLattice(:,2) = realLattice(:,2)
         reduceLattice(:,3) = 0.5_double*(realLattice(:,1)+realLattice(:,3))
      elseif (spaceLattice == 'C') then
         reduceLattice(:,1) = realLattice(:,1)
         reduceLattice(:,2) = 0.5_double*(realLattice(:,1)+realLattice(:,2))
         reduceLattice(:,3) = realLattice(:,3)
      elseif (spaceLattice == 'F') then
         reduceLattice(:,1) = 0.5_double*(realLattice(:,2)+realLattice(:,3))
         reduceLattice(:,2) = 0.5_double*(realLattice(:,1)+realLattice(:,3))
         reduceLattice(:,3) = 0.5_double*(realLattice(:,1)+realLattice(:,2))
      elseif (spaceLattice == 'I') then
         reduceLattice(:,1) =(-0.5_double)*realLattice(:,1) + &
                            &  0.5_double*realLattice(:,2) + &
                            &  0.5_double*realLattice(:,3)
         reduceLattice(:,2) =  0.5_double*realLattice(:,1) + &
                            &(-0.5_double)*realLattice(:,2) + &
                            &  0.5_double*realLattice(:,3)
         reduceLattice(:,3) =  0.5_double*realLattice(:,1) + &
                            &  0.5_double*realLattice(:,2) + &
                            &(-0.5_double)*realLattice(:,3)
      elseif (spaceLattice == 'R') then
         reduceLattice(:,1) =  twoThrd*realLattice(:,1) + &
                            &  oneThrd*realLattice(:,2) + &
                            &  oneThrd*realLattice(:,3)
         reduceLattice(:,2) =(-oneThrd)*realLattice(:,1) + &
                            &  oneThrd*realLattice(:,2) + &
                            &  oneThrd*realLattice(:,3)
         reduceLattice(:,3) =(-oneThrd)*realLattice(:,1) + &
                            &(-twoThrd)*realLattice(:,2) + &
                            &  oneThrd*realLattice(:,3)
      endif

      ! Compute the inverse of the reduced lattice.
      reduceLatticeInv(:,:) = reduceLattice(:,:)
      call dgetrf(3,3,reduceLatticeInv,3,pivotIndices,info)
      if (info /= 0) then
         write (6,*) "Not able to perform inversion of lattice part 1: info=",&
               & info
         stop
      endif
      call dgetri(3,reduceLatticeInv,3,pivotIndices,work,lwork,info)
      if (info /= 0) then
         write (6,*) "Not able to perform inversion of lattice part 2: info=",&
               & info
         stop
      endif

      ! Deallocate data that is no longer needed.
      deallocate (work)

   end subroutine reduceLatticeVectors


   subroutine getFractABCPrime(symmFractABC,symmLattice,&
               & reduceFractABC,reduceLatticeInv)

      ! Import necessary parameter modules.
      use O_Kinds
      use O_Constants

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed dummy parameters.
      real (kind=double), dimension (3)   :: symmFractABC
      real (kind=double), dimension (3,3) :: symmLattice
      real (kind=double), dimension (3)   :: reduceFractABC
      real (kind=double), dimension (3,3) :: reduceLatticeInv

      ! Define local variables.
      integer :: i
      real (kind=double), dimension(3) :: directXYZ

      ! Pxyz  = PxX   + PyY   + PzZ
      ! Pabc  = PaA   + PbB   + PcC
      ! Pabc' = Pa'A' + Pb'B' + Pc'C'

      ! L    = symmLattice
      ! L'   = reduceLattice
      ! L'-1 = reduceLatticeInv

      ! Pxyz  = Pabc  L
      ! Pxyz  = Pabc' L'  ->  Pxyz L'-1 = Pabc' L'L'-1  ->  Pxyz L'-1 = Pabc'
      ! Pabc' = Pxyz L'-1

      ! Convert the fractional abc coordinates from the symmetry-applied set of
      !   atoms to the primitive fractional abc coordinates.  This is done by
      !   converting to xyz cartesian coordinates and then to the reduced abc
      !   coordinate system.  Note the order of the parameters to the matmul
      !   function call.  This is important.  Because fortran is column major
      !   in the way that it accesses matrices we have to order the parameters
      !   this way or else we would have to transpose the matrices before doing
      !   the multiplication.

      directXYZ(:)      = matmul(symmLattice(:,:),symmFractABC(:))
      reduceFractABC(:) = matmul(reduceLatticeInv(:,:),directXYZ(:))

      ! Demand that the coordinates are strictly less than 1.0 to help ensure
      !   that we don't include the same atom twice (e.g. 0.000 0.000 0.000 and
      !   1.000 1.000 1.000).
      do i = 1, 3
         if (abs(reduceFractABC(i) - 1.0_double) < smallThresh) then
            reduceFractABC(i) = reduceFractABC(i) - 1.0_double
         endif
      enddo

   end subroutine getFractABCPrime


   function inPrimitiveCell (fractPosition)

      ! Import necessary parameter modules.
      use O_Kinds
      use O_Constants

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed dummy parameters.
      real (kind=double), dimension (3) :: fractPosition

      ! Define local variables.
      integer :: inPrimitiveCell
      integer :: i

      ! Assume that the atomic position is within the primitive cell.
      inPrimitiveCell = 1

      ! If the fractional position is somewhere between 1+epsilon and 0-epsilon
      !   then the atom is considered to be within the cell for that coordinate
      !   axis.  If it falls outside that range, then it isn't inside the cell.
      do i = 1, 3
         if ((fractPosition(i) >= 1.0_double+smallThresh) .or. &
            &(fractPosition(i) < 0.0_double-smallThresh)) then
            inPrimitiveCell = 0
            exit
         endif
      enddo

   end function inPrimitiveCell


   function isUniqueAtom (numReduceAtom,fractPositions,testPosition)

      ! Import necessary parameter modules.
      use O_Kinds
      use O_Constants

      ! Make sure no funny variables are defined.
      implicit none

      ! Define passed dummy parameters.
      integer :: numReduceAtom
      real (kind=double), dimension (3,numReduceAtom) :: fractPositions
      real (kind=double), dimension (3) :: testPosition

      ! Define local variables.
      integer :: isUniqueAtom
      integer :: i,j

      ! Assume that the atomic position stored in
      !   fractPositions(:,numReduceAtom+1) is a unique atom in the cell.
      isUniqueAtom = 1

      ! Compare the current atom to each other atom one at a time.  If the
      !   current atom and the other atom differ in any coordinate by more than
      !   smallThresh then the atoms are considered to be different.  If they
      !   are no more than smallThresh different for each coordinate, then they
      !   are considered to be the same and we will abort the rest of the
      !   search and return a value of 0 to indicate that the given atom
      !   matches a previously listed one.
      do i = 1, numReduceAtom
         do j = 1, 3

            ! Compare a coordinate.  If it is the same, then the atom may be a
            !   duplicate.  If it is different then we can move on to the next
            !   atom.
            if (abs(testPosition(j)-fractPositions(j,i)) < 100*smallThresh) then
               isUniqueAtom = 0
            else
               isUniqueAtom = 1
               exit
            endif
         enddo

         ! The isUniqueAtom value can only be 0 here if there was a match for
         !   all three coordinates.  In such a case the current atom is a
         !   duplicate and so we abort the rest of the search and return a 0.
         if (isUniqueAtom == 0) then
            exit
         endif
      enddo

   end function isUniqueAtom


   subroutine cleanCrystalSystem

      ! Make sure no funny variables are defined.
      implicit none

      ! Deallocate memory blocks.
      deallocate (atomFractABC)
      deallocate (atomElementID)
      deallocate (atomSpeciesID)
      deallocate (symmFractABC)
      deallocate (symmElementID)
      deallocate (symmSpeciesID)
      
   end subroutine cleanCrystalSystem

end module O_CrystalSystem

program applySpaceGroup

   ! This program is used to apply a set of symmetry operations onto a set of
   !   fractional atomic positions.

   ! All the data are read from standard input and written to standard output.

   ! Import necessary parameter modules.
   use O_Kinds
   use O_Constants

   ! Import necessary object modules.
   use O_SpaceGroupOperations
   use O_CrystalSystem

   ! Make sure no funny variables are defined.
   implicit none

   open(unit=5,file='sginput',form='formatted',status='unknown')
   open(unit=6,file='sgoutput',form='formatted',status='unknown')

   ! Read the symmetry operations and create the necessary data structures to
   !   hold and access that information.
   call readSpaceOps(5) ! Read from standard input.

   ! Read lattice parameters, and request for full cell.
   call readLattice(5) ! Read from standard input.

   ! Read the element data and atomic coordinates in fractional a,b,c units.
   call readAtomicData(5) ! Read from standard input.

   ! Once the space group data and atomic coordinates have been read we can
   !   prepare the data structure for holding the symmetry atoms.
   call prepareSymmAtomicCoords

   ! Apply the symmetry operations to the given atomic positions ignoring
   !   duplicates.
   call applySymmetry

   ! Attempt to make a primitive cell.
   call reduceCell(spaceLattice)

   ! Write the new lattice parameters.
   call writeLattice(6) ! Write to standard output.

   ! Write the element data and symmetry created atomic coordinates in
   !   fractional a,b,c units.
   call writeAtomicData(6) ! Write to standard output.

   ! Clean up allocated memory.
   call cleanCrystalSystem
   call cleanSpaceGroupOperations

end program applySpaceGroup
