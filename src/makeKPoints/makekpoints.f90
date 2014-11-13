!#include "../config.h"
module PointGroupOperations_O

   ! Import necessary modules
   use O_Kinds

   ! Make sure no funny variables are defined.
   implicit none

   ! Begin list of module data.
   integer :: numPointOps
   real (kind=double), allocatable, dimension (:,:)   :: translations
   real (kind=double), allocatable, dimension (:,:,:) :: pointOps
   real (kind=double), allocatable, dimension (:,:,:) :: abcRecipPointOps


   ! Begin list of module functions.
   contains


   ! This subroutine will read data from the given unit.  It will initialize
   !   the number of point group operations and the values of the point group
   !   operations.  In the process it will allocate pointOps.  The point group
   !   operations are derived from the space group operations data base.  They
   !   are simply the same operations without any translational component.
   subroutine readPointOps (fileUnit)

      ! Make sure no funny variables are defined.
      implicit none

      ! Define dummy variables.
      integer, intent (in) :: fileUnit

      ! Define local variables.
      integer :: i,j
      integer :: spaceGroupNum
      integer :: numSpaceOps
      integer :: numShifts

      ! Read past the documentation line that defines the space group.
      read (fileUnit,*)

      ! Read the root space group number for this space group.
      read (fileUnit,*) spaceGroupNum

      ! Read the number of space group operations from the file unit and the
      !   the number of shifted repetitions that exist in the list of
      !   operations.  The number of point group operations will be equal to
      !   the number of space group operations divided by the number of shifted
      !   sets.
      read (fileUnit,*) numSpaceOps, numShifts

      ! Compute the number of point group operations.
      numPointOps = numSpaceOps / numShifts

      ! Initialize storage space for the point group operations and any
      !   included translations.
      allocate (pointOps(3,3,numPointOps))
      allocate (translations(3,numPointOps))
      pointOps(:,:,:) = 0.0_double
      translations(:,:) = 0.0_double

      ! Read the cartesian form of the space group operations from the file
      !   unit.  The first row indicates the contributions to the new a
      !   coordinate, the second for b and the third row for c.  The fourth
      !   contains the translations to be applied to each a,b,c axis.
      do i = 1, numPointOps

         ! Read the empty line between operations.
         read (fileUnit,*)

         ! Read the space group operation.
         read (fileUnit,*) pointOps(:,:,i) ![abc contrib][new abc]

         ! Read the translation line.
         read (fileUnit,*) translations(:,i) ![new abc]

!         ! Make sure that this point group operation does not match any
!         !   previously read in point group operations.  If it does not match
!         !   then we can save it.  Matches might happen if two space group
!         !   operations had the same xyz contribution components for each new
!         !   xyz value but different translation effects.  1=duplicate;0=not
!         if (samePointOps(tempPointOp) == 0) then
!            numPointOps = numPointOps + 1
!            pointOps(:,:,numPointOps) = tempPointOp(:,:)
!         endif

      enddo ! i=1,numPointOps

      ! Read past all the space group operations that are just point group
      !   operations plus extra translation.
      do i = 1, numSpaceOps - numPointOps
         do j = 1,5
            read (fileUnit,*)
         enddo
      enddo

   end subroutine readPointOps


   function samePointOps (tempPointOp)

      ! Import necessary parameter modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define dummy parameters.
      real (kind=double), intent (in), dimension (3,3) :: tempPointOp

      ! Define local variables.
      integer :: i,j,k
      integer :: notAMatch
      integer :: samePointOps

      ! This is done through two assumptions that work opposite to each other.

      ! Assume that this temp point group operation doesn't match previous ones.
      samePointOps = 0

      do i = 1, numPointOps

         ! Assume that the tempPointOp DOES match the current pointOp.
         notAMatch=0

         do j = 1,3
            do k = 1,3
               if (tempPointOp(k,j) /= pointOps(k,j,i)) then
                  ! This operation is not a match.
                  notAMatch = 1
                  exit
               endif
            enddo
            if (notAMatch == 1) exit
         enddo

         ! If it still DOES match, then we found a duplicate and we report it.
         if (notAMatch == 0) then
            samePointOps = 1
            exit
         endif
      enddo

   end function samePointOps


   ! This subroutine will take as input a real space lattice and a reciprocal
   !   space lattice.  It will use them to compute the point group operations
   !   for the non-orthogonal reciprocal space lattice.  In the process it will
   !   allocate space for abcRecipPointOps.  It uses the abc pointOps that
   !   were previously read in.
   subroutine computeABCRecipPointOps (realLattice,recipLattice)

      ! Import necessary parameter modules.
      use O_Kinds
      use O_Constants

      ! Make sure no funny variables are defined.
      implicit none

      ! Define dummy parameters.
      real (kind=double), intent (in), dimension (3,3) :: realLattice
      real (kind=double), intent (in), dimension (3,3) :: recipLattice

      ! Define local variables.
      integer :: i,j,k,l
      real (kind=double), dimension (3,3) :: tempPointOp

      ! Initialize storage space for the abcRecipPointOps.
      allocate (abcRecipPointOps(3,3,numPointOps))

      do i = 1, numPointOps

         ! Initialize the accumulator for this point group operation.
         tempPointOp(:,:) = 0.0_double

         do j = 1,3       ! Loop over x,y,z components of the real  lattice.
            do k = 1,3    ! Loop over x,y,z components of the recip lattice.
               do l = 1,3 ! Loop over a,b,c components of the recip lattice.

                  ! Array-loop over real a,b,c.
                  tempPointOp(:,l) = tempPointOp(:,l) + realLattice(j,:) * &
                        & pointOps(j,k,i) * recipLattice(k,l) / &
                        & (2.0_double * pi)
               enddo  ! l=1,3 (loop over recip a,b,c)
            enddo  ! k=1,3 (recip,x,y,z)
         enddo  ! j=1,3 (real x,y,z)

         ! Save the point group operation.
         abcRecipPointOps(:,:,i) = tempPointOp(:,:)

      enddo  ! i=1,numPointOps
   end subroutine computeABCRecipPointOps

end module PointGroupOperations_O


module Lattice_O

   ! Import necessary modules
   use O_Kinds
   use O_Constants

   ! Make sure no funny variables are defined.
   implicit none

   ! Define the module

   ! Begin list of module data.
   real (kind=double) :: realVolume
   real (kind=double) :: recipVolume
   real (kind=double), dimension (3,3) :: realLattice  ! (xyz,abc)
   real (kind=double), dimension (3,3) :: recipLattice ! (xyz,abc)


   ! Begin list of module functions.
   contains


   ! Read real lattice information from the given file unit.
   subroutine readRealLattice (fileUnit)

      ! Define dummy parameter.
      integer, intent (in) :: fileUnit

      read (fileUnit,*) realLattice(:,:)
   end subroutine readRealLattice


   ! Compute the reciprocal cell vectors and the real and reciprocal lattice
   !   cell volumes.
   subroutine computeLatticeData

      ! Make sure no funny variables are defined.
      implicit none

      ! Define local variables.
      integer :: i,j
      ! The following cycle variables are numbers from 1 to 3.  They are
      !    greater than the loop variable (i or j) by the trailing number, and
      !    they cycle back to 1 when they are greater than the the loop index
      !    m3.  (e.g. i = 1; iCycle1 = 2)  (e.g. i = 3; iCycle1 = 1).
      integer :: iCycle1, iCycle2
      integer :: jCycle1, jCycle2

      ! Given the real space primitive lattice vectors (say a b c) each defined
      !   in terms of the orthogonal coordinate vectors (x y z) (a.k.a i k j) we
      !   can find the reciprocal primitive lattice vectors (say a' b' c') by:
      !   a' = 2*Pi * (b x c)/(a . b x c)
      !   b' = 2*Pi * (c x a)/(a . b x c)
      !   c' = 2*Pi * (a x b)/(a . b x c)
      !   where . = dot product, x = cross product

      do i = 1, 3
         realVolume = 0
         iCycle1 = mod(i,3) + 1
         iCycle2 = mod(i+1,3) + 1
         do j = 1, 3
            jCycle1 = mod(j,3) + 1
            jCycle2 = mod(j+1,3) + 1
            recipLattice(j,i) = &
                  & (realLattice(jCycle1,iCycle1)  * &
                  &  realLattice(jCycle2,iCycle2)) - &
                  & (realLattice(jCycle2,iCycle1)  * &
                  &  realLattice(jCycle1,iCycle2))
            realVolume = realVolume + recipLattice(j,i) * realLattice(j,i)
         enddo

         recipLattice(:,i) = 2.0_double * Pi * recipLattice(:,i) / realVolume
      enddo

      realVolume = abs(realVolume) ! Volume is always positive.
      recipVolume = ((2.0_double*Pi)**3)/realVolume

   end subroutine computeLatticeData

end module Lattice_O



module KPointMesh_O

   ! Import necessary paramter modules.
   use O_Kinds

   ! Make sure no funny variables are defined.
   implicit none

   ! Begin list of module data.
   integer :: doGamma  ! 1=Make 1 gamma kpoint; 0=Make 1 general kpoint.
   integer :: numMeshKPoints
   integer :: numFoldedKPoints
   integer, dimension(3) :: numABCKPoints
   integer, allocatable, dimension(:) :: kPointTracker
   real (kind=double) :: kpThresh    ! Threshhold for which two kpoints are
         ! considered to be symmetricly the same via point group operations.
   real (kind=double) :: weightSum   ! Sum of initial weights of all kpoints.
         ! weightSum = 1 for relativistic; 2 for non-relativistic
   real (kind=double) :: initWeight  ! Initial weight of each mesh kpoint.
   real (kind=double), dimension(3) :: abcShift
   real (kind=double), dimension(3) :: abcDelta
   real (kind=double), allocatable, dimension(:)   :: kPointWeight
   real (kind=double), allocatable, dimension(:,:) :: abcMeshKPoints
   real (kind=double), allocatable, dimension(:,:) :: abcFoldedKPoints


   ! Begin list of module functions.
   contains


   ! Read the number of kpoints in each a,b,c direction and the shift applied
   !   to the initial mesh from the given file unit.
   subroutine readMeshParameters(fileUnit)

      ! Define dummy parameters.
      integer, intent (in) :: fileUnit

      ! Read the data.
      read(fileUnit,*) numABCKPoints(:)
      read(fileUnit,*) abcShift(:)
      read(fileUnit,*) doGamma
   end subroutine readMeshParameters


 
   ! Initialize the kpoint mesh according to the given data.  Allocate space to
   !   hold the reciprocal abc mesh points and track them in the process.
   subroutine initMesh

      ! Import necessary paramter modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define the local parameters.
      integer :: i,j,k
      integer, dimension (3) :: loopIndex

      ! Define the threshhold for when two kpoints are symmetry reduced via
      !   point group operations.  (i.e. how close do two kpoints have to be
      !   when one is folded near the other for them to be considered the same?)
      kpThresh = 0.00001_double  ! A bit arbitrary!

      ! Initialize the weight to 2 always.  (There are 2 electrons per state
      !   by default.  If a spin-polarized calculation is done, then the 1 e-
      !   per state is accounted for in the olcao fortran program.)
      weightSum = 2.0_double

      ! Compute the number of uniform mesh points in order to allocate space to
      !   hold the kpoints.
      numMeshKPoints = product(numABCKPoints(:))

      ! Compute the mesh points separation step size.
      abcDelta(:) = 1.0_double / numABCKPoints(:)

      ! Allocate necessary space.
      allocate (kPointTracker(numMeshKPoints))
      allocate (abcMeshKPoints(3,numMeshKPoints))

      ! Loop over abc kpoints numbers to create the initial uniform mesh in
      !   fractional units of the abc reciprocal space cell.
      numMeshKPoints = 0
      do i = 1, numABCKpoints(1)

         ! Store the i loop index.
         loopIndex(1) = i

         do j = 1, numABCKpoints(2)

            ! Store the j loop index.
            loopIndex(2) = j

            do k = 1, numABCKpoints(3)

               ! Increment the number of uniform mesh kpoints.
               numMeshKPoints = numMeshKPoints + 1

               ! Store the k loop index.
               loopIndex(3) = k

               ! Compute the current mesh kpoint location.
               abcMeshKPoints(:,numMeshKPoints) = &
                     & (loopIndex(:)-1.0_double+abcShift(:)) * abcDelta(:)

            enddo  ! k=1,numABCKpoints(3)
         enddo  ! j=1,numABCKpoints(2)
      enddo  ! i=1,numABCKpoints(1)


      ! Initialize the kpoint tracker.  When a new irreducable kpoint is found
      !   the value is set the negative of the number of irriducable kpoints
      !   found so far.  Then, then kpoint tracker value for all other kpoints
      !   that reduce to this one will get the same value (-# of irriducable
      !   kpoints found so far).
      do i = 1, numMeshKPoints
         kPointTracker(i) = i
      enddo

   end subroutine initMesh


   ! This subroutine will reduce the uniform mesh into a folded mesh within
   !   the first Brillouin zone according to a set of given point group
   !   operations.  The points on the uniform mesh that are equivalent will
   !   be combined together and their weight added.  The weightSum is the
   !   amount that all the kpoint weights must add up to.  In non-relativistic
   !   and non-spin polarized calculations each state has 2 electrons.  The
   !   weightSum is then 2.  In relativistic calculations every state is
   !   explicitly accounted for and so each state has 1 electron.  However,
   !   since we don't want to change the kpoint input file just for a spin
   !   polarized calculation we leave the weight sum as 2 all the time and the
   !   effect is accounted for in the olcao fortran program via the 'spin'
   !   variable that is 1 for non-spinpol and 2 for spinpol calculations.
   subroutine foldMesh(numPointOps,abcRecipPointOps)

      ! numPointOps = # of point group operations.
      ! abcRecipPointOps = Defined set of point group operations in recip space.

      ! Import necessary paramter modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define dummy parameters
      integer, intent (in) :: numPointOps
      real (kind=double), intent (in), dimension(3,3,numPointOps) :: &
            & abcRecipPointOps

      ! Define the local variables.
      integer :: i,j
      real (kind=double), dimension(3) :: foldedKPoint

      ! Allocate space to hold the folded kpoint positions and folded kpoint
      !   weights.  Note that we allocate the same amount of space that is used
      !   to hold the mesh kpoints because we do not yet know how many kpoints
      !   will be reduced and so this is much faster.  The memory usage will
      !   not be significant.
      allocate (abcFoldedKPoints(3,numMeshKPoints))
      allocate (kPointWeight(numMeshKPoints))

      ! Compute the weight assigned to each of the uniform mesh kpoints.
      initWeight = weightSum/real(numMeshKPoints,double)

      ! In the case that only 1 kpoint was requested it must be a general
      !   kpoint and so we explicitly define it.
      if (numMeshKPoints == 1) then

         numFoldedKPoints = 1
         kPointWeight(1)  = weightSum

         if (doGamma == 0) then
            ! Arbitrary general kpoint
            abcFoldedKPoints(1,1) = 0.125_double
            abcFoldedKPoints(2,1) = 0.25_double
            abcFoldedKPoints(3,1) = 1.0_double / 3.0_double
         else
            ! Gamma kpoint
            abcFoldedKPoints(:,1) = 0.0_double
         endif
         return
      endif

      ! Initialize the count for the number of kpoints that exist after the
      !   point group operations have been applied to the uniform mesh.
      numFoldedKPoints = 0

      ! Consider each mesh kpoints in turn and then see which other kpoints can
      !   be reduced to it.
      do i = 1, numMeshKPoints

         ! Has this kpoint not already been reduced to another one?  If it was
         !   already found to be the same as another kpoint due to symmetry
         !   then the value would be -1 * the number of folded kpoints found
         !   so far.
         if (kPointTracker(i) == i) then

            ! We have found a new irreducible point so we mark it with a
            !   negative value in the kpoint tracker and store its value as one
            !   of the reduced kpoints.
            numFoldedKPoints = numFoldedKPoints + 1
            kPointTracker(i) = -numFoldedKPoints
            kPointWeight(numFoldedKPoints) = initWeight

            ! Save this kpoint position as one of the irreducable kpoints.
            abcFoldedKPoints(:,numFoldedKPoints) = abcMeshKPoints(:,i)

            ! In the case that this was the last kpoint from the original
            !   unreduced mesh we will not have to perform any checks for other
            !   symmetric points.  The search for symmetry matching kpoints
            !   begins with the current point and proceeds forward to
            !   numMeshKPoints in the compareKPoints subroutine.  When a match
            !   is found then the current kpoint weight is increased and the
            !   higher index number kpoint is marked negative.
            if (i == numMeshKPoints) cycle

            ! Now we will operate on the irreducible kpoint with the symmetry
            !  operations of this system's point group.  The idea is to find
            !  other kpoints that are symmetric to this one.
            do j = 1, numPointOps

               ! Fold the kpoint for this point group operation.
               call foldKPoint (foldedKPoint,abcRecipPointOps(:,:,j),&
                     & abcMeshKPoints(:,i))

               ! Find and mark the kpoints that are related to the current
               !   irreducable kpoint according to the folding that was just
               !   applied.
               call compareKPoints (foldedKPoint,i)

            enddo  ! j=1, numPointOps

         endif  ! (kPointTracker(i) == i)
      enddo  ! i = 1, numMeshKPoints
   end subroutine foldMesh


   ! This subroutine will apply a given point group operation to a given kpoint
   !   and will return a folded kpoint.
   subroutine foldKPoint (foldedKPoint,abcRecipPointOp,abcMeshKPoint)

      ! Import necessary paramter modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define dummy parameters.  NOTE that the second parameter is singular
      !   "abcRecipPointOp" and not "abcRecipPointOps".
      real (kind=double), intent (inout), dimension (3)   :: foldedKPoint
      real (kind=double), intent (in),    dimension (3,3) :: abcRecipPointOp
      real (kind=double), intent (in),    dimension (3)   :: abcMeshKPoint

      ! Define local variables.
      integer :: i

      do i = 1, 3
         foldedKPoint(i) = sum(abcRecipPointOp(i,:) * abcMeshKPoint(:))

         ! Translate the folded kpoint to the interval 0-1.
         if (foldedKPoint(i) < 0.0_double) then
            foldedKPoint(i) = foldedKPoint(i) + 1.0_double
         elseif (foldedKPoint(i) > 1.0_double) then
            foldedKPoint(i) = foldedKPoint(i) - 1.0_double
         endif
      enddo  ! i=1,3

   end subroutine foldKPoint


   ! This subroutine is used to compare which of the remaining unfolded kpoints
   !   are equivalent to the given folded kpoint.
   subroutine compareKPoints (foldedKPoint,i)

      ! Import necessary paramter modules.
      use O_Kinds

      ! Make sure no funny variables are defined.
      implicit none

      ! Define the dummy parameters
      real (kind=double), intent (in), dimension(3) :: foldedKPoint
      integer, intent (in) :: i

      ! Define local variables.
      integer :: j,k
      integer :: isMatch

      ! Search through the remaining kpoints.
      do j = i+1,numMeshKPoints

         ! Consider only those kpoints that have not yet been matched to an
         !   irriducable kpoint.
         if (kPointTracker(j) == j) then

!            ! There are three separate ways that the two kpoints can be
!            !   considered to match.

            !  (1) Is the vector difference negligable?
            isMatch = 1  ! Assume a match.
            do k = 1,3
               ! If the value in any one direction is sufficiently different
               !   then we don't have a match.
               if (abs(foldedKPoint(k) - abcMeshKPoints(k,j)) > kpThresh) then
                  isMatch = 0
                  exit
               endif
            enddo

            if (isMatch == 1) then
               call saveKPoint(j)
               cycle
            endif


!            ! (2) Is the shifted vector difference negligable?
!            isMatch = 1  ! Assume a match.
!            do k = 1,3
!               ! If the value in any one direction is sufficiently different
!               !   including the lattice shift then we don't have a match.
!               if (abs(1.0_double - foldedKPoint(k) - abcMeshKPoints(k,j)) > &
!                     & kpThresh) then
!                  isMatch = 0
!                  exit
!               endif
!            enddo

!            if (isMatch == 1) then
!write (6,*) "foldedKP=",foldedKPoint(:)
!write (6,*) "abcMeshKP=",abcMeshKPoints(:,j)
!               call saveKPoint(j)
!write (6,*) "Matched #2"
!               cycle
!            endif


!            ! (3) Is the vector sum negligable?
!            isMatch = 1  ! Assume a match.
!            do k = 1,3
!               if (abs(foldedKPoint(k) + abcMeshKPoints(k,j)) > kpThresh) then
!                  isMatch = 0
!                  exit
!               endif
!            enddo

!            if (isMatch == 1) then
!               call saveKPoint(j)
!write (6,*) "Matched #3"
!               cycle
!            endif
         endif
      enddo
   end subroutine compareKPoints


   subroutine saveKPoint (j)

      implicit none

      ! Define dummy parameters.
      integer, intent (in) :: j

      kPointWeight(numFoldedKPoints) = kPointWeight(numFoldedKPoints) + &
            & initWeight
      kPointTracker(j) = -numFoldedKPoints
   end subroutine saveKPoint


   ! This subroutine will print the folded abc kpoint information in a format
   !   suitable for reading by the olcao program.
   subroutine printKPoints (fileUnit)

      implicit none

      ! Define dummy parameters.
      integer, intent (in) :: fileUnit

      ! Define local variables.
      integer :: i

      write (fileUnit,fmt="(a17)") "NUM_BLOCH_VECTORS"
      write (fileUnit,fmt="(i9.9)") numFoldedKPoints

      write (fileUnit,fmt="(a19)") "NUM_WEIGHT_KA_KB_KC"
      do i = 1, numFoldedKPoints
         write (fileUnit,fmt="(i5,2x,f20.16,3(2x,f12.8))") i,kPointWeight(i),&
               & abcFoldedKPoints(:,i)
      enddo
   end subroutine printKPoints

end module KPointMesh_O



program makekpoints

   ! This program is used to generate a set of kpoints for a given crystal and
   !   kpoint mesh request.

   ! Import necessary parameter modules.
   use O_Kinds
   use O_Constants

   ! Import necessary object modules.
   use PointGroupOperations_O
   use Lattice_O
   use KPointMesh_O


   ! Make sure no funny variables are defined.
   implicit none

   open (unit=50,file='kpSpecs.dat',form='formatted')
   open (unit=51,file='kpSpecs.out',form='formatted')

   ! Read the lattice and compute all of its information.
   call readRealLattice(50)
   call computeLatticeData

   ! Read the point group operations and compute their values on the given
   ! a,b,c lattice.
   call readPointOps(50)
   call computeABCRecipPointOps(realLattice,recipLattice)

   ! Read the kpoint mesh parameters, initialize the mesh, fold the mesh, and
   !   print the results.  It should be generally understood that the variables
   !   with abc in the name refer to the reciprocal lattice abc vectors and not
   !   the real space lattice.  This is only given explictly for the variable
   !   abcRecipPointOps though.
   call readMeshParameters(50)
   call initMesh
   call foldMesh (numPointOps,abcRecipPointOps)
   call printKPoints(51)

   close (50)
   close (51)

end program makekpoints
