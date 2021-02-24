!#include "../config.h"

!module O_LAPACKDGESV
!   interface
!      subroutine dgesv(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
!         integer :: N
!         integer :: NRHS
!         real (kind=double), dimension(LDA,N) :: A
!         integer :: LDA
!         integer, dimension(N) :: IPIV
!         real (kind=double), dimension(LDB,N) :: B
!         integer :: LDB
!         integer :: INFO
!      end subroutine dgesv
!   end interface
!
!   contains
!
!subroutine solveDGESV (N, NRHS, A, B, INFO)
!end module O_LAPACKDGESV

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


   ! Begin list of module subroutines.
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

      ! Define return value.
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
   integer, allocatable, dimension (:) :: numBZVertices
   integer, allocatable, dimension (:) :: numBZEdges
   integer, allocatable, dimension (:) :: numBZFacets
   real (kind=double), dimension (4,26) :: planeEquation
   real (kind=double), dimension (3,26) :: latticePointXYZ

   type vertexType
      real (kind=double), dimension (3) :: coord
      integer, dimension (3) :: planeTripleIndex
   end type vertexType
   type (vertexType), dimension (100000) :: fullVertexList
   type (vertexType), dimension (100) :: uniqueVertexList
   integer :: numFullVertices

   type edgeType
      type (vertexType), dimension (2) :: vertex
      integer, dimension (2) :: sharedPlaneIndex ! The two vertices are made
            ! from the same two planes. Here, we hold the plane index numbers
            ! of those two planes.
      integer, dimension (2,2) :: sharedPlaneTripleIndexIndices ! Each of the
            ! plane index numbers is also held in the vertex as one of its
            ! planeTripleIndex numbers. Here, we hold the array index number
            ! (1, 2, or 3) that links the sharePlaneIndex number to their
            ! positions inside the vertex%planeTripleIndex array. (1,:) links
            ! sharedPlaneIndex(1) to vertices 1 and 2. (2,:) links the
            ! sharedPlaneIndex(2) to vertices 1 and 2.
   end type edgeType
   type (edgeType), dimension (1000) :: fullEdgeList
   type (edgeType), dimension (1000) :: uniqueEdgeList
   integer :: numFullEdges

   type facetType
      integer :: numVertices
      type (vertexType), dimension (20) :: vertex
   end type facetType
   type (facetType), dimension (50) :: facetList
   

   ! Begin list of module subroutines.
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


   subroutine computeBrillouinZones (maxBZ)

      ! Make sure no funny variables are defined.
      implicit none

      ! Define dummy variables.
      integer :: maxBZ ! Highest BZ we will compute.

      ! Allocate space to hold BZ data.
      allocate (numBZVertices(maxBZ))
      allocate (numBZEdges(maxBZ))
      allocate (numBZFacets(maxBZ))

      ! FIX At present this will only work for the first BZ. (I.e., maxBZ==1)

      ! Construct a list of all facets involved.
      call makeFacetList(maxBZ)

!      ! Construct a list of all possible vertices. A vertex is created by the
!      !   the intersection of three non-parallel planes. Each plane is defined
!      !   by being 1/2 of a reciprocal lattice vector distant from the origin
!      !   of the reciprocal lattice. Each vertex will be recorded as a set of
!      !   coordinates in reciprocal space along with the set of which three
!      !   planes were used to create the vertex. An additional list of the
!      !   *unique* vertices will be stored.
!      write (52,*) "About to makeFullVertexList"
!      call flush (52)
!      call makeFullVertexList
!
!      ! For each plane, produce a record of the vertices that this plane
!      !   participates in. If it participates in a minimum of three, then it
!      !   defines a facet of the Brillouin zone and it is stored as a list of
!      !   the vertices. The vertices for this facet are stored in a spatially
!      !   ring-sorted order. (I.e., one vertex is selected arbitrarily and then
!      !   then others are stored sequentially according to their proximity to
!      !   the previously stored vertex with ties being inconsequential to the
!      !   end result. This will store the vertices in a ring.) The vertices for
!      !   this facet are also added to the list of unique vertices. Of course,
!      !   when adding, they are compared against the existing vertices in the
!      !   list to prevent double listing in the unique set.
!      write (52,*) "About to makeFacetListUniqueVerticesAndFullEdgeList"
!      call flush (52)
!      call makeFacetListUniqueVerticesAndFullEdgeList
!
!      ! Create a list of edges by traversing the list of facets and
!      !   accumulating each unique sequential pair of vertices plus the final
!      !   first-vertex + last-vertex pair.
!      write (52,*) "About to makeUniqueEdgeList"
!      call flush (52)
!      call makeUniqueEdgeList
!
!
   end subroutine computeBrillouinZones

subroutine makeFacetList (maxBZ)

   implicit none

   ! Define passed dummy parameters.
   integer :: maxBZ

   ! Define local variables.
   integer :: i, j, k
   integer :: onOtherPlane, beyondOtherPlane, oppositeQuadrants
   integer :: numFacets
   integer, dimension(26) :: facetList
   real (kind=double) :: t, d ! See algebraic description below.
   real (kind=double), dimension(3) :: intersectionPoint
   real (kind=double) :: smallThresh

   ! Certain calculations require floating point comparisons that are
   !   inherently not exact on a finite precision computer. Therefore, we
   !   use this threshold to define effective equality.
   smallThresh = 1e-8_double

   ! The goal of the algorithm is to create a list of all the planes that
   !   contribute to the Brillouin zone (and, of course, exclude those planes
   !   that do not contribute). A plane is recognized as a contributor if there
   !   is a line from the origin to the lattice point of the plane that does
   !   not cross through any other plane. Conversely, a plane is recognized to
   !   not contribute to the Brillouin zone if a line from the point on the
   !   plane nearest to the origin passes through any other plane, or if that
   !   same point happens to lie exactly on some other plane.
   ! The algorithm works by first assuming that all 26 planes will contribute.
   !   Then, each plan is reviewed to see if either (a) its nearest point lies
   !   on some other plane or (b) a line from the origin to that same nearest
   !   point passes through another plane.

   ! First, we need to create algebraic expressions for each of the planes
   !   available for us to intersect.
   call makePlaneEquations

   numFacets = 0
   do i = 1, 26

      ! Determine if the current lattice point (that defines some plane) also
      !   happens to lie *on* a plane defined by some other lattice point. If
      !   so, then the plane defined by the current lattice point is
      !   disqualified as a participant in the construction of the Brillouin
      !   zone. Assume that the current lattice point will not sit on some
      !   other plane.
      onOtherPlane = 0
      do j = 1, 26

         ! Do not check against oneself.
         if (i == j) cycle

         if (abs(sum(planeEquation(1:3,j) &
               & * (planeEquation(1:3,i) - planeEquation(1:3,j)))) &
               & < smallThresh) then
            onOtherPlane = 1
write (6,*) "Lattice point i sits on plane j."
            exit
         endif
      enddo

      ! If plane i does not lie on another plane, then we need to next
      !   determine if it lies beyond another plane. (I.e., does the point on
      !   plane i that is nearest to the origin have a line segment between it
      !   and the origin that passes through any other plane?)
      if (onOtherPlane == 0) then
write (6,*) "Checking for position."

         beyondOtherPlane = 0
         do j = 1, 26

            ! Do not check against oneself.
            if (i == j) cycle

            ! The general parametric form for the equation of a line is:
            !   x = x_0 + tu;  y = y_0 + tv;  z = z_0 + tw. Our case is a
            !   special case because the vector <u, v, w> is the same as the
            !   position vector <x_0, y_0, z_0> of the lattice point used to
            !   define plane i. Therefore our line is defined according to:
            !   x = u(1+t);  y = v(1+t);  z = w(1+t).
            ! If we substitute these values for x, y, and z into the definition
            !   for the plane defined by the lattice point j we can solve for
            !   t in terms of u, v, w and l, m, n, which are the Cartesian
            !   coordinates of the position vector used to define plane j.
            ! The equation for plane j is: l(x - l) + m(y - m) + n(z - n) = 0.
            !   Upon substitution we have:
            !   l(u(1+t) - l) + m(v(1+t) - m) + n(w(1+t) - n) = 0.
            ! Solving for t yields:
            !   t = (l(l-u) + m(m-v) + n(n-w)) / (lu + mv + nw)
            !   t = (l^2 - lu + m^2 - mv + n^2 - nw) / (lu + mv + nw)
            !   t = (l^2 + m^2 + n^2 - lu - mv - nw) / (lu + mv + nw)
            !   t = (l^2 + m^2 + n^2 - d) / d; d = (lu + mv + nw)
            !   t = ((l^2 + m^2 + n^2) / d) - 1
            ! Given t we can plug in to the above parametric equations for a
            !   line to get the point at which the line intersects the plane j.
            !   If the origin is between this point and the point that defines
            !   plane i then plane j does not block plane i. On the other hand,
            !   if this point and the point that defines plane i are in the
            !   same octant then we compute the distance from the origin to
            !   each point. If this point is closer to the origin than the
            !   point that defines plane i, then plane j does block plane i and
            !   that ensures that plane i cannot contribute to the Brillouin
            !   zone. If this point is farther from the origin than the point
            !   that defines plane i, then plane j does not block all of plane
            !   i and so does not preclude plane i from contributing to the
            !   Brillouin zone. (At the same time, it also does not ensure that
            !   plane i *does* contribute. That will depend on the other planes
            !   too.)
            d = sum(planeEquation(1:3,i) * planeEquation(1:3,j))
            t = (sum(planeEquation(1:3,i)**2) / d) - 1.0_double
            intersectionPoint(:) = planeEquation(1:3,j) * (1.0_double + t)

            ! Check if the intersectionPoint and the point that defines plane i
            !   are in the same (or opposite) quadrants. We will assume that
            !   the points are in opposite quadrants. If any Cartesian
            !   coordinate component pair has the same sign, then that is a
            !   clear indication that they cannot be in the same quadrant.
            oppositeQuadrants = 1
            do k = 1, 3
               if (((intersectionPoint(k) > 0.0_double) .and. &
                     & (planeEquation(k,i) > 0.0_double)) .or. &
                     & ((intersectionPoint(k) < 0.0_double) .and. &
                     & (planeEquation(k,i) < 0.0_double))) then
                  oppositeQuadrants = 0
                  exit
               endif
            enddo

            ! If the point that defines plane i and the line that goes through
            !   that point and the origin then intersects with plane j in the
            !   opposite quadrant then we can be certain that plane j will not
            !   impede plane i from participating in the construction of the
            !   Brillouin zone.
            if (oppositeQuadrants == 1) then
write (6,*) "Found opposite quadrants", i, j
               cycle
            endif

            ! On the other hand, the intersection point and the point that is
            !   used to define plane i are in the same quadrant, then we need
            !   to determine which one is closer to the origin. If the point
            !   used to define plane i is more distant than the intersection
            !   point, then plane i is "beyond the other plane" and so it
            !   cannot contribute to the Brillouin zone.
            if (sum(planeEquation(1:3,i)**2) > sum(intersectionPoint(:)**2)) &
                  & then
               beyondOtherPlane = 1
               exit
            endif
         enddo
      endif

      ! Plane i has satisfied all the requirements of being included in the
      !   construction of the Brillouin zone if it is not on another plane and
      !   not beyond another plane. Upon satisfaction of the requirements we
      !   increase the number of facets for the Brillouin zone and record the
      !   index number of i.
      if ((onOtherPlane == 0) .and. (beyondOtherPlane == 0)) then
         numFacets = numFacets + 1
         facetList(numFacets) = i
      endif
   enddo

!      ! Determine if a line from the origin and normal to facet i passes
!      !   through any of the facets found so far. If the line does pass through
!      !   any facet found so far, then we turn on the flag to indicate that an
!      !   intersection was found.
!      sameSideAsOrigin = 1
!      do j = 1, numFacets(1)
!
!         ! Get shorthand variables for the coefficients for known facet j.
!         coeffs_j(:) = planeEquation(:,facetList(j,1))
!         latticePointXYZ_j(:) = latticePointXYZ(:,facetList(j,1))
!write (6,fmt="(a3,i3,a4,3d12.4)") "cj(", j, ") = ", coeffs_j(1:3)
!write (6,fmt="(a7,i3,a4,3d12.4)") "lpXYZj(", j, ") = ", latticePointXYZ_j(:)
!
!!write (6,*) coeffs_j(1)*(coeffs_i(1) - coeffs_j(1))
!!write (6,*) coeffs_j(2)*(coeffs_i(2) - coeffs_j(2))
!!write (6,*) coeffs_j(3)*(coeffs_i(3) - coeffs_j(3))
!!write (6,*) sum(coeffs_j(1:3)*(coeffs_i(1:3) - coeffs_j(1:3))), coeffs_j(4)
!         if (abs(sum(coeffs_j(1:3) &
!               & * (latticePointXYZ_i(1:3) - latticePointXYZ_j(1:3))) &
!               & - coeffs_j(4)) < smallThresh) then
!            sameSideAsOrigin = 0
!write (6,*) "Plane i and the origin are on opposite sides of plane j"
!            exit
!         endif
!      enddo
!
!      ! Only if no intersection was found do we append facet i to the list of
!      !   preliminarily acceptable facets.
!      if (sameSideAsOrigin == 1) then
!         numFacets(1) = numFacets(1) + 1
!         facetList(numFacets(1),1) = i
!      endif
!   enddo
!
!
!
!
!   numFacets(2) = 0
!   do i = 1, numFacets(1)
!
!      coeffs_i(:) = planeEquation(:,facetList(i,1))
!      latticePointXYZ_i(:) = latticePointXYZ(:,facetList(i,1))
!write (6,fmt="(a7,i3,a4,3d12.4)") "lpXYZi(", i, ") = ", latticePointXYZ_i(:)
!write (6,fmt="(a,i3,a4,4d12.4)") "ci(", i, ") = ", coeffs_i(1:4)
!
!      sameSideAsOrigin = 1
!      do j = 1, numFacets(1)
!
!         if (i==j) cycle
!
!         coeffs_j(:) = planeEquation(:,facetList(j,1))
!write (6,fmt="(a3,i3,a4,3d12.4)") "cj(", j, ") = ", coeffs_j(1:3)
!         latticePointXYZ_j(:) = latticePointXYZ(:,facetList(j,1))
!write (6,fmt="(a7,i3,a4,3d12.4)") "lpXYZj(", j, ") = ", latticePointXYZ_j(:)
!
!!         if (abs(sum(coeffs_i(1:3)*coeffs_j(1:3)) - coeffs_i(4)) &
!!               & < smallThresh) then
!!            sameSideAsOrigin = 0
!!write (6,*) "Lattice point i sits inside plane j."
!!            exit
!!         endif
!
!!         if (abs(sum(coeffs_i(1:3)*(coeffs_j(1:3) - coeffs_i(1:3))) &
!!               & - coeffs_i(4)) < smallThresh) then
!!            sameSideAsOrigin = 0
!!write (6,*) "Lattice point i sits inside plane j."
!!            exit
!!         endif
!
!         if (abs(sum(coeffs_j(1:3) &
!               & * (latticePointXYZ_i(:) - latticePointXYZ_j(:)))) &
!               & < smallThresh) then
!            sameSideAsOrigin = 0
!write (6,*) "Lattice point i sits inside plane j."
!            exit
!         endif
!      enddo
!
!      if (sameSideAsOrigin == 1) then
!          numFacets(2) = numFacets(2) + 1
!          facetList(numFacets(2),2) = facetList(i,1)
!      endif
!   enddo

   ! At this point we have found all the facets that should be used to make
   !   the Brillouin zone.
   write (6,*) numFacets
   write (6,*) facetList(1:numFacets)
   do i = 1, numFacets
      write (6,*) planeEquation(1:3,facetList(i))
   enddo

end subroutine makeFacetList


!
!
!   subroutine makeFullVertexList
!
!      implicit none
!
!      ! Define local variables.
!      integer :: info
!      integer :: i, j, k
!      integer, dimension (3) :: pivotIndices
!      real (kind=double), dimension (3,3) :: A
!      real (kind=double), dimension (3) :: B, x
!
!      ! First, we need to create algebraic expressions for each of the planes
!      !   available for us to intersect.
!      call makePlaneEquations
!
!      ! Initialize the number of vertices in the full vertex list.
!      numFullVertices = 0
!
!      ! Iterate over all triple sets of planes.
!      do i = 1, 26
!         do j = i+1, 26
!            do k = j+1, 26
!
!               A(:,1) = planeEquation(1:3,i)
!               B(1) = planeEquation(4,i)
!
!               A(:,2) = planeEquation(1:3,j)
!               B(2) = planeEquation(4,j)
!
!               A(:,3) = planeEquation(1:3,k)
!               B(3) = planeEquation(4,k)
!
!               call dgesv(3, 1, A, 3, pivotIndices, B, 3, info)
!               if (info == 0) then
!write (52,*) "Found a valid vertex", i, j, k
!call flush (52)
!
!                  ! We found a valid vertex. Increment our global count.
!                  numFullVertices = numFullVertices + 1
!
!                  ! Store the x, y, z coordinates of this vertex.
!                  fullVertexList(numFullVertices)%coord(:) = B(pivotIndices(:))
!
!                  ! Store the indices that define the set of planes used to
!                  !   make the current vertex.
!                  fullVertexList(numFullVertices)%planeTripleIndex(1) = i
!                  fullVertexList(numFullVertices)%planeTripleIndex(2) = j
!                  fullVertexList(numFullVertices)%planeTripleIndex(3) = k
!               else
!write (52,*) "Found an invalid vertex", i, j, k
!call flush (52)
!               endif
!            enddo
!         enddo
!      enddo
!write (52,*) "numFullVertices=", numFullVertices
!call flush(52)
!   end subroutine makeFullVertexList
!
   ! Take each plane (defined as being perpendicular to the members of the set
   !   of neighboring recipical lattice points) and compute the equation that
   !   defines it. Note: given a Cartesian coordinate vector with magnitudes
   !   (l, m, n) in the x, y, z directions, the equation for the plane that
   !   goes through the point (l, m, n) is: l*(x-l) + m*(y-m) + n*(z-n) = 0.
   !   That is: l*x + m*y + n*z - l^2 - m^2 - n^2 = 0 or also expressed as:
   !   l*x + m*y + n*z = d with d = l^2 + m^2 + n^2. To "store" the equation
   !   for each plane we will simply hold the Cartesian <l, m, n> vector and
   !   the value of d.
   subroutine makePlaneEquations

      implicit none

      ! Define local variables.
      integer :: i, j, k, s, planeIndex

      ! The plane equations must be expressed in Cartesian coordinates.
      !   Therefore, for each of the 26 possible (we exclude the origin zero
      !   vector) vectors we need to obtain its x, y, z coefficients from the
      !   reciprocal space lattice vectors. At the same time, we will make a
      !   single integer index for each plane equation.
      planeIndex = 0
      do i = -1, 1
         do j = -1, 1
            do k = -1, 1

               ! Ignore the zero vector because it does not define a plane.
               if ((i==0) .and. (j==0) .and. (k==0)) then
                  cycle
               endif

               ! Increment the plane index after cycling on 0,0,0.
               planeIndex = planeIndex + 1

               ! Get the s=x, y, z coefficients for the current lattice vector.
               !   Multiply by 1/magnitude because the coefficients are the
               !   coefficients of a normal vector for the plane.
               do s = 1, 3
                  latticePointXYZ(s,planeIndex) = (i*recipLattice(s,1) &
                                                + j*recipLattice(s,2) &
                                                + k*recipLattice(s,3))
               enddo

               ! Store them and the solution (constant) term. Note that the
               !   0.5 factor is applied because the planes are formed at the
               !   midpoint between the lattice point and the origin.
               planeEquation(1:3,planeIndex) = &
                     & latticePointXYZ(1:3,planeIndex) * 0.5_double
               planeEquation(4,planeIndex) = &
                     & sum((latticePointXYZ(:,planeIndex) * 0.5_double)**2)
            enddo
         enddo
      enddo

   end subroutine makePlaneEquations
!
!
!   subroutine makeFacetListUniqueVerticesAndFullEdgeList
!
!      implicit none
!
!      ! Define local variables.
!      integer :: i, j, k, l, m
!      integer :: found
!      integer, dimension (2) :: sharedPlaneArrayIndex
!      integer, dimension (2) :: iPlaneArrayIndex
!      type (vertexType) :: tempVertex
!      type (facetType) :: tempAllVerticesFacet
!      type (facetType) :: tempUniqueVerticesFacet
!
!      ! Initialize the count of the number of facets and edges.
!      numBZFacets(1) = 0
!      numBZEdges(1) = 0
!      numBZVertices(1) = 0
!      numFullEdges = 0
!
!      ! Consider each plane and then determine how many vertices are associated
!      !   with it. If the plane has at least three, then we record it as a
!      !   facet. We also record that vertex into a unique vertex list.
!      do i = 1, 26
!write (52,*) "Iteration i=", i
!call flush (52)
!
!         ! Initialize a count of the number of vertices for the current plane.
!         tempAllVerticesFacet%numVertices = 0
!         tempUniqueVerticesFacet%numVertices = 0
!
!         ! Consider each vertex in turn and determine if the current plane
!         !   shows up in the list of planes that define it (the
!         !   planeTripleIndex). If it does, then we record this vertex in the
!         !   temp facet. It will be copied later if there are enough vertices
!         !   associated with it.
!         do j = 1, numFullVertices
!            do k = 1, 3
!               if (fullVertexList(j)%planeTripleIndex(k) == i) then
!                  tempAllVerticesFacet%numVertices &
!                        & = tempAllVerticesFacet%numVertices + 1
!                  call copyVertex(tempAllVerticesFacet%vertex( &
!                        & tempAllVerticesFacet%numVertices), &
!                        & fullVertexList(j))
!
!                  ! Determine if this is a unique vertex to add to the facet.
!                  found = 0
!                  do l = 1, tempUniqueVerticesFacet%numVertices
!                     if (vertexCoordsEqual(tempUniqueVerticesFacet%vertex(l), &
!                           & fullVertexList(j))) then
!                        found = l
!                        exit
!                     endif
!                  enddo
!                  if (found == 0) then
!                     tempUniqueVerticesFacet%numVertices &
!                           & = tempUniqueVerticesFacet%numVertices + 1
!                     call copyVertex(tempUniqueVerticesFacet%vertex( &
!                           & tempUniqueVerticesFacet%numVertices), &
!                           & fullVertexList(j))
!                  endif
!                  exit
!               endif
!            enddo
!         enddo
!
!         ! If this plane has < three unique vertices, it cannot be a facet.
!         if (tempUniqueVerticesFacet%numVertices < 3) then
!            cycle
!         endif
!
!         ! At this point, we have found a plane with enough vertices.
!write (52,*) "Found a plane with tempUniqueVerticesFacet%numVertices=", &
!   & tempUniqueVerticesFacet%numVertices
!call flush (52)
!
!         ! We know that there are at least three vertices that are unique.
!         !   However, the specific vertices pulled from the list of all
!         !   vertices may not be the "right" ones. Now that we know all of
!         !   the vertices in this facet we will pull out the ones that form
!         !   a ring and put them in the unique vertices list.
!
!         ! Find a pair of vertices that have two shared planes that are the
!         !   same.
!         do j = 1, numFullVertices
!
!            do k = j+1, numFullVertices
!
!               if (verticesConnected(fullVertexList(j),fullVertexList(k))) then
!
!               endif
!            enddo
!         enddo
!
!         ! Increment the count of facets in the system and record all of its
!         !   relevant information.
!         numBZFacets(1) = numBZFacets(1) + 1
!         facetList(numBZFacets(1))%numVertices &
!                  & = tempUniqueVerticesFacet%numVertices
!
!         ! All of the vertices in this facet will have plane "i" as one of
!         !   the planes listed in their planeTripleIndex. What we want to
!         !   do now though is reorder the vertices according to a ring-
!         !   sorting so that they define the perimeter of a polygon. We will
!         !   do that by finding the *other* common plane that pairs of
!         !   vertices share. (I.e., all vertices share one common plane, but
!         !   other planes can be present in *at most* two vertices. Those
!         !   two vertices need to be placed in adjacent array indices.
!         ! So, vertex one is assumed to be settled in the ring. Then, we
!         !   will start adding vertices to the ring at index two and then we
!         !   add index three, etc.
!         call copyVertex(facetList(numBZFacets(1))%vertex(1), &
!               & tempAllVerticesFacet%vertex(1))
!         do j = 2, tempFacet%numVertices
!            
!            ! Find the vertex that has the same (non-"i") plane as the
!            !   vertex before it (j-1).
!            do k = j, tempFacet%numVertices
!               
!               ! Check each of the non-"i" planes in the planeTripleIndex of
!               !   both vertices to find a match. We assume that no match
!               !   will be found.
!               sharedPlaneArrayIndex(:) = 0
!               iPlaneArrayIndex(:) = 0
!               do l = 1, 3 ! planeTripleIndex for the j loop
!
!                  ! Record the iPlaneIndex for the j vertex when we find it.
!                  !   Then cycle, because the iPlaneIndex is not the plane
!                  !   we want to match.
!                  if (tempFacet%vertex(j-1)%planeTripleIndex(l) == i) then
!                     iPlaneArrayIndex(1) = l
!                     cycle
!                  endif
!
!                  do m = 1, 3 ! planeTripleIndex for the k loop
!
!                     ! Record the iPlaneIndex for the k vertex. Then cycle
!                     !   because we will never need to match the iPlaneIndex
!                     if (tempFacet%vertex(k)%planeTripleIndex(m) == i) then
!                        iPlaneArrayIndex(2) = m
!                        cycle
!                     endif
!
!                     ! Look for a match and turn on a flag if we find it.
!                     !   If a match is found, then also record the edge
!                     !   that it defines. Note that we don't exit because
!                     !   we might still need to find the iPlaneIndex
!                     if (tempFacet%vertex(j-1)%planeTripleIndex(l) == &
!                           & tempFacet%vertex(k)%planeTripleIndex(m)) then
!                        sharedPlaneArrayIndex(1) = l
!                        sharedPlaneArrayIndex(2) = m
!                     endif
!                  enddo
!               enddo
!
!               ! If we didn't find a match, then go to the next vertex. We
!               !   only need to test one of the array index values.
!               if (sharedPlaneArrayIndex(1) == 0) then
!                  cycle
!               endif
!
!               ! At this point we can say that we found a match. So, we
!               !   process it.
!
!               ! Record the edge in the fullEdgeList.
!               numFullEdges = numFullEdges + 1
!
!               call copyVertex(fullEdgeList(numFullEdges)%vertex(1), &
!                     & tempFacet%vertex(j-1))
!               call copyVertex(fullEdgeList(numFullEdges)%vertex(2), &
!                     & tempFacet%vertex(k))
!
!               ! Store the plane indices of the two planes that define the
!               !   two vertices of this edge.
!               fullEdgeList(numFullEdges)%sharedPlaneIndex(1) = i
!               fullEdgeList(numFullEdges)%sharedPlaneIndex(2) &
!                     & = tempFacet%vertex(j-1)%planeTripleIndex(l)
!
!               ! Store the array indices of the planeTripleIndex associated
!               !    with the vertices' planes.
!               fullEdgeList(numFullEdges) &
!                     & %sharedPlaneTripleIndexIndices(1,1) &
!                     & = sharedPlaneArrayIndex(1)
!               fullEdgeList(numFullEdges) &
!                     & %sharedPlaneTripleIndexIndices(1,2) &
!                     & = iPlaneArrayIndex(1)
!               fullEdgeList(numFullEdges) &
!                     & %sharedPlaneTripleIndexIndices(2,1) &
!                     & = sharedPlaneArrayIndex(2)
!               fullEdgeList(numFullEdges) &
!                     & %sharedPlaneTripleIndexIndices(2,2) &
!                     & = iPlaneArrayIndex(2)
!
!               ! We can now put the vertex that defines the next point on
!               !   the ring into its correct position in the facet vertex
!               !   array.
!               call copyVertex(facetList(numBZFacets(1))%vertex(j), &
!                     & tempFacet%vertex(k))
!               call copyVertex(tempVertex, tempFacet%vertex(j))
!               call copyVertex(tempFacet%vertex(j), tempFacet%vertex(k))
!               call copyVertex(tempFacet%vertex(k), tempVertex)
!
!               ! Once stored, we can exit this loop.
!               exit
!            enddo
!         enddo
!
!         ! Add each of the vertices to the list of unique vertices.
!         do j = 1, facetList(numBZFacets(1))%numVertices
!            found = 0
!            do k = 1, numBZVertices(1)
!               if (verticesEqual(uniqueVertexList(k), &
!                        & facetList(numBZFacets(1))%vertex(j))) then
!                  found = k
!               endif
!            enddo
!            if (found == 0) then
!               numBZVertices(1) = numBZVertices(1) + 1
!               call copyVertex(uniqueVertexList(numBZVertices(1)), &
!                     & facetList(numBZFacets(1))%vertex(j))
!            endif
!         enddo
!      enddo
!
!   end subroutine makeFacetListUniqueVerticesAndFullEdgeList
!
!
!   function verticesEqual(vertex1, vertex2)
!
!      implicit none
!
!      ! Define dummy parameters.
!      type (vertexType) :: vertex1, vertex2
!
!      ! Define the return variable.
!      logical :: verticesEqual
!
!      ! Define local variables.
!      integer :: i, j
!      real (kind=double) :: threshold
!
!      ! Assume that the two vertices are equal and define the threshold.
!      verticesEqual = .true.
!      threshold = 0.000001
!
!      do i = 1, 3
!         if (abs(vertex1%coord(i) - vertex2%coord(i)) > threshold) then
!            verticesEqual = .false.
!            return
!         endif
!
!         if (vertex1%planeTripleIndex(i) /= vertex2%planeTripleIndex(i)) then
!            verticesEqual = .false.
!            return
!         endif
!      enddo
!
!   end function verticesEqual
!
!
!   function vertexCoordsEqual(vertex1, vertex2)
!
!      implicit none
!
!      ! Define dummy parameters.
!      type (vertexType) :: vertex1, vertex2
!
!      ! Define the return variable.
!      logical :: vertexCoordsEqual
!
!      ! Define local variables.
!      integer :: i, j
!      real (kind=double) :: threshold
!
!      ! Assume that the two vertices are equal and define the threshold.
!      verticesEqual = .true.
!      threshold = 0.000001
!
!      do i = 1, 3
!         if (abs(vertex1%coord(i) - vertex2%coord(i)) > threshold) then
!            verticesEqual = .false.
!            return
!         endif
!      enddo
!
!   end function vertexCoordsEqual
!
!
!   subroutine copyVertex(vertex1, vertex2)
!
!      implicit none
!
!      ! Define dummy parameters.
!      type (vertexType) :: vertex1, vertex2
!
!      ! Perform the copy.
!      vertex1%coord(:) = vertex2%coord(:)
!      vertex1%planeTripleIndex(:) = vertex2%planeTripleIndex(:)
!
!   end subroutine copyVertex
!
!   
!   subroutine makeUniqueEdgeList
!
!      implicit none
!
!      ! Define local variables.
!      integer :: i, j
!      integer :: found
!
!      do i = 1, numFullEdges
!         found = 0
!         do j = 1, numBZEdges(1)
!            if (edgesEqual(fullEdgeList(i), uniqueEdgeList(j))) then
!               found = 1
!               exit
!            endif
!         enddo
!
!         if (found == 0) then
!            numBZEdges(1) = numBZEdges(1) + 1
!            call copyEdge(uniqueEdgeList(numBZEdges(1)), fullEdgeList(i))
!         endif
!      enddo
!   end subroutine makeUniqueEdgeList
!
!
!   function edgesEqual (edge1, edge2)
!
!      implicit none
!
!      ! Define dummy parameters.
!      type (edgeType) :: edge1, edge2
!
!      ! Define the return variable.
!      logical :: edgesEqual
!
!      ! Define local variables.
!      integer :: i, j
!      real (kind=double) :: threshold
!
!      if ((verticesEqual(edge1%vertex(1), edge2%vertex(1)) .and. &
!            & verticesEqual(edge1%vertex(2), edge2%vertex(2))) .or. &
!            & (verticesEqual(edge1%vertex(1), edge2%vertex(2)) .and. &
!            & verticesEqual(edge1%vertex(2), edge2%vertex(1)))) then
!         edgesEqual = .true.
!      else
!         edgesEqual = .false.
!      endif
!
!   end function edgesEqual
!
!
!   subroutine copyEdge(edge1, edge2)
!
!      implicit none
!
!      ! Define dummy parameters.
!      type (edgeType) :: edge1, edge2
!
!      call copyVertex(edge1%vertex(1), edge2%vertex(1))
!      call copyVertex(edge1%vertex(2), edge2%vertex(2))
!
!      edge1%sharedPlaneIndex(:) = edge2%sharedPlaneIndex(:)
!      edge1%sharedPlaneTripleIndexIndices(:,:) &
!            & = edge2%sharedPlaneTripleIndexIndices(:,:)
!
!   end subroutine copyEdge


   subroutine printBrillouinZones (maxBZ)

      ! Make sure no funny variables are defined.
      implicit none

      ! Define dummy variables.
      integer :: maxBZ ! Highest BZ we will print.

      ! Define local variables.
      integer :: i, j, k ! Loop indices.

      ! Step through each of the Brillouin zones we will make.
      do i = 1, maxBZ

         ! Print all of the vertices.
         write (51+i,*) "Vertices"
         do j = 1, numBZVertices(i)
            write (51+i,*) uniqueVertexList(j)%coord(:)
         enddo

         ! Print all of the edges.
         write (51+i,*) "Edges"
         do j = 1, numBZEdges(i)
            write (51+i,*) uniqueEdgeList(j)%vertex(1)%coord(:)
            write (51+i,*) uniqueEdgeList(j)%vertex(2)%coord(:)
         enddo

         ! Print all of the facets.
         write (51+i,*) "Facets"
         do j = 1, numBZFacets(i)
            do k = 1, facetList(numBZFacets(i))%numVertices
               write (51+i,*) facetList(numBZFacets(i))%vertex(k)%coord(:)
            enddo

            ! Write the first vertex again to close the polygon
            write (51+i,*) facetList(numBZFacets(i))%vertex(1)%coord(:)
         enddo
      enddo

   end subroutine printBrillouinZones

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


   ! Begin list of module subroutines.
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
!write (52,*) "foldedKP=",foldedKPoint(:)
!write (52,*) "abcMeshKP=",abcMeshKPoints(:,j)
!               call saveKPoint(j)
!write (52,*) "Matched #2"
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
!write (52,*) "Matched #3"
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


module CommandLine_O

   ! Make sure no funny variables are defined.
   implicit none

   ! Integer to track which number command line argument is to be read in next.
   integer :: nextArg

   ! Flags to output variaous kpoint datasets.
   integer :: doBrillouinZone

   contains

   subroutine parseCommandLine
   
      ! Make sure that there are no accidental variable declarations.
      implicit none
   
      ! Define the local variables used for reading.
      character*25 :: commandBuffer
   
      ! Make sure that the first command line parameter to be read is #1.
      nextArg = 1
   
      ! Initialize the command line parameters
      doBrillouinZone = 0
   
      ! Get the command line argument that flags whether or not to output
      !   descriptive information about the Brillouin zone. The number supplied
      !   is the highest Brillouin zone number to analyze. (If a 1 is given,
      !   then we analyze the first BZ. If a 2 is given, then we analyze the
      !   first and second BZ. If a 3 is given, then we do the 1st, 2nd, and
      !   3rd. Etc.)
      call getarg (nextArg, commandBuffer)
      nextArg = nextArg + 1
      read (commandBuffer,*) doBrillouinZone
   
   end subroutine parseCommandLine

end module CommandLine_O


program makekpoints

   ! This program is used to generate a set of kpoints for a given crystal and
   !   kpoint mesh request.

   ! Import necessary parameter modules.
   use O_Kinds
   use O_Constants

   ! Import necessary modules.
   use PointGroupOperations_O
   use Lattice_O
   use KPointMesh_O
   use CommandLine_O


   ! Make sure no funny variables are defined.
   implicit none

   ! Define local variables.
   integer :: i, ioerr
   character*25 :: filename

   ! Read command line parameters.
   call parseCommandLine

   ! Open files for reading and writing.
   open (unit=50, file='kpSpecs.dat', form='formatted')
   open (unit=51, file='kpSpecs.out', form='formatted')
   ioerr = 0
   do i = 1, doBrillouinZone
      write (filename,fmt="(a3,i1)") "BZ.", i
      open (unit=51+i, file=filename, form='formatted', status='new', &
            & iostat=ioerr)
      if (ioerr /= 0) then
         exit
      endif
   enddo

   ! If this program is called multiple times (which it might be by makeinput)
   !   then this will prevent recalculation of the Brillouin zone figure.
   if (ioerr /= 0) then
      doBrillouinZone = 0
   endif

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

   ! In the case that the user wants information about a the Brillouin zone,
   !   then we compute and print it.
   if (doBrillouinZone > 0) then
write (6,*) "About to compute"
call flush(6)
      call computeBrillouinZones(doBrillouinZone)
write (6,*) "About to print"
call flush(6)
      !call printBrillouinZones(doBrillouinZone)
   endif

   close (50)
   close (51)
   do i = 1, doBrillouinZone
      close (51+i)
   enddo

end program makekpoints
