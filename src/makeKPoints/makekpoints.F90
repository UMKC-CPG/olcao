module CommandLine_O

   use O_Kinds

   ! Make sure no funny variables are defined.
   implicit none

   ! Integer to track which number command line argument is to be read in next.
   integer :: nextArg

   ! Flags to control the supplementary files to output.
   integer :: doBrillouinZone
   real (kind=double) :: scaleFactor

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
      scaleFactor = 10.0_double
   
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

      enddo ! i=1,numPointOps

      ! Read past all the space group operations that are just point group
      !   operations plus extra translation.
      do i = 1, numSpaceOps - numPointOps
         do j = 1,5
            read (fileUnit,*)
         enddo
      enddo

   end subroutine readPointOps


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
   real (kind=double), dimension (3) :: recipMag ! Along reciprocal abc axes.
   real (kind=double), dimension (3) :: realMag ! Along real abc axes
   real (kind=double), dimension (3) :: recipAngle ! Radians
   real (kind=double), dimension (3) :: recipAngleDeg

   type vertexType
      real (kind=double), dimension (3) :: coord
      integer, dimension (3) :: planeTripleIndex
   end type vertexType
   type (vertexType), dimension (100) :: uniqueVertexList

   type edgeType
      type (vertexType), dimension (2) :: vertex
   end type edgeType
   type (edgeType), dimension (1000) :: uniqueEdgeList

   type facetType
      integer :: pointID
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

      ! Define local variables
      integer :: i

      read (fileUnit,*) realLattice(:,:) ! Reads xyz of a, then b, then c

      ! Compute the magnitudes of the real lattice abc vectors.
      do i = 1, 3
         realMag(i) = sqrt(sum(realLattice(:,i)**2))
      enddo

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

      ! Compute the reciprocal cell vector magnitudes.
      do i = 1, 3
         recipMag(i) = sqrt(sum(recipLattice(:,i)**2))
      enddo

      recipAngle(1) = acos((recipLattice(1,2) * recipLattice(1,3) + &
                          & recipLattice(2,2) * recipLattice(2,3) + &
                          & recipLattice(3,2) * recipLattice(3,3)) / &
                          & (recipMag(2)*recipMag(3)))
      recipAngle(2) = acos((recipLattice(1,1) * recipLattice(1,3) + &
                          & recipLattice(2,1) * recipLattice(2,3) + &
                          & recipLattice(3,1) * recipLattice(3,3)) / &
                          & (recipMag(1)*recipMag(3)))
      recipAngle(3) = acos((recipLattice(1,1) * recipLattice(1,2) + &
                          & recipLattice(2,1) * recipLattice(2,2) + &
                          & recipLattice(3,1) * recipLattice(3,2)) / &
                          & (recipMag(1)*recipMag(2)))
      recipAngleDeg(1) = recipAngle(1) * 180.0_double / pi
      recipAngleDeg(2) = recipAngle(2) * 180.0_double / pi
      recipAngleDeg(3) = recipAngle(3) * 180.0_double / pi

   end subroutine computeLatticeData


   subroutine printLatticeData (maxBZ)

      use CommandLine_O

      implicit none

      ! Define passed parameters.
      integer :: maxBZ

      ! Define local variables.
      integer :: h, i, j, k, l
      integer :: vertexCount
      type (vertexType), dimension (8) :: recipLattVertices
      type (edgeType), dimension(12) :: recipLattEdges
      type (facetType), dimension(6) :: recipLattFaces

      ! Create the list of vertices
      vertexCount = 0
      do i = 0, 1
      do j = 0, 1
      do k = 0, 1
         vertexCount = vertexCount + 1
         do l = 1, 3
            recipLattVertices(vertexCount)%coord(l) = &
                  & i * recipLattice(l,1) + &
                  & j * recipLattice(l,2) + &
                  & k * reciplattice(l,3)
         enddo
         recipLattVertices(vertexCount)%coord(:) = &
               & recipLattVertices(vertexCount)%coord(:) * scaleFactor
      enddo
      enddo
      enddo

      ! Create the list of edges.
      call copyVertex(recipLattVertices(1),recipLattEdges(1)%vertex(1))
      call copyVertex(recipLattVertices(2),recipLattEdges(1)%vertex(2))
      call copyVertex(recipLattVertices(1),recipLattEdges(2)%vertex(1))
      call copyVertex(recipLattVertices(3),recipLattEdges(2)%vertex(2))
      call copyVertex(recipLattVertices(2),recipLattEdges(3)%vertex(1))
      call copyVertex(recipLattVertices(4),recipLattEdges(3)%vertex(2))
      call copyVertex(recipLattVertices(3),recipLattEdges(4)%vertex(1))
      call copyVertex(recipLattVertices(4),recipLattEdges(4)%vertex(2))
      call copyVertex(recipLattVertices(5),recipLattEdges(5)%vertex(1))
      call copyVertex(recipLattVertices(6),recipLattEdges(5)%vertex(2))
      call copyVertex(recipLattVertices(5),recipLattEdges(6)%vertex(1))
      call copyVertex(recipLattVertices(7),recipLattEdges(6)%vertex(2))
      call copyVertex(recipLattVertices(6),recipLattEdges(7)%vertex(1))
      call copyVertex(recipLattVertices(8),recipLattEdges(7)%vertex(2))
      call copyVertex(recipLattVertices(7),recipLattEdges(8)%vertex(1))
      call copyVertex(recipLattVertices(8),recipLattEdges(8)%vertex(2))
      call copyVertex(recipLattVertices(2),recipLattEdges(9)%vertex(1))
      call copyVertex(recipLattVertices(6),recipLattEdges(9)%vertex(2))
      call copyVertex(recipLattVertices(4),recipLattEdges(10)%vertex(1))
      call copyVertex(recipLattVertices(8),recipLattEdges(10)%vertex(2))
      call copyVertex(recipLattVertices(1),recipLattEdges(11)%vertex(1))
      call copyVertex(recipLattVertices(5),recipLattEdges(11)%vertex(2))
      call copyVertex(recipLattVertices(3),recipLattEdges(12)%vertex(1))
      call copyVertex(recipLattVertices(7),recipLattEdges(12)%vertex(2))


      ! Create the list of faces.
      recipLattFaces(1)%numVertices = 4
      call copyVertex(recipLattVertices(1),recipLattFaces(1)%vertex(1))
      call copyVertex(recipLattVertices(2),recipLattFaces(1)%vertex(2))
      call copyVertex(recipLattVertices(4),recipLattFaces(1)%vertex(3))
      call copyVertex(recipLattVertices(3),recipLattFaces(1)%vertex(4))
      recipLattFaces(2)%numVertices = 4
      call copyVertex(recipLattVertices(1),recipLattFaces(2)%vertex(1))
      call copyVertex(recipLattVertices(2),recipLattFaces(2)%vertex(2))
      call copyVertex(recipLattVertices(6),recipLattFaces(2)%vertex(3))
      call copyVertex(recipLattVertices(5),recipLattFaces(2)%vertex(4))
      recipLattFaces(3)%numVertices = 4
      call copyVertex(recipLattVertices(1),recipLattFaces(3)%vertex(1))
      call copyVertex(recipLattVertices(3),recipLattFaces(3)%vertex(2))
      call copyVertex(recipLattVertices(7),recipLattFaces(3)%vertex(3))
      call copyVertex(recipLattVertices(5),recipLattFaces(3)%vertex(4))
      recipLattFaces(4)%numVertices = 4
      call copyVertex(recipLattVertices(2),recipLattFaces(4)%vertex(1))
      call copyVertex(recipLattVertices(4),recipLattFaces(4)%vertex(2))
      call copyVertex(recipLattVertices(8),recipLattFaces(4)%vertex(3))
      call copyVertex(recipLattVertices(6),recipLattFaces(4)%vertex(4))
      recipLattFaces(5)%numVertices = 4
      call copyVertex(recipLattVertices(3),recipLattFaces(5)%vertex(1))
      call copyVertex(recipLattVertices(4),recipLattFaces(5)%vertex(2))
      call copyVertex(recipLattVertices(8),recipLattFaces(5)%vertex(3))
      call copyVertex(recipLattVertices(7),recipLattFaces(5)%vertex(4))
      recipLattFaces(6)%numVertices = 4
      call copyVertex(recipLattVertices(5),recipLattFaces(6)%vertex(1))
      call copyVertex(recipLattVertices(6),recipLattFaces(6)%vertex(2))
      call copyVertex(recipLattVertices(8),recipLattFaces(6)%vertex(3))
      call copyVertex(recipLattVertices(7),recipLattFaces(6)%vertex(4))


      ! Print the reciprocal lattice vertices.
      do h = 1, maxBZ
         write (51+h,fmt="(a26)") "recip_lattice_vertices = ["
         do i = 1, 8
            write (51+h,advance="NO",fmt="(a1)") "["
            write (51+h,advance="NO",fmt="(f16.12, a2)") &
                  & recipLattVertices(i)%coord(1), ", "
            write (51+h,advance="NO",fmt="(f16.12, a2)") &
                  & recipLattVertices(i)%coord(2), ", "
            write (51+h,advance="NO",fmt="(f16.12)") &
                  & recipLattVertices(i)%coord(3)
            if (i < 8) then
               write (51+h,fmt="(a3)") "],"
            else
               write (51+h,fmt="(a2)") "]]"
            endif
         enddo
      enddo

      ! Print the reciprocal lattice edges.
      do h = 1, maxBZ
         write (51+h,fmt="(a23)") "recip_lattice_edges = ["
         do i = 1, 12
            write (51+h,advance="NO",fmt="(a2)") "[["
            write (51+h,advance="NO",fmt="(f16.12, a2)") &
               & recipLattEdges(i)%vertex(1)%coord(1), ", "
            write (51+h,advance="NO",fmt="(f16.12, a2)") &
               & recipLattEdges(i)%vertex(1)%coord(2), ", "
            write (51+h,advance="NO",fmt="(f16.12, a3)") &
               & recipLattEdges(i)%vertex(1)%coord(3), "], "
            write (51+h,advance="NO",fmt="(a1)") "["
            write (51+h,advance="NO",fmt="(f16.12, a2)") &
               & recipLattEdges(i)%vertex(2)%coord(1), ", "
            write (51+h,advance="NO",fmt="(f16.12, a2)") &
               & recipLattEdges(i)%vertex(2)%coord(2), ", "
            write (51+h,advance="NO",fmt="(f16.12)") &
               & recipLattEdges(i)%vertex(2)%coord(3)

            if (i < 12) then
               write (51+h,fmt="(a4)") "]], "
            else
               write (51+h,fmt="(a3)") "]]]"
            endif
         enddo
      enddo

      ! Print the reciprocal lattice faces.
      do h = 1, maxBZ
         write (51+h,fmt="(a23)") "recip_lattice_faces = ["
         do i = 1, 6
            write (51+h,advance="NO",fmt="(a1)") "["
            do j = 1, recipLattFaces(i)%numVertices
               write (51+h,advance="NO",fmt="(a1)") "["
               write (51+h,advance="NO",fmt="(f16.12, a2)") &
                     & recipLattFaces(i)%vertex(j)%coord(1), ", "
               write (51+h,advance="NO",fmt="(f16.12, a2)") &
                     & recipLattFaces(i)%vertex(j)%coord(2), ", "
               write (51+h,advance="NO",fmt="(f16.12, a3)") &
                     & recipLattFaces(i)%vertex(j)%coord(3), "], "
            enddo

            ! Write the first vertex again to close the polygon.
            write (51+h,advance="NO",fmt="(a1)") "["
            write (51+h,advance="NO",fmt="(f16.12, a2)") &
                  & recipLattFaces(i)%vertex(1)%coord(1), ", "
            write (51+h,advance="NO",fmt="(f16.12, a2)") &
                  & recipLattFaces(i)%vertex(1)%coord(2), ", "
            write (51+h,advance="NO",fmt="(f16.12, a1)") &
                  & recipLattFaces(i)%vertex(1)%coord(3), "]"

            if (i < 6) then
               write (51+h,fmt="(a2)") "],"
            else
               write (51+h,fmt="(a2)") "]]"
            endif
         enddo
      enddo
   end subroutine printLatticeData


   subroutine copyVertex(vertex1, vertex2)

      implicit none

      ! Define dummy parameters.
      type (vertexType) :: vertex1, vertex2

      ! Perform the copy.
      vertex2%coord(:) = vertex1%coord(:)
      vertex2%planeTripleIndex(:) = vertex1%planeTripleIndex(:)

   end subroutine copyVertex


end module Lattice_O



#ifdef COMPVIS
module VPython_O

   ! Use necessary modules.
   use Lattice_O

   ! Make sure that no unwanted variables are created accidentally.
   implicit none

   ! Define module data.
   real (kind=double) :: sphereRadius
   real (kind=double) :: cameraRadius
   real (kind=double) :: planeThickness
   real (kind=double) :: planeExtent

   contains

   subroutine makeVPythonHeader(maxBZ)

      implicit none

      ! Define passed parameters
      integer :: maxBZ

      ! Define local variables
      integer :: i

      ! Compute or define the sizes of things.
      sphereRadius = minval(recipMag(:)) / 10.0_double
      cameraRadius = maxval(recipMag(:)) * 4.0_double
      planeThickness = 0.03
      planeExtent = maxval(recipMag(:)) * 2.0_double

      do i = 1, maxBZ

         write (61+i,fmt="(a)") "#!/usr/bin/env python3"
         write (61+i,fmt="(a)")
         write (61+i,fmt="(a)") "import vpython as vp"
         write (61+i,fmt="(a)")
         write (61+i,fmt="(a)") "# Define global variables"
         write (61+i,fmt="(a)") "# Time controls"
         write (61+i,fmt="(a)") "delta_t = 0.5"
         write (61+i,fmt="(a)") "frame_rate = 30"
         write (61+i,fmt="(a)")
         write (61+i,fmt="(a)") "# Canvas parameters"
         write (61+i,fmt="(a)") "w = 1200"
         write (61+i,fmt="(a)") "h = 675"
         write (61+i,fmt="(a)")
         write (61+i,fmt="(a)") "# Camera parameters"
         write (61+i,fmt="(a,sp,e12.5)") "camera_radius = ", cameraRadius
         write (61+i,fmt="(a)")
         write (61+i,fmt="(a)") "# Object nature"
         write (61+i,fmt="(a,sp,e12.5)") "plane_thickness = ", planeThickness
         write (61+i,fmt="(a,sp,e12.5)") "sphere_radius = ", sphereRadius
         write (61+i,fmt="(a)")
         write (61+i,fmt="(a)") "# Object lists"
         write (61+i,fmt="(a)") "recip_point = []"
         write (61+i,fmt="(a)") "recip_plane = []"
         write (61+i,fmt="(a)") "bz_point = []"
         write (61+i,fmt="(a)") "bz_plane = []"
         write (61+i,fmt="(a)") "facet_vertex = []"
         write (61+i,fmt="(a)")
         write (61+i,fmt="(a)") "# Create the initial scene"
         write (61+i,fmt="(a)")
         write (61+i,fmt="(a)") "# Create the canvas"
         write (61+i,fmt="(a,i1,a)") &
"scene = vp.canvas(title='Making BZ ", i, "', x=0, y=0, width=w, height=h,"
         write (61+i,fmt="(a,i1,a)") &
"                  center=vp.vector(0, 0, 0), background=vp.vector(0, 0, 0),"
         write (61+i,fmt="(a)") &
"                  autoscale=False)"
         write (61+i,fmt="(a)")
         write (61+i,fmt="(a)")
         write (61+i,fmt="(a)") "# Setup the camera"
         write (61+i,fmt="(a)") &
"scene.camera.pos = vp.vector(camera_radius, camera_radius, camera_radius)"
         write (61+i,fmt="(a)") &
"scene.camera.axis = -scene.camera.pos"
         write (61+i,fmt="(a)")
         write (61+i,fmt="(a)")
      enddo

   end subroutine makeVPythonHeader


   subroutine makeVPythonRecipLattice(maxBZ)

      ! Make things safe
      implicit none

      ! Define passed parameters
      integer :: maxBZ

      ! Define local variables
      integer :: h, i, j, k, s, pointCount
      character*8 :: formatString

      ! Create a sphere at each reciprocal lattice point. Make 0,0,0 separate.
      do h = 1, maxBZ
         pointCount = 0
         do i = -1, 1
            do j = -1, 1
               do k = -1, 1

                  ! Cycle at the origin to make it separately at the end.
                  if ((i==0) .and. (j==0) .and. (k==0)) then
                     cycle
                  endif

                  ! Increment the count of lattice points.
                  pointCount = pointCount + 1
                  if (pointCount <= 10) then
                     write (formatString,fmt="(a)") "(i0.1,a)"
                  else
                     write (formatString,fmt="(a)") "(i0.2,a)"
                  endif

                  ! Print the sphere at this lattice site.
                  write (61+h,fmt="(a)") "recip_point.append("
                  write (61+h,advance="no",fmt="(a)") &
"        vp.sphere(pos=vp.vector("
                  do s = 1, 3
                     write(61+h,advance="no",fmt="(sp,e10.3)") & ! sp causes +
                           (i*recipLattice(s,1) &
                           + j*recipLattice(s,2) &
                           + k*recipLattice(s,3))
                     if (s < 3) then
                        write(61+h,advance="no",fmt="(a)") ", "
                     endif
                  enddo
                  write(61+h,fmt="(a)") "),"
                  write(61+h,advance="no",fmt="(a)") &
"                  radius="
                  write(61+h,advance="no",fmt="(sp,e10.3)") sphereRadius
                  write(61+h,fmt="(a)") "))"

                  ! Print a smaller sphere at the associated BZ point site.
                  write (61+h,fmt="(a)") "bz_point.append("
                  write (61+h,advance="no",fmt="(a)") &
"        vp.sphere(pos=vp.vector("
                  do s = 1, 3
                     write(61+h,advance="no",fmt="(sp,e10.3)") & ! sp causes +
                           (i*recipLattice(s,1) &
                           + j*recipLattice(s,2) &
                           + k*recipLattice(s,3)) / 2.0_double
                     if (s < 3) then
                        write(61+h,advance="no",fmt="(a)") ", "
                     endif
                  enddo
                  write(61+h,fmt="(a)") "),"
                  write(61+h,advance="no",fmt="(a)") "                  radius="
                  write(61+h,advance="no",fmt="(sp,e10.3)") sphereRadius / 3.0
                  write(61+h,fmt="(a)") "))"

               enddo ! k
            enddo ! j
         enddo ! i

         ! Add the 0, 0, 0 reciprocal lattice point.
         write(61+h,advance="no",fmt="(a)") &
"recip_point.append(vp.sphere(pos=vp.vector(0, 0, 0), radius="
         write(61+h,advance="no",fmt="(sp,e10.3)") sphereRadius
         write(61+h,fmt="(a)") "))"
      enddo ! h
      

   end subroutine makeVPythonRecipLattice

end module VPython_O
#endif



module BrillouinZones_O

   ! Use necessary modeuls.
   use Lattice_O
#ifdef COMPVIS
   use VPython_O
#endif

   ! Make sure that no variables are implicitly declared
   implicit none

   ! Define module data.

   contains

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

      ! Add vertices to each facet.
      call addFacetVertices(maxBZ)

      ! Sort the order of the vertices for each facet so that they are in
      !   the order of a ring.
      call ringSortFacetVertices(maxBZ)

      ! Create a list of edges by traversing the list of facets and
      !   accumulating each unique sequential pair of vertices plus the final
      !   first-vertex + last-vertex pair.
      call makeUniqueEdgeList(maxBZ)

   end subroutine computeBrillouinZones



   subroutine makeFacetList (currBZ)

      ! Make sure that no variables are accidentally declared.
      implicit none

      ! Define passed dummy parameters.
      integer :: currBZ

      ! Define local variables.
      integer :: i, j
      integer :: onOtherPlane, beyondOtherPlane
#ifdef COMPVIS
      character*8 :: formatString_i, formatString_j
#endif

      ! The goal of the algorithm is to create a list of all the planes that
      !   contribute to the Brillouin zone (and, of course, exclude those
      !   planes that do not contribute). A plane is recognized as a
      !   contributor if there is a line from the origin to the lattice point
      !   of the plane that does not cross through any other plane. Conversely,
      !   a plane is recognized to not contribute to the Brillouin zone if a
      !   line from the point on the plane nearest to the origin passes through
      !   any other plane, or if that same point happens to lie exactly on some
      !   other plane.
      ! The algorithm works by first assuming that all 26 planes will
      !   contribute. Then, each plane is reviewed to see if either (a) its
      !   nearest point lies on some other plane or (b) a line from the origin
      !   to that same nearest point passes through another plane.

      ! First, we need to create algebraic expressions for each of the planes
      !   available for us to intersect.
      call makePlaneEquations (currBZ)

      ! Now, consider each plane in turn.
      numBZFacets(currBZ) = 0
      do i = 1, 26

#ifdef COMPVIS
         ! Prepare a format string for array indices.
         if (i < 11) then
            write (formatString_i,fmt="(a)") "(i0.1,a)"
         else
            write (formatString_i,fmt="(a)") "(i0.2,a)"
         endif
         write (61+currBZ,advance="no",fmt="(a)") "bz_point["
         write (61+currBZ,fmt=formatString_i) i-1, "].color = vp.color.blue"
#endif
         ! Determine if the current lattice point, i, (that defines some plane)
         !   also happens to lie *on* a plane defined by some other lattice
         !   point, j. If so, then the plane defined by the current lattice
         !   point is disqualified as a participant in the construction of the
         !   Brillouin zone. Assume that the current lattice point will not
         !   sit on some other plane.
         onOtherPlane = 0
         do j = 1, 26

            ! Do not check against oneself.
            if (i == j) cycle
!#ifdef COMPVIS
!               ! Prepare a format string and make plane j visible.
!               if (j < 11) then
!                  write (formatString_j,fmt="(a)") "(i0.1,a)"
!               else
!                  write (formatString_j,fmt="(a)") "(i0.2,a)"
!               endif
!               write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
!               write (61+currBZ,fmt=formatString_j) j-1, "].visible = True"
!#endif

            if (abs(sum(planeEquation(1:3,j) &
                  & * (planeEquation(1:3,i) - planeEquation(1:3,j)))) &
                  & < smallThresh) then
               onOtherPlane = 1
!#ifdef COMPVIS
!               write (61+currBZ,fmt="(a)") "vp.sleep(0.010)"
!               write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
!               write (61+currBZ,fmt=formatString_j) j-1, "].visible = False"
!#endif
               exit
!#ifdef COMPVIS
!            else
!               write (61+currBZ,fmt="(a)") "vp.sleep(0.01)"
!               write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
!               write (61+currBZ,fmt=formatString_j) j-1, &
!"].visible = False"
!#endif
            endif
         enddo ! j loop


         ! If plane i does not lie *on* another plane, then we need to next
         !   determine if it lies *beyond* another plane. (I.e., does the point
         !   on plane i that is nearest to the origin have a line segment
         !   between it and the origin that passes through any other plane?)
         if (onOtherPlane == 0) then

#ifdef COMPVIS
            write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
            write (61+currBZ,fmt=formatString_i) i-1, "].visible = True"
            write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
            write (61+currBZ,fmt=formatString_i) i-1, "].color = vp.color.white"
#endif

            beyondOtherPlane = checkBeyondOtherPlane(planeEquation(1:3,i), i)
         endif


         ! Plane i has satisfied all the requirements of being included in the
         !   construction of the Brillouin zone if it is not on another plane
         !   and not beyond another plane. Upon satisfaction of the
         !   requirements we increase the number of facets for the Brillouin
         !   zone and record the index number of i.
         if ((onOtherPlane == 0) .and. (beyondOtherPlane == 0)) then

#ifdef COMPVIS
            write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
            write (61+currBZ,fmt=formatString_i) i-1, "].opacity = 0.20"
            write (61+currBZ,advance="no",fmt="(a)") "bz_point["
            write (61+currBZ,fmt=formatString_i) i-1, "].color = vp.color.white"
#endif
            numBZFacets(currBZ) = numBZFacets(currBZ) + 1
            facetList(numBZFacets(currBZ))%pointID = i
!#ifdef COMPVIS
!         else
!            write (61+currBZ,advance="no",fmt="(a)") "bz_point["
!            write (61+currBZ,fmt=formatString_i) i-1, "].color = vp.color.white"
!            write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
!            write (61+currBZ,fmt=formatString_i) i-1, "].visible = False"
!#endif
         endif
      enddo ! i loop
   end subroutine makeFacetList


   subroutine addFacetVertices(currBZ)

      implicit none

      ! Define passed parameters.
      integer :: currBZ

      ! Define local variables.
      integer :: info
      integer :: g, h, i, j, k, l
      integer, dimension(3) :: pivotIndices
      real (kind=double), dimension(3,3) :: A
      real (kind=double), dimension(3) :: B
#ifdef COMPVIS
      character*8 :: formatString_i, formatString_j, formatString_k
      integer :: s
#endif


      ! Initialize the number of vertices at each facet.
      do i = 1, numBZFacets(currBZ)
         facetList(i)%numVertices = 0
      enddo

      ! Initialize the number of unique vertices.
      numBZVertices(currBZ) = 0

      ! For each facet, we will iterate over all pairs of other facets to
      !   find the set of vertices that belong to the current facet. That is,
      !   we will evaluate every possible facet triplet.
      do i = 1, numBZFacets(currBZ)
#ifdef COMPVIS
         if (facetList(i)%pointID < 11) then
            write (formatString_i,fmt="(a)") "(i0.1,a)"
         else
            write (formatString_i,fmt="(a)") "(i0.2,a)"
         endif
         write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
         write (61+currBZ,fmt=formatString_i) facetList(i)%pointID - 1, &
               & "].visible = True"
         write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
         write (61+currBZ,fmt=formatString_i) facetList(i)%pointID - 1, &
               & "].opacity = 0.2"
#endif
         do j = i+1, numBZFacets(currBZ)
#ifdef COMPVIS
            if (facetList(j)%pointID < 11) then
               write (formatString_j,fmt="(a)") "(i0.1,a)"
            else
               write (formatString_j,fmt="(a)") "(i0.2,a)"
            endif
            write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
            write (61+currBZ,fmt=formatString_j) facetList(j)%pointID - 1, &
                  & "].visible = True"
            write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
            write (61+currBZ,fmt=formatString_j) facetList(j)%pointID - 1, &
                  & "].opacity = 0.2"
#endif
            do k = j+1, numBZFacets(currBZ)
#ifdef COMPVIS
               if (facetList(k)%pointID < 11) then
                  write (formatString_k,fmt="(a)") "(i0.1,a)"
               else
                  write (formatString_k,fmt="(a)") "(i0.2,a)"
               endif
               write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
               write (61+currBZ,fmt=formatString_k) facetList(k)%pointID - 1, &
                     & "].visible = True"
               write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
               write (61+currBZ,fmt=formatString_k) facetList(k)%pointID - 1, &
                     & "].opacity = 0.2"
#endif

               A(1,:) = planeEquation(1:3, facetList(i)%pointID)
               B(1) = planeEquation(4, facetList(i)%pointID)

               A(2,:) = planeEquation(1:3, facetList(j)%pointID)
               B(2) = planeEquation(4, facetList(j)%pointID)

               A(3,:) = planeEquation(1:3, facetList(k)%pointID)
               B(3) = planeEquation(4, facetList(k)%pointID)

               call dgesv(3, 1, A, 3, pivotIndices, B, 3, info)

               ! Ensure that the vertex position does not have a -0.0 or a
               !   very small rounding error.
               do l = 1, 3
                  if (abs(B(l)) < smallThresh) then
                     B(l) = 0.0_double
                  endif
               enddo

               if (info == 0) then

                  ! We found a possible vertex. We now need to make sure that
                  !   it is not outside of any planes. (It is certainly *on*
                  !   three planes because that is how we generated the point.
                  !   However, aside from those three, we need to be sure that
                  !   it is not outside the BZ.)
                  if(checkBeyondOtherPlane(B(:), 0) == 0) then
#ifdef COMPVIS
                     ! Print a sphere at this vertex site.
                     write (61+currBZ,fmt="(a)") "facet_vertex.append("
                     write (61+currBZ,advance="no",fmt="(a)") &
"        vp.sphere(pos=vp.vector("
                     do s = 1, 3
                        write(61+currBZ,advance="no",fmt="(sp,e10.3)") B(s)
                        if (s < 3) then
                           write(61+currBZ,advance="no",fmt="(a)") ", "
                        endif
                     enddo
                     write(61+currBZ,fmt="(a)") "),"
                     write(61+currBZ,advance="no",fmt="(a)") &
"                  radius="
                     write(61+currBZ,advance="no",fmt="(sp,e10.3)") &
                           & sphereRadius / 1.0
                     write(61+currBZ,fmt="(a)") ", color=vp.color.green))"
                     write(61+currBZ,fmt="(a)") "vp.sleep(5.5)"
#endif

                     ! Add this vertex to each of the associated facets, but
                     !   only if it is unique for the facet.
                     do h = 1, 3
                        select case (h)
                           case (1)
                              g = i
                           case (2)
                              g = j
                           case (3)
                              g = k
                        end select

                        ! We will pull a nasty trick to perform the addition
                        !   and ensure that the vertex is unique for this
                        !   facet. We are going to add the vertex to the list
                        !   of vertices for this facet, *BUT* we will only
                        !   increase the number of vertices if this last one
                        !   added is actually unique. We do this so that we
                        !   can create a vertex that can be used in the
                        !   comparison routine and because we will only ever
                        !   print/use the "numVertices" vertices in the list
                        !   for the facet.

                        facetList(g)%vertex( &
                              & facetList(g)%numVertices+1)%coord(:) = B(:)
                        facetList(g)%vertex(facetList(g)%numVertices+1)&
                              & %planeTripleIndex(1) = i
                        facetList(g)%vertex(facetList(g)%numVertices+1)&
                              & %planeTripleIndex(2) = j
                        facetList(g)%vertex(facetList(g)%numVertices+1)&
                              & %planeTripleIndex(3) = k

                        ! Actually, we will increment the number now and then
                        !   if we find any matches we will decrement the number
                        !   and quit.
                        facetList(g)%numVertices = facetList(g)%numVertices + 1
                        do l = 1, facetList(g)%numVertices - 1
                           if (vertexPosEqual(facetList(g)%vertex(l), &
                                 & facetList(g)%vertex(&
                                 & facetlist(g)%numVertices))) then
                              ! We found an equal vertex. Restore the number
                              !   of vertices and quit.
                              facetList(g)%numVertices = &
                                    & facetList(g)%numVertices - 1
                              exit
                           endif
                        enddo ! l
                     enddo ! h

                     ! Add this vertex to the list of vertices only if it is
                     !   unique. We will do the same trick here as we did
                     !   above.
                     numBZVertices(currBZ) = numBZVertices(currBZ) + 1
                     uniqueVertexList(numBZVertices(currBZ))%coord(:) = B(:)
                     uniqueVertexList(numBZVertices(currBZ))%planeTripleIndex=&
                           (/i, j, k/)
                     do l = 1, numBZVertices(currBZ) - 1
                        if (vertexPosEqual(uniqueVertexList(l), &
                              & uniqueVertexList(numBZVertices(currBZ)))) then
                           ! We found an equal vertex. Restor the number of
                           !   vertices and quit.
                           numBZVertices(currBZ) = numBZVertices(currBZ) - 1
                           exit
                        endif
                     enddo
                  endif ! checkBeyondOtherPlane
               endif ! dsegv success
#ifdef COMPVIS
               write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
               write (61+currBZ,fmt=formatString_k) facetList(k)%pointID - 1, &
                     & "].visible = False"
#endif
            enddo ! k
#ifdef COMPVIS
            write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
            write (61+currBZ,fmt=formatString_j) facetList(j)%pointID - 1, &
                  & "].visible = False"
#endif
         enddo ! j
#ifdef COMPVIS
         write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
         write (61+currBZ,fmt=formatString_i) facetList(i)%pointID - 1, &
               & "].visible = False"
#endif
      enddo ! i
   end subroutine addFacetVertices


   ! Check if the plane defined by the lattice index i OR if the point given
   !   by the coordinates coords happen to lie beyond some other plane.
   function checkBeyondOtherPlane(coords, latticeIndex)

      implicit none

      ! Define passed parameters.
      !integer :: i ! Index of the plane i.
      real (kind=double), dimension(3) :: coords
      integer :: latticeIndex

      ! Define the return value.
      integer :: checkBeyondOtherPlane

      ! Define local variables.
      integer :: j, k
      integer :: beyondOtherPlane, oppositeQuadrants
      real (kind=double) :: t, d ! See algebraic description below.
      real (kind=double), dimension(3) :: intersectionPoint
      real (kind=double), dimension(3) :: testPoint
!#ifdef COMPVIS
!      character*8 :: formatString_j
!#endif

      ! This subroutine needs to work in basically the same way regardless of
      !   whether the point we are checking is given implicitly as a plane
      !   index i or as x, y, z coordinates (coords). So, we first regularize
      !   the inputs.
      ! The expectation is that *if* the lattice index is given, then the
      !   test point is the plane equation given by the index number. Note that
      !   this plane equation will then be compared to other plane equations
      !   and so we should skip the comparison when the plane being compared
      !   is the same as the once defined by the lattice index. (See below.)
      !   On the other hand, if the lattice index is 0, then it expected that
      !   the coords will define the point to be compared AND that this point
      !   should be compared against all plane equations.
      if (latticeIndex > 0) then
         testPoint(1:3) = planeEquation(1:3, latticeIndex)
      else
         testPoint(:) = coords(:)
      endif

      beyondOtherPlane = 0
      do j = 1, 26

         ! Do not check against the plane given by the lattice index. Obviously
         !   if the lattice index is given as zero, then all planes are
         !   checked.
         if (j == latticeIndex) cycle

!#ifdef COMPVIS
!         if (j < 11) then
!            write (formatString_j,fmt="(a)") "(i0.1,a)"
!         else
!            write (formatString_j,fmt="(a)") "(i0.2,a)"
!         endif
!         write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
!         write (61+currBZ,fmt=formatString_j) j-1, "].visible = True"
!         write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
!         write (61+currBZ,fmt=formatString_j) j-1, "].color = vp.color.white"
!#endif

         ! The general parametric form for the equation of a line is:
         !   x = x_0 + tu;  y = y_0 + tv;  z = z_0 + tw. Our case is a
         !   special case because the vector <u, v, w> is the same as the
         !   position vector <x_0, y_0, z_0> of the BZ point used to define
         !   plane i. Therefore our line is defined according to:
         !   x = u(1+t);  y = v(1+t);  z = w(1+t).
         ! If we substitute these values for x, y, and z into the
         !   definition for the plane defined by the BZ point j we
         !   can solve for t in terms of u, v, w and l, m, n. <l, m, n>
         !   are the Cartesian coordinates of the position vector used
         !   to define plane j.
         ! The eqn for plane j is: l(x - l) + m(y - m) + n(z - n) = 0.
         !   Upon substitution we have:
         !   l(u(1+t) - l) + m(v(1+t) - m) + n(w(1+t) - n) = 0.
         ! Solving for t yields:
         !   t = (l(l-u) + m(m-v) + n(n-w)) / (lu + mv + nw)
         !   t = (l^2 - lu + m^2 - mv + n^2 - nw) / (lu + mv + nw)
         !   t = (l^2 + m^2 + n^2 - lu - mv - nw) / (lu + mv + nw)
         !   t = (l^2 + m^2 + n^2 - d) / d; d = (lu + mv + nw)
         !   t = ((l^2 + m^2 + n^2) / d) - 1
         ! Given t we can plug in to the above parametric equations for a
         !   line to get the point at which the line intersects the
         !   plane j. If the origin is between this point and the point
         !   that defines plane i then plane j does not block plane i.
         ! On the other hand, if this point and the point that defines
         !   plane i are in the same octant then we compute the distance
         !   from the origin to each point. If this point is closer to
         !   the origin than the point that defines plane i, then plane j
         !   does block plane i and that ensures that plane i cannot
         !   contribute to the Brillouin zone. If this point is farther
         !   from the origin than the point that defines plane i, then
         !   plane j does not block all of plane i and so does not
         !   preclude plane i from contributing to the Brillouin zone.
         !   (At the same time, it also does not ensure that plane i
         !   *does* contribute. That will depend on the other planes
         !   too.)
         d = sum(testPoint(1:3) * planeEquation(1:3,j))
         if (abs(d) < smallThresh) then
            ! Planes are perpendicular, there can be no intersection of a line
            !   from plane defined by testPoint to origin with any point on
            !   plane defined by planeEquation.
            cycle
         endif

         t = (sum(planeEquation(1:3,j)**2) / d) - 1.0_double

         ! The intersection point is a point *on plane j*, but it is
         !   determined by using plane i coordinates. Recall that
         !   <u, v, w> defines the plane i.
         intersectionPoint(:) = testPoint(:) * (1.0_double + t)

         ! Check if the intersectionPoint and the point that defines
         !   plane i are in the same (or opposite) quadrants. We will
         !   assume that the points are in opposite quadrants. If any
         !   Cartesian coordinate component pair has the same sign, then
         !   that is a clear indication that they must be in the same
         !   quadrant.
         oppositeQuadrants = 1
         do k = 1, 3
            if (((intersectionPoint(k) > 0.0_double) .and. &
                  & (testPoint(k) > 0.0_double)) .or. &
                  & ((intersectionPoint(k) < 0.0_double) .and. &
                  & (testPoint(k) < 0.0_double))) then
               oppositeQuadrants = 0
               exit
            endif
         enddo

         ! If the point that defines plane i and the line that goes
         !   through that point and the origin then intersects with
         !   plane j in the opposite quadrant then we can be certain
         !   that plane j will not impede plane i from participating in
         !   the construction of the Brillouin zone.
         if (oppositeQuadrants == 1) then
!#ifdef COMPVIS
!            write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
!            write (61+currBZ,fmt=formatString_j) j-1, "].color = vp.color.red"
!            write (61+currBZ,fmt="(a)") "vp.sleep(0.5)"
!            write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
!            write (61+currBZ,fmt=formatString_j) j-1, "].color = vp.color.white"
!            write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
!            write (61+currBZ,fmt=formatString_j) j-1, "].visible = False"
!#endif
            cycle
         endif

         ! On the other hand, the intersection point and the point that
         !   is used to define plane i are in the same quadrant, then we
         !   need to determine which one is closer to the origin. If the
         !   point used to define plane i is more distant than the
         !   intersection point, then plane i is "beyond the other
         !   plane" and so it cannot contribute to the Brillouin zone.
         if (sum(testPoint(:)**2) > sum(intersectionPoint(:)**2)) then
            beyondOtherPlane = 1
!#ifdef COMPVIS
!            write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
!            write (61+currBZ,fmt=formatString_j) j-1, &
!                  "].color = vp.color.orange"
!            write (61+currBZ,fmt="(a)") "vp.sleep(10)"
!            write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
!            write (61+currBZ,fmt=formatString_j) j-1, "].color = vp.color.white"
!            write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
!            write (61+currBZ,fmt=formatString_j) j-1, "].visible = False"
!#endif
            exit
!#ifdef COMPVIS
!         else
!            write (61+currBZ,fmt="(a)") "vp.sleep(0.1)"
!            write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
!            write (61+currBZ,fmt=formatString_j) j-1, "].visible = False"
!#endif
         endif
      enddo ! j

      ! Store the result to return.
      checkBeyondOtherPlane = beyondOtherPlane

   end function checkBeyondOtherPlane


   subroutine ringSortFacetVertices(currBZ)

      implicit none

      ! Define the passed parameter.
      integer :: currBZ

      ! Define local variables.
      integer :: i, j, k, l
      integer :: nextVertex
      integer, allocatable, dimension(:) :: usedVertices
      real(kind=double) :: argument
      real(kind=double) :: maxAngle
      real(kind=double) :: currentAngle
      real(kind=double), dimension (3) :: centerVector
      real(kind=double), dimension (3) :: vertexVertexVector
      type (facetType) :: tempFacet
!#ifdef COMPVIS
!      write(61+currBZ,fmt="(a)") "facetVec = vp.arrow(shaftwidth=0.05)"
!      write(61+currBZ,fmt="(a)") &
!            & "centerVec = vp.arrow(shaftwidth=0.05, color=vp.color.green)"
!      write(61+currBZ,fmt="(a)") "centerVec.visible = False"
!      write(61+currBZ,fmt="(a)") &
!            & "vvVec = vp.arrow(shaftwidth=0.05, color=vp.color.orange)"
!      write(61+currBZ,fmt="(a)") "vvVec.visible = False"
!#endif

      do i = 1, numBZFacets(currBZ)

!#ifdef COMPVIS
!         write(61+currBZ,advance="no",fmt="(a)") "facetVec.axis=vp.vector("
!         do j = 1, 3
!            write(61+currBZ,advance="no",fmt=("(sp,e10.3)")) &
!                  & planeEquation(j,facetlist(i)%pointID)
!            if (j < 3) then
!               write (61+currBZ,advance="no",fmt="(a)") ", "
!            endif
!         enddo
!         write(61+currBZ,fmt="(a)") ")"
!         write(61+currBZ,fmt="(a)") "vp.sleep(0.5)"
!#endif

         ! Initialize the temp facet that will hold the ring sorted vertices.
         tempFacet%pointID = facetList(i)%pointID
         tempFacet%numVertices = facetList(i)%numVertices
         call copyVertex(facetList(i)%vertex(1), tempFacet%vertex(1))

         ! Allocate space to hold the list of used vertices as we use them.
         !   Because we always assign the first, it is already used and so
         !   we will mark it as such. Also, we will initialize the list of
         !   other used vertices to zero so that any comparison against their
         !   value will pass.
         allocate(usedVertices(facetList(i)%numVertices))
         usedVertices(:) = 0
         usedVertices(1) = 1

         ! Now, we need to add numVertices to the tempFacet. For each vertex
         !   that we want to add we will choose from all vertices in the
         !   current facet (i). The one we will pick to add to the tempFacet
         !   will be the one with the most positive angle between the position
         !   of the most recently added vertex and the center point of the
         !   facet.
         do j = 2, tempFacet%numVertices

            ! Assume that the max angle is zero.
            maxAngle = 0

            ! Compute the vector from the previously added vertex to the
            !   center point of the plane.
            centerVector(:) = planeEquation(1:3,tempFacet%pointID) &
                  & - tempFacet%vertex(j-1)%coord(:)
!#ifdef COMPVIS
!            write(61+currBZ,fmt="(a)") "centerVec.visible = True"
!            write(61+currBZ,advance="no",fmt="(a)") "centerVec.axis=vp.vector("
!            do k = 1, 3
!               write(61+currBZ,advance="no",fmt=("(sp,e10.3)")) centerVector(k)
!               if (k < 3) then
!                  write (61+currBZ,advance="no",fmt="(a)") ", "
!               endif
!            enddo
!            write(61+currBZ,fmt="(a)") ")"
!            write(61+currBZ,advance="no",fmt="(a)") "centerVec.pos=vp.vector("
!            do k = 1, 3
!               write(61+currBZ,advance="no",fmt="(sp,e10.3)") &
!                     & tempFacet%vertex(j-1)%coord(k)
!               if (k < 3) then
!                  write (61+currBZ,advance="no",fmt="(a)") ", "
!               endif
!            enddo
!            write(61+currBZ,fmt="(a)") ")"
!            write(61+currBZ,fmt="(a)") "vp.sleep(2.5)"
!#endif

            ! Iterate through all other vertices of this facet to find the
            !   one with the largest positive angle A-B-C where A is the
            !   point defined by the centerVector, B is the point of the
            !   previously added vertex, and C is each of the other vertices
            !   in this facet.
            kloop: do k = 2, facetList(i)%numVertices
               do l = 1, facetList(i)%numVertices
                  if (usedVertices(l) == k) then
                     cycle kloop
                  endif
               enddo
               vertexVertexVector(:) = facetList(i)%vertex(k)%coord(:) &
                     & - tempFacet%vertex(j-1)%coord(:)
               do l = 1, 3
                  if (abs(vertexVertexVector(l)) < smallThresh) then
                     vertexVertexVector(l) = 0.0_double
                  endif
               enddo
!#ifdef COMPVIS
!               write(61+currBZ,fmt="(a)") "vvVec.visible = True"
!               write(61+currBZ,advance="no",fmt="(a)") "vvVec.axis=vp.vector("
!               do l = 1, 3
!                  write(61+currBZ,advance="no",fmt=("(sp,e10.3)")) &
!                        & vertexVertexVector(l)
!                  if (l < 3) then
!                     write (61+currBZ,advance="no",fmt="(a)") ", "
!                  endif
!               enddo
!               write(61+currBZ,fmt="(a)") ")"
!               write(61+currBZ,advance="no",fmt="(a)") "vvVec.pos=vp.vector("
!               do l = 1, 3
!                  write(61+currBZ,advance="no",fmt="(sp,e10.3)") &
!                        & tempFacet%vertex(j-1)%coord(l)
!                  if (l < 3) then
!                     write (61+currBZ,advance="no",fmt="(a)") ", "
!                  endif
!               enddo
!               write(61+currBZ,fmt="(a)") ")"
!               write(61+currBZ,fmt="(a)") "vp.sleep(2.5)"
!#endif

               ! Theta = acos ((A dot B) / |A| * |B|)
               argument = sum(centerVector(:) * vertexVertexVector(:)) &
                     & / sqrt(sum(centerVector(:)**2)) &
                     & / sqrt(sum(vertexVertexVector(:)**2))
               if (abs((abs(argument) - 1.0_double)) < smallThresh) then
                  argument = sign(1.0_double, argument)
               endif
               currentAngle = acos(argument)

               if (currentAngle > maxAngle) then
                  maxAngle = currentAngle
                  nextVertex = k
               endif
            enddo kloop

            ! After reviewing all other vertices we know the next vertex.
            call copyVertex(facetList(i)%vertex(nextVertex), &
                  & tempFacet%vertex(j))

            ! Store this vertex in the list of used vertices.
            usedVertices(j) = nextVertex
            
         enddo ! j

         ! Free the list of used vertices.
         deallocate(usedVertices)

         ! Now, all the vertices in the tempFacet are copied back over the
         !   vertices in the facetList(i) facet.
         do j = 1, facetList(i)%numVertices
            call copyVertex(tempFacet%vertex(j), facetList(i)%vertex(j))
         enddo

      enddo ! i
   end subroutine ringSortFacetVertices


   subroutine makeUniqueEdgeList (currBZ)

      implicit none

      ! Define passed dummy variables.
      integer :: currBZ

      ! Define local variables.
      integer :: i, j

      ! Traverse the ring list of each facet and add the unique edges to the
      !   list of unique edges.
      do i = 1, numBZFacets(currBZ)

         do j = 1, facetList(i)%numVertices - 1
            call addUniqueEdge(facetList(i)%vertex(j), &
                  & facetList(i)%vertex(j+1), currBZ)
         enddo

         ! Add the last edge from the final vertex back to the first one.
         call addUniqueEdge(facetList(i)%vertex(facetList(i)%numVertices), &
               & facetList(i)%vertex(1), currBZ)
      enddo
   end subroutine makeUniqueEdgeList


   subroutine addUniqueEdge (vertex1, vertex2, currBZ)

      implicit none

      ! Define passed dummy parameters.
      integer :: currBZ
      type (vertexType) :: vertex1
      type (vertexType) :: vertex2

      ! Define local variables.
      integer :: i, j
      integer :: edgeAlreadyPresent

      ! Traverse the current list of unique edges to see if the proposed
      !   (vertex1, vertex2) pair is already listed. If so, then don't do
      !   anything. If not, then add it!

      ! Assume that the proposed edge is not present.
      edgeAlreadyPresent = 0

      do i = 1, numBZEdges(currBZ)

         ! Check if the proposed edge matches this edge.
         if ((vertexPosEqual(vertex1, uniqueEdgeList(i)%vertex(1)) .and. &
               (vertexPosEqual(vertex2, uniqueEdgeList(i)%vertex(2)))) .or. &
               (vertexPosEqual(vertex1, uniqueEdgeList(i)%vertex(2)) .and. &
               (vertexPosEqual(vertex2, uniqueEdgeList(i)%vertex(1))))) then
            edgeAlreadyPresent = 1
            exit
         endif
      enddo

      if (edgeAlreadyPresent == 0) then
         numBZEdges(currBZ) = numBZEdges(currBZ) + 1
         call copyVertex(vertex1, uniqueEdgeList(numBZEdges(currBZ))%vertex(1))
         call copyVertex(vertex2, uniqueEdgeList(numBZEdges(currBZ))%vertex(2))
      else
         return
      endif
   end subroutine addUniqueEdge


   ! Take each plane (defined as being perpendicular to the members of the set
   !   of neighboring recipical lattice points) and compute the equation that
   !   defines it. Note: given a Cartesian coordinate vector with magnitudes
   !   (l, m, n) in the x, y, z directions, the equation for the plane that
   !   goes through the point (l, m, n) is: l*(x-l) + m*(y-m) + n*(z-n) = 0.
   !   That is: l*x + m*y + n*z - l^2 - m^2 - n^2 = 0 or also expressed as:
   !   l*x + m*y + n*z = d with d = l^2 + m^2 + n^2. To "store" the equation
   !   for each plane we will simply hold the Cartesian <l, m, n> vector and
   !   the value of d.
   subroutine makePlaneEquations (currBZ)

      implicit none

      ! Define passed parameters.
      integer :: currBZ

      ! Define local variables.
      integer :: i, j, k, s, planeIndex
#ifdef COMPVIS
      character*8 :: formatString
#endif

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
                  ! Eliminate machine rounding errors and any possible -0.0s.
                  if (abs(latticePointXYZ(s,planeIndex)) < smallThresh) then
                     latticePointXYZ(s,planeIndex) = 0.0_double
                  endif
               enddo

               ! Store them and the solution (constant) term. Note that the
               !   0.5 factor is applied because the planes are formed at the
               !   midpoint between the lattice point and the origin.
               planeEquation(1:3,planeIndex) = &
                     & latticePointXYZ(1:3,planeIndex) * 0.5_double
               planeEquation(4,planeIndex) = &
                     & sum((latticePointXYZ(:,planeIndex) * 0.5_double)**2)

#ifdef COMPVIS
               ! Prepare a format string for array indices.
               if (planeIndex < 11) then
                  write (formatString,fmt="(a)") "(i0.1,a)"
               else
                  write (formatString,fmt="(a)") "(i0.2,a)"
               endif

               ! Create, pause, and then hide the recip lattice planes and the
               !   Brillouin zone planes.
               write (61+currBZ,fmt="(a)") "recip_plane.append("
               write (61+currBZ,advance="no",fmt="(a)") &
"        vp.cylinder(pos=vp.vector("
               do s = 1, 3
                  write (61+currBZ,advance="no",fmt="(sp,e10.3)") &
                        latticePointXYZ(s,planeIndex)
                  if (s < 3) then
                     write (61+currBZ,advance="no",fmt="(a)") ", "
                  endif
               enddo
               write (61+currBZ,fmt="(a)") "),"
               write (61+currBZ,advance="no",fmt="(a)") &
"                    axis=vp.vector("
               do s = 1, 3
                  write (61+currBZ,advance="no",fmt="(sp,e10.3)") &
                        latticePointXYZ(s,planeIndex) + planeThickness
                  if (s < 3) then
                     write (61+currBZ,advance="no",fmt="(a)") ", "
                  endif
               enddo
               write (61+currBZ,fmt="(a)") "),"
               write (61+currBZ,advance="no",fmt="(a)") &
"                    radius="
               write (61+currBZ,advance="no",fmt="(sp,e10.3)") planeExtent
               write (61+currBZ,fmt="(a)") ", length=0.003))"
               write (61+currBZ,fmt="(a)") "bz_plane.append("
               write (61+currBZ,advance="no",fmt="(a)") &
"        vp.cylinder(pos=vp.vector("
               do s = 1, 3
                  write (61+currBZ,advance="no",fmt="(sp,e10.3)") &
                        planeEquation(s,planeIndex)
                  if (s < 3) then
                     write (61+currBZ,advance="no",fmt="(a)") ", "
                  endif
               enddo
               write (61+currBZ,fmt="(a)") "),"
               write (61+currBZ,advance="no",fmt="(a)") &
"                    axis=vp.vector("
               do s = 1, 3
                  write (61+currBZ,advance="no",fmt="(sp,e10.3)") &
                        planeEquation(s,planeIndex) + planeThickness
                  if (s < 3) then
                     write (61+currBZ,advance="no",fmt="(a)") ", "
                  endif
               enddo
               write (61+currBZ,fmt="(a)") "),"
               write (61+currBZ,advance="no",fmt="(a)") &
"                    radius="
               write (61+currBZ,advance="no",fmt="(sp,e10.3)") planeExtent
               write (61+currBZ,fmt="(a)") ", length=0.003))"
               write (61+currBZ,advance="no",fmt="(a)") "recip_point["
               write (61+currBZ,fmt=formatString) planeIndex-1, &
"].color = vp.color.red"
               !write (61+currBZ,fmt="(a)") "vp.sleep(3)"
               write (61+currBZ,advance="no",fmt="(a)") "bz_point["
               write (61+currBZ,fmt=formatString) planeIndex-1, &
"].color = vp.color.white"
               write (61+currBZ,advance="no",fmt="(a)") "recip_plane["
               write (61+currBZ,fmt=formatString) planeIndex-1, &
"].visible = False"
               write (61+currBZ,advance="no",fmt="(a)") "bz_plane["
               write (61+currBZ,fmt=formatString) planeIndex-1, &
"].visible = False"
#endif
            enddo
         enddo
      enddo

   end subroutine makePlaneEquations

   ! Compare the positions of two vertices and declare them equal even if the
   !   vertices are defined using different planeTripleIndices. (Because not
   !   all vertices need to be compared that way.)
   function vertexPosEqual(vertex1, vertex2)
 
      implicit none
 
      ! Define dummy parameters.
      type (vertexType) :: vertex1, vertex2
 
      ! Define the return variable.
      logical :: vertexPosEqual
 
      ! Define local variables.
      integer :: i, j
      real (kind=double) :: threshold
 
      ! Assume that the two vertices are equal and define the threshold.
      vertexPosEqual = .true.
      threshold = 0.000001
 
      do i = 1, 3
         if (abs(vertex1%coord(i) - vertex2%coord(i)) > threshold) then
            vertexPosEqual = .false.
            return
         endif
      enddo
 
   end function vertexPosEqual


   ! Compare two vertices according to both their positions and the values of
   !   their planeTripleIndices.
   function verticesEqual(vertex1, vertex2)
 
      implicit none
 
      ! Define dummy parameters.
      type (vertexType) :: vertex1, vertex2
 
      ! Define the return variable.
      logical :: verticesEqual
 
      ! Define local variables.
      integer :: i, j
      real (kind=double) :: threshold
 
      ! Assume that the two vertices are equal and define the threshold.
      verticesEqual = .true.
      threshold = 0.000001
 
      do i = 1, 3
         if (abs(vertex1%coord(i) - vertex2%coord(i)) > threshold) then
            verticesEqual = .false.
            return
         endif
 
         if (vertex1%planeTripleIndex(i) /= vertex2%planeTripleIndex(i)) then
            verticesEqual = .false.
            return
         endif
      enddo
 
   end function verticesEqual


   subroutine printBrillouinZones (maxBZ)

      ! Use necessary modules.
      use CommandLine_O

      ! Make sure no funny variables are defined.
      implicit none

      ! Define dummy variables.
      integer :: maxBZ ! Highest BZ we will print.

      ! Define local variables.
      integer :: i, j, k, l, m ! Loop indices.
      real (kind=double), dimension(3) :: fractABC
      character*25 :: filename

      ! Step through each of the Brillouin zones we will make.
      do i = 1, maxBZ

         ! Print all of the facets.
         write (51+i,fmt="(a9)") "faces = ["
         do j = 1, numBZFacets(i)
            write (51+i,advance="NO",fmt="(a1)") "["
            do k = 1, facetList(j)%numVertices
               facetList(j)%vertex(k)%coord(:) = &
                     & facetList(j)%vertex(k)%coord(:) * scaleFactor
               write (51+i,advance="NO",fmt="(a1)") "["
               write (51+i,advance="NO",fmt="(f16.12, a2)") &
                     & facetList(j)%vertex(k)%coord(1), ", "
               write (51+i,advance="NO",fmt="(f16.12, a2)") &
                     & facetList(j)%vertex(k)%coord(2), ", "
               write (51+i,advance="NO",fmt="(f16.12, a3)") &
                     & facetList(j)%vertex(k)%coord(3), "], "
            enddo ! k

            ! Write the first vertex again to close the polygon
            write (51+i,advance="NO",fmt="(a1)") "["
            write (51+i,advance="NO",fmt="(f16.12, a2)") &
                  & facetList(j)%vertex(1)%coord(1), ", "
            write (51+i,advance="NO",fmt="(f16.12, a2)") &
                  & facetList(j)%vertex(1)%coord(2), ", "
            write (51+i,advance="NO",fmt="(f16.12, a1)") &
                  & facetList(j)%vertex(1)%coord(3), "]"

            if (j < numBZFacets(i)) then
               write (51+i,fmt="(a2)") "],"
            else
               write (51+i,fmt="(a2)") "]]"
            endif
         enddo ! j numBZFacets(i)

         ! Print all of the edges.
         write (51+i,fmt="(a9)") "edges = ["
         do j = 1, numBZEdges(i)
            uniqueEdgeList(j)%vertex(1)%coord(:) = &
                  & uniqueEdgeList(j)%vertex(1)%coord(:) * scaleFactor
            uniqueEdgeList(j)%vertex(2)%coord(:) = &
                  & uniqueEdgeList(j)%vertex(2)%coord(:) * scaleFactor

            write (51+i,advance="NO",fmt="(a2)") "[["
            write (51+i, advance="NO",fmt="(f16.12, a2)") &
                  & uniqueEdgeList(j)%vertex(1)%coord(1), ", "
            write (51+i, advance="NO",fmt="(f16.12, a2)") &
                  & uniqueEdgeList(j)%vertex(1)%coord(2), ", "
            write (51+i, advance="NO",fmt="(f16.12, a3)") &
                  & uniqueEdgeList(j)%vertex(1)%coord(3), "], "
            write (51+i,advance="NO",fmt="(a1)") "["
            write (51+i, advance="NO",fmt="(f16.12, a2)") &
                  & uniqueEdgeList(j)%vertex(2)%coord(1), ", "
            write (51+i, advance="NO",fmt="(f16.12, a2)") &
                  & uniqueEdgeList(j)%vertex(2)%coord(2), ", "
            write (51+i, advance="NO",fmt="(f16.12)") &
                  & uniqueEdgeList(j)%vertex(2)%coord(3)

            if (j < numBZEdges(i)) then
               write (51+i,fmt="(a4)") "]], "
            else
               write (51+i,fmt="(a3)") "]]]"
            endif
         enddo ! j numBZEdges(i)

         ! Print all of the vertices
         write (51+i,fmt="(a12)") "vertices = ["
         do j = 1, numBZVertices(i)
            uniqueVertexList(j)%coord(:) = &
                  & uniqueVertexList(j)%coord(:) * scaleFactor
            write (51+i,advance="NO",fmt="(a1)") "["
            write (51+i,advance="NO",fmt="(f16.12, a2)") &
                  & uniqueVertexList(j)%coord(1), ", "
            write (51+i,advance="NO",fmt="(f16.12, a2)") &
                  & uniqueVertexList(j)%coord(2), ", "
            write (51+i,advance="NO",fmt="(f16.12)") &
                  & uniqueVertexList(j)%coord(3)

            if (j < numBZVertices(i)) then
               write (51+i,fmt="(a3)") "], "
            else
               write (51+i,fmt="(a2)") "]]"
            endif
         enddo ! j numBZEdges(i)
      enddo ! i maxBZ

   end subroutine printBrillouinZones

end module BrillouinZones_O


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
   subroutine printKPoints (fileUnit, maxBZ, recipLattice)

      use CommandLine_O

      implicit none

      ! Define dummy parameters.
      integer, intent (in) :: fileUnit
      integer :: maxBZ
      real (kind=double), dimension(3,3) :: recipLattice

      ! Define local variables.
      integer :: i, j, k, l
      real (kind=double), dimension(3) :: xyzKPoint

      write (fileUnit,fmt="(a17)") "NUM_BLOCH_VECTORS"
      write (fileUnit,fmt="(i9.9)") numFoldedKPoints

      write (fileUnit,fmt="(a19)") "NUM_WEIGHT_KA_KB_KC"
      do i = 1, numFoldedKPoints
         write (fileUnit,fmt="(i5,2x,f20.16,3(2x,f12.8))") i,kPointWeight(i),&
               & abcFoldedKPoints(:,i)
      enddo

      ! Print the mesh point positions, the folded point positions, and then
      !   the weights.
      if (maxBZ > 0) then
         do i = 1, maxBZ
            write (fileUnit+i,fmt="(a16)") "mesh_kpoints = ["
            do j = 1, numMeshKPoints

               ! Convert the abc coordinate to xyz.
               do k = 1, 3 ! xyz axes

                  ! Initialize the x, y, or z coordinate.
                  xyzKPoint(k) = 0.0_double

                  do l = 1, 3 ! abc axes
                     
                     xyzKPoint(k) = xyzKPoint(k) &
                           + abcMeshKPoints(l,j) * recipLattice(k,l) &
                           * scaleFactor
                  enddo
               enddo

               write (fileUnit+i,advance="no",fmt="(a1)") "["
               write (fileUnit+i,advance="no",fmt="(f16.12, a2)") &
                     & xyzKPoint(1), ", "
               write (fileUnit+i,advance="no",fmt="(f16.12, a2)") &
                     & xyzKPoint(2), ", "
               write (fileUnit+i,advance="no",fmt="(f16.12, a1)") &
                     & xyzKPoint(3), "]"

               if (j < numMeshKPoints) then
                  write (fileUnit+i,fmt="(a2)") ", "
               else
                  write (fileUnit+i,fmt="(a1)") "]"
               endif
            enddo ! j

            write (fileUnit+i,fmt="(a18)") "folded_kpoints = ["
            do j = 1, numFoldedKPoints

               ! Convert the abc coordinate to xyz.
               do k = 1, 3 ! xyz axes

                  ! Initialize the x, y, or z coordinate.
                  xyzKPoint(k) = 0.0_double

                  do l = 1, 3 ! abc axes
                     
                     xyzKPoint(k) = xyzKPoint(k) &
                           + abcFoldedKPoints(l,j) * recipLattice(k,l) &
                           * scaleFactor
                  enddo
               enddo

               write (fileUnit+i,advance="no",fmt="(a1)") "["
               write (fileUnit+i,advance="no",fmt="(f16.12, a2)") &
                     & xyzKPoint(1), ", "
               write (fileUnit+i,advance="no",fmt="(f16.12, a2)") &
                     & xyzKPoint(2), ", "
               write (fileUnit+i,advance="no",fmt="(f16.12, a1)") &
                     & xyzKPoint(3), "]"

               if (j < numFoldedKPoints) then
                  write (fileUnit+i,fmt="(a2)") ", "
               else
                  write (fileUnit+i,fmt="(a1)") "]"
               endif
            enddo ! j

            ! Write the kpoint weights.
            write (fileUnit+i,fmt="(a18)") "kpoint_weights = ["
            do j = 1, numFoldedKPoints - 1
               write (fileUnit+i,fmt="(a1, f16.12, a2)") "[", &
                     & kPointWeight(j), "],"
            enddo ! j
            write (fileUnit+i,fmt="(a1, f16.12, a2)") "[", &
                  & kPointWeight(j), "]]"
         enddo ! i
      endif
   end subroutine printKPoints

end module KPointMesh_O


!module PrintRecipModel_O
!
!   use CommandLine_O
!
!   implicit none
!
!   contains
!
!   subroutine printVertices
!   end subroutine printVertices
!
!   subroutine printEdges
!   end subroutine printEdges
!
!   subroutine printFaces (maxBZ, numBZFacets, facetList)
!
!      implicit none
!
!      ! Define passed parameters.
!      integer :: maxBZ
!      integer, dimension(:) :: numBZFacets
!      type (facetType), 
!
!   end subroutine printFaces
!
!end module PrintRecipModel_O


program makekpoints

   ! This program is used to generate a set of kpoints for a given crystal and
   !   kpoint mesh request.

   ! Import necessary parameter modules.
   use O_Kinds
   use O_Constants

   ! Import necessary modules.
   use PointGroupOperations_O
   use Lattice_O
   use BrillouinZones_O
   use KPointMesh_O
   use CommandLine_O
#ifdef COMPVIS
   use VPython_O
#endif


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

   ! If this program is called multiple times (which would usually happen
   !   under many common conditions when using makeinput) then this will
   !   prevent recalculation of the supplementary Brillouin zone files.
   if (ioerr /= 0) then
      doBrillouinZone = 0
   endif

#ifdef COMPVIS
   ioerr = 0
   do i = 1, doBrillouinZone
      write (filename,fmt="(a3,i1,a3)") "BZ.", i, ".py"
      open (unit=61+i, file=filename, form='formatted', status='new', &
            iostat=ioerr)
      if (ioerr /= 0) then
         exit
      endif
   enddo

   ! Similar to the above except for the computational visualization. Also,
   !   it is assumed that this program will never be called for any practical
   !   purpose beyond computational visualization if COMPVIS is defined. So,
   !   on any subsequent call to this program (so ioerr == 0) we can return.
   if (ioerr /= 0) then
      stop
   endif
#endif


   ! Read the lattice and compute all of its information.
   call readRealLattice(50)
   call computeLatticeData
   if (doBrillouinZone > 0) then
      call printLatticeData(doBrillouinZone) ! Print if requested
   endif


#ifdef COMPVIS
   ! Create the header for computational visualization if requested.
   call makeVPythonHeader(doBrillouinZone)

   ! Create reciprocal lattice for computational visualization if requested.
   call makeVPythonRecipLattice(doBrillouinZone)
#endif


   ! Read the point group operations and compute their values on the given
   ! a,b,c lattice.
   call readPointOps(50)
   call computeABCRecipPointOps(realLattice,recipLattice)

   ! Read the kpoint mesh parameters, initialize the mesh, fold the mesh, and
   !   print the results.  It should be generally understood that the
   !   variables with abc in the name refer to the reciprocal lattice
   !   abcvectors and not the real space lattice.  This is only given
   !   explictly for the variable abcRecipPointOps though.
   call readMeshParameters(50)
   call initMesh
   call foldMesh (numPointOps,abcRecipPointOps)
   call printKPoints(51,doBrillouinZone,recipLattice)

   ! In the case that the user wants information about a the Brillouin zone,
   !   then we compute and print it.
   if (doBrillouinZone > 0) then
      call computeBrillouinZones(doBrillouinZone)
      call printBrillouinZones(doBrillouinZone)
   endif

   ! Clean up any data allocations.
   call cleanUp

   close (50)
   close (51)
   do i = 1, doBrillouinZone
      close (51+i)
#ifdef COMPVIS
      close (61+i)
#endif
   enddo

   contains

   subroutine cleanUp

      ! Import necessary modules.
      use Lattice_O

      if (doBrillouinZone > 0) then
         deallocate (numBZVertices)
         deallocate (numBZEdges)
         deallocate (numBZFacets)
      endif
   end subroutine cleanUp

end program makekpoints
