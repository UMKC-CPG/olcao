module O_KPoints

   ! Import necessary modules.
   use O_Kinds
   use O_Constants

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer :: kPointStyleCode ! In all cases the kpoint information is read
         !   from a separate file than the olcao.dat file. This separates the
         !   atomic description (basis and potential) in the olcao.dat file
         !   and the kpoint mesh in a separate file. The code value in that
         !   file will cause the kpoints to be set up in different ways.
         !   If 0 then the kpoints will be read as an explicit list.
         !   If 1 then the kpoints will be read as a mesh definition along
         !      with a specific shift.
         !   If 2 then the kpoints will be read as a minimum density that
         !      will be applied along each axis. Additionally a specific
         !      shift is also read in.
   integer :: kPointIntgCode ! For each kpoint mesh it may be possible to use
         !   a variety of different integration methods from simple histograms
         !   to linear analytic tetrahedron (LAT) integration, etc. There are
         !   some constraints to be aware of. If the kpoints are given in the
         !   form of a regular grid (not necessarily square) then any LAT
         !   method can be used. A regular grid will certainly be present for
         !   kPointStyleCode values of 1 and 2, but for kPointStyleCode 0 it
         !   is possible that the mesh will not be regular. It is the
         !   responsibilty of the user to know whether or not a LAT method can
         !   be used if they provide the kPoints as an explicit list with
         !   kPointStyleCode 0.
         !   If kPointIntgCode=0 then the histogram method will be used.
         !   If kPointIntgCode=1 then the LAT method will be used.
   integer :: numKPoints ! The number of kpoints in the system.
   integer :: numPaths ! The number of discontinuous high-sym. KP paths.
   integer :: numPathKP ! Number of kpoints that will be used to create the
         !   path between the high symmetry kpoints.
   integer :: isCartesian ! 1 = yes, 0 = no (Are the high symmetry kpoints
         !   given in cartesian coordinates?)
   integer :: numTotalHighSymKP ! Total number of high symmetry kpoints over
         !   all paths.
   integer, dimension (dim3) :: numAxialKPoints ! Number of kpoints along each
         !   of the a, b, c axes defining the kpoint mesh. This is only used
         !   if the kPointStyleCode equals 1 or 2.
   integer, dimension (dim3) :: numAxialMTOP_KP ! Number of kpoints along
         !   each of the a,b,c axes defining the kpoint mesh. This is only
         !   used if either doMTOP_SCF or doMTOP_PSCF flag is set to 1.
   integer, allocatable, dimension (:) :: numHighSymKP ! Number of high
         !   symmetry kpoints that define the vertices of each path to be taken
         !   for the band diagram.
   integer, allocatable, dimension (:,:) :: mtopKPMap ! Map between the
         !   desired sequence of kpoints for MTOP calculations and the one-
         !   dimensional list of kpoints.
   real (kind=double) :: minKPointDensity ! The minimum number of kpoints per
         !   unit reciprocal-space volume (Bohr^-3). The total kpoint count is
         !   at least minKPointDensity * recipCellVolume, distributed as
         !   uniformly as possible across the three axes. Only used when
         !   kPointStyleCode equals 2.
   real (kind=double), dimension (dim3) :: kPointShift ! Fractional amount to
         !   shift the kpoint mesh by along each of the a,b,c axes. This is
         !   only used if the kPointStyleCode equals 1 or 2.
   real (kind=double), dimension (dim3) :: kPointShiftMTOP ! Fractional amount
         !   to shift the kpoint mesh by along each of the a,b,c axes. This is
         !   only used if either doMTOP_SCF or doMTOP_PSCF flag is set to 1.
   real (kind=double), allocatable, dimension (:) :: kPointWeight ! The
         !   weight assigned to each kpoint.
   real (kind=double), allocatable, dimension (:,:) :: kPoints ! The acutal
         !   kpoints.  The first dimension holds the three cartesian xyz
         !   coordinates.  The second dimension is the index over the number
         !   of kpoints.
   real (kind=double), allocatable, dimension (:,:,:) :: highSymKP ! These
         !   are the high symmetry kpoints of the requested path.  The first
         !   dimension holds the three space coordinates of the kpoints
         !   to be indexed by the second dimension (numHighSymKP). The last
         !   dimension indexes which discrete path these points define.
   real (kind=double), allocatable, dimension (:) :: pathKPointMag ! Distance
         !   magnitudes of kpoints on the path for band diagrams.
   complex (kind=double), allocatable, dimension (:,:) :: phaseFactor ! The
         !   cosine and sine of the dot product between each kpoint and each
         !   real space superlattice vector.
   character*1, allocatable, dimension (:,:) :: highSymKPChar ! The characters
         !   that identify each point on the high symmetry kpoint paths.
   integer :: numPointOps ! Number of point group operations for IBZ
         !   symmetry reduction. These are read from the kpoint input file
         !   when kPointStyleCode equals 2. They are the rotational parts of
         !   the space group operations (translations stripped away) and are
         !   given in fractional (abc) coordinates.
   real (kind=double), allocatable, dimension (:,:,:) :: abcPointOps
         !   The 3x3 rotation matrices for each point group operation, stored
         !   in fractional (abc) coordinates. Dimensions: (3,3,numPointOps).
   real (kind=double), allocatable, dimension (:,:,:) :: abcRecipPointOps
         !   The point group operations converted to reciprocal-space abc
         !   coordinates using the real and reciprocal lattice vectors.
         !   Computed by computeRecipPointOps. Dims: (3,3,numPointOps).
   integer :: numFullMeshKP ! The total number of kpoints in the full
         !   uniform mesh BEFORE IBZ reduction. This is the product of
         !   numAxialKPoints(1) * numAxialKPoints(2) * numAxialKPoints(3).
         !   Only set when the mesh is built internally (style codes 1, 2).
   integer, allocatable, dimension (:) :: fullToIBZMap
         !   Maps each full-mesh kpoint index (1..numFullMeshKP) to its
         !   IBZ representative kpoint index (1..numKPoints). Used by
         !   LAT integration to unfold eigenvalues from the IBZ to the
         !   full mesh: eigenValues(band, fullToIBZMap(k), spin).
         !   Only allocated when IBZ reduction is performed.
   integer :: numTetrahedra ! The total number of tetrahedra tiling the
         !   reciprocal cell. Equal to 6 * nA * nB * nC.
   real (kind=double) :: tetraVol ! The BZ fraction
         !   for each tetrahedron: 1 / numTetrahedra.
   integer, allocatable, dimension (:,:) :: tetrahedra
         !   The four corner kpoint indices for each tetrahedron.
         !   Dimensions: (4, numTetrahedra). Indices reference the full
         !   uniform mesh (1..numFullMeshKP), NOT the IBZ-reduced list.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine readKPoints(readUnit, writeUnit)

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3

   ! Import necessary subroutine modules.
   use O_ReadDataSubs

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from
                                        ! which we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Define the local variables used in this subroutine.
   integer :: i ! Loop index variable.
   integer :: counter ! Dummy variable to read the kpoint index number.

   ! Read the kpoint definition style code. I.e., how the mesh should be
   !   constructed.
   call readData(readUnit,writeUnit,kPointStyleCode,len('KPOINT_STYLE_CODE'),&
         & 'KPOINT_STYLE_CODE')

   ! Read the code that defines the type of integration method to be used.
   !   I.e., Weighted sums, tetrahedral, etc.
   call readData(readUnit,writeUnit,kPointIntgCode,len('KPOINT_INTG_CODE'),&
         & 'KPOINT_INTG_CODE')

   ! Depending on kPointStyleCode we will read the associated relevant data.
   if (kPointStyleCode == 0) then ! Read an explicit list of kpoints.

      ! Read the number of kpoints.
      call readData(readUnit,writeUnit,numKPoints,len('NUM_BLOCH_VECTORS'),&
            & 'NUM_BLOCH_VECTORS')

      ! Read the number of kpoints requested on each axis.
      call readData(readUnit,writeUnit,3,numAxialKPoints,&
            & len('NUM_AXIAL_KPOINTS'),'NUM_AXIAL_KPOINTS')

      ! Allocate space to hold the kpoints and their weighting factors.
      allocate (kPoints(dim3,numKPoints))
      allocate (kPointWeight(numKPoints))

      ! Read past the 'NUM_WEIGHT_KA_KB_KC' header.
      call readAndCheckLabel(readUnit,writeUnit,len('NUM_WEIGHT_KA_KB_KC'),&
            & 'NUM_WEIGHT_KA_KB_KC')

      ! Note that the values read in here are in terms of a,b,c fractional
      !   coordinates.  Later, after the reciprocal lattice has been
      !   initialized these values will be changed into x,y,z cartesian
      !   coordinates.
      do i = 1, numKPoints
         read (readUnit,*) counter, kPointWeight(i), kPoints(1:dim3,i)
         write (writeUnit,fmt="(i5,1x,4f15.8)") counter, kPointWeight(i), &
               & kPoints(1:dim3,i)
      enddo
      call flush (writeUnit)

   elseif (KPointStyleCode == 1) then
      ! Read axial numbers of kpoints and a shift.

      ! Read the number of kpoints along each a,b,c axis.
      call readData(readUnit,writeUnit,3,numAxialKPoints,len('NUM_KP_A_B_C'),&
            & 'NUM_KP_A_B_C')

      ! Read the fractional distance that the uniform mesh should be shifted
      !   away from the origin along each a,b,c axis.
      call readData(readUnit,writeUnit,3,kPointShift,len('KP_SHIFT_A_B_C'),&
            & 'KP_SHIFT_A_B_C')

      ! The kpoint mesh will be constructed later, once implicit information
      !   such as the reciprocal lattice has been obtained.

   elseif (KPointStyleCode == 2) then
      ! Read density-based kpoint specification: a minimum volume density,
      !   a shift, and the point group operations needed for IBZ symmetry
      !   reduction.

      ! Read the minimum volume density of kpoints (kpoints per unit reciprocal-
      !   space volume). Note: the file label says LINE_DENSITY for historical
      !   reasons but the value is actually a volume density.
      call readData(readUnit,writeUnit,minKPointDensity,&
            & len('MIN_KP_LINE_DENSITY'),&
            & 'MIN_KP_LINE_DENSITY')

      ! Read the fractional shift for the mesh along each
      !   a,b,c axis.
      call readData(readUnit,writeUnit,3,kPointShift,&
            & len('KP_SHIFT_A_B_C'),&
            & 'KP_SHIFT_A_B_C')

      ! Read the number of point group operations for IBZ reduction. These
      !   are the rotational parts of the space group symmetry operations
      !   (translations stripped) and are given in fractional (abc) coords.
      call readData(readUnit,writeUnit,numPointOps,&
            & len('NUM_POINT_OPS'),'NUM_POINT_OPS')

      ! Read the point group operation matrices. Each is a 3x3 rotation matrix
      !   in abc coordinates. The file may contain blank lines between
      !   operations for readability, so we read one row (column of the
      !   Fortran array) per line and skip blank lines between operations.
      allocate (abcPointOps(3,3,numPointOps))
      call readAndCheckLabel(readUnit,writeUnit,&
            & len('POINT_OPS'),'POINT_OPS')
      do i = 1, numPointOps
         ! Skip any blank line before this operation. The first operation
         !   has no blank line before it, but subsequent ones do.
         if (i > 1) read (readUnit,*)
         ! Read each row of the matrix from one file line. Each line fills
         !   one column of the Fortran array (column-major storage).
         read (readUnit,*) abcPointOps(:,1,i)
         read (readUnit,*) abcPointOps(:,2,i)
         read (readUnit,*) abcPointOps(:,3,i)
         do counter = 1, 3
            write (writeUnit,fmt="(3f14.8)") &
                  & abcPointOps(:,counter,i)
         enddo
         write (writeUnit,*)
      enddo
      call flush (writeUnit)

      ! The kpoint mesh will be constructed later, once implicit information
      !   such as the reciprocal lattice has been obtained.
   endif


   ! Note: tetrahedra generation for LAT (kPointIntgCode==1) is deferred to
   !   initializeKPoints, after the mesh is built and numAxialKPoints is
   !   fully determined.

end subroutine readKPoints


subroutine readSYBDKPoints(readUnit, writeUnit)

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3
   use O_ReadDataSubs
   use O_WriteDataSubs

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Define the local variables used in this subroutine.
   integer :: i, j ! Loop index variable.

   ! A band structure is often plotted as a contiguous path of points in
   !   k-space connecting specific high-symmetry k-points. However, it is not
   !   a requirement. If one so chooses, the band structure can be plotted
   !   using a series of paths that do not have to be contiguous with
   !   each other. See Setyawan and Curtarolo, "High-throughput electronic
   !   band structure calculations: Challenges and tools", Comp. Mat. Sci.,
   !   vol. 49, pages 299-312, (2010).

   ! Hence, we first read in the number of distcontinuous paths to plot, the
   !   number of kpoints to use across all paths, and whether or not the set
   !   of all kpoints are given in Cartesian coordinates or not.
   call readData(readUnit,writeUnit,numPaths,numPathKP,isCartesian,&
                    len('SYBD_INPUT_DATA'),'SYBD_INPUT_DATA')
   write (writeUnit,*) 'Number of discrete paths         = ',numPaths
   write (writeUnit,*) 'Number of path K-Points          = ',numPathKP
   write (writeUnit,*) '1 = cart, 0 = fract:             = ',isCartesian
   call flush (writeUnit)

   ! Allocate space to hold the number of high symmetry kpoints that will
   !   exist for each path. E.g., there may be 6 kpoints that define the first
   !   path and only 2 kpoints that define the second path.
   allocate (numHighSymKP(numPaths))

   ! Read the number of high symmetry kpoints that define the vertices for
   !   each path.
   call readData(readUnit,writeUnit,numPaths,numHighSymKP(:),0,'')

   ! Compute the total number of high symmetry kpoints over all paths.
   numTotalHighSymKP = sum(numHighSymKP)

   write (writeUnit,*) 'Number of high symmetry K-Points = ',numTotalHighSymKP
   call flush (writeUnit)
   
   ! Allocate space to hold the high symmetry kpoints and their identifying
   !   characters.
   allocate (highSymKP(dim3,maxval(numHighSymKP),numPaths))
   allocate (highSymKPChar(maxval(numHighSymKP),numPaths))

   ! Read the high symmetry kpoints for each path.
   do i = 1, numPaths
   
      ! Read the coordinates of the high symmetry kpoints.
      do j = 1, numHighSymKP(i)
         call readData(readUnit,writeUnit,3,highSymKP(:,j,i),0,'')
      enddo
   enddo

   ! Read the set of characters identifying the high symmetry kpoints.
   do i = 1, numPaths
      do j = 1, numHighSymKP(i)
         call readData(readUnit,writeUnit,highSymKPChar(j,i),0)
      enddo
   enddo
   call writeData(writeUnit)

end subroutine readSYBDKPoints


subroutine readMTOPKPoints(readUnit, writeUnit)

   ! Include the necessary modules
   use O_Kinds
   use O_Constants, only: dim3
   use O_ReadDataSubs

   implicit none

   ! Passed parameters
   integer, intent(in) :: readUnit   ! The unit number of the file to read
   integer, intent(in) :: writeUnit  ! The unit number of the file to write

   ! Read input data
   call readAndCheckLabel(readUnit, writeUnit, len('MTOP_INPUT_DATA'), &
         & 'MTOP_INPUT_DATA')
   call readData(readUnit, writeUnit, 3, numAxialMTOP_KP(:), 0, '')
   call readData(readUnit, writeUnit, 3, kPointShiftMTOP(:), 0, '')
   call flush(writeUnit)

end subroutine readMTOPKPoints


! This converts kpoints from abc fractional (of the recip space cell) to
!   cartesian xyz (of the recip space cell).
subroutine convertKPointsToXYZ

   ! Include the modules we need
   use O_Kinds
   use O_Lattice, only: recipVectors

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the local variables used in this subroutine.
   integer :: i ! Loop index variables.
   real (kind=double) :: kPointX, kPointY, kPointZ

   ! Do the conversion.  We must use a temp variable since each kpoint
   !   component is used to calculate the other components of that kpoint.
   !   I.e., we are *overwriting* the abc fractional coordinates with
   !   xyz Cartesian. The units are inverse bohr radii.
   write (20,*) 'Kpoints in x,y,z cartesian form in recip space cell are:'
   do i = 1, numKPoints
      kPointX = dot_product(kPoints(:,i),recipVectors(1,:))
      kPointY = dot_product(kPoints(:,i),recipVectors(2,:))
      kPointZ = dot_product(kPoints(:,i),recipVectors(3,:))
      kPoints(1,i) = kPointX
      kPoints(2,i) = kPointY
      kPoints(3,i) = kPointZ

      write (20,100) i,kPoints(:,i)
   enddo

   100 format (i5,3f15.8)
end subroutine convertKPointsToXYZ


! This subroutine will compute the k dot r phase factors for each kpoint and
!   real space cell in the superlattice.
subroutine computePhaseFactors

   ! Include the modules we need
   use O_Kinds
   use O_Lattice, only: numCellsReal, cellDimsReal

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the local variables used in this subroutine.
   integer :: i,j
   real (kind=double) :: dotProduct

   ! First we allocate space to hold the matrix
   allocate (phaseFactor(numKPoints,numCellsReal))

   ! Assign the values as the dot product of kpoints and celldims.
!write(22,*) "Real Imaginary"
!write(24,*) "Real Imaginary"
   do i = 1, numCellsReal
      do j = 1, numKPoints
         dotProduct = dot_product(kPoints(:,j),cellDimsReal(:,i))
         phaseFactor(j,i) = cmplx(cos(dotProduct),sin(dotProduct),double)
!if ((j==2) .or. (j==26)) then
!write(20,*) "i,j",i,j,real(phaseFactor(j,i),double),aimag(phaseFactor(j,i))
!write(20,*) "cDR",cellDimsReal(:,i)
!endif
      enddo
   enddo

end subroutine computePhaseFactors


subroutine makePathKPoints

   ! Import necessary modules
   use O_Kinds
   use O_Constants, only: dim3
   use O_Lattice, only: recipVectors
   use O_TimeStamps

   ! Make sure that there are not accidental variable declarations.
   implicit none


   ! Define local parameters used in this subroutine
   integer :: i,j,k ! Loop index variables
   integer :: kPointCounter
   integer :: numSegmentKPoints
   real (kind=double) :: averageDelta
   real (kind=double), dimension (3) :: tempKPoint
   real (kind=double), dimension (3) :: segmentDelta
   real (kind=double), allocatable, dimension (:,:,:) :: symKPDistVect
   real (kind=double), allocatable, dimension (:,:)   :: symKPDistMag

   ! Begin the kpoint path creation.
   call timeStampStart(21)


   ! The main purpose of this subroutine is going to be to replace the kpoints
   !   specified in fort.15 with the kpoints defined by the path in the input
   !   file.  So, we begin by redefining the number of kpoints in the system.
   numKPoints = numPathKP

   ! Decallocate the kpoints and weights that may have been read in.
   if (allocated(kPointWeight)) then
      deallocate (kPointWeight)
      deallocate (kPoints)
   endif

   ! Reallocate the kpoint storage with the new parameters.  Note that the
   !   weights are going to store the magnitudes and not the usual kpoint
   !   weights.
   allocate (kPointWeight (numPathKP))
   allocate (kPoints      (dim3,numPathKP))

   ! Allocate space to hold the magnitudes along the path for each point.
   allocate (pathKPointMag (numPathKP))

   ! Allocate space for details of highly symmetric kpoints.
   allocate (symKPDistVect (3,maxval(numHighSymKP),numPaths))
   allocate (symKPDistMag  (maxval(numHighSymKP),numPaths))


   ! Convert the kpoints from fractional coordinates to cartesian if needed.
   !   Normally, the kpoints are given in fractional coordinates of the
   !   reciprocal space cell. This will convert the coordinates to Cartesian
   !   coordinates.
   if (isCartesian /= 1) then

      do i = 1, numPaths
         do j = 1, numHighSymKP(i)
            tempKPoint(:) = highSymKP(1,j,i) * recipVectors(:,1) + &
                          & highSymKP(2,j,i) * recipVectors(:,2) + &
                          & highSymKP(3,j,i) * recipVectors(:,3)
            highSymKP(:,j,i) = tempKPoint(:)
!write (20,*) "j,i = ",j,i
!write (20,fmt="(a,3f15.8)") "highSymKP(:,j,i)",highSymKP(:,j,i)
         enddo
      enddo
   endif

   ! Identify distances between each consecutive symmetric kpoint on the path.
   !   Note that when transitioning from one path to the next, the distance
   !   between those two kpoints should be zero.
   do i = 1, numPaths
      do j = 1, numHighSymKP(i)
         if ((i == 1) .and. (j == 1)) then ! First point on first path
            symKPDistVect(:,1,i) = 0.0_double
            symKPDistMag(1,i)    = 0.0_double
         elseif (j == 1) then ! First point on any path after the first path
               ! needs to be equal in position and magnitude to the last point
               ! on the previous path. This will signify to later programs and
               ! subroutines that we have started a new path.
            symKPDistVect(:,1,i) = symKPDistVect(:,numHighSymKP(i-1),i-1)
            symKPDistMag(1,i)    = symKPDistMag(numHighSymKP(i-1),i-1)
         else
            symKPDistVect(:,j,i) = highSymKP(:,j,i) - highSymKP(:,j-1,i)
            symKPDistMag(j,i) = symKPDistMag(j-1,i) &
                  + sqrt(sum(symKPDistVect(:,j,i)**2))
         endif
!write (20,*) "j,i = ",j,i
!write (20,fmt='(a,3f15.8)') "symKPDistVect(:,j,i)",symKPDistVect(:,j,i)
!write (20,*) "symKPDistMag(j,i)",symKPDistMag(j,i)
      enddo
   enddo

   ! Consider the special case where the number of requested path kpoints is
   !   equal to the number of highly symmetric kpoints provided.
   if (numPathKP == sum(numHighSymKP(:))) then

      ! Copy the highly symmetric kpoint coordinates to the path kpoints.
      kPointCounter = 0
      do i = 1, numPaths
         do j = 1, numHighSymKP(i)
            kPointCounter = kPointCounter + 1
            kPoints(:,kPointCounter)     = highSymKP(:,j,i)
            pathKPointMag(kPointCounter) = symKPDistMag(j,i)

            ! Record the index number of the other high symmetry k points.
            write (20,*) 'HIGH SYMMETRY K POINT INDEX NUMBER :',i
         enddo
      enddo

   else

      ! Determine the scaling factor for the distance between kpoints in a
      !   given segment between highly symmetric kpoints. This is the most
      !   most distant kpoint magnitude divided by the number of kpoints.
      averageDelta = symKPDistMag(numHighSymKP(numPaths),numPaths) / &
            & (numPathKP - 1)
!write (20,*) "numPaths = ",numPaths
!write (20,*) "numHighSymKP = ",numHighSymKP(:)
!write (20,*) "averageDelta = ",averageDelta

      ! Initialize a kpoint counter
      kPointCounter = 0

      do i = 1, numPaths
         do j = 1, numHighSymKP(i)
            if (j /= 1) then

               ! Determine the number of k points to use for the current
               !   segment between highly symmetric k points (j) and (j-1) of
               !   this path.
               numSegmentKPoints = int((symKPDistMag(j,i) - &
                     & symKPDistMag(j-1,i)) / averageDelta)

               ! Check that this number of segment kpoints will not take us
               !   over the limit. If it does, then reduce the number of
               !   points in this segment.
               if (kPointCounter + numSegmentKPoints > numPathKP) then
                  numSegmentKPoints = numPathKP - kPointCounter
               endif
!write (20,*) "j,i = ",j,i
!write (20,*) "numSegmentKPoints = ",numSegmentKPoints
!write (20,*) "symKPDistMag(j,i) = ",symKPDistMag(j,i)
!write (20,*) "symKPDistMag(j-1,i) = ",symKPDistMag(j-1,i)
!call flush(20)
   
               ! If there are not going to be any segment kpoints for this
               !   segment then cycle to the next symmetric kpoint. This is
               !   most commonly the case when we start a new path.
               if (numSegmentKPoints >= 1) then
   
                  ! Get the size of the x,y,z deltas for this segment
                  segmentDelta(:) = symKPDistVect(:,j,i) / numSegmentKPoints
      
                  ! Loop to assign the k point positions for this segment
                  do k = 1, numSegmentKPoints - 1
      
                     ! Increment the K point counter
                     kPointCounter = kPointCounter + 1
      
                     ! Store the position of the next kpoint on this segment by
                     !   adding the above determined delta to the last known
                     !   path kpoint.
                     kPoints(:,kPointCounter) = &
                           & kPoints(:,kPointCounter - 1) + segmentDelta(:)
      
                     ! Calculate the magnitude of the vector for the above path
                     !   KPoint.
                     pathKPointMag(kPointCounter) = &
                           & pathKPointMag(kPointCounter - 1) + &
                           & sqrt(sum(segmentDelta(:)**2))
                  enddo
               endif
            endif
   
            ! Now, record the j-th high symmetry k point. For the first kpoint
            !   (where i==1 and j==1) this will happen first, before the code
            !   segment above.
   
            ! Increment the K point counter.
            kPointCounter = kPointCounter + 1
   
            ! Record the index number of the high symmetry k point.
            write (20,*) 'HIGH SYMMETRY K POINT INDEX NUMBER :',kPointCounter
   
            ! Store the position of the high symmetry k point in the path array.
            kPoints(:,kPointCounter) = highSymKP(:,j,i)
   
            ! Calculate the magnitude of the vector for this high symm KPoint.
            if ((i==1) .and. (j==1)) then
               pathKPointMag(kPointCounter) = 0.0_double
            elseif (j==1) then
               pathKPointMag(kPointCounter) = pathKPointMag(kPointCounter - 1)
            else
               pathKPointMag(kPointCounter) = &
                     & pathKPointMag(kPointCounter - 1) + &
                     & sqrt(sum(segmentDelta(:)**2))
            endif
         enddo
      enddo
   endif

   ! Assign weighting factors to each kpoint on the path so that they are all
   !   equal.
   kPointWeight(:) = 2.0_double / kPointCounter

   ! Use the newly calculated number of kpoints for all future program parts.
   numPathKP   = kPointCounter
   numKPoints  = kPointCounter

   ! Record the index, magnitude and position of each kpoint on the path.
   do i = 1, kPointCounter
      write (20,fmt="(i5,4f8.3,f18.12)") i,pathKPointMag(i),kPoints(:,i),&
            & kPointWeight(i)
   enddo

   deallocate (symKPDistVect)
   deallocate (symKPDistMag)

   ! Finish the kpoint path creation.
   call timeStampEnd(21)

end subroutine makePathKPoints


subroutine initializeKPointMesh(applySymmetry)

   ! Define used modules.
   use O_Kinds
   use O_Constants

   ! Make sure no funny variables are created.
   implicit none

   ! Define passed parameters.
   integer, intent (in) :: applySymmetry

   ! Define local variables.
   integer :: i, j, k, m, n
   integer :: numMeshKPoints
   integer :: numFoldedKPoints
   integer :: isMatch
   integer, dimension (3) :: loopIndex
   integer, allocatable, dimension(:) :: kPointTracker
   real (kind=double) :: kpThresh
   real (kind=double) :: weightSum
   real (kind=double) :: initWeight
   real (kind=double), dimension(3) :: abcDelta
   real (kind=double), dimension(3) :: foldedKPoint
   real (kind=double), allocatable, dimension(:,:) :: &
         & abcMeshKPoints
   real (kind=double), allocatable, dimension(:) :: &
         & foldedWeight
   real (kind=double), allocatable, dimension(:,:) :: &
         & abcFoldedKPoints

   ! Define the threshold for when two kpoints are considered symmetrically
   !   equivalent under point group operations.
   kpThresh = 0.00001_double

   ! Compute the number of uniform mesh points.
   numMeshKPoints = product(numAxialKPoints(:))

   ! Initialize the total weight to 2. (Two electrons per state by default.
   !   Spin-polarized calculations account for the single electron per state
   !   separately via the 'spin' variable.)
   weightSum = 2.0_double

   ! Compute the mesh step size in fractional units.
   abcDelta(:) = 1.0_double / numAxialKPoints(:)

   ! Allocate space for the uniform mesh.
   allocate (abcMeshKPoints(3,numMeshKPoints))

   ! Create the initial uniform mesh in fractional units
   !   of the abc reciprocal space cell.
   numMeshKPoints = 0
   do i = 1, numAxialKPoints(1)
      loopIndex(1) = i
      do j = 1, numAxialKPoints(2)
         loopIndex(2) = j
         do k = 1, numAxialKPoints(3)
            numMeshKPoints = numMeshKPoints + 1
            loopIndex(3) = k
            abcMeshKPoints(:,numMeshKPoints) = &
                  & -0.5_double + (loopIndex(:) &
                  & - 1.0_double + kPointShift(:)) &
                  & * abcDelta(:)
         enddo
      enddo
   enddo

   ! Deallocate any previously allocated kpoint arrays.
   if (allocated(kPoints)) then
      deallocate(kPoints)
      deallocate(kPointWeight)
   endif
   if (allocated(fullToIBZMap)) then
      deallocate(fullToIBZMap)
   endif

   ! Record the full mesh size for tetrahedra generation
   !   and eigenvalue unfolding.
   numFullMeshKP = numMeshKPoints

   ! If symmetry reduction is not requested, copy the full mesh directly.
   !   Build an identity mapping (every point is its own IBZ representative).
   if (applySymmetry == 0) then
      allocate(kPoints(3,numMeshKPoints))
      allocate(kPointWeight(numMeshKPoints))
      kPoints(:,:) = abcMeshKPoints(:,:)
      kPointWeight(:) = weightSum &
            & / real(numMeshKPoints,double)
      numKPoints = numMeshKPoints
      allocate(fullToIBZMap(numMeshKPoints))
      do i = 1, numMeshKPoints
         fullToIBZMap(i) = i
      enddo
      deallocate(abcMeshKPoints)
      return
   endif

   ! -----------------------------------------------
   ! IBZ symmetry reduction (fold the mesh).
   ! -----------------------------------------------
   ! The algorithm applies each point group operation to each mesh kpoint.
   !   If the rotated kpoint matches another mesh kpoint (within kpThresh),
   !   those two points are symmetry-equivalent. The matching point is
   !   removed from further consideration and its weight is added to the
   !   irreducible representative.

   ! Allocate working arrays for the folding process. We allocate the full
   !   mesh size because we do not yet know how many irreducible points
   !   there will be.
   allocate (kPointTracker(numMeshKPoints))
   allocate (abcFoldedKPoints(3,numMeshKPoints))
   allocate (foldedWeight(numMeshKPoints))

   ! Compute the initial weight of each uniform mesh
   !   kpoint.
   initWeight = weightSum / real(numMeshKPoints,double)

   ! Initialize the kpoint tracker. Each mesh point starts as its own
   !   representative (tracker(i) == i). When a point is found to be equivalent
   !   to an earlier one, its tracker is set to the negative of the IBZ
   !   index of that representative.
   do i = 1, numMeshKPoints
      kPointTracker(i) = i
   enddo

   ! Initialize the count of irreducible kpoints.
   numFoldedKPoints = 0

   ! Consider each mesh kpoint in turn.
   do i = 1, numMeshKPoints

      ! Skip points already matched to an earlier IBZ
      !   representative.
      if (kPointTracker(i) /= i) cycle

      ! This is a new irreducible kpoint.
      numFoldedKPoints = numFoldedKPoints + 1
      kPointTracker(i) = -numFoldedKPoints
      foldedWeight(numFoldedKPoints) = initWeight
      abcFoldedKPoints(:,numFoldedKPoints) = &
            & abcMeshKPoints(:,i)

      ! If this is the last mesh point, no further
      !   comparisons are needed.
      if (i == numMeshKPoints) cycle

      ! Apply each point group operation to this kpoint and look for matches
      !   among the remaining unmatched mesh points.
      do m = 1, numPointOps

         ! Compute the folded (rotated) kpoint by
         !   applying point group operation m.
         do n = 1, 3
            foldedKPoint(n) = &
                  & sum(abcRecipPointOps(n,:,m) &
                  & * abcMeshKPoints(:,i))
         enddo

         ! Compare against remaining unmatched points.
         do j = i + 1, numMeshKPoints
            if (kPointTracker(j) /= j) cycle

            ! Check if the difference is negligible
            !   along all three axes.
            isMatch = 1
            do k = 1, 3
               if (abs(foldedKPoint(k) &
                     & - abcMeshKPoints(k,j)) &
                     & > kpThresh) then
                  isMatch = 0
                  exit
               endif
            enddo

            ! If matched, absorb this point's weight
            !   and mark it as folded.
            if (isMatch == 1) then
               foldedWeight(numFoldedKPoints) = &
                     & foldedWeight(numFoldedKPoints) &
                     & + initWeight
               kPointTracker(j) = -numFoldedKPoints
            endif

         enddo ! j (remaining mesh points)
      enddo ! m (point group operations)
   enddo ! i (mesh points)

   ! Copy the irreducible kpoints into the module-level
   !   kPoints and kPointWeight arrays.
   allocate(kPoints(3,numFoldedKPoints))
   allocate(kPointWeight(numFoldedKPoints))
   kPoints(:,1:numFoldedKPoints) = &
         & abcFoldedKPoints(:,1:numFoldedKPoints)
   kPointWeight(1:numFoldedKPoints) = &
         & foldedWeight(1:numFoldedKPoints)
   numKPoints = numFoldedKPoints

   ! Build the full-to-IBZ mapping from kPointTracker. kPointTracker(i) ==
   !   -ibzIndex for every mesh point. Convert to positive IBZ indices.
   allocate(fullToIBZMap(numFullMeshKP))
   do i = 1, numFullMeshKP
      fullToIBZMap(i) = -kPointTracker(i)
   enddo

   ! Clean up working arrays.
   deallocate(abcMeshKPoints)
   deallocate(kPointTracker)
   deallocate(abcFoldedKPoints)
   deallocate(foldedWeight)

end subroutine initializeKPointMesh


subroutine initializeKPoints (inSCF)

   ! Define used modules.
   use O_Kinds
   use O_Constants, only: dim3
   use O_CommandLine, only: doSYBD_SCF, doSYBD_PSCF, doMTOP_SCF, doMTOP_PSCF

   ! Prevent accidental variable definition.
   implicit none

   ! Declare passed parameters.
   integer, intent(in) :: inSCF

   ! Define local variables.

   ! The considerations for this operation are:

   ! At this point an SCF or PSCF kpoint input file has been read in. The
   !   olcao.dat file has been read in. Implicit information has been
   !   computed. The lattices have been initialized.

   ! If either doSYBD_SCF or doSYBD_PSCF variable was set, then we need to
   !   shift to use the defined SYBD kpoints regardless of any other setting.

   ! If either doMTOP_SCF or doMTOP_PSCF variable was set, then we need to
   !   shift to use the defined MTOP kpoints regardless of any other setting.
   !   The MTOP kpoints are a uniform mesh given as a number of kpoints in
   !   each of the a,b,c directions along with the amount of shift from the
   !   origin.

   ! If none of doSYBD_SCF, doSYBD_PSCF, doMTOP_SCF, and doMTOP_PSCF are set,
   !   then we turn to the kpoints as defined in the kpoint input file.
   !   Otherwise, if one of the above situations was true, then there is
   !   nothing left to do and we can leave this subroutine.

   ! If the kpoint style was set to zero, then the kpoints were explicitly
   !   provided in the kpoint input file, the number of kpoints was set,
   !   the kpoint weights were set, and the kpoint values themselves were
   !   set. Therefore, there is nothing left to do and we can leave this
   !   subroutine.

   ! If the kpoint style was set to one, then the number of kpoints in each
   !   a,b,c direction were provided along with the amount of shift that
   !   should be applied to the mesh. Therefore, it is necessary to create
   !   a mesh from this information and reduce it by symmetry.

   ! If the kpoint style was set to two, then the density of kpoints was
   !   provided along with the amount of shift that should be applied to the
   !   mesh. Therefore, it is necessary to create a mesh from this information
   !   and reduce it by symmetry.

   ! Proceed:

   ! Check for the above sequence of situations.
   if ((doSYBD_SCF == 1) .and. (inSCF == 1)) then
      call makePathKPoints
   elseif ((doSYBD_PSCF == 1) .and. (inSCF == 0)) then
      call makePathKPoints
   elseif ((doMTOP_SCF == 1) .and. (inSCF == 1)) then
      numAxialKPoints(:) = numAxialMTOP_KP(:)
      kPointShift(:) = kPointShiftMTOP(:)
      call initializeKPointMesh(0) ! Do not apply symmetry.
      call convertKPointsToXYZ
      call makeMTOPIndexMap
   elseif ((doMTOP_PSCF == 1) .and. (inSCF == 0)) then
      numAxialKPoints(:) = numAxialMTOP_KP(:)
      kPointShift(:) = kPointShiftMTOP(:)
      call initializeKPointMesh(0) ! Do not apply symmetry.
      call convertKPointsToXYZ
      call makeMTOPIndexMap
   elseif (kPointStyleCode == 0) then
      call convertKPointsToXYZ
   elseif (kPointStyleCode == 1) then
      call initializeKPointMesh(1) ! Apply symmetry to reduce KP.
      call convertKPointsToXYZ
   elseif (kPointStyleCode == 2) then
      call computeAxialKPoints
      call computeRecipPointOps
      call initializeKPointMesh(1) ! Apply symmetry to reduce KP.
      call convertKPointsToXYZ
   endif

   ! If the LAT integration method was requested, generate the tetrahedra
   !   and compute the BZ integration weight per tetrahedron. This must
   !   happen after the mesh is built (so that numAxialKPoints is set).
   !   Tetrahedra reference the full uniform mesh, not the IBZ-reduced
   !   kpoints.
   if (kPointIntgCode == 1) then
      call generateTetrahedra
      call computeTetraVol
   endif

   ! Compute the phase factors.
   call computePhaseFactors

end subroutine initializeKPoints


subroutine computeAxialKPoints

   ! Define used modules.
   use O_Kinds
   use O_Constants, only:dim3
   use O_Lattice, only: recipCellVolume, recipMag

   ! Define local variables.
   integer :: i
   integer :: nearIndex
   real (kind=double) :: stepSize
   real (kind=double) :: numKPointsEstimate
   real (kind=double), dimension (dim3) :: fractionalPoints

   ! The user specifies a volume density (kpoints per unit reciprocal-space
   !   volume in Bohr^-3). Multiply by the reciprocal cell volume to get the
   !   target total number of kpoints in the full mesh.
   numKPointsEstimate = minKPointDensity * recipCellVolume

   ! Distribute the total count across the three axes so that the spacing is
   !   as uniform as possible (equal step size along each axis).

   ! To force the step sizes to be equal for each axis, the ratio of the
   !   axial magnitude and the number of kpoints along that axis should be
   !   equal for all axes. I.e., recipMag(1) / numAxialKPoints(1) = stepSize,
   !   and so forth for 2, and 3. With a fixed stepSize this would yield
   !   real numbers for the number of kpoints. We will assume this for now
   !   and then make integer later.

   ! Now, establish a point density constraint:
   !   recipCellVolume / product(numAxialKPoints(:)) = 1 / minKPointsDensity

   ! Therefore, the step size is:
   stepSize = (product(recipMag(:)) / &
         & (recipCellVolume * minKPointDensity))**(1.0_double/3.0_double)

   ! Now, solve for the number of axial points:
   numAxialKPoints(:) = int(recipMag(:) / stepSize)

   ! Compute the residual fractional part.
   fractionalPoints(:) = (recipMag(:) / stepSize) &
         & - real(numAxialKPoints(:),double)

   ! Finally, adjust for the shift to integer values.

   ! First, increment the number of kpoints along axes that are closest to
   !   the next integer number of points until the density exceeds the target.
   do while (minKPointDensity * recipCellVolume > product(numAxialKPoints(:)))
      nearIndex = maxloc(fractionalPoints(:),1)
      numAxialKPoints(nearIndex) = numAxialKPoints(nearIndex) + 1
      fractionalPoints(nearIndex) = 0.0_double
   enddo

   ! Then, decrement the number of kpoints along the axes that are closest to
   !   the previous integer (floor) until the density is less than the target.
   !   First, fix any fractionalPoints that are zero to be one so that they
   !   are never reduced.
   do i = 1, 3
      if (fractionalPoints(i) == 0.0_double) then
         fractionalPoints(i) = 1.0_double
      endif
   enddo
   do while (minKPointDensity * recipCellVolume < product(numAxialKPoints(:)))
      nearIndex = minloc(fractionalPoints(:),1)
      numAxialKPoints(nearIndex) = numAxialKPoints(nearIndex) - 1
      if (numAxialKPoints(nearIndex) == 0) then
         numAxialKPoints(nearIndex) = 1
      endif
      fractionalPoints(nearIndex) = 1.0_double
   enddo

   ! At this point, the numAxialKPoints should be optimized to be as close as
   !   possible to the target density while keeping inter-point spacing as
   !   constant as possible (comparing axes).

end subroutine computeAxialKPoints


subroutine makeMTOPIndexMap

   ! Use necessary modules

   ! Make sure no variables are accidentally declared.
   implicit none

   ! Declare local variables.
   integer :: kPointCount
   integer :: i,j,k

   ! Allocate space to hold the MTOP index map.
   allocate (mtopKPMap(3,numKPoints))

   ! Build a map between the linear list of kpoints as they were created in
   !   the initializeKPointMesh subroutine and the order in which the c-axis
   !   strings of kpoints are needed for MTOP calculations.
   kPointCount = 0
   do i = 1, numAxialKPoints(1)
      do j = 1, numAxialKPoints(2)
         do k = 1, numAxialKPoints(3)
            kPointCount = kPointCount + 1
            mtopKPMap(3,kPointCount) = kPointCount
         enddo
      enddo
   enddo

   ! Build a map between the linear list of kpoints as they were created in
   !   the initializeKPointMesh subroutine and the order in which the b-axis
   !   strings of kpoints are needed for MTOP calculations.
   kPointCount = 0
   do i = 1, numAxialKPoints(1)
      do j = 1, numAxialKPoints(3)
         do k = 1, numAxialKPoints(2)
            kPointCount = kPointCount + 1
            mtopKPMap(2,kPointCount) = getIndexFromIndices(i,k,j)
         enddo
      enddo
   enddo

   ! Build a map between the linear list of kpoints as they were created in
   !   the initializeKPointMesh subroutine and the order in which the a-axis
   !   strings of kpoints are needed for MTOP calculations.
   kPointCount = 0
   do i = 1, numAxialKPoints(2)
      do j = 1, numAxialKPoints(3)
         do k = 1, numAxialKPoints(1)
            kPointCount = kPointCount + 1
            mtopKPMap(1,kPointCount) = getIndexFromIndices(k,i,j)
         enddo
      enddo
   enddo

!write(20,*) "AXIS1"
!do i = 1,kPointCount
!   write(20,*) "i map kp",i,mtopKPMap(1,i), kPoints(:,i)
!enddo
!write(20,*) "AXIS2"
!do i = 1,kPointCount
!   write(20,*) "i map kp",i,mtopKPMap(2,i), kPoints(:,i)
!enddo
!write(20,*) "AXIS3"
!do i = 1,kPointCount
!   write(20,*) "i map kp",i,mtopKPMap(3,i), kPoints(:,i)
!enddo

end subroutine makeMTOPIndexMap


! Given the a,b,c indices of a kpoint in the full kpoint mesh, this function
!   will return then index of the same kpoint in the linear list of kpoints.
function getIndexFromIndices(a,b,c)

   implicit none

   ! Passed parameters.
   integer, intent(in) :: a, b, c
   integer :: getIndexFromIndices

   ! The linear list of kpoints is created by nested loops over kpoints in
   !   the a, b, c axes in that order. Hence, the c-axis kpoints iterate
   !   fastest through the linear list, followed by the b-axis kpoints. So
   !   the linear list index is the current c index + nAKP(2) times the
   !   number of full lines on the c-axis we passed to get to the current
   !   b-axis point, followed by nAKP(2)*nAKP(3) times the number of full
   !   c,b planes we passed to get to the current a-axis point.
   getIndexFromIndices = c + (b-1)*numAxialKPoints(3) + &
         & (a-1)*numAxialKPoints(2)*numAxialKPoints(3)

end function getIndexFromIndices


! Generate the tetrahedra that tile the reciprocal cell.
!   The uniform Monkhorst-Pack mesh defines a grid of nA x nB x nC
!   parallelepipeds. Each parallelepiped has 8 corners and is decomposed
!   into exactly 6 tetrahedra that share the main diagonal M1-M8 (Bloechl
!   1994). The mesh is periodic, so indices wrap with modular arithmetic.
!
!   Parallelepiped corners at grid position (a, b, c):
!     M1 = (a,   b,   c  )  M5 = (a+1, b+1, c  )
!     M2 = (a+1, b,   c  )  M6 = (a+1, b,   c+1)
!     M3 = (a,   b+1, c  )  M7 = (a,   b+1, c+1)
!     M4 = (a,   b,   c+1)  M8 = (a+1, b+1, c+1)
!
!   Six tetrahedra sharing diagonal M1-M8:
!     T1: M1, M2, M5, M8   T4: M1, M4, M7, M8
!     T2: M1, M3, M5, M8   T5: M1, M4, M6, M8
!     T3: M1, M3, M7, M8   T6: M1, M2, M6, M8
!
!   The tetrahedra indices reference the full uniform
!   mesh (1..numFullMeshKP), not the IBZ-reduced list.
subroutine generateTetrahedra

   implicit none

   ! Define local variables.
   integer :: a, b, c, t
   integer :: nA, nB, nC
   integer :: M1, M2, M3, M4, M5, M6, M7, M8

   ! Local shorthand for the axial kpoint counts.
   nA = numAxialKPoints(1)
   nB = numAxialKPoints(2)
   nC = numAxialKPoints(3)

   ! Compute the total number of tetrahedra: 6 per
   !   parallelepiped, one parallelepiped per grid cell.
   numTetrahedra = 6 * nA * nB * nC

   ! Allocate the tetrahedra array.
   if (allocated(tetrahedra)) deallocate(tetrahedra)
   allocate (tetrahedra(4, numTetrahedra))

   ! Build the tetrahedra by iterating over every grid
   !   cell and decomposing it into 6 tetrahedra.
   t = 0
   do a = 1, nA
      do b = 1, nB
         do c = 1, nC

            ! Compute the 8 corner indices of this
            !   parallelepiped with periodic wrapping.
            M1 = getIndexFromIndices( &
                  & a, b, c)
            M2 = getIndexFromIndices( &
                  & mod(a, nA) + 1, b, c)
            M3 = getIndexFromIndices( &
                  & a, mod(b, nB) + 1, c)
            M4 = getIndexFromIndices( &
                  & a, b, mod(c, nC) + 1)
            M5 = getIndexFromIndices( &
                  & mod(a, nA) + 1, &
                  & mod(b, nB) + 1, c)
            M6 = getIndexFromIndices( &
                  & mod(a, nA) + 1, &
                  & b, mod(c, nC) + 1)
            M7 = getIndexFromIndices( &
                  & a, mod(b, nB) + 1, &
                  & mod(c, nC) + 1)
            M8 = getIndexFromIndices( &
                  & mod(a, nA) + 1, &
                  & mod(b, nB) + 1, &
                  & mod(c, nC) + 1)

            ! Six tetrahedra sharing diagonal M1-M8.
            tetrahedra(:, t+1) = (/M1, M2, M5, M8/)
            tetrahedra(:, t+2) = (/M1, M3, M5, M8/)
            tetrahedra(:, t+3) = (/M1, M3, M7, M8/)
            tetrahedra(:, t+4) = (/M1, M4, M7, M8/)
            tetrahedra(:, t+5) = (/M1, M4, M6, M8/)
            tetrahedra(:, t+6) = (/M1, M2, M6, M8/)
            t = t + 6

         enddo ! c
      enddo ! b
   enddo ! a

end subroutine generateTetrahedra


! Compute the BZ integration weight for each tetrahedron. All tetrahedra tile
!   the BZ uniformly, so each one represents an equal fraction 1/numTetrahedra
!   of the full zone.
subroutine computeTetraVol

   use O_Kinds

   implicit none

   tetraVol = 1.0_double &
         & / real(numTetrahedra, double)

end subroutine computeTetraVol


! Convert the abc-coordinate point group operations into reciprocal-space
!   abc-coordinate operations. This is needed for IBZ symmetry reduction
!   of the kpoint mesh. The transformation uses the real and reciprocal
!   lattice vectors: R_recip = R^T * R_real * R_recip_lattice, where the
!   dot product between real and reciprocal lattice vectors provides the
!   metric. The resulting abcRecipPointOps can be applied directly to kpoints
!   expressed in fractional reciprocal coordinates.
subroutine computeRecipPointOps

   ! Import necessary modules.
   use O_Kinds
   use O_Constants, only: pi
   use O_Lattice, only: realVectors, recipVectors

   ! Make sure no funny variables are defined.
   implicit none

   ! Define local variables.
   integer :: i, j, k, l
   real (kind=double), dimension (3,3) :: tempPointOp

   ! Allocate storage for the reciprocal-space operations.
   allocate (abcRecipPointOps(3,3,numPointOps))

   do i = 1, numPointOps

      ! Initialize the accumulator for this operation.
      tempPointOp(:,:) = 0.0_double

      ! Transform the abc real-space point group operation
      !   into an abc reciprocal-space operation using the
      !   real and reciprocal lattice vectors.
      do j = 1, 3    ! x,y,z of the real lattice
         do k = 1, 3 ! x,y,z of the reciprocal lattice
            do l = 1, 3 ! a,b,c of the reciprocal lattice

               ! Array-loop over real a,b,c.
               tempPointOp(:,l) = tempPointOp(:,l) &
                     & + realVectors(j,:) &
                     & * abcPointOps(j,k,i) &
                     & * recipVectors(k,l) &
                     & / (2.0_double * pi)

            enddo ! l (recip a,b,c)
         enddo ! k (recip x,y,z)
      enddo ! j (real x,y,z)

      ! Save the transformed operation.
      abcRecipPointOps(:,:,i) = tempPointOp(:,:)

   enddo ! i (operations)

end subroutine computeRecipPointOps


subroutine cleanUpKPoints

   implicit none

   deallocate (kPointWeight)
   deallocate (kPoints)
   deallocate (numHighSymKP)
   deallocate (highSymKP)
   deallocate (highSymKPChar)

   if (allocated(pathKPointMag)) then
      deallocate(pathKPointMag)
   endif

   ! Only allocated in setup.
   if (allocated(phaseFactor)) then
      deallocate (phaseFactor)
   endif

   ! Only allocated for MTOP.
   if (allocated(mtopKPMap)) then
      deallocate (mtopKPMap)
   endif

   ! Only allocated for density-based kpoint input
   !   (kPointStyleCode == 2).
   if (allocated(abcPointOps)) then
      deallocate (abcPointOps)
   endif
   if (allocated(abcRecipPointOps)) then
      deallocate (abcRecipPointOps)
   endif

   ! Only allocated when the mesh is built internally
   !   (style codes 1, 2) with or without IBZ reduction.
   if (allocated(fullToIBZMap)) then
      deallocate (fullToIBZMap)
   endif

   ! Only allocated when kPointIntgCode == 1 (LAT).
   if (allocated(tetrahedra)) then
      deallocate (tetrahedra)
   endif

end subroutine cleanUpKPoints


end module O_KPoints
