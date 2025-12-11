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
   real (kind=double) :: minKPointDensity ! The number of kpoints along each
         !   a,b,c axis divided by the magnitude of that axis must be greater
         !   than or equal to the value of minKPointDensity. This is only used
         !   if the kPointStyleCode equals 2.
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
      ! Read fractional minimum kpoint separation and a shift.

      ! Read the minimum line density required of kpoints along all axes.
      call readData(readUnit,writeUnit,minKPointDensity,&
            & len('MIN_KP_LINE_DENSITY'),'MIN_KP_LINE_DENSITY')

      ! Read the fractional distance that the uniform mesh should be shifted
      !   away from the origin along each a,b,c axis.
      call readData(readUnit,writeUnit,3,kPointShift,len('KP_SHIFT_A_B_C'),&
            & 'KP_SHIFT_A_B_C')

      ! The kpoint mesh will be constructed later, once implicit information
      !   such as the reciprocal lattice has been obtained.
   endif


   ! Set up additional information depending on type of integration method.
   if (KPointIntgCode == 0) then
      return
   elseif (KPointIntgCode == 1) then
      call generateTetrahedra ! Linear Analytic Tetrahedral (LAT) method.
   endif

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
   !   xyz Cartesian.
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
   do i = 1, numCellsReal
      do j = 1, numKPoints
         dotProduct = dot_product(kPoints(:,j),cellDimsReal(:,i))
         phaseFactor(j,i) = cmplx(cos(dotProduct),sin(dotProduct),double)
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

   ! Define local variables
   integer :: i, j, k
   integer :: numMeshKPoints
   integer, dimension (3) :: loopIndex
   !integer, allocatable, dimension(:) :: kPointTracker
   real (kind=double) :: kpThresh ! Threshhold for which two kpoints are
         ! considered to be symmetricly the same via point group operations.
   real (kind=double) :: weightSum ! Sum of initial weights of all kpoints.
         ! weightSum = 1 for relativistic; 2 for non-relativistic
   real (kind=double), dimension(3) :: abcDelta
   real (kind=double), allocatable, dimension(:,:) :: abcMeshKPoints
   !real (kind=double), allocatable, dimension(:,:) :: abcFoldedKPoints

   ! Define the threshhold for when two kpoints are symmetry reduced via
   !   point group operations.  (i.e. how close do two kpoints have to be
   !   when one is folded near the other for them to be considered the same?)
   kpThresh = 0.00001_double  ! A bit arbitrary!

   ! Compute the number of uniform mesh points in order to allocate space to
   !   hold the kpoints.
   numMeshKPoints = product(numAxialKPoints(:))

   ! Initialize the weight to 2 always.  (There are 2 electrons per state
   !   by default.  If a spin-polarized calculation is done, then the 1 e-
   !   per state is accounted for in the olcao fortran program itself.)
   weightSum = 2.0_double

   ! Compute the mesh points separation step size in fractional units.
   abcDelta(:) = 1.0_double / numAxialKPoints(:)

   ! Allocate necessary space.
   !allocate (kPointTracker(numMeshKPoints))
   allocate (abcMeshKPoints(3,numMeshKPoints))

   ! Loop over abc axial kpoints numbers to create the initial uniform mesh in
   !   fractional units of the abc reciprocal space cell.
   numMeshKPoints = 0
   do i = 1, numAxialKPoints(1)

      ! Store the i loop index.
      loopIndex(1) = i

      do j = 1, numAxialKPoints(2)

         ! Store the j loop index.
         loopIndex(2) = j

         do k = 1, numAxialKPoints(3)

            ! Increment the number of uniform mesh kpoints.
            numMeshKPoints = numMeshKPoints + 1

            ! Store the k loop index.
            loopIndex(3) = k

            ! Compute the current mesh kpoint location.
            abcMeshKPoints(:,numMeshKPoints) = &
                  & (loopIndex(:)-1.0_double+kPointShift(:)) * abcDelta(:)

         enddo  ! k=1,numAxialKPoints(3)
      enddo  ! j=1,numAxialKPoints(2)
   enddo  ! i=1,numAxialKPoints(1)

   ! There are a variety of ways to traverse the linear list of kPoints:
   !   Shorthand:
   !      nAKP(q) is numAxialKPoints(q) with q = a,b,c.
   !      nAKP(q,s) is numAxialKPoints(q,s) with q=a,b,c and s being
   !         the step size.
   !   Given a loop structure over the a,b,c nAKP, the specific kpoint
   !       in the linear array can be found as loop_c*1 + loop_b*len_c
   !       + loop_a*len_b*len_c.
   !
   !   (1) Iterating a counter using a loop structure along the lines of
   !       nAKP(1,1)->nAKP(2,1)->nAKP(3,1) will visit all mesh points
   !       in a series of strings of length numAxialKPoints(3) along the c
   !       axis. Because all incremental step sizes are one, this sequence
   !       will exactly and sequentially visit all elements in the kpoints
   !       array.
   !   (2) Iterating a counter using a loop structure along the lines of
   !       nAKP(1,1)->nAKP(3,1)->nAKP(2,nAKP(3)) will visit all mesh
   !       points in a series of strings of length numAxialKPoints(2) along
   !       the b axis.
   !   (3) Iterating a counter using a loop structure along the lines of
   !       nAKP(2,1)->nAKP(3,1)->nAKP(1,nAKP(2)*nAKP(3)) will visit all
   !       mesh points in a series of strings of length numAxialKPoins(1)
   !       along the a axis.

   ! Presently, we just make the full mesh here and do not reduce by symmetry.
   !   In the future, we will do symmetry reductions here. Until then, we just
   !   copy the abcMeshKPoints into the list of kPoints.

   ! At this point we have the mesh in abcMeshKPoints and the total number of
   !   kPoints in numMeshKPoints. If we are here it is because we are not using
   !   the kPoints given explicitly in the kPoints input file. So, we need to
   !   manage the allocation of space to hold the weights and kpoints.
   if (allocated(kPoints)) then
      deallocate(kPoints)
      deallocate(kPointWeight)
   endif

   ! Once deallocated (or if they were never allocated to begin with), we can
   !   allocate space to hold the kpoints and their weights.
   allocate(kPoints(3,numMeshKPoints))
   allocate(kPointWeight(numMeshKPoints))

   ! Now, finally, we can copy the kPoints until we add the symmetry part
   !   later. Also, we can initialize the weights.
   kPoints(:,:) = abcMeshKPoints(:,:)
   kPointWeight(:) = weightSum/real(numMeshKPoints,double)
   numKPoints = numMeshKPoints

   deallocate(abcMeshKPoints)

   ! If we don't need symmetry, then return.
   if (applySymmetry == 0) then
      return
   endif

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
      call initializeKPointMesh(1) ! Apply symmetry to reduce KP.
      call convertKPointsToXYZ
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

   ! To determine the number of axial KPoints needed to achieve the desired
   !   kpoint density, we first estimate the total number of kpoints needed.
   numKPointsEstimate = minKPointDensity * recipCellVolume

   ! The number of kpoints along each axis should be such that the spacing
   !   is nearly equal (one axis vs. another) and the product of the number
   !   of kpoints along each axis is as close to the estimate as possible.

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

write(20,*) "AXIS1"
do i = 1,kPointCount
   write(20,*) "i map kp",i,mtopKPMap(1,i), kPoints(:,i)
enddo
write(20,*) "AXIS2"
do i = 1,kPointCount
   write(20,*) "i map kp",i,mtopKPMap(2,i), kPoints(:,i)
enddo
write(20,*) "AXIS3"
do i = 1,kPointCount
   write(20,*) "i map kp",i,mtopKPMap(3,i), kPoints(:,i)
enddo

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


subroutine generateTetrahedra

   ! Include the modules we need
   use O_Kinds

   implicit none

   ! Passed parameters

   ! Define the local variables used in this subroutine.

   ! The list of kpoints is structured as a regular mesh that matches the
   !   reciprocal space lattice (angles and magnitudes). The distances between
   !   kpoints along each a, b, c axis may be different, but, along any given
   !   axis the distances are constant. One could envision the mesh as a set of
   !   exactly equal parallelepipeds that are stacked to form the reciprocal
   !   space cell. The corners (vertices) of the parallelepipeds are the set
   !   of kpoints.

   ! From that set of kpoints we must generate tetrahedra. The tricky thing is
   !   that the kpoints are generated using three nested loops so that indexing
   !   the vertices of a given parallelepiped requires careful counting of
   !   the position in the one dimensional the kpoint array.

   ! The 
   
   !vertices. it as a set of equal sized Each orthorhombic
   !   parallelpiped consists of eight vertices that we will label M1 through
   !   M8. The tricky thing of course is that the kpoints are created by
   !   looping through a, b, c axis steps (in that order). In other words the

   

end subroutine generateTetrahedra


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

end subroutine cleanUpKPoints


end module O_KPoints
