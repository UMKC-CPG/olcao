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
         !   If 0 then the histogram method will be used.
         !   If 1 then LAT as described in will be used.
         !   If 2 then LAT as described in will be used.
         !   If 3 then LAT as described in will be used.
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
   integer, allocatable, dimension (:) :: numHighSymKP ! Number of high
         !   symmetry kpoints that define the vertices of each path to be taken
         !   for the band diagram.
   real (kind=double) :: minKPointDensity ! The number of kpoints along each
         !   a,b,c axis divided by the magnitude of that axis must be greater
         !   than or equal to the value of minKPointDensity. This is only used
         !   if the kPointStyleCode equals 2.
   real (kind=double), dimension (dim3) :: kPointShift ! Fractional amount to
         !   shift the kpoint mesh by along each of the a,b,c axes. This is
         !   only used if the kPointStyleCode equals 1 or 2.
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

   ! Read the kpoint definition style code.
   call readData(readUnit,writeUnit,kPointStyleCode,len('KPOINT_STYLE_CODE'),&
         & 'KPOINT_STYLE_CODE')

   ! Read the code that defines the type of integration method to be used.
   call readData(readUnit,writeUnit,kPointIntgCode,len('KPOINT_INTG_CODE'),&
         & 'KPOINT_INTG_CODE')

   ! Depending on kPointStyleCode we will read the associated relevant data.
   if (kPointStyleCode == 0) then ! Read an explicit list of kpoints.

      ! Read the number of kpoints.
      call readData(readUnit,writeUnit,numKPoints,len('NUM_BLOCH_VECTORS'),&
            & 'NUM_BLOCH_VECTORS')

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
         read (15,*)    counter, kPointWeight(i), kPoints(1:dim3,i)
         write (20,100) counter, kPointWeight(i), kPoints(1:dim3,i)
      enddo
      call flush (20)

   elseif (KPointStyleCode == 1) then
      ! Read axial numbers of kpoints and a shift.

      ! Read the number of kpoints along each a,b,c axis.
      call readData(readUnit,writeUnit,3,numAxialKPoints,len('NUM_KP_A_B_C'),&
            & 'NUM_KP_A_B_C')

      ! Read the fractional distance that the uniform mesh should be shifted
      !   away from the origin along each a,b,c axis.
      call readData(readUnit,writeUnit,3,kPointShift,len('KP_SHIFT_A_B_C'),&
            & 'KP_SHIFT_A_B_C')

      ! Expand the requested kpoints into explicit kpoints.
      !call computeKPointMesh

   elseif (KPointStyleCode == 2) then
      ! Read fractional min. kp sep. and a shift.

      ! Read the minimum line density required of kpoints along all axes.
      call readData(readUnit,writeUnit,minKPointDensity,&
            & len('MIN_KP_LINE_DENSITY'),'MIN_KP_LINE_DENSITY')

      ! Read the fractional distance that the uniform mesh should be shifted
      !   away from the origin along each a,b,c axis.
      call readData(readUnit,writeUnit,3,kPointShift,len('KP_SHIFT_A_B_C'),&
            & 'KP_SHIFT_A_B_C')
   endif


   ! If a LAT method is requested, then we need to make the tetrahedra. If
   !   the KPointIntgCode == 0, then we don't have anything additional to do.
   call generateTetrahedra
!   if (KPointIntgCode == 1) then
!      call generateTetWholeBZ
!   elseif (KPointIntgCode == 2) then
!   endif

   100 format (i5,1x,4f15.8)

end subroutine readKPoints


subroutine readSYBDKPoints(readUnit, writeUnit)

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3
   use O_ReadDataSubs

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
   write (20,*) 'Number of discrete paths         = ',numPaths
   write (20,*) 'Number of path K-Points          = ',numPathKP
   write (20,*) '1 = cart, 0 = fract:             = ',isCartesian
   call flush (20)

   ! Allocate space to hold the number of high symmetry kpoints that will
   !   exist for each path. E.g., there may be 6 kpoints that define the first
   !   path and only 2 kpoints that define the second path.
   allocate (numHighSymKP(numPaths))

   ! Read the number of high symmetry kpoints that define the vertices for
   !   each path.
   call readData(readUnit,writeUnit,numPaths,numHighSymKP(:),0,'')

   ! Compute the total number of high symmetry kpoints over all paths.
   numTotalHighSymKP = sum(numHighSymKP)

   write (20,*) 'Number of high symmetry K-Points = ',numTotalHighSymKP
   call flush (20)
   
   ! Allocate space to hold the high symmetry kpoints.
   allocate (highSymKP(dim3,maxval(numHighSymKP),numPaths))

   ! Read the high symmetry kpoints for each path.
   do i = 1, numPaths
   
      ! Read the coordinates of the high symmetry kpoints.
      do j = 1, numHighSymKP(i)
         call readData(readUnit,writeUnit,3,highSymKP(:,j,i),0,'')
      enddo
   enddo
end subroutine readSYBDKPoints


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

   ! Decallocate the kpoints and weights defined in fort.15.
   deallocate (kPointWeight)
   deallocate (kPoints)

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
   deallocate (highSymKP)

   ! Only allocated in setup.
   if (allocated(phaseFactor)) then
      deallocate (phaseFactor)
   endif

end subroutine cleanUpKPoints


end module O_KPoints
