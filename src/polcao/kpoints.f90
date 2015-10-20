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

   integer :: numKPoints ! The number of kpoints in the system.
   integer :: numHighSymmKP ! Number of high symmetry kpoints that define the
         !   vertices of the path to be taken for the band diagram.
   integer :: numPathKP   ! Number of kpoints that will be used
         !   to create the path between the high symmetry kpoints.
   integer :: isCartesian    ! 1 = yes, 0 = no (Are the high symmetry kpoints
         !   given in cartesian coordinates?)
   real (kind=double), allocatable, dimension (:) :: kPointWeight ! The
         !   weight assigned to each kpoint.
   real (kind=double), allocatable, dimension (:,:) :: kPoints ! The acutal
         !   kpoints.  The first dimension holds the three cartesian xyz
         !   coordinates.  The second dimension is the index over the number
         !   of kpoints.
   real (kind=double), allocatable, dimension (:,:) :: highSymmKP ! These
         !   are the high symmetry kpoints of the requested path.  The first
         !   dimension holds the three space coordinates of the kpoints
         !   to be indexed by the second dimension (numHighSymmKP).
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
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Define the local variables used in this subroutine.
   integer :: i ! Loop index variable.
   integer :: counter ! Dummy variable to read the kpoint index number.

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
   !   coordinates.  Later, after the reciprocal lattice has been initialized
   !   these values will be changed into x,y,z cartesian coordinates.
   do i = 1, numKPoints
      read (15,*)    counter, kPointWeight(i), kPoints(:dim3,i)
      write (20,100) counter, kPointWeight(i), kPoints(:dim3,i)
   enddo
   call flush (20)

   100 format (i5,1x,4f15.8)

end subroutine readKPoints


subroutine readSYBDKPoints(readUnit, writeUnit)

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3
   use O_ReadDataSubs

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Define the local variables used in this subroutine.
   integer :: i ! Loop index variable.

   ! Read the number of highly symmetric kpoints that define the vertices of
   !   the reciprocal space path, the number of kpoints to use on that path,
   !   and whether or not the points are given in cartesian coordinates or not.
   call readData(readUnit,writeUnit,numHighSymmKP,numPathKP,isCartesian,&
                    len('SYBD_INPUT_DATA'),'SYBD_INPUT_DATA')

   write (20,*) 'Number of high symmetry K-Points = ',numHighSymmKP
   write (20,*) 'Number of path K-Points          = ',numPathKP
   write (20,*) '1 = cart, 0 = fract:             = ',isCartesian
   call flush (20)

   ! Allocate space to hold the high symmetry kpoints.
   allocate (highSymmKP(dim3,numHighSymmKP))

   ! Read the coordinates of the high symmetry kpoints.
   do i = 1, numHighSymmKP
      call readData(readUnit,writeUnit,3,highSymmKP(:,i),0,'')
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
   integer :: i,j ! Loop index variables
   integer :: kPointCounter
   integer :: numSegmentKPoints
   real (kind=double) :: averageDelta
   integer, allocatable, dimension (:) :: indexKP
   real (kind=double), dimension (3) :: tempKPoint
   real (kind=double), dimension (3) :: segmentDelta
   real (kind=double), allocatable, dimension (:,:) :: symKPDistVect
   real (kind=double), allocatable, dimension (:)   :: symKPDistMag

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
   allocate (symKPDistVect (3,numHighSymmKP))
   allocate (symKPDistMag  (numHighSymmKP))
   allocate (indexKP       (numHighSymmKP))


   ! Convert the kpoints from fractional coordinates to cartesian if needed.
   if (isCartesian /= 1) then

      do i = 1, numHighSymmKP
         tempKPoint(:) = highSymmKP(1,i) * recipVectors(:,1) + &
                       & highSymmKP(2,i) * recipVectors(:,2) + &
                       & highSymmKP(3,i) * recipVectors(:,3)
         highSymmKP(:,i) = tempKPoint(:)
      enddo
   endif

   ! Identify distances between each consecutive symmetric kpoint on the path.
   symKPDistVect(:,1) = 0.0_double
   symKPDistMag(1)    = 0.0_double
   do i = 2, numHighSymmKP
      symKPDistVect(:,i) = highSymmKP(:,i) - highSymmKP(:,i-1)
      symKPDistMag(i) = symKPDistMag(i-1) + sqrt(sum(symKPDistVect(:,i)**2))
   enddo

   ! Consider the special case where the number of requested path kpoints is
   !   equal to the number of highly symmetric kpoints provided.
   if (numPathKP == numHighSymmKP) then

      ! Copy the highly symmetric kpoint coordinates to the path kpoints.
      do i = 1, numHighSymmKP
         kPoints(:,i)     = highSymmKP(:,i)
         pathKPointMag(i) = symKPDistMag(i)
         indexKP(i) = i

         ! Record the index number of the other high symmetry k points.
         write (20,*) 'HIGH SYMMETRY K POINT INDEX NUMBER :',i
      enddo

      ! Record the value for the kPoint counter
      kPointCounter = numHighSymmKP

   else

      ! Initialize a kpoint counter
      kPointCounter = 1

      ! Initialize the first point.
      kPoints(:,1)     = highSymmKP(:,1)
      pathKPointMag(1) = 0.0_double
      indexKP(1) = 1

      ! Record the index number of the first high symmetry k point.
      write (20,*) 'HIGH SYMMETRY K POINT INDEX NUMBER :',1

      ! Determine the scaling factor for the distance between kpoints in a given
      !   segment between highly symmetric kpoints.
      averageDelta = symKPDistMag(numHighSymmKP) / (numPathKP - 1)

      do i = 2, numHighSymmKP

         ! Determine the number of k points to use for the current segment
         !   between symmetric k points (i) and (i-1).
         numSegmentKPoints = (symKPDistMag(i) - symKPDistMag(i-1))/averageDelta

         ! If there are not going to be any segment kpoints for this segment
         !   then cycle to the next symmetric kpoint.
         if (numSegmentKPoints < 1) cycle

         ! Get the size of the x,y,z deltas for this segment
         segmentDelta(:) = symKPDistVect(:,i) / numSegmentKPoints

         ! Loop to assign the k point positions for this segment
         do j = 1, numSegmentKPoints - 1

            ! Increment the K point counter
            kPointCounter = kPointCounter + 1

            ! Store the position of the next kpoint on this segment by adding
            !   the above determined delta to the last known path kpoint.
            kPoints(:,kPointCounter) = &
                  & kPoints(:,kPointCounter - 1) + segmentDelta(:)

            ! Calculate the magnitude of the vector for the above path KPoint.
            pathKPointMag(kPointCounter) = &
                  & pathKPointMag(kPointCounter - 1) + &
                  & sqrt(sum(segmentDelta(:)**2))
         enddo

         ! Now, record the next high symmetry k point.

         ! Increment the K point counter.
         kPointCounter = kPointCounter + 1

         ! Store the index number of the high symmetry k point.
         indexKP(i) = kPointCounter

         ! Record the index number of the first high symmetry k point.
         write (20,*) 'HIGH SYMMETRY K POINT INDEX NUMBER :',kPointCounter

         ! Store the position of the high symmetry k point in the path array.
         kPoints(:,kPointCounter) = highSymmKP(:,i)

         ! Calculate the magnitude of the vector for the above path KPoint.
         pathKPointMag(kPointCounter) = &
               & pathKPointMag(kPointCounter - 1) + &
               & sqrt(sum(segmentDelta(:)**2))
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
   deallocate (indexKP)

   ! Finish the kpoint path creation.
   call timeStampEnd(21)

end subroutine makePathKPoints


subroutine cleanUpKPoints

   implicit none

   deallocate (kPointWeight)
   deallocate (kPoints)
   deallocate (highSymmKP)

   ! Only allocated in setup.
   if (allocated(phaseFactor)) then
      deallocate (phaseFactor)
   endif

end subroutine cleanUpKPoints


end module O_KPoints
