module O_ExchangeCorrelation

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Exchange Correlation mesh definition.
   integer :: printXCMesh ! 0 = Do not print the mesh; 1 = Print it.
   integer :: numSampleVectors ! How many vectors to use.
   real (kind=double), dimension (2) :: angSampleWeights ! This holds the
         !   weighting factor coefficients that are applied to the vectors.
         !   There are two weighting zones:  The region inside the covalent
         !   radius and the region outside the covalent radius (defined by the
         !   potential type definition).  Presently the default input uses the
         !   same factor for both regions and so the covalent radius has little
         !   meaning or use.  Testing should be applied though.
   real (kind=double), allocatable, dimension (:,:) :: angSampleVectors !
         !   The angular sample vectors are held here.  The first dimension
         !   holds the three space coordinates of the vectors, and the
         !   second dimension is the number of vectors.  Note that it is
         !   assumed that the vectors are already normalized to one in the
         !   input file and that no normalization will occur in the program.
   real (kind=double) :: rSampleIn  ! This is the radius that divides the
         !   inner range of a vector from the middle range of a vector.
   real (kind=double) :: rSampleOut ! This is the radius that divides the
         !   outer range of a vector from the middle range of a vector.
   real (kind=double) :: rSampleSpace ! This is the spacing of the sample
         !   points.

   real (kind=double), allocatable, dimension (:)     :: radialWeight
   real (kind=double), allocatable, dimension (:,:)   :: exchCorrOverlap
   real (kind=double), allocatable, dimension (:,:,:) :: exchRhoOp
   integer :: numRayPoints    ! Number of ray points for a specific potential
                              !   site.  Also used to determine the max.
   integer :: maxNumRayPoints ! Max number of ray points of all the potential
                              !   sites in the system.

   ! A couple of intermediate data sets that are retained from the calculation
   !   of certain mesh parameters so that they don't have to be recalculated
   !   again when the actual mesh and overlap values are calculated.
   real (kind=double), allocatable, dimension (:,:) :: siteSampleMaxRadius !
         !   This stores the maximum distance for each angular sample
         !   vector for each potential site.
   real (kind=double), allocatable, dimension (:,:) :: potLatticeSep ! This
         !   set of vectors stored the difference between the current
         !   potential site and the lattice vector closest to it.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


subroutine getECMeshParameters

   ! Include the modules we need.
   use O_Kinds
   use O_Constants, only: dim3, smallThresh
   use O_Lattice, only: numCellsReal, cellDimsReal, findLatticeVector
   use O_PotSites, only: numPotSites, potSites
   use O_PotTypes, only: minPotAlpha, potTypes
   use O_TimeStamps

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the local loop control variables.
   integer :: i,j,k,l ! Loop variables.  loop1=i, nestedloop2=j ...

   ! Local boundary parameters to be assigned from the passed structures.
   real (kind=double) :: currentCovRadiusHalf

   ! Local locational vectors and distance variables
   real (kind=double), dimension (dim3) :: latticeVector ! This vector points
         ! to the lattice point that is closest to the position vector of each
         ! potential site in the system cell.
   real (kind=double), dimension (dim3) :: potSiteSep ! This vector is
         ! the vector difference of the current potential site's position
         ! vector and the position vector of another potential site in a
         ! replicated cell.
   real (kind=double), dimension (dim3) :: currentAngleVector ! This vector
         ! is the current angular sample vector in the loop over sample vectors.
   real (kind=double), allocatable, dimension (:,:,:) :: exchCorrDist ! This
         ! is the distance between a potential site and all of its replicated
         ! images in neighboring cells.
   real (kind=double) :: siteAngleDot ! This is the dot product of the potential
         ! site position vector, and the angular sample vector of the current
         ! loop iteration.

   ! Local boundary parameters.
   integer :: numRepsNeeded
   real (kind=double) :: maxRadius
   real (kind=double) :: minRadius
   real (kind=double) :: currentMaxRadius
   real (kind=double) :: currentPointRadius

   ! Local variable counters or index reference numbers
   integer, dimension (2) :: currentPotType
   integer :: currentNumRayPoints
   integer :: weightingIndex

   ! Local factors and coefficients updated through loop iterations
   real (kind=double) :: covalentRadiusFactor


   ! Start collecting the exchange correlation mesh point counts.
   call timeStampStart(5)


   ! Determine the number of cell replications that must be considered to find
   !   the exchange correlation matrix parameters.  (i.e. can the real super
   !   lattice be truncated to a smaller set?)
   numRepsNeeded = min (4**dim3,numCellsReal)

   ! Allocate space for local arrays and matrices that will not be used later.
   allocate(exchCorrDist (dim3,numRepsNeeded,numPotSites))

   ! Allocate space for local arrays that will be pointed to later by globally
   !   accessable data structures.
   allocate(potLatticeSep       (dim3,numPotSites))
   allocate(siteSampleMaxRadius (numSampleVectors,numPotSites))

   ! In the event that someone wants the exchange correlation mesh printed
   !   for visualization, then this is where that starts.
   if (printXCMesh == 1) then
      write (200,*) numPotSites
      write (200,*) numSampleVectors
   endif

   ! For every replicated cell, compile a list of the distances from the system
   !   origin to each potential site.
   do i = 1, numPotSites

      ! If requested, print the position of the potential site for plotting of
      !   the XC mesh.
      if (printXCMesh == 1) then
         write (200,*) potSites(i)%cartPos(:)
      endif

      ! Determine if this potential site contributes to the exchange
      !   correlation potential.  If not, then skip this site.
      if (potTypes(potSites(i)%potTypeAssn)%covalentRadius < smallThresh) cycle

      ! Determine the lattice vector closest to the position of this
      !   potential site.
      call findLatticeVector(potSites(i)%cartPos,latticeVector)

      ! Store the difference in distance between the nearest lattice
      !   vector, and the potential site vector.
      potLatticeSep(:,i) = potSites(i)%cartPos(:) - latticeVector

      ! Determine the distance from the origin to each replication of the
      !   current potential site within real super lattice or its truncated
      !   limit.
      do j = 1, numRepsNeeded
         exchCorrDist(:,j,i) = potLatticeSep(:,i) + cellDimsReal(:,j)
      enddo
   enddo

   ! Initialize the current guess for the number of points needed for each ray.
   numRayPoints = 0

   ! Determine the maximum radius that any ray should extend.  Use the smallest
   !   (i.e. most broad) to define this external radius.  Note that all atoms
   !   are treated to initially have the same maximum radius.
   maxRadius = sqrt(rSampleOut/minPotAlpha)


   do i = 1, numPotSites

      ! Get the current potential type from the site number.
      currentPotType(1) = potSites(i)%potTypeAssn

      ! Determine if this potential site contributes to the exchange
      !   correlation potential.  If not, then skip this site.
      if (potTypes(currentPotType(1))%covalentRadius < smallThresh) cycle

      ! Store the current potential site's covalent radius divided by two.
      currentCovRadiusHalf = potTypes(currentPotType(1))%covalentRadius / &
            & 2.0_double

      ! Start an index counter for this potential site for the number of
      !   points for this site's mesh.
      currentNumRayPoints = 0

      ! Determine the minimum radius that any point should be located at for
      !   this potential site.  Use the largest alpha (i.e. the most narrow)
      !   to define this internal radius.
      minRadius = sqrt(rSampleIn / potTypes(currentPotType(1))%alphas&
            & (potTypes(currentPotType(1))%numAlphas))

      ! Begin a loop over each of the angular sample vectors
      do j = 1, numSampleVectors

         ! Store the value of the current angular sample vector.
         currentAngleVector = angSampleVectors(:,j)

         ! Initialize the value for the maximum length of the current ray.
         currentMaxRadius = maxRadius

         ! Begin a loop over all the potential sites in the system cell.
         do k = 1, numPotSites

            ! Get the current potential type from the site number.
            currentPotType(2) = potSites(k)%potTypeAssn

            ! Determine if this potential site contributes to the exchange
            !   correlation potential.  If not, then skip this site. 
            if (potTypes(currentPotType(2))%covalentRadius < &
                  & smallThresh) cycle

            ! Calculate a factor based on the covalent radius of each of the
            !   two potential sites that contribute to the vectors calculated
            !   in the next loop.
            covalentRadiusFactor = potTypes(currentPotType(1))%covalentRadius/&
                  & (potTypes(currentPotType(1))%covalentRadius + &
                  &  potTypes(currentPotType(2))%covalentRadius)

            ! Begin a loop over each of the replicated cells.
            do l = 1, numRepsNeeded

               ! Get the seperation between the current (i) potential site and
               !   the potential site (j) in the current replicated cell.
               !   Including the covalent radius factor.
               potSiteSep(:) = covalentRadiusFactor * &
                  & (exchCorrDist(:,l,k) - exchCorrDist(:,1,i))

               ! Get the dot product of the potential site seperation vector
               !   just obtained, and the current angular sample vector.  This
               !   is the projection of the angular sample vector on the pot
               !   site seperation vector.
               siteAngleDot = sum(potSiteSep(:) * currentAngleVector(:))

               ! If these two vectors are perpendicular then we don't consider
               !   them for helping to estimate the maximum ray length.
               if (siteAngleDot <= smallThresh) cycle

               ! Redetermine the value for the maximum length of the current
               !   ray.
               currentMaxRadius = min(sum(potSiteSep(:)**2)/siteAngleDot,&
                  & currentMaxRadius)
            enddo
         enddo

         ! Save the maximum radius obtained for this angular sample vector
         !   for this potential site.
         siteSampleMaxRadius(j,i) = currentMaxRadius

         ! Initialize the radius of the first sample point on this ray
         currentPointRadius = currentMaxRadius

         ! If requested, print the unit vector describing the current ray
         !   direction and the radius of that ray.  (For XCMesh)
         if (printXCMesh == 1) then
            write (200,*) currentAngleVector(:), currentPointRadius
         endif

         ! Initialize the choice of weighting factor
         weightingIndex = 2

         do while (.true.)

            ! Increment the number of points for this ray.
            currentNumRayPoints = currentNumRayPoints + 1

            ! Determine the zone in which the current point exists.  If it is
            !   in the zone below either test parameter then the weighting
            !   factor is changed.
            if (currentPointRadius < currentCovRadiusHalf .or. &
              & currentPointRadius < minRadius) weightingIndex = 1

            ! If the weighting factor for the current zone is too small then
            !   we exit the loop.
            if (angSampleWeights(weightingIndex) < smallThresh) exit

            ! If the radius is less than the minRadius then we record one
            !   more point and exit.  I am not sure how to make these three
            !   if statements work better or be more understandable.  I will
            !   have to spend some time to puzzle it out.
            if (currentPointRadius < minRadius) then
               currentNumRayPoints = currentNumRayPoints + 1
               exit
            endif

            ! Multiply the radius by a factor to reduce its distance.
            currentPointRadius = currentPointRadius * rSampleSpace

            ! If requested, print the radius of the current point on the ray.
            if (printXCMesh == 1) then
               write (200,*) currentPointRadius
            endif
         enddo

         ! If the request to print the XCMesh is true, then we are done doing
         !   it for this ray and we print "END" to signify that.
         if (printXCMesh == 1) then
            write (200,*) "END"
         endif
      enddo

      ! Determine the maximum number of ray points required so far.  This will
      !   later determine the dimension of the exchange correlation mesh
      !   data structure.
      numRayPoints = max (currentNumRayPoints,numRayPoints)

      ! Record that this loop has finished
      if (mod(i,10) .eq. 0) then
         write (20,ADVANCE="NO",FMT="(a1)") "|"
      else
         write (20,ADVANCE="NO",FMT="(a1)") "."
      endif
      if (mod(i,50) .eq. 0) then
         write (20,*) " ",i
      endif
      call flush (20)
   enddo

   ! At this point, the numRayPoints is the maximum possible number of points
   !   for any ray and we can rename it for use as the dimension of the
   !   exchange correlation data structures.
   maxNumRayPoints = numRayPoints

   ! Deallocate arrays and matrices that will not be used later.
   deallocate(exchCorrDist)

   ! Log the date and time we end.
   call timeStampEnd(5)

end subroutine getECMeshParameters


subroutine makeECMeshAndOverlap

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3, smallThresh
   use O_Lattice, only: logElecThresh, numCellsReal, cellSizesReal, &
         & cellDimsReal, findLatticeVector
   use O_PotSites, only: potSites, numPotSites
   use O_PotTypes, only: potTypes, maxNumPotAlphas
   use O_Potential, only: potDim, GGA
   use O_TimeStamps

   ! Import the necessary HDF modules
   use HDF5
   use O_SetupExchCorrHDF5, only: numPoints_did, numPoints, radialWeight_did, &
         & points, exchRhoOp_did, potPoints, exchCorrOverlap_did, potPot

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the local loop control variables.
   integer :: i,j,k,l,m,n ! Loop variables.  loop1=i, nestedloop2=j ...
   integer :: hdferr

   ! Local boundary parameters to be assigned from the passed structures.
   real (kind=double) :: rSampleSpaceSqrt
   real (kind=double), allocatable, dimension (:) :: currentPotAlphas

   ! Local locational vectors and distance variables
   real (kind=double), dimension (dim3) :: latticeVector ! This vector points
         ! to the lattice point that is closest to the position vector of each
         ! potential site in the system cell.
   real (kind=double), dimension (dim3) :: currentPotLatticeSep ! This stores
         ! the above value for the current loop iteration.
   real (kind=double), dimension (dim3) :: currentAngleVector ! This vector
         ! is the current angular sample vector in the loop over sample vectors.
   real (kind=double), allocatable, dimension (:,:) :: exchangePointRadius !
         ! This is the stored radius vector from each potential site to each of
         ! its mesh points.  The dimensions will be dim3, maxNumRayPoints.
   real (kind=double), dimension (dim3) :: exchangePointSiteSep ! This
         ! holds a vector indicating the seperation between the current
         ! exchange point and the current potential site.
   real (kind=double), dimension (dim3) :: currentRayPoint ! This points to
         ! an exchange point on a ray for each iteration of the point, and ray
         ! loops.
   real (kind=double), dimension (dim3) :: latticeOffset ! This vector gives
         ! the relation between each potential site and each exchange point on
         ! a lattice.

   ! Local boundary parameters.
   real (kind=double) :: minRadius
   real (kind=double) :: currentPointRadius

   ! Local variable counters and indices
   integer :: currentNumRayPoints
   integer :: weightingIndex
   integer :: initPotAlphaIndex
   integer :: finPotAlphaIndex
   integer :: currentPotAlphaIndex
   integer :: currentNumPotAlphas
   integer, dimension (2) :: currentPotType


   ! Local factors and coefficients updated through loop iterations
   real (kind=double) :: covalentRadiusHalf
   real (kind=double) :: currentNegligLimit
   real (kind=double) :: weightFactor
   real (kind=double) :: endPointFactor
   real (kind=double) :: latticeOffsetMagSqrd
   real (kind=double) :: currentNegligLimitPlus
   real (kind=double) :: radialmagnitude
   real (kind=double) :: currentAlphaMagnitude

   ! These are the different distances that will be needed for the
   ! first and second derivatives of rho.
   real (kind=double) :: xdistance
   real (kind=double) :: ydistance
   real (kind=double) :: zdistance
   real (kind=double) :: xxdistance
   real (kind=double) :: xydistance
   real (kind=double) :: xzdistance
   real (kind=double) :: yydistance
   real (kind=double) :: yzdistance
   real (kind=double) :: zzdistance

   ! We need to know the size of the last dimension in the exchRhoOp matrix.
   !   This value is determined by whether or not we are doing an LDA (1) or a
   !   GGA (10) based calculation.
   integer :: numOpValues

   ! Start making the exchange correlation matrices and mesh points.
   call timeStampStart(7)

   ! Initialize variables to avoid compiler warnings.
   xdistance = 0.0_double
   ydistance = 0.0_double
   zdistance = 0.0_double
   xxdistance = 0.0_double
   xydistance = 0.0_double
   xzdistance = 0.0_double
   yydistance = 0.0_double
   yzdistance = 0.0_double
   zzdistance = 0.0_double

   ! Pull variables out of the large data structures and assign local values
   !   to them.
   rSampleSpaceSqrt = sqrt(rSampleSpace)

   ! Allocate space for matrices that will be pointed to by global data
   !   structures later.
   allocate (radialWeight    (maxNumRayPoints))
   ! When GGA=0 only space for the rho operator is allocated. When GGA=1 space
   !   is also allocated for first and second derivatives of the rho operator.
   if (GGA == 0) then ! Doing LDA
      numOpValues = 1
   allocate (exchRhoOp       (potDim,maxNumRayPoints,1))
   else ! Doing GGA
      numOpValues = 10
   allocate (exchRhoOp       (potDim,maxNumRayPoints,10))
   endif
   allocate (exchCorrOverlap (potDim,potDim))


   ! Allocate space for matrices and arrays that are used only locally
   allocate (currentPotAlphas    (maxNumPotAlphas))
   allocate (exchangePointRadius (dim3,maxNumRayPoints))

   ! Initialize the exchCorrOverlap matrix to zero since it will later be
   !   created through cumulative summation.
   exchCorrOverlap(:,:) = 0.0_double

   ! Initialize other variables.
   currentNegligLimit  = 0.0_double
   currentNumPotAlphas = 0

   do i = 1, numPotSites

      ! Initialize the radialWeight to zero for each potential iteration.
      radialWeight(:) = 0.0_double

      ! Get the current potential type from the site number.
      currentPotType(1) = potSites(i)%potTypeAssn

      ! Determine if this potential site contributes to the exchange
      !   correlation potential.  If not, then skip this site.
      if (potTypes(currentPotType(1))%covalentRadius < smallThresh) cycle

      ! Store the value of half the current potential site's covalent radius.
      covalentRadiusHalf = potTypes(currentPotType(1))%covalentRadius / &
            & 2.0_double

      ! Initialize the exchange Rho Operator matrix for this potential site.
      exchRhoOp(:,:,:) = 0.0_double

      ! Initialize a counter to index the total number of points.  Note that
      !   this should work out to the same number that was calculated in the
      !   previous subroutine.  The previous subroutine just obtained this
      !   dimension variable for dynamic allocation, along with some other
      !   minor matrices.
      numRayPoints = 0

      ! Initialize the weighting factor
      weightFactor = (1.0_double - rSampleSpace) / rSampleSpaceSqrt

      ! Record the potential lattice seperation for this point.
      currentPotLatticeSep(:) = potLatticeSep(:,i)

      ! Determine the minimum radius that any point should be located at for
      !   this potential site.
      minRadius = sqrt(rSampleIn / potTypes(currentPotType(1))% &
            & alphas(potTypes(currentPotType(1))%numAlphas))


      do j = 1, numSampleVectors

         ! Store the value of the current angular sample vector.
         currentAngleVector(:) = angSampleVectors(:,j)

         ! Store a local variable of the current maximum distance for the
         !   current ray of the current potential site.
         currentPointRadius = siteSampleMaxRadius(j,i)

         ! Store the current weighting index number
         weightingIndex = 2

         ! Initialize a counter to index the total number of point ON THIS RAY.
         currentNumRayPoints = 0

         ! Initialize a factor that will reduce the weight of the first point
         !   treated (that is the point farthest from the potential site) by
         !   a factor of two.  After the first point is processed, this factor
         !   is changed to 1.0 so it has no later effect.
         endPointFactor = 0.5_double


         do while (.true.)

            ! Increment the total number of ray points, and the number of
            !   ray points for this ray.
            numRayPoints = numRayPoints + 1
            currentNumRayPoints = currentNumRayPoints + 1

            ! Store the position of this ray point plus the exchange position
            !   constant.
            exchangePointRadius(:,numRayPoints) = currentPointRadius * &
                  & currentAngleVector(:) + currentPotLatticeSep(:)

            ! Store the weighting factor for this point.
            radialWeight(numRayPoints) = angSampleWeights(weightingIndex) * &
                  & weightFactor * endPointFactor * currentPointRadius**3

            ! Switch the value of the endPointFactor to 1.0 so it has no more
            !   effect on the weighting of points.
            endPointFactor = 1.0_double

            ! Determine the zone in which the current point exists.  If it is
            !   in the zone below either test parameter then the weighting
            !   factor is changed.
            if (currentPointRadius < covalentRadiusHalf .or. &
                  & currentPointRadius < minRadius) weightingIndex = 1

            ! If the weighting factor for the current zone is too small then
            !   we exit the loop.
            if (angSampleWeights(weightingIndex) < smallThresh) exit

            ! If the radius is less than the minRadius then we record one
            !   more point and exit.  I am not sure how to make these three
            !   'if' statements work better or be more understandable.  I will
            !   have to spend some time to puzzle it out later.
            if (currentPointRadius < minRadius) then

               ! Increment the total number of ray points, and the number of
               !   ray points for this ray.
               numRayPoints = numRayPoints + 1
               currentNumRayPoints = currentNumRayPoints + 1

               ! Update the position of the current point radius as usual,
               !   except with a slightly larger spacing due to the sqrt.
               currentPointRadius = currentPointRadius * rSampleSpaceSqrt

               ! Store the position of this ray point plus the exchange position
               !   constant.  Note the addition of a factor 0.75 in front of
               !   the vector definition.
               exchangePointRadius(:,numRayPoints) = 0.75_double * &
                     & currentPointRadius * currentAngleVector(:) + &
                     & currentPotLatticeSep(:)

               ! Store the weighting factor for this point in a bit different
               !   way than before.
               radialWeight(numRayPoints) = angSampleWeights(1) * &
                     & (currentPointRadius**3) / 3.0_double

               exit
            endif

            ! Multiply the radius by a spacing factor to reduce its distance.
            currentPointRadius = currentPointRadius * rSampleSpace
         enddo


         ! Initiate a loop over the points for this ray.
         do k = numRayPoints - currentNumRayPoints + 1, numRayPoints

            ! Store the position vector of the current point on the current ray.
            currentRayPoint(:) = exchangePointRadius(:,k)

            ! Initialize two cumulative index counters on the range of potential
            !   alphas being considered out of all the potential alphas in the
            !   whole system.
            initPotAlphaIndex = 0
            finPotAlphaIndex  = 0

            do l = 1, numPotSites

               ! Get the current potential type from the site number.
               currentPotType(2) = potSites(l)%potTypeAssn

               ! Get the seperation vector from the current exchange point on
               !   the current ray to the current potential site.
               exchangePointSiteSep(:) = currentRayPoint(:) - &
                     & potSites(l)%cartPos(:)

               ! If this potential site is the first in a set of equivalent
               !   potential sites then we update certain parameters and
               !   indices.
               if (potSites(l)%firstPotType == 1) then

                  ! Store the number of potential alphas for this site
                  currentNumPotAlphas = potTypes(currentPotType(2))%numAlphas

                  ! Store local values for the current potential alphas.
                  currentPotAlphas(:currentNumPotAlphas) = &
                        & potTypes(currentPotType(2))%alphas(:&
                        & currentNumPotAlphas)

                  ! Determine the maximum radial distance beyond which values
                  !   are considered negligible.
                  currentNegligLimit = logElecThresh/currentPotAlphas(1)

                  ! Update the alpha indices
                  initPotAlphaIndex = finPotAlphaIndex
                  finPotAlphaIndex  = initPotAlphaIndex + currentNumPotAlphas
               endif

               ! Get the lattice point vector that is closest to the current
               !   exchangePointSiteSep vector.
               call findLatticeVector(exchangePointSiteSep,latticeVector)

               ! Determine the vector which is the exchange point vector - 
               !   the position vector of the current potential site - 
               !   the lattice point vector for the point closest to the
               !   difference between the exhcange point vector and the
               !   potential site vector.  This is used to define the relation
               !   between the potential site and the exchange point over all
               !   the lattice sites in a later loop.
               latticeOffset(:) = exchangePointSiteSep(:) - latticeVector(:)

               ! Get the square of the magnitude of the latticeOffset vector
               latticeOffsetMagSqrd = sum(latticeOffset(:)**2)

               ! Determine the distance beyond which a replicated lattice cell
               !   will not be considered.
               currentNegligLimitPlus = currentNegligLimit + &
                     & latticeOffsetMagSqrd + 2.0_double * &
                     & sqrt(currentNegligLimit * latticeOffsetMagSqrd)

               ! Begin a loop over all the replicated cells.
               do m = 1, numCellsReal

                  ! Once the cells are beyond the required range for this
                  !   potential site we exit the cell loop.
                  if (cellSizesReal(m) > currentNegligLimitPlus) exit

                  ! Determine the magnitude of the radial component.
                  radialMagnitude = sum((latticeOffset(:)-cellDimsReal(:,m))**2)

                  ! If this particular potential site is beyond the required
                  !   range then we cycle on this loop.
                  if (radialMagnitude > currentNegligLimit) cycle

                  ! These are the distances used when calculating the first
                  ! derivatives (xdistance,ydistance,zdistance) and the second
                  ! derivatives (xxdistance,xydistance,xzdistance,yydistance
                  ! yzdistance,zzdistane) for the exchange rho op value as
                  ! a gaussian function.
                  if (GGA == 1) then
                     xdistance = latticeOffset(1)-cellDimsReal(1,m)
                     ydistance = latticeOffset(2)-cellDimsReal(2,m)
                     zdistance = latticeOffset(3)-cellDimsReal(3,m)
                     xxdistance = xdistance * xdistance
                     xydistance = xdistance * ydistance
                     xzdistance = xdistance * zdistance
                     yydistance = ydistance * ydistance
                     yzdistance = ydistance * zdistance
                     zzdistance = zdistance * zdistance
                  endif

                  do n = 1, currentNumPotAlphas

                     ! Determin the range magnitude for the current alpha
                     currentAlphaMagnitude = currentPotAlphas(n) * &
                           & radialMagnitude

                     ! If this range is beyond the negligability limit then
                     !   we exit this loop. and continue with the next cell.
                     if (currentAlphaMagnitude > logElecThresh) exit

                     ! Increment the index for which alpha we are dealing with
                     !   out of all alphas in the whole system
                     currentPotAlphaIndex = initPotAlphaIndex + n

                     ! Save the exchange rho op value for this point and alpha.
                     !   If we are doing a GGA calculation then we also need to
                     !   save the first (x,y,z) derivatives and the second (xx,
                     !   xy, xz, yy, yz, zz) derivatives.
                     if (GGA == 0) then
                        ! Save just this point.
                        exchRhoOp(currentPotAlphaIndex,k,1) = &
                              & exchRhoOp(currentPotAlphaIndex,k,1) + &
                              & exp(-currentAlphaMagnitude)
                     else
                        exchRhoOp(currentPotAlphaIndex,k,1) = &
                              & exchRhoOp(currentPotAlphaIndex,k,1) + &
                              & exp(-currentAlphaMagnitude)
                        exchRhoOp(currentPotAlphaIndex,k,2) = &
                              & currentPotAlphas(n) * &
                              & exchRhoOp(currentPotAlphaIndex,k,2) + &
                              & xdistance*exp(-currentAlphaMagnitude)
                        exchRhoOp(currentPotAlphaIndex,k,3) = &
                              & currentPotAlphas(n) * &
                              & exchRhoOp(currentPotAlphaIndex,k,3) + &
                              & ydistance*exp(-currentAlphaMagnitude)
                        exchRhoOp(currentPotAlphaIndex,k,4) = &
                              & currentPotAlphas(n) * &
                              & exchRhoOp(currentPotAlphaIndex,k,4) + &
                              & zdistance*exp(-currentAlphaMagnitude)
                        exchRhoOp(currentPotAlphaIndex,k,5) = &
                              & (currentPotAlphas(n)**2) * &
                              & exchRhoOp(currentPotAlphaIndex,k,5) + &
                              & xxdistance*exp(-currentAlphaMagnitude)
                        exchRhoOp(currentPotAlphaIndex,k,6) = &
                              & (currentPotAlphas(n)**2) * &
                              & exchRhoOp(currentPotAlphaIndex,k,6) + &
                              & xydistance*exp(-currentAlphaMagnitude)
                        exchRhoOp(currentPotAlphaIndex,k,7) = &
                              & (currentPotAlphas(n)**2) * &
                              & exchRhoOp(currentPotAlphaIndex,k,7) + &
                              & xzdistance*exp(-currentAlphaMagnitude)
                        exchRhoOp(currentPotAlphaIndex,k,8) = &
                              & (currentPotAlphas(n)**2) * &
                              & exchRhoOp(currentPotAlphaIndex,k,8) + &
                              & yydistance*exp(-currentAlphaMagnitude)
                        exchRhoOp(currentPotAlphaIndex,k,9) = &
                              & (currentPotAlphas(n)**2) * &
                              & exchRhoOp(currentPotAlphaIndex,k,9) + &
                              & yzdistance*exp(-currentAlphaMagnitude)
                        exchRhoOp(currentPotAlphaIndex,k,10) = &
                              & (currentPotAlphas(n)**2) * &
                              & exchRhoOp(currentPotAlphaIndex,k,10) + &
                              & zzdistance*exp(-currentAlphaMagnitude)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo


      ! Accumulate values in the exchange correlation overlap matrix here
      do j = 1, numRayPoints
         do k = 1, potDim
            exchCorrOverlap(1:k,k) = exchCorrOverlap(1:k,k) + &
            exchRhoOp(k,j,1) * exchRhoOp(1:k,j,1) * radialWeight(j)
         enddo
      enddo


      ! Save values calculated from this potential site loop.

      ! Write the values for the radialWeight and exchRhoOp in HDF5 format.
      call h5dwrite_f(numPoints_did(i),H5T_NATIVE_INTEGER, &
            & numRayPoints,numPoints,hdferr)
      if (hdferr /= 0) stop 'Failed to write num ray points'
      call h5dwrite_f(radialWeight_did(i),H5T_NATIVE_DOUBLE, &
            & radialWeight(:),points,hdferr)
      if (hdferr /= 0) stop 'Failed to write radial weights'
      call h5dwrite_f(exchRhoOp_did(i),H5T_NATIVE_DOUBLE, &
            & exchRhoOp(:,:,:),potPoints,hdferr)
      if (hdferr /= 0) stop 'Failed to write exchange rho operator'


      ! Record the progress so far.
      if (mod(i,10) .eq. 0) then
         write (20,ADVANCE="NO",FMT="(a1)") "|"
      else
         write (20,ADVANCE="NO",FMT="(a1)") "."
      endif
      if (mod(i,50) .eq. 0) then
         write (20,*) " ",i
      endif
      call flush (20)
   enddo

   ! Save to disk the overlap that was accumulated through each loop

   ! Write the exchCorrOverlap to disk in HDF5 format.
   call h5dwrite_f(exchCorrOverlap_did,H5T_NATIVE_DOUBLE, &
         & exchCorrOverlap(:,:),potPot,hdferr)
   if (hdferr /= 0) stop 'Failed to write exch corr overlap'

   ! Deallocate space for matrices and arrays that are used only locally.
   deallocate (currentPotAlphas)
   deallocate (exchangePointRadius)
   deallocate (potLatticeSep)
   deallocate (siteSampleMaxRadius)
   deallocate (radialWeight)
   deallocate (exchRhoOp)
   deallocate (exchCorrOverlap)

   ! Log the date and time we end.
   call timeStampEnd(7)

end subroutine makeECMeshAndOverlap


subroutine readExchCorrMeshParameters(readUnit, writeUnit)

   ! Import necessary modules.
   use O_Kinds
   use O_ReadDataSubs

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: readUnit ! The unit number of the file from which we
                                   ! are reading.
   integer, intent(in) :: writeUnit ! The unit number of the file to which we
                                    ! are writing.

   ! Define local variables.
   real (kind=double), dimension(3) :: tempData

   ! Read the flag requesting that the XC mesh be printed or not.
   call readData(readUnit,writeUnit,printXCMesh,len('PRINT_XC_MESH'),&
         & 'PRINT_XC_MESH')

   ! Read the number of angular sample vectors.
   call readData(readUnit,writeUnit,numSampleVectors,&
         & len('NUM_ANGULAR_SAMPLE_VECTORS'),'NUM_ANGULAR_SAMPLE_VECTORS')

   ! Read the inner and outer weighting factors.
   call readData(readUnit,writeUnit,2,angSampleWeights(:),&
         & len('WTIN_WTOUT'),'WTIN_WTOUT')

   ! Read the parameters for sampling.
   call readData(readUnit,writeUnit,3,tempData(:),&
         & len('RADIAL_SAMPLE-IN_OUT_SPACING'),'RADIAL_SAMPLE-IN_OUT_SPACING')
   rSampleIn    = tempData(1)
   rSampleOut   = tempData(2)
   rSampleSpace = tempData(3)

end subroutine readExchCorrMeshParameters



! Use the generalized spiral set on S2 as described by E.B. Saff and A.B.J.
!   Kuijlaars, "Distributing Many Points on a Sphere", The Mathematical
!   Intelligencer, Vol. 19, #1, 5-11, (1997).  (Last page of text.)
subroutine makeSampleVectors

   ! Import necessary modules.
   use O_Kinds
   use O_Constants, only: pi, dim3

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define local variables.
   integer :: k
   real (kind=double), allocatable, dimension (:) :: thetaK;
   real (kind=double), allocatable, dimension (:) :: phiK;
   real (kind=double), allocatable, dimension (:) :: hK;

   ! Allocate space to hold the angular sample vectors.
   allocate (angSampleVectors(dim3,numSampleVectors))

   ! Allocate space to hold thetaK, phiK and hK.
   allocate (thetaK(numSampleVectors))
   allocate (phiK(numSampleVectors))
   allocate (hK(numSampleVectors))

   ! Initialize the values of hK.
   do k = 1, numSampleVectors
      hK(k) = -1.0_double + 2.0_double * &
            & (k-1.0_double) / (numSampleVectors-1.0_double)
   enddo

   ! Obtain the values of thetaK.
   do k = 1, numSampleVectors
      thetaK(k) = acos(hK(k))
   enddo

   ! Obtain the values of phiK.
   phiK(1) = 0.0_double
   phiK(numSampleVectors) = 0.0_double
   do k = 2, numSampleVectors-1
      phiK(k) = phiK(k-1) + 3.6_double/sqrt(real(numSampleVectors,double)) * &
            & 1.0_double / sqrt(1.0_double - hK(k)*hK(k))
      do while (phiK(k) > 2.0_double*pi)
         phiK(k) = phiK(k) - 2.0_double*pi
      enddo
   enddo

   ! Convert the spherical coordinates to xyz vectors and store them.
   do k = 1, numSampleVectors
      angSampleVectors(1,k) = sin(thetaK(k)) * cos(phiK(k))
      angSampleVectors(2,k) = sin(thetaK(k)) * sin(phiK(k))
      angSampleVectors(3,k) = cos(thetaK(k))
   enddo

   ! Deallocate the space that held thetaK, phiK and hK.
   deallocate (thetaK)
   deallocate (phiK)
   deallocate (hK)

end subroutine makeSampleVectors


subroutine getMaxNumRayPoints(numPoints_did,numPoints)

   ! Import necessary library modules.
   use HDF5

   ! Import necessary object modules.
   use O_PotSites, only: numPotSites

   ! Make sure no funny variables are defined.
   implicit none

   ! Define passed dummy variables.
   integer(hid_t),   dimension (:) :: numPoints_did
   integer(hsize_t), dimension (1) :: numPoints

   ! Define local variables.
   integer :: i
   integer :: hdferr

   maxNumRayPoints = 0
   do i = 1, numPotSites

      ! Num ray points for this potSite.
      call h5dread_f (numPoints_did(i),H5T_NATIVE_INTEGER,&
            & numRayPoints,numPoints,hdferr)
      if (hdferr /= 0) stop 'Failed to read numRay Points'

      ! Max of all num ray points read so far.
      maxNumRayPoints = max(maxNumRayPoints,numRayPoints)
   enddo

end subroutine getMaxNumRayPoints


subroutine cleanUpExchCorr

   implicit none

   deallocate (angSampleVectors)

end subroutine cleanUpExchCorr


end module O_ExchangeCorrelation
