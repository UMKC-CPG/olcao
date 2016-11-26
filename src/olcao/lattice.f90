module O_Lattice  ! Lattice and superlattice object.

   ! Import necessary modules.
   use O_Kinds
   use O_Constants, only: dim3

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public
   private :: getNumCells, makeLattice, cellBorderCheck, createUnitVectors,&
            & dotProductMatrix, intersectionTest, minimalBorderCheck, &
            & minimalBorderRecord, sharedVolumeCheck, defineOctant, &
            & combineWignerAndOctants, compareOctantsAndWignerCells

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Real space cell information.
   real (kind=double), dimension (dim3,dim3) :: realVectors ! Cell vectors
         !   of the real space lattice for a, b, c, given in x, y, z
         !   components.  The first index is x,y,z and the second index is for
         !   a,b,c.
   real (kind=double), dimension (dim3,dim3) :: realUnitVectors ! Unit vectors
         !   of the real space lattice for a, b, c, given in x, y, z
         !   components.  The first index is x,y,z and the second index is for
         !   a,b,c.
   real (kind=double), dimension (dim3,dim3) :: realUnitNormal ! Unit vectors
         !   that are normal to the bc, ac, ab planes.  They define the
         !   orientation of each plane.
   real (kind=double), dimension (dim3) :: realMag ! Magnitudes of the real
         !   space cell (a,b,c).
   real (kind=double), dimension (dim3) :: realAngles ! Angles of the real
         !   space cell alpha = angle between b,c; beta = angle between a,c;
         !   gamma = angle between a,b.
   real (kind=double), dimension (dim3) :: realPlaneAngles ! Angles of the real
         !   space cell that are between the realUnitNormal vectors and the
         !   corresponding realUnitVector.  Angles are (bc,a ; ac,b ; ab,c)
         !   where the first is the plane defined by the normal and the second
         !   is the lattice unit vector.
   integer, dimension(dim3) :: primRepsReal ! The number of repetitions
         !   needed for each vector direction (a, b, c) to enclose the
         !   negligability limit.
   integer :: numCellsReal ! Number of combinations of primitive vectors for
         !   which the size of the super lattice contains the negligability
         !   limit.  This is basically the product of the values in the
         !   repetitions array.
   real (kind=double) :: negligLimitReal ! This is the radius beyond which
         !   integrations can be neglected.  This sphere should fit inside
         !   the superlattice rather tightly.
   real (kind=double) :: realCellVolume ! Just what the name says.
   real (kind=double), allocatable, dimension (:) :: cellSizesReal ! The sum of
         !   the squares of the dimensions for each primitive lattice
         !   combination used to form the superlattice.
   real (kind=double), allocatable, dimension (:,:) :: cellDimsReal ! The
         !   dimensions for each primitive lattice combination used to form
         !   the superlattice.  The first dimension hold the three space
         !   coordinates.  The second dimension is the index of numCells.

   ! Real space uniform mesh information.
   integer :: numTotalMeshPoints
   integer, dimension (3) :: numMeshPoints
   real (kind=double), dimension(3) :: realFractStrideLength
   real (kind=double), dimension(3) :: realFractCrossArea
   real (kind=double), allocatable, dimension (:,:,:) :: meshValues

   ! Reciprocal space cell information.
   real (kind=double), dimension (dim3,dim3) :: recipVectors ! These
         !   are the unit vectors of the reciprocal lattice a', b', c'
         !   given in x, y, z components.
   integer, dimension(dim3) :: primRepsRecip ! The number of repetitions
         !   needed for each vector direction (a, b, c) to enclose the
         !   negligability limit.
   integer :: numCellsRecip ! Number of combinations of primitive vectors for
         !   which the size of the super lattice contains the negligability
         !   limit.  This is basically the product of the values in the
         !   repetitions array.
   real (kind=double), allocatable, dimension (:) :: cellSizesRecip ! The sum of
         !   the squares of the dimensions for each primitive lattice
         !   combination used to form the superlattice.
   real (kind=double), allocatable, dimension (:,:) :: cellDimsRecip ! The
         !   dimensions for each primitive lattice combination used to form
         !   the superlattice.  The first dimension hold the three space
         !   coordinates.  The second dimension is the index of numCells.
   real (kind=double) :: negligLimitRecip ! This is the radius beyond which
         !   integrations can be neglected.  This sphere should fit inside
         !   the superlattice rather tightly.
   real (kind=double) :: recipCellVolume ! Just what the name says.


   ! Data to help find the super lattice vector closest to a given vector.
   integer, dimension (8) :: numWignerCellsNeeded ! This counts the number
         !   of Wigner-Seitz cells needed to cover a given octant that was
         !   created when the given cell was centered at the origin and
         !   divided into 8 pieces.
   integer, allocatable, dimension (:,:) :: overlappingWignerCells ! For
         !   each of the octants this will record which Wigner-Sietz cells
         !   are needed to completely cover it.  Each of the lattice points
         !   that the Wigner-Seitz cells is centered at is used later in the
         !   determination of the closest lattice vector to any arbitrary
         !   vector.
   real (kind=double), dimension (dim3,dim3) :: invRealVectors ! This matrix
         !   is just the reciprocal matrix / 2 / Pi (aka the inverse of the
         !   real space lattice vectors matrix.  It is useful to have this
         !   matrix precalculated to reduce computation.  The matrix is
         !   used to extract factors to find the closest lattice vector to
         !   any arbitrary vector.

   ! Define values that determine interaction cutoffs for this lattice.
   real (kind=double) :: logBasisFnThresh ! Threshold for basis function
         !   interaction given in terms of the exponent (log) of the number.
   real (kind=double) :: logElecThresh ! Threshold for electroStatic
         !   interaction given in terms of the exponent (log) of the number.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

! This subroutine will read the real space lattice vectors.
subroutine readRealCellVectors(readUnit,writeUnit)

   ! Use necessary modules.
   use O_ReadDataSubs

   ! Make sure no funny variables are defined.
   implicit none 

   ! passed parameters.
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   call readData(readUnit,writeUnit,3,3,realVectors(:,:),&
      & len('CELL_VECTORS'),'CELL_VECTORS')

end subroutine readRealCellVectors

! This subroutine will read the definition for a 3D mesh that can be used to
!   evaluate and display the charge density or specific wave function states.
subroutine readNumMeshPoints(readUnit,writeUnit)

   ! Use necessary modules.
   use O_ReadDataSubs

   ! Make sure no funny variables are defined.
   implicit none
   
   ! passed parameters.
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   call readData(readUnit,writeUnit,3,numMeshPoints(:),0,'')

end subroutine readNumMeshPoints


! This subroutne will determine the reciprocal lattice vectors from the
!   real space lattice vectors.  It will also determine the real space crystal
!   volume as a side effect.
subroutine getRecipCellVectors

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3, pi

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the local variables used in this subroutine.
   integer :: i,j ! Loop index variable.  i = col#; j = row#
   ! The following cycle variables are numbers from 1 to dim3.  They are
   !    greater than the loop variable (i or j) by the trailing number, and
   !    they cycle back to 1 when they are greater than the the loop index and
   !    dim3.  (e.g. i = 1; iCycle1 = 2)  (e.g. i = 3; iCycle1 = 1)
   integer :: iCycle1, iCycle2
   integer :: jCycle1, jCycle2

   ! Given the real space primitive lattice vectors (say a b c) each defined in
   !   terms of the orthogonal coordinate vectors (x y z) (a.k.a i k j) we can
   !   find the reciprocal primitive lattice vectors (say a' b' c') by:
   !   a' = 2*Pi * (b x c)/(a . b x c); where . = dot product, x = cross product
   !   b' = 2*Pi * (c x a)/(a . b x c)
   !   c' = 2*Pi * (a x b)/(a . b x c)

   do i = 1, dim3  ! Reciprocal a,b,c
      realCellVolume = 0
      iCycle1 = mod(i,dim3) + 1
      iCycle2 = mod(i+1,dim3) + 1
      do j = 1, dim3  ! x,y,z
         jCycle1 = mod(j,dim3) + 1
         jCycle2 = mod(j+1,dim3) + 1
         recipVectors(j,i) = realVectors(jCycle1,iCycle1) * &
               & realVectors(jCycle2,iCycle2) - realVectors(jCycle2,iCycle1) * &
               & realVectors(jCycle1,iCycle2)
         realCellVolume = realCellVolume + recipVectors(j,i) * realVectors(j,i)
      enddo

      recipVectors(:,i) = 2.0_double * pi * recipVectors(:,i) / realCellVolume
   enddo

   realCellVolume = abs(realCellVolume) ! Volume is always positive
   recipCellVolume = ((2.0_double*pi)**3)/realCellVolume

   ! Also, since we're here, let's compute the inverse of the real lattice
   !   vector matrix.
   invRealVectors(:,:) = recipVectors(:,:) / (2.0_double * pi)

   write (20,100) 'The real  cell volume is: ',realCellVolume
   write (20,100) 'The recip cell volume is: ',recipCellVolume
   write (20,*)   'The recip lattice vectors are: '
   write (20,fmt="(3d18.8)") recipVectors(1:dim3,1)
   write (20,fmt="(3d18.8)") recipVectors(1:dim3,2)
   write (20,fmt="(3d18.8)") recipVectors(1:dim3,3)
   100 format (a,e12.5)

end subroutine getRecipCellVectors


subroutine setCutoffThresh(basisFnConvgTemp,electroConvgTemp)

   ! Include necessary modules.
   use O_Kinds ! Variable precision defined for intrinsic types

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   real (kind=double) :: basisFnConvgTemp
   real (kind=double) :: electroConvgTemp

   ! The meaning of these convergence numbers can be understood as follows:
   ! Consider two atoms at different sites.  Consider an atomic orbital from
   !   each atom.  Each atomic orbital is composed of a set of Gaussian
   !   functions.  Consider the most broad Gaussian from both sets.  When the
   !   overlap between the two Gaussians is considered, it can be represented
   !   by another Gaussian function at some point between the atomic sites.
   !   This is from the Gaussian product theorem.  The idea then is that the
   !   exact integral of this new Gaussian for all orbital types (QN_nlm)
   !   should only be performed *if* the value of the Gaussian at a distance
   !   of *at least* the separation distance between the two atoms is greater
   !   than basisFnConvgTemp.  Call the square of the separation distance
   !   atomSiteSepSqrd.
   ! The value of basisFnConvgTemp is used to compute the Gaussian exponent
   !   value.  For example, if basisFnConvgTemp = 1e-16 then we would write
   !   something like:  1e-16 = exp(-alpha*r^2) = exp(-x).  The value of x is
   !   -ln(1e-16) = 36.48.
   ! The value of alpha in the above expression is known only when the two
   !   particular atoms are known.  The alpha for the product of two Gaussians
   !   is a1 * a2 / (a1 + a2) = a3.  So, if we multiply our value of x by 1/a3
   !   we will have the value of r^2 that will produce 1e-16 = exp(-a3 * r^2).
   ! So, the logBasisFnThresh that is computed below from the basisFnConvgTemp
   !   represents the product of the alpha and r^2 terms in the Gaussian
   !   exponent.  By multiplying by 1/alpha we get a so-called negligability
   !   distance (currentNegligLimit) that is the r^2 part.
   ! If currentNegligLimit > atomSiteSepSqrd then the value of the Gaussian at
   !   a distance of sqrt(atomSiteSepSqrd) is *at least* 1e-16.  On the other
   !   hand, if currentNegligLimit < atomSiteSepSqrd then the value of the
   !   Gaussian at a distance of sqrt(atomSiteSepSqrd) is *less than* 1e-16.
   !   In such a case we say that there is no need to do the integral.
   ! Note that while this does define a hard rule there is not so much physical
   !   intuition as to why this should be the method to discard "too small"
   !   integrals.  The nature of the smallest alpha (a1) from atom #1 and the
   !   smallest alpha (a2) from atom 2 is that the Gaussian produced by their
   !   product will likely be near the midpoint between them.  The
   !   negligability distance is the same for a pair of atoms regardless of the
   !   separation distance so it does make sense to compare this to the
   !   atom separation distance, but the exact number for the cutoff does not
   !   (I think) really give an accurate way to quantify the error in the
   !   calculation.  It feels more ad hoc and "we do it because it works" in
   !   nature.  That's fine but not totally satisfying.
   ! Note also, the electronic convergence follows an identical argument.

   logBasisFnThresh = -log(basisFnConvgTemp)
   logElecThresh   = -log(electroConvgTemp)

end subroutine setCutoffThresh



! This subroutine will initialize two super lattices.  One is real space and
!   the other is in reciprocal space.  Each is constructed from the lattice of
!   their respective spaces that were given in input.  The goal is create a
!   periodic lattice in each space that encompases a sphere whose radius
!   represents the negligability limit.  In theory, all the integrals done in
!   OLCAO suite should be over all space to be absolutely accurate.  Due to
!   finite computing power we make a limit of negligability, beyond which the
!   values of the integerals are not considered.  This limit is controlled by
!   the atomic basis function alphas, and the potential function alphas that
!   govern the decaying exponential terms in their respective descriptions.
subroutine initializeLattice (doRecip)

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3
   use O_AtomicTypes, only: minAtomicAlpha
   use O_PotTypes, only: minPotAlpha
   use O_TimeStamps

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer, intent (in) :: doRecip ! If 1 then compute the recip super lattice.

   ! Define the local variables used in this subroutine.
   integer :: i ! Loop variables

   real (kind=double), dimension(dim3) :: sumSquare
   real (kind=double) :: maxSquareReal

   ! Log the beginning of the lattice initialization.
   call timeStampStart(3)

   ! The determination of the size of the negligability sphere is not an
   !   obvious process.  Some of the code that is here is taken from the
   !   original program without a true understanding of its meaning.  Where
   !   that occures I will try to make a note in the form of a question.

   ! The size of the negligability sphere in real space is determined by the
   !   smallest atomic or potential alpha.  This is the term that dominates
   !   at long distances, so the edge of importance for terms associated with
   !   the smallest alpha will be defined as the negligability boundary.

   ! So, we should be looking for the smallest alpha, and using that to find the
   !   boundary.  However, in the old code we look for the max of 2/atomicAlpha,
   !   and 1/potAlpha + 0.5/atomicAlpha.  Why do we do this?  Lizhi Ouyang says
   !   that the atomicAlphas appear in the integreal in a pair and that is the
   !   reason for the 2 in the numerator, and that the potAlphas only appear in
   !   solitary form so we only need a 1 in the numerator.  I have yet to
   !   confirm this.  Why then do we need the + 0.5/atomicAlpha in the max()
   !   side associated with the potential alphas?
   ! realNegLimit = threshold * max (2/atomicAlpha,1/potAlpha+0.5/atomicAlpha)
   negligLimitReal = logBasisFnThresh * max(2.0_double/minAtomicAlpha, &
         & 1.0_double/minPotAlpha + 0.5_double/minAtomicAlpha)


   ! The boundary is then adjusted by the size of the primitive cell in another
   !   unexplain way.  First we look at the primitive vectors and determine the
   !   largest one from the squares.  x^2 + y^2 + z^2, find the max of a,b,c.
   maxSquareReal = 0.0_double
   do i = 1, dim3
      sumSquare(i) = sum(realVectors(:,i)**2)
      maxSquareReal = max(sumSquare(i),maxSquareReal)
   enddo

   ! limit = limit + 0.75*maxSquareReal+sqrt(3*maxSquareReal*limit)
   negligLimitReal = negligLimitReal + 0.75_double * &
         & maxSquareReal + sqrt(3.0_double * maxSquareReal * negligLimitReal)

   ! What I don't understand is why the limit has to be so large if it is only
   !   controlled by the alpha values in the decaying exponential.  Can we set
   !   the limit just by the alphas?  Why not?  What assumptions about the
   !   structure were made to arrive at the above determination?

   ! Lizhi Ouyang may have addressed this issue in his version of the program.
   !   There he compares the above calculated real negligability limit to the
   !   following expression:  2 * logBasisFnThresh / minAtomicAlpha.  I will do
   !   the same here, but I don't yet know the exact reason for this expression.
   negligLimitReal = max(negligLimitReal, 2.0_double * logBasisFnThresh / &
         & minAtomicAlpha)

   ! Because some electrostatic long range interactions are computed in
   !   reciprocal space we must also make a supercell in recprocal space that
   !   encompases a sphere where the radius is defined again by the minimum
   !   alphas in the decaying exponential.  This section was introduced in
   !   Lizhi Ouyang's version of the program set.
   negligLimitRecip = 4.0 * logElecThresh * minPotAlpha


   ! Make an estimate of the size of each dimension of the superlattice for real
   !   space and for reciprocal space.  The size is given in terms
   !   of the number of primitive cells that need to be replicated in a given
   !   direction to encompass the limit of negligability.  This is done for
   !   both real space, and reciprocal space.
   do i = 1, dim3
      primRepsReal(i)=2 * ceiling(sqrt(negligLimitReal) / &
             & sqrt(sum(realVectors(:,i)**2))) + 1
      primRepsRecip(i)=2 * ceiling(sqrt(negligLimitRecip) / &
             & sqrt(sum(recipVectors(:,i)**2))) + 1
   enddo

   ! Construct the super lattice.  Each combination of lattice direction
   !   repetitions that results in a superlattice where the size is less than
   !   the negligability limit will be counted, and the lattice size information
   !   will be recorded.  It is important that the lattices be arranged in an
   !   array in nondecreasing order by the sum of the squares of the x,y,z
   !   components.  This job is passed off to a subroutine because it must be
   !   done for both the real space and reciprocal space lattices.

   call getNumCells(primRepsReal,numCellsReal,negligLimitReal,realVectors)

   ! Now that we know the array sizes we can allocate space for them.
   allocate(cellSizesReal(numCellsReal))
   allocate(cellDimsReal(dim3,numCellsReal))


   call makeLattice(primRepsReal,numCellsReal,cellSizesReal,cellDimsReal,&
         & negligLimitReal,realVectors)


   if (doRecip .eq. 1) then

      call getNumCells(primRepsRecip,numCellsRecip,negligLimitRecip,&
            & recipVectors)

      ! Now that we know the array sizes we can allocate space for them.
      allocate(cellSizesRecip(numCellsRecip))
      allocate(cellDimsRecip(dim3,numCellsRecip))


      call makeLattice(primRepsRecip,numCellsRecip,cellSizesRecip,&
            & cellDimsRecip,negligLimitRecip,recipVectors)

   endif

   ! Note the number of lattice points in real and reciprocal space.
   write (20,*) 'The number of real  lattice points is:  ',numCellsReal
   write (20,*) 'The number of recip lattice points is:  ',numCellsRecip


   ! Log the date and time we end.
   call timeStampEnd(3)

end subroutine initializeLattice


subroutine initialize3DMesh

   implicit none

   call getRealMagnitudes

   call getRealUnitVectors

   call getRealFractStrideLength

   call getRealAngles

   call getRealFractCrossArea

   call getRealNormalVectors

   call getRealPlaneAngles

end subroutine initialize3DMesh


subroutine getNumCells (primReps,numCells,negligLimit,vectors)

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3, bigThresh

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer, dimension(dim3) :: primReps
   integer :: numCells
   real (kind=double) :: negligLimit
   real (kind=double), dimension(dim3,dim3) :: vectors

   ! Define the local variables used in this subroutine.
   integer :: i,j,k,l,m ! Loop variables.  loop1=i, nestedloop2=j ...
   integer :: isWithinNegligLimit ! Initialized to zero within each call to
                                  !   cellBorderCheck.
   integer, dimension(dim3) :: currentRep
   real (kind=double) :: tempSumOfSquares
   real (kind=double), dimension(dim3) :: tempDimensions
   integer :: smallestOverCounter
   real (kind=double) :: smallestOverNeglig


   ! Initialize the number of cells that will fit in the negligability limit.
   numCells = 0

   ! This is a bit of an annoying way to do this.  We don't know the size of
   !   dimension for the lattice Sizes array and the latticeDims matrix.  So
   !   we go through the whole procedure here twice.  The first time to
   !   determine the size of the array we need, and then second time to fill
   !   the arrays after we allocate space for them.  A better way to do this
   !   would be to make a dynamic static allocation (if that makes sense) that
   !   will certainly hold all the results.  Then fill it is as far as needed,
   !   copy the valid part of the big array to a permanent data structure, and
   !   deallocate the big array.

   ! Initialize a variable that tracks the smallest distance squared
   !   that is greater than the negligability limit.  There is no need to
   !   initialize a counter to track the number of times that this sumOfSquares
   !   distance has been found since this is initialized to bigThresh and will
   !   certainly never be the real smallestOverNeglig.
   smallestOverNeglig  = bigThresh
   smallestOverCounter = 0

   ! Start three nested loops for the three dimensions.
   do i = 1, primReps(1)
      currentRep(1) = i - (1 + primReps(1))/2
      do j = 1, primReps(2)
         currentRep(2) = j - (1 + primReps(2))/2
         do k = 1, primReps(3)
            currentRep(3) = k - (1 + primReps(3))/2

            ! Calculate the vector to the current cell.
            tempDimensions(:)=0.0_double
            do l = 1, dim3  ! Picks a, b, c vector
               do m = 1, dim3  ! Picks x, y, z vector
                  tempDimensions(m) = tempDimensions(m) + &
                                    & vectors(m,l) * currentRep(l)
               enddo
            enddo

            ! Calculate the distance squared to that lattice point.
            tempSumOfSquares = sum(tempDimensions(:)**2)

            ! Check to see if that point is outside the negligability limit.
            if (tempSumOfSquares > negligLimit) then

               ! If the lattice point is outside the negligability limit, then
               !   we next must check to see if any part of that lattice
               !   point's cell is possibly inside the negligability limit.

               ! If, at any time, some part of this cell is within the
               !   negligability limit then we must include it in the count.
               !   Otherwise we cycle without counting this cell since it is
               !   entirly outside the negligability sphere.

               ! This is accomplished by doing a similar task as the above
               !   calculation of the vector to the seven untested vertices
               !   of the current cell.
               call cellBorderCheck (vectors(:,:), tempDimensions(:), &
                     & negligLimit, isWithinNegligLimit)

               ! If it is determined that the cell is not within the
               !   negligability limit then we cycle to avoid including
               !   this cell in the count.
               if (isWithinNegligLimit == 0) then
                  cycle
               endif
            endif

            ! Increment the array size counter by 1 since this lattice
            !   point is within the negligability limit.
            numCells = numCells + 1
         enddo
      enddo
   enddo
end subroutine getNumCells


subroutine makeLattice (primReps,numCells,cellSizes,cellDims,negligLimit,&
      & vectors)

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3
   use O_SortSubs ! Subroutines for mergesort sorting.

   ! Define the dummy variables passed to this subroutine.
   integer, dimension(dim3) :: primReps
   integer :: numCells
   real (kind=double), dimension (:) :: cellSizes
   real (kind=double), dimension (:,:) :: cellDims
   real (kind=double) :: negligLimit
   real (kind=double), dimension(dim3,dim3) :: vectors

   ! Define the local variables used in this subroutine.
   integer :: i,j,k,l,m ! Loop variables.  loop1=i, nestedloop2=j ...
   integer :: isWithinNegligLimit ! Initialized to zero within each call to
                                  !   cellBorderCheck.
   integer, dimension(dim3) :: currentRep
   real (kind=double) :: tempSumOfSquares
   real (kind=double), dimension(dim3) :: tempDimensions

   ! Define the local variables used for sorting.
   real (kind=double), allocatable, dimension (:)   :: cellSizesTemp
   real (kind=double), allocatable, dimension (:,:) :: cellDimsTemp
   integer,            allocatable, dimension (:)   :: sortOrder
   integer,            allocatable, dimension (:)   :: segmentBorders

   ! Allocate space to temporarily hold the directly computed data.
   allocate (cellSizesTemp(numCells))
   allocate (cellDimsTemp(dim3,numCells))
   allocate (sortOrder(numCells))
   allocate (segmentBorders(numCells+1))

   ! Now, we repeat the procedure except that we record the values this time.
   !   Another alternative way to do this would be to allocate a rediciously
   !   HUGE static array based on the given cell size, the negligability limit,
   !   and the alphas provided.  Then, fill out that array, determining the
   !   actual precise size at the same time.  Then copy the results from the
   !   HUGE array into the exact sized array and deallocate the HUGE one.  That
   !   would probably be faster.
   numCells = 0
   do i = 1, primReps(1)
      currentRep(1) = i - (1 + primReps(1))/2
      do j = 1, primReps(2)
         currentRep(2) = j - (1 + primReps(2))/2
         do k = 1, primReps(3)
            currentRep(3) = k - (1 + primReps(3))/2
            tempDimensions=0.0_double
            do l = 1, dim3
               do m = 1, dim3
                  tempDimensions(m) = tempDimensions(m) + &
                                    & vectors(m,l) * currentRep(l)
               enddo
            enddo
            tempSumOfSquares = sum(tempDimensions(1:dim3)**2)
            if (tempSumOfSquares > negligLimit) then

               ! If the lattice point is outside the negligability limit, then
               !   we next must check to see if any part of that lattice
               !   point's cell is possibly inside the negligability limit.

               ! If, at any time, some part of this cell is within the
               !   negligability limit then we must include it in the count.
               !   Otherwise we cycle without counting this cell since it is
               !   entirly outside the negligability sphere.

               ! This is accomplished by doing a similar task as the above
               !   calculation of the vector to the seven untested vertices
               !   of the current cell.
               call cellBorderCheck (vectors(:,:), tempDimensions(:), &
                     & negligLimit, isWithinNegligLimit)

               ! If it is determined that the cell is not within the
               !   negligability limit then we cycle to avoid including
               !   this cell in the count.
               if (isWithinNegligLimit == 0) then
                  cycle
               endif
            endif
            numCells = numCells + 1
            cellSizesTemp(numCells) = tempSumOfSquares
            cellDimsTemp(:dim3,numCells) = tempDimensions
         enddo
      enddo
   enddo

   ! Sort the results in order of increasing sum of squares using a merge sort.

   ! Initialize the segment borders so that every number is a border.
   segmentBorders(1) = 0
   do i = 1, numCells
      segmentBorders(i+1) = i
   enddo

   call mergeSort (cellSizesTemp,cellSizes,sortOrder,segmentBorders,numCells)

   ! Copy the results to the cellDims array
   do i = 1, numCells
      cellDims(:dim3,i) = cellDimsTemp(:dim3,sortOrder(i))
   enddo

   ! Deallocate the unnecessary arrays now.
   deallocate (cellSizesTemp)
   deallocate (cellDimsTemp)
   deallocate (sortOrder)
   deallocate (segmentBorders)

end subroutine makeLattice


subroutine cellBorderCheck (vectors, tempDimensions, negligLimit, &
      & isWithinNegligLimit)

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension(dim3,dim3) :: vectors
   real (kind=double), dimension(dim3) :: tempDimensions
   real (kind=double) :: negligLimit
   integer :: isWithinNegligLimit

   ! Define the local variables used in this subroutine.
   integer :: i ! Loop variables
   real (kind=double), dimension (dim3) :: cellShifts
   real (kind=double), dimension (dim3) :: currentTest

   ! Initialize the isWithinNegligLimit return value to zero to assume that
   !   this cell is not within the negligability limit.
   isWithinNegligLimit = 0

   ! Initialize the distances that must be shifted in each direction.
   cellShifts(:) = 0.0_double
   do i = 1, dim3
      cellShifts(:) = cellShifts(:) + vectors(:,i)
   enddo

   ! Begin testing the 7 possible locations that would indicate whether
   !   this cell has some part within the negligability distance.  This is
   !   done by following a path from one point to the next via a sequence of
   !   cellShifts.

   currentTest(:) = tempDimensions(:)
   currentTest(1) = currentTest(1) + cellShifts(1)
   if (sum(currentTest(:)**2) < negligLimit) then
      isWithinNegligLimit = 1
      return
   endif

   currentTest(2) = currentTest(2) + cellShifts(2)
   if (sum(currentTest(:)**2) < negligLimit) then
      isWithinNegligLimit = 1
      return
   endif

   currentTest(1) = currentTest(1) - cellShifts(1)
   if (sum(currentTest(:)**2) < negligLimit) then
      isWithinNegligLimit = 1
      return
   endif

   currentTest(3) = currentTest(3) + cellShifts(3)
   if (sum(currentTest(:)**2) < negligLimit) then
      isWithinNegligLimit = 1
      return
   endif

   currentTest(1) = currentTest(1) + cellShifts(1)
   if (sum(currentTest(:)**2) < negligLimit) then
      isWithinNegligLimit = 1
      return
   endif

   currentTest(2) = currentTest(2) - cellShifts(2)
   if (sum(currentTest(:)**2) < negligLimit) then
      isWithinNegligLimit = 1
      return
   endif

   currentTest(1) = currentTest(1) - cellShifts(1)
   if (sum(currentTest(:)**2) < negligLimit) then
      isWithinNegligLimit = 1
      return
   endif

end subroutine cellBorderCheck


subroutine findLatticeVector (arbitraryVector, latticeVector)

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3, bigThresh

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension (:), intent(in)  :: arbitraryVector
   real (kind=double), dimension (:), intent(out) :: latticeVector

   ! Define the local variables used in this subroutine.
   integer :: i ! Loop variable.
   integer :: currentWignerCell
   integer :: minWignerCell
   integer :: octantNumber
   real (kind=double) :: exactAxisCell
   real (kind=double) :: roundedAxisCell
   real (kind=double) :: minDistance
   real (kind=double) :: currentDistance
   real (kind=double), dimension (dim3) :: remainderVector

   ! Initialize the minimum wigner cell to a large number.
   minWignerCell = 100000000

   ! Initialize the return lattice vector to 0.
   latticeVector(:) = 0.0_double

   ! Initialize the octant number tracker.  This will keep track of which
   !   octant the resultant vector goes through.  For an understanding of
   !   the term "octant" check "initializeFindVec" below.
   octantNumber = 1
   do i = 1, dim3

      ! Determine the number of cells in a specific axial direction that
      !   are of a distance less than that of the incoming arbitrary vector.
      exactAxisCell = sum(arbitraryVector(:)*invRealVectors(:,i))
      roundedAxisCell = real(int(abs(exactAxisCell) + 0.5_double),double)
      if (exactAxisCell < 0.0_double) then
         roundedAxisCell = -roundedAxisCell
      endif

      ! Increment the lattice vector to that point
      latticeVector(:) = latticeVector(:) + realVectors(:,i) * roundedAxisCell

      ! Identify the current octant of the lattice vector
      if (exactAxisCell - roundedAxisCell < 0.0_double) then
         octantNumber = octantNumber + 2**(i - 1)
      endif
   enddo

   ! Determine the number of Wigner-Seitz cells that are used to enclose the
   !   current octant of the lattice vector.  If there are more than one
   !   Wigner-Seitz cells for this octant, then the lattice vector must be
   !   checked for each of them to see which is the closest.
   if (numWignerCellsNeeded(octantNumber) > 1) then

      ! Define the current remainder vector.
      remainderVector = arbitraryVector - latticeVector

      ! Initialize the min distance
      minDistance = bigThresh

      do i = 1, numWignerCellsNeeded(octantNumber)

         ! Get the cell number of the current Wigner-Seitz cell
         currentWignerCell = overlappingWignerCells(i,octantNumber)

         ! Get the distance to that lattice point
         currentDistance = sum((remainderVector(:) - &
                          & cellDimsReal(:,currentWignerCell))**2)

         ! Determine if this distance is shorter than the current min
         !   distance, and update the necessary variables if so.
         if (currentDistance < minDistance) then
            minDistance   = currentDistance
            minWignerCell = currentWignerCell
         endif
      enddo

      ! If the minimum distance lattice point is not 1 then we must
      !   adjust the value of the nearest lattice point.
      if (minWignerCell > 1) then
         latticeVector = latticeVector + cellDimsReal(:,minWignerCell)
      endif
   endif
end subroutine findLatticeVector



! The basic purpose of this subroutine is to create a convenient system
!   to rapidly search for the superlattice lattice vector that is nearest
!   to an arbitrary vector.  The procedure to do so can be summed up as
!   follows:
! 1)  Create a Wigner-Seitz cell.  (Cell in real space of minimal volume that
!     can be replicated to fill all space.)
! 2)  Use the system's given cell parameters to define eight cells around the
!     origin.  These cells are called octants and each is a 1/8 sized
!     replication of the system's given cell that has been translated.  Imagine
!     the system's given cell centered on the origin and divided into 8 smaller
!     (but same shaped) units.
! 3)  Compare a Wigner-Seitz cell centered at every relavent (close to the
!     origin) lattice point to each octant mini cell.  When the Wigner-Seitz
!     cell and the octant have a shared volume then the lattice point of the
!     Wigner-Seitz cell is recorded as relavent for the later search for a
!     lattice vector that is nearest to an arbitrary vector.
subroutine initializeFindVec

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3, pi
   use O_TimeStamps

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define local variables
   real (kind=double), allocatable, dimension (:) :: rootCellSizes ! The sqrt
         ! of the values in cellSizes * 0.5.  This is used to help find
         ! the values of the wigner-seitz cell (real space smallest volume
         ! cell).
   real (kind=double), allocatable, dimension (:,:) :: cellUnitDims ! The values
         ! in cellDims divided by the corresponding rootCellSize value
         ! before the 0.5 multiplication.  These are basically unit vectors
         ! that point in the direction of the origin of a cell in the
         ! superlattice.
   real (kind=double), allocatable, dimension (:) :: rootWignerSizes ! Copies
         ! of the values from rootCellSizes that contribute to the Wigner-
         ! Seitz cell.
   real (kind=double), allocatable, dimension (:,:) :: wignerUnitDims ! Copies
         ! of the values from cellUnitDims that contribute to the Wigner-
         ! Seitz cell.
   real (kind=double), allocatable, dimension (:) :: rootOctantSizes ! The
         ! distance to the plane that forms one border of an octant of the
         ! given cell.
   real (kind=double), allocatable, dimension (:,:) :: octantUnitDims ! The
         ! direction to the plane that forms one border of an octant of the
         ! given cell.
   real (kind=double), allocatable, dimension (:,:) :: unitDotMatrix ! This is
         ! a dot product matrix for every unit vector of cellUnitDims with
         ! every other cellUnitDims vector.
   integer, allocatable, dimension (:) :: cellSegmentCounter ! This array will
         ! be used to count the number of times a given cellDims plane
         ! contributes to a segment that is contained within all other
         ! cellDims planes.
   integer :: maxWignerCell ! This is simply an array index of the cell in
         ! probe%realSuperLattice%cellSizes which has the greatest distance
         ! from the origin.  It is used to limit the search for volume sharing
         ! Wigner-Seitz cells.
   integer :: numWignerSegments ! This is a counter that is incremented in a
         ! subroutine for each border segment that has a score of 3 or more in
         ! the search for borders of the minimal volume cell (Wigner-Seitz).

   ! Start setting up structures to find lattice vectors.
   call timeStampStart(4)

   ! Allocate space to hold the square roots of the cellSizes array, and the
   !   values for each cell's unit vectors.  Also allocate space for the
   !   counter tracking each cell's contribution to the Wigner-Seitz cell
   !   NOTE that we do not store the values for the cell at 0,0,0 in any case.
   allocate(rootCellSizes (numCellsReal - 1))
   allocate(cellUnitDims  (dim3,numCellsReal - 1))

   ! Initialize the variable that contains the index number for the cell that
   !   is a maximum distance from the origin.
   maxWignerCell = 0

   ! Create a matrix that is the inverse of the real space lattice vectors
   !   matrix and is related to the reciprocal space lattice vectors by a
   !   factor of 2*pi.
   invRealVectors(:,:) = recipVectors(:,:) / (2.0_double * pi)

   ! Determine the unit vectors for each of the SuperLattice cell vectors
   !   except the cell at the center with the dim values of 0,0,0 and size
   !   value of one cell.
   call createUnitVectors(rootCellSizes,cellUnitDims)

   ! Now we determine the wigner-seitz cell based on the superlattice
   !   that was constructed earlier.  The wigner-seitz cell is
   !   the smallest cell in real space that can be replicated to form the
   !   complete superlattice.  It is akin to the first bruillion zone in
   !   reciprocal space except that this is in real space.  The procedure
   !   to the best of my understanding is as follows:

   ! 1)  Create a dot product matrix of the unit vectors for each cell.  This
   !     is used as a reference later when calculating certain other vectors.
   ! 2)  Consider each pair of cellDims vectors.  Each vector indicates the
   !     normal of a plane.  Then consider the line formed by the intersection
   !     of the pair of planes.  (Clearly, the pairs of cellDims vectors that
   !     produce parallel planes will be skipped.)
   ! 3)  Now we consider all the other cellDims vectors in turn.  Each of
   !     these cellDims vectors produces a third plane.  This plane may or
   !     may not intersect the line segment formed in (2).
   !   i)   If it does intersect the segment then the bounds of the segment
   !        are reduced.
   !      a)  If the segment bounds are reduced such that the segment has a
   !          negative length, then this segment is no longer considered for
   !          inclusion in the Wigner-Seitz cell, and the next pair of
   !          cellDims vectors is considered from (2).
   !   ii)  If it does not intersect the segment then we consider if the
   !        segment is "inside" or "outside" the plane.  "Inside" means that
   !        the segment is closer to the system origin than 1/2 the distance
   !        to the third plane.  "Outside" means the opposite.  The 1/2 is
   !        used in the determination of the Wigner-Seitz cell.
   !      a)   If the segment is "inside" then the segment can possibly be a
   !           component in the Wigner-Seitz cell.  Therefor, the comparison
   !           of this segment continues.
   !      b)   If the segment is "outside" then the segment can never be a
   !           component in the Wigner-Seitz cell.  Therefor all consideration
   !           of this pair from (2) is stopped and the next pair is
   !           considered.
   ! 4)  If a segment completes all the comparisons of (3) without being
   !     rejected then both of the cellDims vectors that contributed to the
   !     formation of that segment have their score incremented by one
   !     (scores all start at 0).
   ! 5)  cellDims vectors with a score of three or more are included in the
   !     Wigner-Seitz cell.  At the moment I have not puzzled out exactly why
   !     it must be three.  I am sure that there is a simple geometric
   !     argument.  If someone in the future wants to detail it here that
   !     would be most welcome.

   ! NOTE that the above explaination is the best that I have been able to
   !   understand from the original (poorly documented) code.  My docs may
   !   have some conceptual errors but I have not found logic errors in the
   !   old code so I will replicate its functionality exactly.  Therefor if
   !   there is some important concept that I am not understanding from the
   !   old code it will still be included even though my documentation does
   !   not mention it.

   ! Create a dot product matrix of all the unit vectors for each cell.
   !   This matrix is used to find a vector to a point on the line created
   !   when two planes intersect.  This step corresponds to step 1 above.

   ! Allocate space for the resultant matrix first.  This has to be done
   !   outside the procedure call since the matrix has an allocatable attribute.
   allocate (unitDotMatrix(numCellsReal - 1,numCellsReal - 1))

   ! Then compute the matrix.
   call dotProductMatrix(cellUnitDims,cellUnitDims,unitDotMatrix)

   ! Compare the line formed by the intersection of the planes defined by
   !   any given pair of cellDims vectors (normals) to the planes formed by
   !   every other cellDims vector.  If a segment of the line exists within
   !   every other cellDims vector, then that pair of cellDims vectors will
   !   have their score incremented by one.  The scores are used later.  This
   !   step corresponds to steps 2, 3, and 4 above.
   !   NOTE that the integer (0) at the end is used to indicate that
   !   identical vectors should not be searched for in the cellUnitDims
   !   list of vectors.  In this execution of the intersection test there
   !   will be no identical vectors.

   ! First allocate space to store the counter for each cell and initialize it.
   allocate(cellSegmentCounter(numCellsReal - 1))
   cellSegmentCounter(:) = 0

   ! Then perform the intersection test.
   call intersectionTest(cellUnitDims,unitDotMatrix,rootCellSizes,&
                       & cellSegmentCounter,numCellsReal - 1,0)

   ! Check the results stored in the cellSegmentCounter.  Any cell that has a
   !   score of 3 or more will be included in the Wigner-Seitz cell.
   call minimalBorderCheck(cellSegmentCounter,numWignerSegments,maxWignerCell)

   ! Allocate sufficient space to hold the number of Wigner-Seitz cell borders.
   allocate(rootWignerSizes(numWignerSegments))
   allocate(wignerUnitDims(dim3,numWignerSegments))

   ! Record the unit vectors and origin distances for the cells that contribute
   !   to the borders for the Wigner-Seitz minimal border cell.
   call minimalBorderRecord(cellSegmentCounter,cellUnitDims,rootCellSizes,&
                          & wignerUnitDims,rootWignerSizes)

   ! Now that the Wigner-Seitz cell has been created we will define an
   !   octant of the given cell.  This is the given cell centerd on the
   !   origin and divided into eight parts.  Only one octant of that cell
   !   needs to be defined here since it can be easily translated to
   !   represent any octant of the given cell.  In the old program this was
   !   termed the "prismatic cell", but I was unable to find a consistant
   !   precise deminition of "prismatic cell" that matched what this process
   !   does so I think that term should not be used here.

   ! First we have to allocate space to hold the values for the octant.
   allocate (rootOctantSizes(6)) ! Always 6 because only 6 sides are needed to
                                 !   define any bravis lattice.
   allocate (octantUnitDims(dim3,6)) ! Six for the same reason as above.

   ! Then we can define the octant.
   call defineOctant(invRealVectors,octantUnitDims,rootOctantSizes,realVectors)

   ! Now the Wigner-Seitz cell and the octants will be compared.  Each octant
   !   will be considered in turn.  For each octant it will be determined
   !   how many, and at what lattice point a set of Wigner-Seitz cells will
   !   have to be positioned to completely enclose that octant.  The number
   !   of Wigner-Seitz cells necessary will be stored in:
   !   numWignerCellsNeeded.  The lattice point index of the Wigner-
   !   Seitz cells that are necessary to cover a lattice point will be stored
   !   in:  overlappingWignerCells.
   call compareOctantsAndWignerCells(wignerUnitDims,rootWignerSizes, &
         & octantUnitDims,rootOctantSizes,realVectors,maxWignerCell)

   ! That's all folks.  The number of Wigner-Seitz cells needed and their
   !   locations (lattice points) are all that was needed.  All we have to do
   !   now is deallocate the stuff that won't be used ever again.
   deallocate (rootCellSizes)
   deallocate (cellUnitDims)
   deallocate (unitDotMatrix)
   deallocate (cellSegmentCounter)
   deallocate (rootWignerSizes)
   deallocate (wignerUnitDims)
   deallocate (rootOctantSizes)
   deallocate (octantUnitDims)

   ! Log the date and time we finish.
   call timeStampEnd(4)

end subroutine initializeFindVec

subroutine createUnitVectors (rootCellSizes,cellUnitDims)

   ! Include the modules we need.
   use O_Kinds
   use O_Constants, only: dim3


   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension (:)   :: rootCellSizes
   real (kind=double), dimension (:,:) :: cellUnitDims
   

   ! Define the local variables used in this subroutine.
   integer :: i ! Loop index variable.
   
   ! Assign the values for all but the 0,0,0 cell.
   do i = 1, numCellsReal - 1

      ! Get the square root of the current lattice size.  Note that the first
      !   lattice point at 0,0,0 is skipped.
      rootCellSizes(i) = sqrt(cellSizesReal(i+1))

      ! Define a unit vector of the current lattice point's origin vector
      cellUnitDims(:dim3,i) = cellDimsReal(:dim3,i+1) / rootCellSizes(i)

      ! Multiply by 0.5 to form boundary distance between two
      !   adjacent minimal sized cells.  (wigner-seitz cells).
      rootCellSizes(i) = 0.5 * rootCellSizes(i)
   enddo
end subroutine createUnitVectors

subroutine dotProductMatrix (vectorArray1,vectorArray2,matrix)

   ! Include the necessary modules.
   use O_Kinds
   use O_Constants, only: dim3

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension (:,:) :: vectorArray1
   real (kind=double), dimension (:,:) :: vectorArray2
   real (kind=double), dimension (:,:) :: matrix

   ! Define the local variables used in this subroutine.
   integer :: i,j ! Loop variables.  loop1=i, nestedloop2=j ...
   integer :: n   ! Size of the incoming array.

   n = size (vectorArray1,2)

   do i = 1,n
      do j = i,n
         matrix(i,j) = dot_product(vectorArray1(:dim3,i),vectorArray2(:dim3,j))
         matrix(j,i) = matrix(i,j) ! This may be unnecessary
      enddo
   enddo
end subroutine dotProductMatrix

subroutine intersectionTest (unitDims, unitDotMatrix, rootSizes, cellCounter,&
                            & numVectors, identTrap)

   ! Include the necessary modules.
   use O_Kinds
   use O_Constants, only: smallThresh

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension (:,:) :: unitDims
   real (kind=double), dimension (:,:) :: unitDotMatrix
   real (kind=double), dimension (numVectors)   :: rootSizes
   integer, dimension (numVectors) :: cellCounter
   integer :: numVectors
   integer :: identTrap

   ! Define the local variables used in this subroutine.
   integer :: i,j,k ! Loop variables.  loop1=i, nestedloop2=j ...
   integer, allocatable, dimension (:) :: trapVector
   integer :: failure ! Simple flag to detect failure conditions
   real (kind=double), dimension (3) :: unitSegment
   real (kind=double), dimension (3) :: toSegmentVector
   real (kind=double) :: factor1, factor2 ! Simple multiplication factors that
                                          !   are saved to reduce computation.
   real (kind=double) :: upperBound, lowerBound ! Bounds on the length of the
                                                !   segment.
   real (kind=double) :: orientation, inclusion ! Calculated numbers that
                                                !   affect the choice of wigner
                                                !   seitz cell and the lattice
                                                !   points for searching.

   ! Initialize the counter that tracks the number of times that cell vector
   !   contributes to a segment that is contained within all the other
   !   segments.
   cellCounter(:) = 0

   ! Initialize the trapVector array to all zeros.  If trapping of identical
   !   vectors was not requested, then this will allow vectors to be compared.
   allocate (trapVector(numVectors))
   trapVector(:) = 0

   ! If trapping of identical vectors was requested then this will mark the
   !   second vector in every pair of identical vectors.  This will allow the
   !   first to be considered, but the second will be skipped.
   if (identTrap == 1) then
      do i = 1,numVectors-1

         ! Skip any vectors that have already been marked.
         if (trapVector(i) == 1) cycle

         do j = i+1,numVectors

            ! Are vector 1 and vector 2 at the same orientation (parallel)?
            !   If NOT, then cycle.  (i.e. only proceed if parallel)
            if (abs(abs(unitDotMatrix(i,j))-1.0) > smallThresh) cycle

            ! Are vector 1 and vector 2 at the same distance?  If NOT, then
            !   cycle.  (i.e. only proceed if at the same distance).
            if (abs(rootSizes(i) - unitDotMatrix(i,j) * rootSizes(j)) > &
               &smallThresh) cycle

            ! At this point vectors 1 and 2 are sufficiently parallel and
            !   at sufficiently similar distances to be regarded as the
            !   same vectors.  Now, we only mark vector #2.  In this way,
            !   when the search is done below, vector 1 will be considered
            !   as a border that could help define a shared volume between
            !   the current Wigner-Seitz cell and the current octet, but its
            !   identical twin (vector 2) will not be considered.  This
            !   prevents double counting of a single border.
            trapVector(j) = 1
         enddo
      enddo
   endif

   do i = 1,numVectors-1

      ! Skip any vectors that have been marked for being the second part of
      !  an identical pair.
      if (trapVector(i) == 1) cycle

      do j = i+1,numVectors

         ! Skip any vectors that have been marked for being the second part
         !  of an identical pair.
         if (trapVector(j) == 1) cycle

         ! If the current two planes are parallel then they must be skipped
         !   since they do not form a line segment via intersection.
         if (abs(abs(unitDotMatrix(i,j))-1.0) < smallThresh) cycle

         ! Identify a vector that points to a point on the line.
         factor1 = (rootSizes(i) - rootSizes(j) * unitDotMatrix(i,j))/&
                   (1.0 - unitDotMatrix(i,j)**2)
         factor2 = (rootSizes(j) - rootSizes(i) * unitDotMatrix(i,j))/&
                   (1.0 - unitDotMatrix(i,j)**2)
         toSegmentVector(:) = factor1 * unitDims(:,i) + factor2 * unitDims(:,j)

         ! Identify a unit vector that points in the direction of the segment.
         !   This is simply a cross product.
         unitSegment(1) = unitDims(2,i) * unitDims(3,j) - &
                         &unitDims(3,i) * unitDims(2,j)
         unitSegment(2) = unitDims(3,i) * unitDims(1,j) - &
                         &unitDims(1,i) * unitDims(3,j)
         unitSegment(3) = unitDims(1,i) * unitDims(2,j) - &
                         &unitDims(2,i) * unitDims(1,j)

         ! Initialize the bounds for the length of the segment to large numbers.
         upperBound =  10000000.0_double
         lowerBound = -10000000.0_double

         ! Initialize the failure flag to success.
         failure = 0

         ! Compare this segment to the planes formed by every other vector in
         !   unitDims and rootSizes.
         do k = 1, numVectors

            ! Do not include the currently selected vectors in the comparison.
            if ((k == i) .or. (k == j)) cycle

            ! Calculate a value that relates the degree to which the k loop
            !   vector is perpendicular (0) or parallel (1) to the unitSegment
            !   vector.
            orientation = sum(unitSegment(:) * unitDims(:,k))

            ! Calculate a value that applies to cases where the k loop vector
            !   is sufficiently perpendicular to the unitSegment vector.  We
            !   have to check that if the segment is very far away from the
            !   plane formed by the k loop vector then does the k loop vector
            !   include or exclude the segment?
            inclusion = rootSizes(k) - sum(toSegmentVector(:) * unitDims(:,k))

            ! This part is tricky.  First determine if the orientation of the
            !   unitSegment vector is parallel (close to 1) or perpendicular
            !   (close to 0) with the current k loop vector.
            if (abs(orientation) < smallThresh) then

               ! Here they are oriented close to perpendicular so we must find
               !   if it is possible for the unitSegment to be a part of the
               !   minimal volume.

               !   This would happen if the k loop vector (even though it is
               !   close to perpendicular) had a large magnitude.  If it could
               !   be included, then we just cycle to the next k loop vector.
               !   Otherwise it cannot be included, and we must choose a new
               !   unitSegment since this choice of unitSegment has failed.
               if (inclusion > -smallThresh) then
                  cycle
               else
                  ! The unitSegment choice is outside the plane formed by the
                  !   current k loop vector.  (The k loop vector is a normal)
                  failure = 1
                  exit
               endif
            else

               ! Here the unitSegment and the k loop plane are sufficiently
               !   non-perpendicular so that we can calculate where the plane
               !   will intersect with the line segment.  This will cut the
               !   upper or lower bounds for the end points of the unitSegment.
               if (orientation < 0.0 .and.&
                  &inclusion/orientation > lowerBound) then
                  lowerBound = inclusion/orientation
               endif
               if (orientation > 0.0 .and.&
                  &inclusion/orientation < upperBound) then
                  upperBound = inclusion/orientation
               endif

               ! Did we cut it so that the segment still has a positive length?
               if (upperBound > lowerBound + smallThresh) then
                  ! Try the next k loop vector to cut it even smaller.
                  cycle
               else
                  ! It was cut to a negative length and cannot exist as a
                  !   border for the minimal volume.
                  failure = 1
                  exit
               endif
            endif
         enddo

         ! The segment still has a positive length and has not been dis-
         !   included by any perpendicular vectors so we add a score point to
         !   each of the vectors that contributed to the creation of the
         !   segment.
         if (failure == 0) then
            cellCounter(i) = cellCounter(i) + 1
            cellCounter(j) = cellCounter(j) + 1
         endif
      enddo
   enddo

   ! Deallocate anything that was allocated and will not be used later.
   deallocate (trapVector)

end subroutine intersectionTest

subroutine minimalBorderCheck (cellCounter,numWignerSegments,maxWignerCell)

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer, dimension (:) :: cellCounter
   integer :: numWignerSegments
   integer :: maxWignerCell

   ! Define the local variables used in this subroutine.
   integer :: i  ! i = loop index
   integer :: numCells

   ! Establish the number of possible segments defined by the cells
   numCells = size (cellCounter)

   ! Determine how many of those cells contribute to the Wigner-Seitz cell so
   !   that we can allocate the wigner data structures
   numWignerSegments = 0
   do i = 1,numCells
      if (cellCounter(i) >= 3) then
         numWignerSegments = numWignerSegments + 1

         ! Determine the maximun index (and hence maximum distance) Wigner-
         !   Seitz cell.
         if (i > maxWignerCell) then
            maxWignerCell = i
         endif
      endif
   enddo

end subroutine minimalBorderCheck

subroutine minimalBorderRecord (cellCounter,cellUnitDims,rootCellSizes,&
                               &wignerUnitDims,rootWignerSizes)

   ! Include the modules we need
   use O_Kinds

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer, dimension (:) :: cellCounter
   real (kind=double), dimension (:,:) :: cellUnitDims
   real (kind=double), dimension (:)   :: rootCellSizes
   real (kind=double), dimension (:,:) :: wignerUnitDims
   real (kind=double), dimension (:)   :: rootWignerSizes

   ! Define the local variables used in this subroutine.
   integer :: i  ! i = loop index
   integer :: numCells
   integer :: numWignerSegments

   ! Establish the number of possible segments defined by the cells
   numCells = size (cellCounter)

   ! Initialize the counter for the number of segments that contribute
   !   to the Wigner-Seitz cell.
   numWignerSegments = 0

   ! Check each index in cellCounter(:).  Any cell with a score of 3 or more
   !   is to be included as a contributer to the Wigner-Seitz cell.  If
   !   included, then record rootCellSize (distance to cell) and cellUnitDims
   !   (direction to cell) in the corresponding Wigner-Seitz variables.
   do i = 1,numCells
      if (cellCounter(i) >= 3) then
         numWignerSegments = numWignerSegments + 1
         rootWignerSizes(numWignerSegments) = rootCellSizes(i)
         wignerUnitDims(:dim3,numWignerSegments) = cellUnitDims(:dim3,i)
      endif
   enddo

end subroutine minimalBorderRecord

subroutine sharedVolumeCheck (cellCounter,latticePointIndex,octantIndex,&
                             &numWignerCellsNeeded,overlappingWignerCells)

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer, dimension (:) :: cellCounter
   integer :: latticePointIndex
   integer :: octantIndex
   integer, dimension (8) :: numWignerCellsNeeded
   integer, dimension (:,:) :: overlappingWignerCells

   ! Define the local variables used in this subroutine.
   integer :: numPlanes

   ! Initialize the plane counter
   numPlanes = 0

   ! Count the number of planes that have a score of 3 or more
   numPlanes = numPlanes + count(cellCounter(:)>=3)

   ! If the number of planes with a score of 3 or more is at least 4, then
   !   this indicates a shared volume between the current octant and the
   !   current Wigner-Seitz cell.
   if (numPlanes >= 4) then

      ! Increment the number of Wigner-Seitz cells that are needed to
      !   completely enclose the current octant.
      numWignerCellsNeeded(octantIndex) = numWignerCellsNeeded(octantIndex) + 1

      ! Record the lattice point of the current Wigner-Seitz cell for use later
      !   in determining the closest lattice point to any arbitrary vector.
      overlappingWignerCells(numWignerCellsNeeded(octantIndex),octantIndex) = &
         & latticePointIndex
   endif

end subroutine sharedVolumeCheck

subroutine defineOctant (invRealVectors,octantUnitDims,rootOctantSizes,&
      & realVectors)

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension (dim3,dim3) :: invRealVectors
   real (kind=double), dimension (:,:)       :: octantUnitDims
   real (kind=double), dimension (:)         :: rootOctantSizes
   real (kind=double), dimension (dim3,dim3) :: realVectors

   ! Define the local variables used in this subroutine.
   integer :: j ! j = loop counter
   integer :: currentIndex ! There are six sides to the octant but we only have
                           !   to do three operations to determine all six
                           !   sides.  So we have this counter to allow simple
                           !   back references to indices.
   real (kind=double) :: unitDotProduct ! This is a dot product over each
                           ! column of the real space cell vectors
                           !   (realVectors), and the inverse of the real space
                           !   cell vectors (invRealVectors).
   real (kind=double) :: unitDistance ! The square root of the sum of squares
                           !   for each column in the invRealVectors.

   do j = 1,dim3

      ! Identify the current index
      currentIndex = 2*j

      ! Obtain the dot product and distance
      unitDotProduct = sum(realVectors(:,j) * invRealVectors (:,j))
      unitDistance = sqrt(sum(invRealVectors(:,j)**2))

      ! Assign unit vector directions.  These vectors point in the direction
      !   of normals to planes that form the borders of the octant being
      !   defined.
      octantUnitDims(:,currentIndex)   = invRealVectors(:,j)/unitDistance
      octantUnitDims(:,currentIndex-1) = -octantUnitDims(:,currentIndex)

      ! Assign the distances.  Remember that these vectors are normals to
      !   planes that define the octant.
      rootOctantSizes(currentIndex)   = 0.5_double*unitDotProduct/unitDistance
      rootOctantSizes(currentIndex-1) = 0.0_double
   enddo

end subroutine defineOctant

subroutine combineWignerAndOctants(wignerUnitDims,rootWignerSizes,&
      & octantUnitDims,rootOctantSizes,wignerOctantUnitDims,&
      & rootWignerOctantSizes)

   ! Include the modules we need
   use O_Kinds

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension (:,:) :: wignerUnitDims
   real (kind=double), dimension (:)   :: rootWignerSizes
   real (kind=double), dimension (:,:) :: octantUnitDims
   real (kind=double), dimension (:)   :: rootOctantSizes
   real (kind=double), dimension (:,:) :: wignerOctantUnitDims
   real (kind=double), dimension (:)   :: rootWignerOctantSizes

   ! Define the local variables used in this subroutine.
   integer :: i ! loop counter
   integer :: numWignerBorders, numOctantBorders ! The number of borders needed
              !  to define the Wigner-Seitz cell and the octant of the given
              !   cell.
   integer :: currentIndex ! This is a simple counter to track where the
              ! values in the array and matrix should be stored.

   ! Determine the number of borders for each type of cell.
   numWignerBorders = size (rootWignerSizes)
   numOctantBorders = size (rootOctantSizes)

   ! Initialize the currentIndex marker
   currentIndex = 0

   do i = 1, numWignerBorders
      currentIndex = currentIndex + 1
      rootWignerOctantSizes(currentIndex)  = rootWignerSizes(i)
      wignerOctantUnitDims(:,currentIndex) = wignerUnitDims(:,i)
   enddo

   do i = 1, numOctantBorders
      currentIndex = currentIndex + 1
      rootWignerOctantSizes(currentIndex)  = rootOctantSizes(i)
      wignerOctantUnitDims(:,currentIndex) = octantUnitDims(:,i)
   enddo

end subroutine combineWignerAndOctants

subroutine compareOctantsAndWignerCells(wignerUnitDims,rootWignerSizes,&
      & octantUnitDims,rootOctantSizes,realVectors,maxWignerCell)

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension (:,:) :: wignerUnitDims
   real (kind=double), dimension (:)   :: rootWignerSizes
   real (kind=double), dimension (:,:) :: octantUnitDims
   real (kind=double), dimension (:)   :: rootOctantSizes
   real (kind=double), dimension (dim3,dim3) :: realVectors
   integer :: maxWignerCell

   ! Define the local variables used in this subroutine.
   integer :: i,j,k ! i = loop counter, j = 1st nested loop counter, ...
   integer :: numWignerBorders, numOctantBorders, numTotalBorders
   integer, allocatable, dimension (:) :: cellSegmentCounterLocal ! This
         ! will store the contributions of each border segment from the
         ! current combination of Wigner-Seitz cell and octant.
   real (kind=double), allocatable, dimension (:) :: rootWignerOctantSizes !
         ! The combination of the values from the rootWignerSizes and
         ! rootOctantSizes arrays with the Wigner-Seitz cells first, and the
         ! octants cells last.
   real (kind=double), allocatable, dimension (:,:) :: wignerOctantUnitDims !
         ! The combination of the values from the wignerUnitDims and
         ! octantUnitDims matrices with the Wigner-Seitz cells first, and the
         ! octants cells last.
   real (kind=double), allocatable, dimension (:)   :: currentSizes ! The
         ! rootWignerOctantSizes for the current choice of octant and the
         ! current choice of Wigner-Seitz cell.
   real (kind=double), dimension (dim3) :: octantTranslation ! This is a vector
         ! that is used to indicate the direction that the original octant
         ! should be translated to produce the current octant to study.
   real (kind=double), allocatable, dimension (:,:) :: unitDotMatrix ! This is
         ! a dot product matrix for every unit vector of wignerOctantUnitDims
         ! with every other wignerOctantUnitDims vector.  NOTE that this is
         ! NOT the same data structure as the unitDotMatrix referenced in the
         ! subroutine that calls this one.  That unitDotMatrix was not passed
         ! here and does not exist here due to scoping rules.

   ! Initialize the number of wigner cells needed.
   numWignerCellsNeeded(:) = 0

   ! Store the length of the rootWignerSizes and rootOctantSizes arrays for
   !   faster reference later on.
   numWignerBorders = size (rootWignerSizes)
   numOctantBorders = size (rootOctantSizes)
   numTotalBorders  = numWignerBorders + numOctantBorders

   ! The Wigner-Seitz cell and the octant of the given cell will be joined
   !   into the same array.  This will make some computations easier.  It is
   !   important to note that some of the vectors describing the two cells may
   !   be identical.  When the octant and the Wigner-Seitz cell are compared
   !   the second of the two identical vectors will be skipped so that a
   !   single border is not counted twice.

   ! We must first allocate sufficient space to hold the wignerOctant
   !   combination array and matrix
   allocate (wignerOctantUnitDims(dim3,numTotalBorders))
   allocate (rootWignerOctantSizes(numTotalBorders))

   ! Then we can combine the two arrays into one.
   call combineWignerAndOctants(wignerUnitDims,rootWignerSizes,&
                               &octantUnitDims,rootOctantSizes,&
                               &wignerOctantUnitDims,rootWignerOctantSizes)

   ! Recreate the unitDotMatrix from the wignerOctantUnitDims matrix.

   ! Allocate space for the resultant matrix first.  This has to be done
   !   outside the procedure call since the matrix has an allocatable attribute.
   allocate (unitDotMatrix(numTotalBorders,numTotalBorders))

   ! Then compute the matrix
   call dotProductMatrix(wignerOctantUnitDims,wignerOctantUnitDims,&
                        &unitDotMatrix)

   ! Allocate space for the currentSizes that will be modified through each
   !   iteration of the loops below.
   allocate (currentSizes(numTotalBorders))

   ! Allocate space to store the result of the next calculation.  That is the
   !   index numbers of the lattice points occupied by a Wigner-Seitz cell that
   !   are required to totally enclose each octant.  For now this allocation
   !   will be static since I can't image needing more than 50 for a single
   !   octant I think that this will be enough.  (Famous last words.)
   ! In the future I would like to dynamically allocate this.  That will
   !   require finding out the maximum number of Wigner-Seitz cells that any
   !   octant will need.  This will force some duplicate calculations.  (i.e.
   !   finding out how much space is needed is just as difficult (I think) as
   !   finding out the answers to the problem itself.)
   allocate (overlappingWignerCells(50,8))

   ! Begin a loop to test each of the eight octants.
   do i = 1, 8
      ! Define the translation amount for the origin this test octant.  The
      !   process here is to use a binary representation of the i loop index
      !   variable to define how much the vectors giving the octant borders
      !   should have their lengths modified.  The "how much" is either zero
      !   or an amount given by the system's given real space cell vectors.
      !   It is very simple to work through this with a simple cubic cell.
      do j = 1, dim3
         octantTranslation(j) = 0.0_double
         do k = 1, dim3
            octantTranslation(j) = octantTranslation(j) - 0.5_double * &
               & real((mod(i - 1,2**k) / 2**(k-1)) * realVectors(j,k),double)
         enddo
      enddo

      ! Shift the original octet by the amount determined above and save it in
      !   the currentSizes array.  This defines the lengths of the vectors to
      !   the current octant borders.
      do j = numWignerBorders + 1, numWignerBorders + numOctantBorders
         currentSizes(j) = rootWignerOctantSizes(j) + sum(octantTranslation(:)*&
                         & wignerOctantUnitDims(:,j))
      enddo

      ! Initiate another loop that will position a Wigner-Seitz cell at
      !   every lattice point (including the 0,0,0 point) and test for a
      !   shared volume between that Wigner-Seitz cell and the current octet
      !   of the given cell.
      do j = 1, numCellsReal

         ! Do not consider lattice points that are more than three maximum
         !   length Wigner-Seitz cells away since they will certainly not have
         !   a shared volume with any octet.
         if (cellSizesReal(j) > 3.0_double * cellSizesReal(maxWignerCell)) exit

         ! Shift the Wigner-Seitz cell to the current lattice point.
         do k = 1, numWignerBorders
            currentSizes(k) = rootWignerOctantSizes(k) + &
                  & sum(cellDimsReal(:,j) * wignerOctantUnitDims(:,k))
         enddo

         ! Compare the line formed by the intersection of the planes
         !   defined by any given pair of Wigner-Seitz or octet 
         !   wignerOctetDims vectors (normals) to the planes formed by
         !   every other wignerOctetDims vector.  If a segment of the
         !   line exists within every other wignerOctetDims vector, then
         !   that pair of either Wigner-Seitz or octet vectors will have
         !   their score incremented by one.  The scores are used later.
         !   NOTE that the integer (1) at the end is used to indicate that
         !   identical vectors should be searched for in the wignerOctetDims
         !   list of vectors.  In this execution of the intersection test
         !   there may be some identical vectors.

         ! First allocate space to store the counter for each cell.
         if (.not.allocated(cellSegmentCounterLocal)) then
            allocate(cellSegmentCounterLocal(size(currentSizes)))
            cellSegmentCounterLocal(:) = 0
         endif

         ! Then perform the intersection test.
         call intersectionTest(wignerOctantUnitDims,unitDotMatrix,currentSizes,&
               & cellSegmentCounterLocal,size(currentSizes),1)

         ! Check the results stored in the cellSegmentCounterLocal.  Any plane
         !   formed from the wignerOctetDims that has a score of 3 or more
         !   will increment a counter for the current Wigner-Seitz cell.
         !   If the counter is greater than four, then the Wigner-Seitz
         !   cell and current octet share a common volume.  Because of this
         !   the lattice point of the current Wigner-Seitz cell will be
         !   recorded for later searches for lattice vectors closest to an
         !   arbitrary vector.
         call sharedVolumeCheck(cellSegmentCounterLocal,j,i,&
               & numWignerCellsNeeded,overlappingWignerCells)

      enddo
   enddo

   ! Deallocate all the data structures that were defined here, but will not
   !   be used later.
   deallocate (wignerOctantUnitDims)
   deallocate (rootWignerOctantSizes)
   deallocate (unitDotMatrix)
   deallocate (currentSizes)
   deallocate (cellSegmentCounterLocal)

end subroutine compareOctantsAndWignerCells


subroutine getRealMagnitudes

   implicit none

   ! Define local variables.
   integer :: i

   ! Compute the magnitudes of a, b, and c.
   do i = 1, 3
      realMag(i) = sqrt(sum(realVectors(:,i)**2))
   enddo
end subroutine getRealMagnitudes


subroutine getRealUnitVectors

   implicit none

   ! Define local variables.
   integer :: i

   ! Compute the unit vectors of the real lattice.
   do i = 1, 3
      realUnitVectors(:,i) = realVectors(:,i) / realMag(i)
   enddo
end subroutine getRealUnitVectors


subroutine getRealAngles

   implicit none

   ! Define local variables.
   integer :: i
   integer :: index1, index2

   ! Compute the angles between the vectors (alpha,beta,gamma:  bc,ac,ab).
   do i = 1,3
      index1 = mod(i,3)+1
      index2 = mod(i+1,3)+1
      realAngles(i) = acos(dot_product(realVectors(:,index1),&
            & realVectors(:,index2)) / realMag(index1) / &
            & realMag(index2))
   enddo
end subroutine getRealAngles


subroutine getRealNormalVectors

   ! Import necessary modules.
   use O_MathSubs

   implicit none

   ! Define local variables.
   integer :: i
   integer :: index1, index2

   ! Compute the unit normal vectors of the real space cell that define the
   !   orientation of each plane (bc, ac, ab).  (This is used to help find the
   !   volume of space represented by one data point of a 3D mesh.)
   do i = 1,3
      index1 = mod(i,3)+1
      index2 = mod(i+1,3)+1
      call crossProduct(realUnitNormal(:,i),realUnitVectors(:,index1),&
            realUnitVectors(:,index2))
      realUnitNormal(:,i) = realUnitNormal(:,i) / sin(realAngles(i))
   enddo
end subroutine getRealNormalVectors


subroutine getRealFractStrideLength

   implicit none

   ! Define local variables.
   integer :: i

   ! Based on the number of mesh points in the 3D mesh, determine the
   !   fractional stride length.
   do i = 1, 3
      realFractStrideLength(i) = 1.0_double / real(numMeshPoints(i),double)
   enddo
end subroutine getRealFractStrideLength


subroutine getRealFractCrossArea

   implicit none

   ! Define local variables.
   integer :: i
   integer :: index1, index2

   ! Compute the fractional cross sectional areas (bc, ac, ab).  This is the
   !   cross sectional area of a unit in the mesh and essentially represents
   !   the amount of area that one data point represents.  NOTE:  The cross
   !   sectional area of ab is not influenced by the c vector's orientation.
   !   If a and b were perpendicular then the fractCrossArea(3) would simply be
   !   a*b.  However, they may be non-orthogonal with an angle between them
   !   (gamma) /= 90 degrees.  In this case the area is computed as the area
   !   for a parallelogram = a*b*sin(gamma).
   do i = 1,3
      index1 = mod(i,3)+1
      index2 = mod(i+1,3)+1
      realFractCrossArea(i) = &
            & realMag(index1) * realFractStrideLength(index1) * &
            & realMag(index2) * realFractStrideLength(index2) * &
            & sin(realAngles(i))
   enddo
end subroutine getRealFractCrossArea


subroutine getRealPlaneAngles

   implicit none

   ! Define local variables.
   integer :: i
   integer :: index1, index2

   ! Compute the angles between the unit normal vectors of the planes and the
   !   lattice unit vector not in the plane.  Angles are (bc,a ; ac,b ; ab,c)
   !   where the first is the plane (normal) and the second is the other unit
   !   lattice vector.  This angle is necessary because (for example) if the
   !   c lattice vector is *not* perpendicular to the ab plane then the volume
   !   represented by a mesh data point will *not* be area(ab)*mag(c).  Instead
   !   it will be area(ab)*mag(c)*sin(this angle).
   do i = 1,3
      index1 = mod(i,3)+1
      index2 = mod(i+1,3)+1
      realPlaneAngles(i) = &
            asin(dot_product(realUnitNormal(:,i),realUnitVectors(:,i)))
   enddo
!write (6,*) "planeAngles bc,a ac,b ab,c = ",planeAngles(:)*180.0_double/pi



end subroutine getRealPlaneAngles


subroutine cleanUpLattice

   implicit none

   deallocate (cellDimsReal)
   deallocate (cellSizesReal)
   if (allocated(cellDimsRecip)) then
      deallocate (cellDimsRecip)
      deallocate (cellSizesRecip)
   endif

   deallocate (overlappingWignerCells)

end subroutine cleanUpLattice


end module O_Lattice
