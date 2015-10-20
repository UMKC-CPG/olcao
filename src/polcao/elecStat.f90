module O_ElectroStatics

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public
   private :: neutralAndNuclearQPot, residualQ

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


subroutine makeElectrostatics()

   ! Use necessary modules.
   use O_TimeStamps

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Start making the matrices
   call timeStampStart(14)

   ! Create the real space matrices
   call neutralAndNuclearQPot()

   ! Create the reciprocal space matrices
   call residualQ()

   ! Log the date and time we end.
   call timeStampEnd(14)

end subroutine makeElectrostatics


subroutine neutralAndNuclearQPot()

   ! Include the definition modules we need.
   use HDF5
   use MPI
   use O_Kinds
   use O_Constants, only: dim3, smallThresh, pi
   use O_SetupElecStatHDF5, only: potAlphaOverlap_did, nonLocalNeutQPot_did, &
         & nonLocalNucQPot_did, localNeutQPot_did, localNucQPot_did, potPot, pot
   use O_Lattice, only: numCellsReal, cellDimsReal, logElecThresh, &
         & findLatticeVector
   use O_PotTypes, only: potTypes, maxNumPotAlphas
   use O_PotSites, only: numPotSites, potSites
   use O_Potential, only: potDim

   ! Import external subroutines.
   use O_MathSubs, only: erf
   use O_ParallelSubs

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the local loop control variables.
   integer :: i,j,k,l,m ! Loop variables.  loop1=i, nestedloop2=j ...
   integer :: hdferr

   real (kind=double), allocatable, dimension (:,:) :: nonLocalNeutQPot
   real (kind=double), allocatable, dimension (:)   :: nonLocalNucQPot
   real (kind=double), allocatable, dimension (:)   :: localNucQPot
   real (kind=double), allocatable, dimension (:,:) :: localNeutQPot
   real (kind=double), allocatable, dimension (:,:) :: potAlphaOverlap ! This
      ! matrix is potDim x potDim in size and is used for multiple purposes.
      ! It represents the overlap between all pairs of the Gaussian functions
      ! that are used for the electronic potential and valence charge density.
      ! The all pairs includes neighboring cells.  Note that the coefficients
      ! in front of the Gaussian functions are not included because they are
      ! still unknown for both the potential and valence charge density.


   ! Iteration dependent variables.  For each pair of potential sites these
   !   values will change when the potential type of either site changes.
   integer, dimension (2) :: currentType ! This stores the type number of the
         ! current i-loop and j-loop atoms.  (1) is i, and (2) is j.
   integer, dimension (2) :: currentNumAlphas ! This stores the number of
         ! alphas for potential type of the current i-loop and j-loop atoms.
         ! (1) is i, and (2) is j.
   integer :: currentMinNumAlpha ! This is the minimum of the two values in
         ! currentNumAlphas(:)
   integer :: currentMaxNumAlpha ! This is the maximum of the two values in
         ! currentNumAlphas(:)
   real (kind=double), allocatable, dimension (:,:) :: currentAlphas
   real (kind=double), allocatable, dimension (:,:) :: currentAlphasSqrt

   ! Local index record keeping variables.
   integer :: initAlphaIndex(2), finAlphaIndex(2) ! Indices that track which
         ! number alphas are being dealt with out of all alphas in the system.
         ! This is in reference to potDim.
   integer :: matrixIndex(2) ! Index record keeping variables.
   integer :: loopIndex(2) ! Index record keeping variables.
   integer :: arrayIndex ! Index record keeping variables.

   ! Flag variables to indicate some state in the calculation and control the
   !   program flow.
   integer :: aboveThreshold ! A simple 0,1 flag to test whether or not we have
         ! exceeded the requested level of accuracy of a given calculation set.

   ! Local intermediate matrices.  These are used to speed calculation of the
   !   above matrices.
   real (kind=double), allocatable, dimension (:,:) :: potOverlapCoeff ! This
         ! matrix is used primarily in the calculation of the potential
         ! overlap matrix.  It is (Pi/atom1,alpha1 + atom2,alpha2)**1.5.
   real (kind=double), allocatable, dimension (:,:) :: reducedPotAlpha ! This is
         ! a matrix describing a relation between the potential alphas of the
         ! current two potential sites.  It has a form analagous to the
         ! reduced mass of the two body central potential problem in mechanics:
         ! m1*m2/(m1+m2).  Here m1 is alpha1 of atom1 and m2 is alpha2 of atom2.
   real (kind=double), allocatable, dimension (:,:) :: reducedPotAlphaSqrt !
         ! This is just the square root of the above matrix.
   real (kind=double), allocatable, dimension (:,:) :: neutralQPotCoeff ! This
         ! matrix is Pi**3 / (atom1,alpha1 * atom2,alpha2)**1.5.  It is used
         ! for the summation of the local neutral valence potential-charge
         ! matrix and the nonlocal neutral valence potential-charge matrix.

   ! Local locational vectors and distance variables
   real (kind=double), dimension (dim3) :: potSiteSep ! This is simply the
         ! vector difference between two potential site location vectors.
   real (kind=double), dimension (dim3) :: latticeVector ! This vector points
         ! to the lattice point that is closest to the potSiteSep vector.
   real (kind=double), dimension (dim3) :: minSepVector  ! This is the vector
         ! difference between the potential site seperation vector (potSiteSep)
         ! and the closest lattice point vector (latticeVector).  It represents
         ! the minimal separation vector between two sites.
   real (kind=double) :: potSiteOffset ! This is the distance from potential
         ! site one to potential site two offset by the lattice Vector.
   real (kind=double) :: potSiteOffsetSqrd ! This is the square of the above
         ! value.

   ! Misc. factors that are not defined directly in the constants module.
   real (kind=double) :: piXX3       ! Pi**3
   real (kind=double) :: piXX32      ! Pi**1.5 = Pi**(3/2)
   real (kind=double) :: sqrtPi      ! Sqrt(Pi)
   real (kind=double) :: x2InvSqrtPi ! 2.0 / Sqrt(Pi)
   real (kind=double) :: x2pi        ! 2.0 * Pi

   ! Misc. intermediate calculations of different interactions.
   real (kind=double) :: nuclearFactor ! This can be made into a temp matrix
         ! like some other parameters too.
   real (kind=double) :: screenMagnitude  ! This depends on the potential site
         ! seperation and so cannot easily be made into a simple matrix.
         ! However, some of the multiplications can be removed so that only
         ! the effect of the superlattice potential site seperation is left.
   real (kind=double) :: coeff ! The same as the above applies to this.
   real (kind=double) :: potMagnitude

   ! Parallel variables
   integer :: minSite, maxSite ! The min and max potential site for this
                               ! process.
   integer :: mpiSize, mpiErr, mpiRank

   ! Allocate space for the coulomb potential matrices that are created here.
   allocate (nonLocalNeutQPot (potDim,potDim))
   allocate (nonLocalNucQPot  (potDim))
   allocate (localNucQPot     (potDim))
   allocate (localNeutQPot    (potDim,potDim))
   allocate (potAlphaOverlap  (potDim,potDim))

   ! Allocate space for the temporary intermediate matrices that will be used
   !   only during this subroutine.
   allocate (potOverlapCoeff     (maxNumPotAlphas,maxNumPotAlphas))
   allocate (reducedPotAlpha     (maxNumPotAlphas,maxNumPotAlphas))
   allocate (reducedPotAlphaSqrt (maxNumPotAlphas,maxNumPotAlphas))
   allocate (neutralQPotCoeff    (maxNumPotAlphas,maxNumPotAlphas))
   allocate (currentAlphas       (maxNumPotAlphas,2))
   allocate (currentAlphasSqrt   (maxNumPotAlphas,2))

   ! Initialize any factors not taken from the constants module.
   x2Pi        = 2.0_double * pi
   sqrtPi      = sqrt(pi)
   piXX3       = pi**3
   piXX32      = pi**1.5
   x2InvSqrtPi = 2.0_double / sqrtPi 

   ! Initialize all the matrices that will be filled in a cumulative or partial
   !   way.
   nonLocalNeutQPot (:,:) = 0.0_double
   nonLocalNucQPot  (:)   = 0.0_double
   localNucQPot     (:)   = 0.0_double
   localNeutQPot    (:,:) = 0.0_double
   potAlphaOverlap  (:,:) = 0.0_double
   reducedPotAlpha  (:,:) = 0.0_double

   ! Initialize the counters that track what alphas we are currently
   !   dealing with out of all the alphas in the whole system (potDim).
   initAlphaIndex(1) = 0
   finAlphaIndex(1) = 0

   ! Initialize other variables.
   currentNumAlphas(:) = 0
   currentType(:)      = 0
   currentMinNumAlpha  = 0
   currentMaxNumAlpha  = 0

   ! Many different things are computed in this subroutine.
   !
   ! A key aspect of OLCAO is that the charge density is cast into the form of
   !   a sum of atom centered Gaussian functions.  This form will allow for
   !   easy evaluation of the Coulomb electrostatic potential energy and the
   !   exchange-correlation energy.  Casting the valence charge density into
   !   the form of atom centered Gaussian functions follows a very similar
   !   procedure to the one used for the core charge density.  A set of
   !   equations is produced that relates the overlap of every Gaussian pair
   !   and the overlap between each Gaussian and the exact charge density
   !   represenstation (square of the solid state wave function).  In this
   !   subroutine, the overlap between each Gaussian pair including the
   !   Gaussians from image cells in the lattice are accumulated.  Note, that
   !   this matrix is real and symmetric and hence we really only need the
   !   upper triangle.  (This is the potAlphaOverlap matrix.) 

   ! To make the coulomb potential matrices we first create a smaller set of
   !   matrices.  This set of matrices is created for each possible combination
   !   of non-equivalent potential types.  The matrices are essentially factors
   !   that are constant during the lattice vector loop.  Each pair of
   !   potential sites runs through the hundred or so real space lattice
   !   vectors.

   ! It is also useful that these matrices only have to be calculated each time
   !   the potential sites change to a different potential type.  However,
   !   (although still useful for the lattice vector loop) this becomes less
   !   efficient when the number of potential types grows, especially with
   !   respect to the number of potential sites.

   ! In this light, there are a couple of ways that this can be sped up.
   !1) The matrices will be symmetric in the case that the two sets of
   !   potential alphas currently under consideration are the same. This
   !   would occur in the case of more than one potential type representing
   !   the same element (e.g. ca1 and ca2 from fluorapatite).
   !2) The matrices for every identical pair of elements will be identical.
   !   (e.g. The matrix formed from o1 and o2 from fluorapatite will be the
   !   exact same as the matrix formed by (o1,o3), and (o2,o3).  Also, the
   !   matrix formed by (ca1,o1) will be the same as the matrix set formed by
   !   (ca1,o2), (ca1,o3), (ca2,o1), (ca2,o2), and (ca3,o3).  etc.

   ! So, a little more complexity can reduce the number of times these must be
   !   calculated.  For small systems with only a few dozen types this is
   !   clearly of little consequence and not a prioriety yet.  However, if we
   !   later wish to do much larger systems with hundreds of types and
   !   thousands of sites such an improvement may be useful.

   ! The basic algorithm here is to consider the interaction between the
   !   Gaussians on each potential site (i) and the Gaussians on each other
   !   potential site (j) possibly in another replicated cell (k).  The
   !   exponential term in the Gaussian is used as a measure for the degree of
   !   interaction such that when this value reaches some threshold the level
   !   of interaction is considered to be negligable and any further
   !   interaction integrals can be disregarded.
   ! Some details are as follows:  The exponential alphas for the Gaussians
   !   at each site are arranged in monotonically increasing order.  The
   !   interactions between any pair of sites proceeds by consideration of the
   !   two longest ranged Gaussians first (one from each site).  If these two
   !   Gaussians are determined to have a non-negligable interaction then the
   !   interaction intergrals for this pair are computed.  Then the longest
   !   ranged Gaussian from site i and all the other Gaussians from the other
   !   site j are used for calculation.  Then the longest ranged Gaussian from
   !   site j and all the other Gaussians from the other site i are used for
   !   calculation.  (This is the first column, and then the first row of the
   !   interaction matrix between these two sites.)  Then the second longest
   !   ranged Gaussians of both sites are considered and the same procedure
   !   is used.  The second column and second row are computed (excluding the
   !   interaction with the longest ranged Gaussian because that was already
   !   done).  This proceeds down the diagonal of the matrix until at some
   !   point two Gaussians on the diagonal pass the threshold and have a
   !   negligable interaction.  At this point the algorithm stops for this
   !   pair of sites and the next site (j) in replicated cell (k) is chosen.
   ! In most cases the negligability limit will be reached before all the
   !   Gaussian pairs are considered.  This is because the diagonal elements
   !   (which are used to determine the negligability limit) represent shorter
   !   and shorter ranged Gaussians as the calculation proceeds.  Indeed, the
   !   the chance that the negligability limit will be reached before the final
   !   pair of Gaussians for a particular pair of sites is almost 100%.  The
   !   later Gaussians are very narrow and not broad at all.  The only time
   !   that *all* the Gaussians from a particular site pair will have non-
   !   negligable overlap is when the two sites are actually the same site in
   !   the original cell.  In this case all Gaussians will be the same and even
   !   the most sharp will overlap.
   ! In the past the program allowed for the use of auxiliary potential sites
   !   that were not associated with atomic sites.  This method has not been
   !   used in a very long time due to the arbitrary nature of when, how many,
   !   and where to add the auxiliary functions.  However, if they were added,
   !   then it would be possible to not reach the negligability limit before
   !   they were almost all considered (note the tricky aspect that the first
   !   alpha loop only goes to numAlphas - 1 or so).  Imagine only 3 or so very
   !   broad Gaussians at some auxiliary site and you could see how this might
   !   happen.  Now, in the future I may add the ability to have auxiliary
   !   functions back in to the OLCAO method, but at the moment that will be
   !   a lot more trouble than it is worth because it was rarely used, rather
   !   arbitrary, and there are some considerable computational efficiencies
   !   that have been developed *assuming* that all potential (charge) sites
   !   map to atomic sites too.  Still, we will keep this current algorithm
   !   because it isn't all that bad even for our now simpler case of no
   !   auxiliary functions.

   ! Initialize variables and data structures for parallel computation.
   call MPI_Comm_Size(MPI_COMM_WORLD,mpiSize,mpiErr)
   call MPI_Comm_Rank(MPI_COMM_WORLD,mpiRank,mpiErr)
   call loadBalMPI(numPotSites,minSite,maxSite,mpiRank,mpiSize)

   ! For parallel calculation, the processes that have a minSite>1 will need
   !   to compute the current status of a number of variables *as if* the
   !   algorithm had been running in serial. In the future, this should be
   !   replaced with a direct and intelligent assignment.
   if (minSite > 1) then
      do i = 1, minSite-1
         if (potSites(i)%firstPotType == 1) then
            ! Assign local copies of the potential type based variables.
            currentType(1)      = potSites(i)%potTypeAssn
            currentNumAlphas(1) = potTypes(currentType(1))%numAlphas
            currentAlphas(:currentNumAlphas(1),1) = &
               & potTypes(currentType(1))%alphas(:currentNumAlphas(1))
            currentAlphasSqrt(:currentNumAlphas(1),1) = &
               & sqrt(currentAlphas(:currentNumAlphas(1),1))

            ! Establish the beginning and ending indices for the set of
            !   alphas (potential terms) that we are dealing with out of
            !   all of the alphas in the whole system.
            initAlphaIndex(1)   = finAlphaIndex(1)
            finAlphaIndex(1)    = finAlphaIndex(1) + currentNumAlphas(1)
         endif
      enddo
   endif

   ! Begin computing over the assigned range of sites for this process.
   do i = minSite, maxSite
      if (potSites(i)%firstPotType == 1) then

         ! Assign local copies of the potential type assignment of the current
         !   potential site, and the number of alphas for that type.
         currentType(1)         = potSites(i)%potTypeAssn
         currentNumAlphas(1)    = potTypes(currentType(1))%numAlphas
         currentAlphas(:currentNumAlphas(1),1) = &
               & potTypes(currentType(1))%alphas(:currentNumAlphas(1))
         currentAlphasSqrt(:currentNumAlphas(1),1) = &
               & sqrt(currentAlphas(:currentNumAlphas(1),1))

         ! Update the indices for which alphas we are dealing with out of all
         !   the alphas in the whole system.
         initAlphaIndex(1) = finAlphaIndex(1)
         finAlphaIndex(1)  = finAlphaIndex(1) + currentNumAlphas(1)
      endif


      ! Initialize a pair of counters to track what alphas we are currently
      !   dealing with out of all the alphas in the whole system for the second
      !   nested loop of all alphas in the system.
      initAlphaIndex(2) = 0
      finAlphaIndex(2) = 0

      do j = 1, numPotSites
         if (potSites(j)%firstPotType == 1) then

            ! Assign local copies of the potential type assignment of the
            !   current potential site, and the number of alphas for that type.
            currentType(2)         = potSites(j)%potTypeAssn
            currentNumAlphas(2)    = potTypes(currentType(2))%numAlphas
            currentAlphas(:currentNumAlphas(2),2) = &
                  & potTypes(currentType(2))%alphas(:currentNumAlphas(2))
            currentAlphasSqrt(:currentNumAlphas(2),2) = &
                  & sqrt(currentAlphas(:currentNumAlphas(2),2))

            ! Update the indices for which alphas we are dealing with out of
            !   all the alphas in the whole system.
            initAlphaIndex(2) = finAlphaIndex(2)
            finAlphaIndex(2)  = finAlphaIndex(2) + currentNumAlphas(2)


            ! Identify some useful values for this loop iteration that will be
            !   used to control later loops and summations.
            currentMinNumAlpha = min(currentNumAlphas(1),currentNumAlphas(2))
            currentMaxNumAlpha = max(currentNumAlphas(1),currentNumAlphas(2))

            do k = 1, currentNumAlphas(2)

               ! Assign values for the matrix that contributes the alphas of
               !   exponential functions used to describe the potential
               !   overlap matrix, and the real space electrostatic non-local
               !   matrix, etc.  The equation is:
               !   (alpha1 * alpha2) / (alpha1 + alpha2) for each alpha in the
               !   list of alphas for the current types.
               reducedPotAlpha(:currentNumAlphas(1),k) = &
                     & currentAlphas(:currentNumAlphas(1),1) * &
                     & currentAlphas(k,2) / &
                     &(currentAlphas(:currentNumAlphas(1),1) + &
                     & currentAlphas(k,2))

               ! Assign values for the matrix that contributes to the local
               !  (bare-nuclear/neutral-valence) potential charge matrices.
               neutralQPotCoeff(:currentNumAlphas(1),k) = &
                     & piXX3 / (currentAlphas(:currentNumAlphas(1),1) * &
                     & currentAlphasSqrt(:currentNumAlphas(1),1) * &
                     & currentAlphas(k,2) * currentAlphasSqrt(k,2))

               ! Assign values for the matrix that defines coefficients that
               !   determine the values in the potential overlap matrix later.
               !   The equation is (pi/(alpha1 + alpha2)**1.5 for each alpha
               !   in the list of alphas for the current types.
               potOverlapCoeff(:currentNumAlphas(1),k) = &
                     & (pi/(currentAlphas(:currentNumAlphas(1),1) + &
                     & currentAlphas(k,2)))**1.5
            enddo

            ! Compute here the square root of the above reducedPotAlpha matrix.
            reducedPotAlphaSqrt(:currentNumAlphas(1),:currentNumAlphas(2)) = &
                  & sqrt(reducedPotAlpha(:currentNumAlphas(1),&
                  & :currentNumAlphas(2)))

         endif

         ! Determine the seperation between the two current potential sites.
         !   Note that since the two potential sites are both within the same
         !   system cell (unit cell whatever), the difference between them will
         !   always be within one unit cell of the octants of the 3-D cartesean
         !   coordinate system.
         potSiteSep(:) = potSites(i)%cartPos(:) - potSites(j)%cartPos(:)

         ! Find the origin (vector) of the superlattice site closest to the
         !   difference between the position vectors of two potential sites.
         !   Since the potSiteSep is always in one of the first octants of the
         !   3-D cartesean coordinate system the closest superlattice site will
         !   always be either the system origin, or the first lattice point
         !   encountered in any of the 26 possible directions.

         call findLatticeVector (potSiteSep,latticeVector)


         ! Determine the distance between the latticeVector calculated above,
         !   and the potSiteSep vector also calculated above.  This is the
         !   minimum separation distance between the two sites accounting for
         !   periodic boundary conditions.
         minSepVector(:) = potSiteSep(:) - latticeVector(:)


         do k = 1, numCellsReal
            potSiteOffsetSqrd = sum((minSepVector(:) - cellDimsReal(:,k))**2)
            potSiteOffset = sqrt(potSiteOffsetSqrd)

            ! Initialize a flag to check whether or not we are beyond the
            !   requested accuracy level.  This will only remain zero when the
            !   site i and site j are the same site and k = 1 (site j is in the
            !   original cell, not a replicated cell).
            aboveThreshold = 0

            ! Initiate a loop to fill in the non-local part and the gaussian
            !   local part. 
            do l = 1, currentMinNumAlpha - 1

               ! Define the strength of the overlap for the most distant
               !   reaching alphas.  (All later alpha comparisons will be
               !   against values larger than those in reducedPotAlpha(l,l).
               potMagnitude = reducedPotAlpha(l,l) * potSiteOffsetSqrd

               ! If this magnitude is outside the range of our requested
               !   accuracy level, then we don't need to calculate anymore.
               if (potMagnitude > logElecThresh) then
                  aboveThreshold = 1
                  exit
               endif
               do m = l, currentNumAlphas(2)

                  ! Identify the current matrix indices.
                  matrixIndex(1) = initAlphaIndex(1) + l
                  matrixIndex(2) = initAlphaIndex(2) + m

                  ! Determine if this current set will apply to the local
                  !   or the non-local part of the electroStatic matrix set.
                  if (abs(potSiteOffsetSqrd) >= smallThresh) then
                     nonLocalNeutQPot(matrixIndex(2),matrixIndex(1)) = &
                        & nonLocalNeutQPot(matrixIndex(2),matrixIndex(1)) + &
                        & neutralQPotCoeff(l,m) / potSiteOffset * &
                        & (- erf(currentAlphasSqrt(l,1) * potSiteOffset) + &
                        & (+ erf(reducedPotAlphaSqrt(l,m) * potSiteOffset)))
                  else
                     localNeutQPot(matrixIndex(2),matrixIndex(1)) = &
                        & localNeutQPot(matrixIndex(2),matrixIndex(1)) + &
                        & neutralQPotCoeff(l,m) * x2InvSqrtPi * &
                        & reducedPotAlphaSqrt(l,m)
                  endif

                  ! Calculate the electroStatic potential Overlap matrix ??
                  potAlphaOverlap(matrixIndex(2),matrixIndex(1)) = &
                     & potAlphaOverlap(matrixIndex(2),matrixIndex(1)) + &
                     & potOverlapCoeff(l,m) * exp(-reducedPotAlpha(l,m) * &
                     & potSiteOffsetSqrd)
               enddo

               ! This next loop is the same as the previous one except that
               !   the indices on the matrices and arrays are reversed.
               do m = l+1, currentNumAlphas(1)

                  ! Identify the current matrix indices.
                  matrixIndex(1) = initAlphaIndex(1) + m
                  matrixIndex(2) = initAlphaIndex(2) + l

                  ! Determine if this current set will apply to the local
                  !   or the non-local part of the electroStatic matrix set.
                  if (abs(potSiteOffsetSqrd) >= smallThresh) then
                     nonLocalNeutQPot(matrixIndex(2),matrixIndex(1)) = &
                        & nonLocalNeutQPot(matrixIndex(2),matrixIndex(1)) + &
                        & neutralQPotCoeff(m,l) / potSiteOffset * &
                        & (- erf(currentAlphasSqrt(m,1) * potSiteOffset) + &
                        & (+ erf(reducedPotAlphaSqrt(m,l) * potSiteOffset)))
                  else
                     localNeutQPot(matrixIndex(2),matrixIndex(1)) = &
                        & localNeutQPot(matrixIndex(2),matrixIndex(1)) + &
                        & neutralQPotCoeff(m,l) * x2InvSqrtPi * &
                        & reducedPotAlphaSqrt(m,l)
                  endif

                  ! Calculate the electroStatic potential Overlap matrix ??
                  potAlphaOverlap(matrixIndex(2),matrixIndex(1)) = &
                     & potAlphaOverlap(matrixIndex(2),matrixIndex(1)) + &
                     & potOverlapCoeff(m,l) * exp(-reducedPotAlpha(m,l) * &
                     & potSiteOffsetSqrd)
               enddo
            enddo

            ! Check if the above loop set was aborted or not.  If not, it means
            !   that we have to
            if (aboveThreshold == 0) then
               do m = l,currentMaxNumAlpha
                  if (currentNumAlphas(1) <= currentNumAlphas(2)) then

                     ! Identify the current loop indices.
                     loopIndex(1) = l
                     loopIndex(2) = m

                     ! Identify the current matrix indices.
                     matrixIndex(1) = initAlphaIndex(1) + l
                     matrixIndex(2) = initAlphaIndex(2) + m
                  else

                     ! Identify the current loop indices.
                     loopIndex(1) = m
                     loopIndex(2) = l

                     ! Identify the current matrix indices.
                     matrixIndex(1) = initAlphaIndex(1) + m
                     matrixIndex(2) = initAlphaIndex(2) + l
                  endif

                  ! Determine if this current set will apply to the local
                  !   or the non-local part of the electroStatic matrix set.
                  if (abs(potSiteOffsetSqrd) >= smallThresh) then
                     nonLocalNeutQPot(matrixIndex(2),matrixIndex(1)) = &
                        & nonLocalNeutQPot(matrixIndex(2),matrixIndex(1)) + &
                        & neutralQPotCoeff(loopIndex(1),loopIndex(2)) / &
                        & potSiteOffset * &
                        & (- erf(currentAlphasSqrt(loopIndex(1),1) * &
                        & potSiteOffset) + &
                        & (erf(reducedPotAlphaSqrt(loopIndex(1),loopIndex(2)) *&
                        & potSiteOffset)))
                  else
                     localNeutQPot(matrixIndex(2),matrixIndex(1)) = &
                        & localNeutQPot(matrixIndex(2),matrixIndex(1)) +&
                        & neutralQPotCoeff(loopIndex(1),loopIndex(2)) * &
                        & x2InvSqrtPi * &
                        & reducedPotAlphaSqrt(loopIndex(1),loopIndex(2))
                  endif

                  ! Calculate the electroStatic potential Overlap matrix ??
                  potAlphaOverlap(matrixIndex(2),matrixIndex(1)) = &
                     & potAlphaOverlap(matrixIndex(2),matrixIndex(1)) + &
                     & potOverlapCoeff(loopIndex(1),loopIndex(2)) * &
                     & exp(-reducedPotAlpha(loopIndex(1),loopIndex(2)) * &
                     & potSiteOffsetSqrd)
               enddo
            endif

            ! Begin another loop to calculate the gaussian screened nuclear
            !   charge integral vector.
            do m = 1, currentNumAlphas(1)

               ! Identify the index of the charge array to be filled
               arrayIndex = initAlphaIndex(1) + m

               ! Calculate the nuclear potential factor
               nuclearFactor = currentAlphas(m,1) + &
                     & potTypes(currentType(2))%nucAlpha

               ! Calculate the magnitude of the screening vector for this index?
               screenMagnitude = currentAlphas(m,1) * potSiteOffsetSqrd * &
                     & potTypes(currentType(2))%nucAlpha / nuclearFactor

               if (screenMagnitude < logElecThresh) then
                  if (abs(screenMagnitude) >= smallThresh) then

                     ! Caluclate the coefficient factor for the evaluation
                     !   that follows
                     coeff = currentAlphas(m,1) * potSiteOffset / &
                        & sqrt(nuclearFactor)

                     nonLocalNucQPot(arrayIndex) = &
                        & nonLocalNucQPot(arrayIndex) - &
                        & potTypes(currentType(2))%nucCharge * &
                        & piXX32 * exp(-screenMagnitude) * erf(coeff) / &
                        & coeff / nuclearFactor
                  else
                     nonLocalNucQPot(arrayIndex) = &
                        & nonLocalNucQPot(arrayIndex) - &
                        & potTypes(currentType(2))%nucCharge * &
                        & exp(-screenMagnitude) * x2pi / nuclearFactor
                  endif
               else
                  exit ! Exit with the screenmagnitude is sufficiently small.
               endif
            enddo
         enddo

         ! Determine the local bare nuclear potential charge matrix (array
         !   really since only the diagonal values are non-zero).

         ! If the seperation is very small then the indices are referring to
         !   the same atom and so we include this contribution.
         if (dot_product(potSiteSep(:),potSiteSep(:)) < smallThresh) then

            ! This assignment is simply -2 * pi * (current nuclear charge) / 
            !   this atom's potential alphas.  The +1 is used in the range
            !   since the init and fin AlphaIndex values are low by one since
            !   that make life much clearer up above.

            localNucQPot(initAlphaIndex(1)+1:finAlphaIndex(1)) = &
               & localNucQPot(initAlphaIndex(1)+1:finAlphaIndex(1)) - &
               & x2pi * potTypes(currentType(1))%nucCharge / &
               & potTypes(currentType(1))%alphas(:currentNumAlphas(1))
         endif
      enddo
   enddo

   ! Deallocate the variables that will not be used later.
   deallocate (potOverlapCoeff)
   deallocate (reducedPotAlpha)
   deallocate (reducedPotAlphaSqrt)
   deallocate (currentAlphas)
   deallocate (currentAlphasSqrt)
   deallocate (neutralQPotCoeff)

   ! Accumulate the results using MPI_REDUCE. In the future the bounds of the
   !   local buffer that were actually used should be sent explicitly so that
   !   we don't spend extra time sending and accumulating a bunch of zeros.
   call MPI_BARRIER(MPI_COMM_WORLD,mpiErr)
   call MPI_REDUCE(nonLocalNeutQPot(:,:),nonLocalNeutQPot(:,:),potDim*potDim,&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpiErr)
   call MPI_BARRIER(MPI_COMM_WORLD,mpiErr)
   call MPI_REDUCE(localNeutQPot(:,:),localNeutQPot(:,:),potDim*potDim,&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpiErr)
   call MPI_BARRIER(MPI_COMM_WORLD,mpiErr)
   call MPI_REDUCE(nonLocalNucQPot(:),nonLocalNucQPot(:),potDim,&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpiErr)
   call MPI_BARRIER(MPI_COMM_WORLD,mpiErr)
   call MPI_REDUCE(localNucQPot(:),localNucQPot(:),potDim,&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpiErr)
   call MPI_BARRIER(MPI_COMM_WORLD,mpiErr)
   call MPI_REDUCE(potAlphaOverlap(:,:),potAlphaOverlap(:,:),potDim*potDim,&
         & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpiErr)

   call MPI_BARRIER(MPI_COMM_WORLD,mpiErr)
   if (mpiRank==0) then
     call h5dwrite_f(potAlphaOverlap_did,H5T_NATIVE_DOUBLE, &
           & potAlphaOverlap(:,:),potDims2,hdferr)
     if (hdferr /= 0) stop 'Failed to write pot alpha overlap'
     call h5dwrite_f(nonLocalNeutQPot_did,H5T_NATIVE_DOUBLE, &
           & nonLocalNeutQPot(:,:),potDims2,hdferr)
     if (hdferr /= 0) stop 'Failed to write non local neut q pot'
     call h5dwrite_f(nonLocalNucQPot_did,H5T_NATIVE_DOUBLE, &
           & nonLocalNucQPot(:),potDims1,hdferr)
     if (hdferr /= 0) stop 'Failed to write non local nuc q pot'
     call h5dwrite_f(localNeutQPot_did,H5T_NATIVE_DOUBLE, &
           & localNeutQPot(:,:),potDims2,hdferr)
     if (hdferr /= 0) stop 'Failed to write local neut q pot'
     call h5dwrite_f(localNucQPot_did,H5T_NATIVE_DOUBLE, &
           & localNucQPot(:),potDims1,hdferr)
     if (hdferr /= 0) stop 'Failed to write local nuc q pot'
   endif
   ! Dellocate the variables from the global data structures that will not
   !   be used later.
   deallocate (nonLocalNeutQPot)
   deallocate (nonLocalNucQPot)
   deallocate (localNucQPot)
   deallocate (localNeutQPot)
   deallocate (potAlphaOverlap)

end subroutine neutralAndNuclearQPot

subroutine residualQ()

   ! Include the modules we need
   use HDF5
   use MPI
   use O_Kinds
   use O_Constants, only: dim3, pi, smallThresh
   use O_SetupElecStatHDF5, only: nonLocalResidualQ_did, potTypesPot
   use O_Lattice, only: numCellsRecip, cellSizesRecip, cellDimsRecip, &
         & logElecThresh, realCellVolume, numCellsReal, cellDimsReal, &
         & findLatticeVector
   use O_PotSites, only: numPotSites, potSites
   use O_PotTypes, only: numPotTypes, potTypes, maxNumPotAlphas, minPotAlpha
   use O_Potential, only: potDim

   ! Import external subroutines.
   use O_MathSubs, only: erf
   use O_ParallelSubs

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the local loop control variables.
   integer :: i,j,k,l ! Loop variables.  loop1=i, nestedloop2=j ...
   integer :: hdferr

   ! Primary data structure for the results of this subroutine.
   real (kind=double), allocatable, dimension (:,:) :: nonLocalResidualQ

   ! Iteration dependent variables.  For each pair of potential sites these
   !   values will change when the potential type of either site changes.
   integer, dimension (2) :: currentType ! This stores the type number of the
         ! current i-loop and j-loop atoms.  (1) is i, and (2) is j.
   integer, dimension (2) :: currentNumAlphas ! This stores the number of
         ! alphas for potential type of the current i-loop and j-loop atoms.
         ! (1) is i, and (2) is j.
   integer :: currentElement
   integer :: lastElement
   real (kind=double), allocatable, dimension (:,:) :: currentAlphas
   real (kind=double), allocatable, dimension (:,:) :: currentAlphasSqrt
   real (kind=double), allocatable, dimension (:,:) :: currentAlphasFactor

   ! Local index record keeping variables.
   integer :: initAlphaIndex, finAlphaIndex ! Indices that track which
         ! number alphas are being dealt with out of all alphas in the system.
   integer :: matrixIndex(2) ! Index record keeping variables.

   ! Local intermediate matrices.  These are used to speed calculation of the
   !   above matrices.
   real (kind=double), allocatable, dimension (:,:) :: reducedPotAlphaSqrt

   ! Local locational vectors and distance variables
   real (kind=double), dimension (dim3) :: potSiteSep ! This is simply the
         ! vector difference between two potential site location vectors.
   real (kind=double), dimension (dim3) :: latticeVector ! This vector points
         ! to the lattice point that is closest to the potSiteSep vector.
   real (kind=double), dimension (dim3) :: minSepVector ! This is the vector
         ! difference between the potential site seperation vector (potSiteSep)
         ! and the closest lattice point vector (latticeVector).  It represents
         ! the minimum separation vector between two sites.
   real (kind=double) :: potSiteOffset ! This is the distance from potential
         ! site 1 to potential site 2 offset by the lattice Vector.
   real (kind=double) :: potSiteOffsetSqrd ! This is the square of the above
         ! value.

   ! Misc. intermediate factors and coefficients.
   real (kind=double) :: minReducedAlpha ! This is the relation between the
         ! current alpha being calculated and the minimal alpha of the whole
         ! system.  It is the same type of calculation as the
         ! reducedPotAlphaSqrt described above (except without the sqrt).
   real (kind=double) :: inv4MinReducedAlpha ! 1/4 * 1/minReducedAlpha
   real (kind=double) :: alphaFactor ! Like the minReducedAlpha except that
         ! this is just a simple multiplication of the two alphas.
   real (kind=double) :: alphaSelfFactor ! This is the above alphaFactor *
         ! x2InvSqrtPi * (reducedPotAlphaSqrt for the current alpha - the
         ! minReducedAlpha from above)
   real (kind=double) :: alphaVolumeFactor ! This is the above alphaFactor *
         ! 4 * pi / real space cell volume.
   real (kind=double) :: logElecThreshSqrt
   real (kind=double), allocatable, dimension (:) :: phase
   real (kind=double), allocatable, dimension (:) :: minReducedAlphaSqrt ! Sqrt
         ! of the above minReducedAlpha for each potential alpha of a type.
   real (kind=double), allocatable, dimension (:,:) :: recipCellExp

   ! Misc. factors that are derived from values in the constants module.
   real (kind=double) :: piXX3       ! Pi**3
   real (kind=double) :: piXX32      ! Pi**1.5 = Pi**(3/2)
   real (kind=double) :: sqrtPi      ! Sqrt(Pi)
   real (kind=double) :: x2InvSqrtPi ! 2.0 / sqrt(pi)

   ! Parallel Variables
   integer :: minSite,maxSite
   integer :: ga_nonLocResidQ
   integer :: mpiSize,mpiRank,mpiErr

   ! Allocate space for the coulomb potential matrices that are created here.
   allocate (nonLocalResidualQ (numPotTypes,potDim))

   ! Allocate space for the temporary matrices that will be used only during
   !   this subroutine.
   allocate (reducedPotAlphaSqrt (maxNumPotAlphas,maxNumPotAlphas))
   allocate (currentAlphas       (maxNumPotAlphas,2))
   allocate (currentAlphasSqrt   (maxNumPotAlphas,2))
   allocate (currentAlphasFactor (maxNumPotAlphas,2))
   allocate (recipCellExp        (numCellsRecip,maxNumPotAlphas))
   allocate (minReducedAlphaSqrt (maxNumPotAlphas))
   allocate (phase               (numCellsRecip))

   ! Initialize the matrix to zero since the process to fill it is a
   !   cumulative summation.
   nonLocalResidualQ (:,:) = 0.0_double

   ! Initialize the pi factor and logElecThreshSqrt
   piXX3             = pi**3
   piXX32            = pi**1.5
   sqrtPi            = sqrt(Pi)
   x2InvSqrtPi       = 2.0_double / sqrtPi
   logElecThreshSqrt = sqrt(logElecThresh)

   ! Initialize other variables.
   currentNumAlphas(:) = 0
   currentType(:)      = 0

   ! Initialize a pair of counters to track what alphas we are currently
   !   dealing with out of all the alphas in the whole system.
   initAlphaIndex = 0
   finAlphaIndex  = 0

   ! Initialize the current element to zero.
   currentElement = 0

   ! Initialize the parallel data structures and variables.
   call MPI_Comm_Size(MPI_COMM_WORLD,mpiSize,mpiErr)
   call MPI_Comm_Rank(MPI_COMM_WORLD,mpiRank,mpiErr)
   call loadBalMPI(numPotSites,minSite,maxSite,mpiRank,mpiSize)

   ! For parallel calculation, the processes that have a minSite>1 will need
   !   to compute the current status of a number of variables *as if* the
   !   algorithm had been running in serial. In the future, this should be
   !   replaced with a direct and intelligent assignment.
   if (minSite > 1) then
      do i = 1, minSite-1
         if (potSites(i)%firstPotType == 1) then
            ! Assign local copies of the potential type based variables.
            currentType(1)      = potSites(i)%potTypeAssn
            currentNumAlphas(1) = potTypes(currentType(1))%numAlphas
            currentAlphasSqrt(:currentNumAlphas(1),1) = &
                  & sqrt(currentAlphas(:currentNumAlphas(1),1))
            currentAlphasFactor(:currentNumAlphas(1),1) = piXX32 / &
                  & (currentAlphas(:currentNumAlphas(1),1) * &
                  & currentAlphasSqrt(:currentNumAlphas(1),1))
            lastElement    = currentElement
            currentElement = potTypes(currentType(1))%elementID

            ! Establish the beginning and ending indices for the set of
            !   alphas (potential terms) that we are dealing with out of
            !   all of the alphas in the whole system.
            initAlphaIndex(1)   = finAlphaIndex(1)
            finAlphaIndex(1)    = finAlphaIndex(1) + currentNumAlphas(1)

            if (lastElement /= currentElement) then
               do j = 1, currentNumAlphas(1)
                  ! Get the potential alpha overlap between the current alpha
                  !   and the minimal potential alpha of the whole system.
                  minReducedAlpha = currentAlphas(j,1) * minPotAlpha / &
                        & (currentAlphas(j,1) + minPotAlpha)

                  ! Square root of the above
                  minReducedAlphaSqrt(j) = sqrt(minReducedAlpha)

                  ! One fourth of the inverse of minReducedAlpha
                  inv4MinReducedAlpha = 0.25_double / minReducedAlpha


                  ! Compute the exp decay factor for the reciprocal space cells.
                  recipCellExp(2:,j) = exp(-cellSizesRecip(2:) * &
                        & inv4MinReducedAlpha) / cellSizesRecip(2:)
               enddo
            endif
         endif
      enddo
   endif

   ! Begin computing over the assigned range of sites for this process.
   do i = minSite, maxSite

      if (potSites(i)%firstPotType == 1) then

         ! Assign local copies of the potential type assignment of the current
         !   potential site, and the number of alphas for that type.
         currentType(1)           = potSites(i)%potTypeAssn
         currentNumAlphas(1)      = potTypes(currentType(1))%numAlphas
         currentAlphas(:currentNumAlphas(1),1) = &
               & potTypes(currentType(1))%alphas(:currentNumAlphas(1))
         currentAlphasSqrt(:currentNumAlphas(1),1) = &
               & sqrt(currentAlphas(:currentNumAlphas(1),1))
         currentAlphasFactor(:currentNumAlphas(1),1) = piXX32 / &
               & (currentAlphas(:currentNumAlphas(1),1) * &
               & currentAlphasSqrt(:currentNumAlphas(1),1))
         lastElement    = currentElement
         currentElement = potTypes(currentType(1))%elementID

         ! Update the indices for which alphas we are dealing with out of all
         !   the alphas in the whole system.
         initAlphaIndex = finAlphaIndex
         finAlphaIndex  = finAlphaIndex + currentNumAlphas(1)

         ! Compute key variables for the alphas of this potential type.  This
         !   only needs to be done when the alphas actually change (which can
         !   only happen when the element changes).

         if (lastElement /= currentElement) then
            do j = 1, currentNumAlphas(1)

               ! Get the potential alpha overlap between the current alpha and
               !   the minimal potential alpha of the whole system.
               minReducedAlpha = currentAlphas(j,1) * minPotAlpha / &
                            &(currentAlphas(j,1) + minPotAlpha)

               ! Square root of the above
               minReducedAlphaSqrt(j) = sqrt(minReducedAlpha)

               ! One fourth of the inverse of minReducedAlpha
               inv4MinReducedAlpha = 0.25_double / minReducedAlpha


               ! Compute the exp decay factor for the reciprocal space cells.
               recipCellExp(2:,j) = exp(-cellSizesRecip(2:) * &
                     & inv4MinReducedAlpha) / cellSizesRecip(2:)
            enddo
         endif

      endif

      ! Initialize the matrix index counter for the second nested loop.
      matrixIndex(2) = 0

      do j = 1, numPotSites
         if (potSites(j)%firstPotType == 1) then

            ! Increment the matrixIndex(2) value
            matrixIndex(2) = matrixIndex(2) + 1

            ! Assign local copies of the potential type assignment of the
            !   current potential site, and the number of alphas for that type.
            currentType(2)           = potSites(j)%potTypeAssn
            currentNumAlphas(2)      = potTypes(currentType(2))%numAlphas
            currentAlphas(:currentNumAlphas(2),2) = &
                  & potTypes(currentType(2))%alphas(:currentNumAlphas(2))
            currentAlphasSqrt(:currentNumAlphas(2),2) = &
                  & sqrt(currentAlphas(:currentNumAlphas(2),2))
            currentAlphasFactor(:currentNumAlphas(2),2) = piXX32 / &
                  & (currentAlphas(:currentNumAlphas(2),2) * &
                  & currentAlphasSqrt(:currentNumAlphas(2),2))

            do k = 1, currentNumAlphas(2)

               ! Assign values for the matrix that contributes the alphas of
               !   exponential functions used to describe the potential
               !   overlap matrix, and the real space electrostatic non-local
               !   matrix, etc.  The equation is:
               !   (alpha1 * alpha2) / (alpha1 + alpha2) for each alpha in the
               !   list of alphas for the current types.
               reducedPotAlphaSqrt(:currentNumAlphas(1),k) = sqrt( &
                  & currentAlphas(:currentNumAlphas(1),1) * &
                  & currentAlphas(k,2) / &
                  &(currentAlphas(:currentNumAlphas(1),1) + &
                  & currentAlphas(k,2)))
            enddo
         endif

         ! Determine the seperation between the two current potential sites.
         !   Note that since the two potential sites are both within the same
         !   system cell (unit cell whatever), the difference between them will
         !   always be within one of unit cell the octants of the 3-D cartesean
         !   coordinate system.
         potSiteSep(:) = potSites(i)%cartPos(:) - potSites(j)%cartPos(:)

         ! Compute the phase factors for the reciprocal space lattice.
         do k = 1, numCellsRecip
            phase(k) = cos(sum(cellDimsRecip(:,k) * potSiteSep(:)))
         enddo

         ! Find the origin (vector) of the superlattice site closest to the
         !   difference between the position vectors of two potential sites.
         !   Since the potSiteSep is always in one of the first octants of the
         !   3-D cartesean coordinate system the closest superlattice site will
         !   always be either the system origin, or the first lattice point
         !   encountered in any of the 26 possible directions.
         call findLatticeVector (potSiteSep,latticeVector)

         ! Determine the distance between the latticeVector calculated above,
         !   and the potSiteSep vector also calculated above.
         minSepVector(:) = potSiteSep(:) - latticeVector(:)

         do k = 1, currentNumAlphas(1)

            ! Set the value for the first matrix index of the reciprocal
            !   electrostatic matrix.
            matrixIndex(1) = initAlphaIndex + k

            ! The next set of little assignment statements is designed to make
            !   the loops over real space and reciprocal space cells as fast
            !   as possible by eliminating any factors or coefficients that
            !   are independent of the potential site seperation over all
            !   potential sites in all the replicated cells.

            ! Calculate an alpha interaction factor with the current second
            !   nested loop minimal alpha.
            alphaFactor = currentAlphasFactor(k,1) * &
                        & currentAlphasFactor(1,2)

            ! Calculate the coefficient used for the case of self reference in
            !   the real part of the ewald summation.
            alphaSelfFactor = x2InvSqrtPi * alphaFactor* &
                            & (reducedPotAlphaSqrt(k,1)-minReducedAlphaSqrt(k))

            ! Calculate a relation between the above factor and the cell
            !   volume in real space.
            alphaVolumeFactor = 4.0_double * pi * alphaFactor / realCellVolume

            ! Real part of Ewald Summation
            do l = 1, numCellsReal
               potSiteOffsetSqrd = sum((minSepVector(:) - cellDimsReal(:,l))**2)

               if (potSiteOffsetSqrd < smallThresh) then
                  nonLocalResidualQ(matrixIndex(2),matrixIndex(1)) = &
                     & nonLocalResidualQ(matrixIndex(2),matrixIndex(1)) + &
                     & alphaSelfFactor
               else
                  potSiteOffset = sqrt(potSiteOffsetSqrd)
                  if (minReducedAlphaSqrt(k) * potSiteOffset > &
                        & logElecThreshSqrt .and. reducedPotAlphaSqrt(k,1) * &
                        & potSiteOffset > logElecThreshSqrt) then
                     exit
                  else
                     nonLocalResidualQ(matrixIndex(2),matrixIndex(1)) = &
                        & nonLocalResidualQ(matrixIndex(2),matrixIndex(1)) + &
                        & alphaFactor * (- erf(minReducedAlphaSqrt(k) * &
                        & potSiteOffset) + (erf(reducedPotAlphaSqrt(k,1) * &
                        & potSiteOffset))) / potSiteOffset
                  endif
               endif
            enddo

!            ! Zero point value for reciprocal summation
!            nonLocalResidualQ(matrixIndex(2),matrixIndex(1)) = &
!               & nonLocalResidualQ(matrixIndex(2),matrixIndex(1)) - &
!               & alphaVolumeFactor * (0.25 / minPotAlpha)

            ! The second term (in the parenthesis and after the sum) is the
            !   zero point value for the reciprocal summation.
            nonLocalResidualQ(matrixIndex(2),matrixIndex(1)) = &
                  & nonLocalResidualQ(matrixIndex(2),matrixIndex(1)) + &
                  & alphaVolumeFactor * (sum(recipCellExp(2:numCellsRecip,k) * &
                  & phase(2:numCellsRecip)) - (0.25_double / minPotAlpha))


!            do l = 2, numCellsRecip
!               nonLocalResidualQ(matrixIndex(2),matrixIndex(1)) = &
!                     & nonLocalResidualQ(matrixIndex(2),matrixIndex(1)) + &
!                     & recipCellExp(l,k) * alphaVolumeFactor * phase(l)
!            enddo
         enddo
      enddo

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
   write (20,*)

   ! Deallocate the variables that will not be used later.
   deallocate (reducedPotAlphaSqrt)
   deallocate (currentAlphas)
   deallocate (currentAlphasSqrt)
   deallocate (currentAlphasFactor)
   deallocate (cellDimsRecip)
   deallocate (cellSizesRecip)
   deallocate (phase)
   deallocate (recipCellExp)
   deallocate (minReducedAlphaSqrt)

   ! Accumulate the results into the distributed global array. In the future
   !   the bounds of the local buffer that were actually used should be sent
   !   to the ga_accumulation subroutine so that we don't spend extra time
   !   sending and accumulating a bunch of zeros.
   call MPI_BARRIER(MPI_COMM_WORLD,mpiErr)
   call MPI_REDUCE(nonLocalResidualQ(:,:),nonLocalResidualQ(:,:),&
         & numPotTypes*potDim,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,&
         & mpiErr)

   call MPI_BARRIER(MPI_COMM_WORLD,mpiErr)
   if (mpiRank==0) then
      ! Write to disk the matrices that will be used by main later.
      call h5dwrite_f(nonLocalResidualQ_did,H5T_NATIVE_DOUBLE, &
            & nonLocalResidualQ(:,:),potTypesPot,hdferr, &
            & xfer_prp=pot_xferpid)
      if (hdferr /= 0) stop 'Failed to write non local residual charge'
   endif

   ! Deallocate matrix that will not be used later.
   deallocate (nonLocalResidualQ)

end subroutine residualQ

end module O_ElectroStatics
