module O_LOEN

   ! Import necessary modules.
   use O_KINDS

   ! Make sure that no variables are accidentally declared.
   implicit none

   ! Define module data.
   complex (kind=double), allocatable, dimension(:,:) :: bsComp ! Bispec comp.
   complex (kind=double), allocatable, dimension(:) :: bsCompSum
   real (kind=double), allocatable, dimension(:,:) :: neighDist
   real (kind=double), allocatable, dimension(:,:) :: neighSwitchFn ! f_c
   real (kind=double), allocatable, dimension(:,:,:) :: neighxyzCoords
   real (kind=double), allocatable, dimension(:,:,:) :: neighttpCoords
   integer, allocatable, dimension(:) :: neighCount ! Neighbors of each pot.
   integer, allocatable, dimension(:,:) :: neighSite ! Neigh site idx #s.
   integer, allocatable, dimension(:,:) :: neighCell ! Neigh cell idx #s.

   contains

subroutine createNeighborList

   ! Import necessary modules.
   use O_Kinds
   use O_Constants
   use O_Input
   use O_Lattice
   use O_PotSites, only : numPotSites
   use O_Basis, only : initializePotSite


   ! Define local variables for logging and loop control.
   integer :: i, j, k

   ! Site specific variables that change with each site pair loop iteration.
   integer,            dimension (2) :: currentPotTypes
   integer,            dimension (2) :: currentElements
   real (kind=double), dimension (dim3,2) :: currentPositions
   real (kind=double), dimension (2) :: currentZFactors  ! Nuclei charges.
   real (kind=double), dimension (2) :: currentNucAlphas ! Exp alphas.

   ! Local position and direction vectors and radii.
   real (kind=double), dimension (dim3) :: latticeVector ! Vector to lattice
         ! point closest to the difference between the unit cell positions for
         ! site 1 and site 2.
   real (kind=double), dimension (dim3) :: shiftedPotPos ! The position of
         ! site 2 shifted to each relevant lattice point.
   real (kind=double) :: potSiteSepSqrd ! The square of the minimum distance
         ! seperating site 1 and site 2 according to their unit cell positions
         ! shifted by the lattice point closest to their difference.
   real (kind=double) :: shiftedPotSiteSep ! The seperation distance between
         ! pot site 1 and the shifted position of pot site 2.
   real (kind=double) :: shiftedPotSiteSepSqrd ! The seperation distance
         ! between pot site 1 and the shifted position of pot site 2 squared.
   real (kind=double) :: maxLatticeRadius ! Maximum radius beyond which no
         ! lattice points will be considered for integration.

   ! Additional local variables.
   real (kind=double) :: cutoffLoEnSqrd

   cutoffLoEnSqrd = cutoffLoEn**2

   ! Allocate space to hold data.
   allocate(neighCount(numPotSites))
   allocate(neighSite(maxNumNeigh,numPotSites))
   allocate(neighCell(maxNumNeigh,numPotSites))
   allocate(neighDist(maxNumNeigh,numPotSites))
   allocate(neighSwitchFn(maxNumNeigh,numPotSites))
   allocate(neighxyzCoords(3,maxNumNeigh,numPotSites))
   allocate(neighttpCoords(3,maxNumNeigh,numPotSites))

   ! Initialize the count of the number of neighbors within the cutoff for
   !   each possible central site.
   neighCount(:) = 0

   do i = 1, numPotSites

      ! Obtain local copies of key data from larger global data structures for
      !   the first looped atom.
      call initializePotSite(i,currentPotTypes(1),currentElements(1),&
            & currentZFactors(1),currentNucAlphas(1),currentPositions(:,1))

      ! Now, we need to compute a list of the sites within the cutoff radius.
      !   We will loop over all other sites in all replicated cells (until
      !   the cutoff has been reached).

      ! Begin a loop over all sites in the system. Note that we *do*
      !   double count because we need to completely quantify the local
      !   environment of each site.
      do j = 1, numPotSites

         ! Obtain local copies of key data from larger global data structures
         !   for the second looped atom.
         call initializePotSite(j,currentPotTypes(2),currentElements(2),&
               & currentZFactors(2),currentNucAlphas(2),currentPositions(:,2))

         ! Find the lattice point closest to the difference between the two
         !   atom sites. (This will allow us to get the smallest distance
         !   between the two sites so we can easily check if we need to check
         !   more cells of the periodic superlattice.)
         call findLatticeVector((currentPositions(:,1) &
               & - currentPositions(:,2)), latticeVector)

         ! Determine the square of the minimum seperation distance between the
         !   two sites.
         potSiteSepSqrd = sum((currentPositions(:,1) &
               & - currentPositions(:,2) - latticeVector(:))**2)

         ! Determine if this atom is within the current cutoff. Cycle if not.
         !   If we pass this "if" statement, that means this atom will
         !   certainly need to be included in the list of neighbor atoms.
         if (potSiteSepSqrd > cutoffLoEnSqrd) cycle

         ! Determine the maximum radius beyond which no lattice point will be
         !   considered to contribute to the overlap integral for this site 
         !   pair.  (This is the law of cosines c^2 = a^2 + b^2 + 2ab*cos(g)
         !   with g = gamma = angle between a and b.  Here potSiteSepSqrd=a^2,
         !   cutoffLoEnSqrd = b^2, g=0.)
         maxLatticeRadius = potSiteSepSqrd + cutoffLoEnSqrd + 2.0_double*&
                          & sqrt(potSiteSepSqrd * cutoffLoEnSqrd)

         ! Check if we have to check more lattice points than were determined
         !   to be needed by comparing the maxlatticeRadius to the maximum
         !   radius of the lattice points defined earlier.
         if (maxLatticeRadius > cellSizesReal(numCellsReal)) then
            write (20,*) 'More lattice points needed for this bispectrum comp.'
            write (20,*) 'maxLatticeRadius requested=',maxLatticeRadius
            write (20,*) 'max available=',cellSizesReal(numCellsReal)
            stop
         endif

         ! ! Begin a loop over all the lattice points to shift the position of
         !   atom number 2 to all the replicated cells.
         do k = 1, numCellsReal

            ! Exit the loop when we have exceeded the necessary number of
            !   lattice points based on distance.
            if (cellSizesReal(k) > maxLatticeRadius) exit

            ! Obtain the position of site #2 shifted by the current lattice.
            shiftedPotPos(:) = currentPositions(:,2) + latticeVector(:) + &
                  & cellDimsReal(:,k)

            ! Obtain the square of the magnitude of the seperation vector
            !   between site 1 and the shifted position of site 2.
            shiftedPotSiteSepSqrd = sum ((currentPositions(:,1) - &
                  & shiftedPotPos(:))**2)
!write(6,*) "i,j,k,sepSqrd=",i,j,k,shiftedPotSiteSepSqrd
            ! Determine if this shifted pot position puts it outside of the
            !   above determined negligability limit for this site pair.
            if (shiftedPotSiteSepSqrd > cutoffLoEnSqrd) cycle

            ! Ignore cases where the shifted pot position and the current
            !   pot position are the same. (I.e., there is no need to include
            !   self-interaction for the local environment analysis.)
            if (shiftedPotSiteSepSqrd < smallThresh) cycle
!write(6,*) "Accepted"

            ! At this point, site is within the cutoff radius so we count it.

            ! Increment the counter of the number of sites in the neighbor
            !   list for site i.
            neighCount(i) = neighCount(i) + 1

            ! Record the index of the j site into the i neighbor list of
            !   sites. Similarly, record the specific replicated cell index k.
            neighSite(neighCount(i),i) = j
            neighCell(neighCount(i),i) = k
            
            ! Record the relative xyz coordinates of the j site in cell k into
            !   the list of the i neighbors.
            neighxyzCoords(:,neighCount(i),i) = shiftedPotPos(:) &
                  & - currentPositions(:,1)

            ! Record the distance between site i and site j in cell k.
            shiftedPotSiteSep = sqrt(shiftedPotSiteSepSqrd)
            neighDist(neighCount(i),i) = shiftedPotSiteSep

            ! Compute the value of the "switch function" to smoothly connect
            !   the region inside the cutoff with the region outside the
            !   cutoff.
            neighSwitchFn(neighCount(i),i) = 0.5_double &
                  & * (cos(pi * shiftedPotSiteSep/cutoffLoEn) + 1.0_double)

            ! For the next part, we will often refer to "Quantum Theory of
            !   Angular Momentum: Irreducible Tensors, Spherical Harmonics,
            !   Vector Coupling Coefficients, 3nj Symbols" by Varshalovich DA,
            !   Moskalev AN, and Khersonski VK.; Singapore; Teaneck, NJ, USA:
            !   World Scientific Pub; 1988.
            ! For short: QTAM.

            ! The goal now is to convert the relative xyz coordinates of this
            !   atom inside the ball defined by the cutoffLoEn into a set of
            !   angular coordinates theta_0, theta, and phi that sit on the
            !   surface of a three-sphere (i.e., the three-dimensional surface
            !   of a four-dimensional solid object enclosed by a constant 4D
            !   radius).
            ! The reason is that if we take the specific distribution of atoms
            !   inside the 3D ball and cast them into spherical coordinates on
            !   the surface of a three-sphere, then we can used the spherical
            !   harmonics for 4D space as a natural basis to expand the
            !   function that defines the distribution of atoms. Note that
            !   the spherical harmonics for 4D systems requires only three
            !   parameters in much the same way that the spherical harmonics
            !   used for 3D systems takes only two parameters.
            ! The practical upshot is that we can use *one* basis set of
            !   4D spherical harmonic functions instead of needing two basis
            !   sets: one of 3D spherical harmonics and one of 1D radial
            !   functions.
            ! The outline of the details is as follows (noting that this
            !   outline does not exactly match the step-by-step algorithmic
            !   process because many steps may be omitted or reordered for
            !   computational efficiency):
            ! (1) Consider the xyz and spherical polar coordinates.
            ! (2) Compute the theta_0, theta, phi coordinates where the
            !     theta and phi spherical polar coodinates define the
            !     direction *to* the atom from the current target atom *and*
            !     they also define an axis about which we could imagine
            !     rotating by an amount omega. (Compare Fig. 1.2 (page 4) and
            !     Fig. 1.6 (page 23) of QTAM.)
            ! (3) The spherical polar rho coordinate is transformed into the
            !     theta_0 (omega in QTAM) coordinate. That is, we abstracly
            !     map the radial distance *along* the direction defined by
            !     theta and phi onto a rotation *about* the same axis defined
            !     by theta and phi.
            ! (4) The maximum rotation is pi (or a bit less, see below). This
            !     is because omega (0-pi), theta (0-pi), phi (0-2pi) can
            !     describe *any* rotation uniquely. If we permitted omega to
            !     go from 0-2pi, then it would be possible to express each
            !     rotation two different ways by also selecting different
            !     theta and phi values. (Specifically: theta->(pi-theta) and
            !     phi->(pi+phi).) This constraint is expressed in the first
            !     paragraph of 1.4.2 of QTAM (page 23).
            ! Note that "3D" spherical harmonics take two variables (theta and
            !   phi) and are independent of the third dimension of radius.
            !   Similarly, the "4D" spherical harmonics take three variables
            !   (say theta_0, theta, and phi) and are independent of the third
            !   dimension of radius.
            ! Regarding creation of the angular coordinates let's skip over
            !   the creation of theta_0 and return to it later. Instead, we
            !   first discuss the easier variables of theta and phi. These
            !   coordinates are exactly the same as the coordinates used in
            !   3D spherical polar coordinates.
            ! Theta takes values from 0 (when the direction to the atom is
            !   aligned along the +z axis) to Pi (when the direction to the
            !   atom is aligned along the -z axis).
            ! Phi takes values from 0 (+x axis) to 2Pi (+x axis again) with
            !   increasing values of phi shifting counter-clockwise toward the
            !   +y axis. See Figure 1.2 in section 1.1.2 of QTAM, page 4.
            ! The determination of theta and phi defines a vector direction
            !   toward a specific atom from the current target atom. The
            !   distance to that atom is traditionally thought of as the
            !   radial distance rho. However, we will consider that distance
            !   as being the measure of an angle, theta_0, defined by:
            !   theta_0 = theta_0_max * r/R_cut. We can choose any value for
            !   theta_0_max, but there are certain more convenient choices.
            !   (See below where the computation is actually performed.)
            ! First, realize that theta_0 plays a role in the three-sphere
            !   much as theta plays in the two sphere. I.e., it selects among
            !   the possible "latitudes" that the embedded two-sphere
            !   (parameterized by theta and phi) will sit on.
            ! Conveniently, the mathematical form of the 4D hyperspherical
            !   harmonics are closely related to the elements of unitary
            !   transformation matrices that perform a rotation about an
            !   axis defined by theta and phi through an angle defined by
            !   omega = 2*theta_0. Specifically, this is expressed in section
            !   4.5.3 of QTAM (page 82).
            ! However, as mentioned in section 4.5.4 of QTAM (page 82)
            !   equation 14, the range of variables that define the 4D
            !   spherical harmonics is 0 <= theta <= pi; 0 <= phi <= 2pi;
            !   0 <= omega <= 2pi. That is different from the range of the
            !   same variables that define a rotation in the case of omega
            !   only. Repeating, for defining unique rotations we have omega
            !   between 0 and pi, but for the 4D spherical harmonics we have
            !   omega between 0 and 2pi.
            ! Rotations may be understood by studying section 1.4.2 and
            !   Figure 1.6 in QTAM (page 23). Then, (as expressed in section
            !   1.4.3), if we understand that rather than using a traditional
            !   3-vector it is possible to define an arbitrary point using
            !   Equation 8 in section 1.4.3 of QTAM (page 24), we can later
            !   see that a rotation matrix can be expressed as a 2x2 matrix:
            !   Equation 11. (I.e., we may equate a rotation and a point.)
            ! Then, in sub-sections (a), (b) and (c) of section 1.4.4 the
            !   different sets of variables that can be used to define a
            !   rotation are related to each other. That is useful because
            !   from Equation 33 in section 1.4.5 (page 28) we obtain a
            !   simple expression for the rotation operator in terms of
            !   omega, theta, and phi. Because the different sets of variables
            !   for defining a rotation are equivalent, we should be able to
            !   express a rotation in either omega (theta_0), theta, phi or
            !   in Euler angles (alpha, beta, gamma) to achieve the same
            !   rotation.
            ! Crucially, in equation 35, we see that the matrix elements of
            !   the rotation operator when expressed in terms of Euler angles
            !   (alpha, beta, gamma) specifically associated with eigenstates
            !   of angular momentum can be expressed in terms of so-called
            !   Wigner-D functions.
            ! The practical upshot of that last paragraphs is the logical
            !   connection between rotation matrices in terms omega, theta,
            !   phi and the very useful spherical harmonics (eigenstates of
            !   angular momentum). That is, we can probably find a way to
            !   relate omega, theta, phi rotations of spherical harmonics
            !   and Wigner-D functions. That is precisely what we find in
            !   Equation 3 of section 4.5.2 (page 81) of QTAM.
            ! There, we learn that eigen functions of arbitrary omega, theta
            !   phi rotations (i.e. the 4D spherical harmonics which serve
            !   as a basis for such rotations) can be expressed in terms of
            !   Wigner-D functions.
            ! Then, we just need the explicit form of the Wigner-D functions
            !   which can be found in Equation (1) of section 4.3 (page 76)
            !   of QTAM followed by the explicit form of the "d" functions
            !   which can be any one of equations 2-5 in the same section.
            !   Somewhat arbitrarily, we use equation 4.

            ! Compute the theta_0, theta, and phi coordinates of this atom
            !   with respect to the current target atom.

            ! Compute theta_0 from r: theta_0 = angleSqueeze * pi * r/R_cutoff
            !   theta_0 nominally goes from 0 to pi to define a unique
            !   rotation. However, we are mapping radial distance onto an
            !   angle and the "early" part from 0 to some modestly small angle
            !   and the "later" part from pi to pi-(some modestly small angle)
            !   are "compressed". A variation in radial position maps onto a
            !   variation in angle. Consider a variation in radial position
            !   equal to 0.1. If the variation occurs close to r=0 or r=R_cut
            !   then it will be mapped into a fairly small range of theta_0
            !   values. Fortunately, there are unlikely to be any other atoms
            !   near the central atom. So, we don't have to worry about the
            !   r=0 side of things. To avoid the compression problem with the
            !   r=R_cut side of things we will only map the radial values into
            !   a limited range. Specifically between 0 and angleSqueeze*pi.
            !   So, if r=0 then theta_0, but if r=R_cut then
            !   theta_0=angleSqueeze*pi.
            ! Note additionally, that we will only be using a subset of the
            !   possible input values for the 4D spherical harmonics. The 4D
            !   spherical harmonics accomodate 0<=omega<=2pi but we will only
            !   use 0<=theta_0<=angleSqueeze*pi.
            neighttpCoords(1,neighCount(i),i) = &
                  & angleSqueeze * pi * shiftedPotSiteSep / cutoffLoEn

            ! Compute theta: theta = arccos(z/r) where z can take on values
            !   between +R_cutoff and -R_cutoff. Hence, theta will naturally
            !   take on values between 0 and pi. 
            neighttpCoords(2,neighCount(i),i) = &
                  & acos(neighxyzCoords(3,neighCount(i),i) / shiftedPotSiteSep)

            ! Compute phi:
            neighttpCoords(3,neighCount(i),i) = atan2( &
                  & neighxyzCoords(2,neighCount(i),i), &
                  & neighxyzCoords(1,neighCount(i),i))
            ! Phi in range 0,2pi.
            if (neighttpCoords(3,neighCount(i),i) < 0.0_double) then
               neighttpCoords(3,neighCount(i),i) = &
                     & neighttpCoords(3,neighCount(i),i) + 2.0_double * pi
            endif
            if (neighttpCoords(3,neighCount(i),i) >= 2.0_double * pi) then
               neighttpCoords(3,neighCount(i),i) = &
                     & neighttpCoords(3,neighCount(i),i) - 2.0_double * pi
            endif

!write (6,*) neighCount(i), i
!write (6,*) neighxyzCoords(:,neighCount(i),i)
!write (6,*) neighttpCoords(:,neighCount(i),i)

            ! Correct any numerical errors.
            if (neighttpCoords(1,neighCount(i),i) > pi) then
               neighttpCoords(1,neighCount(i),i) = pi
            endif
            if (neighttpCoords(1,neighCount(i),i) < 0.0_double) then
               neighttpCoords(1,neighCount(i),i) = 0.0_double
            endif
            if (neighttpCoords(2,neighCount(i),i) > pi) then
               neighttpCoords(2,neighCount(i),i) = pi
            endif
            if (neighttpCoords(2,neighCount(i),i) < 0.0_double) then
               neighttpCoords(2,neighCount(i),i) = 0.0_double
            endif
            if (neighttpCoords(3,neighCount(i),i) > 2.0_double * pi) then
               neighttpCoords(3,neighCount(i),i) = 2.0_double * pi
            endif
            if (neighttpCoords(3,neighCount(i),i) < 0.0_double) then
               neighttpCoords(3,neighCount(i),i) = 0.0_double
            endif
!write (6,*) neighttpCoords(:,neighCount(i),i)

         enddo ! k numCellsReal
      enddo ! j atom2
   enddo ! i atom1

end subroutine createNeighborList


function u_densityCoeff (twoj, twom, twomp, siteIndex)

   ! Use necessary modules
   use O_Kinds
   use O_MathSubs

   ! Define function return variable.
   complex(kind=double) :: u_densityCoeff

   ! Define passed parameters.
   integer :: twoj, twom, twomp
   integer :: siteIndex


   ! Define local variables.
   integer :: i

   ! Initialize the u_densityCoeff
!   u_densityCoeff = hypersphericalHarmonic4D(twoj,twom,twomp, &
!         & (/0.0d0,0.0d0,0.0d0/))
   u_densityCoeff = cmplx(0.0_double, 0.0_double, double)
!write (6,fmt="(a,2e12.3)") "u = ", u_densityCoeff
!write (6,fmt="(5i)"), neighCount(siteIndex),twoj,twom,twomp,siteIndex

   do i = 1, neighCount(siteIndex)
!write (6,fmt="(a,i3,e14.5)") "neighSwitchFn = ", i, neighSwitchFn(i,siteIndex)
      u_densityCoeff = u_densityCoeff &
            & + neighSwitchFn(i,siteIndex) &
            & * hypersphericalHarmonic4D(twoj,twom,twomp,&
            & neighttpCoords(:,i,siteIndex))
!write (6,fmt="(i2,4e15.5)") i, neighSwitchFn(i,siteIndex), &
!      & neighttpCoords(:,i,siteIndex)
!write (6,fmt="(a,2e15.5)") "  ",u_densityCoeff
   enddo
!write (6,fmt="(a,2e15.5)") "u = ", u_densityCoeff

end function u_densityCoeff


subroutine computeBispectrumComponent

   ! Use necessary modules
   use O_Kinds
   use O_Input, only : twoj1, twoj2
   use O_PotSites, only : numPotSites
   use O_MathSubs, only : clebschGordan

   ! Make sure no variables are accidentally declared.
   implicit none

   ! Define local variables.
   integer :: a, b, c, d, i, j
   integer :: twom, twomp, twoj
   integer :: twom1Index, twom2Index
   integer :: twom1pIndex, twom2pIndex
   integer :: twojIndex
   integer :: twoMaxjInt, twoMinjInt, twoMaxjHOInt, twoMinjHOInt
   real (kind=double), allocatable, dimension (:,:,:) :: cgc
   complex (kind=double), allocatable, dimension (:,:,:) :: uPreComp
complex (kind=double) :: b_temp
integer :: termCount

   ! For a given twoj1 and twoj2 we will compute and sum together the bsComp
   !   for all possible j values in steps of 1 between (twoj1 + twoj2)/2 and
   !   abs(twoj1 - twoj2)/2.

   ! The form is B_j1,j2,j =
   !   SUM_m1,m1p=-j1->j1 (
   !   SUM_m2,m2p=-j2->j2 (
   !   SUM_m,mp=-j->j ( (u_m,mp^j)* C_j1,m1,j2,m2^jm C_j1,m1p,j2,m2p^jmp 
   !      u_m1,m1p^j1 u_m2,m2p^j2 )))
   ! Note that C_j1,m1,j2,m2^jm C_j1,m1p,j2,m2p^jmp is written as:
   !   H_j1,m1,m1p,j2,m2,m2p,j,m,mp in the reference below.
   ! Following: Thompson AP, Swiler LP, Trott CR, Foiles SM, Tucker GJ.,
   !   "Spectral neighbor analysis method for automated generation of
   !   quantum-accurate interatomic potentials.", Journal of Computational
   !   Physics, Vol.285, p.31630. (2015), Available from:
   !   https://www.sciencedirect.com/science/article/pii/S0021999114008353

   ! Note that, for a fixed specific environment of atoms indexed by i and
   !   at positions theta_0, theta, phi, the u values depend only on the input
   !   of one j value and two m values. Hence, we only need to compute the u
   !   values once for each possible set of inputs.
   ! The tricky issue is that the inputs for u can be either integers or half
   !   odd integers. If j1 and j2 are integers, then m1, m1p, m2, m2p, j, m, mp
   !   are all integers. However, if either j1 or j2 is a half odd integer,
   !   then that set (either j1, m1, m1p OR j2, m2, m2p) will be half odd
   !   integer while the other will be integer; AND the j1, m, mp set will be
   !   half odd integer. Alternatively, if j1 AND j2 are half odd integer, then
   !   the associated u values will be given half odd integer input, BUT the
   !   j, m, mp inputs will be integer valued.
   ! So, we should just plan on constructing arrays of integer and half odd
   !   integer solutions for u given the j1 and j2 inputs. If there are no
   !   half odd integer solutions, then the length of that array will be zero.

   ! Figure out the sizes of the integer and half-integer arrays.

   if ((modulo(twoj1,2) == 0) .and. (modulo(twoj2,2) == 0)) then
!write (6,*) "Choice 1"
      ! Both j1 and j2 are integer and thus j values are integer.
      twoMaxjInt = twoj1 + twoj2
      twoMinjInt = twoj1 - twoj2 ! Known: twoj1 >= twoj2 when reading input.
      twoMaxjHOInt = 0
      twoMinjHOInt = 2 ! Prevents computation of _any_ HO u terms.

   elseif ((modulo(twoj1,2) == 1) .and. (modulo(twoj2,2) == 1)) then
!write (6,*) "Choice 2"
      ! Both j1 and j2 are half odd integer and thus j values are integer.
      twoMaxjInt = twoj1 + twoj2
      twoMinjInt = twoj1 - twoj2 ! Known: twoj1 >= twoj2 when reading input.
      twoMaxjHOInt = twoj1 ! We guaranteed: twoj1 >= twoj2 when reading input.
      twoMinjHOInt = twoj2 ! We guaranteed: twoj1 >= twoj2 when reading input.

   elseif (modulo(twoj1,2) == 1) then
!write (6,*) "Choice 3"
      ! Only j1 is half odd integer and thus j values are half odd integer.
      !   This is because we already determined that only one of j1 and j2
      !   are half odd integer through the first "if" statement.
      twoMaxjInt = twoj2
      twoMinjInt = twoj2
      twoMaxjHOInt = twoj1 + twoj2
      twoMinjHOInt = twoj1 - twoj2 ! Known: twoj1 >= twoj2 when reading input.

   else
!write (6,*) "Choice 4"
      ! Only j2 is half odd integer and thus j values are half odd integer.
      !   This is through the process of elimination.
      twoMaxjInt = twoj1
      twoMinjInt = twoj1
      twoMaxjHOInt = twoj1 + twoj2
      twoMinjHOInt = twoj1 - twoj2 ! Known: twoj1 >= twoj2 when reading input.
   endif

!write (6,*) "twoMinjInt twoMaxjInt = ", twoMinjInt, twoMaxjInt
!write (6,*) "twoMinjHOInt twoMaxjHOInt = ", twoMinjHOInt, twoMaxjHOInt

   ! Allocate space to hold the cgc values and initialize to zero.
   allocate(cgc(twoj2+1, twoj1+1, twoj2+1))
   cgc(:,:,:) = cmplx(0.0d0,0.0d0,double)

   ! Go through all pairs of m1 m2. For each pair determine which of the
   !   possible j values the current m1+m2 sum is less than or equal to.
   !   Record the
   ! Create the "bottom" cgc j value.
   twom1Index = 0
   do a = -twoj1, twoj1, 2 ! 2m1 values.
      twom1Index = twom1Index + 1
      twom2Index = 0
      do b = -twoj2, twoj2, 2 ! 2m2 values.
         twom2Index = twom2Index + 1
         twom = a+b
         twoj = twoj1+twoj2
         twojIndex = 0
         do while ((abs(twom) <= twoj) .and. (abs(twoj1-twoj2) <= twoj)) !twoj
            twojIndex = twojIndex + 1
            ! Store the fact that this m1+m2 pair has a cgc for this j value.
            cgc(twojIndex,twom1Index,twom2Index) = &
                  & clebschGordan(twoj1, twoj2, twoj, a, b, twom)
!write (6,fmt="(6i3)") twoj1, twoj2, twoj, a, b, twom
!write (6,fmt="(a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a)") "cg = CG(j1=S(",twoj1, &
!      & ")/2, m1=S(",a,")/2, j2=S(",twoj2,")/2, m2=S(",b,")/2, j3=S(",&
!      &twoj,")/2, m3=S(",twom,")/2)"
!write (6,fmt="(a)") "print(f'{cg.doit()}, {cg.doit().evalf()}')"
!write (6,fmt="(a,i1,a,i1,a,i1,a,e12.3)") "cgc(",twojIndex,",",twom1Index,&
!      &",",twom2Index,")=",cgc(twojIndex,twom1Index,twom2Index)
            twoj = twoj - 2
         enddo
      enddo
   enddo

   ! Allocate space to hold the pre-computed u density coefficients.
   allocate(uPreComp(0:max(twoMaxjInt,twoMaxjHOInt),&
         & -max(twoMaxjInt,twoMaxjHOInt):max(twoMaxjInt,twoMaxjHOInt),&
         & -max(twoMaxjInt,twoMaxjHoInt):max(twoMaxjInt,twoMaxjHOInt)))
!   allocate(uPreComp(0:max(twoMaxjInt,twoMaxjHOInt),&
!         & 0:max(twoMaxjInt+1,twoMaxjHOInt),&
!         & 0:max(twoMaxjInt+1,twoMaxjHOInt)))

   ! Allocate space to hold the bispectrum component result for each allowed j
   !   for each atom and then initialize it to zero. Note that:
   !   ((2j_1 + 2j_2) - (2j_1 - 2j_2))/2 + 1 = 2j_2 + 1
   allocate(bsComp((twoj2 + 1), numPotSites))
   allocate(bsCompSum(numPotSites))
   bsComp(:,:) = cmplx(0.0_double,0.0_double,double)
   bsCompSum(:) = cmplx(0.0_double,0.0_double,double)

   ! Visit each potential site and compute its bispectrum component.
   do i = 1, numPotSites

      ! Pre-compute the integer based u values.
!write (6,*) "Starting integer based u values."
      do a = twoMinjInt, twoMaxjInt, 2
         do b = -a, a, 2
            do c = -a, a, 2
               uPreComp(a, b, c) = u_densityCoeff(a, b, c, i)
            enddo
         enddo
      enddo
!write (6,*) "Ending integer based u values."

      ! Pre-compute the half-odd integer based u values.
!write (6,*) "Starting half-odd integer based u values."
      do a = twoMinjHOInt, twoMaxjHOInt, 2
         do b = -a, a, 2
            do c = -a, a, 2
               uPreComp(a, b, c) = u_densityCoeff(a, b, c, i)
            enddo
         enddo
      enddo
!write (6,*) "Ending half-odd integer based u values."


      ! Accumulate the bispectrum component.
!!termCount = 0
!      twom2Index = 0
!      do a = -twoj2, twoj2, 2 ! m2
!         twom2Index = twom2Index + 1
!         twom1Index = 0
!         do b = -twoj1, twoj1, 2 ! m1
!            twom1Index = twom1Index + 1
!            twom1pIndex = 0
!            twom = a+b
!            do c = -twoj1, twoj1, 2 ! m1 prime
!               twom1pIndex = twom1pIndex + 1
!               twom2pIndex = 0
!               do d = -twoj2, twoj2, 2 ! m2 prime
!                  twom2pIndex = twom2pIndex + 1
!                  twomp = c+d
!                  twoj = twoj1+twoj2
!                  twojIndex = 0
!                  do while ((abs(twom) <= twoj) .and. (abs(twomp) <= twoj) &
!                        & .and. (abs(twoj1-twoj2) <= twoj))
!                     twojIndex = twojIndex + 1
!!termCount = termCount + 1
!                     b_temp = &
!                           & + conjg(uPreComp(twoj,twom,twomp)) &
!                           & * cgc(twojIndex,twom1Index,twom2Index) &
!                           & * cgc(twojIndex,twom1pIndex,twom2pIndex) &
!                           & * uPreComp(twoj1,b,c) &
!                           & * uPreComp(twoj2,a,d)
!                     bsComp(twojIndex,i) = bsComp(twojIndex,i) + b_temp
!!                     bsComp(twojIndex,i) = bsComp(twojIndex,i) &
!!                           & + conjg(uPreComp(twoj,twom,twomp)) &
!!                           & * cgc(twojIndex,twom1Index,twom2Index) &
!!                           & * cgc(twojIndex,twom1pIndex,twom2pIndex) &
!!                           & * uPreComp(twoj1,a,b) &
!!                           & * uPreComp(twoj2,c,d)
!!if (twojIndex == 1) then
!   !write (6,*) "twojIndex=",twojIndex
!   !write(6,fmt="(a,7i3)") "abcd twoj twom twomp",a,b,c,d,twoj,twom,twomp
!   !write(6,fmt="(a,9i3)") "count m1 m2 m m1p m2p mp twoj twojIdx ", termCount, b, a, a+b, c, d, c+d, twoj, twojIndex
!   !write(6,fmt="(a,2e15.6)") "b_temp = ", b_temp
!   !write(6,fmt="(a,2e15.6)") "b_Accum(twojIndex,i)", bsComp(twojIndex,i)
!   !write(6,fmt="(a,2e15.6)") "conjg(uPreComp(twoj,twom,twomp))",&
!   !   & conjg(uPreComp(twoj,twom,twomp))
!   !write(6,fmt="(a,2e15.6)") "cgc(twojIndex,twom1Index,twom2Index)",&
!   !   & cgc(twojIndex,twom1Index,twom2Index)
!   !write(6,fmt="(a,2e15.6)") "cgc(twojIndex,twom1pIndex,twom2pIndex)",&
!   !   & cgc(twojIndex,twom1pIndex,twom2pIndex)
!   !write(6,fmt="(a,2e15.6)") "uPreComp(twoj1,a,b)",uPreComp(twoj1,a,c)
!   !write(6,fmt="(a,2e15.6)") "uPreComp(twoj2,c,d)",uPreComp(twoj2,b,d)
!!endif
!                     twoj = twoj - 2
!                  enddo ! while, j
!               enddo ! d, m2 prime
!            enddo ! c, m2
!         enddo ! b, m1 prime
!      enddo ! a, m1
!!write (6,*) "Term count = ", termCount
      ! Accumulate the bispectrum component.
      twom1Index = 0
      do a = -twoj1, twoj1, 2 ! m1
         twom1Index = twom1Index + 1
         twom1pIndex = 0
         do b = -twoj1, twoj1, 2 ! m1 prime
            twom1pIndex = twom1pIndex + 1
            twom2Index = 0
            do c = -twoj2, twoj2, 2 ! m2
               twom2Index = twom2Index + 1
               twom2pIndex = 0
               twom = a+c
               do d = -twoj2, twoj2, 2 ! m2 prime
                  twom2pIndex = twom2pIndex + 1
                  twomp = b+d
                  twoj = twoj1+twoj2
                  twojIndex = 0
                  do while ((abs(twom) <= twoj) .and. (abs(twomp) <= twoj) &
                        & .and. (abs(twoj1-twoj2) <= twoj))
                     twojIndex = twojIndex + 1
                     bsComp(twojIndex,i) = bsComp(twojIndex,i) &
                           & + conjg(uPreComp(twoj,twom,twomp)) &
                           & * cgc(twojIndex,twom1Index,twom2Index) &
                           & * cgc(twojIndex,twom1pIndex,twom2pIndex) &
                           & * uPreComp(twoj1,a,b) &
                           & * uPreComp(twoj2,c,d)
                     twoj = twoj - 2
                  enddo ! while, j
               enddo ! d, m2 prime
            enddo ! c, m2
         enddo ! b, m1 prime
      enddo ! a, m1
   enddo ! i, potential site index

   ! Print out the results in a human readable format to the std. output file.
   do i = 1, numPotSites
      write (20,fmt="(a,i10)") "Potential Site: ", i
      bsCompSum = cmplx(0.0_double,0.0_double,double)
      do j = 1, twoj2+1
         bsCompSum = bsCompSum + bsComp(j,i)
         write(20,fmt="(a,i4,a,e25.16)") "Bispectrum Component for 2j = ", &
            & (twoj1+twoj2) - (j-1)*2, " is: ", real(bsComp(j,i),double)
      enddo

      write (20,fmt="(a,e25.16)") "Accumulated Bispectrum Component is: ", &
            & real(bsCompSum(i),double)
   enddo



   ! Open the output file for machine readable results.
   open (unit=21,file='fort.21',status='new',form='formatted')

   ! Create the header for the columns of data.
   write (21,fmt="(a)") "        site          2j         BiSpecComp"//&
         & "              Total_BSC"

   ! Print out the results in a machine readable format.
   do i = 1, numPotSites
      do j = 1, twoj2+1
         write (21,fmt="(2i12,2e25.16)") i, (twoj1+twoj2) - (j-1)*2, &
               & real(bsComp(j,i),double), real(bsCompSum(i),double)
      enddo
   enddo

   ! Close the output file.
   close(21)

end subroutine computeBispectrumComponent

end module O_LOEN
