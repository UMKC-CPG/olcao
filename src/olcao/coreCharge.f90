module O_CoreCharge

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Core charge density cast into gaussians used for electronic potential.
   real (kind=double), allocatable, dimension (:)   :: coreChargeDensity

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

! Compute the core charge fitted with functions centered on each atom with
!   the constraint of having exact charge.  Note that the potential functions
!   are the same ones used here to fit the core charge.
subroutine makeCoreRho

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: pi, maxOrbitals
   use O_AtomicTypes, only: atomTypes, maxNumAtomAlphas, maxNumCoreRadialFns
   use O_AtomicSites, only: coreDim
   use O_PotTypes, only: potTypes, maxNumPotAlphas
   use O_PotSites, only: potSites, numPotSites
   use O_Potential, only: potDim
   use O_TimeStamps

   ! Import the necessary HDF modules
   use HDF5
   use O_SetupElecStatHDF5, only: coreChargeDensity_did, pot

   ! Import the necessary external interfaces
   use O_LAPACKDPOSVX

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the local variables used in this subroutine.
   integer :: i,j,k,l,m ! Loop variables.  loop1=i, nestedloop2=j ...
   integer :: hdferr
   integer :: info ! LAPACK error

   ! Iteration dependent variables.  For each potential site iteration these
   !   values will change when the potential type of the site changes.
   integer :: currType
   integer :: currNumPotAlphas
   integer :: currNumAtomAlphas
   integer :: currNumCoreRadialFns
   integer :: currAngMomIndex
   integer, allocatable, dimension (:) :: currNumOrbAlphas
   integer, allocatable, dimension (:) :: currCoreQN_lList
   real (kind=double), allocatable, dimension (:) :: currPotAlphas
   real (kind=double), allocatable, dimension (:) :: currAtomAlphas
   real (kind=double), allocatable, dimension (:) :: coreIntegrals ! This is
      ! used to store integrals of each Gaussian (i).  The index is over i.
      ! This is the I' vector described below.
   real (kind=double), allocatable, dimension (:) :: coreCoeffs ! This is
      ! used to store the initially determined coefficients that fit the
      ! Gaussian functions to the core atomic orbitals.  The index is over i.
      ! This is the G vector described below.
   real (kind=double), allocatable, dimension (:) :: coreCoeffsCorrected ! This
      ! is used to store the coefficients after they have been corrected to
      ! produce the exact core charge.
   real (kind=double), allocatable, dimension (:,:) :: coreBasisOverlap ! This
      ! is used to store integrals representing the overlap between all core
      ! atomic orbitals and each Gaussian (i).  The index is over i.  This is
      ! the I vector described below.  (Upon calling DPOSVX the contents of
      ! this vector will be overwritten with the solution.  This solution is
      ! then copied to the coreCoeffs for easier understanding.)  NOTE
      ! that the second index is always 1.  This is because the DPOSVX
      ! interface requires a multidimensional array.
   real (kind=double), allocatable, dimension (:,:) :: gaussSelfOverlap ! This
      ! stores the overlap between pairs of Gaussian functions used to
      ! represent the core charge.
   real (kind=double) :: maxAlphaRange
   real (kind=double) :: delta ! Lagrange multiplier
   real (kind=double) :: currNumCoreElectrons ! Current core charge exactly.
   real (kind=double) :: currFittedCharge ! Core charge using initial fit. 
   real (kind=double), dimension (4) :: factor ! Integration constants for
      ! integrals of two Gaussians from a basis function with their associated
      ! spherical harmonics.

   ! Local intermediate matrix used to speed calculation of other matrices
   !   and arrays.  (Essentially the above factor multiplied by the term from
   !   the integration that involves the exponential alphas.)
   real (kind=double), allocatable, dimension (:,:,:) :: alphaFactor

   ! Misc. factors that are derived from values in the constants module.
   real (kind=double) :: piXX32

   ! Initialize the pi factor
   piXX32  = pi**1.5_double

   ! Notes on integration:
   ! Regular integration of spherical Gaussian functions (note that the
   !   r^2*sin(theta) is just necessary for integration in spherical
   !   coordinates and it does not represent any angular character to the
   !   function):  int(r^0*exp(-ar^2)*r^2*sin(theta)) = pi^(3/2) / a^(3/2)

   ! Integration of two Gaussian functions from a particular basis function:
   !   The two Gaussians are characterized by a1 and a2 for the coefficients
   !   in the exponent and a factor of r^l where l is the QN_l value for the
   !   basis function that these Gaussians help to form.  This will look like:
   !   int(r^2*exp(-(a1+a2)*r^2) * r^2*sin(theta)) for two p-type Gaussians.

   ! This is not the whole story though because each basis function is also
   !   multiplied by a real spherical harmonic function so the integration of
   !   the two Gaussians should also include a factor from the integration of
   !   the square of the real spherical harmoic (because each basis function
   !   has one).  This will create the following integral for px-type:
   !   int(r^2*exp(-(a1+a2)*r^2) * r^2*sin(theta) * 
   !   (sin(theta)*cos(phi)*sqrt(3/4Pi))^2).  The 1st term is the sum of two
   !   p-type Gaussians, the 2nd term is for integration in spherical coords.,
   !   and the third term is the square of the px=Y1-1 spherical harmonic.
   !   This will work out to be sqrt(Pi)*3/8 * 1/(a1+a2)^(5/2) for px.  This
   !   is also the same for py and pz because the different spherical harmonics
   !   for p-type will all integrate to the same value
   !   (int(Y11) = int(Y10) = int(Y1-1)).

   ! All the s, p, d, and f-type spherical harmonics will each integrate to one
   !   value regardless of the QN_m for the p, d, or f orbitals.  The following
   !   descriptive equations lump all alphas to one term, use a generic Ylm,
   !   and include the r^2*sin(theta) for integration in spherical coords.
   !   Note the exponent on the alpha terms.  The alpha terms will be factored
   !   in during the actual computational loop.
   !   int(r^0*exp(-ar^2)*Y0m^2*r^2*sin(theta)) = pi^(1/2) *   (1/ 4) / a^(3/2)
   !   int(r^2*exp(-ar^2)*Y1m^2*r^2*sin(theta)) = pi^(1/2) *   (3/ 8) / a^(5/2)
   !   int(r^4*exp(-ar^2)*Y2m^2*r^2*sin(theta)) = pi^(1/2) *  (15/16) / a^(7/2)
   !   int(r^6*exp(-ar^2)*Y3m^2*r^2*sin(theta)) = pi^(1/2) * (105/32) / a^(9/2)
   factor(1) = sqrt(pi) *   1.0_double /  4.0_double
   factor(2) = sqrt(pi) *   3.0_double /  8.0_double
   factor(3) = sqrt(pi) *  15.0_double / 16.0_double
   factor(4) = sqrt(pi) * 105.0_double / 32.0_double

   ! It is then possible to simply multiply each integration by a factor of
   !   2*l+1 to get the accumulated contribution from each QN_m orbital of a
   !   particular QN_l angular momentum.  (e.g. the dxy, dxz, dyz, d(x^2-y^2),
   !   and d(3z^2-r^2) cases are all the same so just multiply by 5 to get the
   !   accumulated integral of all 5.)
   factor(1) = factor(1) * 1.0_double
   factor(2) = factor(2) * 3.0_double
   factor(3) = factor(3) * 5.0_double
   factor(4) = factor(4) * 7.0_double

   ! Then, because these basis functions are normalized to 1 and we want to
   !   contemplate the charge in each orbital we need to multiply by a factor
   !   of two because each orbital has two electrons in it.
   factor(:) = factor(:) * 2.0_double

   ! Finally, another factor of 2 is applied for a not so obvious reason.
   !   So, we will explain it here.  When computing the square of the
   !   atomic orbital function |u_jA|^2 down some ways in this subroutine you
   !   will recall that the atomic orbital contains a summation of Gaussian
   !   functions (times Ylm) so u_jA = [C1*exp(-a1*r^2) + C2*exp(-a2*r^2) +
   !   ... + CN*exp(-aN*r^2)]*Ylm.  When squared you have |u_jA|^2 =
   !   [C1*C1*exp(-(a1+a1)*r^2) + 2*C1*C2*exp(-(a1+a2)*r^2) +
   !   2*C1*C3*exp(-(a1+a3)*r^2) + ... CN*CN*exp(-(aN+aN)*r^2)]*Ylm^2.  Note
   !   the factor of 2 on all the terms where the indices are different.  This
   !   comes from the same situation as (x+y)^2 = x^2 + xy + yx + y^2 =
   !   x^2 + 2xy + y^2.  We are avoiding the double counting and have included
   !   the factor of 2 way up here so that it does not have to be multiplied
   !   repeatedly.  A factor of 0.5 is multiplied onto all the terms with the
   !   same indices to effectively remove the factor of 2 given here.
   factor(:) = factor(:) * 2.0_double

!Old values.  These match the above numbers including both factors of 2, but
!  the above description is a lot more clean.  I just left this here as a
!  reminder to people on how *not* to document a program.  (i.e. This had
!  absolutely zero documentation to show why these numbers were picked.)  (I
!  had added documentation for sqrtPi4 to indicate that it was equal to the
!  square root of Pi divided by 4.  Previously it was just written as
!  0.443113462 without *any* indication as to where that came from.  Insane.)
!   factor(1) = 2.0_double    * sqrtPi4 * 2.0_double
!   factor(2) = 9.0_double    * sqrtPi4 * 2.0_double
!   factor(3) = 37.5_double   * sqrtPi4 * 2.0_double
!   factor(4) = 183.75_double * sqrtPi4 * 2.0_double


   ! Allocate space for the global array containing the core charge density.
   allocate (coreChargeDensity(potDim))

   ! Allocate space for the local arrays and matrices
   allocate (coreBasisOverlap       (maxNumPotAlphas,2))
   allocate (coreIntegrals          (maxNumPotAlphas))
   allocate (coreCoeffs             (maxNumPotAlphas))
   allocate (coreCoeffsCorrected    (maxNumPotAlphas))
   allocate (gaussSelfOverlap       (maxNumPotAlphas,maxNumPotAlphas))
   allocate (currPotAlphas          (maxNumPotAlphas))
   allocate (currAtomAlphas         (maxNumAtomAlphas))
   allocate (currNumOrbAlphas       (maxOrbitals))
   allocate (currCoreQN_lList       (maxNumCoreRadialFns))
   allocate (alphaFactor (maxNumAtomAlphas,maxNumAtomAlphas,maxOrbitals))

   ! Start making the intermediate Core Charge Density
   call timeStampStart(13)

   ! If there are no core wave functions then we can skip this whole thing.
   if (coreDim == 0) return

   ! Initialize the core basis overlap, integrals, initial core coefficients,
   !   and corrected core coefficients.
   coreBasisOverlap    (:,:) = 0.0_double
   coreIntegrals       (:)   = 0.0_double
   coreCoeffs          (:)   = 0.0_double
   coreCoeffsCorrected (:)   = 0.0_double

   ! Initialize the coreChargeDensity array.
   coreChargeDensity (:) = 0.0_double

   ! Initialize the gaussSelfOverlap.
   gaussSelfOverlap (:,:) = 0.0_double

! OKAY, at present you can just ignore all of the following comments (in this
!   section).  This does not appear to be the way that things should work, but
!   I don't yet understand why.....UGH! The evidence is that the total energy
!   calculations are off by a bit now and this messes with the alignment of
!   XANES/ELNES spectra from different types (of the same element) in a given
!   system (e.g. alpha-Boron) such that the combination of the two spectra do
!   not line up to match experiment as well. Therefore, the code now has gone
!   back to the old implementation of the Lagrange Multiplier method.  I will
!   detail the math when I have time.  ^_^  (It isn't 100% wrong at all, but
!   there is some subtle detail that needs to be fleshed out. Proceed with some
!   care.)

   ! Our overall goal is to compute the charge density distribution of the core
   !   orbitals of the atomic basis functions and then use a sum of Gaussian
   !   functions to attempt to describe the same radial charge distribution.
   !   We can then say that the core charge distribution can be given by two
   !   different equations.  (One based on the atomic basis functions and the
   !   other based on a linear combination of Gaussians.)
   ! Rho_cA = 2 * SUM_orb(|u_orbA|^2)   That is, the core charge density
   !   (Rho_c) of an atom (A) is two times the sum of the squares of
   !   the orbital wave functions.  (Recall that these basis functions are all
   !   real at the present moment.)  The factor of 2 here is due to the fact
   !   that each orbital is used to represent spin degenerate states (i.e.
   !   there are 2 electrons per orbital).  Also, note here that the condition
   !   of normalization must hold for the orbital wave functions so that the
   !   sum of their contributions equals the total number of electrons in the
   !   core.  (This should have been made true during their construction.)
   !~Rho_cA = SUM_termsA(G_iA * exp(-b_iA * r^2))   That is, the core charge
   !   density when cast into Gaussians (~Rho_c) of an atom (A) is a sum of a
   !   set of Gaussian functions centered on the atomic site.  The number of
   !   Gaussian terms (i) depends on the atom being represented.  The values of
   !   the coefficients (G_iA) are unknown at this point and need to be solved
   !   for. The values of the coefficients (b_iA) are known and depend on the
   !   atom being represented (b=greek beta).  We also require that the
   !   normalization condition to hold for this expression so that the
   !   summation will equal the total number of electrons in the core.
   ! We set Rho_cA = ~Rho_cA so that we can solve for the G_iA coefficients.
   !   If we take this the naive approach we quickly see a problem:
   !
   !   2*SUM_orbA(|u_orbA|^2) = SUM_iA(G_iA*exp(-b_iA*r^2))
   !
   !   Integrate both sides over r.  (u_orbA is a function of r of course).
   !
   !   2*SUM_orbA(int(|u_orbA|^2)) = SUM_iA(G_iA*int(exp(-b_iA*r^2)))
   !
   !   2[int(|u_1A|^2) + int(|u_2A|^2) + ... + int(|u_NorbA|^2)] = 
   !      G_1A*int(exp(-b_1A*r^2)) + G_2A*int(exp(-b_2A*r^2)) + ... +
   !      G_NtermsA*int(exp(-b_NtermsA*r^2))
   !
   !   The left hand side simply evaluates to a single number (say X) while the
   !   right hand side has Nterms variables in it.  Clearly, one equation and
   !   Nterms unknowns is a problem.  To resolve this we will instead first
   !   multiply Rho_cA and ~Rho_cA by exp(-b_jA*r^2):
   !
   !   2*[SUM_orbA(|u_orbA|^2)*exp(-b_jA*r^2)] =
   !      SUM_iA(G_iA*exp(-b_iA*r^2)*exp(-b_jA*r^2)).
   !
   !   Then, we integrate and get an equation for each index 'j'.
   !
   !   2*SUM_orbA[int((|u_orbA|^2)*exp(-b_jA*r^2))] =
   !      SUM_iA[G_iA * int(exp(-b_iA*r^2)*exp(-b_jA*r^2))]
   !
   !   The left and the right sides can still be evaluated analytically and
   !   this will produce a set of Nterms equations with Nterms unknowns.  One
   !   for each b_jA.  This set of linear equations can be easily solved via
   !   the LAPACK routine dposvx so that all G_iA are known within the error of
   !   the LAPACK's fitting routine (which is ~1e-16).
   ! The meaning of the G_iA is SUM_termsA(G_iA*int(exp(-b_iA*r^2))) = Num
   !   core electrons.  The G_iA are coefficients for each Gaussian term used
   !   to describe the core charge density.  Well, almost.  The dposvx routine
   !   will not be able to produce coefficients G_iA that make the sum given
   !   here = Num core electrons because it is an unconstrained fitting.  There
   !   will be some error.  But, this is unacceptable because we need *exact*
   !   charge in the Gaussian description to make the long range Coulomb sums
   !   that appear later in the program converge.  If we do not have exact
   !   charge the sums will not converge.  Bad.
   ! So, let us define some things to make discussion easier.
   !   . = dot product
   !   N_c = Number of core electrons.
   !   N_t = Number of terms (indices of j).
   !   l = lambda = Lagrangian multiplier.
   !   vector_j is a column vector over the index j.
   !   matrix_ij is a matrix over rows i and columns j.
   !   I  = vector_j(int(|u_orb|^2 * exp(-b_jA*r^2)))   [Vector of integrals]
   !        Overlap between all core atomic orbitals and each j Gaussian fn.
   !        This talks about how much charge the jth Gaussian could represent.
   !
   !   I' = vector_i(int(exp(-b_jA*r^2)))               [Vector of integrals]
   !        Total volume of each Gaussian function.
   !
   !   S  = matrix_ij(int(exp(-(b_iA+b_jA)*r^2)))       [Matrix of integrals]
   !        Overlap between Gaussian function pairs.
   !       
   !   G  = vector_j(Coefficients of I' from LAPACK to give N_c (almost)).
   !  ~G  = vector_j(Corrected coefficients of I' to give N_c (exactly)).
   !
   ! Describing the earlier ideas in a linear algebra way we have:
   !   SG  = I
   !    G  = S^(-1)I   [This is solved for the G coefficients.]
   !
   ! We want to have G.I' - N_c = 0, but instead we get:
   !    G.I' -N_c = d /= 0
   !
   ! So we need to have a corrected G (~G) instead, then:
   !   ~G.I' - N_c = 0   (exactly).
   !
   ! We will use the method of Lagrange multipliers to solve this problem.  We
   !   want to extremize the function f(I') = ~G.I' - G.I' subject to the
   !   constraint h(I') = ~G.I'-N_c = 0.  The objective function is f and the
   !   constraint function is h so our Lagrangian is F(I') = f(I') + l*h(I') =
   !   ~G.I' - G.I' + l(~G.I' - N_c).
   !   The l = lambda = our lagrangian multiplier.
   !
   ! Then we set the gradient of F equal to zero and obtain a set of equations,
   !   one for each element in I' plus one for l.   (p=partial derivative)
   !   pF/pI'1 = partial of F wrt I' index 1 = ~G1 - G1 + l*~G1 = 0
   !   pF/pI'2 = partial of F wrt I' index 2 = ~G2 - G2 + l*~G2 = 0
   !   pF/pI'3 = partial of F wrt I' index 3 = ~G3 - G3 + l*~G3 = 0
   !   .
   !   .
   !   .
   !   pF/pI'N_t = partial of F wrt I' index N_t = ~GN_t - GN_t + l*~GN_t = 0
   !   pF/pl = partial wrt l = ~G1*I'1 + ~G2*I'2 + ... + ~GN_t*I'N_t - N_c = 0
   !
   ! We can write this set of equations in vector notation as:
   !   ~G - G + l*~G = 0  and  ~G.I'=N_c
   !
   ! Rewrite the first equation in terms of ~G:
   !   ~G + l*~G = G
   !   ~G(1+l)   = G
   !   ~G        = G/(1+l)
   !
   ! Substitute this into the second equation to solve for lambda (l).
   !   [G/(1+l)].I'   = N_c
   !   G.I'           = N_c * (1+l)
   !   G.I' / N_c     = 1+l
   !   G.I' / N_c - 1 = l
   !
   ! Substitute lambda (l) back in to the first equation to find the corrected
   !   coefficients that will give exact charge.
   !   ~G = G / (1+l)
   !   ~G = G / (1+[G.I' / N_c - 1])
   !   ~G = G / (G.I' / N_c)
   !   ~G = G * N_c / G.I'

   ! Clearly the corrected coefficients that minimize the objective function
   !   under the requirements of the constraint function are the old
   !   coefficients times the ratio of the real number of core electrons for
   !   the particular atom divided by the number of core electrons obtained
   !   from the original fitting.  (A simple scaling.)

   ! Initiate a loop over all the charge sites.
   do i = 1, numPotSites

      ! If this potential site is equivalent to another site (i.e. it is not
      !   the first atom of this potential type), then we can skip this site.
      if (potSites(i)%firstPotType == 0) cycle

      ! Note that it is assumed that the position of the atom at index i is the
      !   position as the potential at index i.

      ! Determine the type number for this atom/potential site.
      currType = potSites(i)%potTypeAssn

      ! If there are no core states for this particular atom type, then we can
      !   skip this potential site too.
      if (atomTypes(currType)%numCoreStates == 0) cycle

      ! Determine the number of potential alphas for this site.
      currNumPotAlphas = potTypes(currType)%numAlphas

      ! Determine the number of atomic alphas for this site.
      currNumAtomAlphas = atomTypes(currType)%numOrbAlphas(1)

      ! Determine the number of core wave functions for the current atom type.
      currNumCoreRadialFns = atomTypes(currType)%numCoreRadialFns

      ! Make a local copy of the current potential alphas for this type.
      currPotAlphas(:currNumPotAlphas) = &
            & potTypes(currType)%alphas(:currNumPotAlphas)

      ! Make a local copy of the current atomic alphas for this type.
      currAtomAlphas(:currNumAtomAlphas) = &
            & atomTypes(currType)%alphas(:currNumAtomAlphas)

      ! Make a local copy of the current number of alphas for each orbital type.
      currNumOrbAlphas(:) = atomTypes(currType)%numOrbAlphas(:)

      ! Make a local copy of the core orbital types for the current atom.
      currCoreQN_lList(:currNumCoreRadialFns) = &
            & atomTypes(currType)%coreQN_lList(:currNumCoreRadialFns)

      do j = 1, currNumPotAlphas

         ! Initialize the data structure for calculating the core charge
         !   atomic orbital basis overlap.  This is used to make (I).
         coreBasisOverlap(j,1) = 0.0_double
         coreBasisOverlap(j,2) = piXX32 / currPotAlphas(j)**1.5_double
         coreIntegrals(j)      = coreBasisOverlap(j,2)

         ! Fill out a two dimensional matrix of coefficients that will be
         !    used to detemine the core basis overlap below.  The efficiency
         !    of making the matrix ahead of time is that dimensions 1 (m) and
         !    2 (l) are the same for all orbital wave functions of the same
         !    orbital angular momentum.  So if there are 4 s core states, then
         !    this part of the calculation only has to be done once, not four
         !    times.  Note also that each m,l matrix is symmetric so we only
         !    have to calculate the upper triangle.
         do k = 1, atomTypes(currType)%maxCoreQN_l + 1
            do l = 1, atomTypes(currType)%numOrbAlphas(k)
               do m = l, atomTypes(currType)%numOrbAlphas(k)

                  ! Compute the inverse of the sum of the three alphas that
                  !   are from the different Gaussians (two atomic and one
                  !   potential).
                  maxAlphaRange = 1.0_double / (currPotAlphas(j) + &
                     & currAtomAlphas(l) + currAtomAlphas(m))

                  ! Almost finalize the integral of the product of the three
                  !   Gaussians.  We have yet to include the coefficients of
                  !   the atomic Gaussians.  This part here will account for
                  !   the alpha component and the other constants described
                  !   above.
                  alphaFactor(m,l,k) = factor(k) * sqrt(maxAlphaRange) * &
                        & maxAlphaRange**(k)
               enddo
            enddo
         enddo

         ! Using the above calculated values we can get the coreBasisOverlap.
         !   The basic idea is to square the atomic orbital functions.  These
         !   functions are actually a summation of Gaussian functions
         !   themselves and so we will get a summation of the product of each
         !   term with every other term. (a+b)*(a+b) = a^2+2ab+b^2.  The
         !   factor of 2 is included for all cross terms (e.g. 2ab) and the
         !   factor of 2 needs to be removed for the squared terms (e.g. a^2
         !   and b^2).  NOTE that a Gaussian function for representing the
         !   core charge has already been included in the alphaFactor above
         !   via the maxAlphaRange variable.
         do k = 1, currNumCoreRadialFns

            ! Get the array index number for the current orbital angular mom.
            currAngMomIndex = currCoreQN_lList(k) + 1

            ! Accumulate the cross product terms.
            do l = 1, atomTypes(currType)%numOrbAlphas(currAngMomIndex)
               do m = l, atomTypes(currType)%numOrbAlphas(currAngMomIndex)
                  coreBasisOverlap(j,1) = coreBasisOverlap(j,1) + &
                        & alphaFactor(m,l,currAngMomIndex) * &
                        & atomTypes(currType)%coreRadialFns(l,k) * &
                        & atomTypes(currType)%coreRadialFns(m,k)
               enddo

               ! Subtract off the factor of two that had been preincluded in
               !   "factor(k)" above.
               coreBasisOverlap(j,1) = coreBasisOverlap(j,1) - &
                     & alphaFactor(l,l,currAngMomIndex) * &
                     & atomTypes(currType)%coreRadialFns(l,k) * &
                     & atomTypes(currType)%coreRadialFns(l,k) * 0.5_double
            enddo
         enddo

         ! Calculate the upper triangle of the Gaussian function self overlap
         !   matrix.  Only the upper triangle is needed because the matrix is
         !   symmetric.
         do k = j, currNumPotAlphas
            gaussSelfOverlap(j,k) = piXX32 / (currPotAlphas(j) + &
                  & currPotAlphas(k))**1.5_double
            gaussSelfOverlap(k,j) = gaussSelfOverlap(j,k)
         enddo
      enddo


      ! Use the linear equation solver from lapack to get a least squares fit.
      !   This is a very precise fitting but it will not produce coefficients
      !   that will reproduce the total core charge of this site.
      call solveDPOSVX (currNumPotAlphas,2,&
            & gaussSelfOverlap(:currNumPotAlphas,:currNumPotAlphas),&
            & currNumPotAlphas,coreBasisOverlap(:currNumPotAlphas,:),info)

      if (info /= 0) then
         write (20, *) 'Core charge dposvx failed. INFO= ', info
         stop
      endif

      ! Copy the solution to a better named array.  (This is (G).)  Note that
      !   the above call to solveDPOSVX destroyed the prior contents of
      !   coreBasisOverlap.
      coreCoeffs(:currNumPotAlphas) = coreBasisOverlap(:currNumPotAlphas,1)


      ! Compute the total core charge from the fitted parameters.
      currFittedCharge = dot_product(coreCoeffs(:currNumPotAlphas),&
            & coreIntegrals(:currNumPotAlphas))


      ! Report the integrated core charge density as represented by Gaussian
      !   functions before normalization.
      write (20,*) "Core rho before correction for type ",currType, " = ", &
            & currFittedCharge

      ! Compute the number of core electrons for the atom at the current
      !   potential site (index i).
      currNumCoreElectrons = &
            & 2.0_double * real(atomTypes(currType)%numCoreStates,double)

      ! Compute the corrected coefficients.
      coreCoeffsCorrected(:currNumPotAlphas) = &
            & coreCoeffs(:currNumPotAlphas) * &
            & currNumCoreElectrons / currFittedCharge

      write (20,*) "Fitted core rho after correction1 for type ",currType, &
            & " = ", dot_product( &
            & coreCoeffsCorrected(:currNumPotAlphas), &
            & coreIntegrals(:currNumPotAlphas))

      ! Find the lagrange multiplier for exact charge
      delta = (currNumCoreElectrons - &
         & sum(coreBasisOverlap(:currNumPotAlphas,1) * &
         & coreIntegrals(:currNumPotAlphas))) / &
         & sum(coreBasisOverlap(:currNumPotAlphas,2) * &
         & coreIntegrals(:currNumPotAlphas))

      ! Save the exactly normalized charge density.
      coreChargeDensity&
         & (potTypes(currType)%cumulAlphaSum + 1: &
         &  potTypes(currType)%cumulAlphaSum + currNumPotAlphas) = &
         & coreBasisOverlap(:currNumPotAlphas,1) + delta * &
         & coreBasisOverlap(:currNumPotAlphas,2)
   enddo ! Potential site index i

   ! Deallocate the local data structures used in this subroutine
   deallocate (coreBasisOverlap)
   deallocate (coreIntegrals)
   deallocate (coreCoeffs)
   deallocate (coreCoeffsCorrected)
   deallocate (gaussSelfOverlap)
   deallocate (currPotAlphas)
   deallocate (currAtomAlphas)
   deallocate (currNumOrbAlphas)
   deallocate (currCoreQN_lList)
   deallocate (alphaFactor)

   ! Write the core charge density to disk for use later.
   call h5dwrite_f(coreChargeDensity_did,H5T_NATIVE_DOUBLE, &
         & coreChargeDensity(:),pot,hdferr)
   if (hdferr /= 0) stop 'Failed to write core charge density'

   ! Deallocate the core charge density matrix since it is no longer needed.
   deallocate (coreChargeDensity)

   ! Log the date and time we end.
   call timeStampEnd(13)

end subroutine makeCoreRho

end module O_CoreCharge
