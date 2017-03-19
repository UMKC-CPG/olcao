module O_PotentialUpdate

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define variables used to update the solid state potential. These arrays are
   !   retained across all SCF iterations.
   real (kind=double), allocatable, dimension (:) :: convergenceRecord
   real (kind=double), allocatable, dimension (:,:) :: totalEnergyRecord
   real (kind=double), allocatable, dimension (:,:,:) :: usedPotCoeffs ! Holds
         ! the currently used potential coefficients along with previous sets
         ! of used coefficients equal to the given feedbackLevel+1. The first
         ! index is potDim, the second is feedbackLevel+1 and the last is spin.
         ! Note that the +1 arises because we need to hold the *current* plus
         ! feedbackLevel number previous iterations.
   real (kind=double), allocatable, dimension (:,:,:) :: guessedPotCoeffs !
         ! Holds the currently guessed potential coefficients along with
         ! previous sets of guessed coefficients equal to the given
         ! feedbackLevel+1. The first index is potDim, the second is
         ! feedbackLevel+1 and the last is spin. Note that the +1 arises because
         ! we need to hold the *current* plus feedbackLevel number previous
         ! iterations.

   real (kind=double), allocatable, dimension (:,:)   :: tempPotCoeffs
   real (kind=double), allocatable, dimension (:,:,:) :: tempGuessedPotCoeffs
   real (kind=double), allocatable, dimension (:,:,:) :: tempUsedPotCoeffs


!   ! Define variables used to update the solid state potential.
!   real (kind=double), allocatable, dimension (:,:) :: xl0,xl1,xl2 ! Holds the
!         !   actually used potential coefficients from the current and previous
!         !   two iterations.
!   real (kind=double), allocatable, dimension (:,:) :: yl0,yl1,yl2 ! Holds the
!         !   guesses for the next set of potential coefficients from the
!         !   current and previous two iterations.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


subroutine makeSCFPot (totalEnergy)

   ! Import necessary modules.
   use HDF5
   use O_Kinds
   use O_LAPACKDPOSVX
   use O_TimeStamps
   use O_Input, only: numElectrons
   use O_Populate, only: occupiedEnergy
   use O_KPoints, only: numKPoints
   use O_Constants, only: smallThresh
   use O_AtomicTypes, only: atomTypes
   use O_AtomicSites, only: coreDim
   use O_PotSites, only: numPotSites, potSites
   use O_PotTypes, only: numPotTypes, potTypes
   use O_ExchangeCorrelation, only: maxNumRayPoints, radialWeight, exchRhoOp, &
         & exchCorrOverlap, numRayPoints
   use O_ElectroStatics, only: nonLocalNeutQPot, localNeutQPot, localNucQPot, &
         & nonLocalNucQPot, nonLocalResidualQ, potAlphaOverlap
   use O_SetupExchCorrHDF5, only: numPoints_did, radialWeight_did, &
         & exchRhoOp_did, exchCorrOverlap_did, numPoints, points, potPoints
   use O_Potential, only: spin, potDim, intgConsts, spinSplitFactor, &
         & potAlphas, potCoeffs, currIteration, lastIteration, feedbackLevel, &
         & relaxFactor, xcCode, converged, convgTest, GGA, numPlusUJAtoms, &
         & plusUJAtomID, plusUJ
   use O_SetupElecStatHDF5, only: potAlphaOverlap_did, coreChargeDensity_did, &
         & nonLocalNeutQPot_did, localNeutQPot_did, localNucQPot_did, &
         & nonLocalNucQPot_did, nonLocalResidualQ_did, pot, potPot, potTypesPot
   use O_ValeCharge, only: potRho, chargeDensityTrace, nucPotTrace, &
         & kineticEnergyTrace
   use O_MainPotRhoHDF5, only: potCoeffs_did, totalRhoCoeffs_did, &
         & valeRhoCoeffs_did, spinDiffRhoCoeffs_did, alphas_did, terms

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed dummy variables.
   real (kind=double) :: totalEnergy

   ! Define the local variables used in this subroutine.
   integer :: i,j,k ! Loop index variables
   integer :: hdferr
   integer :: errorCount
   integer :: info
   integer :: numOpValues
   integer :: currentType
   integer :: currentNumAlphas
   integer :: currentCumulAlphaSum
   integer :: potTypeInitIndex
   integer :: potTypeFinIndex
   integer :: siteIndex
   integer :: potTermCount
   integer :: currNumAlphas
   integer :: maxFeedback
   integer :: totalEnergyImprovement
   real (kind=double) :: fittedCharge
   real (kind=double) :: fittedSpinDiffCharge
   real (kind=double) :: spinDiffDifference
   real (kind=double) :: electronDiff
   real (kind=double) :: currentNucCharge
   real (kind=double) :: sumIntegratedCharge
   real (kind=double) :: kineticEnergy
   real (kind=double) :: exchCorrEnergy
   real (kind=double) :: exchCorrEnergyDiff
   real (kind=double) :: elecStatEnergy
   real (kind=double) :: elecStatEnergyDiff
   real (kind=double) :: coreSum
   real (kind=double) :: totalSum
   real (kind=double) :: spinDiffSum
   real (kind=double) :: SX 
   real (kind=double) :: SY
   real (kind=double) :: SZ
   real (kind=double) :: SXX
   real (kind=double) :: SXY
   real (kind=double) :: SXZ
   real (kind=double) :: SYY 
   real (kind=double) :: SYZ
   real (kind=double) :: SZZ
   real (kind=double) :: SXC 
   real (kind=double) :: SYC
   real (kind=double) :: SZC
   real (kind=double) :: SXXC
   real (kind=double) :: SXYC
   real (kind=double) :: SXZC
   real (kind=double) :: SYYC
   real (kind=double) :: SYZC
   real (kind=double) :: SZZC
   real (kind=double) :: SXS
   real (kind=double) :: SYS
   real (kind=double) :: SZS

   real (kind=double) :: blendingFactor
   real (kind=double) :: testableDelta
   real (kind=double) :: radialWeightSum
   real (kind=double) :: weightedPotDiff
   real (kind=double) :: totalMagneticMoment
   real (kind=double), dimension (4) :: currentExchCorrPot
   real (kind=double), dimension (2) :: totalEnergySpin
   real (kind=double), dimension (2) :: kineticEnergySpin
   real (kind=double), dimension (2) :: exchCorrEnergySpin
   real (kind=double), allocatable, dimension (:,:) :: generalRho
   real (kind=double), allocatable, dimension (:,:) :: tempOverlap
   real (kind=double), allocatable, dimension (:,:) :: elecStatPot
   real (kind=double), allocatable, dimension (:,:) :: exchCorrPot
   real (kind=double), allocatable, dimension (:,:) :: exchCorrRho
   real (kind=double), allocatable, dimension (:,:) :: exchCorrRhoCore
   real (kind=double), allocatable, dimension (:,:) :: exchCorrRhoSpin
   real (kind=double), allocatable, dimension (:,:) :: outputPot
   real (kind=double), allocatable, dimension (:,:) :: realSpacePotDiff
   real (kind=double), allocatable, dimension (:)   :: averageDelta
   real (kind=double), allocatable, dimension (:)   :: maxDelta
   real (kind=double), allocatable, dimension (:)   :: typesMagneticMoment
   real (kind=double), allocatable, dimension (:,:) :: potDifference ! Holds
         ! the difference between the currently used potential coefficients and
         ! the current guess for the next set of potential coefficients.

   real (kind=double) :: th1
   real (kind=double) :: th2
   real (kind=double) :: t1
   real (kind=double) :: t2
   real (kind=double) :: s11
   real (kind=double) :: s12
   real (kind=double) :: s22
   real (kind=double) :: det
   real (kind=double), allocatable, dimension (:) :: drl1,drl2
   real (kind=double), allocatable, dimension (:) :: rl0,rl1,rl2

   ! Log the date and time we start.
   call timeStampStart (18)

   ! Allocate space for a temporary array that will hold various incarnations
   !   of the various overlap matrices.
   allocate (tempOverlap    (potDim,potDim))

   ! Allocate space for local matrix that store various charge density values.
   allocate (generalRho     (potDim,6+spin))

   ! Allocate space for the magnetic moments for each potential type.
   allocate (typesMagneticMoment (numPotTypes)) 

   ! Copy the electrostatic potential constants of integration and the 
   !   charge density potential (potRho) into a data structure
   generalRho(:,1)     = potRho(:,1)  ! (Spin up + Spin down) or total
   generalRho(:,2)     = intgConsts(:)
   if (spin == 2) then
      generalRho(:,3)  = potRho(:,2)  !  Spin up - Spin down
   endif

!write (20,*) "potRho 1 iter",currIteration
!write (20,*) potRho(:,1)
!write (20,*) "potRho 2 iter",currIteration
!write (20,*) potRho(:,2)
!write (20,*) "potRhoMagDiff"
!do i = 1, potDim
!   write (20,*) abs(potRho(i,1))-abs(potRho(i,2))
!enddo
!call flush(20)

   ! Deallocate the potRho since it will not be needed until the
   !   next scf iteration calls makeValenceRho.
   deallocate (potRho)

   ! Allocate space for the electrostatic potential overlap matrix.
   allocate (potAlphaOverlap (potDim,potDim))

   ! Read the electrostatic potential overlap matrix
   call h5dread_f (potAlphaOverlap_did,H5T_NATIVE_DOUBLE,&
         & potAlphaOverlap(:,:),potPot,hdferr)
   if (hdferr /= 0) stop 'Failed to read pot alpha overlap'

   ! Copy the upper triangle of the potAlphaOverlap matrix to prevent its
   !   destruction.
   do i = 1, potDim
      tempOverlap(1:i,i) = potAlphaOverlap(1:i,i)
   enddo

   ! Solve the linear set of equations with LAPACK to get the coefficients for
   !   total valence charge and (for spin polarized calcs.) the spin difference
   !   of the valence charge.
   call solveDPOSVX (potDim,spin+1,tempOverlap,potDim,generalRho(:,:spin+1),&
         & info)

   if (info /= 0) then
      write (20, *) 'Fitting dposvx failed. INFO= ', info
      stop
   endif

   ! Now we need to make sure that the spin difference is scaled in the same
   !   way as the total charge was.
   if (spin == 2) then
      ! Compute the fitted valence charge difference.  (Up - Down)
      fittedSpinDiffCharge = dot_product(intgConsts(:potDim),&
            & generalRho(:potDim,3))

      ! Compare the fitted charge difference to the charge difference computed
      !   from the charge density matrix (in valeCharge.f90). We will treat the
      !   earlier charge spin difference as the "real one".
      spinDiffDifference = chargeDensityTrace(2) - fittedSpinDiffCharge

      ! Preserve the spin (UP-DOWN) charge density coefficients in 
      !   generalRho(:,8).
      generalRho(:potDim,8) = generalRho(:potDim,3)

      ! Scale the spin difference to match that computed by the charge density
      !   trace.
      generalRho(:potDim,8) = generalRho(:potDim,8) + generalRho(:potDim,2) * &
            & spinDiffDifference/sum(intgConsts(:potDim)*generalRho(:potDim,2))

      ! Show that the newly scaled fitted spin difference is consistent with the
      !   charge density trace.
      write (20,fmt="(a30,e16.8,a8,e16.8)") "Spin Difference Matching Error",&
            & spinDiffDifference," Out Of ", chargeDensityTrace(2)

      ! Demonstrate that the constrained fit will produce an exact spin diff.
      write (20,fmt="(a33,e16.8)") "Constrained spin difference e- = ", &
            & dot_product(intgConsts(:potDim),generalRho(:potDim,8))
   endif
   
   ! Compute the fitted valence charge.  (Up + Down for spin polarized calcs.)
   fittedCharge = dot_product(intgConsts(:potDim),generalRho(:potDim,1))

   ! Calculate the difference between the electron number and the fitted charge
   !   determined above.  This represents the difficulty of the Gaussians to
   !   accurately represent the true charge density.
   electronDiff = numElectrons - fittedCharge

   ! Record the electron fitting error.
   write (20,fmt="(a37,e16.8,a8,i8)") &
         & 'Unconstrained Vale Electron Fit Error',&
         & electronDiff, ' Out Of ',numElectrons
   call flush (20)

   ! Correct and store the valence charge coefficients in generalRho(:,2).
   !   This would overwrite the spin difference for spin polarized
   !   calculations, but that was already saved in generalRho(:,8).
!generalRho(:potDim,2) = generalRho(:potDim,1)*numElectrons/fittedCharge
   generalRho(:potDim,2) = generalRho(:potDim,1) + generalRho(:potDim,2) * &
         & electronDiff / sum(intgConsts(:potDim)*generalRho(:potDim,2))

   ! Similarly scale and store the spin difference.

   ! Demonstrate that the constrained fit will produce an exact valence charge.
   write (20,*) "Constrained fit num e- = ",dot_product(intgConsts(:potDim), &
         & generalRho(:potDim,2))

   ! At this point we will apply a kick to each of the coefficients for the
   !   fitted charge for each potential type.  This is done only on the first
   !   iteration of spin polarized calculations.  On the first iteration, the
   !   spin difference is always equal to zero.  We will kick each coefficient
   !   by an amount proportional to the amount of charge represented by the
   !   Gaussian (the spinSplitFactor for each type).
   if ((currIteration == 1) .and. (spin == 2)) then
      do i = 1, numPotTypes
         do j = potTypes(i)%cumulAlphaSum+1, &
              & potTypes(i)%cumulAlphaSum+potTypes(i)%numAlphas
            generalRho(j,8) = generalRho(j,2) * spinSplitFactor(j)
         enddo
      enddo
   endif

   ! Read the core charge density into generalRho 3.
   if (coreDim /= 0) then
      call h5dread_f (coreChargeDensity_did,H5T_NATIVE_DOUBLE,&
            & generalRho(:,3),pot,hdferr)
      if (hdferr /= 0) stop 'Failed to read core charge density'
   else
      generalRho(:,3) = 0.0_double
   endif

   ! At present generalRho(:,2) holds the valence coefficients of the up+down
   !   charge for both the spin polarized and non spin polarized case.  The
   !   generalRho(:,3) holds the core charge coefficients.  We will now compute
   !   the total (up+down, core+valence) and store it in generalRho(:,1).
   generalRho(:,1) = generalRho(:,2) + generalRho(:,3)

   ! Now that all of the important charge coefficients have been computed, we
   !   can store them for this iteration.

   ! Store the potential coefficients that were actually used in this
   !   iteration.
   do i = 1, spin
      call h5dwrite_f (potCoeffs_did(currIteration,i),&
            & H5T_NATIVE_DOUBLE,potCoeffs(:,i),terms,hdferr)
      if (hdferr /= 0) then
         write (20,*) 'Cannot write total (or spin up) potential coefficients.'
         stop
      endif
   enddo

   ! Store the total (valence+core (up+down)) charge density coefficients.
   call h5dwrite_f (totalRhoCoeffs_did(currIteration),H5T_NATIVE_DOUBLE,&
         & generalRho(:,1),terms,hdferr)
   if (hdferr /= 0) stop 'Cannot write total Rho coefficients.'

   ! Store the valence (up+down) charge density coefficients.
   call h5dwrite_f (valeRhoCoeffs_did(currIteration),H5T_NATIVE_DOUBLE,&
         & generalRho(:,2),terms,hdferr)
   if (hdferr /= 0) stop 'Cannot write vale Rho coefficients.'

   ! Store the spin difference charge density coefficients if necessary.
   if (spin == 2) then
      call h5dwrite_f (spinDiffRhoCoeffs_did(currIteration),H5T_NATIVE_DOUBLE,&
            & generalRho(:,8),terms,hdferr)
      if (hdferr /= 0) stop 'Cannot write spin difference Rho coefficients.'
   endif

   ! Store the exponential alphas (that are shared by all the above coefficient
   !   sets) on the first iteration only.
   if (currIteration == 1) then
      call h5dwrite_f (alphas_did,H5T_NATIVE_DOUBLE,potAlphas(:),&
            & terms,hdferr)
      if (hdferr /= 0) stop 'Cannot write Gaussian exponential alpahs.'
   endif

   ! Divide the charge so far into neutral spherically symmetric pieces and a
   !   residual part.  The previous sentence is the original documentation, but
   !   what appears to actually be done is just sum together the core charge
   !   density (3) with the adjusted valence charge density (2) and store the
   !   result in (1).  In the case of spin polarized calculations, (2) is also
   !   the total valence charge density, and is considered to be the up+down.
   !   generalRho(:,1) = TOTAL      CHARGE (CORE + VALENCE) and (UP + DOWN)
   !   generalRho(:,2) = VALENCE    CHARGE
   !   generalRho(:,3) = CORE       CHARGE
   !   generalRho(:,4) = NEUTRAL    CHARGE
   !   generalRho(:,5) = RESIDUAL   CHARGE (Held in first index of each type.)
   !   generalRho(:,6) = RESIDUAL   CHARGE (Held in first numPotTypes indices.)
   !   generalRho(:,7) = NUCLEAR    CHARGE (Held in last index of each type.)
   !   generalRho(:,8) = DIFFERENCE CHARGE (SPIN DIFFERENCE or ABSENT)
   ! As a (possibly obvious) note: All charges are from e- except for NUCLEAR.
   generalRho(:,4) = generalRho(:,1)
   generalRho(:,5) = 0.0_double
   generalRho(:,6) = 0.0_double
   generalRho(:,7) = 0.0_double
!write (20,*) "generalRho1 iter",currIteration
!write (20,*) generalRho(:,1)
!write (20,*) "generalRho2 iter",currIteration
!write (20,*) generalRho(:,2)
!write (20,*) "generalRho3 iter",currIteration
!write (20,*) generalRho(:,3)
!write (20,*) "generalRho4 iter",currIteration
!write (20,*) generalRho(:,4)
!call flush(20)

   ! This part is not well documented.  But later it is mentioned that the
   !   goal is to simulate the nuclear charges by affecting the shortest
   !   range potential in a way by the nuclear charge.
   ! Compute the difference between the nuclear charge and the total (up+down,
   !   core+valence) charge that is at the site associated with this type.
   !   Then, modify the first coefficient of the total (u+d,c+v) charge so that
   !   the broadest Gaussian reflects that deviation from neutral.
   do i = 1,numPotTypes

      ! Initialize variables for the current potential type
      currentNumAlphas = potTypes(i)%numAlphas
      currentNucCharge = potTypes(i)%nucCharge * potTypes(i)%multiplicity
      currentCumulAlphaSum = potTypes(i)%cumulAlphaSum

      ! Identify the beginning and ending indices of the list of potential
      !   and fitted charge density terms for this potential type (from the
      !   full list of all potential terms for the whole system).
      potTypeInitIndex = currentCumulAlphaSum + 1
      potTypeFinIndex  = currentCumulAlphaSum + currentNumAlphas

      ! Compute the total fitted charge associated with the terms of this
      !   potential type. This does not include the multiplicity of all
      !   contributions from all sites of this type. This is just one "generic"
      !   site of this current type. Note also that this is pulled from
      !   generalRho(:,1) which is up+down and core+valence.
      sumIntegratedCharge = sum(intgConsts(potTypeInitIndex:potTypeFinIndex) * &
            & generalRho(potTypeInitIndex:potTypeFinIndex,1))

      ! Assign values to generalRho (:,4:7)

      ! Store the difference between the fitted charge for this type and the
      !   nuclear charge associated with sites of this type. Note that all
      !   other terms in generalRho(:,5) except the initial index number for
      !   each type are held at zero. By dividing by the integration constant
      !   for the associated term we recover a proper coefficient (like the
      !   other coefficients for the fitted charge density) that will be used
      !   subsequently.
      generalRho(potTypeInitIndex,5) = (sumIntegratedCharge-currentNucCharge) &
            & / intgConsts(potTypeInitIndex)

      ! Recall that generalRho(:,4) holds a copy of the total (up+down,
      !   core+valence) fitted charge density. Now, we modify the first
      !   coefficient of the current type (which is also the smallest
      !   coefficient of the current type and therefore also the broadest
      !   gaussian). The idea is that generalRho(:,4) will hold the total
      !   electronic charge at the site associated with this type with a slight
      !   perturbation that represents the degree to which the charge deviates
      !   from neutral thus producing a "neutral" charge. (Recall that there is
      !   no rigid connection implying that the fitted charge at this site
      !   "belongs" to the atom at the current site.) Another way to think of
      !   it is that the generalRho(:,4) charge starts out describing whatever
      !   charge happens to be at a site and then the difference between that
      !   and the neutral charge is subtracted off, producing a net neutral
      !   charge. 
      generalRho(potTypeInitIndex,4) = generalRho(potTypeInitIndex,1) - &
            & generalRho(potTypeInitIndex,5)

      ! This is a bit awkward, but it is somewhat useful later on. Here, we
      !   keep track of the same deviation values as computed above and stored
      !   in generalRho(:,5) [first index of each type] and instead put the
      !   quantities all together in a row starting with index #1 and going up
      !   to index #numPotTypes. (That is, for generalRho six in particular,
      !   the indices do not align with the other indices of the other
      !   generalRho values, the intgConsts, or the alphas.)
      generalRho(i,6) = generalRho(potTypeInitIndex,5)

      ! Compute a proper coefficient to describe the nuclear charge using the
      !   sharpest Gaussian of the fitted electronic charge for the current
      !   type. It may be worth noting that the Gaussian used here may not be
      !   (i.e. almost certainly isn't) the same as the Gaussian used to
      !   describe the nucleus for the electron-nucleus interaction integrals.
      !   Presently, the nuclear Guassian has an exponential coefficient that
      !   is conventionally set to 20.0 for interaction integral calculation
      !   while the sharpest Gaussian for the charge density is on the order of
      !   1e5 to 1e9. Thus, there may be some important difference that is not
      !   being accounted for.
      generalRho(potTypeFinIndex,7) = -currentNucCharge / &
            & intgConsts(potTypeFinIndex)
   enddo
!write (20,*) "generalRho5 iter",currIteration
!write (20,*) generalRho(:,5)
!write (20,*) "generalRho6 iter",currIteration
!write (20,*) generalRho(:,6)
!write (20,*) "generalRho7 iter",currIteration
!write (20,*) generalRho(:,7)
!write (20,*) "generalRho8 iter",currIteration
!write (20,*) generalRho(:,8)
!call flush(20)


   ! Allocate space for holding specific components of the electrostatic
   !   potential. These electrostatic potential components will be used to
   !   produce related electrostatic contributions to the total energy of the
   !   system. The expressions for the contributions to the total energy have
   !   the form: <rho|op|rho> where the <rho| and |rho> are some kinds of
   !   fitted charge representation and the op is a real, symmetric operator.
   !   When the operator operates on a given |rho> it will produce a potential
   !   that will then be "felt" by the <rho| on the left to produce an
   !   electrostatic energy. Naturally though, because the operator is
   !   symmetric, there is no reason to not think that the operator could act
   !   toward the left <rho| creating a potential that is "felt" by the right
   !   <rho| to produce the same energetic result.

   ! Note importantly that as usual, the quantities stored here are just
   !   coefficients for the same Gaussian functions that are used in other
   !   contexts (e.g. charge).

   ! The reason to make a particular choice is that we choose to gather only
   !   specific charge-charge interactions when constructing the effective
   !   potential. For example, the potential in the Hamiltonian will already
   !   include a nuclear-electron interaction that is computed directly and thus
   !   it is not necessary to include it in the effective potential. Further,
   !   by being selective with the choice of op|rho> we can reduce the need for
   !   repetitive computations. (I.e. we can re-use one op|rho> with different
   !   values for <rho|.)

   ! There are six sets of quantities that are filled as follows:
   !   elecStatPot(:,1) = Nonlocal neutralized total electronic sites.
   !   elecStatPot(:,2) = Local neutralized valence electronic sites.
   !   elecStatPot(:,3) = Local core electronic sites.
   !   elecStatPot(:,4) = Nonlocal nuclear sites.
   !   elecStatPot(:,5) = Local nuclear sites.
   !   elecStatPot(:,6) = Nonlocal residual electronic sites.

   ! Recall that "neutralized" means that the broadest Gaussian has been
   !   modified to subtract the difference between the neutral (isolated atom)
   !   charge at the potential site and the actual fitted charge at the
   !   potential site.

   ! Recall that each of these will be paired with one or more other charge
   !   descriptions to produce an electrostatic energy contribution. The
   !   nonlocal terms will describe contributions from charges on other sites
   !   affecting charges on a given site. The local terms will describe
   !   contributions from charges on the same site as charge on a given site.

   allocate (elecStatPot (potDim,6))


   ! Allocate space for the non-local real-space electrostatic operator matrix.
   allocate (nonLocalNeutQPot (potDim,potDim))


   ! Read the non-local real-space electrostatic operator.
   call h5dread_f (nonLocalNeutQPot_did,H5T_NATIVE_DOUBLE,&
         & nonLocalNeutQPot(:,:),potPot,hdferr)
   if (hdferr /= 0) stop 'Failed to read non local neut q pot'


   ! Operate on the spherically symmetric neutral charge. After this, the
   !   elecStatPot(:,1) will hold quantities that describe the electrostatic
   !   potential due to the electronic (Q) influence of all other sites (hence
   !   "nonLocal") with residual components that cause deviations from neutral
   !   removed (hence "NeutQ") that will ultimately be "felt" by other fitted
   !   charge terms of some kind <rho|. (Side note, recall that generalRho(:,4)
   !   holds the neutral total charge (core+valence and up+down).)
   do i = 1,potDim
      elecStatPot(i,1) = sum (nonLocalNeutQPot(:,i) * generalRho(:,4))
   enddo

   deallocate (nonLocalNeutQPot)


   ! Allocate space for the local neutral charge potential operator matrix.
   allocate (localNeutQPot (potDim,potDim))

   ! Read the local real-space electrostatic operator.
   call h5dread_f (localNeutQPot_did,H5T_NATIVE_DOUBLE,&
         & localNeutQPot(:,:),potPot,hdferr)
   if (hdferr /= 0) stop 'Failed to read local neut q pot'


   ! Operate on the (valence - residual) and core charges seperately. Recall
   !   that generalRho(:,5) is all zeros except for the first term of each
   !   type. That first term is used to represent the deviation from neutral
   !   charge at a given site. Also recall that generalRho(:,2) is the up+down
   !   fitted valence charge (which is not typically neutral with its nucleus
   !   at a given site). Finally, recall that generalRho(:,3) is the core
   !   fitted charge without any residual component removed, (obviously it is
   !   also up+down).
   ! When complete, elecStatPot(:,2) will hold terms describing the potential
   !   due to the electrostatic influence of neutral valence electronic charge
   !   *at the same site* (hence "local") as a given other charge term (to be
   !   identified later). Note that the difference (generalRho(:,2) -
   !   generalRho(:,5)) will produce a neutral fitted valence charge as the
   !   source of the electronic charge generating the potential.
   ! Also, when complete, elecStatPot(:,3) will hold terms describing the
   !   potential due to the electrostatic influence of core electronic charge
   !   *at the same site* (hence "local") as a given other charge term (to be
   !   identified later).
   do i = 1, potDim
      elecStatPot(i,2) = sum(localNeutQPot(:,i) * &
            & (generalRho(:,2) - generalRho(:,5)))
      elecStatPot(i,3) = sum(localNeutQPot(:,i) * generalRho(:,3))
   enddo

   deallocate (localNeutQPot)


   ! Read in the nonlocal nuclear integral vector (gaussian screened). This
   !   vector is basically just like the below vector except that it describes
   !   the potential at a site due to nuclei at other sites.
   call h5dread_f (nonLocalNucQPot_did,H5T_NATIVE_DOUBLE,&
         & elecStatPot(:,4),pot,hdferr)
   if (hdferr /= 0) stop 'Failed to read non local nuc q pot'


   ! Read in the local nuclear integral vector. This vector (elecStatPot(:,5))
   !   already holds the potential at a given site due to nuclear charges at
   !   the same site as the given type (hence "NucQ" and "local"). This
   !   is a simple "read" process because the nuclei are static in position
   !   (unless forces are introduced in the future) and constant in their
   !   magnitude (unlike the electronic charge distriubutions). Because they
   !   never change during the SCF cycle, they were computed once already in
   !   setup and do not need to be computed again.
   call h5dread_f (localNucQPot_did,H5T_NATIVE_DOUBLE,&
         & elecStatPot(:,5),pot,hdferr)
   if (hdferr /= 0) stop 'Failed to read local nuc q pot'


   ! Allocate space for the reciprocal-space potential operator.
   allocate (nonLocalResidualQ (numPotTypes,potDim))

   ! Read in the reciprocal-space potential operator.
   call h5dread_f (nonLocalResidualQ_did,H5T_NATIVE_DOUBLE,&
         & nonLocalResidualQ(:,:),potTypesPot,hdferr)
   if (hdferr /= 0) stop 'Failed to read local residual q pot'



   ! Operate on the residual charge. The residual (a.k.a. polar) charge
   !   influence converges much faster in reciprocal space than it does in real
   !   space. This is why it was separated out. Read about Ewald summation.
   ! Note that elecStatPot(:,6) = Electrostatic potential from the residual
   !   charge located on other sites.
   do i = 1,potDim
      elecStatPot(i,6) = sum(nonLocalResidualQ(:numPotTypes,i) * &
            & generalRho(:numPotTypes,6))
   enddo

   ! Deallocate space for the reciprocal-space potential operator.
   deallocate (nonLocalResidualQ)


   ! Calculate the total electrostatic energy minus that of the isolated cores.
   !   The total charge, with the nuclear charge simulated by the shortest
   !   range gaussian on each site, is integrated against all the potentials
   !   produced by the non-local and reciprocal space operators.  Only the
   !   valence charge is integrated against the potentials produced by the
   !   local operators.

   ! 

   ! Note:
   ! TOTAL  = core+valence and up+down e- charge (inc. deviations from neut)
   ! CORE   = core and up+down e- charge (fully neutral always)
   ! VALE   = valence and up+down e- charge (not neutral)
   ! NUC    = nuclear (protonic) charge
   ! NEUT_T = neutralized modification of TOTAL e- charge
   ! NEUT_V = neutralized modification of VALE e- charge
   ! RES    = residual e- charge (i.e. the difference between TOTAL and NEUT_T)

   ! Additional IMPORTANT reminder: The electrostatic energy that is computed
   !   here excludes any contribution from nuclear-nuclear, nuclear-core, and
   !   core-core charge interactions. Those contributions are constants.

   ! The expanded formula looks like (using a <||> term #comment form):
   ! E = <TOTAL+NUC|nonlocal_neut|NEUT_T> # Neut part of e- --> other sites
   !   + <TOTAL+NUC|nonlocal_res|RES> # Residual part of e- --> other sites
   !   + <VALE|local_neut|NEUT_V> # Neut part of vale e- --> local site vale
   !   + <VALE|local_neut|CORE> # Neut core e- --> local site vale
   !   + <VALE|local_nuc|NUC> # Nuclei p --> local site vale
   !   + <NEUT_V|local_neut|CORE> # Neut core e- --> local site neut-ized vale
   !   + <NEUT_V|local_nuc|NUC> # Nuclei p --> local site neutralized vale
   elecStatEnergy = sum( &
         & 0.5_double * (generalRho(:,1) + generalRho(:,7)) * & ! TOTAL+NUC
         & (elecStatPot(:,1) + elecStatPot(:,6)) & ! NEUT_T + RES (non-local)
         & + 0.5_double * generalRho(:,2) * & ! VALE, next: NEUT_V, CORE, NUC
         & (elecStatPot(:,2) + elecStatPot(:,3) + elecStatPot(:,5)) & ! local
         & + 0.5_double * (generalRho(:,2) - generalRho(:,5)) * & ! NEUT_V
         & (elecStatPot(:,3) + elecStatPot(:,5))) ! CORE + NUC

   ! Interpretation follows: The electrostatic energy is a sum of different
   !   components where each component represents the response of some part of
   !   the total charge to the potential created by some other part of the
   !   total charge (where "total" includes nuclear and electronic parts).
   ! The "first level parts" are the nuclei and the electrons.
   ! The "second level parts" are the core and valence electrons.
   ! The "third level parts" are the valence electrons that (with the help of
   !   the core electrons) would exactly cancel out the nuclear charges and
   !   the residual valence electrons representing the bit of valence charge
   !   that "goes beyond" the isolated atom neutral condition.
   ! The interactions take the form of "local" and "nonlocal". Essentially, we
   !   need to have certain parts (nuclei, core, neutral vale, residual vale)
   !   interact with other parts in a local and nonlocal way *excluding*
   !   and self interactions. Note that non-local means that charge centered
   !   on one site experiences a potential due to charges on other sites and
   !   that local refers to charge centered on one site experiencing a
   !   potential due to charges on the same site (just not with itself
   !   obviously).
   ! To get the useful electrostatic energy, we need to get the energy of the
   !   nuclei, core, neutral valence, and non-neutral (regular) valence in the
   !   potentials created by certain charge distributions. The charge "feeling"
   !   a potential is on the "left" and the charge creating a potential is on
   !   the "right" of the above set of terms in <||> format. All terms are
   !   detailed below.
   ! NUC in the potential due to ANYTHING
   !   nuclear-nonlocal_nuclear : NOT INCLUDED CONSTANT
   !   nuclear-local_nuclear: NOT INCLUDED SELF INTERACTION
   !   nuclear-nonlocal_core : The first term includes a NUC on the left and a
   !        NEUT_T on the right which itself includes the CORE and neutralized
   !        valence (NEUT_V). The operator is nonlocal. SHOULD NOT BE INCLUDED
   !        BECAUSE THIS IS CONSTANT.
   !   nuclear-local_core : NOT INCLUDED CONSTANT
   !   nuclear-nonlocal_neut_vale : The first term includes a NUC on the left
   !        and a NEUT_T on the right which itself includes the neutral
   !        valence. The operator is nonlocal.
   !   nuclear-nonlocal_res_vale : The second term includes a NUC on the left
   !        and a RES on the right. The operator is nonlocal.
   ! CORE in the potential due to ANYTHING
   !   core-nonlocal_core : The first term includes a TOTAL on the left which
   !        itself includes the core and a NEUT_T on the right which itself
   !        includes the CORE. The operator is nonlocal. SHOULD NOT BE INCLUDED
   !        BECAUSE THIS IS CONSTANT.
   !   core-local_core : NOT INCLUDED CONSTANT
   !   core-nonlocal_neut_vale : The first term includes a TOTAL on the left
   !        which itself includes the CORE and a NEUT_T on the right which
   !        itself includes the neutral valence NEUT_V. The operator is
   !        nonlocal.
   !   core-nonlocal_res_vale : The second term includes a TOTAL on the left
   !        which itself includes the CORE and a RES on the right. The operator
   !        is nonlocal.
   !   core-local_neut_vale : The sixth term includes a NEUT_V on the left and a
   !        CORE on the right. The operator is local.
   !   core-local_res_vale : The fourth term includes VALE on the left which
   !        combines NEUT_V and RES and a CORE on the right. The operator is
   !        local.
   !   


   ! Correct for the difference between the actual and fitted valence charge
   !   integrated against the gaussian screened nucleus.

   elecStatEnergyDiff = elecStatEnergy
   elecStatEnergy = elecStatEnergy - sum(generalRho(:,2) * elecStatPot(:,4)) + &
         & nucPotTrace(1) ! UP+DOWN or Total
   elecStatEnergyDiff = elecStatEnergyDiff - elecStatEnergy


   write (20,*) 'Difference between actual and fitted charge:', &
         & elecStatEnergyDiff
   call flush (20)

   ! Now, we compute the effective electrostatic potential that a given
   !   electron would feel for the current fitted charge distribution. This
   !   potential will be comprised of influences from other electrons divided
   !   into local and nonlocal contributions from core and valence electrons.
   !   The valence electrons are further subdivided into neutralized and
   !   residual contributions. The nuclei into local and nonlocal contributions
   !   as well. The nonlocal contributions are
   ! The local core contribution comes from elecStatPot(:,3).
   ! The nonlocal core contribution comes from a portion of elecStatPot(:,1).
   ! The nonlocal neutralized valence come from the other portion of
   !   elecStatPot(:,1).
   ! The local neutralized valence comes from elecStatPot(:,2).
   ! The local nuclear contribution comes from elecStatPot(:,5).
   ! The nonlocal nuclear contribution comes from elecStatPot(:,4).
   ! The residual valence part comes from elecStatPot(:,6). (I am suspecting
   !   that this is done separately for numerical accuracy reasons.)


   ! NOW, do a least squares fit for the potential with the overlap matrix.
   !   charge located on other sites. Thus: sum(elecStat(:,1,2,3,5,6)) = Total
   !   electrostatic potential. If we then subtract the gaussian screened
   !   nuclear charge (elecStat(:,4)) from the above sum we will obtain the
   !   electrostatic contribution to the potential needed for SCF calculations.
   elecStatPot(:,1) = elecStatPot(:,1) + elecStatPot(:,2) + elecStatPot(:,3) - &
         &            elecStatPot(:,4) + elecStatPot(:,5)
   elecStatPot(:,2) = elecStatPot(:,6)

   ! Copy the upper triangle of the potAlphaOverlap matrix to prevent its
   !   destruction.
   do i = 1, potDim
      tempOverlap(1:i,i) = potAlphaOverlap(1:i,i)
   enddo

   ! Solve with LAPACK dposvx
   call solveDPOSVX (potDim,2,tempOverlap,potDim,elecStatPot(:,:2),info)

   if (info /= 0) then
      write (20, *) 'Electrostatics dposvx failed. INFO= ', info
      stop
   endif

   ! Begin the exchange-correlation potential fitting

   ! The maxNumRayPoints was determined already in the access to the setup HDF5
   !   data.  Here we just copy the value from points(1).
   maxNumRayPoints = points(1)

   ! Allocate space to hold the exchange-correlation potential, exchange
   !   correlation radial weights, Rho Matrix Operator, and the resultant Rho.
   allocate (exchCorrPot  (potDim,2+spin))
   ! If running GGA calculation space is allocated for the first and second
   ! derivatives of the exchange rho operator.
   if (GGA == 0) then
      numOpValues = 1
   else
      numOpValues = 10
   endif
   allocate (exchCorrRho     (numOpValues,maxNumRayPoints))
   allocate (exchCorrRhoCore (numOpValues,maxNumRayPoints))
   allocate (exchRhoOp       (potDim,maxNumRayPoints,numOpValues))

   if (spin == 1) then
      allocate (exchCorrRhoSpin (1,maxNumRayPoints))
   else
      allocate (exchCorrRhoSpin (4,maxNumRayPoints))
   endif

   allocate (radialWeight (maxNumRayPoints))


   ! Initialize the exchange-correlation potential accumulator
   exchCorrPot(:,:) = 0.0_double

   do i = 1, numPotSites
!errorCount = 0
      currentType = potSites(i)%potTypeAssn

      ! Cycle to the next potential site if the covalent radius of the
      !   potential type for this site is too small.
      if (potTypes(currentType)%covalentRadius < smallThresh) cycle

      call h5dread_f (numPoints_did(i),H5T_NATIVE_INTEGER,&
            & numRayPoints,numPoints,hdferr)
      if (hdferr /= 0) stop 'Failed to read num ray points'

      call h5dread_f (radialWeight_did(i),H5T_NATIVE_DOUBLE,&
            & radialWeight(:),points,hdferr)
      if (hdferr /= 0) stop 'Failed to read radial weights'

      call h5dread_f (exchRhoOp_did(i),H5T_NATIVE_DOUBLE,&
            & exchRhoOp(:,:,:),potPoints,hdferr)
      if (hdferr /= 0) stop 'Failed to read exch rho operator'


      ! The exchange correlation matrix times the charge vector produces the
      !   the real space charge
      do j = 1, numRayPoints
         coreSum  = sum(exchRhoOp(:potDim,j,1) * generalRho(:potDim,3))
         totalSum = sum(exchRhoOp(:potDim,j,1) * generalRho(:potDim,1))

         if (coreSum < 0.0_double) then
            coreSum = 0.0_double
         endif
         exchCorrRhoCore (1,j) = coreSum

         if (totalSum < 0.0_double) then
            totalSum = 0.0_double
         endif
         exchCorrRho (1,j) = totalSum

         if (spin == 2) then
            spinDiffSum = sum(exchRhoOp(:potDim,j,1) * generalRho(:potDim,8))
            if (abs(spinDiffSum) > totalSum) then
!errorCount = errorCount + 1
!               write (20,*) "A charge difference between up and down spins ",&
!                     & "has been found that is greater"
!               write (20,*) "than the total charge. This is a numerical error",&
!                     & " that should only occur"
!               write (20,*) "for very small charges. In this case:"
!               write (20,*) "The total charge is:",totalSum
!               write (20,*) "The spin difference is:",spinDiffSum
               spinDiffSum = sign(totalSum,spinDiffSum)/10.0_double
!write (20,*) "exchRhoOp generalRho1 generalRho8 i,j=",i,j
!do k = 1, potDim
!write (20,*) exchRhoOp(k,j,1),generalRho(k,1),generalRho(k,8)
!enddo
            endif
            exchCorrRhoSpin(1,j) = spinDiffSum  ! Valence spin difference
         else
            exchCorrRhoSpin(1,j) = 0.0_double
         endif

         if (GGA == 1) then
            SX = -1.0_double * sum(2.0_double * exchRhoOp(:potDim,j,2) * &
                  & generalRho(:potDim,1))
            SY = -1.0_double * sum(2.0_double * exchRhoOp(:potDim,j,3) * &
                  & generalRho(:potDim,1))
            SZ = -1.0_double * sum(2.0_double * exchRhoOp(:potDim,j,4) * &
                  & generalRho(:potDim,1))
            SXX = sum(4.0_double * exchRhoOp(:potDim,j,5) &
                  & - 2.0_double * exchRhoOp(:potDim,j,1) * &
                  & generalRho(:potDim,1))
            SXY = sum(4.0_double * exchRhoOp(:potDim,j,6))
            SXZ = sum(4.0_double * exchRhoOp(:potDim,j,7))
            SYY = sum(4.0_double * exchRhoOp(:potDim,j,8) &
                  &- 2.0_double * exchRhoOp(:potDim,j,1) * &
                  & generalRho(:potDim,1))
            SYZ = sum(4.0_double * exchRhoOp(:potDim,j,9))
            SZZ = sum(4.0_double * exchRhoOp(:potDim,j,10) &
                  &- 2.0_double * exchRhoOp(:potDim,j,1) * &
                  & generalRho(:potDim,1))


            SXC = -1.0_double * sum(2.0_double * exchRhoOp(:potDim,j,2) * &
                  & generalRho(:potDim,3))
            SYC = -1.0_double * sum(2.0_double * exchRhoOp(:potDim,j,3) * &
                  & generalRho(:potDim,3))
            SZC = -1.0_double * sum(2.0_double * exchRhoOp(:potDim,j,4) * &
                  & generalRho(:potDim,3))
            SXXC = sum(4.0_double * exchRhoOp(:potDim,j,5) &
                  & - 2.0_double * exchRhoOp(:potDim,j,1) * &
                  & generalRho(:potDim,3))
            SXYC = sum(4.0_double * exchRhoOp(:potDim,j,6))
            SXZC = sum(4.0_double * exchRhoOp(:potDim,j,7))
            SYYC = sum(4.0_double * exchRhoOp(:potDim,j,8) &
                  &- 2.0_double * exchRhoOp(:potDim,j,1) * &
                  & generalRho(:potDim,3))
            SYZC = sum(4.0_double * exchRhoOp(:potDim,j,9))
            SZZC = sum(4.0_double * exchRhoOp(:potDim,j,10) &
                  &- 2.0_double * exchRhoOp(:potDim,j,1) * &
                  & generalRho(:potDim,3))

            if (spin == 2) then
               SXS = -1.0_double * sum(2.0_double * exchRhoOp(:potDim,j,2) * &
                     & generalRho(:potDim,8))
               SYS = -1.0_double * sum(2.0_double * exchRhoOp(:potDim,j,3) * &
                     & generalRho(:potDim,8))
               SZS = -1.0_double * sum(2.0_double * exchRhoOp(:potDim,j,4) * &
                     & generalRho(:potDim,8))
            endif
         endif

         if (GGA == 1) then
            if (SXC < smallThresh) then
               SXC = smallThresh
            endif
            if (SYC < smallThresh) then
               SYC = smallThresh
            endif
            if (SZC < smallThresh) then
               SZC = smallThresh
            endif
             if (SXS < smallThresh) then
               SXS = smallThresh
            endif
            if (SYS < smallThresh) then
               SYS = smallThresh
            endif
            if (SZS < smallThresh) then
               SZS = smallThresh
            endif
         endif

         ! Note that the indexing scheme has changed from the one used
         ! in exchRhoOp(:,:,10) in order to incorporate coreSum and spinDiffSum
         if (GGA == 1) then
            exchCorrRho (2,j) = SX
            exchCorrRho (3,j) = SY
            exchCorrRho (4,j) = SZ
            exchCorrRho (5,j) = SXX
            exchCorrRho (6,j) = SXY
            exchCorrRho (7,j) = SXZ
            exchCorrRho (8,j) = SYY
            exchCorrRho (9,j) = SYZ
            exchCorrRho (10,j) = SZZ

            exchCorrRhoCore (2,j) = SXC
            exchCorrRhoCore (3,j) = SYC
            exchCorrRhoCore (4,j) = SZC
            exchCorrRhoCore (5,j) = SXXC
            exchCorrRhoCore (6,j) = SXYC
            exchCorrRhoCore (7,j) = SXZC
            exchCorrRhoCore (8,j) = SYYC
            exchCorrRhoCore (9,j) = SYZC
            exchCorrRhoCore (10,j) = SZZC
            if (spin == 2) then
               exchCorrRhoSpin (2,j) = SXS
               exchCorrRhoSpin (3,j) = SYS
               exchCorrRhoSpin (4,j) = SZS
            endif
         endif
      enddo ! j = 1, numRayPoints

      ! Define the exchange correlation functions
      ! 100 = Wigner
      ! 101 = Ceperley-Alder
      ! 102 = Hedin-Lundqvist
      ! 150 = Ceperley-Alder
      ! 151 = von Barth-Hedin
      ! 152 = unknown
      ! 200 = PBE96 (Perdew, Burke, and Enzerhof)

      if (xcCode == 100) then
            do j = 1, numRayPoints
               ! These functions have been converted to subroutines because the
               !   return value of the second array index value was 0.0 on the
               !   ia64 HP-UX machine sirius.  This seems to work.
               call wignerXC(exchCorrRho(1,j),currentExchCorrPot(1:2))
               call wignerXCEnergy(exchCorrRhoCore(1,j),currentExchCorrPot(3))
               do k = 1,3
                  exchCorrPot(:,k) = exchCorrPot(:,k) + &
                        & radialWeight(j) * currentExchCorrPot(k) * &
                        & exchRhoOp(:potDim,j,1)
               enddo
            enddo
      elseif (xcCode == 101) then
            do j = 1, numRayPoints
               call ceperleyAlderXC(exchCorrRho(1,j),currentExchCorrPot(1:2))
               call ceperleyAlderXCEnergy(exchCorrRhoCore(1,j),&
                     & currentExchCorrPot(3))
               do k = 1,3
                  exchCorrPot(:,k) = exchCorrPot(:,k) + &
                        & radialWeight(j) * currentExchCorrPot(k) * &
                        & exchRhoOp(:potDim,j,1)
               enddo
            enddo
      elseif (xcCode == 102) then
            do j = 1, numRayPoints
               call hedinLundqvistXC(exchCorrRho(1,j),&
                     & currentExchCorrPot(1:2))
               call hedinLundqvistXCEnergy(exchCorrRhoCore(1,j),&
                     & currentExchCorrPot(3))
               do k = 1,3
                  exchCorrPot(:,k) = exchCorrPot(:,k) + &
                        & radialWeight(j) * currentExchCorrPot(k) * &
                        & exchRhoOp(:potDim,j,1)
               enddo
            enddo
      elseif (xcCode == 150) then
            ! Ceperley and Alder Exchange-Correlation (LSDA)
            do j = 1, numRayPoints
               call ceperleyAlderSP(exchCorrRho(1,j),exchCorrRhoSpin(1,j),&
                     & exchCorrRhoCore(1,j),currentExchCorrPot(:))
!write (20,*) "i,j,pts",i,j,numRayPoints
!write (20,*) exchCorrRho(1,j),exchCorrRhoSpin(1,j),&
!   & exchCorrRhoCore(1,j),currentExchCorrPot(:)
               do k = 1,4
                  exchCorrPot(:,k) = exchCorrPot(:,k) + &
                        & radialWeight(j) * currentExchCorrPot(k) * &
                        & exchRhoOp(:,j,1)
               enddo
            enddo
!write (20,*) "exchCorrPot(:,1) i=",i
!write (20,*) exchCorrPot(:,1)
!write (20,*) "exchCorrPot(:,2) i=",i
!write (20,*) exchCorrPot(:,2)
!write (20,*) "exchCorrPot(:,3) i=",i
!write (20,*) exchCorrPot(:,3)
!write (20,*) "exchCorrPot(:,4) i=",i
!write (20,*) exchCorrPot(:,4)

      elseif (xcCode == 151) then
            ! von Barth and Hedin Exchange-Correlation (LSDA)
            do j = 1, numRayPoints
               call vonBarthHedin(exchCorrRho(1,j),exchCorrRhoSpin(1,j),&
                     & exchCorrRhoCore(1,j),currentExchCorrPot(:))
               do k = 1,4
                  exchCorrPot(:,k) = exchCorrPot(:,k) + &
                        & radialWeight(j) * currentExchCorrPot(k) * &
                        & exchRhoOp(:potDim,j,1)
               enddo
            enddo
      elseif (xcCode == 152) then
            ! Old SPmain.for Exchange Correlation (LSDA) (Should be like the
            !   von barth and Hedin above.)
            do j = 1, numRayPoints
               call oldEXCORR(exchCorrRho(1,j),exchCorrRhoSpin(1,j),&
                     & exchCorrRhoCore(1,j),currentExchCorrPot(:))
               do k = 1,4
                  exchCorrPot(:,k) = exchCorrPot(:,k) + &
                        & radialWeight(j) * currentExchCorrPot(k) * &
                        & exchRhoOp(:potDim,j,1)
               enddo
            enddo
      elseif (xcCode == 200) then
            ! PBE96 (Perdew, Burke, and Enzerhof), should be like
            ! GGA in old ggamain.for
            do j = 1, numRayPoints

               call pbe96(exchCorrRho(1,j),exchCorrRho(2,j), &
                  & exchCorrRho(3,j),exchCorrRho(4,j), &
                  & exchCorrRho(5,j),exchCorrRho(6,j),exchCorrRho(7,j),&
                  & exchCorrRho(8,j), &
                  & exchCorrRho(9,j),exchCorrRho(10,j),currentExchCorrPot(1:2))
               call pbe96(exchCorrRhoCore(1,j),exchCorrRhoCore(2,j), &
                  & exchCorrRhoCore(3,j),exchCorrRhoCore(4,j), &
                  & exchCorrRhoCore(5,j),exchCorrRhoCore(6,j), &
                  & exchCorrRhoCore(7,j),exchCorrRhoCore(8,j), &
                  & exchCorrRhoCore(9,j),exchCorrRho(10,j), &
                  & currentExchCorrPot(3:4))
               exchCorrPot(:,1) = exchCorrPot(:,1) + &
                  & radialWeight(j) * currentExchCorrPot(1) * &
                  & exchRhoOp(:potDim,j,1)
               exchCorrPot(:,2) = exchCorrPot(:,2) + &
                  & radialWeight(j) * currentExchCorrPot(2) * &
                  & exchRhoOp(:potDim,j,1)
               exchCorrPot(:,3) = exchCorrPot(:,3) + &
                  & radialWeight(j) * currentExchCorrPot(4) * &
                  & exchRhoOp(:potDim,j,1)
            enddo
!      elseif (XC_CODE == 250) then
!           ! pbe96spin (Perdew, Burke, and Enzerhof), should be like
!           ! GGA in old ggasmain.for
!           do j = 1, numRayPoints
!
!              call pbe96spin(exchCorrRho(1,j),exchCorrRho(3,j),exchCorrRho(4,j), &
!                 & exchCorrRho(5,j),exchCorrRho(6,j), &
!                 & exchCorrRho(7,j),exchCorrRho(8,j),exchCorrRho(9,j),exchCorrRho(10,j), &
!                 & exchCorrRho(11,j),exchCorrRho(12,j),exchCorrRho(22,j),exchCorrRho(23,j),exchCorrRho(24,j), &
!                 & currentExchCorrPot(1:2))
!              call pbe96spin(exchCorrRho(2,j),exchCorrRho(3,j),exchCorrRho(13,j), &
!                 & exchCorrRho(14,j),exchCorrRho(15,j), &
!                 & exchCorrRho(16,j),exchCorrRho(17,j),exchCorrRho(18,j),exchCorrRho(19,j), &
!                 & exchCorrRho(20,j),exchCorrRho(21,j),exchCorrRho(22,j),exchCorrRho(23,j),exchCorrRho(24,j), &
!                 & currentExchCorrPot(3:4))
!              exchCorrPot(:,1) = exchCorrPot(:,1) + &
!                 & radialWeight(j) * currentExchCorrPot(1) * &
!                 & exchRhoOp(:potDim,j,1)
!              exchCorrPot(:,2) = exchCorrPot(:,2) + &
!                 & radialWeight(j) * currentExchCorrPot(2) * &
!                 & exchRhoOp(:potDim,j,1)
!              exchCorrPot(:,3) = exchCorrPot(:,3) + &
!                 & radialWeight(j) * currentExchCorrPot(4) * &
!                 & exchRhoOp(:potDim,j,1)
!              exchCorrPot(:,4) = exchCorrPot(:,3) + &
!                 & radialWeight(j) * currentExchCorrPot(4) * &
!                 & exchRhoOp(:potDim,j,1)
!           enddo
      endif
!write (20,*) "i,error count=",i,errorCount
!call flush (20)
   enddo ! i = 1, numPotSites

   ! Allocate space to hold the exchCorrOverlap
   allocate (exchCorrOverlap(potDim,potDim))

   ! Read the exchange correlation overlap.
   call h5dread_f (exchCorrOverlap_did,H5T_NATIVE_DOUBLE,&
         & exchCorrOverlap(:,:),potPot,hdferr)
   if (hdferr /= 0) stop 'Failed to read exch corr overlap'


   ! Solve the linear set of equations with LAPACK
   call solveDPOSVX (potDim,2+spin,exchCorrOverlap,potDim,&
         & exchCorrPot(:,:2+spin),info)

   if (info /= 0) then
      write (20, *) 'Exchange-correlation dposvx failed. INFO= ', info
      stop
   endif

   ! Calculate the total exchange correlation energy minus the core exchange
   !   correlation energy to get the valence exchange correlation energy.
   exchCorrEnergy = 0.0_double
   do i = 1, potDim
      exchCorrEnergy = exchCorrEnergy + &
         & sum(potAlphaOverlap(:,i)*exchCorrPot(:,1+spin)) * generalRho(i,1) - &
         & sum(potAlphaOverlap(:,i)*exchCorrPot(:,2+spin)) * generalRho(i,3)
   enddo

   ! Compute the spin dependent exchange correlation energy. Note that because
   !   generalRho(:,8) holds a charge spin difference the contribution from the
   !   core is identically zero. Therefore, nothing needs to be subtracted off.
   if (spin == 2) then
      exchCorrEnergyDiff = 0.0_double
      do i = 1, potDim
         exchCorrEnergyDiff = exchCorrEnergyDiff &
            & + sum(potAlphaOverlap(:,i)*exchCorrPot(:,1+spin))*generalRho(i,8)
      enddo

      ! Compute the up and then down contributions to the total exchange
      !   correlation energy. Up = (Total + Diff)/2; Down = (Total - Diff)/2
      exchCorrEnergySpin(1) = (exchCorrEnergy+exchCorrEnergyDiff) / 2.0_double
      exchCorrEnergySpin(2) = (exchCorrEnergy-exchCorrEnergyDiff) / 2.0_double
   endif


   ! Compute (copy) the kinetic energy and the spin polarized contributions to
   !   the kinetic energy.
   kineticEnergy = kineticEnergyTrace(1) ! UP + DOWN or TOTAL
   if (spin == 2) then
      kineticEnergySpin(1) = &
            & (kineticEnergyTrace(1) + kineticEnergyTrace(2)) / 2.0_double
      kineticEnergySpin(2) = &
            & (kineticEnergyTrace(1) - kineticEnergyTrace(2)) / 2.0_double
   endif

   ! Calculate the total energy and the spin polarized contributions to the
   !   total energy.
   totalEnergy = kineticEnergy + elecStatEnergy + exchCorrEnergy
   if (spin == 2) then
      totalEnergySpin(1) = kineticEnergySpin(1) + elecStatEnergy &
            & + exchCorrEnergySpin(1)
      totalEnergySpin(2) = kineticEnergySpin(2) + elecStatEnergy &
            & + exchCorrEnergySpin(2)
   endif

   if (totalEnergy > 0.0_double) then
      write (20,*) "Total energy is positive!"
      write (20,*) "Your structure might have a serious problem. Aborting."
      stop
   endif


   ! Deallocate exchange correlation arrays now that the total energy is known.
   deallocate (exchCorrRho)
   deallocate (exchCorrRhoCore)
   deallocate (exchCorrRhoSpin)
   deallocate (exchCorrOverlap)

   ! Allocate space to hold data structures for determining the next set of
   !   potential coefficients and for measuring the degree of convergence that
   !   we have achieved since the previous iteration. Note that we only need to
   !   allocate these entities once because they are used in every SCF iteration
   !   and they don't consume too much space so they will not "get in the way"
   !   of other much larger data structures that are needed in other parts of
   !   the SCF process (e.g. diagonalization).
   ! Note importantly, that the indices for usedPotCoeffs and guessedPotCoeffs
   !   in the second array index are one higher than the numbers expressed in
   !   the Anderson paper. That is, index 1 is for x^l+0 and index 2 is for
   !   x^l-1 and index 3 is for x^l-2, etc.
   if (.not. allocated(usedPotCoeffs)) then
      allocate (convergenceRecord(feedbackLevel+1))
      allocate (totalEnergyRecord(feedbackLevel+2,spin))
      convergenceRecord(:) = 100000.0_double ! Some large number
      totalEnergyRecord(:,:) = 0.0_double

      allocate (usedPotCoeffs(potDim,feedbackLevel+1,spin)) ! Anderson's x
      allocate (guessedPotCoeffs(potDim,feedbackLevel+1,spin)) ! Anderson's y

      allocate (tempPotCoeffs(potDim,spin))
      allocate (tempGuessedPotCoeffs(potDim,feedbackLevel+1,spin))
      allocate (tempUsedPotCoeffs(potDim,feedbackLevel+1,spin))

      do i = 1, spin
         do j = 1, feedbackLevel+1
            usedPotCoeffs(:,j,i) = potCoeffs(:,i)
            guessedPotCoeffs(:,j,i) = elecStatPot(:,1) + elecStatPot(:,2) &
                  & + exchCorrPot(:,i)
         enddo
      enddo
   endif

   ! Allocate space to hold the output potentials
!   allocate (outputPot (potDim,spin))
!   if (.not.allocated(xl0)) then
!      allocate (xl0 (potDim,spin))
!      allocate (xl1 (potDim,spin))
!      allocate (xl2 (potDim,spin))
!      allocate (yl0 (potDim,spin))
!      allocate (yl1 (potDim,spin))
!      allocate (yl2 (potDim,spin))
!
!      do i = 1, spin
!         xl1(:,i) = potCoeffs(:,i)
!         xl2(:,i) = potCoeffs(:,i)
!         yl1(:,i) = elecStatPot(:,1) + elecStatPot(:,2) + exchCorrPot(:,i)
!         yl2(:,i) = elecStatPot(:,1) + elecStatPot(:,2) + exchCorrPot(:,i)
!      enddo
!   endif


   allocate (potDifference(potDim,spin)) ! For measuring convergence.
   do i = 1, spin

      ! NOT BEING USED NOW, DIDN'T WORK SO GREAT. That may be due to not
      !   matching the changes in energy with the "correct" pair of used vs.
      !   guessed coefficients. Not sure so I'm leaving it for now ...
      ! If the total energy of the current iteration is worse (higher) than the
      !   total energy of each of the previously saved iterations, then we will
      !   reverse the assignment of used and guessed potential coefficients for
      !   this iteration. This will have the effect of "pushing against" the
      !   particular blending that produced this slightly worse energy.
!write (20,*) totalEnergy, totalEnergyRecord(:)
if (spin == 1) then
   totalEnergyRecord(1,:) = totalEnergy
else
   totalEnergyRecord(1,:) = totalEnergySpin(:)
endif
!if (totalEnergy < maxval(totalEnergyRecord(:))) then
!   totalEnergyImprovement = 1
!else
!   totalEnergyImprovement = 0
!endif

!      if (totalEnergyImprovement == 1) then
         usedPotCoeffs(:,1,i) = potCoeffs(:,i) ! Anderson's x^l
         guessedPotCoeffs(:,1,i) = elecStatPot(:,1) + elecStatPot(:,2) &
               & + exchCorrPot(:,i) ! Anderson's y^l
!      else
!         guessedPotCoeffs(:,1,i) = potCoeffs(:,i) ! Anderson's x^l
!         usedPotCoeffs(:,1,i) = elecStatPot(:,1) + elecStatPot(:,2) &
!               & + exchCorrPot(:,i) ! Anderson's y^l
!      endif

      ! Form the difference between the potential coefficients that have been
      !   guessed for the next iteration and those that were used in the current
      !   iteration.
      potDifference(:,i) = guessedPotCoeffs(:,1,i) - usedPotCoeffs(:,1,i)
   enddo
!do i = 1, potDim
!write (20,*) "exchCorr",exchCorrPot(i,1),exchCorrPot(i,2)
!enddo
!do i = 1, potDim
!write (20,*) "elecStat1,2",elecStatPot(i,1),elecStatPot(i,2)
!enddo

!   do i = 1, spin
!      xl0(:,i) = potCoeffs(:,i)
!      yl0(:,i) = elecStatPot(:,1) + elecStatPot(:,2) + exchCorrPot(:,i)
!
!      ! Form the difference between last iteration and the current new
!      !   calculation of the potential.  (Note that this calculation will be
!      !   mixed later with the potential of last iteration to make the
!      !   convergence more smooth.)
!      outputPot(:,i) = yl0(:,i) - xl0(:,i)
!   enddo


   ! Deallocate matrices that have been copied to temporary arrays.
   deallocate (elecStatPot)
   deallocate (exchCorrPot)


   ! Begin real-space self consistancy analysis


   ! Initialize a site index number to track the number of sites that are
   !   included.  Unless there are any sites that have a very small covalent
   !   radius this will track the same values as the loop index i.
   siteIndex = 0

   ! Initialize the value of the convergence delta that will be tested to see
   !   if the potentials have converged.
   testableDelta = 0.0_double

   ! Allocate matrices to obtain the testable delta
   allocate (realSpacePotDiff (maxNumRayPoints,spin))
   allocate (averageDelta     (numPotSites))
   allocate (maxDelta         (numPotSites))


   do i = 1, numPotSites
      currentType = potSites(i)%potTypeAssn

      ! Cycle if the covalent radius for the potential type of this potential
      !   site is sufficiently small.
      if (potTypes(currentType)%covalentRadius < smallThresh) cycle

      siteIndex = siteIndex + 1

      ! Read the mesh information for the current potential site.
      call h5dread_f (numPoints_did(i),H5T_NATIVE_INTEGER,&
            & numRayPoints,numPoints,hdferr)
      if (hdferr /= 0) stop 'Failed to read num ray points #2'
      call h5dread_f (radialWeight_did(i),H5T_NATIVE_DOUBLE,&
            & radialWeight(:),points,hdferr)
      if (hdferr /= 0) stop 'Failed to read radial weights #2'
      call h5dread_f (exchRhoOp_did(i),H5T_NATIVE_DOUBLE,&
            & exchRhoOp(:,:,:),potPoints,hdferr)
      if (hdferr /= 0) stop 'Failed to read exch rho operator #2'


      ! Calculate the exchange-correlation rho matrix operator times the
      !   potential vector difference to get the difference in real space.
      do j = 1, spin
         do k = 1, numRayPoints
!            realSpacePotDiff(k,j) = sum(exchRhoOp(:,k,1) * outputPot(:,j))
            realSpacePotDiff(k,j) = sum(exchRhoOp(:,k,1) * potDifference(:,j))
         enddo

         ! Determine the delta to be tested.
         maxDelta(siteIndex)     = maxval(realSpacePotDiff(:numRayPoints,j))
         radialWeightSum         = sum (radialWeight(:numRayPoints))
         weightedPotDiff         = sum (radialWeight(:numRayPoints) * &
               & realSpacePotDiff(:numRayPoints,j)**2)
         averageDelta(siteIndex) = sqrt(weightedPotDiff / radialWeightSum)
!write (20,*) "averageDelta,i,j",averageDelta(siteIndex),i,j
!write (20,*) "maxDelta",maxDelta(siteIndex)
         testableDelta           = max (averageDelta(siteIndex),testableDelta)
      enddo
   enddo

   ! Deallocate arrays to obtain the testableDelta
   deallocate (potDifference)
   deallocate (realSpacePotDiff)
   deallocate (averageDelta)
   deallocate (maxDelta)
   deallocate (exchRhoOp)
   deallocate (radialWeight)
!   deallocate (outputPot)

   ! Generate the potentials for the next iteration using the method D.G.
   !   Anderson, J. Assoc. Comp. Mach. 12, 547 (1965).

   ! The parameter here referred to as the feedbackLevel is called M by
   !   Anderson.  If M=0, the straight 'relaxed' iteration is used.  If M=2,
   !   the 2x2 determinant in anderson's method is checked.  If it is too
   !   small, the first order M=1 extrapolation is used instead.  M is
   !   effectively reduced to 0 or 1 if the iteration is 0, 1, or 2.

!   allocate ( rl0 (potDim))
!   allocate ( rl1 (potDim))
!   allocate ( rl2 (potDim))
!   allocate (drl1 (potDim))
!   allocate (drl2 (potDim))


   ! At this point we must determine what potential coefficients to use during
   !   the next iteration. The D.G. Anderson method implemented in the
   !   blendPotentials subroutine may be used.

   ! Blend the currently used potential, the current guess, and (possibly) a set
   !   number of previous potential functions in an effort to accelerate the
   !   convergence while retaining smooth character to the convergence process.
!write (20,*) "Iteration: ",currIteration
!write (20,*) "guessed11 used11 guessed21 used21 guessed31 used31"
!do i = 1, potDim
!write (20,*) guessedPotCoeffs(i,1,1), usedPotCoeffs(i,1,1)
!enddo
!do i = 1, potDim
!write (20,*) guessedPotCoeffs(i,2,1), usedPotCoeffs(i,2,1)
!enddo
!do i = 1, potDim
!write (20,*) guessedPotCoeffs(i,3,1), usedPotCoeffs(i,3,1)
!enddo
!write (20,*) "guessed12 used12 guessed22 used22 guessed32 used32"
!do i = 1, potDim
!write (20,*) guessedPotCoeffs(i,1,2), usedPotCoeffs(i,1,2)
!enddo
!do i = 1, potDim
!write (20,*) guessedPotCoeffs(i,2,2), usedPotCoeffs(i,2,2)
!enddo
!do i = 1, potDim
!write (20,*) guessedPotCoeffs(i,3,2), usedPotCoeffs(i,3,2)
!enddo
!do i = 1, potDim
!write (20,*) yl0(i,1), xl0(i,1)
!enddo
!do i = 1, potDim
!write (20,*) yl1(i,1), xl1(i,1)
!enddo
!do i = 1, potDim
!write (20,*) yl2(i,1), xl2(i,1)
!enddo
!write (20,*) "guessed12 used12 guessed22 used22 guessed32 used32"
!do i = 1, potDim
!write (20,*) yl0(i,2), xl0(i,2)
!enddo
!do i = 1, potDim
!write (20,*) yl1(i,2), xl1(i,2)
!enddo
!do i = 1, potDim
!write (20,*) yl2(i,2), xl2(i,2)
!enddo

!   do i = 1, spin
!      call blendPotentialsSCF(1,potDim,potCoeffs(:,i),&
!            & guessedPotCoeffs(:,:,i),usedPotCoeffs(:,:,i),totalEnergyRecord)
!   enddo
   do i = 1, spin
      call blendPotentialsSCF(1,potDim,potCoeffs(:,i),&
            & guessedPotCoeffs(:,:,i),usedPotCoeffs(:,:,i),&
            & totalEnergyRecord(:,i))
   enddo
!   if (currIteration > feedbackLevel+2) then
!      do i = 1, spin
!         call blendPotentialsTE(1,potDim,potCoeffs(:,i),&
!               & guessedPotCoeffs(:,:,i),usedPotCoeffs(:,:,i),&
!               & totalEnergyRecord(:,i))
!      enddo
!   endif

   call shiftPotentials(guessedPotCoeffs,usedPotCoeffs,totalEnergyRecord)

!   if (spin == 1) then ! Do the simple spin non-polarized case.
!      call blendPotentialsSCF(1,potDim,potCoeffs(:,1),&
!            & guessedPotCoeffs(:,:,1), usedPotCoeffs(:,:,1),totalEnergyRecord)
!   else
!
!      tempGuessedPotCoeffs(:,:,1) = guessedPotCoeffs(:,:,1) &
!            & + guessedPotCoeffs(:,:,2)
!      tempGuessedPotCoeffs(:,:,2) = guessedPotCoeffs(:,:,1) &
!            & - guessedPotCoeffs(:,:,2)
!      tempUsedPotCoeffs(:,:,1) = usedPotCoeffs(:,:,1) &
!            & + usedPotCoeffs(:,:,2)
!      tempUsedPotCoeffs(:,:,2) = usedPotCoeffs(:,:,1) &
!            & - usedPotCoeffs(:,:,2)
!
!      do i = 1, spin
!         call blendPotentialsSCF(1,potDim,tempPotCoeffs(:,i),&
!               & tempGuessedPotCoeffs(:,:,i), tempUsedPotCoeffs(:,:,i),&
!               & totalEnergyRecord)
!      enddo
!
!      potCoeffs(:,1) = (tempPotCoeffs(:,1) + tempPotCoeffs(:,2))/2.0_double
!      potCoeffs(:,2) = (tempPotCoeffs(:,1) - tempPotCoeffs(:,2))/2.0_double
!      guessedPotCoeffs(:,:,1) = (tempGuessedPotCoeffs(:,:,1) &
!            & + tempGuessedPotCoeffs(:,:,2))/2.0_double
!      guessedPotCoeffs(:,:,2) = (tempGuessedPotCoeffs(:,:,1) &
!            & - tempGuessedPotCoeffs(:,:,2))/2.0_double
!      usedPotCoeffs(:,:,1) = (tempUsedPotCoeffs(:,:,1) &
!            & + tempUsedPotCoeffs(:,:,2))/2.0_double
!      usedPotCoeffs(:,:,2) = (tempUsedPotCoeffs(:,:,1) &
!            & - tempUsedPotCoeffs(:,:,2))/2.0_double
!   endif
!   jointPotCoeffs(1:potDim)                   = potCoeffs(:,1)
!   jointPotCoeffs(potDim+1:potDim*2)          = potCoeffs(:,2)
!   jointGuessedPotCoeffs(1:potDim,:)          = guessedPotCoeffs(:,:,1)
!   jointGuessedPotCoeffs(potDim+1:potDim*2,:) = guessedPotCoeffs(:,:,2)
!   jointUsedPotCoeffs(1:potDim,:)             = usedPotCoeffs(:,:,1)
!   jointUsedPotCoeffs(potDim+1:potDim*2,:)    = usedPotCoeffs(:,:,2)
!
!   call blendJointPotentials(1,potDim*2,jointPotCoeffs(:),&
!         & jointGuessedPotCoeffs(:,:),jointUsedPotCoeffs(:,:),totalEnergyRecord)
!
!   ! Disassemble both spins into separate arrays.
!   potCoeffs(:,1)          = jointPotCoeffs(1:potDim)
!   potCoeffs(:,2)          = jointPotCoeffs(potDim+1:potDim*2)
!   guessedPotCoeffs(:,:,1) = jointGuessedPotCoeffs(1:potDim,:)
!   guessedPotCoeffs(:,:,2) = jointGuessedPotCoeffs(potDim+1:potDim*2,:)
!   usedPotCoeffs(:,:,1)    = jointUsedPotCoeffs(1:potDim,:)
!   usedPotCoeffs(:,:,2)    = jointUsedPotCoeffs(potDim+1:potDim*2,:)


!   ! This is basically copied from the old code with little modification for
!   !   the names since I don't know what they mean or do in most cases.
!   do i = 1, spin
!      th1 = 0.0_double
!      th2 = 0.0_double
!      if (feedbackLevel /= 0) then
!         if (currIteration > 1) then
!            rl0(:)  = yl0(:,i) - xl0(:,i)
!            rl1(:)  = yl1(:,i) - xl1(:,i)
!            rl2(:)  = yl2(:,i) - xl2(:,i)
!            drl1(:) = rl0(:) - rl1(:)
!            drl2(:) = rl0(:) - rl2(:)
!
!            ! This part is new.  Make temp matrices holding the potAlphaOverlap
!            !   times the drl1, and drl2 vectors.
!
!            ! Initialize the summation variables
!            s11 = 0.0_double
!            s12 = 0.0_double
!            s22 = 0.0_double
!            t1  = 0.0_double
!            t2  = 0.0_double
!
!            ! The cases for the drl1 vector are done first. 
!            do j = 1, potDim
!               tempOverlap(:,j) = potAlphaOverlap(:,j) * drl1(j)
!               s11 = s11 + sum(tempOverlap(:,j) * drl1(:))
!               t1  = t1  + sum(tempOverlap(:,j) * rl0(:))
!            enddo
!
!            ! The cases for the drl2 vector are done second.
!            do j = 1, potDim
!               tempOverlap(:,j) = potAlphaOverlap(:,j) * drl2(j)
!               s12 = s12 + sum(tempOverlap(:,j) * drl1(:))
!               s22 = s22 + sum(tempOverlap(:,j) * drl2(:))
!               t2  = t2  + sum(tempOverlap(:,j) * rl0(:))
!            enddo
!            th1 = t1/s11
!            if (feedbackLevel /= 1) then
!               if (currIteration > 2) then
!
!                  ! Calculate the determinate.
!                  det = s11*s22 - s12*s12
!                  if (det/(s11*s22) >= 0.00000001_double) then
!                     th1 = ( s22*t1 - s12*t2)/det
!                     th2 = (-s12*t1 + s11*t2)/det
!                  endif
!               endif
!            endif
!         endif
!      endif
!
!      ! Form the next iteration and save the previous two.
!      potCoeffs(:,i) = &
!            & (1.0_double - relaxFactor) * &
!            & ((1.0_double-th1-th2)*xl0(:,i) + th1*xl1(:,i) + th2*xl2(:,i)) + &
!            & relaxFactor * &
!            & ((1.0_double-th1-th2)*yl0(:,i) + th1*yl1(:,i) + th2*yl2(:,i))
!
!      xl2(:,i) = xl1(:,i)
!      xl1(:,i) = xl0(:,i)
!      yl2(:,i) = yl1(:,i)
!      yl1(:,i) = yl0(:,i)
!   enddo


   ! Record the potential terms (alphas), potential coefficients, total charge
   !   density (core+valence), and valence charge density (up+down), and for
   !   spin polarized cases the difference (up-down).
   rewind (8)

   ! Write the file header that says the total number of types.
   write (8,fmt="(a9,i5)") "NUM_TYPES",numPotTypes

   ! Normally, one would expect this type of loop to go from 1 to spin.
   !   However, in this case we need to go from 1 to 2. For the spin==2 case,
   !   the reason is obvious: we need to write the spin up and spin down
   !   potential coefficients and other data. For the spin==1 case, we should
   !   only need to write one set of coefficients because there is no spin up
   !   or spin down concept here, the system is spin degenerate. However, to
   !   maintain consistency in the structure of the file we will repeat the
   !   same potential coefficients and they will have a "SPIN_DN" label. Do not
   !   be confused though, the coefficients for the SPIN_DN section are not used
   !   by a spin non-polarized calculation at all. They are just cruft that is
   !   kept around to maintain file format consistency.
   do i = 1, 2

      ! Initialize the counter for the number of terms.
      potTermCount = 0

      if (i == 1) then
         write (8,fmt="(a18)") "TOTAL__OR__SPIN_UP"
      else
         write (8,fmt="(a7)") "SPIN_DN"
      endif

      do j = 1, numPotTypes

         ! Write the type header that says the number of terms for this type.
         write (8,fmt="(i5)") potTypes(j)%numAlphas

         do k = 1, potTypes(j)%numAlphas

            ! Increment the counter.
            potTermCount = potTermCount + 1

            ! Write the term. Note that the spin==1 case always writes a
            !   potCoeffs term from index 1 while the spin==2 case can write a
            !   potCoeffs term from index i==1 or i==2.
            if (spin == 1) then
               write (8,fmt="(2(1x,e17.10),3(1x,e13.6))") &
                     & potCoeffs(potTermCount,1),potAlphas(potTermCount),&
                     & generalRho(potTermCount,1),generalRho(potTermCount,2),&
                     & 0.0_double
            else
               write (8,fmt="(2(1x,e17.10),3(1x,e13.6))") &
                     & potCoeffs(potTermCount,i),potAlphas(potTermCount),&
                     & generalRho(potTermCount,1),generalRho(potTermCount,2),&
                     & generalRho(potTermCount,8)
            endif
         enddo
      enddo
   enddo
   call flush (8)

   ! Record the potential information for any plusUJ terms.
   write (8,fmt="(a18)") "NUM__PLUSUJ__TERMS"
   write (8,fmt="(i5)") numPlusUJAtoms
   do i = 1, 2 ! Spin up and down.

      ! Write the header tag for this block of data.
      if (i == 1) then
         write (8,fmt="(a18)") "TOTAL__OR__SPIN_UP"
      else
         write (8,fmt="(a7)") "SPIN_DN"
      endif

      do j = 1, numKPoints

         ! Write the atom ID numbers followed by the plusUJ terms for each of
         !   the numPlusUJAtoms. Note again that the spin==1 case will produce
         !   redundent data in the i==2 iteration.
         do k = 1, numPlusUJAtoms
            !do k = 1, 7 ! Temp. comment until I decide to do 5x5 or 7x7 blocks.
            if (spin == 1) then
               write (8,fmt="(i5)") plusUJAtomID(k)
               write (8,fmt="(4e18.10)") plusUJ(1:4,k,1,j)
               write (8,fmt="(3e18.10)") plusUJ(5:7,k,1,j)
            else
               write (8,fmt="(i5)") plusUJAtomID(k)
               write (8,fmt="(4e18.10)") plusUJ(1:4,k,i,j)
               write (8,fmt="(3e18.10)") plusUJ(5:7,k,i,j)
            endif
         enddo
      enddo
   enddo
   call flush (8)

   ! Compute the total magnetic moment of the system and of individual
   !   potential types.  (Only done for spin polarized calculations.)
   if (spin == 2) then

      ! Initialize the accumlators for the total system magnetic moment and the
      !    individual type's magnetic moment.
      totalMagneticMoment = 0.0_double
      typesMagneticMoment(:) = 0.0_double

      do i = 1, numPotTypes

         ! Initialize variables for the current potential type
         currentNumAlphas = potTypes(i)%numAlphas
         currentCumulAlphaSum = potTypes(i)%cumulAlphaSum

         ! Identify the indices of the summations.
         potTypeInitIndex = currentCumulAlphaSum + 1
         potTypeFinIndex  = currentCumulAlphaSum + currentNumAlphas

         typesMagneticMoment(i) = &
               & sum(generalRho(potTypeInitIndex:potTypeFinIndex,8) &
               & * intgConsts(potTypeInitIndex:potTypeFinIndex))

         totalMagneticMoment = totalMagneticMoment + typesMagneticMoment(i) * &
               & potTypes(i)%multiplicity
      enddo

      ! Record the total magnetic moment for this iteration.
      write (20,*) "Total magnetic moment (au) = ",totalMagneticMoment

      ! Record the magnetic moments for each individual type.
      rewind (13)
      do i = 1, numPotTypes
         write (13,*) i,typesMagneticMoment(i),atomTypes(i)%elementName,&
               & atomTypes(i)%speciesID,atomTypes(i)%typeID
      enddo
   endif

   ! Record the key iteration information to a seperate file.
   if (spin == 1) then
      write (7,fmt="(i5,3f13.8,f17.8)") currIteration, occupiedEnergy, &
            & electronDiff, testableDelta, totalEnergy
   else
      write (7,fmt="(i5,3f13.8,2f17.8)") currIteration, occupiedEnergy, &
            & electronDiff, testableDelta, totalEnergy, totalMagneticMoment
   endif
   call flush (7)


   ! Record the key iteration information to the main output file and a
   !   separate file for easier plotting.
   write (20,*) 'KINETIC ENERGY:        ',kineticEnergy
   write (20,*) 'ELECTROSTATIC ENERGY:  ',elecStatEnergy
   write (20,*) 'EXCH-CORR ENERGY:      ',exchCorrEnergy
   write (20,*) 'TOTAL ENERGY:          ',totalEnergy

   write (14,fmt="(i5,4e18.10)") currIteration, kineticEnergy, elecStatEnergy, &
         & exchCorrEnergy, totalEnergy


   ! At this point we can check for the convergence of the system.  If it has
   !   converged, then we mark it so, otherwise we simply increment the
   !   iteration counter.
   if (testableDelta < convgTest) then
      converged = 1
   else
      currIteration = currIteration + 1
   endif


   ! Deallocate arrays and matrices that are not needed now.
   deallocate (typesMagneticMoment)
   deallocate (generalRho)
   deallocate (chargeDensityTrace)
   deallocate (nucPotTrace)
   deallocate (kineticEnergyTrace)
   deallocate (tempOverlap)
   deallocate (potAlphaOverlap)
!   deallocate (rl0)
!   deallocate (rl1)
!   deallocate (rl2)
!   deallocate (drl1)
!   deallocate (drl2)

   ! Log the date and time we end.
   call timeStampEnd (18)

end subroutine makeSCFPot

! Use the method of D. G. Anderson (J. Assoc. Comp. Mach. 12, 547 (1965)) to
!   compute the potential coefficients for the next SCF iteration. The problem
!   is that we have a non-linear optimization problem but that we have a
!   good(-ish) initial guess for the solution.

! In a generalized secant method one would form a secant hyperplane (in N
!   dimensions, where N is the number of terms in the potential) through two
!   points (possible solutions). In the method used here we will only make a
!   secant hyperline through those same possible solutions. I.e. we will reduce
!   the data in the representation of the coefficients and apply the secant
!   method to that so as to make a simpler overall expression for accelerated
!   convergence with a few practical benefits and some remaining issues. The
!   key idea is that the "most correct" hyperplane approach has a high
!   computational cost (because of the need for large matrix inversion (which
!   I have not yet attempted to apply to this problem so I don't actually know
!   how big of an issue it really it), but a "linearized" form that does not
!   necessarily intesect the subspace defining the solution (as a generalized
!   secant method would) would be easy to compute and it would get close to the
!   right direction. All we have to do then is make an algorithm that will get
!   us as close to the solution subspace as possible and then we should
!   converge to it. (Some issues remain, and they will be discussed later.)

! In the D. G. Anderson paper, the basic iteration is z^l+1 = Gz^l. (Equation
!   4.1). The "l+1" and "l" superscripts are not exponents. Instead they
!   indicate the next iteration and current iteration respectively. The G is an
!   abstract iteration operator, which in our case represents all of the rest of
!   the work in the SCF cycle that is required to take a guess for the
!   potential terms and produce a new guess. The z values represent the
!   potential coefficients from one iteration to the next. The iteration
!   sequence is broken into a coupled pair of two iterative sequences called x
!   and y. The original guess that was actually used is called x^l and the new
!   guess for the next iteration is called y^l. The y^l = Gx^l relationship is
!   defined in equation 4.2. The x and y are each manifest as a 1D array of
!   potential function coefficients. (More explicit details will be givin below
!   about how x, y, and other arrays appear in the program code.)

! Now, there is some difference between what was actually used (x) and the guess
!   for the current iteration (y). This is called a residual and it is defined
!   by equation 4.3: r^l = y^l - x^l.

! At this point the concept of a N-vector inner product is defined as (u,v) = 
!   SUM(u_i * v_i * w_i; i=1..N) where the u and v are arbitrary vectors and
!   the w_i are weighting factors, (Equation 4.4). In our scenario, the
!   inner product takes the form <u|O|v> where u is a row vector, v is a column
!   vector, and O is a positive definite matrix (operator) that is manifest as
!   the potAlphaOverlap matrix. (That is, the operator is the matrix that
!   represents the degree of overlap between Gaussians of different potential
!   terms at different potential sites.) In this sense, the weighting factors
!   are defined by the overlap. The overlap carries meaning because a larger
!   overlap indicates a larger potential in that region (or, for charge, a
!   greater degree of electronic repulsion). However, we could in the future
!   consider weighting factors that are additinally designed to favor
!   convergence for certain terms instead of other terms in the array. The
!   rational for pursuing a weighting factor is given in the paragraph
!   beginning with "Variants of these algorithms..." on page 554. Serious
!   consideration should be given to this line of thinking in the future.

! The simplest algorithm is then given in equations 4.5 though 4.9 and it is
!   called the "extrapolation" method. The extrapolation method is in contrast
!   to the "relaxation" method that is introduced in Equation 3.3 and discussed
!   in the subsequent paragraphs. The idea is that there is some operator that
!   acts on z to produce a zero vector (Fz=0, Equation 3.1) only when the
!   "right" z vector has been found. Or, equivalently, there is an operator G
!   that leaves z unchanged (z=Gz, Equation 3.2) only when the "right" z vector
!   has been found. These expressions can be combined generally to give
!   Equation 3.3, Gz=z-HFz, for a regular homogeneous operator H. That equation
!   will be true for some optimal choice of z. The relaxation approach to
!   finding the optimal z starts with successive substitution of z into
!   z^l+1=Gz^l (Equation 3.4). For a best choice of H (often just some
!   empirically discovered multiplicative constant or array of multiplicative
!   constants) the iterative process will drive toward a solution. Indeed, if
!   the value of H is defined as the inverse of the Jacobian matrix of F with
!   respect to z, then a generalized Newton-Raphson method is created. However,
!   the problem is that the relaxation approach is too slow to be tractable 
!   for our type of problem. Thus, we use the extrapolation method that will
!   accelerate the convergence. (Again, the simplest extrapolation algorithm is
!   given in equations 4.5 through 4.9. We will further modify that form with
!   equations 4.15 through 4.18 for our work.)

! Ultimately, we will want the potential function coefficients for the next
!   iteration to be defined according to equation 4.9 with B^l being equal to
!   the relaxation factor given in the olcao.dat input file and held here as
!   "relaxFactor" (usually 0.2). The purpose of Eqn 4.9 is to mix the currently
!   used potential coefficients x with the terms of the next guess y. Actually,
!   this equation will mix u and v which themselves are constructed from the x
!   and y values of the current and previous iterations. The current u is mixed
!   with a difference between u and v (which is like a residual (r) vector that
!   represents the difference between used and guessed coefficients). The u and
!   v include terms from the current and a number of previous iterations. The u
!   and v are each built using a linearized secant-method approach. The goal of
!   the secant method is to find the root of a function. In this case, the
!   function that we want to find the root of is a type of residual function.
!   It is a real-space difference between the potential that we actually used
!   and the potential that we directly computed (before any blending) for the
!   next iteration.. The potential difference is computed on the
!   exchange-correlation mesh. When that residual function is equal to zero
!   then we have achieved self-consistency. (This is the testableDelta variable
!   that appears earlier in the potential update subroutine.)

! Looking just at u (with v being built in the same way) we consider three
!   points 1, 2, and 3. Point 3 is the new point (i.e. u) while point 2 is the
!   current point (x^l) and point 1 (x^l-1) is the previous point. We have that
!   3 = 2 - f(2)/f'(2) which is basically a Newton's method formula. Then we
!   approximate the derivative f'(2) with f'(2) ~= [f(2)-f(1)] / [2-1].
!   Combining we get the expression 3 = 2 - f(2) * [2-1] / [f(2)-f(1)] which
!   when put into the form of the Anderson paper is Equation 4.5: u^l = x^l +
!   theta^l * (x^l-1 -x^l). Theta will be determined next, but it represents
!   the f(2)/[f(2)-f(1)] part of the method. In actuality, it will not
!   explicitly be computed as that part of the method but will instead be
!   selected to try to zero out the residual. Perhaps a better way to say it is
!   that the [2-1] is a difference in the used coefficients from the previous
!   iteration to the current one (the run) while the 1/[f(2)-f(1)] describes
!   the extent to which the residual has changed (the rise).

! The u vector is the base for the next iteration of coefficients and it is
!   blended with some v. The u vector is built from the current x and a
!   difference between the current and previous values of x. The difference is
!   multiplied by a coefficient (theta) that is specifically selected so as to
!   zero out the linearized residual R. We define R as R=0.5*(v-u,v-u). The
!   expression that yields this condition is dR/dTheta = (r^l+1 - r^l,v - u)=0.
!   The linearized residual is designed as such to represent a measure of the
!   difference between our potential function coefficient guesses and the
!   potential function coefficient values that are actually used. The idea being
!   that the more we are able to make the y values we guess look just like the
!   x values we just used, the closer we are to a converged solution. Thus,
!   taking the derivative of R with respect to the undetermined parameters theta
!   and setting the solution equal to zero should point to theta values that
!   will make the residual smaller. With only 3, 2, and 1 points, the value of
!   theta is given by Equation 4.8.

! Now, it is discussed on pages 554-556 (starting near the bottom of 554) that
!   a generalization of the method to higher degrees (M) can be done and that
!   it is useful for M<~5. Beyond that, the earlier iterations are of limited
!   value when determining new iterations. The usefulness of incorporating more
!   previous iterations at all is that it improves the accuracy of the next
!   step. The down side is that as convergence is reached the iterates become
!   more and more alike and thus the system of equations becomes ill
!   conditioned and possibly singular.

! The implementation of this generalized approach is a straightforward extension
!   of the previous discussion with the following exception. The dR/dTheta = 0
!   requirement will lead to a system of M equations (4.18) instead of a single
!   equation for theta. The system of equations must be solved via linear
!   algebra methods, but this is fairly simple. The only caveat is that as the
!   system moves toward convergence, the system of equations will become more
!   ill conditioned. If the determinant of the matrix formulation of the problem
!   becomes close to zero, then the program will need to fall back to a lower
!   degree solution. Similarly, when the number of iterations performed so far
!   is less than the desired degree, then a lower degree solution will need to
!   be used until the higher degree option is available (i.e. when there are
!   enough previous iterations).

! This subroutine will take a parameter called the feedbackLevel, which is given
!   in the olcao.dat input. The feedbackLevel represents the number M from
!   Anderson's paper. Also note, that usedPotCoeffs is equal to x and
!   guessedPotCoeffs is equal to y. Further, the indices are off by one. That
!   is, the 1 in usedPotCoeffs(:,1,:) refers to x^l, not x^l+1 or x^l-1.
!   Further, if a 2 were there, that would be understood as x^l-1. A 3 index
!   would be for x^l-2 etc.
subroutine blendPotentialsSCF(firstTerm, numTerms, outCoeffs, inGuessedCoeffs,&
      & inUsedCoeffs, totalEnergyRecord)

   use O_Kinds
   use O_Potential, only: spin,potDim,currIteration,feedbackLevel,relaxFactor
   use O_ElectroStatics, only: potAlphaOverlap
   use O_LAPACKDPOSVX

   implicit none

   ! Define the dummy variables that are passed to this function.
   integer, intent(in) :: firstTerm
   integer, intent(in) :: numTerms
   real (kind=double), dimension(:), intent(out) :: outCoeffs
   real (kind=double), dimension(:,:), intent(inout) :: inGuessedCoeffs
   real (kind=double), dimension(:,:), intent(inout) :: inUsedCoeffs
   real (kind=double), dimension(:), intent(inout) :: totalEnergyRecord
!   real (kind=double), dimension(:), intent(inout) :: convergenceRecord

   ! Define the local variables.
   integer :: info
   integer :: i,j,k
   integer :: maxFeedback
   integer, dimension(1) :: worstEnergyIndex
   real (kind=double), allocatable, dimension(:) :: theta ! The x from Ax=B.
   real (kind=double), allocatable, dimension(:) :: tempArray ! An intermediate
         ! array for constructing the final set of potential coefficients.
   real (kind=double), allocatable, dimension(:) :: solutions ! The B from Ax=B
   real (kind=double), allocatable, dimension(:,:) :: solutionsTemp ! The B as
         ! it must be passed to the LAPACK DPOSVX subroutine. (I.e. as a matrix
         ! with multiple right hand sides. In this case though, there is only
         ! one right hand side.) Note also that the dimension for solutionsTemp
         ! is determined by a loop index from maxFeedback down to one. It is
         ! not a constant feedbackLevel or maxFeedback. The idea is that if a
         ! higher degree solution fails, we will shift down to a lower degree
         ! one automatically.
   real (kind=double), allocatable, dimension(:,:) :: matrix ! The A from Ax=B
   real (kind=double), allocatable, dimension(:,:) :: matrixTemp ! The A as it
         ! must be passed to the LAPACK DPOSVX subroutine. It is basically the
         ! same as "matrix" except that the dimension is determined as for
         ! solutionsTemp above and for the same reasons.
   real (kind=double), allocatable, dimension(:,:) :: rl ! Difference between
         ! the actually used and the guessed potential coefficients for each of
         ! the feedbackLevel sets that are retained, plus the current iteration.
         ! Thus, the first index is numTerms and the second is feedbackLevel+1.
   real (kind=double), allocatable, dimension(:,:) :: drl ! Differences
         ! between the first rl and the other rl arrays. The first index is
         ! numTerms and the second is feedbackLevel. (The +1 is not needed here
         ! because we don't need a difference between the current set (first rl)
         ! and itself.)

   ! Determine the maximum amount of feedback that is possible (a function of
   !   the currIteration number) or desired (a function of feedbackLevel). A
   !   certain number of iterations must pass before a particular feebackLevel
   !   can be used. For example is the desired feedbackLevel is 2, then at least
   !   two iterations must complete before we will have 2 prior sets of data to
   !   use for feedback. (Thus we must *be on* iteration #3 before a
   !   feedbackLevel of two can become active.)
   maxFeedback = min(feedbackLevel,currIteration-1)
!write (20,*) "feedback",feedbackLevel,currIteration,maxFeedback

!   ! Find the energetically worst iteration and copy the last iteration over it
!   !   so as to save the iteration in the last index from being shifted out if
!   !   it is better than the worst. Obviously, if the last is also the worst,
!   !   then it will just be copied over itself and then shifted out.
!   ! Don't forget to copy the energy as well.
!   if (currIteration > feedbackLevel) then
!      worstEnergyIndex = maxloc(totalEnergyRecord(1:maxFeedback+2))
!      if (worstEnergyIndex(1) == 1) then
!         do i = 1, potDim
!            inUsedCoeffs(i,1) = sum(inUsedCoeffs(i,2:maxFeedback+1)) &
!                  & / (real(maxFeedback+1,double))
!            inGuessedCoeffs(i,1) = sum(inGuessedCoeffs(i,2:maxFeedback+1)) &
!                  & / (real(maxFeedback+1,double))
!            totalEnergyRecord(1) = sum(totalEnergyRecord(2:maxFeedback+2)) &
!                  & / (real(maxFeedback+1,double))
!         enddo
!      endif
!   endif

   ! Allocate space for the operating data structures.
   allocate (tempArray(numTerms))
   allocate (theta(maxFeedback+1)) ! The first index is special (1-sum(others)).
   allocate (rl(numTerms,maxFeedback+1)) ! +1 to hold the current iteration.
   if (maxFeedback > 0) then
      allocate (drl(numTerms,maxFeedback)) ! No +1 needed.
   endif

   ! Initialize all of the operating parameters that are needed for the
   !   algorithm.

   ! In the event that we don't do any feedback we need to initialize theta
   !   to 1. (If there is no feedback then theta in an array of length 1.)
   theta(:) = 0.0_double
   theta(1) = 1.0_double ! Initialized to 1.

   ! Compute the difference between the guessed and actually used potential
   !   coefficients for the current iteration and every level of feedback. We
   !   can call these deltas.
   do i = 1, maxFeedback+1
      rl(:,i) = (inGuessedCoeffs(:,i) - inUsedCoeffs(:,i))
!write (20,*) "rl i=",i
!write (20,*) rl(:,i)
   enddo

   ! Compute the difference between the deltas from different iterations.
   !   Specifically, get the difference between the current delta (1) and
   !   each other delta.
   if (maxFeedback > 0) then
      do i = 1, maxFeedback
         drl(:,i) = rl(:,1) - rl(:,i+1)
      enddo
   endif

   ! If we are going to bother to do any accelerated convergence then we
   !   proceed. Otherwise we just skip to the step where we form the next
   !   iteration.
   if (maxFeedback /= 0) then

      ! Allocate space to hold the A matrix and B solutions in the largest
      !   case scenario.
      allocate (solutions(maxFeedback))
      allocate (solutionsTemp(maxFeedback,1))
      allocate (matrix(maxFeedback,maxFeedback))
      allocate (matrixTemp(maxFeedback,maxFeedback))

      ! Initialize the matrix and solutions in preparation for accumulation.
      matrix(:,:) = 0.0_double
      solutions(:) = 0.0_double
      tempArray(:) = 0.0_double

      ! Assemble the system of linear equations for Ax=B. For A, we only need
      !   to compute the upper triangle. (The variables we have are A=matrix,
      !   B=solutions.)
!write (20,*) "potAlphaOverlap"
!write (20,*) potAlphaOverlap
      do i = 1, maxFeedback
!write (20,*) "drl i=",i
!write (20,*) drl(:,i)
         do j = 1, i
            do k = firstTerm, firstTerm+numTerms-1
               tempArray(1:numTerms) = potAlphaOverlap(firstTerm:numTerms,k) &
                     & * drl(k-firstTerm+1,i)
               matrix(j,i) = matrix(j,i) + sum(tempArray(1:numTerms) &
                     & * drl(1:numTerms,j))
               if (j == 1) then
                  solutions(i) = solutions(i) + sum(tempArray(1:numTerms) &
                        & * rl(1:numTerms,1))
               endif
            enddo
         enddo
      enddo

      ! Work backwards from the most feedback to the least. With each
      !   iteration we will assemble the A and B from Ax=B. On the first
      !   iteration we will compute the maximum possible feedback. However,
      !   we recognize from the D.G. Anderson paper (554-556) that as the
      !   feedbackLevel is increased there is a greater chance that the
      !   system of linear equations will become ill conditioned (and thus
      !   singular). When we compute the DPOSVX solution we check if there
      !   are any errors
      do i = maxFeedback, 1, -1

         ! Copy the actual matrix A and solutions B into temporary data
         !   structures because the originals will be destroyed within the
         !   dposvx subroutine.
         matrixTemp(:,:) = matrix(:,:)
         solutionsTemp(:,1) = solutions(:)

         ! Solve the system of linear equations Ax=B.
         call solveDPOSVX(i,1,matrixTemp(1:i,1:i),i,solutionsTemp(1:i,1),&
               & info)

         ! Determine if the solution is acceptable. If not, then we cycle to
         !   a lower degree attempt and try again. If we are at the last
         !   chance and we still get an error, then we die and complain.
         if (info /= 0) then
            if (i == 1) then
               write (20,*) "Failed to solve DPOSVX for potential blending."
               write (20,*) "This was the last chance. Stopping."
               write (20,*) "matrix:",matrix(1,1)
               write (20,*) "solutions:",solutions(1)
               stop
            else
               write (20,*) "Failed to solve DPOSVX for potential blending."
               write (20,*) "Trying again with fewer feedback terms."
               write (20,*) "Current feedback term was:",i
               write (20,*) "matrix:",matrix(1:i,1:i)
               write (20,*) "solutions:",solutions(1:i)
            endif
            cycle
         endif

         ! Copy the results into the theta array taking note that the first
         !   index of theta is special. It will be used for
         !   1.0-sum(other_theta_values).
         theta(:) = 0.0_double
         theta(1) = 1.0_double - sum(solutionsTemp(1:i,1))
         theta(2:i+1) = solutionsTemp(1:i,1)

         ! If any solution was successfully obtained that is all we need and
         !   we can skip out of the rest of this loop.
         exit
      enddo

      ! Deallocate the data structures for the linear equation problem.
      deallocate (matrix)
      deallocate (matrixTemp)
      deallocate (solutions)
      deallocate (solutionsTemp)
   endif ! if (maxFeedback /= 0)

   ! Form the potential coefficients for the next iteration.

   ! Use a temp array to collect the contribution from all usedPotCoeffs of
   !   each feedback level and the current iteration. Note that theta(1) is
   !   special and that the other theta values were computed from the system
   !   of linear equations.
   tempArray(:) = 0.0_double
   do i = 1, maxFeedback+1
      tempArray(:) = tempArray(:) + theta(i)*inUsedCoeffs(:,i)
   enddo
   outCoeffs(:) = (1.0_double - relaxFactor) * tempArray(:)

   tempArray(:) = 0.0_double
   do i = 1, maxFeedback+1
      tempArray(:) = tempArray(:) + theta(i)*inGuessedCoeffs(:,i)
   enddo
   outCoeffs(:) = outCoeffs(:) + relaxFactor * tempArray(:)

!   ! Find the energetically worst iteration and copy the last iteration over it
!   !   so as to save the iteration in the last index from being shifted out if
!   !   it is better than the worst. Obviously, if the last is also the worst,
!   !   then it will just be copied over itself and then shifted out. Don't
!   !   forget to copy the energy as well.
!   if (currIteration > feedbackLevel) then
!      inUsedCoeffs(:,worstEnergyIndex(1)) = inUsedCoeffs(:,maxFeedback+1)
!      inGuessedCoeffs(:,worstEnergyIndex(1)) = inGuessedCoeffs(:,maxFeedback+1)
!      totalEnergyRecord(worstEnergyIndex(1)) = totalEnergyRecord(maxFeedback+2)
!   endif

   ! Progressively shift the oldest iterations out by replacing them with an
   !   earlier iteration. (Note that this will do a bit more work than needed
   !   during the first few SCF cycles, but then it will work fine.
   do i = feedbackLevel, 1, -1
      inUsedCoeffs(:,i+1) = inUsedCoeffs(:,i)
      inGuessedCoeffs(:,i+1) = inGuessedCoeffs(:,i)
   enddo
!   do i = feedbackLevel, 1, -1
!      inUsedCoeffs(:,i+1) = inUsedCoeffs(:,i)
!      inGuessedCoeffs(:,i+1) = inGuessedCoeffs(:,i)
!      totalEnergyRecord(i+2) = totalEnergyRecord(i+1)
!   enddo
!   totalEnergyRecord(2) = totalEnergyRecord(1)

!   ! Progressively shift the energetically worse iterations out by replacing
!   !   them with energetically more favorable iterations.
!   do i = 1, feedbackLevel+1
!      if (totalEnergy < totalEnergyRecord(i)) then
!         do j = feedbackLevel, i, -1
!            totalEnergyRecord(j+1) = totalEnergyRecord(j)
!            inUsedCoeffs(:,j+1)    = inUsedCoeffs(:,j)
!            inGuessedCoeffs(:,j+1) = inGuessedCoeffs(:,j)
!         enddo
!         totalEnergyRecord(i) = totalEnergy
!         exit
!      endif
!   enddo

   ! Deallocate space for the operating data structures of the algorithm.
   if (maxFeedback > 0) then
      deallocate (drl)
   endif
   deallocate (rl)
   deallocate (theta)
   deallocate (tempArray)
!   deallocate (totalEnergyPercentage)
!   deallocate (convergencePercentage)
!   deallocate (percentageAverage)
!   deallocate (weightFactor)

end subroutine blendPotentialsSCF


! Just like the blendPotentialsSCF subroutine except that instead of trying to
!   bring the residual to zero we will try to minimize the total energy (TE).
subroutine blendPotentialsTE(firstTerm, numTerms, outCoeffs, inGuessedCoeffs, &
      & inUsedCoeffs, totalEnergyRecord)

   use O_Kinds
   use O_Potential, only: spin, currIteration, feedbackLevel, relaxFactor
   use O_ElectroStatics, only: potAlphaOverlap
   use O_LAPACKDPOSVX

   implicit none

   ! Define the dummy variables that are passed to this function.
   integer, intent(in) :: firstTerm
   integer, intent(in) :: numTerms
   real (kind=double), dimension(:), intent(inout) :: outCoeffs
   real (kind=double), dimension(:,:), intent(inout) :: inGuessedCoeffs
   real (kind=double), dimension(:,:), intent(inout) :: inUsedCoeffs
   real (kind=double), dimension(:), intent(inout) :: totalEnergyRecord

   ! Define the local variables.
   integer :: info
   integer :: i,j,k
   integer :: maxFeedback
   real (kind=double) :: SCF_TE_weight
   real (kind=double), allocatable, dimension(:) :: theta ! The x from Ax=B.
   real (kind=double), allocatable, dimension(:) :: tempArray ! An intermediate
         ! array for constructing the final set of potential coefficients.
   real (kind=double), allocatable, dimension(:) :: tempArray2 ! An intermediate
         ! array for constructing the final set of potential coefficients.
   real (kind=double), allocatable, dimension(:) :: solutions ! The B from Ax=B
   real (kind=double), allocatable, dimension(:,:) :: solutionsTemp ! The B as
         ! it must be passed to the LAPACK DPOSVX subroutine. (I.e. as a matrix
         ! with multiple right hand sides. In this case though, there is only
         ! one right hand side.) Note also that the dimension for solutionsTemp
         ! is determined by a loop index from maxFeedback down to one. It is
         ! not a constant feedbackLevel or maxFeedback. The idea is that if a
         ! higher degree solution fails, we will shift down to a lower degree
         ! one automatically.
   real (kind=double), allocatable, dimension(:,:) :: matrix ! The A from Ax=B
   real (kind=double), allocatable, dimension(:,:) :: matrixTemp ! The A as it
         ! must be passed to the LAPACK DPOSVX subroutine. It is basically the
         ! same as "matrix" except that the dimension is determined as for
         ! solutionsTemp above and for the same reasons.
   real (kind=double), allocatable, dimension(:) :: rlEnergy
!   real (kind=double), allocatable, dimension(:,:) :: rl ! Difference between
         ! the actually used and the guessed potential coefficients for each of
         ! the feedbackLevel sets that are retained, plus the current iteration.
         ! Thus, the first index is numTerms and the second is feedbackLevel+1.
   real (kind=double), allocatable, dimension(:) :: drlEnergy
!   real (kind=double), allocatable, dimension(:,:) :: drl ! Differences
         ! between the first rl and the other rl arrays. The first index is
         ! numTerms and the second is feedbackLevel. (The +1 is not needed here
         ! because we don't need a difference between the current set (first rl)
         ! and itself.)
! CAN BE USED TO BLEND SCF and TE APPROACHES. (Recall that the outCoeff that
!   comes into this subroutine contains the result of the blendPotentialSCF
!   subroutine that ran just before this one. Thus, we can choose to keep some,
!   none, or all of that potential. A value of 1 eliminates all of it and uses
!   only the TE.
   SCF_TE_weight = 1.0_double

   ! Determine the maximum amount of feedback that is possible (a function of
   !   the currIteration number) or desired (a function of feedbackLevel). A
   !   certain number of iterations must pass before a particular
   !   feedbackLevel can be used. For example is the desired feedbackLevel
   !   is 2, then at least two iterations must complete before we will have
   !   2 prior sets of data to use for feedback. (Thus we must *be on*
   !   iteration #3 before a feedbackLevel of two can become active.)
   maxFeedback = min(feedbackLevel,currIteration-2)
write (20,*) "maxFeedback",maxFeedback
call flush (20)

   ! Allocate space for the operating data structures.
   allocate (tempArray(numTerms))
   allocate (tempArray2(numTerms))
!   allocate (rl(numTerms,maxFeedback+1)) ! +1 to hold the current iteration.
!   if (maxFeedback > 0) then
!      allocate (drl(numTerms,maxFeedback)) ! No +1 needed.
!   endif
   allocate (theta(maxFeedback+1)) ! The first index is special (1-sum(others)).
   allocate (rlEnergy(maxFeedback+1))
   allocate (drlEnergy(maxFeedback))

   ! Initialize all of the operating parameters that are needed for the
   !   algorithm.

   ! In the event that we don't do any feedback we need to initialize theta
   !   to 1. (If there is no feedback then theta in an array of length 1.)
   theta(:) = 0.0_double
   theta(1) = 1.0_double ! Initialized to 1.

   ! Compute the difference between the guessed and actually used potential
   !   coefficients for the current iteration and every level of feedback. We
   !   can call these deltas.
!write (20,*) "TERecord(:)",totalEnergyRecord(:)
   do i = 1, maxFeedback+1
      rlEnergy(i) = totalEnergyRecord(i) - totalEnergyRecord(i+1)
!write (20,*) "i TEi+1 TEi",i,totalEnergyRecord(i),totalEnergyRecord(i+1)
   enddo

   ! Compute the difference between the deltas from different iterations.
   !   Specifically, get the difference between the current delta (1) and
   !   each other delta.
   do i = 1, maxFeedback
      drlEnergy(i) = rlEnergy(1) - rlEnergy(i+1)
!write (20,*) "i rlE1 rlEi+1",i,rlEnergy(1),rlEnergy(i+1)
   enddo

   ! Allocate space to hold the A matrix and B solutions in the largest
   !   case scenario.
   allocate (solutions(maxFeedback))
   allocate (solutionsTemp(maxFeedback,1))
   allocate (matrix(maxFeedback,maxFeedback))
   allocate (matrixTemp(maxFeedback,maxFeedback))

   ! Initialize the matrix and solutions in preparation for accumulation.
   matrix(:,:) = 0.0_double
   solutions(:) = 0.0_double
   tempArray(:) = 0.0_double

   ! Assemble the system of linear equations for Ax=B. For A, we only need
   !   to compute the upper triangle. (The variables we have are A=matrix,
   !   B=solutions.)
   do i = 1, maxFeedback
      solutions(i) = drlEnergy(i) * rlEnergy(1)
      do j = 1, i
         matrix(j,i) = drlEnergy(i) * drlEnergy(j)
      enddo
!      do j = firstTerm, firstTerm+numTerms-1
!         tempArray(1:numTerms) = potAlphaOverlap(firstTerm:numTerms,j) &
!               & * drl(j-firstTerm+1,i)
!         solutions(i) = solutions(i) + sum(tempArray(1:numTerms) &
!               & * rl(1:numTerms,1))
!         do k = 1, i
!            matrix(k,i) = matrix(k,i) + sum(tempArray(1:numTerms) &
!                  & * drl(1:numTerms,k))
!         enddo
!      enddo
   enddo

   ! Work backwards from the most feedback to the least. With each
   !   iteration we will assemble the A and B from Ax=B. On the first
   !   iteration we will compute the maximum possible feedback. However,
   !   we recognize from the D.G. Anderson paper (554-556) that as the
   !   feedbackLevel is increased there is a greater chance that the
   !   system of linear equations will become ill conditioned (and thus
   !   singular). When we compute the DPOSVX solution we check if there
   !   are any errors
   do i = maxFeedback, 1, -1

      ! Copy the actual matrix A and solutions B into temporary data
      !   structures because the originals will be destroyed within the
      !   dposvx subroutine.
      matrixTemp(:,:) = matrix(:,:)
      solutionsTemp(:,1) = solutions(:)

      ! Solve the system of linear equations Ax=B.
      call solveDPOSVX(i,1,matrixTemp(1:i,1:i),i,solutionsTemp(1:i,1),&
            & info)

      ! Determine if the solution is acceptable. If not, then we cycle to
      !   a lower degree attempt and try again. If we are at the last
      !   chance and we still get an error, then we die and complain.
      if (info /= 0) then
         if (i == 1) then
            write (20,*) "Failed to solve DPOSVX for potential blending."
            write (20,*) "This was the last chance. Stopping."
            write (20,*) "matrix:",matrix(1,1)
            write (20,*) "solutions:",solutions(1)
            write (20,*) "DPOSVX info = ",info
            stop
         else
            write (20,*) "Failed to solve DPOSVX for potential blending."
            write (20,*) "Trying again with fewer feedback terms."
            write (20,*) "Current feedback term was:",i
            write (20,*) "matrix:",matrix(1:i,1:i)
            write (20,*) "solutions:",solutions(1:i)
            write (20,*) "DPOSVX info = ",info
         endif
         cycle
      endif

      ! Copy the results into the theta array taking note that the first
      !   index of theta is special. It will be used for
      !   1.0-sum(other_theta_values).
      theta(:) = 0.0_double
      theta(1) = 1.0_double - sum(solutionsTemp(1:i,1))
      theta(2:i+1) = solutionsTemp(1:i,1)

      ! If any solution was successfully obtained that is all we need and
      !   we can skip out of the rest of this loop.
      exit
   enddo

   ! Deallocate the data structures for the linear equation problem.
   deallocate (matrix)
   deallocate (matrixTemp)
   deallocate (solutions)
   deallocate (solutionsTemp)

   ! Form the potential coefficients for the next iteration.

   ! Use a temp array to collect the contribution from all usedPotCoeffs of
   !   each feedback level and the current iteration. Note that theta(1) is
   !   special and that the other theta values were computed from the system
   !   of linear equations.
!   outCoeffs(:) = 0.0_double
!   do i = 1, maxFeedback+1
!      outCoeffs(:) = outCoeffs(:) + theta(i)*inUsedCoeffs(:,i)
!   enddo

!   tempArray(:) = 0.0_double
!   do i = 1, maxFeedback+1
!      tempArray(:) = tempArray(:) + theta(i)*inUsedCoeffs(:,i)
!   enddo

   tempArray(:) = 0.0_double
   do i = 1, maxFeedback+1
      tempArray(:) = tempArray(:) + theta(i)*inUsedCoeffs(:,i)
   enddo
   tempArray2(:) = (1.0_double - relaxFactor) * tempArray(:)

   tempArray(:) = 0.0_double
   do i = 1, maxFeedback+1
      tempArray(:) = tempArray(:) + theta(i)*inGuessedCoeffs(:,i)
   enddo
   tempArray2(:) = tempArray2(:) + relaxFactor * tempArray(:)

   ! Average together the contributions from the SCF and TE approaches.
   outCoeffs(:) = (1.0_double - SCF_TE_weight)*outCoeffs(:) &
         & + SCF_TE_weight*tempArray2(:)

   ! Deallocate space for the operating data structures of the algorithm.
   if (currIteration > 2) then
      deallocate (theta)
      deallocate (rlEnergy)
      deallocate (drlEnergy)
      deallocate (tempArray)
      deallocate (tempArray2)
   endif

end subroutine blendPotentialsTE

subroutine shiftPotentials(inGuessedCoeffs,inUsedCoeffs,totalEnergyRecord)

   use O_Kinds
   use O_Potential, only: spin,feedbackLevel

   ! Define the dummy variables that are passed to this function.
   real (kind=double), dimension(:,:,:), intent(inout) :: inGuessedCoeffs
   real (kind=double), dimension(:,:,:), intent(inout) :: inUsedCoeffs
   real (kind=double), dimension(:,:), intent(inout) :: totalEnergyRecord

   ! Define local variables.
   integer :: i, j

   ! Progressively shift the oldest iterations out by replacing them with an
   !   earlier iteration. (Note that this will do a bit more work than needed
   !   during the first few SCF cycles, but then it will work fine.
   do i = 1, spin
      do j = feedbackLevel, 1, -1
         inUsedCoeffs(:,j+1,i) = inUsedCoeffs(:,j,i)
         inGuessedCoeffs(:,j+1,i) = inGuessedCoeffs(:,j,i)
         totalEnergyRecord(j+2,i) = totalEnergyRecord(j+1,i)
      enddo
      totalEnergyRecord(2,i) = totalEnergyRecord(1,i)
   enddo

!   ! Progressively shift the energetically worse iterations out by replacing
!   !   them with energetically more favorable iterations.
!   do i = 1, feedbackLevel+1
!      if (totalEnergy < totalEnergyRecord(i)) then
!         do j = feedbackLevel, i, -1
!            totalEnergyRecord(j+1) = totalEnergyRecord(j)
!            inUsedCoeffs(:,j+1)    = inUsedCoeffs(:,j)
!            inGuessedCoeffs(:,j+1) = inGuessedCoeffs(:,j)
!         enddo
!         totalEnergyRecord(i) = totalEnergy
!         exit
!      endif
!   enddo

end subroutine shiftPotentials

subroutine blendJointPotentials(firstTerm, numJointTerms, outCoeffs,&
!      & inGuessedCoeffs, inUsedCoeffs)
      & inGuessedCoeffs, inUsedCoeffs, totalEnergyRecord)

   use O_Kinds
   use O_Potential, only: spin, currIteration, feedbackLevel, relaxFactor
   use O_ElectroStatics, only: potAlphaOverlap
   use O_LAPACKDPOSVX

   implicit none

   ! Define the dummy variables that are passed to this function.
   integer, intent(in) :: firstTerm
   integer, intent(in) :: numJointTerms
   real (kind=double), dimension(:), intent(out) :: outCoeffs
   real (kind=double), dimension(:,:), intent(inout) :: inGuessedCoeffs
   real (kind=double), dimension(:,:), intent(inout) :: inUsedCoeffs
   real (kind=double), dimension(:), intent(inout) :: totalEnergyRecord

   ! Define the local variables.
   integer :: info
   integer :: i,j,k,l
   integer :: maxFeedback
   real (kind=double), allocatable, dimension(:) :: theta ! The x from Ax=B.
   real (kind=double), allocatable, dimension(:) :: tempArray ! An intermediate
         ! array for constructing the set of linear equations.
   real (kind=double), allocatable, dimension(:) :: tempJointArray ! An
         ! intermediate array for constructing the final set of potential
         ! coefficients.
   real (kind=double), allocatable, dimension(:) :: solutions ! The B from Ax=B
   real (kind=double), allocatable, dimension(:,:) :: solutionsTemp ! The B as
         ! it must be passed to the LAPACK DPOSVX subroutine. (I.e. as a matrix
         ! with multiple right hand sides. In this case though, there is only
         ! one right hand side.) Note also that the dimension for solutionsTemp
         ! is determined by a loop index from maxFeedback down to one. It is
         ! not a constant feedbackLevel or maxFeedback. The idea is that if a
         ! higher degree solution fails, we will shift down to a lower degree
         ! one automatically.
   real (kind=double), allocatable, dimension(:,:) :: matrix ! The A from Ax=B
   real (kind=double), allocatable, dimension(:,:) :: matrixTemp ! The A as it
         ! must be passed to the LAPACK DPOSVX subroutine. It is basically the
         ! same as "matrix" except that the dimension is determined as for
         ! solutionsTemp above and for the same reasons.
   real (kind=double), allocatable, dimension(:,:) :: rl ! Difference between
         ! the actually used and the guessed potential coefficients for each of
         ! the feedbackLevel sets that are retained, plus the current iteration.
         ! Thus, the first index is numJointTerms and the second is
         ! feedbackLevel+1.
   real (kind=double), allocatable, dimension(:,:) :: drl ! Differences
         ! between the first rl and the other rl arrays. The first index is
         ! numJointTerms and the second is feedbackLevel. (The +1 is not needed
         ! here because we don't need a difference between the current set
         ! (first rl) and itself.)

   ! Determine the maximum amount of feedback that is possible (a function of
   !   the currIteration number) or desired (a function of feedbackLevel). A
   !   certain number of iterations must pass before a particular feebackLevel
   !   can be used. For example is the desired feedbackLevel is 2, then at least
   !   two iterations must complete before we will have 2 prior sets of data to
   !   use for feedback. (Thus we must *be on* iteration #3 before a
   !   feedbackLevel of two can become active.)
   maxFeedback = min(feedbackLevel,currIteration-1)

   ! Allocate space for the operating data structures.
   allocate (tempArray(numJointTerms/2))
   allocate (tempJointArray(numJointTerms))
   allocate (theta(maxFeedback+1)) ! The first index is special (1-sum(others)).
   allocate (rl(numJointTerms,maxFeedback+1)) ! +1 to hold current iteration.
   if (maxFeedback > 0) then
      allocate (drl(numJointTerms,maxFeedback)) ! No +1 needed.
   endif

   ! Initialize all of the operating parameters that are needed for the
   !   algorithm.

   ! In the event that we don't do any feedback we need to initialize theta
   !   to 1. (If there is no feedback then theta in an array of length 1.)
   theta(:) = 0.0_double
   theta(1) = 1.0_double ! Initialized to 1.

   ! Compute the difference between the guessed and actually used potential
   !   coefficients for the current iteration and every level of feedback. We
   !   can call these deltas.
   do j = 1, maxFeedback+1
      rl(:,j) = inGuessedCoeffs(:,j) - inUsedCoeffs(:,j)
   enddo

   ! Compute the difference between the deltas from different iterations.
   !   Specifically, get the difference between the current delta (1) and
   !   each other delta.
   if (maxFeedback > 0) then
      do j = 1, maxFeedback
         drl(:,j) = rl(:,1) - rl(:,j+1)
      enddo
   endif

   ! If we are going to bother to do any accelerated convergence then we
   !   proceed. Otherwise we just skip to the step where we form the next
   !   iteration.
   if (maxFeedback /= 0) then

      ! Allocate space to hold the A matrix and B solutions in the largest
      !   case scenario.
      allocate (solutions(maxFeedback))
      allocate (solutionsTemp(maxFeedback,1))
      allocate (matrix(maxFeedback,maxFeedback))
      allocate (matrixTemp(maxFeedback,maxFeedback))

      ! Initialize the matrix and solutions in preparation for accumulation.
      matrix(:,:) = 0.0_double
      solutions(:) = 0.0_double

      ! Assemble the system of linear equations for Ax=B. For A, we only need
      !   to compute the upper triangle.
      do j = 1, maxFeedback
         do k = 1, j
            do l = firstTerm, firstTerm+numJointTerms/2-1
               tempArray(:) = potAlphaOverlap(:,l)*drl(l-firstTerm+1,j)
               matrix(k,j) = matrix(k,j) + sum(tempArray(:) * &
                     & drl(firstTerm:numJointTerms/2,k))
               if (k == 1) then
                  solutions(j) = solutions(j) + sum(tempArray(:) * &
                        & rl(firstTerm:numJointTerms/2,1))
               endif
            enddo
            do l = firstTerm, firstTerm+numJointTerms/2-1
               tempArray(:) = potAlphaOverlap(:,l)*drl(l+numJointTerms/2,j)
               matrix(k,j) = matrix(k,j) + sum(tempArray(:) * &
                     & drl(numJointTerms/2+1:numJointTerms,k))
               if (k == 1) then
                  solutions(j) = solutions(j) + sum(tempArray(:) * &
                        & rl(numJointTerms/2+1:numJointTerms,1))
               endif
            enddo
         enddo
      enddo

      ! Work backwards from the most feedback to the least. With each
      !   iteration we will assemble the A and B from Ax=B. On the first
      !   iteration we will compute the maximum possible feedback. However,
      !   we recognize from the D.G. Anderson paper (554-556) that as the
      !   feedbackLevel is increased there is a greater chance that the
      !   system of linear equations will become ill conditioned (and thus
      !   singular). When we compute the DPOSVX solution we check if there
      !   are any errors.
      do j = maxFeedback, 1, -1

         ! Copy the actual matrix A and solutions B into temporary data
         !   structures because the originals will be destroyed within the
         !   dposvx subroutine.
         matrixTemp(:,:) = matrix(:,:)
         solutionsTemp(:,1) = solutions(:)
!write (20,*) "matrixTemp"
!write (20,*) matrixTemp(1:j,1:j)
!write (20,*) "solutionsTemp"
!write (20,*) solutionsTemp(1:j,1)

         ! Solve the system of linear equations Ax=B.
         call solveDPOSVX(j,1,matrixTemp(1:j,1:j),j,solutionsTemp(1:j,1),&
               & info)

         ! Determine if the solution is acceptable. If not, then we cycle to
         !   a lower degree attempt and try again. If we are at the last
         !   chance and we still get an error, then we die and complain.
         if (info /= 0) then
            if (j == 1) then
               write (20,*) "Failed to solve DPOSVX for potential blending."
               write (20,*) "This was the last chance. Stopping."
               write (20,*) "matrix:",matrix(1,1)
               write (20,*) "solutions:",solutions(1)
               stop
            else
               write (20,*) "Failed to solve DPOSVX for potential blending."
               write (20,*) "Trying again with fewer feedback terms."
               write (20,*) "Current feedback term was:",j
               write (20,*) "matrix:",matrix(1:j,1:j)
               write (20,*) "solutions:",solutions(1:j)
            endif
            cycle
         endif

         ! Copy the results into the theta array taking note that the first
         !   index of theta is special. It will be used for
         !   1.0-sum(other_theta_values).
         theta(:) = 0.0_double
         theta(1) = 1.0_double - sum(solutionsTemp(1:j,1))
         theta(2:j+1) = solutionsTemp(1:j,1)

         ! If any solution was successfully obtained that is all we need and
         !   we can skip out of the rest of this loop.
         exit
      enddo

      ! Deallocate the data structures for the linear equation problem.
      deallocate (matrix)
      deallocate (matrixTemp)
      deallocate (solutions)
      deallocate (solutionsTemp)
   endif ! if (maxFeedback /= 0)

   ! Form the potential coefficients for the next iteration.

   ! Use a temp array to collect the contribution from all usedPotCoeffs of
   !   each feedback level and the current iteration. Note that theta(1) is
   !   special and that the other theta values were computed from the system
   !   of linear equations.
   tempJointArray(:) = 0.0_double
   do j = 1, maxFeedback+1
      tempJointArray(:) = tempJointArray(:) + theta(j)*inUsedCoeffs(:,j)
   enddo
   outCoeffs(:) = (1.0_double - relaxFactor) * tempJointArray(:)

   tempJointArray(:) = 0.0_double
   do j = 1, maxFeedback+1
      tempJointArray(:) = tempJointArray(:) + theta(j)*inGuessedCoeffs(:,j)
   enddo
   outCoeffs(:) = outCoeffs(:) + relaxFactor * tempJointArray(:)

   ! Progressively shift the oldest iterations out by replacing them with an
   !   earlier (lower middle index) iteration.
   do j = feedbackLevel, 1, -1
      inUsedCoeffs(:,j+1) = inUsedCoeffs(:,j)
      inGuessedCoeffs(:,j+1) = inGuessedCoeffs(:,j)
      totalEnergyRecord(j+1) = totalEnergyRecord(j)
   enddo

   ! Deallocate space for the operating data structures of the algorithm.
   if (maxFeedback > 0) then
      deallocate (drl)
   endif
   deallocate (rl)
   deallocate (theta)
   deallocate (tempArray)
   deallocate (tempJointArray)

end subroutine blendJointPotentials



subroutine wignerXC (rho,answer)

   ! Include kind definitions
   use O_Kinds

   implicit none

   ! Define dummy variable passed to this function and the return value.
   real (kind=double) :: rho
   real (kind=double), dimension (2) :: answer

   ! Define local variables used in this function
   real (kind=double) :: rho13
   real (kind=double) :: rho13Factor

! 0.738 = 3/4 * 2/(pi*alpha*rs) * 0.5  (Hartree).
! 0.984 = 2/(pi*alpha*rs) * 0.5  (Hartree).

   ! Assign the local variables used in this routine
!   rho13 = rho**Wtemp(3) ! rho^(1/3)
   rho13 = rho**0.3333333333333333333_double ! rho^(1/3)
   rho13Factor = (1.0_double + 12.57_double * rho13)

   ! Determine the XC potential(1) and energy(2) via the Wigner method.
   answer(1) = -rho13 * (0.984_double + &
         & (0.943656_double + 8.8963_double * rho13) / (rho13Factor)**2)
   answer(2) = -0.738_double * rho13 * (1.0_double + 0.959_double / rho13Factor)

!   answer(1) = -rho13 * (Wtemp(1) + &
!         & (0.943656_double + 8.8963_double * rho13) / (rho13Factor)**2)
!   answer(2) = Wtemp(2) * rho13 * (1.0_double + 0.959_double / rho13Factor)

end subroutine wignerXC

subroutine wignerXCEnergy (rho,answer)

   ! Include kind definitions
   use O_Kinds

   implicit none

   ! Define dummy variable passed to this function and the return value.
   real (kind=double) :: rho ! Core rho.
   real (kind=double) :: answer

   ! Define local variables used in this function
   real (kind=double) :: rho13 ! rho^(1/3)

! 0.738 = 3/4 * 2/(pi*alpha*rs) * 0.5  (Hartree).

   ! Determine the XC energy via the Wigner method.
   rho13 = rho**0.3333333333333333333_double ! rho^(1/3)
   answer = -0.738_double * rho13 * (1.0_double + 0.959_double / &
         & (1.0_double + 12.57_double * rho13))

!   rho13 = rho**Wtemp(3)
!   answer = Wtemp(2) * rho13 * (1.0_double + 0.959_double / &
!         & (1.0_double + 12.57_double * rho13))

end subroutine wignerXCEnergy



subroutine ceperleyAlderXC (rho,answer)

   ! Include kind definitions
   use O_Kinds
   use O_Constants, only: pi, smallThresh

   implicit none

   ! Define dummy variable passed to this function and the return value.
   real (kind=double) :: rho
   real (kind=double),dimension (2) :: answer

   ! Define local variables used in this function
   real (kind=double) :: alpha  ! alpha sub 0
   real (kind=double) :: rs     ! r sub s
   real (kind=double) :: g      ! gamma
   real (kind=double) :: B1
   real (kind=double) :: B2
   real (kind=double) :: ec     ! e sub c
   real (kind=double) :: ex     ! e sub x
   real (kind=double) :: ex0    ! e sub x sup 0
   real (kind=double) :: sqrtrs ! sqrt(rs)
   real (kind=double) :: lnrs   ! ln(rs)
   real (kind=double) :: AF     ! A sub P (Paramagnetic) was AP
   real (kind=double) :: BF     ! B sub P (Paramagnetic) was BP
   real (kind=double) :: CF     ! C sub P (Paramagnetic) was CP
   real (kind=double) :: DF     ! D sub P (Paramagnetic) was DP

   if (rho < smallThresh) then
      answer(:) = 0.0_double
      return
   endif

   alpha = (4.0_double/(9.0_double*pi))**0.33333333
   rs = (3.0_double/(4.0_double*pi*rho))**0.33333333
   ex0 = 3.0_double/(2.0_double * pi * alpha)

   if (rs >= 1.0) then
      g  = -0.1423
      B1 = 1.0529
      B2 = 0.3334
      ! answer(2) = ex + ec.
      sqrtrs = sqrt(rs)
      ex = -ex0/rs
      ec = g / (1 + B1*sqrtrs + B2*rs)
      answer(2) = (ex + ec)

      ! answer(1) = ux + uc
      answer(1) = (4.0/3.0*ex + ec * ec/g * &
            & (1 + 7.0/6.0*B1*sqrtrs + 4.0/3.0*B2*rs))
   else
      AF =  0.0311
      BF = -0.0480
      CF =  0.0020
      DF = -0.0116

      lnrs = log(rs) ! log is the natural log (ln) in fortran.

      ! answer(2) = ex + ec.
      ex = -ex0/rs
      answer(2) = (ex + AF*lnrs + BF + CF*rs*lnrs + DF*rs)

      ! answer(1) = ux + uc
      answer(1) = (4.0/3.0*ex + AF*lnrs + &
            & (BF - 1.0/3.0*AF) + 2.0/3.0*CF*rs*lnrs + &
            & 1.0/3.0*(2*DF - CF)*rs)
   endif

   answer(:) = answer(:) * 0.5

end subroutine ceperleyAlderXC

subroutine ceperleyAlderXCEnergy (rho,answer)

   ! Include kind definitions
   use O_Kinds
   use O_Constants, only: pi, smallThresh

   implicit none

   ! Define dummy variable passed to this function and the return value.
   real (kind=double) :: rho ! Core Rho
   real (kind=double) :: answer

   ! Define local variables used in this function
   real (kind=double) :: alpha  ! alpha sub 0
   real (kind=double) :: rs     ! r sub s
   real (kind=double) :: g      ! gamma
   real (kind=double) :: B1
   real (kind=double) :: B2
   real (kind=double) :: ec     ! e sub c
   real (kind=double) :: ex     ! e sub x
   real (kind=double) :: ex0    ! e sub x sup 0
   real (kind=double) :: sqrtrs ! sqrt(rs)
   real (kind=double) :: lnrs   ! ln(rs)
   real (kind=double) :: AF     ! A sub P (Paramagnetic) was AP
   real (kind=double) :: BF     ! B sub P (Paramagnetic) was BP
   real (kind=double) :: CF     ! C sub P (Paramagnetic) was CP
   real (kind=double) :: DF     ! D sub P (Paramagnetic) was DP

   if (rho < smallThresh) then
      answer = 0.0_double
      return
   endif

   alpha = (4.0_double/(9.0_double*pi))**0.33333333
   rs = (3.0_double/(4.0_double*pi*rho))**0.33333333
   ex0 = 3.0_double/(2.0_double * pi * alpha)

   if (rs >= 1.0) then
      g  = -0.1423
      B1 = 1.0529
      B2 = 0.3334
      ! answer = ex + ec.
      sqrtrs = sqrt(rs)
      ex = -ex0/rs
      ec = g / (1 + B1*sqrtrs + B2*rs)
      answer = (ex + ec)
   else
      AF =  0.0311
      BF = -0.0480
      CF =  0.0020
      DF = -0.0116

      lnrs = log(rs) ! log is the natural log (ln) in fortran.

      ! answer = ex + ec.
      ex = -ex0/rs
      answer = (ex + AF*lnrs + BF + CF*rs*lnrs + DF*rs)
   endif

   answer = answer * 0.5

end subroutine ceperleyAlderXCEnergy

subroutine hedinLundqvistXC (rho,answer)

   ! Include kind definitions and the HL constants.
   use O_Kinds
   use O_Constants, only: pi, smallThresh

   implicit none

   ! Define dummy variable passed to this function and the return value.
   real (kind=double) :: rho ! Charge density at a specific point.
   real (kind=double), dimension (2) :: answer

   ! Define local variables used in this function
!   real (kind=double) :: invRho13 ! (1/rho)^(1/3)

   real (kind=double) :: alpha
   real (kind=double) :: rs
   real (kind=double) :: x

   if (rho < smallThresh) then
      answer(:) = 0.0_double
      return
   endif

   alpha = (4.0_double/(9.0_double*pi))**0.33333333
   rs = (3.0_double/(4.0_double*pi*rho))**0.33333333
   x = rs/21.0_double

   answer(1) = (-2.0_double/(pi * alpha * rs) - 0.045_double * &
         & log(1.0_double + 1.0_double/x))

   answer(2) = (-0.75_double * 2.0_double / (pi * alpha * rs) - &
         & 0.045_double * ((1+x**3.0_double) * log(1.0_double+1.0_double/x) + &
         & x/2.0_double - x**2.0_double - 0.33333333))

   answer(:) = answer(:) * 0.5

!   invRho13 = (1.0_double/rho)**(1.0_double/3.0_double)

   ! answer(1) = ux + uc.
!   answer(1) = (HLtemp(6)/invRho13 - 0.045_double*log(1+HLtemp(2)/invRho13))* &
!         & 0.5_double ! Hartree

   ! answer(2) = ex + ec.
!   answer(2) = (-0.045_double * ((1+HLtemp(1)/rho)*log(1+HLtemp(2)/invRho13) + &
!         & HLtemp(3)*invRho13 - HLtemp(4)*invRho13*invRho13 - &
!         & 1.0_double/3.0_double) - HLtemp(5)/invRho13) * 0.5_double ! Hartree

end subroutine hedinLundqvistXC

subroutine hedinLundqvistXCEnergy (rho,answer)

   ! Include kind definitions and the HL constants.
   use O_Kinds
   use O_Constants, only: pi, smallThresh

   implicit none

   ! Define dummy variable passed to this function and the return value.
   real (kind=double) :: rho ! Charge density at a specific core point.
   real (kind=double) :: answer

   ! Define local variables used in this function
!   real (kind=double) :: invRho13 ! (1/rho)^(1/3)

   real (kind=double) :: alpha
   real (kind=double) :: rs
   real (kind=double) :: x

   if (rho < smallThresh) then
      answer = 0.0_double
      return
   endif

   alpha = (4.0_double/(9.0_double*pi))**0.33333333
   rs = (3.0_double/(4.0_double*pi*rho))**0.33333333
   x = rs/21.0_double

   answer = (-0.75_double * 2.0_double / (pi * alpha * rs) - &
         & 0.045_double * ((1+x**3.0_double) * log(1.0_double+1.0_double/x) + &
         & x/2.0_double - x**2.0_double - 0.33333333))

   answer = answer * 0.5

!   invRho13 = (1.0_double/rho)**(1.0_double/3.0_double)

!   answer = (-0.045_double * ((1+HLtemp(1)/rho) * log(1+HLtemp(2)/invRho13) + &
!         & HLtemp(3)*invRho13 - HLtemp(4)*invRho13*invRho13 - &
!         & 1.0_double/3.0_double) - HLtemp(5)/invRho13) * 0.5_double ! Hartree

end subroutine hedinLundqvistXCEnergy

subroutine vonBarthHedin (totalRho,spinDiffRho,coreRho,answer)

   ! Include kind definitions and the spin polarized vBH constants
   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! Define dummy variables passed to this function and the return value.
   real (kind=double) :: totalRho,spinDiffRho,coreRho
   real (kind=double), dimension (4) :: answer

   ! Define local variables used in this function
   real (kind=double) :: x ! Rho(+) / RhoTotal
 !  real (kind=double) :: rho13
   real (kind=double) :: epsxP
   real (kind=double) :: epsxF
   real (kind=double) :: muxP
   real (kind=double) :: fofx
   real (kind=double) :: epsx
!   real (kind=double) :: rsrP
!   real (kind=double) :: rsrF
   real (kind=double) :: mucP
   real (kind=double) :: mucF
   real (kind=double) :: epscP
   real (kind=double) :: epscF
 !  real (kind=double) :: epscDiff
   real (kind=double) :: epsc
 !  real (kind=double) :: coefficient
!   real (kind=double) :: factor

   real (kind=double) :: alpha
   real (kind=double) :: rs
   real (kind=double) :: gamma
   real (kind=double) :: a
   real (kind=double) :: z
   real (kind=double) :: two13
   real (kind=double) :: rP
   real (kind=double) :: rF
   real (kind=double) :: cP
   real (kind=double) :: cF
   real (kind=double) :: tauc

   ! Initialize the answer
   answer(:) = 0.0_double

   ! Initialize some variables
   alpha = (4.0_double/(9.0_double*pi))**0.33333333_double
   rs = (3.0_double/(4.0_double*pi*totalRho))**0.33333333_double
   a = 2.0_double**(-1.0_double/3.0_double)
   gamma = 4.0_double/3.0_double * a / (1.0_double - a)
   two13 = 2.0**(0.3333333333_double)

   rP = 21.0_double
   rF = 2.0_double**(4.0_double/3.0_double)*rP
   cP = 0.045_double
   cF = cP/2.0_double

   ! Abort if the incoming total charge is sufficiently large, then compute the
   !   exchange correlation for it.
   if (totalRho > 1.0d-25) then


      ! Reminder of the origin of the vBHtemp values.
      !   vBHtemp(1)  = 1/3  ! Exponent value
      !   vBHtemp(2)  = -3/(2pi (4/9 pi)^(1/3) ((3/(4 pi))^(1/3)))*0.5 (Hartree)
      !   vBHtemp(3)  = 4/3 vBHtemp(2)
      !   vBHtemp(4)  = 1/(1-a); a = 2**(-1/3)
      !   vBHtemp(5)  = 4/3  ! Exponent value
      !   vBHtemp(6)  = a = 2**(-1/3)
      !   vBHtemp(7)  = gamma = 4/3  *  a/(1-a)
      !   vBHtemp(8)  = ((3/(4 pi))**(1/3)) / 21
      !   vBHtemp(9)  = ((3/(4 pi))**(1/3)) / (2^(4/3) * 21)
      !   vBHtemp(10) = 1/(1-a) * (2*x^4/3 - a); x=1/2

      ! Prepare the necessary parameters for the calculation of epsilon sub x.
      !   Note that some of the intermediate values calculated here will be
      !   reused later.

      ! Compute the value of x = (Rho(+) / Rho)
      x = (totalRho + spinDiffRho) / (2.0_double * totalRho)
      if (x < 0.0_double) then
         x = 0.0_double
      elseif (x > 1.0_double) then
         x = 1.0_double
      endif

      epsxP = -3.0_double/(2.0_double * pi * alpha) / rs

      epsxF = two13 * epsxP

      muxP = 4.0_double/3.0_double * epsxP

!      muxF = two13 * muxP

      fofx = 1.0/(1.0 - a) * (x**(4.0/3.0) + &
            & (1.0 - x)**(4.0/3.0) - a)

      epsx = epsxP + (epsxF - epsxP) * fofx

      z = rs/rP

      epscP = -cP*((1 + (z)**3.0) * log(1.0+1.0/z) + z/2.0 - &
            & (z)**2.0 - 1.0/3.0)

      z = rs/rF

      epscF = -cF*((1 + (z)**3.0) * log(1.0+1.0/z) + z/2.0 - &
            & (z)**2.0 - 1.0/3.0)

      epsc = epscP + (epscF - epscP)*fofx

      answer(3) = epsx + epsc

      mucP = -cP * log(1.0 + rP/rs)

      mucF = -cF * log(1.0 + rF/rs)

      tauc = mucF - mucP - 4.0/3.0 * (epscF - epscP)

      ! Compute the value of mu sub xc (+).
      answer(1) = (muxP + gamma*(epscF-epscP)) * (2.0*x)**(1.0/3.0) + mucP - &
            & gamma*(epscF - epscP) + tauc * fofx

      ! Compute the value of mu sub xc (-).
      answer(2) = (muxP + gamma*(epscF-epscP)) * (2.0*(1.0-x))**(1.0/3.0)+mucP-&
            & gamma*(epscF - epscP) + tauc * fofx
   endif



   if (coreRho <= 1.0E-25) then
      answer(:) = answer(:) * 0.5
      return
   endif


   rs = (3.0_double/(4.0_double*pi*coreRho))**0.33333333
   x = 0.5_double

   epsxP = -3.0/(2.0 * pi * alpha) / rs

   epsxF = two13 * epsxP

   fofx = 1.0/(1.0 - a) * (x**(4.0/3.0) + &
         & (1.0 - x)**(4.0/3.0) - a)

   epsx = epsxP + (epsxF - epsxP) * fofx

   z = rs/rP

   epscP = -cP*((1 + (z)**3.0) * log(1.0+1.0/z) + z/2.0 - &
         & (z)**2.0 - 1.0/3.0)

   z = rs/rF

   epscF = -cF*((1 + (z)**3.0) * log(1.0+1.0/z) + z/2.0 - &
         & (z)**2.0 - 1.0/3.0)

   epsc = epscP + (epscF - epscP)*fofx


   answer(4) = epsx + epsc

   answer(:) = answer(:) * 0.5

   return

!   ! Compute Total Rho to the 1/3 power
!   rho13 = totalRho**(vBHtemp(1))
!
!   ! Compute value of epsilon sub x sup P
!   epsxP = vBHtemp(2) * rho13
!
!   ! Compute the value of mu sub x sup P
!   muxP = vBHtemp(3) * rho13
!
!   ! Compute the value of fofx.  [Read f of x, or f(x)]
!   fofx = vBHtemp(4)*(x**vBHtemp(5) + (1-x)**vBHtemp(5) - vBHtemp(6))
!
!   ! Compute the value of epsilon sub x for the total charge
!   epsx = epsxP + muxP * fofx / vBHtemp(7)
!
!
!   ! Prepare the necessary parameters for the calculation of epsilon sub c.
!   !   Note again that some of the intermediate values calculated here will be
!   !   reused later.
!
!   ! Compute the value of (r sub s) / r sup P
!   rsrP = vBHtemp(8) / rho13
!
!   ! Compute the value of (r sub s) / r sup F
!   rsrF = vBHtemp(9) / rho13
!
!   ! Compute the value of mu sub c sup P.
!   !   Using 0.045/2 (0.0225) will convert to Hartree.
!   mucF = -0.0225_double * log (1.0_double + 1.0_double/rsrP)
!
!   ! Compute the value of mu sub c sup F.
!   !   Using 0.0225/2 (0.01125) will convert to Hartree.
!   mucF = -0.01125_double * log (1.0_double + 1.0_double/rsrF)
!
!   ! Compute the value of epsilon sub c sup P
!   epscP = mucF*(1+rsrP**3.0_double) - 0.0225_double * &
!         & (rsrP/2.0_double - rsrP**2.0_double - vBHtemp(1)) ! Hartree
!
!   ! Compute the value of epsilon sub c sup F
!   epscF = mucF*(1+rsrF**3.0_double) - 0.01125_double * &
!         & (rsrF/2.0_double - rsrF**2.0_double - vBHtemp(1)) ! Hartree
!
!   ! Compute the difference between epscF and epscP with the gamma factor
!   epscDiff = vBHtemp(7) * (epscF - epscP)
!
!   ! Compute the value of epsilon sub c for the total charge
!   epsc = epscP + epscDiff * fofx / vBHtemp(7)
!
!   ! Compute the value of epsilon sub xc for the total charge
!   answer(3) = epsx + epsc
!
!
!   ! Prepare the necessary parameters for the calculation of v sub xc.  The
!   !   computation of the (+) and (-) values is very similar and so the
!   !   parameters for one will be reused again for the second.
!
!   ! Compute the coefficient for the (2x)**1/3 or (2(1-x))**1/3 terms.
!   coefficient = muxP + epscDiff
!
!   ! Compute the additive factor for the potential terms.
!   factor = mucF - epscDiff + &
!         & (mucF - mucF - epscDiff / (vBHtemp(4) * vBHtemp(6))) * fofx
!
!   ! Compute the value of mu sub xc (+).
!   answer(1) = coefficient * ((2.0_double*x)**vBHtemp(1)) + factor
!
!   ! Compute the value of mu sub xc (-).
!   answer(2) = coefficient * ((2.0_double*(1.0_double-x))**vBHtemp(1)) + factor
!
!
!
!   ! Check to see if the core charge contribution is negligable.
!   if (coreRho < 1.0d-25) then
!      return
!   endif
!
!   ! Prepare the necessary parameters for the calculation of epsilon sub c
!   !   for the core charge.  Note very importantly that in this case it is
!   !   assumed that the core charge has no spin polarization so that the spin
!   !   up charge is exactly equal to the spin down charge.  In such a case
!   !   the value of x = 1/2.
!
!
!   ! Compute Core Rho to the 1/3 power
!   rho13 = coreRho**(vBHtemp(1))
!
!   ! Compute value of epsilon sub x sup P
!   epsxP = vBHtemp(2) * rho13
!
!   ! Compute the value of mu sub x sup P
!   muxP = vBHtemp(3) * rho13
!
!   ! Compute the value of fofx.  [Read f of x, or f(x)].  In this case the
!   !   value of x=1/2 and so this is computed entirely ahead of time as
!   !   vBHtemp(10).
!   fofx = vBHtemp(4)*(2*x**vBHtemp(5) - vBHtemp(6))
!
!   ! Compute the value of epsilon sub x for the core charge
!   epsx = epsxP + muxP * vBHtemp(10) / vBHtemp(7)
!
!
!
!   ! Prepare the necessary parameters for the calculation of epsilon sub c.
!   !   Note again that some of the intermediate values calculated here will be
!   !   reused later.
!
!   ! Compute the value of (r sub s) / r sup P
!   rsrP = vBHtemp(8) / rho13
!
!   ! Compute the value of (r sub s) / r sup F
!   rsrF = vBHtemp(9) / rho13
!
!   ! Compute the value of mu sub c sup P
!   !   Using 0.045/2 (0.0225) will convert to Hartree.
!   mucF = -0.0225_double * log (1 + 1/rsrP)
!
!   ! Compute the value of mu sub c sup F
!   !   Using 0.0225/2 (0.01125) will convert to Hartree.
!   mucF = -0.01125_double * log (1 + 1/rsrF)
!
!   ! Compute the value of epsilon sub c sup P
!   epscP = mucF*(1+rsrP**3) - 0.0225_double * &
!         & (rsrP/2 - rsrP**2 - vBHtemp(1)) ! Hartree as above.
!
!   ! Compute the value of epsilon sub c sup P
!   epscF = mucF*(1+rsrF**3) - 0.01125_double * &
!         & (rsrF/2 - rsrF**2 - vBHtemp(1)) ! Hartree as above.
!
!   ! Compute the difference between epscF and epscP with the gamma factor
!   epscDiff = vBHtemp(7) * (epscF - epscP)
!
!   ! Compute the value of epsilon sub c for the total charge
!   epsc = epscP + epscDiff * fofx / vBHtemp(7)
!
!   ! Compute the value of epsilon sub xc for the core charge
!   answer(4) = epsx + epsc
!

end subroutine vonBarthHedin


subroutine ceperleyAlderSP (totalRho,spinDiffRho,coreRho,answer)

   ! Include kind definitions and the spin polarized vBH constants
   use O_Kinds
   use O_Constants, only: pi, smallThresh

   implicit none

   ! Define dummy variables passed to this function and the return value.
   real (kind=double), intent(in) :: totalRho,spinDiffRho,coreRho
   real (kind=double), intent(out), dimension (4) :: answer

   ! Define local varibles
   real (kind=double) :: zeta
   real (kind=double) :: fofzeta
   real (kind=double) :: dfofzeta
   real (kind=double) :: rs
   real (kind=double) :: g
   real (kind=double) :: sqrtrs ! sqrt(rs)
   real (kind=double) :: lnrs   ! ln(rs)
   real (kind=double) :: alpha
   real (kind=double) :: ecP
   real (kind=double) :: ecF
   real (kind=double) :: ec
   real (kind=double) :: exPara
   real (kind=double) :: exF ! Dont want to use exP since it is reserved.
   real (kind=double) :: ex
   real (kind=double) :: ex0
   real (kind=double) :: fourThirds  ! 4/3
   real (kind=double) :: twoThirds   ! 2/3
   real (kind=double) :: oneThird    ! 1/3
   real (kind=double) :: B1P
   real (kind=double) :: B1F
   real (kind=double) :: B2P
   real (kind=double) :: B2F
   real (kind=double) :: AP     ! A sub P
   real (kind=double) :: BP     ! B sub P
   real (kind=double) :: CP     ! C sub P
   real (kind=double) :: DP     ! D sub P
   real (kind=double) :: AF     ! A sub F
   real (kind=double) :: BF     ! B sub F
   real (kind=double) :: CF     ! C sub F
   real (kind=double) :: DF     ! D sub F
   real (kind=double) :: mucF   ! Ferromagnetic
   real (kind=double) :: mucP   ! Paramagnetic
   real (kind=double) :: muxF   ! Ferromagnetic
   real (kind=double) :: muxP   ! Paramagnetic
   real (kind=double) :: mucUp
   real (kind=double) :: mucDn
   real (kind=double) :: muxUp
   real (kind=double) :: muxDn

   answer(:) = 0.0_double
   fourThirds = 4.0_double/3.0_double
   twoThirds  = 2.0_double/3.0_double
   oneThird   = 1.0_double/3.0_double


   if (totalRho >= smallThresh) then

      alpha = (4.0_double/(9.0_double*pi))**0.33333333_double
      rs = (3.0_double/(4.0_double*pi*totalRho))**0.33333333_double
      zeta = spinDiffRho/totalRho
      ex0 = 3.0_double/(2.0_double * pi * alpha)

      if (rs >= 1.0_double) then
         g   = -0.1423_double
         B1P =  1.0529_double
         B1F =  1.3981_double
         B2P =  0.3334_double
         B2F =  0.2611_double
         sqrtrs = sqrt(rs)
         fofzeta = ((1.0_double+zeta)**fourThirds &
               & + (1.0_double-zeta)**fourThirds - 2.0_double) &
               & / (2.0_double**fourThirds - 2.0_double)
         dfofzeta = (fourThirds * (1.0_double+zeta)**oneThird - &
                  &  fourThirds * (1.0_double-zeta)**oneThird)/ &
                  &  (2.0_double**fourThirds - 2.0_double)
         ! answer(3) = ex + ec.
         ecP = g / (1.0_double + B1P*sqrtrs + B2P*rs)
         ecF = g / (1.0_double + B1F*sqrtrs + B2F*rs)
         ec  = ecP + fofzeta*(ecF-ecP)
         exPara = -ex0/rs
         exF = exPara * 2.0_double**oneThird
         ex  = exPara + fofzeta*(exF-exPara)
         answer(3) = ex + ec

         ! answer(1,2) = ux + uc
         mucP = ecP * ecP / g * (1.0_double + 7.0_double/6.0_double*B1P*sqrtrs &
               & + 4.0_double/3.0_double*B2P*rs)
         mucF = ecF * ecF / g * (1.0_double + 7.0_double/6.0_double*B1F*sqrtrs &
               & + 4.0_double/3.0_double*B2F*rs)
         muxP  = fourThirds * exPara 
         muxF  = fourThirds * exF
         mucUp = mucP + fofzeta*(mucF-mucP) + (ecF-ecP) &
            & *( 1.0_double-zeta)*dfofzeta
         mucDn = mucP + fofzeta*(mucF-mucP) + (ecF-ecP) &
            & *(-1.0_double-zeta)*dfofzeta
         muxUp = muxP + fofzeta*(muxF-muxP) + (exF-exPara) &
            & *( 1.0_double-zeta)*dfofzeta
         muxDn = muxP + fofzeta*(muxF-muxP) + (exF-exPara) &
            & *(-1.0_double-zeta)*dfofzeta
         answer(1) = (muxUp + mucUp)
         answer(2) = (muxDn + mucDn)
      else
         AP =  0.03110_double
         AF =  0.01555_double
         BP = -0.04800_double
         BF = -0.02690_double
         CP =  0.00200_double
         CF =  0.00070_double
         DP = -0.01160_double
         DF = -0.00480_double

         lnrs = log(rs) ! log is the natural log (ln) in fortran.
         fofzeta = ((1.0_double+zeta)**fourThirds &
               & + (1.0_double-zeta)**fourThirds - 2.0_double) &
               & / (2.0_double**fourThirds - 2.0_double)
         dfofzeta = (fourThirds * (1.0_double+zeta)**oneThird - &
                  &  fourThirds * (1.0_double-zeta)**oneThird)/ &
                  &  (2.0_double**fourThirds - 2.0_double)

         ! answer(3) = ex + ec.
         ecP = AP*lnrs + BP + CP*rs*lnrs + DP*rs
         ecF = AF*lnrs + BF + CF*rs*lnrs + DF*rs
         ec = ecP + fofzeta*(ecF-ecP)
         exPara = -ex0/rs
         exF = exPara * 2.0_double**oneThird
         ex = exPara + fofzeta*(exF-exPara)
         answer(3) = ex + ec

         ! answer(1,2) = ux + uc
         mucP  = AP*lnrs + (BP - oneThird*AP) + twoThirds*CP*rs*lnrs + &
               & oneThird*(2.0_double*DP - CP)*rs
         mucF  = AF*lnrs + (BF - oneThird*AF) + twoThirds*CF*rs*lnrs + &
               & oneThird*(2.0_double*DF - CF)*rs
         muxP  = fourThirds * exPara 
         muxF  = fourThirds * exF
         mucUp = mucP + fofzeta*(mucF-mucP) + (ecF-ecP)   *( 1.0-zeta)*dfofzeta
         mucDn = mucP + fofzeta*(mucF-mucP) + (ecF-ecP)   *(-1.0-zeta)*dfofzeta
         muxUp = muxP + fofzeta*(muxF-muxP) + (exF-exPara)*( 1.0-zeta)*dfofzeta
         muxDn = muxP + fofzeta*(muxF-muxP) + (exF-exPara)*(-1.0-zeta)*dfofzeta
         answer(1) = (muxUp + mucUp)
         answer(2) = (muxDn + mucDn)
      endif
      answer(1) = 0.5_double * answer(1)
      answer(2) = 0.5_double * answer(2)
      answer(3) = 0.5_double * answer(3)
   endif

   ! Consider the core charge where zeta = 0 so f(zeta) = 0.
   if (coreRho >= smallThresh) then
      alpha = (4.0_double/(9.0_double*pi))**0.33333333_double
      rs = (3.0_double/(4.0_double*pi*coreRho))**0.33333333_double
      sqrtrs = sqrt(rs)
      ex0 = 3.0_double/(2.0_double * pi * alpha)

      if (rs >= 1.0) then
         g   = -0.1423_double
         B1P =  1.0529_double
         B2P =  0.3334_double

         ! answer(4) = ex + ec (core charge) (zeta = 0)
         ec = g / (1 + B1P*sqrtrs + B2P*rs)
         answer(4) = (-ex0 / rs + ec)
      else
         AP =  0.03110_double
         BP = -0.04800_double
         CP =  0.00200_double
         DP = -0.01160_double
         lnrs = log(rs)

         ! answer(4) = ex + ec (core charge) (zeta = 0)
         ec = AP*lnrs + BP + CP*lnrs + DP*rs
         answer(4) = (-ex0 / rs + ec)
      endif

      answer(4) = 0.5_double * answer(4)
   endif
end subroutine ceperleyAlderSP


subroutine oldEXCORR(rh,sold,rhc,answer)

   ! Include kind definitions and the spin polarized vBH constants
   use O_Kinds

!   implicit none
   implicit real*8(a-h,o-z)

   ! Define dummy variables passed to this function and the return value.
   ! rh=+ totalRho, sold = spinDiffRho,rhc = coreRho
   real (kind=double), intent(in) :: rh,sold,rhc
   real (kind=double), intent(out), dimension (4) :: answer

! rh = total rho (up + down)
! sold = rho difference (up - down)
! rhc = total rho from the core orbitals (up + down)

! a= 1/(2^(0.3333333333))
! t43 = 2^(4/3)
! alpha = (4/9/pi)^0.3333333333
! topia = 2/pi/alpha
! topia34 = topia * 3/4
! conv = (3/4/pi)^0.3333333333

   real (kind=double) :: vxcup,vxcdn,excrs,excrc

   data a/0.793700525984d0/,third/0.3333333333333d0/,b/0.7734d0/, &
         & c21/21.d0/,alpha/0.521061762d0/,topia/1.221774115d0/, &
         & gm/5.129762808d0/,t43/2.519842099d0/, &
         & topia34/0.916330586d0/,fthd/1.333333333333333d0/
   data conv/0.6203504d0/
   data small/1.d-25/

   answer(:) = 0.0_double

   vxcup=0.0d0
   vxcdn=0.0d0
   excrs=0.0d0
   excrc=0.0d0

   if (rh < small) return
   x=(rh + sold)/(rh*2.d0)
   if (x.lt.0.d0) x=0.d0
   if(x.gt.1.d0) x=1.d0
   xs=(2.d0*x)**third-1.d0  ! Looks like x sub 1, but not quite.
   xsd=(2.d0*(1.d0-x))**third  ! Looks like 1/(x sub 2)
   rs=conv/(rh**third)
   z=c21/rs  ! Appears to be defined in hedin-Lundqvist x = rs/A,A=21. So z=1/x
   zp=t43*z
!.....
!.....Set up correlation functionals (will get facotor of topia/rs later.)
!.....
   ecp=g1(z)
   ecf=g1(zp)/a
!.....exchange functionals
   exp=-topia34/rs
   exf=exp/a
!.....
!.....functionals for the potentials
!.....
   ucp=g2(z)
   ucf=g2(zp)/a
   exfac=1.d0+ucp
   vc=gm*(ecf-ecp)
   exfac=exfac +xs*(1.d0+vc)
   tc=ucf-ucp-fthd*(ecf-ecp)
   fx=(x**fthd + (1.d0-x)**fthd-a)/(1.d0-a)
   exfac=exfac+tc*fx
!.....
!.....potentials
!.....vxcup=(uxp+vc)*(2x)**(1/3) + uxp-vc+tc*fx
!..... Von Barth-Hedin, Equation 6.2
!.....
   vxcup=-exfac*topia/rs
   vxcdn=(-topia/rs)*((1.d0+vc)*xsd + ucp + Tc*fx-vc)
!..... convert to hartree units
   vxcup=vxcup/2.d0
   vxcdn=vxcdn/2.d0
!.....
!.....energy FUNCTIONal
!.....
   excp=exp-topia*ecp/rs
   excf=exf-topia*ecf/rs
   excrs=excp+ (excf-excp)*fx
!..... to hartrees
   excrs=excrs/2.d0
   if(rhc.lt.small) then
      answer(1) = vxcup
      answer(2) = vxcdn
      answer(3) = excrs
      return
   endif
!..... core
   rsc=conv/(rhc**third)
   zc=c21/rsc
   zpc=t43*zc
   ecpc=g1(zc)
   ecfc=g1(zpc)/a
!.....exchange functionals
   expc=-topia34/rsc
   exfc=expc/a
   fxc=(2.0*0.5**fthd-a)/(1.d0-a)
   excpc=expc-topia*ecpc/rsc
   excfc=exfc-topia*ecfc/rsc
   excrc=excpc + (excfc-excpc)*fxc
!..... to hartrees
   excrc=excrc/2.d0

   answer(1) = vxcup
   answer(2) = vxcdn
   answer(3) = excrs
   answer(4) = excrc

end subroutine oldEXCORR

subroutine pbe96(rho,rhox,rhoy,rhoz,rhoxx,rhoxy,rhoxz,rhoyy,rhoyz,rhozz,answer)
   use O_Kinds
   use O_Constants

   ! Referenced papers

   ! PhysRevLett 77.3865 "Generalized Gradient Approximation Made Simple"

   ! PhysRevB 75.195108 "Functional form of the generalized gradient 
   ! approximation for exchange: The PBEalpha functional"


   implicit none

   ! Subroutine inputs
   real (kind=double) :: rho,rhox,rhoy,rhoz,rhoxx,rhoxy,rhoxz,rhoyy,rhoyz,rhozz

   ! Subroutine outputs
   real (kind=double), dimension (2) :: answer

   ! Variables used in both exchange and correlation calculations
   real (kind=double) :: seitzRadius

   ! Variables declared for exchange energy and potential
   real (kind=double) :: AX,Mu,Kappa,u,v,EXUNIF,FXPBE,s,kf,ks,gradN,Fs,Fss
   real (kind=double) :: exchangeEnergy,exchangePotential

   ! Variables used for correlation energy and potential
   real (kind=double) :: FZZ,F,EC,ECRS,FZ,ECZET,COMM,COMM2,COMM3,VCUP
   real (kind=double) :: VCDN,DVCUP,DVCDN,G,PON,B,Q4,Q5,H,T
   real (kind=double) :: Beta,Gamm,DELT,ETA,GAM
   real (kind=double) :: correlationEnergy,correlationPotential
   real (kind=double) :: Q0,Q2,Q3,BG,BEC,FAC,FACT2,Q9
   real (kind=double) :: FACT1,FACT0,FACT5,FACT3,HRST,GZ
   real (kind=double) :: HB,HBT,HRS,HT,HZT,HTT,HZ,Q8,PREF
   real (kind=double) :: VV, UU, WW
   ! GCOR2 inputs and outputs
   real (kind=double) :: EU,EURS,RTRS,EP,EPRS,ALFM,ALFRSM


   ! The local Seitz radius.
   seitzRadius = (3D0/4D0*pi*rho)**(1D0/3D0)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!                                                                !!
   !!                             EXCHANGE                           !!
   !!                                                                !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   ! The local Fermi wave vector.
   kf = (3D0*(pi**2D0)*rho)**(1D0/3D0)

   ! The Thomas-Fermi screening wave number.
   ks = sqrt(4D0*kf/pi)

   gradN = sqrt(rhox*rhox+rhoy*rhoy+rhoz*rhoz)
   
   ! Dimensionless density gradient.
   s = gradN/(2D0*kf*rho)

   AX = -0.7385587663820224_double 
   Mu = 0.2195149727645171_double
   Kappa = 0.8040_double

   !
   EXUNIF = AX*rho**(1D0/3D0)

   ! Compare to eqn 14 from PhysRevLett 77.3865
   FXPBE = 1D0 + Kappa - Kappa/(1D0+(Mu/Kappa)*(s**2D0))

   ! Exchange energy
   exchangeEnergy = EXUNIF*FXPBE
   
   !!!!!!!!!!!!!!!    EXCHANGE POTENTIAL     !!!!!!!!!!!!!!!!!!!!!!

   u = ((rhox*(rhox*rhoxx+rhoy*rhoxy+rhoz*rhoxz)/gradN)+ &
      & (rhoy*(rhox*rhoxy+rhoy*rhoyy+rhoz*rhoyz)/gradN)+ &
      & (rhoz*(rhox*rhoxz+rhoy*rhoyz+rhoz*rhozz)/gradN))/ &
      & (8D0*(rho**2D0)*(kf**3D0))

   v = (rhoxx +rhoyy + rhozz)/(4D0*rho*(kf**2D0))

   ! Find first and second derivatives of FX w.r.t s.

   ! Compare to PhysRevB 75.195108 equation 7
   Fs=2D0*Kappa*(Mu/Kappa)/((1D0+(Mu/Kappa)*(s**2D0))**2D0)
   Fss=-4D0*(Mu/Kappa)*s*Fs/(1D0+(Mu/Kappa)*(s**2D0))
   exchangePotential = EXUNIF*((4D0/3D0)*FXPBE-(u-(4D0/3D0)*s**3D0)*FSS-v*FS)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!                                                                !!
   !!                       CORRELATION                              !!
   !!                                                                !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   RTRS=sqrt(seitzRadius)

   CALL GCOR2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0, &
      & 0.49294D0,RTRS,EU,EURS)
   CALL GCOR2(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0, &
      & 0.62517D0,RTRS,EP,EPRS)
   CALL GCOR2(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0, &
      & 0.49671D0,RTRS,ALFM,ALFRSM)   

   !!!!!!!!!!!!!!!!    LSD POTENTIAL and ENERGY   !!!!!!!!!!!!!!!!!!!!!!

   GAM = 0.5198420997897463_double
   FZZ = 8D0/(9D0*GAM)
 

   F=((1D0**(4D0/3D0))+(1D0**(4D0/3D0))-2D0)/GAM
   EC=EU-ALFM*F/FZZ
   ! ECRS=dEC/dRS
   ECRS=EURS-ALFRSM*F/FZZ

   ! ECZET=dEC/dZETA
   ECZET= -(ALFM/FZZ)
 
   ! FZ=dF/dZETA
   FZ=(4D0/3D0)*((1D0**(1D0/3D0))-(1D0**(1D0/3D0)))/GAM
   COMM=EC-seitzRadius*ECRS/3
   VCUP=COMM+ECZET
   VCDN=COMM-ECZET
   DVCUP=0.0_double
   DVCDN=0.0_double

   !!!!!!!!!!!!!!!!!!!!   PBE CORRELATION ENERGY      !!!!!!!!!!!!!!!!!!!!!
     
   Beta = 0.0667245506031492_double
   Gamm = 0.0310906908696549_double
   DELT = Beta/Gamm
   T = gradN/(2D0*rho*ks)

   ! G = PHI(ZETA) from 77.3865
   G=((1D0**(2D0/3D0))+(1D0**(2D0/3D0)))/2D0

   ! Compare to terms in equation 8 from 77.3865
   PON=-EC/((G**3D0)*Gamm)

   ! Compare to equation 8 from 77.3865. A = B
   B=DELT/(DEXP(PON)-1D0)
   Q4=1D0+B*(T**2D0)
   Q5=1D0+B*(T**2D0)+(B**2D0)*(T**4D0)

   ! Compare to equation 7 from 77.3865
   H=(G**3D0)*(Beta/DELT)*DLOG(1D0+DELT*Q4*(T**2D0)/Q5)

   ! Compare to equation 3 from 77.3865
   correlationEnergy=EC+H

   !!!!!!!!!!!!!!!!!    PBE CORRELATION POTENTIAL      !!!!!!!!!!!!!!

   UU= ((rhox*(rhox*rhoxx+rhoy*rhoxy+rhoz*rhoxz)/gradN)+ &
      & (rhoy*(rhox*rhoxy+rhoy*rhoyy+rhoz*rhoyz)/gradN)+ &
      & (rhoz*(rhox*rhoxz+rhoy*rhoyz+rhoz*rhozz)/gradN))/ &
      & (8D0*(rho**2D0)*(ks**3D0))

   VV= (rhoxx+rhoyy+rhozz)/(4D0*rho*(ks**2D0))

   WW= 0D0

   GZ=((1D0+ETA)**(-1D0/6D0)- &
      & (1D0+ETA)**(-1/6D0))/3D0
   FAC=DELT/B+1D0
   BG=-3D0*(B**2D0)*EC*FAC/(Beta*(G**4D0))
   BEC=(B**2D0)*FAC/(Beta*(G**3D0))
   Q8=(Q5**2D0)+DELT*Q4*Q5*(T**2D0)
   Q9=1D0+2D0*B*(T**2D0)
   HB=-Beta*(G**3D0)*B*(T**6D0)*(2D0+B*(T**2D0))/Q8
   HRS=-(seitzRadius/3D0)*HB*BEC*ECRS
   FACT0=2D0*DELT-6D0*B
   FACT1=Q5*Q9+Q4*(Q9**2D0)
   HBT=2D0*Beta*(G**3D0)*(T**4D0)*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
   HRST=(seitzRadius/3D0)*(T**2D0)*HBT*BEC*ECRS
   HZ=3D0*GZ*H/G+HB*(BG*GZ+BEC*ECZET)
   HT=2D0*Beta*(G**3D0)*Q9/Q8
   HZT=3D0*GZ*HT/G+HBT*(BG*GZ+BEC*ECZET)
   FACT2=Q4*Q5+B*(T**2D0)*(Q4*Q9+Q5)
   FACT3=2D0*B*Q5*Q9+DELT*FACT2
   HTT=4D0*Beta*(G**3D0)*T*(2D0*B/Q8-(Q9*FACT3/Q8)/Q8)
   COMM2=H+HRS+HRST+(T**2D0)*HT/6D0+7D0*(T**3D0)*HTT/6D0
   PREF=HZ-GZ*(T**2D0)*HT/G
   FACT5=GZ*(2D0*HT+T*HTT)/G
   COMM3=COMM2-UU*HTT-VV*HT-WW*(HZT-FACT5)
   DVCUP=COMM3+PREF
   DVCDN=COMM3-PREF

   correlationPotential = VCUP+DVCUP

   !!!!!!!!!!!!!!    Combining Correlation and Potential    !!!!!!!!!!!!!!

   ! XC_potential.
   answer(1) = exchangePotential + correlationPotential
   ! Total XC_energy.
   answer(2) = exchangeEnergy + correlationEnergy

end subroutine pbe96

subroutine GCOR2(A,A1,B1,B2,B3,B4,RTRS,GG,GGRS)
   
   implicit none

   real (kind=double) :: A,A1,B1,B2,B3,B4,RTRS,GG,GGRS,Q0,Q1,Q2,Q3
 
   Q0=-2*A*(1+A1*(RTRS**2))
   Q1=2*A*RTRS*(B1+RTRS*(B2+RTRS*(B3+B4*RTRS)))
   Q2=DLOG(1+1/Q1)
   GG=Q0*Q2
   Q3=A*(B1/RTRS+2*B2+RTRS*(3*B3+4*B4*RTRS))
   GGRS=-2*A*A1*Q2-Q0*Q3/(Q1*(1+Q1))
   
end subroutine GCOR2

!
!*************************************************************
!
function g1(x)
   implicit real*8(a-h,o-z)
   data one/1.d0/,thd/0.33333333333333d0/,hf/0.5d0/, &
         & c0/18.d0/,c1/0.075d0/,c2/-0.107142857d0/, &
         & c3/0.166666666666666d0/, &
         & c4/-0.3d0/,c5/0.75d0/,b/0.7734d0/
   if(x.ge.0.1d0) then
      x1=one/x
      g1=(one+x1*x1*x1)*log(one+one/x1)+hf*x1-x1*x1-thd
      g1=g1*b*x1
   else
      g1=((((-x/c0+c1)*x+c2)*x+c3)*x+c4)*x+c5
      g1=g1*b
   endif
end function g1
!
!*************************************************************
!
function g2(x)
   implicit real*8(a-h,o-z)
   data one/1.d0/,b/0.7734d0/,c0/6.d0/,c1/0.2d0/,c2/-0.25d0/, &
         & c3/0.333333333333333d0/,c4/-0.5d0/
   if(x.ge.0.1d0) then
      x1=one/x
      g2=b*x1*log(one+one/x1)
   else
      g2=((((-x/c0+c1)*x+c2)*x+c3)*x+c4)*x+one
      g2=g2*b
   endif
end function g2


subroutine cleanUpPotentialUpdate

   implicit none

   if (allocated (usedPotCoeffs)) then
      deallocate (usedPotCoeffs)
      deallocate (guessedPotCoeffs)
      deallocate (tempUsedPotCoeffs)
      deallocate (tempGuessedPotCoeffs)
      deallocate (tempPotCoeffs)
   endif

!   if (allocated (xl0)) then
!      deallocate (xl0)
!      deallocate (xl1)
!      deallocate (xl2)
!      deallocate (yl0)
!      deallocate (yl1)
!      deallocate (yl2)
!   endif

end subroutine cleanUpPotentialUpdate


end module O_PotentialUpdate
