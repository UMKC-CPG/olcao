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
   real (kind=double), allocatable, dimension (:,:) :: potDifference ! Holds
         ! the difference between the currently used potential coefficients and
         ! the current guess for the next set of potential coefficients.
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

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed dummy variables.
   real (kind=double) :: totalEnergy

   ! Define the local variables used in this subroutine.
   integer :: i,j,k ! Loop index variables
   integer :: hdferr
   integer :: info
   integer :: numOpValues
   integer :: currentType
   integer :: currentNumAlphas
   integer :: currentCumulAlphaSum
   integer :: potTypeInitIndex
   integer :: potTypeFinIndex
   integer :: siteIndex
   integer :: potTermCount
   real (kind=double) :: fittedCharge
   real (kind=double) :: fittedSpinDiffCharge
   real (kind=double) :: electronDiff
   real (kind=double) :: currentNucCharge
   real (kind=double) :: sumIntegratedCharge
   real (kind=double) :: kineticEnergy
   real (kind=double) :: exchCorrEnergy
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

   real (kind=double) :: testableDelta
   real (kind=double) :: radialWeightSum
   real (kind=double) :: weightedPotDiff
   real (kind=double) :: totalMagneticMoment
   real (kind=double), dimension (4) :: currentExchCorrPot
   real (kind=double), allocatable, dimension (:,:) :: generalRho
   real (kind=double), allocatable, dimension (:,:) :: tempOverlap
   real (kind=double), allocatable, dimension (:,:) :: elecStatPot
   real (kind=double), allocatable, dimension (:,:) :: exchCorrPot
   real (kind=double), allocatable, dimension (:,:) :: exchCorrRho
   real (kind=double), allocatable, dimension (:,:) :: exchCorrRhoCore
   real (kind=double), allocatable, dimension (:,:) :: exchCorrRhoSpin
   real (kind=double), allocatable, dimension (:,:) :: realSpacePotDiff
   real (kind=double), allocatable, dimension (:)   :: averageDelta
   real (kind=double), allocatable, dimension (:)   :: maxDelta
   real (kind=double), allocatable, dimension (:)   :: typesMagneticMoment

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
      write (20, *) 'dposvx failed. INFO= ', info
      stop
   endif

   ! Compute the fitted valence charge.  (Up + Down for spin polarized calcs.)
   fittedCharge = dot_product(intgConsts(:potDim),generalRho(:potDim,1))
   if (spin == 2) then
      ! Compute the fitted valence charge difference.  (Up - Down)
      fittedSpinDiffCharge = dot_product(intgConsts(:potDim),&
            & generalRho(:potDim,3))

      ! Preserve the spin (UP-DOWN) charge density coefficients in 
      !   generalRho(:,8).
      generalRho(:potDim,8) = generalRho(:potDim,3)
   endif

   ! Calculate the difference between the electron number and the fitted charge
   !   determined above.  This represents the difficulty of the Gaussians to
   !   accurately represent the true charge density.
   electronDiff = numElectrons - fittedCharge

   ! Record the electron fitting error.
   write (20,fmt="(a37,f12.8,a8,i12)") &
         & 'Unconstrained Vale Electron Fit Error',&
         & electronDiff, ' Out Of ',numElectrons
   call flush (20)

   ! Correct and store the valence charge coefficients in generalRho(:,2).
   !   This would overwrite the spin difference for spin polarized
   !   calculations, but that was already saved in generalRho(:,8).
!generalRho(:potDim,2) = generalRho(:potDim,1)*numElectrons/fittedCharge
   generalRho(:potDim,2) = generalRho(:potDim,1) + generalRho(:potDim,2) * &
         & electronDiff / sum(intgConsts(:potDim)*generalRho(:potDim,2))

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
            & generalRho(:potDim,3),pot,hdferr)
      if (hdferr /= 0) stop 'Failed to read core charge density'
   else
      generalRho(:potDim,3) = 0.0_double
   endif

   ! At present generalRho(:,2) holds the valence coefficients for both the
   !   spin polarized and non spin polarized case.  The generalRho(:,3) holds
   !   the core charge coefficients.  We will now compute the total and store
   !   it in generalRho(:,1).
   generalRho(:potDim,1) = generalRho(:potDim,2) + generalRho(:potDim,3)

   

   ! Divide the charge so far into neutral spherically symmetric pieces and a
   !   residual part.  The previous sentence is the original documentation, but
   !   what appears to actually be done is just sum together the core charge
   !   density (3) with the adjusted valence charge density (2) and store the
   !   result in (1).  In the case of spin polarized calculations, (2) is also
   !   the total valence charge density, and is considered to be the up+down.
   !   generalRho(:,1) = TOTAL      CHARGE  (CORE + VALENCE)
   !   generalRho(:,2) = VALENCE    CHARGE
   !   generalRho(:,3) = CORE       CHARGE
   !   generalRho(:,4) = NEUTRAL    CHARGE
   !   generalRho(:,5) = RESIDUAL   CHARGE 1?
   !   generalRho(:,6) = RESIDUAL   CHARGE 2?
   !   generalRho(:,7) = NUCLEAR    CHARGE
   !   generalRho(:,8) = DIFFERENCE CHARGE (SPIN DIFFERENCE or ABSENT)
   generalRho(:,4) = generalRho(:,1)
   generalRho(:,5) = 0.0_double
   generalRho(:,6) = 0.0_double
   generalRho(:,7) = 0.0_double

   ! This part is not well documented.  But later it is mentioned that the
   !   goal is to simulate the nuclear charges by affecting the shortest
   !   range potential in a way by the nuclear charge.
   do i = 1,numPotTypes

      ! Initialize variables for the current potential type
      currentNumAlphas = potTypes(i)%numAlphas
      currentNucCharge = potTypes(i)%nucCharge * potTypes(i)%multiplicity
      currentCumulAlphaSum = potTypes(i)%cumulAlphaSum

      ! Identify the indices of the summations.
      potTypeInitIndex = currentCumulAlphaSum + 1
      potTypeFinIndex  = currentCumulAlphaSum + currentNumAlphas

      ! Calculate the sum of the coefficients of integration times the total
      !   charge determined for this potential type.
      sumIntegratedCharge = sum(intgConsts(potTypeInitIndex:potTypeFinIndex) * &
            & generalRho(potTypeInitIndex:potTypeFinIndex,1))

      ! Assign values to generalRho (:,4:7)
      generalRho(potTypeInitIndex,5)=(sumIntegratedCharge-currentNucCharge) / &
            & intgConsts(potTypeInitIndex)

      generalRho(potTypeInitIndex,4) = generalRho(potTypeInitIndex,1) - &
            & generalRho(potTypeInitIndex,5)

      generalRho(i,6) = generalRho(potTypeInitIndex,5)

      generalRho(potTypeFinIndex,7) = -currentNucCharge / &
            & intgConsts(potTypeFinIndex)
   enddo


   ! Allocate space for the electrostatic potential
   allocate (elecStatPot (potDim,6))


   ! Allocate space for the non-local real-space electrostatic operator matrix.
   allocate (nonLocalNeutQPot (potDim,potDim))


   ! Read the non-local real-space electrostatic operator.
   call h5dread_f (nonLocalNeutQPot_did,H5T_NATIVE_DOUBLE,&
         & nonLocalNeutQPot(:,:),potPot,hdferr)
   if (hdferr /= 0) stop 'Failed to read non local neut q pot'


   ! Operate on the spherically symmetric neutral charge.
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


   ! Operate on the (valence - residual) and core charges seperately.
   do i = 1, potDim
      elecStatPot(i,2) = sum(localNeutQPot(:,i) * &
            & (generalRho(:,2) - generalRho(:,5)))
      elecStatPot(i,3) = sum(localNeutQPot(:,i) * generalRho(:,3))
   enddo

   deallocate (localNeutQPot)


   ! Read in the local nuclear integral vector
   call h5dread_f (localNucQPot_did,H5T_NATIVE_DOUBLE,&
         & elecStatPot(:,5),pot,hdferr)
   if (hdferr /= 0) stop 'Failed to read local nuc q pot'



   ! Read in the non local nuclear integral vector (gaussian screeneed)
   call h5dread_f (nonLocalNucQPot_did,H5T_NATIVE_DOUBLE,&
         & elecStatPot(:,4),pot,hdferr)
   if (hdferr /= 0) stop 'Failed to read non local nuc q pot'



   ! Allocate space for the reciprocal-space potential operator.
   allocate (nonLocalResidualQ (numPotTypes,potDim))

   ! Read in the reciprocal-space potential operator
   call h5dread_f (nonLocalResidualQ_did,H5T_NATIVE_DOUBLE,&
         & nonLocalResidualQ(:,:),potTypesPot,hdferr)
   if (hdferr /= 0) stop 'Failed to read local residual q pot'



   ! Operate on the residual charge.
   ! elecStat(:,6) = Electrostatic potential from the residual charge.
   ! sum(elecStat(:,1,2,3,5,6)) = Total electrostatic potential.
   ! Then subtract the gaussian screened nuclear charge (elecStat(:,4)) from
   !   the above sum to obtain the electrostatic potential needed for the
   !   SCF calculations.
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
   elecStatEnergy = sum(0.5_double  * (generalRho(:,1) + generalRho(:,7)) * &
         & (elecStatPot(:,1) + elecStatPot(:,6)) + generalRho(:,2) * &
         & 0.5_double * (elecStatPot(:,2) + elecStatPot(:,3) + &
         & elecStatPot(:,5)) + 0.5_double * (generalRho(:,2) - &
         & generalRho(:,5))  * (elecStatPot(:,3) + elecStatPot(:,5)))


   ! Correct for the difference between the actual and fitted valence charge
   !   integrated against the gaussian screened nucleus.

   elecStatEnergyDiff = elecStatEnergy
   elecStatEnergy = elecStatEnergy - sum(generalRho(:,2) * elecStatPot(:,4)) + &
         & nucPotTrace(1) ! UP+DOWN or Total
   elecStatEnergyDiff = elecStatEnergyDiff - elecStatEnergy


   write (20,*) 'Actual and fitted charge-screened nucleus difference',&
         & elecStatEnergyDiff
   call flush (20)


   ! NOW, do a least squares fit for the potential with the overlap matrix.
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
      write (20, *) 'dposvx failed. INFO= ', info
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

         if (spin == 2) then
            spinDiffSum = sum(exchRhoOp(:potDim,j,1) * generalRho(:potDim,8))
            exchCorrRhoSpin(1,j) = spinDiffSum  ! Valence spin difference
         else
            exchCorrRhoSpin(1,j) = 0.0_double  ! Valence spin difference
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

         if (GGA == 0) then
            if (coreSum < 0.0_double) then
               coreSum = 0.0_double
            endif
         else
            if (coreSum < smallThresh) then
               coreSum = smallThresh
            endif
  
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

         if (totalSum < smallThresh) then
! This is totally ridiculous.  If the following statement is not present, then
!   the program will not work.  The values for totalSum will not be
!   calculated correctly.  For the first potential all the values will be less
!   than smallThresh.  For the later potential sites (i) the values will be
!   correct.  I can only assume that the error lies in a bad allocation
!   statement somewhere, or an access outside of some array boundry, but I have
!   been unable to determine where.  If these write statements are present,
!   then the program will function correctly even though the statements are
!   never executed.  This tells me that the compiler is optimizing differently
!   for the two cases on the HP-UX ia64 machine sirius.  Further, on the Tru64
!   machine hilbert this works fine without these write statements.  I tried
!   to compile both cases with the -C flag to check array bounds.  The Tru64
!   compiler gave no errors, while the HP-UX one wouldn't even get past the
!   parseInput subroutine that has been used flawlessly for years.  Annoying.
!write (20,*) "totalSum=",totalSum
!call flush (20)
            totalSum = smallThresh
!write (20,*) "totalSum=",totalSum
!call flush (20)
         endif

         exchCorrRho (1,j) = totalSum
         exchCorrRhoCore (1,j) = coreSum

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
      enddo

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
                  exchCorrPot(:potDim,k) = exchCorrPot(:potDim,k) + &
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
                  exchCorrPot(:potDim,k) = exchCorrPot(:potDim,k) + &
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
                  exchCorrPot(:potDim,k) = exchCorrPot(:potDim,k) + &
                        & radialWeight(j) * currentExchCorrPot(k) * &
                        & exchRhoOp(:potDim,j,1)
               enddo
            enddo
      elseif (xcCode == 150) then
            ! Ceperley and Alder Exchange-Correlation (LSDA)
            do j = 1, numRayPoints
               call ceperleyAlderSP(exchCorrRho(1,j),exchCorrRhoSpin(1,j),&
                     & exchCorrRhoCore(1,j),currentExchCorrPot(1:4))
               do k = 1,4
                  exchCorrPot(:potDim,k) = exchCorrPot(:potDim,k) + &
                        & radialWeight(j) * currentExchCorrPot(k) * &
                        & exchRhoOp(:potDim,j,1)
               enddo
            enddo

      elseif (xcCode == 151) then
            ! von Barth and Hedin Exchange-Correlation (LSDA)
            do j = 1, numRayPoints
               call vonBarthHedin(exchCorrRho(1,j),exchCorrRhoSpin(1,j),&
                     & exchCorrRhoCore(1,j),currentExchCorrPot(1:4))
               do k = 1,4
                  exchCorrPot(:potDim,k) = exchCorrPot(:potDim,k) + &
                        & radialWeight(j) * currentExchCorrPot(k) * &
                        & exchRhoOp(:potDim,j,1)
               enddo
            enddo
      elseif (xcCode == 152) then
            ! Old SPmain.for Exchange Correlation (LSDA) (Should be like the
            !   von barth and Hedin above.)
            do j = 1, numRayPoints
               call oldEXCORR(exchCorrRho(1,j),exchCorrRhoSpin(1,j),&
                     & exchCorrRhoCore(1,j),currentExchCorrPot(1:4))
               do k = 1,4
                  exchCorrPot(:potDim,k) = exchCorrPot(:potDim,k) + &
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
               exchCorrPot(:potDim,1) = exchCorrPot(:potDim,1) + &
                  & radialWeight(j) * currentExchCorrPot(1) * &
                  & exchRhoOp(:potDim,j,1)
               exchCorrPot(:potDim,2) = exchCorrPot(:potDim,2) + &
                  & radialWeight(j) * currentExchCorrPot(2) * &
                  & exchRhoOp(:potDim,j,1)
               exchCorrPot(:potDim,3) = exchCorrPot(:potDim,3) + &
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
!              exchCorrPot(:potDim,1) = exchCorrPot(:potDim,1) + &
!                 & radialWeight(j) * currentExchCorrPot(1) * &
!                 & exchRhoOp(:potDim,j,1)
!              exchCorrPot(:potDim,2) = exchCorrPot(:potDim,2) + &
!                 & radialWeight(j) * currentExchCorrPot(2) * &
!                 & exchRhoOp(:potDim,j,1)
!              exchCorrPot(:potDim,3) = exchCorrPot(:potDim,3) + &
!                 & radialWeight(j) * currentExchCorrPot(4) * &
!                 & exchRhoOp(:potDim,j,1)
!              exchCorrPot(:potDim,4) = exchCorrPot(:potDim,3) + &
!                 & radialWeight(j) * currentExchCorrPot(4) * &
!                 & exchRhoOp(:potDim,j,1)
!           enddo
      endif
   enddo

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
      write (20, *) 'dposvx failed. INFO= ', info
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



   kineticEnergy = kineticEnergyTrace(1) ! UP + DOWN or TOTAL

   ! Calculate the total energy
   totalEnergy = kineticEnergy + elecStatEnergy + exchCorrEnergy


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
   if (.not. allocated(potDifference)) then
      allocate (potDifference(potDim,spin)) ! For measuring convergence.
      allocate (usedPotCoeffs(potDim,feedbackLevel+1,spin)) ! Anderson's x
      allocate (guessedPotCoeffs(potDim,feedbackLevel+1,spin)) ! Anderson's y

      do i = 1, spin
         do j = 1, feedbackLevel+1
            usedPotCoeffs(:,j,i) = potCoeffs(:,i)
            guessedPotCoeffs(:,j,i) = elecStatPot(:,1) + elecStatPot(:,2) + &
                  & exchCorrPot(:,i)
         enddo
      enddo
   endif


   do i = 1, spin
      usedPotCoeffs(:,1,i) = potCoeffs(:,i) ! Anderson's x^l
      guessedPotCoeffs(:,1,i) = elecStatPot(:,1) + elecStatPot(:,2) &
            & + exchCorrPot(:,i) ! Anderson's y^l

      ! Form the difference between the potential coefficients that have been
      !   guessed for the next iteration and those that were used in the current
      !   iteration.
      potDifference(:,i) = guessedPotCoeffs(:,1,i) - usedPotCoeffs(:,1,i)

   enddo


   ! Deallocate matrices that have been copied to temporary arrays.
   deallocate (elecStatPot)
   deallocate (exchCorrPot)

   if ((currIteration == 1) .or. (currIteration == lastIteration)) then
      write (20,*) 'OUTPUT INTERESTING STUFF HERE'
   endif


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

      ! Read the mesh information.
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
            realSpacePotDiff(k,j) = sum(exchRhoOp(:,k,1) * potDifference(:,j))
         enddo

         ! Determine the delta to be tested.
         maxDelta(siteIndex)     = maxval(realSpacePotDiff(:numRayPoints,j))
         radialWeightSum         = sum (radialWeight(:numRayPoints))
         weightedPotDiff         = sum (radialWeight(:numRayPoints) * &
               & realSpacePotDiff(:numRayPoints,j)**2)
         averageDelta(siteIndex) = sqrt(weightedPotDiff / radialWeightSum)
         testableDelta           = max (averageDelta(siteIndex),testableDelta)
      enddo
   enddo

   ! Deallocate arrays to obtain the testableDelta
   deallocate (realSpacePotDiff)
   deallocate (averageDelta)
   deallocate (maxDelta)
   deallocate (exchRhoOp)
   deallocate (radialWeight)

   ! Blend the currently used potential, the current guess, and (possibly) a set
   !   number of previous potential functions in an effort to accelerate the
   !   convergence while retaining smooth character to the convergence process.
   call blendPotentials()

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
         write (8,fmt="(a20)") "TOTAL OR SPIN_UP"
      else
         write (8,fmt="(a20)") "SPIN_DN"
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
   write (8,fmt="(a20)") "NUM PLUSUJ TERMS"
   write (8,fmt="(i5)") numPlusUJAtoms
   do i = 1, 2 ! Spin up and down.

      ! Write the header tag for this block of data.
      if (i == 1) then
         write (8,fmt="(a20)") "TOTAL OR SPIN_UP"
      else
         write (8,fmt="(a20)") "SPIN_DN"
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

         typesMagneticMoment(i) = sum(&
               & generalRho(potTypeInitIndex:potTypeFinIndex,8))

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

   write (14,*) currIteration, kineticEnergy, elecStatEnergy, &
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
   deallocate (nucPotTrace)
   deallocate (kineticEnergyTrace)
   deallocate (tempOverlap)
   deallocate (potAlphaOverlap)

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
!   secant method would be easy to compute and it would get close to the right
!   direction. All we have to do then is make an algorithm that will get us as
!   close to the solution subspace as possible and then we should converge to
!   it. (Some issues remain, and they will be discussed later.)

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
!   are defined by the overlap. However, we could in the future consider
!   weighting factors that are additinally designed to favor convergence for
!   certain terms over convergence in other terms in the potential function
!   coefficient array. The rational for pursuing a weighting factor is given in
!   the paragraph beginning with "Variants of these algorithms..." on page 554.
!   Serious consideration should be given to this line of thinking in the
!   future.

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
!   "relaxFactor". The purpose of Equation 4.9 is to blend or mix the currently
!   used potential coefficients x with the terms of the next guess y. Actually,
!   this equation will mix u and v which themselves are constructed from the x
!   and y values of the current and previous iterations. The u and v are each
!   built using a linearized secant-method approach. That is (looking just at u
!   with v being built in the same way) we consider three points 1, 2, and 3.
!   Point 3 is the new point (i.e. u) while point 2 is the current point (x^l)
!   and point 1 (x^l-1) is the previous point. We have that 3 = 2 - f(2)/f'(2)
!   which is basically a Newton's method formula. Then we approximate the
!   derivative f'(2) with f'(2) ~= [f(2)-f(1)] / [2-1]. Combining we get the
!   expression 3 = 2 - f(2)* [2-1] / [f(2)-f(1)] which when put into the form of
!   the Anderson paper is Equation 4.5: u^l = x^l + theta^l * (x^l-1 -x^l).
!   Theta will be determined next, but it represents the f(2)/[f(2)-f(1)] part
!   of the method. In actuality, it will not explicitly be equal to that part of
!   the method but will instead be selected to try to minimize another quantity.
!   Note that the signs are correct.

! The u vector is the base for the next iteration of coefficients and it is
!   blended with some v. The u vector is built from the current x and a
!   difference between the current and previous values of x. The difference is
!   multiplied by a coefficient (theta) that is specifically selected so as to
!   minimize the linearized residual R. We define R as R=0.5*(v-u,v-u). The
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
subroutine blendPotentials()

   use O_Kinds
   use O_Potential, only: spin, potDim, currIteration, feedbackLevel, &
         & relaxFactor, potCoeffs
   use O_ElectroStatics, only: potAlphaOverlap
   use O_LAPACKDPOSVX

   implicit none

   ! Define the dummy variables that are passed to this function.

   ! Define the local variables.
   integer :: info
   integer :: i,j,k,l
   integer :: maxFeedback
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
         ! Thus, the first index is potDim and the second is feedbackLevel+1.
   real (kind=double), allocatable, dimension(:,:) :: drl ! Differences
         ! between the first rl and the other rl arrays. The first index is
         ! potDim and the second is feedbackLevel. (The +1 is not needed here
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

   ! Allocate space for the operating data structures.
   allocate (tempArray(potDim))
   allocate (theta(maxFeedback+1)) ! The first index is special (1-sum(others)).
   allocate (rl(potDim,maxFeedback+1)) ! +1 to also hold the current iteration.
   if (maxFeedback > 0) then
      allocate (drl(potDim,maxFeedback)) ! No +1 needed.
   endif

   ! Presently we are treating the two spin directions independently. However,
   !   it might be a good idea in the future to try to integrate them together.
   do i = 1, spin

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
         rl(:,j) = guessedPotCoeffs(:,j,i) - usedPotCoeffs(:,j,i)
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
               do l = 1, potDim
                  tempArray(:) = potAlphaOverlap(:,l)*drl(l,j)
                  matrix(k,j) = matrix(k,j) + sum(tempArray(:)*drl(:,k))
                  if (k == 1) then
                     solutions(j) = solutions(j) + sum(tempArray(:)*rl(:,1))
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
         do j = maxFeedback, 1, -1

            ! Copy the actual matrix A and solutions B into temporary data
            !   structures because the originals will be destroyed within the
            !   dposvx subroutine.
            matrixTemp(:,:) = matrix(:,:)
            solutionsTemp(:,1) = solutions(:)

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
                  stop
               else
                  write (20,*) "Failed to solve DPOSVX for potential blending."
                  write (20,*) "Trying again with fewer feedback terms."
                  write (20,*) "Current feedback term was",j
               endif
               cycle
            endif

            ! Copy the results into the theta array taking note that the first
            !   index of theta is special. It will be used for
            !   1.0-sum(other_theta_values).
            theta(1) = 1.0_double - sum(solutionsTemp(1:maxFeedback,1))
            theta(2:maxFeedback+1) = solutionsTemp(1:maxFeedback,1)

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
      do j = 1, maxFeedback+1
         tempArray(:) = tempArray(:) + theta(j)*usedPotCoeffs(:,j,i)
      enddo
      potCoeffs(:,i) = (1.0_double - relaxFactor) * tempArray(:)

      tempArray(:) = 0.0_double
      do j = 1, maxFeedback+1
         tempArray(:) = tempArray(:) + theta(j)*guessedPotCoeffs(:,j,i)
      enddo
      potCoeffs(:,i) = potCoeffs(:,i) + relaxFactor * tempArray(:)

      ! Progressively shift the oldest iterations out by replacing them with an
      !   earlier (lower middle index) iteration.
      do j = feedbackLevel, 1, -1
         usedPotCoeffs(:,j+1,i) = usedPotCoeffs(:,j,i)
         guessedPotCoeffs(:,j+1,i) = guessedPotCoeffs(:,j,i)
      enddo
   enddo

   ! Deallocate space for the operating data structures of the algorithm.
   if (maxFeedback > 0) then
      deallocate (drl)
   endif
   deallocate (rl)
   deallocate (theta)
   deallocate (tempArray)

end subroutine blendPotentials


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
   real (kind=double) :: totalRho,spinDiffRho,coreRho
   real (kind=double), dimension (4) :: answer

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
   fourThirds = 4.0/3.0
   twoThirds  = 2.0/3.0
   oneThird   = 1.0/3.0


   if (totalRho > smallThresh) then

      alpha = (4.0_double/(9.0_double*pi))**0.33333333
      rs = (3.0_double/(4.0_double*pi*totalRho))**0.33333333
      zeta = spinDiffRho/totalRho
      ex0 = 3.0_double/(2.0_double * pi * alpha)

      if (rs >= 1.0) then
         g   = -0.1423
         B1P =  1.0529
         B1F =  1.3981
         B2P =  0.3334
         B2F =  0.2611
         sqrtrs = sqrt(rs)
         fofzeta = ((1.0+zeta)**fourThirds + (1.0-zeta)**fourThirds - 2.0) / &
               & (2.0**fourThirds - 2.0)
         dfofzeta = (fourThirds * (1.0+zeta)**oneThird - &
                  &  fourThirds * (1.0-zeta)**oneThird)/ &
                  &  (2.0**fourThirds - 2.0)
         ! answer(3) = ex + ec.
         ecP = g / (1 + B1P*sqrtrs + B2P*rs)
         ecF = g / (1 + B1F*sqrtrs + B2F*rs)
         ec  = ecP + fofzeta*(ecF-ecP)
         exPara = -ex0/rs
         exF = exPara * 2.0_double**oneThird
         ex  = exPara + fofzeta*(exF-exPara)
         answer(3) = ex + ec

         ! answer(1,2) = ux + uc
         mucP = ecP * ecP / g * (1 + 7.0/6.0*B1P*sqrtrs + 4.0/3.0*B2P*rs)
         mucF = ecF * ecF / g * (1 + 7.0/6.0*B1F*sqrtrs + 4.0/3.0*B2F*rs)
         muxP  = fourThirds * exPara 
         muxF  = fourThirds * exF
         mucUp = mucP + fofzeta*(mucF-mucP) + (ecF-ecP)   *( 1.0-zeta)*dfofzeta
         mucDn = mucP + fofzeta*(mucF-mucP) + (ecF-ecP)   *(-1.0-zeta)*dfofzeta
         muxUp = muxP + fofzeta*(muxF-muxP) + (exF-exPara)*( 1.0-zeta)*dfofzeta
         muxDn = muxP + fofzeta*(muxF-muxP) + (exF-exPara)*(-1.0-zeta)*dfofzeta
         answer(1) = (muxUp + mucUp)
         answer(2) = (muxDn + mucDn)
      else
         AP =  0.03110
         AF =  0.01555
         BP = -0.04800
         BF = -0.02690
         CP =  0.00200
         CF =  0.00070
         DP = -0.01160
         DF = -0.00480

         lnrs = log(rs) ! log is the natural log (ln) in fortran.
         fofzeta = ((1.0+zeta)**fourThirds + (1.0-zeta)**fourThirds - 2.0) / &
               & (2.0**fourThirds - 2.0)
         dfofzeta = (fourThirds * (1.0+zeta)**oneThird - &
                  &  fourThirds * (1.0-zeta)**oneThird)/ &
                  &  (2.0**fourThirds - 2.0)

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
               & oneThird*(2.0*DP - CP)*rs
         mucF  = AF*lnrs + (BF - oneThird*AF) + twoThirds*CF*rs*lnrs + &
               & oneThird*(2.0*DF - CF)*rs
         muxP  = fourThirds * exPara 
         muxF  = fourThirds * exF
         mucUp = mucP + fofzeta*(mucF-mucP) + (ecF-ecP)   *( 1.0-zeta)*dfofzeta
         mucDn = mucP + fofzeta*(mucF-mucP) + (ecF-ecP)   *(-1.0-zeta)*dfofzeta
         muxUp = muxP + fofzeta*(muxF-muxP) + (exF-exPara)*( 1.0-zeta)*dfofzeta
         muxDn = muxP + fofzeta*(muxF-muxP) + (exF-exPara)*(-1.0-zeta)*dfofzeta
         answer(1) = (muxUp + mucUp)
         answer(2) = (muxDn + mucDn)
      endif
      answer(1) = 0.5 * answer(1)
      answer(2) = 0.5 * answer(2)
      answer(3) = 0.5 * answer(3)
   endif

   ! Consider the core charge where zeta = 0 so f(zeta) = 0.
   if (coreRho > smallThresh) then
      rs = (3.0_double/(4.0_double*pi*coreRho))**0.33333333
      sqrtrs = sqrt(rs)

      if (rs >= 1.0) then

         ! answer(4) = ex + ec (core charge) (zeta = 0)
         ec = g / (1 + B1P*sqrtrs + B2P*rs)
         answer(4) = (-ex0 / rs + ec)
      else
         AP =  0.03110
         BP = -0.04800
         CP =  0.00200
         DP = -0.01160
         lnrs = log(rs)

         ! answer(4) = ex + ec (core charge) (zeta = 0)
         ec = AP*lnrs + BP + CP*lnrs + DP*rs
         answer(4) = (-ex0 / rs + ec)
      endif

      answer(4) = 0.5 * answer(4)
   endif
end subroutine ceperleyAlderSP


subroutine oldEXCORR(rh,sold,rhc,answer)

   ! Include kind definitions and the spin polarized vBH constants
   use O_Kinds

!   implicit none
   implicit real*8(a-h,o-z)

   ! Define dummy variables passed to this function and the return value.
   ! rh=+ totalRho, sold = spinDiffRho,rhc = coreRho
   real (kind=double) :: rh,sold,rhc
   real (kind=double), dimension (4) :: answer

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
      deallocate (potDifference)
   endif

end subroutine cleanUpPotentialUpdate


end module O_PotentialUpdate
