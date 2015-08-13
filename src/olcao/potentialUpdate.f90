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

   ! Define variables used to update the solid state potential.
   real (kind=double), allocatable, dimension (:,:) :: xl0,xl1,xl2 ! Holds the
         !   actually used potential coefficients from the current and previous
         !   two iterations.
   real (kind=double), allocatable, dimension (:,:) :: yl0,yl1,yl2 ! Holds the
         !   guesses for the next set of potential coefficients from the
         !   current and previous two iterations.

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
         & relaxFactor, xcCode, converged, convgTest
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
   real (kind=double) :: spinDiffSum
   real (kind=double) :: nonCoreTotalSum
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
   real (kind=double), allocatable, dimension (:,:) :: outputPot
   real (kind=double), allocatable, dimension (:,:) :: realSpacePotDiff
   real (kind=double), allocatable, dimension (:)   :: averageDelta
   real (kind=double), allocatable, dimension (:)   :: maxDelta
   real (kind=double), allocatable, dimension (:)   :: typesMagneticMoment

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
   call solveDPOSVX (potDim,spin+1,tempOverlap,potDim,generalRho(:,:spin+1))

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
   call solveDPOSVX (potDim,2,tempOverlap,potDim,elecStatPot(:,:2))


   ! Begin the exchange-correlation potential fitting

   ! The maxNumRayPoints was determined already in the access to the setup HDF5
   !   data.  Here we just copy the value from points(1).
   maxNumRayPoints = points(1)

   ! Allocate space to hold the exchange-correlation potential, exchange
   !   correlation radial weights, Rho Matrix Operator, and the resultant Rho.
   allocate (exchCorrPot  (potDim,2+spin))
   allocate (exchCorrRho  (1+spin,maxNumRayPoints))
   allocate (radialWeight (maxNumRayPoints))
   allocate (exchRhoOp    (potDim,maxNumRayPoints))


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
            & exchRhoOp(:,:),potPoints,hdferr)
      if (hdferr /= 0) stop 'Failed to read exch rho operator'


      ! The exchange correlation matrix times the charge vector produces the
      !   the real space charge
      do j = 1, numRayPoints
         coreSum         = sum(exchRhoOp(:potDim,j) * generalRho(:potDim,3))
         nonCoreTotalSum = sum(exchRhoOp(:potDim,j) * generalRho(:potDim,1))

         if (spin == 2) then
            spinDiffSum = sum(exchRhoOp(:potDim,j) * generalRho(:potDim,8))
            exchCorrRho (3,j) = spinDiffSum  ! Valence spin difference
         endif

         if (coreSum < 0.0_double) then
            coreSum = 0.0_double
         endif
         if (nonCoreTotalSum < smallThresh) then
! This is totally ridiculous.  If the following statement is not present, then
!   the program will not work.  The values for nonCoreTotalSum will not be
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
!write (20,*) "nonCoreTotalSum=",nonCoreTotalSum
!call flush (20)
            nonCoreTotalSum = smallThresh
!write (20,*) "nonCoreTotalSum=",nonCoreTotalSum
!call flush (20)
         endif

         exchCorrRho (1,j) = nonCoreTotalSum
         exchCorrRho (2,j) = coreSum
      enddo

      ! Define the exchange correlation functions
      ! 100 = Wigner
      ! 101 = Ceperley-Alder
      ! 102 = Hedin-Lundqvist
      ! 150 = Ceperley-Alder
      ! 151 = von Barth-Hedin
      ! 152 = unknown

      if (xcCode == 100) then
            do j = 1, numRayPoints
               ! These functions have been converted to subroutines because the
               !   return value of the second array index value was 0.0 on the
               !   ia64 HP-UX machine sirius.  This seems to work.
               call wignerXC(exchCorrRho(1,j),currentExchCorrPot(1:2))
               call wignerXCEnergy(exchCorrRho(2,j),currentExchCorrPot(3))
               do k = 1,3
                  exchCorrPot(:potDim,k) = exchCorrPot(:potDim,k) + &
                        & radialWeight(j) * currentExchCorrPot(k) * &
                        & exchRhoOp(:potDim,j)
               enddo
            enddo
      elseif (xcCode == 101) then
            do j = 1, numRayPoints
               call ceperleyAlderXC(exchCorrRho(1,j),currentExchCorrPot(1:2))
               call ceperleyAlderXCEnergy(exchCorrRho(2,j),&
                     & currentExchCorrPot(3))
               do k = 1,3
                  exchCorrPot(:potDim,k) = exchCorrPot(:potDim,k) + &
                        & radialWeight(j) * currentExchCorrPot(k) * &
                        & exchRhoOp(:potDim,j)
               enddo
            enddo
      elseif (xcCode == 102) then
            do j = 1, numRayPoints
               call hedinLundqvistXC(exchCorrRho(1,j),&
                     & currentExchCorrPot(1:2))
               call hedinLundqvistXCEnergy(exchCorrRho(2,j),&
                     & currentExchCorrPot(3))
               do k = 1,3
                  exchCorrPot(:potDim,k) = exchCorrPot(:potDim,k) + &
                        & radialWeight(j) * currentExchCorrPot(k) * &
                        & exchRhoOp(:potDim,j)
               enddo
            enddo
      elseif (xcCode == 150) then
            ! Ceperley and Alder Exchange-Correlation (LSDA)
            do j = 1, numRayPoints
               call ceperleyAlderSP(exchCorrRho(1,j),exchCorrRho(3,j),&
                     & exchCorrRho(2,j),currentExchCorrPot(1:4))
               do k = 1,4
                  exchCorrPot(:potDim,k) = exchCorrPot(:potDim,k) + &
                        & radialWeight(j) * currentExchCorrPot(k) * &
                        & exchRhoOp(:potDim,j)
               enddo
            enddo

      elseif (xcCode == 151) then
            ! von Barth and Hedin Exchange-Correlation (LSDA)
            do j = 1, numRayPoints
               call vonBarthHedin(exchCorrRho(1,j),exchCorrRho(3,j),&
                     & exchCorrRho(2,j),currentExchCorrPot(1:4))
               do k = 1,4
                  exchCorrPot(:potDim,k) = exchCorrPot(:potDim,k) + &
                        & radialWeight(j) * currentExchCorrPot(k) * &
                        & exchRhoOp(:potDim,j)
               enddo
            enddo
      elseif (xcCode == 152) then
            ! Old SPmain.for Exchange Correlation (LSDA) (Should be like the
            !   von barth and Hedin above.)
            do j = 1, numRayPoints
               call oldEXCORR(exchCorrRho(1,j),exchCorrRho(3,j),&
                     & exchCorrRho(2,j),currentExchCorrPot(1:4))
               do k = 1,4
                  exchCorrPot(:potDim,k) = exchCorrPot(:potDim,k) + &
                        & radialWeight(j) * currentExchCorrPot(k) * &
                        & exchRhoOp(:potDim,j)
               enddo
            enddo
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
         & exchCorrPot(:,:2+spin))


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
   deallocate (exchCorrOverlap)

   ! Allocate space to hold the output potentials
   allocate (outputPot (potDim,spin))
   if (.not.allocated(xl0)) then
      allocate (xl0 (potDim,spin))
      allocate (xl1 (potDim,spin))
      allocate (xl2 (potDim,spin))
      allocate (yl0 (potDim,spin))
      allocate (yl1 (potDim,spin))
      allocate (yl2 (potDim,spin))

      do i = 1, spin
         xl1(:,i) = potCoeffs(:,i)
         xl2(:,i) = potCoeffs(:,i)
         yl1(:,i) = elecStatPot(:,1) + elecStatPot(:,2) + exchCorrPot(:,i)
         yl2(:,i) = elecStatPot(:,1) + elecStatPot(:,2) + exchCorrPot(:,i)
      enddo
   endif


   do i = 1, spin
      xl0(:,i) = potCoeffs(:,i)
      yl0(:,i) = elecStatPot(:,1) + elecStatPot(:,2) + exchCorrPot(:,i)

      ! Form the difference between last iteration and the current new
      !   calculation of the potential.  (Note that this calculation will be
      !   mixed later with the potential of last iteration to make the
      !   convergence more smooth.)
      outputPot(:,i) = yl0(:,i) - xl0(:,i)

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
            & exchRhoOp(:,:),potPoints,hdferr)
      if (hdferr /= 0) stop 'Failed to read exch rho operator #2'


      ! Calculate the exchange-correlation rho matrix operator times the
      !   potential vector difference to get the difference in real space.
      do j = 1, spin
         do k = 1, numRayPoints
            realSpacePotDiff(k,j) = sum(exchRhoOp(:,k) * outputPot(:,j))
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
   deallocate (outputPot)

   ! Generate the potentials for the next iteration using the method D.G.
   !   Anderson, J. Assoc. Comp. Mach. 12, 547 (1965).

   ! The parameter here referred to as the feedbackLevel is called M by
   !   Anderson.  If M=0, the straight 'relaxed' iteration is used.  If M=2,
   !   the 2x2 determinant in anderson's method is checked.  If it is too
   !   small, the first order M=1 extrapolation is used instead.  M is
   !   effectively reduced to 0 or 1 if the iteration is 0, 1, or 2.

   allocate ( rl0 (potDim))
   allocate ( rl1 (potDim))
   allocate ( rl2 (potDim))
   allocate (drl1 (potDim))
   allocate (drl2 (potDim))

   ! This is basically copied from the old code with little modification for
   !   the names since I don't know what they mean or do in most cases.
   do i = 1, spin
      th1 = 0.0_double
      th2 = 0.0_double
      if (feedbackLevel /= 0) then
         if (currIteration > 1) then
            rl0(:)  = yl0(:,i) - xl0(:,i)
            rl1(:)  = yl1(:,i) - xl1(:,i)
            rl2(:)  = yl2(:,i) - xl2(:,i)
            drl1(:) = rl0(:) - rl1(:)
            drl2(:) = rl0(:) - rl2(:)

            ! This part is new.  Make temp matrices holding the potAlphaOverlap
            !   times the drl1, and drl2 vectors.

            ! Initialize the summation variables
            s11 = 0.0_double
            s12 = 0.0_double
            s22 = 0.0_double
            t1  = 0.0_double
            t2  = 0.0_double

            ! The cases for the drl1 vector are done first. 
            do j = 1, potDim
               tempOverlap(:,j) = potAlphaOverlap(:,j) * drl1(j)
               s11 = s11 + sum(tempOverlap(:,j) * drl1(:))
               t1  = t1  + sum(tempOverlap(:,j) * rl0(:))
            enddo

            ! The cases for the drl2 vector are done second.
            do j = 1, potDim
               tempOverlap(:,j) = potAlphaOverlap(:,j) * drl2(j)
               s12 = s12 + sum(tempOverlap(:,j) * drl1(:))
               s22 = s22 + sum(tempOverlap(:,j) * drl2(:))
               t2  = t2  + sum(tempOverlap(:,j) * rl0(:))
            enddo
            th1 = t1/s11
            if (feedbackLevel /= 1) then
               if (currIteration > 2) then

                  ! Calculate the determinate.
                  det = s11*s22 - s12*s12
                  if (det/(s11*s22) >= 0.01_double) then
                     th1 = ( s22*t1 - s12*t2)/det
                     th2 = (-s12*t1 + s11*t2)/det
                  endif
               endif
            endif
         endif
      endif

      ! Form the next iteration and save the previous two.
      potCoeffs(:,i) = &
            & (1.0_double - relaxFactor) * &
            & ((1.0_double-th1-th2)*xl0(:,i) + th1*xl1(:,i) + th2*xl2(:,i)) + &
            & relaxFactor * &
            & ((1.0_double-th1-th2)*yl0(:,i) + th1*yl1(:,i) + th2*yl2(:,i))

      xl2(:,i) = xl1(:,i)
      xl1(:,i) = xl0(:,i)
      yl2(:,i) = yl1(:,i)
      yl1(:,i) = yl0(:,i)
   enddo


   ! Record the potential terms (alphas), potential coefficients, total charge
   !   density (core+valence), and valence charge density (up+down), and for
   !   spin polarized cases the difference (up-down).
   rewind (8)

   ! Write the file header that says the total number of types.
   write (8,fmt="(a9,i5)") "NUM_TYPES",numPotTypes

   do i = 1, spin

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

            ! Write the term.
            if (spin == 1) then
               write (8,fmt="(2(1x,e16.10),3(1x,e12.6))") &
                     & potCoeffs(potTermCount,i),potAlphas(potTermCount),&
                     & generalRho(potTermCount,1),generalRho(potTermCount,2),&
                     & 0.0_double
            else
               write (8,fmt="(2(1x,e16.10),3(1x,e12.6))") &
                     & potCoeffs(potTermCount,i),potAlphas(potTermCount),&
                     & generalRho(potTermCount,1),generalRho(potTermCount,2),&
                     & generalRho(potTermCount,8)
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
   deallocate (rl0)
   deallocate (rl1)
   deallocate (rl2)
   deallocate (drl1)
   deallocate (drl2)

   ! Log the date and time we end.
   call timeStampEnd (18)

end subroutine makeSCFPot


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
   alpha = (4.0_double/(9.0_double*pi))**0.33333333
   rs = (3.0_double/(4.0_double*pi*totalRho))**0.33333333
   a = 2.0_double**(-1.0_double/3.0_double)
   gamma = 4.0_double/3.0_double * a / (1.0_double - a)
   two13 = 2.0**(0.3333333333)

   rP = 21.0
   rF = 2.0**(4.0/3.0)*rP
   cP = 0.045
   cF = cP/2.0

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

      epsxP = -3.0/(2.0 * pi * alpha) / rs

      epsxF = two13 * epsxP

      muxP = 4.0/3.0 * epsxP

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

   if (allocated (xl0)) then
      deallocate (xl0)
      deallocate (xl1)
      deallocate (xl2)
      deallocate (yl0)
      deallocate (yl1)
      deallocate (yl2)
   endif

end subroutine cleanUpPotentialUpdate


end module O_PotentialUpdate
