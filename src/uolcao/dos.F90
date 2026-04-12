module O_DOS

   ! Import the necessary modules.
   use O_Kinds

   ! Make sure that no variables are declared accidentally.
   implicit none

   ! Define module data specifically for TDOS calculation for each SCF cycle.
   real (kind=double), allocatable, dimension (:,:,:) :: tdos
   real (kind=double), allocatable, dimension (:)     :: currentEnergyValues

   ! Define module data for the general PDOS/TDOS calculation.
   integer :: numEnergyPoints
   real (kind=double), allocatable, dimension (:)     :: energyScale

   contains

subroutine computeIterationTDOS

   ! Import the necessary modules.
   use O_Kinds
   use O_Constants,       only: pi, hartree
   use O_Populate,        only: occupiedEnergy
   use O_SecularEquation, only: energyEigenValues
   use O_KPoints,         only: numKPoints, kPointWeight
   use O_Potential,       only: spin, lastIteration, currIteration
   use O_Input, only: sigmaDOS, eminDOS, emaxDOS, deltaDOS, numStates

   ! Make sure that no variables are declared accidentally.
   implicit none

   ! Define local variables
   integer :: i,j,k,l
   real (kind=double) :: sigmaSqrtPi
   real (kind=double) :: expTerm
   real (kind=double) :: expFactor

   ! Normalization for Gaussian broadening.
   sigmaSqrtPi = sqrt(pi) * sigmaDOS

   if (.not. allocated(tdos)) then

      ! Determine the number of energy buckets to be computed for.
      numEnergyPoints = int((emaxDOS - eminDOS ) / deltaDOS)

      ! Allocate memory to hold the data for those points. The funny-stuff
      !   with last iteration: If lastIteration is zero, that is a special
      !   case. We still just compute one iteration, but it is also a signal
      !   to not re-compute the wave function later.
      if (lastIteration > 0) then
         allocate (tdos(numEnergyPoints,spin,lastIteration))
      else
         allocate (tdos(numEnergyPoints,spin,1))
      endif
      allocate (energyScale(numEnergyPoints))
      allocate (currentEnergyValues(numStates))

      ! Initialize the data for accumulation.
      tdos(:,:,:) = 0.0_double

      ! Compute the values for the energy scale.
      do i = 1, numEnergyPoints
         energyScale(i) = eminDOS + (i-1) * deltaDOS
      enddo
   endif


   do i = 1, spin
      do j = 1, numKPoints

         ! Obtain a copy of the current energy values for this spin and kpoint,
         !   and shift them in accordance  with the highest occupied energy
         !   level.
         currentEnergyValues(:) = (energyEigenValues(:,j,i) - occupiedEnergy)

         do k = 1, numStates

            ! Apply broadening to each state.
            do l = 1, numEnergyPoints

               ! Compute the exponential term for the broadening of this point.
               !   Note that from the usual Gaussian term of exp(-alpha*x^2)
               !   we are have (eV-eS)^2 = x^2 and (1/sigma)^2 = alpha.
               expTerm = ((currentEnergyValues(k)-energyScale(l))/sigmaDOS)**2

               ! If the exponential term is less than 50 we apply the
               !   broadening.
               if (expTerm < 50.0_double) then

                  ! Compute the exponential factor.  It is at this point that
                  !   the kpoint weighting factor is applied to make the energy
                  !   bucket we are considering be filled with the right number
                  !   of electrons.  Note that we must account for the value
                  !   of "spin".  This will put either 1 or 2 electrons in each
                  !   state.
                  expFactor = exp(-expTerm) / sigmaSqrtPi / hartree * &
                        & kPointWeight(j) / real(spin,double)

                  ! Store the broadened TDOS.
                  tdos(l,i,currIteration) = tdos(l,i,currIteration) + expFactor
               endif
            enddo
         enddo
      enddo
   enddo

end subroutine computeIterationTDOS

subroutine printIterationTDOS

   ! Import the necessary modules.
   use O_Kinds
   use O_Constants, only: hartree
   use O_Potential, only: spin, currIteration

   ! Make sure that no variables are declared accidentally.
   implicit none

   ! Define local variables
   integer :: i,j
!   integer :: pointCount
   integer :: fileID
   character*9 :: fileName

   ! Note that the current iteration was already incremented by one when this
   !   subroutine is called so we must refer back to the previous iteration
   !   when accessing data.

   do i = 1, spin

      ! Define the file ID to open.
      fileID = currIteration+i*1000

      ! Define the file name.
      write (fileName,fmt="(a5,i4)") "fort.",fileID

      ! Open the file for creating the opendx output.
      open (unit=fileID,file=fileName,status='unknown',form='formatted')


do j = 1, numEnergyPoints
! Make sure to convert to eV upon output.
write (fileID,fmt="(1x,3f12.8)") energyScale(j)*hartree,tdos(j,1,1),&
      & tdos(j,1,currIteration-1)

enddo

   enddo


   ! Deallocate unused memory.
   deallocate (tdos)
   deallocate (energyScale)
   deallocate (currentEnergyValues)

end subroutine printIterationTDOS


subroutine computeDOS(inSCF)

   ! Import the necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Potential,   only: spin
   use O_Populate,    only: electronPopulation
   use O_KPoints, only: numKPoints, kPointWeight, &
         & kPointIntgCode, numPointOps
   use O_Constants, only: pi, hartree, lAngMomCount
   use O_AtomicSites, only: valeDim, &
         & numAtomSites, atomSites, invAtomPerm
   use O_AtomicTypes, only: numAtomTypes, &
         & atomTypes, maxNumValeStates
   use O_Input, only: numStates, sigmaDOS, &
         & eminDOS, emaxDOS, deltaDOS, &
         & detailCodePDOS
#ifndef GAMMA
   use O_SecularEquation, only: valeValeOL, &
         & valeVale, energyEigenValues, &
         & readDataSCF, readDataPSCF
#else
   use O_SecularEquation, only: &
         & valeValeOLGamma, valeValeGamma, &
         & energyEigenValues, &
         & readDataSCF, readDataPSCF
#endif


   ! Make sure that no funny variables are used.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: inSCF

   ! Define local variables.
   integer :: h, i, j, k, l  ! Loop index variables.
   character*17 :: formatString
   character*1, dimension (lAngMomCount) :: &
         & QN_lLetter
   character*14, dimension (4,7) :: QN_mLetter
   integer, allocatable, dimension (:) :: cumulNumDOS
   integer, allocatable, dimension (:) :: pdosIndex
         !   Given a valence state, this stores the place within pdosAccum that
         !   this state should be stored.
   integer, allocatable, dimension (:) :: &
         & numAtomStates
   integer :: cumulDOSTotal
   integer :: numSQN_l
   integer :: numPQN_l
   integer :: numDQN_l
   integer :: numFQN_l
   integer :: initSIndex
   integer :: initPIndex
   integer :: initDIndex
   integer :: initFIndex
   integer :: currentType
   integer :: valeDimIndex
   integer :: initIndex
   integer :: finIndex
   integer :: energyLevelCounter
   integer :: numCols
   integer :: stateSpinKPointIndex
   real (kind=double) :: occupancyNumber
   real (kind=double) :: oneValeRealAccum
   real (kind=double) :: expTerm
   real (kind=double) :: expFactor
   real (kind=double) :: sigmaSqrtPi
   real (kind=double) :: integratedArea
   real (kind=double) :: numStatesInRange
   real (kind=double) :: totalElectronsComputed
   real (kind=double) :: electronFactor
   real (kind=double) :: currentPopulation
   real (kind=double), allocatable, &
         & dimension (:)   :: pdosAccum
   real (kind=double), allocatable, &
         & dimension (:)   :: localizationIndex
   real (kind=double), allocatable, &
         & dimension (:)   :: totalSystemDos
   real (kind=double), allocatable, &
         & dimension (:)   :: energyValuesAvg
   real (kind=double), allocatable, &
         & dimension (:,:) :: pdosComplete
   real (kind=double), allocatable, &
         & dimension (:,:) :: electronNumber
#ifndef GAMMA
   complex (kind=double), allocatable, &
         & dimension (:) :: waveFnSqrd
#else
   real (kind=double), allocatable, &
         & dimension (:) :: waveFnSqrdGamma
#endif

   ! LAT PDOS variables (used only when kPointIntgCode == 1).
   real (kind=double), allocatable, &
         & dimension(:,:,:) :: projArray
   integer, allocatable, &
         & dimension(:,:) :: channelPermTbl

   ! Log the date and time we start.
   call timeStampStart (19)

   ! Mode 3 restriction: per-atom per-lm PDOS with LAT requires D^l rotation
   !   matrices for the channel permutation, which are not available. Stop with
   !   a clear error (DESIGN 1.4).
   if (kPointIntgCode == 1 .and. &
         & detailCodePDOS == 3) then
      write (20, *) 'ERROR: LAT PDOS (kPointIntg'
      write (20, *) '  Code=1) with per-atom per-'
      write (20, *) '  lm detail (detailCodePDOS='
      write (20, *) '  3) requires D^l rotation'
      write (20, *) '  matrices that are not'
      write (20, *) '  available. Use mode 0, 1,'
      write (20, *) '  or 2 with LAT integration.'
      stop 'computeDOS: mode 3 + LAT unsupported'
   endif

   ! Normalization for Gaussian broadening.
   sigmaSqrtPi = sqrt(pi) * sigmaDOS

   ! Allocate arrays and matrices for this computation.
   allocate (pdosIndex     (valeDim))
   allocate (numAtomStates (numAtomSites))
   ! Store DOS for each type's orbital sum.
   if     (detailCodePDOS == 0) then
      allocate (cumulNumDOS (numAtomTypes + 1))
   ! Store DOS for each atom's TDOS.
   elseif (detailCodePDOS == 1) then
      allocate (cumulNumDOS (1))
   ! Store DOS for each QN_nl resolved atom.
   elseif (detailCodePDOS == 2) then
      allocate (cumulNumDOS (numAtomSites + 1))
   ! Store DOS for each QN_nlm resolved atom.
   elseif (detailCodePDOS == 3) then
      allocate (cumulNumDOS (numAtomSites + 1))
   endif

   ! The pdosAccum array and the waveFnSqrd/overlap arrays are only needed for
   !   the Gaussian path. The LAT path handles its own allocations inside
   !   computeProjections_LAT.
   if (kPointIntgCode /= 1) then
      if     (detailCodePDOS == 0) then
         allocate (pdosAccum (valeDim))
      elseif (detailCodePDOS == 1) then
         allocate (pdosAccum (numAtomSites))
      elseif (detailCodePDOS == 2) then
         allocate (pdosAccum (valeDim))
      elseif (detailCodePDOS == 3) then
         allocate (pdosAccum (valeDim))
      endif
#ifndef GAMMA
      allocate (waveFnSqrd(valeDim))
      allocate (valeValeOL(valeDim, valeDim))
      if (inSCF == 0) then
         allocate (valeVale(valeDim, numStates, 1))
      endif
#else
      allocate (waveFnSqrdGamma(valeDim))
      allocate (valeValeOLGamma(valeDim, valeDim))
      if (inSCF == 0) then
         allocate (valeValeGamma( &
               & valeDim, numStates, 1))
      endif
#endif
   endif

   ! Define the QN_l letters.
   QN_lLetter(1) = 's'
   QN_lLetter(2) = 'p'
   QN_lLetter(3) = 'd'
   QN_lLetter(4) = 'f'

   ! Define the QN_m resolved letters.
   QN_mLetter(1,1) = 'r'
   QN_mLetter(2,1) = 'x'
   QN_mLetter(2,2) = 'y'
   QN_mLetter(2,3) = 'z'
   QN_mLetter(3,1) = 'xy'
   QN_mLetter(3,2) = 'xz'
   QN_mLetter(3,3) = 'yz'
   QN_mLetter(3,4) = 'xx~yy'
   QN_mLetter(3,5) = '2zz~xx~yy'
   QN_mLetter(4,1) = 'xyz'
   QN_mLetter(4,2) = 'xxz~yyz'
   QN_mLetter(4,3) = 'xxx~3yyx'
   QN_mLetter(4,4) = '3xxy~yyy'
   QN_mLetter(4,5) = '2zzz~3xxz~3yyz'
   QN_mLetter(4,6) = '4zzx~xxx~yyx'
   QN_mLetter(4,7) = '4zzy~xxy~yyy'

   ! Initialize other variables.
   cumulNumDOS(:) = 0
   cumulDOSTotal  = 0

   if (detailCodePDOS == 0) then

      ! Initialize counter to index the cumulative sum of QN_l orbitals for all
      !   atomic types.  (This PDOS will give a DOS for each QN_nl pair of each
      !   atomic type.)
      cumulNumDOS(1) = 0

      ! Loop to record the number of orbitals that each type contributes.  (An
      !   orbital is just a QN_nl pair.)
      do i = 1, numAtomTypes
         cumulNumDOS(i+1) = cumulNumDOS(i) + &
               & sum(atomTypes(i)%numQN_lValeRadialFns(:))
      enddo

      ! Record the total number of orbitals summed over all types.
      cumulDOSTotal = cumulNumDOS(numAtomTypes+1)

   elseif (detailCodePDOS == 1) then

      ! There will be one DOS curve for each atom and that is it.
      cumulDOSTotal = numAtomSites

   elseif (detailCodePDOS == 2) then

      ! Initialize counter to index the cumulative sum of QN_l orbitals for all
      !   atoms.  (This PDOS will give a DOS for each QN_nl pair of each atom.)
      cumulNumDOS(1) = 0

      ! Loop to record the number of orbitals that for each atom contributes.
      !   (An orbital is just a QN_nl pair.)
      do i = 1, numAtomSites
         cumulNumDOS(i+1) = cumulNumDOS(i) + sum( &
               & atomTypes(atomSites(i)%atomTypeAssn)%numQN_lValeRadialFns(:))
      enddo

      ! Record the total number of orbitals summed over all atoms.
      cumulDOSTotal = cumulNumDOS(numAtomSites+1)

   elseif (detailCodePDOS == 3) then

      ! Initialize counter to index the cumulative sum of QN_l orbitals for
      !   all atoms.  (This PDOS will give a DOS for each QN_nlm set for each
      !   atom.)
      cumulNumDOS(1) = 0

      ! Loop to record the index number for each atom's orbitals.
      do i = 1, numAtomSites
         cumulNumDOS(i+1) = cumulNumDOS(i) + &
            & atomTypes(atomSites(i)%atomTypeAssn)%numQN_lValeRadialFns(1)*1 + &
            & atomTypes(atomSites(i)%atomTypeAssn)%numQN_lValeRadialFns(2)*3 + &
            & atomTypes(atomSites(i)%atomTypeAssn)%numQN_lValeRadialFns(3)*5 + &
            & atomTypes(atomSites(i)%atomTypeAssn)%numQN_lValeRadialFns(4)*7 
      enddo

      ! Record the total number of QN_m resolved orbitals summed over all atoms.
      cumulDOSTotal = cumulNumDOS(numAtomSites+1)

   endif

   ! Initialize valeDimIndex to record the index number in valeDim that each
   !   pdos state is at.
   valeDimIndex = 0

   ! Loop over every atom in the system to index where the pdos values for
   !   each atom should be stored.
   do i = 1, numAtomSites

      ! Obtain the type of the current atom.
      currentType = atomSites(i)%atomTypeAssn

      ! Identify and store the number of valence states for this atom.
      numAtomStates(i) = atomTypes(currentType)%numValeStates

      ! In the case where the PDOS should be collected by types we loop
      !   through each QN_nl pair for this atom and record the index where its
      !   PDOS should be recorded according to its type.  In the case where
      !   the PDOS is collected by atoms we link each QN_nl pair of this atom
      !   with the index of this atom.  In the case where the PDOS is collected
      !   for each QN_nl pair of each atom the index is as in the first case
      !   that now we index for each atom as opposed to each type.
      if (detailCodePDOS == 0) then

         numSQN_l = atomTypes(currentType)%numQN_lValeRadialFns(1)
         numPQN_l = atomTypes(currentType)%numQN_lValeRadialFns(2)
         numDQN_l = atomTypes(currentType)%numQN_lValeRadialFns(3)
         numFQN_l = atomTypes(currentType)%numQN_lValeRadialFns(4)

         initSIndex = cumulNumDOS(currentType)
         initPIndex = cumulNumDOS(currentType) + numSQN_l
         initDIndex = cumulNumDOS(currentType) + numSQN_l + numPQN_l
         initFIndex = cumulNumDOS(currentType) + numSQN_l + numPQN_l + numDQN_l

         do j = 1, numSQN_l
            valeDimIndex = valeDimIndex + 1
            pdosIndex(valeDimIndex) = initSIndex + j
         enddo
         do j = 1, numPQN_l
            do k = 1,3
               valeDimIndex = valeDimIndex + 1
               pdosIndex(valeDimIndex) = initPIndex + j
            enddo
         enddo
         do j = 1, numDQN_l
            do k = 1,5
               valeDimIndex = valeDimIndex + 1
               pdosIndex(valeDimIndex) = initDIndex + j
         enddo
         enddo
         do j = 1, numFQN_l
            do k = 1,7
               valeDimIndex = valeDimIndex + 1
               pdosIndex(valeDimIndex) = initFIndex + j
            enddo
         enddo
      elseif (detailCodePDOS == 1) then
         do j = 1, numAtomStates(i)
            valeDimIndex = valeDimIndex + 1
            pdosIndex(valeDimIndex) = i  ! NOTE THAT THIS IS 'i', NOT 'j'.
         enddo
      elseif (detailCodePDOS == 2) then ! Consider spdf for each atom too.

         numSQN_l = atomTypes(currentType)%numQN_lValeRadialFns(1)
         numPQN_l = atomTypes(currentType)%numQN_lValeRadialFns(2)
         numDQN_l = atomTypes(currentType)%numQN_lValeRadialFns(3)
         numFQN_l = atomTypes(currentType)%numQN_lValeRadialFns(4)

         initSIndex = cumulNumDOS(i)
         initPIndex = cumulNumDOS(i) + numSQN_l
         initDIndex = cumulNumDOS(i) + numSQN_l + numPQN_l
         initFIndex = cumulNumDOS(i) + numSQN_l + numPQN_l + numDQN_l

         do j = 1, numSQN_l
            valeDimIndex = valeDimIndex + 1
            pdosIndex(valeDimIndex) = initSIndex + j
         enddo
         do j = 1, numPQN_l
            do k = 1,3
               valeDimIndex = valeDimIndex + 1
               pdosIndex(valeDimIndex) = initPIndex + j
            enddo
         enddo
         do j = 1, numDQN_l
            do k = 1,5
               valeDimIndex = valeDimIndex + 1
               pdosIndex(valeDimIndex) = initDIndex + j
         enddo
         enddo
         do j = 1, numFQN_l
            do k = 1,7
               valeDimIndex = valeDimIndex + 1
               pdosIndex(valeDimIndex) = initFIndex + j
            enddo
         enddo
      elseif (detailCodePDOS == 3) then
         do j = 1, numAtomStates(i)
            valeDimIndex = valeDimIndex + 1
            pdosIndex(valeDimIndex) = valeDimIndex ! Each QN_nlm is saved.
         enddo
      endif
   enddo

   ! Determine the number of energy buckets to be computed for.
   numEnergyPoints = int((emaxDOS - eminDOS ) / deltaDOS)

   ! Allocate space to hold the pdos and localization index results
   allocate (localizationIndex (numStates))
   allocate (energyScale       (numEnergyPoints))
   allocate (totalSystemDos    (numEnergyPoints))
   allocate (electronNumber    (maxNumValeStates,numAtomSites))
   allocate (energyValuesAvg   (numStates))
   allocate (pdosComplete      (cumulDOSTotal,numEnergyPoints))

   ! Assign values to the energy scale.
   do i = 1, numEnergyPoints
      energyScale(i) = eminDOS + (i-1) * deltaDOS
   enddo

   ! For the LAT path, build the channel permutation table before the spin loop
   !   (spin-independent). This maps PDOS channel indices through the inverse
   !   atom permutation for IBZ unfolding during tetrahedron corner assembly
   !   (DESIGN 1.4, PSEUDOCODE 8.1).
   if (kPointIntgCode == 1) then
      call buildChannelPermTable(detailCodePDOS, &
            & numPointOps, cumulDOSTotal, &
            & cumulNumDOS, numAtomSites, &
            & invAtomPerm, channelPermTbl)
   endif

   do h = 1, spin

      ! Track the stateSpinKPoint index number.
      stateSpinKPointIndex = (h-1) * numStates

      ! Record which calculation is being done.
      if (spin == 2) then
         if (h == 1) then
            write (20,*) "Computing spin up DOS."
         else
            write (20,*) "Computing spin down DOS."
         endif
      endif

      ! Initialize various arrays and matrices.
      electronNumber  (:,:) = 0.0_double
      pdosComplete    (:,:) = 0.0_double
      localizationIndex (:) = 0.0_double

      ! Branch on integration method. The LAT path uses a two-pass design
      !   (project then integrate via tetrahedra); the Gaussian path uses a
      !   single-pass design (project and broaden simultaneously). Both fill
      !   pdosComplete, electronNumber, and localizationIndex. The output phase
      !   that follows is shared.
      if (kPointIntgCode == 1) then

         ! ----------------------------------------- LAT two-pass PDOS (DESIGN
         ! 1.4). -----------------------------------------

         ! Pass 1: stream IBZ k-points, compute Mulliken projections, store in
         !   projArray. Also accumulates electronNumber and localizationIndex
         !   for the diagnostic.
         write (20, *) "LAT PDOS Pass 1: projections"
         call computeProjections_LAT(inSCF, h, &
               & numKPoints, numStates, &
               & numAtomSites, numAtomStates, &
               & pdosIndex, valeDim, &
               & cumulDOSTotal, spin, projArray, &
               & electronNumber, localizationIndex)

         ! Pass 2: tetrahedron integration with Bloechl corner weights and
         !   channel permutation for IBZ unfolding.
         write (20, *) "LAT PDOS Pass 2: integration"
         call integratePDOS_LAT(projArray, &
               & channelPermTbl, pdosComplete, h, &
               & numStates, cumulDOSTotal, &
               & numEnergyPoints, energyScale)

         ! Free the projection array (large memory).
         deallocate (projArray)

      else

      ! ----------------------------------------- Gaussian single-pass PDOS
      ! (original code). -----------------------------------------

      ! Begin accumulating the DOS values
      do i = 1, numKPoints

         ! Track the stateSpinKPoint index number.
         if ((spin == 2) .and. (i /= 1)) then
            stateSpinKPointIndex = stateSpinKPointIndex + numStates
         endif

         ! Determine if we are doing the DOS in a post-SCF calculation, or
         !   within an SCF calculation. 
         if (inSCF == 1) then
            ! Read necessary data from SCF (setup,main) data structures.
            call readDataSCF(h,i,numStates,1) ! 1 = Overlap matrixCode
         else
            ! Read necessary data from post SCF (intg,band) data structures.
            call readDataPSCF(h,i,numStates,1) ! 1 = Overlap matrixCode
         endif


         do j = 1, numStates

            ! Track the stateSpinKPoint index number.
            stateSpinKPointIndex = stateSpinKPointIndex + 1

!            occupancyNumber = electronPopulation(stateSpinKPointIndex)

            ! Determine the occupancy number for this state.  The default
            !   assumption is that each state contains 2 electrons.  This is
            !   because the sum of all kpoint weighting factors is equal to 2.
            !   This will make spin polarized states half occupied and spin
            !   non-polarized states fully occupied.
            if (energyEigenValues(j,i,h) <= 0.0_double) then
               occupancyNumber = 1.0_double/real(spin,double)
            else
               occupancyNumber = 0.0_double
            endif

            ! Initialize the accumlator for the wave function--overlap
            !   interaction.
            pdosAccum(:) = 0.0_double

            ! Initialize a counter that tracks the index of valeDim (total
            !   number of states used in the system).
            valeDimIndex = 0

            ! Begin a loop over all atoms
            do k = 1, numAtomSites

               ! Loop over all the valence states for this atom
               do l = 1, numAtomStates(k)

                  ! Increment the valeDimIndex
                  valeDimIndex = valeDimIndex + 1

#ifndef GAMMA
                  ! Compute the square of the wave function for each element.
                  waveFnSqrd(:valeDim) = &
                        & conjg(valeVale(valeDimIndex,j,1)) * &
                        & valeVale(:valeDim,j,1)

                  ! Compute the effects of overlap for the real part only
                  !   (real*real) + (imag*imag).
                  oneValeRealAccum = sum( &
                       & real(waveFnSqrd(:valeDim),double) * &
                       & real(valeValeOL(:valeDim,valeDimIndex),double) + &
                       & aimag(waveFnSqrd(:valeDim)) * &
                       & aimag(valeValeOL(:valeDim,valeDimIndex)))
!write (22, *) "h, i, j, k, l"
!write (22, *) h, i, j, k, l
!write (22, *) "oneValeRealAccum"
!write (22, *) oneValeRealAccum
!!write (22, *) "waveFnSqrd, valeValeOL"
!!do m = 1, valeDim
!!write (22, *) waveFnSqrd(m), valeValeOL
!!enddo
#else
                  ! Compute the square of the wave function for each element.
                  waveFnSqrdGamma(:valeDim) = &
                        & valeValeGamma(valeDimIndex,j,1) * &
                        & valeValeGamma(:valeDim,j,1)

                  ! Compute the effects of overlap.
                  oneValeRealAccum = sum(waveFnSqrdGamma(:valeDim) * &
                        & valeValeOLGamma(:valeDim,valeDimIndex))
#endif

                  ! Store the current electron number assignment.
                  electronNumber(l,k) = electronNumber(l,k) + &
                        & oneValeRealAccum * kPointWeight(i) * occupancyNumber
!                  electronNumber(l,k) = electronNumber(l,k) + &
!                        & oneValeRealAccum * occupancyNumber


                  pdosAccum(pdosIndex(valeDimIndex)) = &
                        & pdosAccum(pdosIndex(valeDimIndex)) + &
                        & oneValeRealAccum / real(spin,double)


                  ! Store the square of the accumulation as the localization
                  !   index for this state.
                  localizationIndex(j) = localizationIndex(j) + &
                        & oneValeRealAccum * oneValeRealAccum * &
                        & kPointWeight(i) / real(spin,double)
!                  localizationIndex(j) = localizationIndex(j) + &
!                        & oneValeRealAccum * oneValeRealAccum * & &
!                        occupancyNumber
               enddo
            enddo


            ! Apply broadening to the pdos values according to the given
            !   pdos sigma factor.
            do k = 1, numEnergyPoints

               ! Compute the exponential term for the broadening of this point.
               !   Note that from the usual Gaussian term of exp(-alpha*x^2)
               !   we are have (eV-eS)^2 = x^2 and (1/sigma)^2 = alpha.
               expTerm = ((energyEigenValues(j,i,h) - &
                     & energyScale(k))/sigmaDOS)**2

               ! If the exponential term is less than 50 we apply the
               !   broadening.
               if (expTerm < 50.0_double) then

                  ! Compute the exponential factor.  It is at this point that
                  !   the kpoint weighting factor is applied to make the energy
                  !   bucket we are considering be filled with the right number
                  !   of electrons.  The pdosAccum below already has the
                  !   1/spin factor included.  Note that the hartree conversion
                  !   must be included here because the sigmaDOS (that was used
                  !   to make sigmaSqrtPi) is in units of hartree and the
                  !   exponential and kPointWeights are unitless. Thus, if we
                  !   want to get States / [eV Cell] then we need to convert
                  !   the sigmaSqrtPi back to eV. Because it is in the
                  !   denominator we need to divide the whole expression by
                  !   eV/hartree to get a final energy unit of eV. Recall that
                  !   the hartree variable equals 27.211... eV/hartree.
                  expFactor = exp(-expTerm) / sigmaSqrtPi / hartree * &
                        & kpointWeight(i)

                  ! Broaden and store the complete pdos
                  pdosComplete(1:cumulDOSTotal,k) = &
                        & pdosComplete(1:cumulDOSTotal,k) &
                        & + pdosAccum(1:cumulDOSTotal) * expFactor
               endif
            enddo
         enddo ! (j numStates)

         ! Record that this kpoint has been finished.
         if (mod(i,10) .eq. 0) then
            write (20,ADVANCE="NO",FMT="(a1)") "|"
         else
            write (20,ADVANCE="NO",FMT="(a1)") "."
         endif
         if (mod(i,50) .eq. 0) then
            write (20,*) " ",i
         endif
         call flush (20)

      enddo ! (i numKPoints)

      ! Reset the line pointer to the beginning of the next line if necessary.
      if (mod(numKPoints,10) .ne. 0) then
         write (20,*)
         call flush (20)
      endif

      endif ! (kPointIntgCode branch)


      ! Compute the total DOS for the whole system.
      totalSystemDos(:) = 0.0_double
      do i = 1, numEnergyPoints
         totalSystemDos(i) = totalSystemDos(i) + sum(pdosComplete(:,i))
      enddo


      ! Compute the number of electron states present in the range of eminDOS
      !   to emaxDOS for this spin orientation or combination.
      numStatesInRange = 0.0_double
      do i = 1, numKPoints
         do j = 1, numStates
            if ((energyEigenValues(j,i,h) > eminDOS) .and. &
                  & (energyEigenValues(j,i,h) < emaxDOS)) then
               numStatesInRange = numStatesInRange + kPointWeight(i) / &
                     & real(spin,double)
            endif
         enddo
      enddo


      ! Compute the total number of electrons determined from the pdos
      !   calculation.
      totalElectronsComputed = sum(electronNumber(:,:))


      ! Compute the total number of electrons in each spin direction (up,down).
      !   This number is then be compared to the value for the
      !   totalElectronsComputed to normalize the spectra.
      currentPopulation = 0.0_double
      energyLevelCounter = 0
      do i = 1, numKPoints
         do j = 1, spin
            do k = 1, numStates
               energyLevelCounter = energyLevelCounter + 1
               ! Accumulate only electrons for the current spin orientation.
               !   In the case of spin non-polarized calculations we will
               !   accumulate all populations.
               if (j == h) then
                  currentPopulation = currentPopulation + &
                        & electronPopulation(energyLevelCounter)
               endif
            enddo
         enddo
      enddo


      ! Compute the ratio of the exact electron count to the computed electron
      !   count. For the Gaussian path this corrects broadening-tail truncation.
      !   For LAT, tetrahedron corner weights provide exact BZ integration so
      !   electronFactor should be ~1.0 (DESIGN 1.4).
      electronFactor = currentPopulation &
            & / totalElectronsComputed

      if (kPointIntgCode == 1) then
         ! LAT path: log the ratio as a diagnostic but do NOT apply it. A ratio
         !   far from 1.0 signals an integration bug.
         write (20, *) 'LAT electronFactor: ', &
               & electronFactor
      else
         ! Gaussian path: apply the correction to electronNumber, TDOS, and
         !   PDOS.
         electronNumber(:,:) = &
               & electronNumber(:,:) * electronFactor
         totalSystemDos(:) = &
               & totalSystemDos(:) * electronFactor
         pdosComplete(:,:) = &
               & pdosComplete(:,:) * electronFactor
      endif


      ! Find the total area under the computed TDOS. This should equal the
      !   number of electron spin states available in the system (occupied +
      !   unoccupied) OVER THE REQUESTED ENERGY RANGE ONLY.
      integratedArea = sum( &
            & totalSystemDos(1:numEnergyPoints - 1) &
            & + totalSystemDos(2:numEnergyPoints)) &
            & * deltaDOS * 0.5_double * hartree

      ! Record the exact and calculated values for electrons and states.
      write (20, *) 'Electrons Calculated:    ', &
            & totalElectronsComputed
      write (20, *) 'Electrons Expected:      ', &
            & currentPopulation
      write (20, *) 'Spin States Calculated:  ', &
            & integratedArea
      write (20, *) 'Spin States Expected:    ', &
            & numStatesInRange


      ! Begin recording the results to disk.

      ! Record the total system DOS, converting the scale to eV. For the LAT
      !   path, TDOS has already been written by computeTDOS_LAT to fort.60/61,
      !   so we skip the TDOS write here.
      if (kPointIntgCode /= 1) then
         write (59+h, fmt="(a14,a14)") &
               & "    ENERGY(eV)", "          TDOS"
         do i = 1, numEnergyPoints
            write (59+h, fmt="(f14.4,f14.6)") &
                  & energyScale(i) * hartree, &
                  & totalSystemDos(i)
         enddo
      endif

      ! Loop over types for the types-based DOS
      if (detailCodePDOS == 0) then

         ! Print the key bits of information for the PDOS output.
         write (69+h,fmt="(a7)") 'STYLE 1'
         write (69+h,fmt="(a10,i6)") 'NUM_UNITS ',numAtomTypes
         write (69+h,fmt="(a11,i9)") 'NUM_POINTS ',numEnergyPoints

         ! Print the energy scale used by all types, converting to eV.
         do i = 1, numEnergyPoints
            write (69+h,fmt="(f16.8)") energyScale(i) * hartree
         enddo

         do i = 1, numAtomTypes

            ! Determine the indices of the pdosComplete to use for this type,
            !   beginning with the first orbital and ending with the last for
            !   this type.
            initIndex = cumulNumDOS(i) + 1
            finIndex  = cumulNumDOS(i+1)


            ! Write the label and orbital information for this type.
            write (69+h,fmt="(a13,i5)") 'SEQUENCE_NUM ',i
            write (69+h,fmt="(a13,a3)") 'ELEMENT_NAME ',atomTypes(i)%elementName
            write (69+h,fmt="(a11,i5)") 'SPECIES_ID ',atomTypes(i)%speciesID
            write (69+h,fmt="(a8,i5)") 'TYPE_ID ',atomTypes(i)%typeID
            write (69+h,fmt="(a10,1x,i2)") 'COL_LABELS',finIndex-initIndex+2
            write (69+h,ADVANCE='NO',fmt="(a6)") 'TOTAL '
            numCols = 1
            do j = 1, lAngMomCount  ! 1=s; 2=p; 3=d; 4=f
               do k = 1, atomTypes(i)%numQN_lValeRadialFns(j)

                  numCols = numCols + 1

                  ! Write the orbital definition in the standard 1s,2s,3s,2p,3p
                  !   notation.  The number is the n quantum number and we
                  !   exclude the orbitals that are in the core and begin
                  !   counting with the orbitals in the valence.  Note that the
                  !   s orbitals ascend like 1s, 2s, 3s while the p orbitals
                  !   ascend like 2p, 3p.  The (+j-1) accounts for this need.
                  write (69+h,ADVANCE='NO',fmt="(i1,a1,1x)") &
                        & atomTypes(i)%numQN_lCoreRadialFns(j)+k+j-1, &
                        & QN_lLetter(j)

                  if (mod(numCols,6) == 0) then
                     write (69+h,*)
                  endif
               enddo
            enddo

            ! Advance the file pointer if the data fit was not exact.
            if (mod(numCols,6) /= 0) then
               write (69+h,*)
            endif

            ! Write the raw data for this unit.  The first value is the sum of
            !   all orbitals for this type.  (TDOS for this type.)  The next
            !   values are the QN_nl resolved PDOS for this type.
            do j = 1, numEnergyPoints
               write (69+h,ADVANCE="NO",fmt="(e13.6,1x)") &
                  & sum(pdosComplete(initIndex:finIndex,j))
               numCols = 1
               do k = initIndex,finIndex
                  numCols = numCols + 1
                  write (69+h,ADVANCE="NO",fmt="(e13.6,1x)") pdosComplete(k,j)
                  if (mod(numCols,6) == 0) then
                    write (69+h,*)
                  endif
               enddo

               ! Advance the file pointer if the data fit was not exact.
               if (mod(numCols,6) /= 0) then
                  write (69+h,*)
               endif
            enddo
         enddo
      ! Loop over atoms for the atom-based DOS
      elseif (detailCodePDOS == 1) then

         ! Print the key bits of information for the PDOS output.
         write (69+h,fmt="(a7)") 'STYLE 1'
         write (69+h,fmt="(a10,i6)") 'NUM_UNITS ', numAtomSites
         write (69+h,fmt="(a11,i9)") 'NUM_POINTS ', numEnergyPoints

         ! Print the energy scale used by all atoms, converting to eV.
         do i = 1, numEnergyPoints
            write (69+h,fmt="(f16.8)") energyScale(i) * hartree
         enddo

         do i = 1, numAtomSites

            ! Obtain the type of the current atom.
            currentType = atomSites(i)%atomTypeAssn

            ! Write the label information for this atom.
            write (69+h,fmt="(a13,i5,1x)") 'SEQUENCE_NUM ',i
            write (69+h,fmt="(a13,a3)") 'ELEMENT_NAME ',&
                  & atomTypes(currentType)%elementName
            write (69+h,fmt="(a11,i5)") 'SPECIES_ID ',&
                  & atomTypes(currentType)%speciesID
            write (69+h,fmt="(a8,i5)") 'TYPE_ID ',&
                  & atomTypes(currentType)%typeID
            write (69+h,fmt="(a10,1x,i2)") 'COL_LABELS',1
            write (69+h,fmt="(a6)") 'TOTAL '

            ! Write the raw data for this unit.
            do j = 1, numEnergyPoints
               write (69+h,fmt="(f16.8)") pdosComplete(i,j)
            enddo
         enddo
      ! Loop over atoms & orbitals (valeDim).
      elseif (detailCodePDOS == 2) then

         ! Print the key bits of information for the PDOS output.
         write (69+h,fmt="(a7)") 'STYLE 1'
         write (69+h,fmt="(a10,i6)") 'NUM_UNITS ',numAtomSites
         write (69+h,fmt="(a11,i9)") 'NUM_POINTS ', numEnergyPoints

         ! Print the energy scale used by all atoms, converting to eV.
         do i = 1, numEnergyPoints
            write (69+h,fmt="(f16.8)") energyScale(i) * hartree
         enddo

         ! Initialize a record of which valence index the loop is on.
         valeDimIndex = 0

         do i = 1, numAtomSites

            ! Obtain the type of the current atom.
            currentType = atomSites(i)%atomTypeAssn

            ! Determine the indices of the pdosComplete to use for this atom,
            !   beginning with the first orbital and ending with the last for
            !   this atom.
            initIndex = cumulNumDOS(i)+1
            finIndex  = cumulNumDOS(i+1)


            ! Write the label information for this atom.
            write (69+h,fmt="(a13,i5,1x)") 'SEQUENCE_NUM ',i
            write (69+h,fmt="(a13,a3)") 'ELEMENT_NAME ',&
                  & atomTypes(currentType)%elementName
            write (69+h,fmt="(a11,i5)") 'SPECIES_ID ',&
                  & atomTypes(currentType)%speciesID
            write (69+h,fmt="(a8,i5)") 'TYPE_ID ',&
                  & atomTypes(currentType)%typeID
            write (69+h,fmt="(a10,1x,i2)") 'COL_LABELS',finIndex-initIndex+2
            write (69+h,ADVANCE='NO',fmt="(a6)") 'TOTAL '
            numCols = 1
            do j = 1, lAngMomCount
               do k = 1, atomTypes(currentType)%numQN_lValeRadialFns(j)

                  numCols = numCols + 1

                  ! Write the orbital definition in the standard 1s,2s,3s,2p,3p
                  !   notation.  The number is the n quantum number and we
                  !   exclude the orbitals that are in the core and begin
                  !   counting with the orbitals in the valence.  Note that the
                  !   s orbitals ascend like 1s, 2s, 3s while the p orbitals
                  !   ascend like 2p, 3p.  The (+j-1) accounts for this need.
                  write (69+h,ADVANCE='NO',fmt="(i1,a1,1x)") &
                        & atomTypes(currentType)%numQN_lCoreRadialFns(j)+k+j-1,&
                        & QN_lLetter(j)

                  if (mod(numCols,6) == 0) then
                     write (69+h,*)
                  endif
               enddo
            enddo

            ! Advance the file pointer if the data fit was not exact.
            if (mod(numCols,6) /= 0) then
               write (69+h,*)
            endif

            ! Write the raw data for this unit.  The first value is the TDOS
            !   for this atom.  The other values are the QN_nl resolved PDOS
            !   for this atom.
            do j = 1, numEnergyPoints
               write (69+h,ADVANCE="NO",fmt="(e13.6,1x)") &
                     & sum(pdosComplete(initIndex:finIndex,j))
               numCols = 1
               do k = initIndex,finIndex
                  numCols = numCols + 1
                  write (69+h,ADVANCE="NO",fmt="(e13.6,1x)") pdosComplete(k,j)
                  if (mod(numCols,6) == 0) then
                    write (69+h,*)
                  endif
               enddo

               ! Advance the file pointer if the data fit was not exact.
               if (mod(numCols,6) /= 0) then
                  write (69+h,*)
               endif
            enddo
         enddo
      ! Loop over atoms and QN_nlm orbitals.
      elseif (detailCodePDOS == 3) then

         ! Print the key bits of information for the PDOS output.
         write (69+h,fmt="(a7)") 'STYLE 1'
         write (69+h,fmt="(a10,i6)") 'NUM_UNITS ',numAtomSites
         write (69+h,fmt="(a11,i9)") 'NUM_POINTS ', numEnergyPoints

         ! Print the energy scale used by all atoms, converting to eV.
         do i = 1, numEnergyPoints
            write (69+h,fmt="(f16.8)") energyScale(i) * hartree
         enddo

         ! Initialize a record of which valence index the loop is on.
         valeDimIndex = 0

         do i = 1, numAtomSites

            ! Obtain the type of the current atom.
            currentType = atomSites(i)%atomTypeAssn

            ! Determine the indices of the pdosComplete to use for this atom,
            !   beginning with the first orbital and ending with the last for
            !   this atom.
            initIndex = cumulNumDOS(i)+1
            finIndex  = cumulNumDOS(i+1)


            ! Write the label information for this atom.
            write (69+h,fmt="(a13,i5,1x)") 'SEQUENCE_NUM ',i
            write (69+h,fmt="(a13,a3)") 'ELEMENT_NAME ',&
                  & atomTypes(currentType)%elementName
            write (69+h,fmt="(a11,i5)") 'SPECIES_ID ',&
                  & atomTypes(currentType)%speciesID
            write (69+h,fmt="(a8,i5)") 'TYPE_ID ',&
                  & atomTypes(currentType)%typeID
            write (69+h,fmt="(a10,1x,i2)") 'COL_LABELS',finIndex-initIndex+2
            write (69+h,ADVANCE='NO',fmt="(a6)") 'TOTAL '
            numCols = 1
            do j = 1, lAngMomCount
               do k = 1, atomTypes(currentType)%numQN_lValeRadialFns(j)
                  do l = 1, (j-1)*2+1

                     numCols = numCols + 1

                     ! Write the QN_nlm orbital definition in the following
                     !   notation:  1s,2s,3s,2px,2py,2pz,3px,3py,3pz,...  The
                     !   first number is the n QN (excluding all orbitals from
                     !   the core).  The first letter is the l QN, and the
                     !   string after that is the m QN in xyz notation.  Note
                     !   that the s orbitals ascend like 1s, 2s, 3s, while the
                     !   p orbitals ascend like 2p, 3p.  The (+j-1) accounts
                     !   for this need.
                     write (formatString,fmt="(a11,i2.2,a4)") "(i1,a1,a1,a", &
                           & len_trim(QN_mLetter(j,l)),",1x)"

                     write (69+h,ADVANCE='NO',fmt=formatString) &
                           & atomTypes(currentType)%&
                           & numQN_lCoreRadialFns(j)+k+j-1, &
                           & QN_lLetter(j),"_",QN_mLetter(j,l)

                     if (mod(numCols,6) == 0) then
                        write (69+h,*)
                     endif
                  enddo
               enddo
            enddo

            ! Advance the file pointer if the data fit was not exact.
            if (mod(numCols,6) /= 0) then
               write (69+h,*)
            endif

            ! Write the raw data for this unit.  The first value is the TDOS
            !   for this atom.  The other values are the QN_nlm resolved PDOS
            !   for this atom.
            do j = 1, numEnergyPoints
               write (69+h,ADVANCE="NO",fmt="(e13.6,1x)") &
                     & sum(pdosComplete(initIndex:finIndex,j))
               numCols = 1
               do k = initIndex,finIndex
                  numCols = numCols + 1
                  write (69+h,ADVANCE="NO",fmt="(e13.6,1x)") pdosComplete(k,j)
                  if (mod(numCols,6) == 0) then
                    write (69+h,*)
                  endif
               enddo

               ! Advance the file pointer if the data fit was not exact.
               if (mod(numCols,6) /= 0) then
                  write (69+h,*)
               endif
            enddo
         enddo
      endif


      ! Compute the weighted average of the energy of a given band across all
      !   kpoints. This energy value is used as the energy for the localization
      !   index.

      ! Initialize the averageEnergy (stored in energyValuesAvg).
      energyValuesAvg(:) = 0.0_double

      do i = 1, numKPoints

         ! Collect weighted sum of the energy for this kpoint over all states.
         energyValuesAvg(:) = energyValuesAvg(:) + &
               & energyEigenValues(:numStates,i,h) * kPointWeight(i)
      enddo

      ! Divide the result by two since the sum of kpoint weights is two and not
      !   one.
      energyValuesAvg(:) = energyValuesAvg(:) / 2.0_double
      localizationIndex(:) = localizationIndex(:) / 2.0_double

      ! Record the results to disk as the localization index for each state,
      !   converting the energy scale to eV.
      do i = 1, numStates
         write (79+h,fmt="(1x,f24.8,2x,f24.8)") energyValuesAvg(i) * hartree,&
               & localizationIndex(i)
      enddo
   enddo ! (h spin)

   ! Deallocate all the unnecessary matrices and arrays. The Gaussian-only
   !   arrays (pdosAccum, waveFnSqrd, overlap) are guarded because the LAT path
   !   handles its own allocations.
   deallocate (cumulNumDOS)
   deallocate (pdosIndex)
   deallocate (numAtomStates)
   deallocate (energyScale)
   deallocate (localizationIndex)
   deallocate (totalSystemDos)
   deallocate (pdosComplete)
   deallocate (electronNumber)

   ! Gaussian-path-only deallocations.
   if (kPointIntgCode /= 1) then
      deallocate (pdosAccum)
#ifndef GAMMA
      deallocate (waveFnSqrd)
      if (inSCF == 0) then
         deallocate (valeValeOL)
      endif
#else
      deallocate (waveFnSqrdGamma)
      if (inSCF == 0) then
         deallocate (valeValeOLGamma)
      endif
#endif
   endif ! (kPointIntgCode /= 1)

   ! LAT-specific deallocation: channel permutation table was built before the
   !   spin loop.
   if (kPointIntgCode == 1) then
      deallocate (channelPermTbl)
   endif

   ! Log the date and time we end.
   call timeStampEnd (19)

end subroutine computeDOS


! Compute the total density of states using the Linear Analytic Tetrahedron
!   (LAT) method (Bloechl, Jepsen, & Andersen, PRB 49, 16223, 1994). This
!   approach decomposes the Brillouin zone into tetrahedra and integrates the
!   DOS analytically within each one, eliminating the need for a broadening
!   parameter.
!
!   The algorithm loops over bands and tetrahedra. For each tetrahedron, the
!   four corner eigenvalues are looked up via the fullKPToIBZKPMap (which maps
!   full-mesh kpoint indices to IBZ kpoint indices), sorted, and the Bloechl
!   analytic formulas are applied to each energy grid point.
!
!   Output is written to the same TDOS file (unit 60/61) in the same format as
!   the Gaussian broadening path. PDOS is not computed in this subroutine
!   (future work).
subroutine computeTDOS_LAT

   ! Import the necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Potential, only: spin
   use O_KPoints, only: numKPoints, numTetrahedra, &
         & tetraVol, tetrahedra, &
         & fullKPToIBZKPMap, kPointWeight
   use O_Constants, only: hartree
   use O_Input, only: numStates, eminDOS, emaxDOS, &
         & deltaDOS
   use O_SecularEquation, only: energyEigenValues

   ! Make sure that no funny variables are used.
   implicit none

   ! Define local variables.
   integer :: h         ! Spin loop index.
   integer :: n         ! Band (state) loop index.
   integer :: t         ! Tetrahedron loop index.
   integer :: iE        ! Energy grid loop index.
   integer :: corner    ! Corner loop index.
   integer :: ibzK      ! IBZ kpoint index for corner.
   integer :: minIdx    ! Index of minimum for sort.
   real (kind=double) :: energy ! Current grid point.
   real (kind=double) :: integratedArea
   real (kind=double) :: tempVal
   real (kind=double) :: kpWtSum ! sum(kPointWeight)
   real (kind=double), allocatable, dimension(:) :: &
         & totalSystemDos
   real (kind=double), dimension(4) :: eps
   real (kind=double), dimension(4) :: cornerDOSWt

   ! Log the date and time we start.
   call timeStampStart (19)

   ! Determine the number of energy grid points.
   numEnergyPoints = int((emaxDOS - eminDOS) / deltaDOS)

   ! Allocate and fill the energy scale.
   if (allocated(energyScale)) deallocate(energyScale)
   allocate (energyScale(numEnergyPoints))
   do iE = 1, numEnergyPoints
      energyScale(iE) = eminDOS + (iE - 1) * deltaDOS
   enddo

   ! Allocate the TDOS accumulation array.
   allocate (totalSystemDos(numEnergyPoints))

   ! Normalization factor: the tetrahedron BZ integration sums tetraVol to 1.0
   !   over the full BZ, but the Gaussian path uses kPointWeight (which sums to
   !   2.0) as its BZ integration weight. The factor of 2 in kPointWeight
   !   accounts for the two electron spin states per band in a spin-unpolarized
   !   calculation (spin=1). Both paths divide by spin separately, so to match
   !   conventions the LAT path must include this same weight sum. For spin=1:
   !   2/1 = 2 spin states per band. For spin=2: 2/2 = 1 spin state per band.
   kpWtSum = sum(kPointWeight(1:numKPoints))

   ! Loop over spin orientations.
   do h = 1, spin

      ! Initialize the TDOS array.
      totalSystemDos(:) = 0.0_double

      ! Record which calculation is being done.
      if (spin == 2) then
         if (h == 1) then
            write (20,*) "Computing spin up LAT TDOS."
         else
            write (20,*) "Computing spin dn LAT TDOS."
         endif
      endif

      ! Loop over bands (states).
      do n = 1, numStates

         ! Loop over tetrahedra.
         do t = 1, numTetrahedra

            ! Look up the four corner kpoint indices from the tetrahedra array
            !   (full mesh indices) and map them to IBZ kpoint indices.
            do corner = 1, 4
               ibzK = fullKPToIBZKPMap(tetrahedra(corner, t))
               eps(corner) = &
                     & energyEigenValues(n, ibzK, h)
            enddo

            ! Sort the four eigenvalues in ascending order using a simple
            !   selection sort (4 elements).
            do corner = 1, 3
               minIdx = corner
               do ibzK = corner + 1, 4
                  if (eps(ibzK) < eps(minIdx)) then
                     minIdx = ibzK
                  endif
               enddo
               if (minIdx /= corner) then
                  tempVal = eps(corner)
                  eps(corner) = eps(minIdx)
                  eps(minIdx) = tempVal
               endif
            enddo

            ! Loop over the energy grid and accumulate the DOS contribution from
            !   this tet.
            do iE = 1, numEnergyPoints
               energy = energyScale(iE)

               ! Skip if outside eigenvalue range.
               if (energy < eps(1) .or. &
                     & energy >= eps(4)) cycle

               ! Compute per-corner DOS weights. The TDOS uses only their sum.
               call bloechlCornerDOSWt( &
                     & energy, eps, cornerDOSWt)

               ! Accumulate into the TDOS. tetraVol is the BZ fraction for one
               !   tet. kpWtSum matches the kPointWeight convention (sums to 2).
               !   1/spin accounts for spin degeneracy. 1/hartree converts
               !   states/Hartree to states/eV.
               totalSystemDos(iE) = &
                     & totalSystemDos(iE) &
                     & + sum(cornerDOSWt) &
                     & * tetraVol * kpWtSum &
                     & / real(spin, double) &
                     & / hartree

            enddo ! iE (energy grid)
         enddo ! t (tetrahedra)

         ! Progress indicator.
         if (mod(n, 50) == 0) then
            write (20, ADVANCE="NO", FMT="(a1)") "."
            call flush(20)
         endif

      enddo ! n (states)

      write (20, *) ""

      ! Compute the integrated area under the TDOS curve using the trapezoidal
      !   rule.
      integratedArea = sum( &
            & totalSystemDos(1:numEnergyPoints - 1) &
            & + totalSystemDos(2:numEnergyPoints)) &
            & * deltaDOS * 0.5_double * hartree
      write (20, *) "LAT TDOS integrated area: ", &
            & integratedArea

      ! Write the TDOS to file (same format as Gaussian).
      write (59+h, fmt="(a14,a14)") &
            & "    ENERGY(eV)", "          TDOS"
      do iE = 1, numEnergyPoints
         write (59+h, fmt="(f14.4,f14.6)") &
               & energyScale(iE) * hartree, &
               & totalSystemDos(iE)
      enddo

   enddo ! h (spin)

   ! Clean up. Deallocate both the local TDOS array and the module-level
   !   energyScale so that computeDOS can re-allocate it independently.
   deallocate (totalSystemDos)
   deallocate (energyScale)

   ! Log the date and time we end.
   call timeStampEnd (19)

end subroutine computeTDOS_LAT


! Compute the four Bloechl corner integration weights for an arbitrary energy E
!   and sorted eigenvalues e1 <= e2 <= e3 <= e4. These weights determine how
!   much of each corner's partial-DOS projection to include at energy E. They
!   are the LAT replacement for Gaussian broadening in the energy-resolved PDOS
!   (DESIGN 1.4, PSEUDOCODE 3a).
!
!   The formulas follow from the vertex-averaging property of linear functions
!   over tetrahedra: the integral of barycentric coordinate lambda_i over any
!   sub-tetrahedron equals V_sub/4 times the sum of lambda_i at the sub-tet's
!   four vertices.
!
!   Four cases arise depending on where E falls relative to the sorted corner
!   eigenvalues:
!     Case 0: E < e1 or E >= e4 (trivial bounds) Case 1: e1 <= E < e2 (small
!     sub-tet near e1) Case 2: e2 <= E < e3 (pentahedral middle) Case 3: e3 <= E
!     < e4 (complement of Case 1)
!
!   Same physics as populate.F90's corner weight computation, but parameterized
!   on an arbitrary energy E rather than E_Fermi = 0.
subroutine bloechlCornerWeights(energy, sortedEps, &
      & cornerWt)

   use O_Kinds

   implicit none

   ! Passed parameters.
   real (kind=double), intent(in) :: energy
   real (kind=double), dimension(4), intent(in) :: &
         & sortedEps
   real (kind=double), dimension(4), intent(out) :: &
         & cornerWt

   ! Local variables.
   real (kind=double) :: e1, e2, e3, e4
   real (kind=double) :: e31, e41, e32, e42
   real (kind=double) :: t2, t3, t4, f
   real (kind=double) :: a, b, c, d
   real (kind=double) :: v_I, v_II, v_III
   real (kind=double) :: s1, s2, s3, f_un
   real (kind=double) :: denom
   real (kind=double), parameter :: tol = 1.0d-12

   e1 = sortedEps(1)
   e2 = sortedEps(2)
   e3 = sortedEps(3)
   e4 = sortedEps(4)

   ! Case 0a: energy below all corners. No spectral weight from this tetrahedron
   !   at this energy.
   if (energy < e1) then
      cornerWt(:) = 0.0_double
      return
   endif

   ! Case 0b: energy above all corners. By vertex averaging over the full
   !   tetrahedron, each corner gets an equal share of 1/4.
   if (energy >= e4) then
      cornerWt(:) = 0.25_double
      return
   endif

   ! Case 1: e1 <= E < e2. The iso-energy surface cuts the three edges from
   !   corner 1, forming a small sub-tetrahedron with apex at corner 1.
   if (energy < e2) then
      denom = (e2-e1) * (e3-e1) * (e4-e1)
      if (abs(denom) < tol) then
         cornerWt(:) = 0.0_double
         return
      endif
      t2 = (energy - e1) / (e2 - e1)
      t3 = (energy - e1) / (e3 - e1)
      t4 = (energy - e1) / (e4 - e1)
      f = t2 * t3 * t4
      cornerWt(2) = f * t2 / 4.0_double
      cornerWt(3) = f * t3 / 4.0_double
      cornerWt(4) = f * t4 / 4.0_double
      cornerWt(1) = f - cornerWt(2) &
            & - cornerWt(3) - cornerWt(4)
      return
   endif

   ! Case 2: e2 <= E < e3 (middle range). Corners 1 and 2 lie below; corners 3
   !   and 4 lie above. The occupied region is a pentahedron decomposed into
   !   three sub-tetrahedra (T_I, T_II, T_III).
   if (energy < e3) then
      e31 = e3 - e1
      e41 = e4 - e1
      e32 = e3 - e2
      e42 = e4 - e2
      if (e31 * e41 < tol .or. &
            & e32 * e42 < tol) then
         cornerWt(:) = 0.0_double
         return
      endif

      ! Intersection parameters: fractional positions where the iso-energy
      !   surface cuts each edge.
      a = (energy - e1) / e31
      b = (energy - e1) / e41
      c = (energy - e2) / e32
      d = (energy - e2) / e42

      ! Sub-tetrahedra volume ratios (as fractions of the full tetrahedron
      !   volume).
      v_I   = a * b
      v_II  = a * d * (1.0_double - b)
      v_III = (1.0_double - a) * c * d

      ! Corner weights from vertex averaging over the three sub-tetrahedra.
      cornerWt(1) = ( &
            & v_I * (3.0_double - a - b) &
            & + v_II * (2.0_double - a - b) &
            & + v_III * (1.0_double - a) &
            & ) / 4.0_double
      cornerWt(2) = ( &
            & v_I &
            & + v_II * (2.0_double - d) &
            & + v_III * (3.0_double - c - d) &
            & ) / 4.0_double
      cornerWt(3) = ( &
            & v_I * a + v_II * a &
            & + v_III * (a + c) &
            & ) / 4.0_double
      cornerWt(4) = ( &
            & v_I * b &
            & + v_II * (b + d) &
            & + v_III * d &
            & ) / 4.0_double
      return
   endif

   ! Case 3: e3 <= E < e4. Only corner 4 lies above the energy. The unoccupied
   !   region is a small sub-tet near corner 4 (complement of Case 1).
   denom = (e4-e1) * (e4-e2) * (e4-e3)
   if (abs(denom) < tol) then
      cornerWt(:) = 0.25_double
      return
   endif
   s1 = (e4 - energy) / (e4 - e1)
   s2 = (e4 - energy) / (e4 - e2)
   s3 = (e4 - energy) / (e4 - e3)
   f_un = s1 * s2 * s3
   cornerWt(1) = 0.25_double &
         & - f_un * s1 / 4.0_double
   cornerWt(2) = 0.25_double &
         & - f_un * s2 / 4.0_double
   cornerWt(3) = 0.25_double &
         & - f_un * s3 / 4.0_double
   cornerWt(4) = (1.0_double - f_un) &
         & - cornerWt(1) - cornerWt(2) &
         & - cornerWt(3)

end subroutine bloechlCornerWeights


! Compute the four per-corner DOS density weights (cornerDOSWt_LAT) for one
!   tetrahedron at a given energy. These are the energy derivatives of the
!   cumulative corner integration weights returned by bloechlCornerWeights:
!
!     cornerDOSWt_LAT(c) = d/dE [cornerIntgWt(c)]
!
!   Units: 1/energy (same as eigenvalue units). Their sum equals the total
!   per-tetrahedron DOS (the dosContrib from the original TDOS code).
!
!   Used by both the TDOS (sum only) and the PDOS (per-corner, to weight
!   Mulliken projections). See DESIGN 1.3 and PSEUDOCODE 2a.
subroutine bloechlCornerDOSWt(energy, sortedEps, &
      & cornerDOSWt)

   use O_Kinds

   implicit none

   ! Passed parameters.
   real (kind=double), intent(in) :: energy
   real (kind=double), dimension(4), intent(in) :: &
         & sortedEps
   real (kind=double), dimension(4), intent(out) :: &
         & cornerDOSWt

   ! Local variables.
   real (kind=double) :: e1, e2, e3, e4
   real (kind=double) :: e31, e41, e32, e42
   real (kind=double) :: t2, t3, t4, f, gTotal
   real (kind=double) :: a, b, cv, dv
   real (kind=double) :: da, db, dc, dd
   real (kind=double) :: v_I, v_II, v_III
   real (kind=double) :: dv_I, dv_II, dv_III
   real (kind=double) :: s1, s2, s3, f_un
   real (kind=double) :: denom
   real (kind=double), parameter :: tol = 1.0d-12

   e1 = sortedEps(1)
   e2 = sortedEps(2)
   e3 = sortedEps(3)
   e4 = sortedEps(4)

   ! Case 0: energy outside eigenvalue range. No spectral density from this
   !   tetrahedron.
   if (energy < e1 .or. energy >= e4) then
      cornerDOSWt(:) = 0.0_double
      return
   endif

   ! ------------------------------------------------- Case 1: e1 <= E < e2
   ! ------------------------------------------------- The cumulative weights
   ! are w(j) = f*t_j/4
   !   for j=2,3,4, and w(1) = f - w(2) - w(3) - w(4), where f = t2*t3*t4 and
   !   t_j = (E-e1) / (e_j - e1). Apply the product rule:
   !     d(f*t_j)/dE = gTotal*t_j + f/(e_j - e1)
   !   where gTotal = df/dE = 3*(E-e1)^2 / denom.
   if (energy < e2) then
      denom = (e2-e1) * (e3-e1) * (e4-e1)
      if (abs(denom) < tol) then
         cornerDOSWt(:) = 0.0_double
         return
      endif
      t2 = (energy - e1) / (e2 - e1)
      t3 = (energy - e1) / (e3 - e1)
      t4 = (energy - e1) / (e4 - e1)
      f = t2 * t3 * t4
      gTotal = 3.0_double &
            & * (energy - e1)**2 / denom

      cornerDOSWt(2) = (gTotal * t2 &
            & + f / (e2 - e1)) / 4.0_double
      cornerDOSWt(3) = (gTotal * t3 &
            & + f / (e3 - e1)) / 4.0_double
      cornerDOSWt(4) = (gTotal * t4 &
            & + f / (e4 - e1)) / 4.0_double
      cornerDOSWt(1) = gTotal &
            & - cornerDOSWt(2) &
            & - cornerDOSWt(3) &
            & - cornerDOSWt(4)
      return
   endif

   ! ------------------------------------------------- Case 2: e2 <= E < e3
   ! (middle range) ------------------------------------------------- The
   ! cumulative weights use three sub-tetrahedra
   !   volumes v_I, v_II, v_III with intersection parameters a, b, c, d. We
   !   compute volume derivatives dv_I, dv_II, dv_III and parameter derivatives
   !   da, db, dc, dd, then apply the product rule to each corner weight
   !   expression.
   if (energy < e3) then
      e31 = e3 - e1
      e41 = e4 - e1
      e32 = e3 - e2
      e42 = e4 - e2
      if (e31 * e41 < tol .or. &
            & e32 * e42 < tol) then
         cornerDOSWt(:) = 0.0_double
         return
      endif

      ! Intersection parameters.
      a  = (energy - e1) / e31
      b  = (energy - e1) / e41
      cv = (energy - e2) / e32
      dv = (energy - e2) / e42

      ! Sub-tetrahedra volume ratios.
      v_I   = a * b
      v_II  = a * dv * (1.0_double - b)
      v_III = (1.0_double - a) * cv * dv

      ! Parameter derivatives (d/dE).
      da = 1.0_double / e31
      db = 1.0_double / e41
      dc = 1.0_double / e32
      dd = 1.0_double / e42

      ! Volume derivatives (d/dE).
      dv_I   = b / e31 + a / e41
      dv_II  = dv * (1.0_double - b) / e31 &
            & + a * (1.0_double - b) / e42 &
            & - a * dv / e41
      dv_III = -cv * dv / e31 &
            & + (1.0_double - a) * dv / e32 &
            & + (1.0_double - a) * cv / e42

      ! Corner 1: w(1) = [v_I*(3-a-b) + v_II*(2-a-b) + v_III*(1-a)] / 4
      cornerDOSWt(1) = ( &
            & dv_I * (3.0_double - a - b) &
            & + v_I * (-da - db) &
            & + dv_II * (2.0_double - a - b) &
            & + v_II * (-da - db) &
            & + dv_III * (1.0_double - a) &
            & + v_III * (-da) &
            & ) / 4.0_double

      ! Corner 2: w(2) = [v_I + v_II*(2-d) + v_III*(3-c-d)] / 4
      cornerDOSWt(2) = ( &
            & dv_I &
            & + dv_II * (2.0_double - dv) &
            & + v_II * (-dd) &
            & + dv_III &
            & * (3.0_double - cv - dv) &
            & + v_III * (-dc - dd) &
            & ) / 4.0_double

      ! Corner 3: w(3) = [v_I*a + v_II*a + v_III*(a+c)] / 4
      cornerDOSWt(3) = ( &
            & dv_I * a + v_I * da &
            & + dv_II * a + v_II * da &
            & + dv_III * (a + cv) &
            & + v_III * (da + dc) &
            & ) / 4.0_double

      ! Corner 4: w(4) = [v_I*b + v_II*(b+d) + v_III*d] / 4
      cornerDOSWt(4) = ( &
            & dv_I * b + v_I * db &
            & + dv_II * (b + dv) &
            & + v_II * (db + dd) &
            & + dv_III * dv &
            & + v_III * dd &
            & ) / 4.0_double
      return
   endif

   ! ------------------------------------------------- Case 3: e3 <= E < e4
   ! ------------------------------------------------- The unoccupied region is
   ! a small sub-tet near
   !   corner 4 with fraction f_un = s1*s2*s3 where s_j = (e4-E)/(e4-e_j).
   !   Derivatives:
   !     ds_j/dE = -1/(e4-e_j) df_un/dE = -gTotal
   !   For j=1,2,3: dw(j)/dE = (gTotal*s_j
   !       + f_un/(e4-e_j)) / 4
   denom = (e4-e1) * (e4-e2) * (e4-e3)
   if (abs(denom) < tol) then
      cornerDOSWt(:) = 0.0_double
      return
   endif
   s1 = (e4 - energy) / (e4 - e1)
   s2 = (e4 - energy) / (e4 - e2)
   s3 = (e4 - energy) / (e4 - e3)
   f_un = s1 * s2 * s3
   gTotal = 3.0_double &
         & * (e4 - energy)**2 / denom

   cornerDOSWt(1) = (gTotal * s1 &
         & + f_un / (e4 - e1)) / 4.0_double
   cornerDOSWt(2) = (gTotal * s2 &
         & + f_un / (e4 - e2)) / 4.0_double
   cornerDOSWt(3) = (gTotal * s3 &
         & + f_un / (e4 - e3)) / 4.0_double
   cornerDOSWt(4) = gTotal &
         & - cornerDOSWt(1) &
         & - cornerDOSWt(2) &
         & - cornerDOSWt(3)

end subroutine bloechlCornerDOSWt


! Build the channel permutation lookup table for LAT PDOS IBZ unfolding
!   (PSEUDOCODE 8.1). For each point group operation R and PDOS channel alpha,
!   channelPermTbl(R, alpha) gives the channel index at the IBZ k-point whose
!   projection should be used for channel alpha at the full-mesh k-point related
!   by R.
!
!   The permutation depends on the PDOS detail mode: Mode 0 (per-type, per-l):
!     identity. Type-level
!       sums are invariant under R because R permutes atoms within each type.
!     Mode 1 (per-atom total): channel = atom index, so permute via invAtomPerm.
!     Mode 2 (per-atom, per-l): remap the atom index via invAtomPerm while
!       preserving the l-shell offset within the atom (same species => same
!       orbital structure).
!     Mode 3 (per-atom, per-lm): not supported with LAT (requires D^l rotation
!       matrices).
!
!   Built once before the spin loop since channel permutation is
!   spin-independent.
subroutine buildChannelPermTable(detailCode, &
      & numOps, totalChannels, cumulDOS, &
      & numSites, invPerm, channelPermTbl)

   use O_Kinds

   implicit none

   ! Passed parameters.
   integer, intent(in) :: detailCode
   integer, intent(in) :: numOps
   integer, intent(in) :: totalChannels
   integer, dimension(:), intent(in) :: cumulDOS
   integer, intent(in) :: numSites
   integer, dimension(:,:), intent(in) :: invPerm
   integer, allocatable, dimension(:,:), &
         & intent(out) :: channelPermTbl

   ! Local variables.
   integer :: opR       ! Point group operation index.
   integer :: alpha     ! Channel loop index.
   integer :: atomA     ! Atom index decoded from alpha.
   integer :: permAtom  ! Permuted atom index.
   integer :: baseOld   ! Cumulative offset for atomA.
   integer :: baseNew   ! Cumulative offset for permAtom.
   integer :: nOrbitals ! Number of l-shells for atomA.
   integer :: orbOff    ! Orbital offset loop index.

   allocate (channelPermTbl(numOps, totalChannels))

   ! Mode 0: per-type, per-l. Type-level sums are invariant under R (R permutes
   !   atoms within each type, so the type sum is unchanged). Channel
   !   permutation is the identity.
   if (detailCode == 0) then
      do opR = 1, numOps
         do alpha = 1, totalChannels
            channelPermTbl(opR, alpha) = alpha
         enddo
      enddo
      return
   endif

   ! Mode 1: per-atom total. Channel index = atom index. Permute directly via
   !   invAtomPerm.
   if (detailCode == 1) then
      do opR = 1, numOps
         do atomA = 1, numSites
            channelPermTbl(opR, atomA) = &
                  & invPerm(opR, atomA)
         enddo
      enddo
      return
   endif

   ! Mode 2: per-atom, per-l. Decode alpha into (atom, l-offset), permute the
   !   atom via invAtomPerm, re-encode using the permuted atom's cumulative
   !   offset. The l-shell offset is unchanged because same species implies
   !   identical orbital structure.
   if (detailCode == 2) then
      do opR = 1, numOps
         do atomA = 1, numSites
            permAtom = invPerm(opR, atomA)
            baseOld = cumulDOS(atomA)
            baseNew = cumulDOS(permAtom)
            nOrbitals = cumulDOS(atomA + 1) &
                  & - cumulDOS(atomA)
            do orbOff = 1, nOrbitals
               channelPermTbl(opR, &
                     & baseOld + orbOff) = &
                     & baseNew + orbOff
            enddo
         enddo
      enddo
      return
   endif

end subroutine buildChannelPermTable


! Pass 1 of the LAT PDOS computation (PSEUDOCODE 8.2). Stream through IBZ
!   k-points, read eigenvectors and overlap from HDF5, compute Mulliken
!   projections for all bands and channels, and store into projArray.
!
!   The Mulliken decomposition is identical to the existing Gaussian path: for
!   each basis function, waveFnSqrd * valeValeOL gives the Mulliken weight of
!   that function in a given band. These are summed by channel according to
!   pdosIndex.
!
!   In addition, electronNumber and localizationIndex are accumulated using the
!   same formulas as the Gaussian path (kPointWeight * step-function occupancy)
!   to provide the normalization diagnostic.
!
!   projArray is allocated here and must be deallocated by the caller after
!   integratePDOS_LAT completes.
subroutine computeProjections_LAT(inSCF, spinIdx, &
      & numKP, numSt, numSites, numAtmSt, pIndex, &
      & vDim, cumDOSTotal, spinCount, projArr, &
      & elecNum, locIdx)

   use O_Kinds
   use O_Potential, only: spin
   use O_KPoints, only: kPointWeight
   use O_Constants, only: pi, hartree
   use O_AtomicSites, only: valeDim
   use O_SecularEquation, only: energyEigenValues, &
         & readDataSCF, readDataPSCF
#ifndef GAMMA
   use O_SecularEquation, only: valeValeOL, valeVale
#else
   use O_SecularEquation, only: valeValeOLGamma, &
         & valeValeGamma
#endif

   implicit none

   ! Passed parameters.
   integer, intent(in) :: inSCF
   integer, intent(in) :: spinIdx
   integer, intent(in) :: numKP
   integer, intent(in) :: numSt
   integer, intent(in) :: numSites
   integer, dimension(:), intent(in) :: numAtmSt
   integer, dimension(:), intent(in) :: pIndex
   integer, intent(in) :: vDim
   integer, intent(in) :: cumDOSTotal
   integer, intent(in) :: spinCount
   real (kind=double), allocatable, &
         & dimension(:,:,:), intent(out) :: projArr
   real (kind=double), dimension(:,:), &
         & intent(inout) :: elecNum
   real (kind=double), dimension(:), &
         & intent(inout) :: locIdx

   ! Local variables.
   integer :: i, j, k, l       ! Loop indices.
   integer :: valeDimIdx        ! Tracks position in the full valence dimension.
   real (kind=double) :: occupNum  ! Step-function occupancy for normalization
         !   diagnostic.
   real (kind=double) :: oneValeRA ! Mulliken projection for one basis function.
#ifndef GAMMA
   complex (kind=double), allocatable, &
         & dimension (:) :: waveFnSqrd
#else
   real (kind=double), allocatable, &
         & dimension (:) :: waveFnSqrdGamma
#endif

   ! Allocate the projection array: (channel, band, IBZ kpoint). This is the
   !   main memory cost of the LAT PDOS computation.
   allocate (projArr(cumDOSTotal, numSt, numKP))
   projArr(:,:,:) = 0.0_double

   ! Allocate work arrays for the Mulliken product.
#ifndef GAMMA
   allocate (waveFnSqrd(vDim))
   allocate (valeValeOL(vDim, vDim))
   if (inSCF == 0) then
      allocate (valeVale(vDim, numSt, 1))
   endif
#else
   allocate (waveFnSqrdGamma(vDim))
   allocate (valeValeOLGamma(vDim, vDim))
   if (inSCF == 0) then
      allocate (valeValeGamma(vDim, numSt, 1))
   endif
#endif

   ! Stream through IBZ k-points, reading eigenvectors and overlap one k-point
   !   at a time.
   do i = 1, numKP

      ! Read eigenvectors + overlap for this IBZ kpoint and spin orientation.
      if (inSCF == 1) then
         call readDataSCF(spinIdx, i, numSt, 1)
      else
         call readDataPSCF(spinIdx, i, numSt, 1)
      endif

      do j = 1, numSt

         ! Step-function occupancy for the electron number diagnostic. States
         !   below the Fermi level (shifted to 0) get 1/spin; above get 0.
         if (energyEigenValues(j,i,spinIdx) &
               & <= 0.0_double) then
            occupNum = 1.0_double &
                  & / real(spinCount, double)
         else
            occupNum = 0.0_double
         endif

         ! Reset the valence dimension tracker.
         valeDimIdx = 0

         ! Loop over all atoms and their valence states to compute Mulliken
         !   projections.
         do k = 1, numSites
            do l = 1, numAtmSt(k)

               valeDimIdx = valeDimIdx + 1

#ifndef GAMMA
               ! Complex case: compute waveFnSqrd = conjg(C_mu) * C_nu for all
               !   nu, then dot with overlap.
               waveFnSqrd(:vDim) = &
                     & conjg(valeVale( &
                     & valeDimIdx, j, 1)) &
                     & * valeVale(:vDim, j, 1)

               oneValeRA = sum( &
                     & real(waveFnSqrd(:vDim), &
                     & double) &
                     & * real(valeValeOL(:vDim, &
                     & valeDimIdx), double) &
                     & + aimag(waveFnSqrd( &
                     & :vDim)) &
                     & * aimag(valeValeOL(:vDim,&
                     & valeDimIdx)))
#else
               ! Gamma case: real eigenvectors.
               waveFnSqrdGamma(:vDim) = &
                     & valeValeGamma( &
                     & valeDimIdx, j, 1) &
                     & * valeValeGamma( &
                     & :vDim, j, 1)

               oneValeRA = sum( &
                     & waveFnSqrdGamma(:vDim) &
                     & * valeValeOLGamma(:vDim,&
                     & valeDimIdx))
#endif

               ! Accumulate into the channel determined by pdosIndex. The 1/spin
               !   factor ensures that the projection sums correctly for
               !   spin-polarized calculations.
               projArr(pIndex(valeDimIdx), j, i) &
                     & = projArr(pIndex( &
                     & valeDimIdx), j, i) &
                     & + oneValeRA &
                     & / real(spinCount, double)

               ! Accumulate electron number for the normalization diagnostic
               !   (same formula as the Gaussian path).
               elecNum(l, k) = elecNum(l, k) &
                     & + oneValeRA &
                     & * kPointWeight(i) * occupNum

               ! Accumulate localization index (same formula as Gaussian path).
               locIdx(j) = locIdx(j) &
                     & + oneValeRA * oneValeRA &
                     & * kPointWeight(i) &
                     & / real(spinCount, double)
            enddo
         enddo
      enddo ! j (states)

      ! Progress indicator.
      if (mod(i, 10) == 0) then
         write (20, ADVANCE="NO", &
               & FMT="(a1)") "|"
      else
         write (20, ADVANCE="NO", &
               & FMT="(a1)") "."
      endif
      if (mod(i, 50) == 0) then
         write (20, *) " ", i
      endif
      call flush(20)

   enddo ! i (kpoints)

   ! Newline after progress if needed.
   if (mod(numKP, 50) /= 0) then
      write (20, *)
      call flush(20)
   endif

   ! Clean up work arrays. projArray is kept alive for the caller to use in Pass
   !   2.
#ifndef GAMMA
   deallocate (waveFnSqrd)
   if (inSCF == 0) then
      deallocate (valeValeOL)
   endif
#else
   deallocate (waveFnSqrdGamma)
   if (inSCF == 0) then
      deallocate (valeValeOLGamma)
   endif
#endif

end subroutine computeProjections_LAT


! Pass 2 of the LAT PDOS computation (PSEUDOCODE 8.3). Loop over bands and
!   tetrahedra. For each tetrahedron, look up the four corner eigenvalues via
!   the IBZ map, sort them with a tracked permutation, compute Bloechl corner
!   weights at each energy grid point, and accumulate weighted projections into
!   pdosComplete using the channel permutation table for IBZ unfolding.
!
!   The corner weight at each energy E distributes spectral density among the
!   four tetrahedron corners. The channel permutation maps the projection stored
!   at the IBZ k-point to the correct channel at the full-mesh k-point.
!
!   pdosComplete must be initialized to zero by the caller before this
!   subroutine is called.
subroutine integratePDOS_LAT(projArr, &
      & channelPermTbl, pdosComp, spinIdx, &
      & numSt, cumDOSTotal, numEPts, eScale)

   use O_Kinds
   use O_Constants, only: hartree
   use O_KPoints, only: numKPoints, numTetrahedra, &
         & tetraVol, tetrahedra, &
         & fullKPToIBZKPMap, fullKPToIBZOpMap, &
         & kPointWeight
   use O_SecularEquation, only: energyEigenValues

   implicit none

   ! Passed parameters.
   real (kind=double), dimension(:,:,:), &
         & intent(in) :: projArr
   integer, dimension(:,:), intent(in) :: &
         & channelPermTbl
   real (kind=double), dimension(:,:), &
         & intent(inout) :: pdosComp
   integer, intent(in) :: spinIdx
   integer, intent(in) :: numSt
   integer, intent(in) :: cumDOSTotal
   integer, intent(in) :: numEPts
   real (kind=double), dimension(:), &
         & intent(in) :: eScale

   ! Local variables.
   integer :: n         ! Band (state) loop index.
   integer :: t         ! Tetrahedron loop index.
   integer :: c         ! Corner loop index (1-4).
   integer :: iE        ! Energy grid loop index.
   integer :: alpha     ! Channel loop index.
   integer :: orig      ! Original corner before sort.
   integer :: minIdx    ! Sort: index of minimum.
   integer :: opR       ! Point group op for corner.
   integer :: kIBZc     ! IBZ kpoint for corner.
   integer :: permAlpha ! Permuted channel index.
   real (kind=double) :: energy   ! Current grid pt.
   real (kind=double) :: tempVal  ! Sort swap temp.
   real (kind=double) :: kpWtSum  ! sum(kPointWeight)
   integer :: tempInt             ! Sort swap temp.

   ! Per-tetrahedron arrays.
   real (kind=double), dimension(4) :: eps
   real (kind=double), dimension(4) :: cornerDOSWt
   integer, dimension(4) :: kFull  ! Full-mesh corners.
   integer, dimension(4) :: kIBZ   ! IBZ corners.
   integer, dimension(4) :: opIdx  ! Op index per corner.
   integer, dimension(4) :: sigma  ! Sort permutation.

   ! Normalization: match the kPointWeight convention used by the Gaussian path.
   !   See the comment in computeTDOS_LAT for the full explanation.
   kpWtSum = sum(kPointWeight(1:numKPoints))

   ! Loop over bands (states).
   do n = 1, numSt

      ! Loop over tetrahedra.
      do t = 1, numTetrahedra

         ! Look up the four corner k-point indices from the tetrahedra array
         !   (full mesh) and map them to IBZ k-point indices and operation
         !   indices.
         do c = 1, 4
            kFull(c) = tetrahedra(c, t)
            kIBZ(c) = &
                  & fullKPToIBZKPMap(kFull(c))
            opIdx(c) = &
                  & fullKPToIBZOpMap(kFull(c))
            eps(c) = energyEigenValues( &
                  & n, kIBZ(c), spinIdx)
            sigma(c) = c
         enddo

         ! Sort the four eigenvalues in ascending order, tracking the
         !   permutation sigma so we can map sorted corners back to their
         !   original IBZ kpoint and operation.
         do c = 1, 3
            minIdx = c
            do iE = c + 1, 4
               if (eps(iE) < eps(minIdx)) then
                  minIdx = iE
               endif
            enddo
            if (minIdx /= c) then
               tempVal = eps(c)
               eps(c) = eps(minIdx)
               eps(minIdx) = tempVal
               tempInt = sigma(c)
               sigma(c) = sigma(minIdx)
               sigma(minIdx) = tempInt
            endif
         enddo

         ! Loop over the energy grid and accumulate weighted projections into
         !   pdosComplete.
         do iE = 1, numEPts
            energy = eScale(iE)

            ! Skip if outside eigenvalue range.
            if (energy < eps(1) .or. &
                  & energy >= eps(4)) cycle

            ! Compute per-corner DOS density weights (not the cumulative corner
            !   weights from bloechlCornerWeights, which are for integrated
            !   properties only).
            call bloechlCornerDOSWt( &
                  & energy, eps, cornerDOSWt)

            ! Accumulate weighted projections. Each sorted corner c maps back to
            !   its original corner sigma(c), whose IBZ kpoint and operation
            !   index determine the projection lookup. The channel permutation
            !   table handles the IBZ unfolding of the channel index.
            do c = 1, 4
               if (abs(cornerDOSWt(c)) &
                     & < 1.0d-30) cycle
               orig = sigma(c)
               opR = opIdx(orig)
               kIBZc = kIBZ(orig)

               do alpha = 1, cumDOSTotal
                  permAlpha = channelPermTbl( &
                        & opR, alpha)
                  pdosComp(alpha, iE) = &
                        & pdosComp(alpha, iE) &
                        & + cornerDOSWt(c) &
                        & * tetraVol * kpWtSum &
                        & / hartree &
                        & * projArr(permAlpha, &
                        & n, kIBZc)
               enddo
            enddo
         enddo ! iE (energy grid)
      enddo ! t (tetrahedra)

      ! Progress indicator for long computations.
      if (mod(n, 50) == 0) then
         write (20, ADVANCE="NO", &
               & FMT="(a1)") "."
         call flush(20)
      endif

   enddo ! n (states)

   write (20, *) ""

end subroutine integratePDOS_LAT


end module O_DOS
