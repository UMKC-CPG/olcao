module O_OptcTransitions

   ! Import necessary modules
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define module data.
   real (kind=double) :: orbitalDiff ! Energy difference between highest
         !   occupied state and lowest unoccupied state. (Normal energy onset.)
   real (kind=double) :: energyCutoff ! The energy cut off value determined
         !   from the input data and the type of calculation.  This represents
         !   the highest energy level above 0 included in the calc.
   real (kind=double) :: energyMin ! The minimum energy in the energy window
         !   to be computed for.
   real (kind=double) :: maxTransEnergy ! The maximum transition energy
         !   determined from the input data and the type of calculation.  This
         !   represents the largest transition energy included in the calc.
   real (kind=double), allocatable, dimension (:) :: energyScale ! The fine
         !   scale used for the energy axis.

   integer :: maxPairs ! The largest value in the transCounter array defined
         !   below.
   integer, allocatable, dimension (:,:) :: transCounter ! A count of the
         !   number of transitions present for each kpoint of both spin dirs.

   integer, allocatable, dimension (:,:) :: firstOccupiedState ! Array with the
         !   index number of the first occupied state for each kpoint and spin.
   integer, allocatable, dimension (:,:) :: lastOccupiedState ! Array with the
         !   index number of the last occupied state for each kpoint and spin.
   integer, allocatable, dimension (:,:) :: firstUnoccupiedState ! Array with
         !   the index number of the first unoccupied state for each
         !   kpoint and spin.
   integer, allocatable, dimension (:,:) :: lastUnoccupiedState ! Array with the
         !   index number of the last unoccupied state for each kpoint and spin.

   real (kind=double), allocatable, dimension (:) :: indirectGap ! This is the
         !   smallest energy difference between the highest occupied state and
         !   the lowest unoccupied state for all kpoints.
   real (kind=double), allocatable, dimension (:) :: directGap ! This is the
         !   smallest energy difference between the highest occupied state and
         !   the lowest unoccupied state for any single kpoint.

   real (kind=double), allocatable, dimension (:,:,:,:) :: energyMom
   real (kind=double), allocatable, dimension (:,:,:)   :: energyDiff

#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:,:) :: valeValeMom
   complex (kind=double), allocatable, dimension (:,:,:)   :: coreValeOL
#else
   real (kind=double), allocatable, dimension (:,:,:) :: valeValeMomGamma
   real (kind=double), allocatable, dimension (:,:)   :: coreValeOLGamma
#endif

   real (kind=double), allocatable, dimension (:,:) :: sigmaEAccumulator

contains

subroutine getEnergyStatistics

   ! Import necessary modules.
   use O_Kinds
   use O_Constants
   use O_KPoints
   use O_Input
   use O_CommandLine
   use O_SecularEquation
   use O_Potential

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables.
   integer :: h,i,j,k ! Loop index variables.
   real (kind=double) :: maxOccupiedEnergy
   real (kind=double) :: minUnoccupiedEnergy
   real (kind=double) :: currentEnergyDiff
   real (kind=double) :: currentGap
   integer :: firstInit
   integer :: lastInit
   integer :: firstFin
   integer :: lastFin


   ! Pull variables out of imported modules.
   if (stateSet == 0) then      ! Standard optical properties calculation.
      ! The energy onset for standard optical properties calculations is the
      !   band gap width in eV.  The energy scale for evaluating the
      !   accumulated (and broadened) transitions will begin as close to 0 eV
      !   as conveniently possible.
      energyCutoff   = cutoffEnOPTC
      maxTransEnergy = maxTransEnOPTC
   elseif (stateSet == 1) then  ! PACS type calculation
      ! PACS calculations have the interesting feature that the energy onset is
      !   at some very high energy that is dependent on the particular target
      !   atom being excited (elemental dependency).  This is because the
      !   transitions we are considering are from a deep core state to the
      !   conduction band.  We want that energy onset to be defined by the
      !   difference in total energy between the ground state and the excited
      !   state instead of the orbital energy difference.  (WHY?)

      !   To make that happen we must first subtract away the original energy
      !   onset which is defined by the orbital energy difference between the
      !   initial state and the lowest energy state that the electron is
      !   excited into (bottom of the conduction band).
      energyCutoff   = bigThresh

      ! Compute the max transition energy as the total energy difference
      !   between ground and excited states.  Since we want the first value to
      !   be a nice round number (multiple of 5) we subtract out so-called
      !   "slack" which is the remainder of an integer division by 5.  This is
      !   the lowest energy in the output data set (no transitions at this
      !   energy though).  Then, we add the energy window we want to compute
      !   for to get the maximum transition energy.
      energyMin      = totalEnergyDiffPACS - mod(totalEnergyDiffPACS,&
            & onsetEnergySlackPACS)
      maxTransEnergy = energyMin + energyWindowPACS
   else                         ! Sigma(E) type calculation
      energyCutoff   = cutoffEnSIGE
      maxTransEnergy = maxTransEnSIGE
   endif


   ! Allocate arrays
   allocate (transCounter         (numKPoints,spin))
   allocate (firstOccupiedState   (numKPoints,spin))
   allocate (lastOccupiedState    (numKPoints,spin))
   allocate (firstUnoccupiedState (numKPoints,spin))
   allocate (lastUnoccupiedState  (numKPoints,spin))
   allocate (directGap   (spin))
   allocate (indirectGap (spin))

   ! The purpose of this subroutine is to gather important statistics and
   !   indices for use later on.  The important values that will be determined
   !   are:  1) The first and last occupied state and the first and last 
   !   unoccupied state for each kpoint and each spin that each different type
   !   of calculation cares about (e.g. PACS transitions are from one core
   !   state to numerous CB states; OPTC transitions are from many VB states to
   !   numerous CB states; and SIGE considers a finite range near the fermi
   !   level); 2) The direct band gap for each spin; 3) The minimum indirect
   !   band gap for each spin; 4) The number of transitions for each kpoint and
   !   spin so that we can allocate memory easily; 5) whatever ...

   do h = 1, spin

      ! Initialize the values used to find the indirect gap and the direct gap.
      maxOccupiedEnergy = -bigThresh
      minUnoccupiedEnergy = bigThresh
      directGap(h) = bigThresh

      ! Pacs calculations need to have the resultant spectra shifted according
      !   to the difference in orbital energies from the ground and excited
      !   states.  To find the amount of the shift we initialize the search
      !   number.
      if (stateSet == 1) then  ! Doing PACS calculation.
         orbitalDiff = bigThresh
      endif

      do i = 1, numKPoints
         do j = 1, numStates

            ! Find the last occupied state for this KPoint. (This works because
            !   the energy values have already been shifted such that the
            !   highest occupied state of all kpoints is at 0.0 eV.)  Also note
            !   that this will deal with degenerate highest occupied states.
            if (energyEigenValues(j,i,h) > 0.0_double) then

               ! Store the index values for the last occupied and first
               !   unoccupied states for this kpoint.
               lastOccupiedState(i,h) = j-1
               firstUnoccupiedState(i,h) = j

               ! Compute the size of the direct gap for this kpoint.
               currentGap = abs(energyEigenValues(j,i,h) - &
                     & energyEigenValues(j-1,i,h))

               ! Compare the smallest overall gap to the gap for this kpoint to
               !   get the new smallest overall direct gap so far.
               if (currentGap < directGap(h)) then
                  directGap(h) = currentGap
               endif

               ! Compute the energy values of these states so that we can later
               !   determine the indirect band gap.
               if (energyEigenValues(j-1,i,h) > maxOccupiedEnergy) then
                  maxOccupiedEnergy = energyEigenValues(j-1,i,h)
               endif
               if (energyEigenValues(j,i,h) < minUnoccupiedEnergy) then
                  minUnoccupiedEnergy = energyEigenValues(j,i,h)
               endif

               ! If we are doing a sige calculation then we are only concerned
               !   with states that are a few eV from the Fermi energy.  In this
               !   case we will search for the beginning and ending states.
               if (stateSet == 2) then

                  ! Loop higher than the current j-loop state to find the
                  !   lowest unoccupied state that is greater than the maximum
                  !   transition energy from the Fermi level.
                  do k = j, numStates
                     if (energyEigenValues(k,i,h) > maxTransEnergy) then
                         lastUnoccupiedState(i,h) = k
                         exit
                     endif
                     if (k == numStates) then
                        lastUnoccupiedState(i,h) = numStates
                     endif
                  enddo


                  ! Loop lower than the current j-loop state to find the
                  !   highest occupied state that is less than the maximum
                  !   transition energy from the Fermi level.
                  do k = 1, j-1
                     if (abs(energyEigenValues(j-k,i,h)) > maxTransEnergy) then
                        firstOccupiedState(i,h) = j-k
                        exit
                     endif
                     if (k == j-1) then
                        firstOccupiedState(i,h) = 1
                     endif
                  enddo
               endif

               ! No need to search through any higher states in the outer j
               !   loop.  All relevant information has been obtained.
               exit
            endif
         enddo

         ! Initialize the counter for the number of transitions for this kpoint.
         transCounter(i,h) = 0

         ! For normal optical properties calculations the last unoccupied state
         !   we care about is the last (highest) state in the calculation, and
         !   the first occupied state is always the first (lowest) state in the
         !   calculation.  For PACS, these values depend on which core state
         !   has been excited.  For Sigma(E) these values depend

         if (stateSet == 0) then ! Doing normal optical properties calculation.
            firstOccupiedState(i,h) = 1
            ! lastOccupiedState determined above.
            ! firstUnoccupiedState determined above.
            lastUnoccupiedState(i,h) = numStates
         elseif (stateSet == 1) then ! Doing a PACS calculation.
            firstOccupiedState(i,h) = firstInitStatePACS ! From O_Input
            lastOccupiedState(i,h)  = lastInitStatePACS  ! From O_Input
            ! firstUnoccupiedState determined above.
            lastUnoccupiedState(i,h) = numStates
         elseif (stateSet == 2) then ! Doing a Sigma(E) calculation.
            ! firstOccupiedState determined above.
            ! lastOccupiedState determined above.
            ! firstUnoccupiedState determined above.
            ! lastUnoccupiedState determined above.
         else
            ! Error, no other options.
            stop "Check optical properties command line parameter:  stateSet"
         endif

         ! Store the state variables for temporary use as loop indices.
         firstInit = firstOccupiedState(i,h)
         lastInit  = lastOccupiedState(i,h)
         firstFin  = firstUnoccupiedState(i,h)
         lastFin   = lastUnoccupiedState(i,h)

         if (stateSet == 1) then ! Doing PACS calculation.
            ! Determine the orbital energy difference for this kpoint and
            !   compare it to the smallest difference yet obtained.
            orbitalDiff = min(orbitalDiff,&
                  & energyEigenValues(firstFin,i,h) - &
                  & energyEigenValues(lastInit,i,h))
         endif

         ! Loop over all the possible transitions for this kpoint to determine
         !   the number of accepted transitions for this kpoint.  We will also
         !   refine the value for the lastUnoccupied state to consider.
         do j = firstInit, lastInit
            do k = firstFin, lastFin

               ! If the energy of the final state is higher than the requested
               !   cut-off then we adjust the record for the last unoccupied
               !   state and go to the next initial state because all the
               !   remaining final states will be greater.
               if (energyEigenValues(k,i,h) > energyCutoff) then
                  lastUnoccupiedState(i,h) = k-1
                  exit
               endif

               ! Compute the energy of the transition from the current states.
               currentEnergyDiff = energyEigenValues(k,i,h) - &
                     & energyEigenValues(j,i,h)

               ! Check if the energy difference is less than the maximum
               !   transition energy that the input file requested computation
               !   for.  If it fails, then we adjust the record for the last
               !   unoccupied state and go to the next initial state because
               !   all the remaining final states for this energy will be
               !   greater.
               if (currentEnergyDiff > maxTransEnergy) then
                  lastUnoccupiedState(i,h) = k-1
                  exit
               endif

               ! At this point the transition is valid and one we would want to
               !   compute.  Unfortunately there isn't a good way to save the
               !   energyDiff computation that we did above for later use, and
               !   it will have to be done again.  (BOO)  (Unless we do some
               !   static memory allocation. (BOO))
               transCounter(i,h) = transCounter(i,h) + 1
            enddo
         enddo
      enddo

      ! Determine the band gap.
      indirectGap(h) = minUnoccupiedEnergy - maxOccupiedEnergy

      if (spin == 1) then
         write (20,*) "Indirect Band Gap(eV) = ",indirectGap(h)*hartree
         write (20,*) "Direct Band Gap(eV)   = ",directGap(h)*hartree
         call flush (20)
      elseif (h == 1) then
         write (20,*) "(Up) Indirect Band Gap(eV) = ",indirectGap(h)*hartree
         write (20,*) "(Up) Direct Band Gap(eV)   = ",directGap(h)*hartree
         call flush (20)
      else ! spin == 2 and h == 2
         write (20,*) "(Dn) Indirect Band Gap(eV) = ",indirectGap(h)*hartree
         write (20,*) "(Dn) Direct Band Gap(eV)   = ",directGap(h)*hartree
         call flush (20)
      endif

      ! Determine the maximum number of transition pairs of all the kpoints and
      !   for both spins.
      maxPairs = max(maxPairs,maxval(transCounter(:,h)))
   enddo ! (h spin)

   ! Now that the number of transitions for each kpoint are known we can 
   !   allocate space to hold information based on the number transitions.
   if (stateSet /= 2) then  ! Not doing a Sigma(E) calculation.
      allocate (energyDiff (maxPairs,numKPoints,spin))
      allocate (energyMom  (dim3,maxPairs,numKPoints,spin))

      ! Initialize these arrays.
      energyDiff(:maxPairs,:numKPoints,:) = 0.0_double
      energyMom(:dim3,:maxPairs,:numKPoints,:) = 0.0_double
   endif

end subroutine getEnergyStatistics



subroutine computeTransitions


   ! Import the necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Input
   use O_KPoints
   use O_CommandLine
   use O_MatrixSubs
   use O_SecularEquation
   use O_AtomicSites
   use O_IntegralsPSCF
   use O_Potential

   ! Import the necessary HDF data.
   use HDF5
   use O_PSCFBandHDF5

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define local variables.
   integer :: h,i,j ! Loop index variables
   real (kind=double) :: energyShift
   real (kind=double), allocatable, dimension (:,:) :: tempRealValeVale
#ifndef GAMMA
   real (kind=double), allocatable, dimension (:,:) :: tempImagValeVale
#endif


   ! Log the date and time we start.
   call timeStampStart (23)

   ! Allocate the matrix to hold the wave functions and initialize it.
#ifndef GAMMA
   allocate (valeVale (valeDim,numStates,1,1))
   valeVale(:,:,1,1) = cmplx(0.0_double,0.0_double,double)
#else
   allocate (valeValeGamma(valeDim,numStates,1))
   valeValeGamma(:,:,1) = 0.0_double
#endif


   do h = 1, spin

      ! Begin a loop over the number of kpoints
      do i = 1, numKPoints


         ! Allocate temporary reading matrices.
#ifndef GAMMA
         allocate (tempRealValeVale(valeDim,numStates))
         allocate (tempImagValeVale(valeDim,numStates))
#else
         allocate (tempRealValeVale(valeDim,numStates))
#endif

         if (stateSet /= 1) then  ! Not doing a PACS calculation.

            ! Read the datasets for this kpoint.
#ifndef GAMMA
            call readMatrix(eigenVectorsBand_did(:,i,h),valeVale(:,:,1,1),&
                  & tempRealValeVale(:,:),tempImagValeVale(:,:),&
                  & valeStatesBand,valeDim,numStates)
#else
            call readMatrixGamma(eigenVectorsBand_did(1,i,h),&
                  & valeValeGamma(:,:,1),valeStatesBand,valeDim,numStates)
#endif
         else

#ifndef GAMMA
            ! Read the data for the ground state for this kpoint.
            call readPartialWaveFns(eigenVectorsBand_did(:,i,h),&
                  & valeVale(:,:,1,1),tempRealValeVale(:,:),&
                  & tempImagValeVale(:,:),valeStatesBand,&
                  & firstOccupiedState(i,h),lastOccupiedState(i,h),&
                  & valeDim,numStates)

            ! Read the data for the excited state for this kpoint.
            call readPartialWaveFns(eigenVectorsBand2_did(:,i,h),&
                  & valeVale(:,:,1,1),tempRealValeVale(:,:),&
                  & tempImagValeVale(:,:),valeStatesBand,&
                  & firstUnoccupiedState(i,h),lastUnoccupiedState(i,h),&
                  & valeDim,numStates)
#else
            ! Read the data for the ground state for this kpoint.
            call readPartialWaveFnsGamma(eigenVectorsBand_did(1,i,h),&
                  & valeValeGamma(:,:,1),tempRealValeVale(:,:),valeStatesBand,&
                  & firstOccupiedState(i,h),lastOccupiedState(i,h),&
                  & valeDim,numStates)

            ! Read the data for the excited state for this kpoint.
            call readPartialWaveFnsGamma(eigenVectorsBand2_did(1,i,h),&
                  & valeValeGamma(:,:,1),tempRealValeVale(:,:),valeStatesBand,&
                  & firstUnoccupiedState(i,h),lastUnoccupiedState(i,h),&
                  & valeDim,numStates)
#endif
         endif


         ! Read the orthogonalization coefficients after allocating space to
         !   hold them (and the temp reading matrix).
#ifndef GAMMA
         deallocate (tempRealValeVale)
         deallocate (tempImagValeVale)
         allocate   (tempRealValeVale (coreDim,valeDim))
         allocate   (tempImagValeVale (coreDim,valeDim))
         allocate   (coreValeOL (coreDim,valeDim,1))
         if (coreDim /= 0) then
            call readMatrix(coreValeBand_did(:,i),coreValeOL(:,:,1),&
                  & tempRealValeVale(:,:),tempImagValeVale(:,:),&
                  & coreValeBand,coreDim,valeDim)
         endif
         deallocate (tempRealValeVale)
         deallocate (tempImagValeVale)
#else
         deallocate (tempRealValeVale)
         allocate   (coreValeOLGamma  (coreDim,valeDim))
         if (coreDim /= 0) then
            call readMatrixGamma(coreValeBand_did(1,i),coreValeOLGamma(:,:),&
                  & coreValeBand,coreDim,valeDim)
         endif
#endif


         ! Perform the computations in serial or all together.
         if (serialXYZ == 0) then
#ifndef GAMMA

            allocate   (valeValeMom (valeDim,valeDim,1,3))
            ! Get the integral results for the x, y, z momentum matrices.
            ! Runcode:  3 = XMom, 4 = YMom, 5 = ZMom; (j+2)
            do j = 1, 3
               call getIntgResults (valeValeMom(:,:,:,j),coreValeOL,&
                     & i,j+2,valeValeBand_did(i),valeValeBand,0,1)
            enddo
            ! Deallocate matrices that are no longer needed in this iteration to
            !   make space for those that are needed.
            deallocate (coreValeOL)
#else
            allocate   (valeValeMomGamma (valeDim,valeDim,3))
            ! Get the integral results for the x, y, z momentum matrices.
            ! Runcode:  3 = XMom, 4 = YMom, 5 = ZMom; (j+2)
            do j = 1, 3
               call getIntgResults (valeValeMomGamma(:,:,j),coreValeOLGamma,&
                     & j+2,valeValeBand_did(i),valeValeBand,0,1)
            enddo
            ! Deallocate matrices that are no longer needed in this iteration to
            !   make space for those that are needed.
            deallocate (coreValeOLGamma)
#endif

            if (stateSet /= 2) then  ! Not doing a Sigma(E) calculation.
               call computePairs (i,0,h)
            else
               call computeSigmaE (i,0,h)
            endif

         else
            do j = 1, 3
#ifndef GAMMA
               allocate   (valeValeMom (valeDim,valeDim,1,1))
               ! Runcode:  3 = XMom, 4 = YMom, 5 = ZMom; (j+2)
               call getIntgResults (valeValeMom(:,:,:,1),coreValeOL,&
                     & i,j+2,valeValeBand_did(i),valeValeBand,0,1)
#else
               allocate   (valeValeMomGamma (valeDim,valeDim,1))
               ! Runcode:  3 = XMom, 4 = YMom, 5 = ZMom; (j+2)
               call getIntgResults (valeValeMomGamma(:,:,1),coreValeOLGamma,&
                     & j+2,valeValeBand_did(i),valeValeBand,0,1)
#endif

               if (stateSet /= 2) then  ! Not doing a Sigma(E) calc.
                  call computePairs (i,j,h)
               else
                  call computeSigmaE (i,j,h)
               endif
            enddo

#ifndef GAMMA
            deallocate (coreValeOL)
#else
            deallocate (coreValeOLGamma)
#endif
         endif


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

      enddo

      ! Add a final return if one was not already made.
      if (mod(numKPoints,10) /= 0) then
         write (20,*)
      endif


      ! In the case of XANES calculations we need to shift the calculated
      !   energy values.  The amount of the shift is equal to the difference
      !   between (the total energy difference between the ground and excited
      !   states) and (the calculated LUMO,core difference).
      if (stateSet == 1) then  ! Doing PACS calculation.
         energyShift = totalEnergyDiffPACS - orbitalDiff
         energyDiff(:,:,h) = energyDiff(:,:,h) + energyShift
      endif
   enddo ! (h spin)

   ! Deallocate unnecessary matrices and arrays
#ifndef GAMMA
   deallocate (valeVale)
#else
   deallocate (valeValeGamma)
#endif
   deallocate (firstOccupiedState)
   deallocate (lastOccupiedState)
   deallocate (firstUnoccupiedState)
   deallocate (lastUnoccupiedState)
   deallocate (directGap)
   deallocate (indirectGap)

   ! Log the date and time we end.
   call timeStampEnd (23)

end subroutine computeTransitions



subroutine computePairs (currentKPoint,xyzComponents,spinDirection)

   ! Import the necessary modules.
   use O_Kinds
   use O_KPoints
   use O_SortSubs
   use O_Input
   use O_CommandLine
   use O_AtomicSites
   use O_SecularEquation

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer :: currentKPoint
   integer :: xyzComponents ! 0=all, 1=x, 2=y, 3=z
   integer :: spinDirection

   ! Define local variables.
   integer :: i,j,k ! Loop index variables
   integer :: initComponent
   integer :: finComponent
   integer :: transPairCount
   integer :: firstInit
   integer :: lastInit
   integer :: firstFin
   integer :: lastFin
   integer :: finalStateIndex
   integer, allocatable, dimension (:) :: sortOrder
   integer, allocatable, dimension (:) :: segmentBorders
   real    (kind=double) :: currentEnergyDiff
   real    (kind=double), allocatable, dimension (:)     :: energyDiffTemp
   real    (kind=double), allocatable, dimension (:,:)   :: energyMomTemp
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: conjWaveMomSum
   complex (kind=double),              dimension (dim3)  :: valeValeXMom
#else
   real    (kind=double), allocatable, dimension (:,:,:) :: conjWaveMomSumGamma
   real    (kind=double),              dimension (dim3)  :: valeValeXMomGamma
#endif

   ! Initialize a counter for the current number of transition pairs
   transPairCount = 0

   ! Make shorthand for the state indices.
   firstInit = firstOccupiedState(currentKPoint,spinDirection)
   lastInit  = lastOccupiedState(currentKPoint,spinDirection)
   firstFin  = firstUnoccupiedState(currentKPoint,spinDirection)
   lastFin   = lastUnoccupiedState(currentKPoint,spinDirection)

   ! Determine the range of components (xyz) that should be considered.
   if (xyzComponents == 0) then
      initComponent = 1
      finComponent = 3
   else
      initComponent = 1
      finComponent = 1
   endif


#ifndef GAMMA

   ! Allocate space to hold the sum(conjg(valeVale(:,j)) * valeVale_Mom(:,k,1))
   !   for each of the possible final states.  This is done since the values
   !   are independent of the initial states.  The 1 for the valeVale is
   !   for the 1 kpoint.  The finComponent is 3 for all three at once, and
   !   1 for when X, Y, Z are done separately.
   allocate (conjWaveMomSum (valeDim,lastFin-firstFin+1,finComponent))

   ! Compute the sum over the final states.
   do i = initComponent, finComponent
      finalStateIndex = 0
      do j = firstFin, lastFin
         ! Define the final index for conjWaveMomSum
         finalStateIndex = finalStateIndex + 1
         do k = 1, valeDim
            conjWaveMomSum(k,finalStateIndex,i) = &
                  & sum(conjg(valeVale(:,j,1,1)) * valeValeMom(:,k,1,i))
         enddo
      enddo
   enddo

   ! Deallocate the momentum matrix elements.
   deallocate (valeValeMom)

#else
   ! Documentation similar to the above non-gamma case.
   allocate (conjWaveMomSumGamma (valeDim,lastFin-firstFin+1,finComponent))


   ! Compute the sum over the final states.
   do i = initComponent, finComponent

      ! Make the upper triangle correct for Hermiticity.  Recall that for
      !   the Gamma K Point all the matrices are real (except the momentum
      !   matrix which was multiplied by a -i and is hence imaginary).
      !   Since it must be Hermitian we need to apply that now.
      do j = 1, valeDim
         valeValeMomGamma(1:j,j,i) = -valeValeMomGamma(1:j,j,i)
      enddo

      finalStateIndex = 0
      do j = firstFin, lastFin

         ! Increment the finalStateIndex for conjWaveMomSum
         finalStateIndex = finalStateIndex + 1

         do k = 1, valeDim
            conjWaveMomSumGamma(k,finalStateIndex,i) = &
                  & sum(valeValeGamma(:,j,1) * valeValeMomGamma(:,k,i))
         enddo
      enddo
   enddo

   ! Deallocate the momentum matrix elements.
   deallocate (valeValeMomGamma)
#endif


   ! Allocate space for the energy difference.
   allocate (energyDiffTemp (maxPairs))
   allocate (energyMomTemp  (finComponent,maxPairs))

   ! Initialize the temporary energy transition array.
   energyDiffTemp(:) = 0.0_double

   ! Allocate space to hold the indices for each segment of the energyDiff
   !   array.
   allocate (segmentBorders (lastInit-firstInit+2))

   ! Initialize the first index since it will always be 0.
   segmentBorders(1) = 0

write (20,*) "firstInit, lastInit=",firstInit,lastInit
write (20,*) "firstFin, lastFin=",firstFin,lastFin
call flush (20)
   ! Begin the double loop to determine the transition energies.
   do i = firstInit, lastInit
      finalStateIndex = 0
      do j = firstFin, lastFin

         ! If the energy of the final state is higher than the requested
         !   cut-off we go to the next initial state.
         if (energyEigenValues(j,currentKPoint,spinDirection) > &
               & energyCutoff) exit

         ! Compute the energy of the transition from the current states.
         currentEnergyDiff = energyEigenValues(j,currentKPoint,spinDirection)-&
               & energyEigenValues(i,currentKPoint,spinDirection)

         ! Check if the energy difference is less than the maximum
         !   transition energy that the input file requested computation
         !   for.  If it fails, then we go to the next initial state because
         !   all the remaining final states for this energy will be greater.
         if (currentEnergyDiff > maxTransEnergy) exit

         ! Increment the number of transition pairs counted so far.
         transPairCount = transPairCount + 1

         ! Store the transition energy for the current pair.
         energyDiffTemp(transPairCount) = currentEnergyDiff

         ! Increment the final state index for the conjWaveMomSum
         finalStateIndex = finalStateIndex + 1

#ifndef GAMMA

         ! Loop to obtain the wave function times the momentum integral.
         do k = initComponent,finComponent
             valeValeXMom(k) = sum(valeVale(:,i,1,1) * &
                   & conjWaveMomSum(:,finalStateIndex,k))

            ! Compute the imaginary component of the square of the valeValeXMom
            !   to obtain the momentum energy????  Note that the reason it
            !   is imaginary is because of the negative sign included in the
            !   getIntgResults subroutine for the momentum matrix.
            energyMomTemp(k,transPairCount) = &
                  &  real(valeValeXMom(k),double)**2 + aimag(valeValeXMom(k))**2
         enddo
#else

         ! Loop to get the wave function times the momentum matrix element.
         do k = initComponent,finComponent
            valeValeXMomGamma(k) = sum(valeValeGamma(:,i,1) * &
                  & conjWaveMomSumGamma(:,finalStateIndex,k))

            ! Compute the real component of the square of the valeValeXMom to
            !   obtain the momentum energy????
            energyMomTemp(k,transPairCount) = valeValeXMomGamma(k)**2
         enddo
#endif
      enddo

      ! Save the index for the end border of this segment.
      segmentBorders(i - firstInit + 2) = transPairCount
   enddo


   ! Deallocate unnecessary matrix
#ifndef GAMMA
   deallocate (conjWaveMomSum)
#else
   deallocate (conjWaveMomSumGamma)
#endif

   ! Determine if there was only one segment.  In this case we don't have to
   !   sort anything.

   ! Sort energyDiffTemp into energyDiff, and obtain the indices for the
   !   correct sorted order of energyDiff so that we can copy the energy
   !   momentum directly.

   allocate (sortOrder (transPairCount))

   call mergeSort (energyDiffTemp,energyDiff(:,currentKPoint,spinDirection),&
         & sortOrder,segmentBorders,transPairCount)

   ! Copy energyMomTemp to the real energyMom data structure using the sorting
   !   order determined in the mergeSort subroutine.
   if (xyzComponents == 0) then
      do i = 1, transPairCount
         energyMom(:,i,currentKPoint,spinDirection) = &
               & energyMomTemp(:,sortOrder(i))
      enddo
   else
      do i = 1, transPairCount
         energyMom(xyzComponents,i,currentKPoint,spinDirection) = &
               & energyMomTemp(1,sortOrder(i))
      enddo
   endif

   ! Deallocate unnecessary arrays and matrices
   deallocate (energyDiffTemp)
   deallocate (energyMomTemp)
   deallocate (segmentBorders)
   deallocate (sortOrder)

end subroutine computePairs



subroutine computeSigmaE (currentKPoint,xyzComponents,spinDirection)

   ! Import the necessary data modules.
   use O_Kinds
   use O_Constants
   use O_KPoints
   use O_SortSubs
   use O_Input
   use O_Lattice
   use O_AtomicSites
   use O_SecularEquation

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer :: currentKPoint
   integer :: xyzComponents
   integer :: spinDirection

   ! Define local variables.
   integer :: i,j,k ! Loop index variables
   integer :: initComponent
   integer :: finComponent
   integer :: firstInit
   integer :: lastInit
   integer :: firstFin
   integer :: lastFin
   integer :: numEnergyStates
   integer :: finalStateIndex
   integer :: initialStateIndex
   integer :: numEnergyPoints
   real (kind=double) :: expAlpha
   real (kind=double) :: broadenEnergyDiff
   real (kind=double) :: kPointFactor
   real (kind=double) :: sigmaSqrtPi
   real (kind=double) :: conversionFactor
   real    (kind=double), allocatable, dimension (:,:)   :: energyRangeAll
   real    (kind=double), allocatable, dimension (:,:,:) :: momMatrix
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: conjWaveMomSum
   complex (kind=double)                                 :: valeValeXMom
#else
   real    (kind=double), allocatable, dimension (:,:,:) :: conjWaveMomSumGamma
   real    (kind=double)                                 :: valeValeXMomGamma
#endif

   ! Compute the number of energy points to evaluate.
   numEnergyPoints = int((maxTransEnSIGE * 2) / deltaSIGE + 1)

   ! Define constants for normalizing the broadening gaussian.
   sigmaSqrtPi = sqrt(pi) * sigmaSIGE

   ! Obtain the state bounds for this kpoint.
   firstInit =   firstOccupiedState(currentKPoint,spinDirection)
   lastInit  =    lastOccupiedState(currentKPoint,spinDirection)
   firstFin  = firstUnoccupiedState(currentKPoint,spinDirection)
   lastFin   =  lastUnoccupiedState(currentKPoint,spinDirection)

   ! Get the number of energy states in the system.
   numEnergyStates = lastFin - firstInit + 1

   ! Allocate space for the resulting sigmaE and initialize it to zero during
   !   the first kpoint iteration.  Also create the energy scale and initialize
   !   it.  NOTE:  The sigmaEAccumulator
   if (.not. allocated(sigmaEAccumulator)) then
      allocate (sigmaEAccumulator (numEnergyPoints,dim3))
      sigmaEAccumulator (:,:) = 0.0_double

      allocate (energyScale (numEnergyPoints))
      do i = 1, numEnergyPoints
         energyScale(i) = -maxTransEnSIGE + deltaSIGE * (i-1)
      enddo
   endif

   ! Determine the range of components (xyz) that should be considered.
   if (xyzComponents == 0) then
      initComponent = 1
      finComponent = 3
   else
      initComponent = 1
      finComponent = 1
   endif

#ifndef GAMMA
   ! Allocate space to hold the sum(conjg(valeVale(:,j)) * valeVale_Mom(:,k,1))
   !   for each of the possible energy states.
   allocate (conjWaveMomSum (valeDim,numEnergyStates,finComponent))

   ! Compute the sum over all the energy states.
   do i = initComponent,finComponent
      finalStateIndex = 0
      do j = firstInit, lastFin
         ! Increment the final state index for conjWaveMomSum
         finalStateIndex = finalStateIndex + 1
         do k = 1, valeDim
            conjWaveMomSum(k,finalStateIndex,i) = &
                  & sum(conjg(valeVale(:,j,1,1) * valeValeMom(:,k,1,i)))
         enddo
      enddo
   enddo

   ! Deallocate the momentum matrix elements
   deallocate (valeValeMom)

#else

   allocate (conjWaveMomSumGamma (valeDim,numEnergyStates,finComponent))

   ! Compute the sum over all the energy states.
   do i = initComponent,finComponent

      finalStateIndex = 0
      do j = firstInit, lastFin

         ! Increment the final state index for conjWaveMomSum
         finalStateIndex = finalStateIndex + 1

         do k = 1, valeDim
            conjWaveMomSumGamma(k,finalStateIndex,i) = &
                  & sum(valeValeGamma(:,j,1) * valeValeMomGamma(:,k,i))
         enddo
      enddo
   enddo

   ! Deallocate the momentum matrix elements
   deallocate (valeValeMomGamma)
#endif


   ! Allocate space for the fully computed momentum matrix.
   allocate (momMatrix (finComponent,numEnergyStates,numEnergyStates))

   ! Begin the double loop to determine the momentum matrix elements between
   !   each possible state.
   initialStateIndex = 0
   do i = firstInit, lastFin
      finalStateIndex = 0
      do j = firstInit, lastFin

         ! We don't care about self interactions so we skip the i=j case.
         if (i == j) then
            cycle
         endif

         ! Define the indices for the conjWaveMomSum and the momentum matrix
         !   elements.
         initialStateIndex = initialStateIndex + 1
         finalStateIndex   = finalStateIndex + 1

#ifndef GAMMA
         ! Loop to obtain the wave function times the momentum integral.
         do k = initComponent,finComponent
             valeValeXMom = sum(valeVale(:,i,1,1) * &
                   & conjWaveMomSum(:,finalStateIndex,k))

            ! Compute the imaginary component of the square of the
            !   valeValeXMom to obtain the momentum matrix element.  Note
            !   that the reason it is imaginary is because of the negative
            !   sign included in the getIntgResults subroutine for the
            !   momentum matrix.
            momMatrix(k,initialStateIndex,finalStateIndex) = &
                  &  real(valeValeXMom,double)**2 + &
                  & aimag(valeValeXMom)**2
         enddo
#else

         ! Loop to obtain the wave function times the momentum integral.
         do k = initComponent, finComponent
            valeValeXMomGamma = sum(valeValeGamma(:,i,1) * &
                  & conjWaveMomSumGamma(:,finalStateIndex,k))

            ! Compute the real component of the square of the valeValeXMom to
            !   obtain the momentum energy????
            momMatrix(k,initialStateIndex,finalStateIndex) = &
                  & valeValeXMomGamma**2
         enddo
#endif
      enddo
   enddo

   ! Allocate space for one energy range for each initial state and one energy
   !   range for each final state.
   allocate (energyRangeAll (numEnergyPoints, numEnergyStates))

   ! Determine the weighting effect of this kpoint and include the normalizaion
   !   factor for the gaussian.
   kPointFactor = kPointWeight(currentKPoint) * 0.5_double / sigmaSqrtPi 


   ! Define a gaussian of unit area for each energy state and compute its value
   !   for each energy point.  Note that the energy range is in hartree atomic
   !   units.
   do i = 1, numEnergyStates
      do j = 1, numEnergyPoints
         ! Determine the energy difference between the current state and the
         !   current point on the energy scale.
         broadenEnergyDiff = energyEigenValues(firstInit+i-1,currentKPoint,&
               & spinDirection) - energyScale(j)

         ! Determine the exponential alpha for the broadening factor.
         expAlpha = broadenEnergyDiff * broadenEnergyDiff / &
               & (sigmaSIGE * sigmaSIGE)

         ! Compute the value of the gaussian broadened energy state at this
         !   point on the energy scale.
         energyRangeAll(j,i) = exp(-expAlpha) * kPointFactor
      enddo
   enddo

   ! Multiply the gaussian for each state by the gaussian for each other state
   !   and the momentum matrix factor between the two states.
   if (xyzComponents == 0) then
      do i = 1, numEnergyStates
         do j = 1, numEnergyStates

            ! Ignore self interactions
            if (i == j) then
               cycle
            endif

            do k = 1,dim3
               sigmaEAccumulator(:,k) = sigmaEAccumulator(:,k) + &
                  & energyRangeAll(:,i) * energyRangeAll(:,j) * &
                  & momMatrix(k,i,j)
            enddo
         enddo
      enddo
   else
      do i = 1, numEnergyStates
         do j = 1, numEnergyStates

            ! Ignore self interactions
            if (i == j) then
               cycle
            endif

            sigmaEAccumulator(:,xyzComponents) = &
               & sigmaEAccumulator(:,xyzComponents) + &
               & energyRangeAll(:,i) * energyRangeAll(:,j) * &
               & momMatrix(1,i,j)
         enddo
      enddo
   endif

   ! Print the results during the last KPoint iteration for the last component
   !   or set of components.
   if ((currentKPoint == numKPoints) .and. ((xyzComponents == 0) .or. &
         & (xyzComponents == 3))) then

      ! This is a bit tricky and so it will help a lot to know something
      !   about cgs esu units before you start.  The Kubo-Greenwood formula
      !   (KGF) given in this reference:  Greenwood, Proc. Phys. Soc. London,
      !   71,585, (1958) is in cgs esu units.  We did the calculation in a.u.
      !   and the result must be presented in SI.  Arrrg.  Note some deviations
      !   from the formula as presented and as we calculated.  (1)  The
      !   formula used the velocity operator while this calculation used the
      !   momentum operator, so our equations have the electron mass squared in
      !   the denominator to compensate.  (2)  The formula should be considered
      !   as a vector outer product that produces a rank-2 tensor.  Our
      !   computed result is an x, y, z vector that is then averaged when it is
      !   printed out.  That is where the 1/3 factor comes from in our equation.

      ! The equation used is as follows:
      !    2 * Pi * hbar * e^2 / (3 * m^2 * Omega) *
      !    Sum ( |Pnm|^2 * delta(E-En) * delta(E-Em))

      ! In a.u. hbar = 1, e = 1, and m = 1 (They still have units though)

      ! In a.u. we have E=energy, T=time, M=mass, L=length.

      ! The first step is to apply the terms in the coefficient that are not
      !   equal to one.  (Note that we do not divide by three here since we
      !   will perform the averaging later.  We also do not change the units of
      !   the cell volume since that will be accounted for next.)
      conversionFactor = 2 * pi / realCellVolume


      ! The units of the KGF in a.u. are seen as:
      ! SigmaE = (ET EL) / (M^2 L^3) * (M^2 L^2)/(T^2) * 1/(E^2)
      ! SigmaE = 1/T
      ! hBar is ET, charge^2 (Q^2) is ML^3/T^2 = EL, we divide by the volume
      !   explicitly unlike in the KGF.  The tricky part is the charge because
      !   it is not a separate fundamental unit in cgs.  In SI it is in a
      !   separate fundamental unit called Coulombs, but in cgs it is not.  In
      !   cgs the units of Q are sqrt(g cm^3 / s^2).

      ! The second step is then to convert the 1/T in a.u. to 1/sec in cgs.
      !   This can be done by the conversion factor of 2.418884326505x10-17 s
      !   equals one atomic unit of time taken from NIST (auTime).  Note that
      !   the factor of 1e-17 is not included in the constant auTime and it
      !   must be accounted for later.

      ! This will convert the result from (a.u. 1/T) to (cgs emu 1/s).  Note
      !   that it still needs to be multiplied by 1e+17.
      conversionFactor = conversionFactor / auTime

      ! The third step is then to convert the 1/sec in cgs to 1/(ohm m) in
      !   SI.  This can be done with the knowledge that ohm = ET/Q^2 = Js/C^2.
      !   Writing that in cgs esu we have:  (g cm^2 / s^2) s s^2/(g cm^3) = 
      !   s / cm.  Then 1/(ohm cm) = cm / (s cm) = 1/s in cgs units.  This
      !   shows that the two units are equivalent.  Now we just need the
      !   conversion factor between them.
      !   We have from Jackson:  1 / (ohm m) = (1e-16 * c^2) * 1e9 1/s
      !   Therefor: 1/s (cgs emu) = 1/(1e-16 * c^2) * 1e-9 1/(ohm m)
      !   The value of 1e-16 * c^2 is DEFINED as:  8.9875517873681764.

      ! This will convert the result from (cgs emu 1/s) to (SI 1/(ohm m)).
      !   Note that it still needs to be multiplied by an additional 1e-9.
      !   This creates a total exponent multiplication factor of 1e+8 so far.
      conversionFactor = conversionFactor / lightFactor

      ! Finally, we want to express the result in (micro ohm cm)^-1 which is
      !   equal to (1e-8 ohm m)^-1 = 1e8(ohm m)^-1.  This creates a total
      !   exponent multiplication factor of 1 so we don't have to make any more
      !   modifications.

      ! Adjust the sigmaE to have the correct units.
      sigmaEAccumulator(:,:) = sigmaEAccumulator(:,:) * conversionFactor

      do i = 1, numEnergyPoints

         ! Write the averaged total, and the x, y, z components, making sure to
         !   convert the energy scale into eV.
         write (49+spinDirection,fmt="(5e20.8e4)") energyScale(i)*hartree,&
            & sum(sigmaEAccumulator(i,:))/3.0_double,sigmaEAccumulator(i,:)
      enddo

      deallocate (sigmaEAccumulator)
      deallocate (energyScale)
   endif

   ! Deallocate unnecessary matrices.
#ifndef GAMMA
   deallocate (conjWaveMomSum)
#else
   deallocate (conjWaveMomSumGamma)
#endif
   deallocate (momMatrix)
   deallocate (energyRangeAll)

end subroutine computeSigmaE


end module O_OptcTransitions
