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

      ! Allocate memory to hold the data for those points.
      allocate (tdos(numEnergyPoints,spin,lastIteration))
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
               expTerm = ((currentEnergyValues(k) - &
                     & energyScale(l))/sigmaDOS)**2

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
                  tdos(l,i,currIteration) = tdos(l,i,currIteration) + &
                        & expFactor
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
      fileID = 999+i

      ! Define the file name.
      write (fileName,fmt="(a5,i4)") "fort.",fileID

      ! Open the file for creating the opendx output.
      open (unit=fileID,file=fileName,status='unknown',form='formatted')


      ! Define the mesh
!      write (fileID,*) "object 1 class gridpositions counts ",numEnergyPoints,&
!            & currIteration
!      write (fileID,*) "origin -30 0"
!      write (fileID,*) "delta ",dosInput%deltaDOS,0.0_double
!      write (fileID,*) "delta ",0.0_double,1.0_double
!      write (fileID,*)


      ! Define the connections
!      write (fileID,*) "object 2 class gridconnections counts ",&
!            & numEnergyPoints,currIteration
!      write (fileID,*)


      ! Define the data
!      write (fileID,*) "object 3 class array type float rank 0 items ",&
!            & numEnergyPoints*currIteration," data follows"

do j = 1, numEnergyPoints
! Make sure to convert to eV upon output.
write (fileID,fmt="(1x,3f12.8)") energyScale(j)*hartree,tdos(j,1,1),&
      & tdos(j,1,currIteration-1)

enddo

      ! Write the data.
!      pointCount = 0
!      do j = 1, numEnergyPoints
!         do k = 1, currIteration
!
!            ! Record the point
!            pointCount = pointCount + 1
!            write (fileID,ADVANCE="NO",fmt="(1x,f12.8)") tdos(j,i,k)
!
!            ! Write a newline after every 5 points
!            if (pointCount == 5) then
!               write (fileID,*)
!               pointCount = 0
!            endif
!         enddo
!      enddo

      ! Add a final newline if necessary.
!      if (pointCount /= 0) then
!         write (fileID,*)
!      endif

!      ! Tack on the ending information.
!      write (fileID,*) 'attribute "dep" string "positions"'
!      write (fileID,*)   
!      write (fileID,*) 'object "tdos" class field'
!      write (fileID,*) 'component "positions" value 1'
!      write (fileID,*) 'component "connections" value 2'
!      write (fileID,*) 'component "data" value 3'
!      write (fileID,*) 'end'
!      write (fileID,*)   
   enddo


   ! Deallocate unused memory.
   deallocate (tdos)
   deallocate (energyScale)
   deallocate (currentEnergyValues)

end subroutine printIterationTDOS


subroutine computeDOS

   ! Import the necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Potential,   only: spin
   use O_CommandLine, only: doDOS
   use O_Populate,    only: electronPopulation
   use O_KPoints,     only: numKPoints, kPointWeight
   use O_Constants,   only: pi, hartree, maxOrbitals
   use O_AtomicSites, only: valeDim, numAtomSites, atomSites
   use O_AtomicTypes, only: numAtomTypes, atomTypes, maxNumValeStates
   use O_Input,       only: numStates, sigmaDOS, eminDOS, emaxDOS, deltaDOS, &
         & detailCodePDOS
#ifndef GAMMA
   use O_SecularEquation, only: valeValeOL, valeVale, energyEigenValues, &
         & readDataSCF, readDataPSCF
#else
   use O_SecularEquation, only: valeValeOLGamma, valeValeGamma, &
         & energyEigenValues, readDataSCF, readDataPSCF
#endif


   ! Make sure that no funny variables are used.
   implicit none

   ! Define local variables.
   integer :: h,i,j,k,l ! Loop index variables
   character*17 :: formatString
   character*1, dimension (maxOrbitals) :: QN_lLetter
   character*14, dimension (4,7) :: QN_mLetter
   integer, allocatable, dimension (:) :: cumulNumDOS
   integer, allocatable, dimension (:) :: pdosIndex ! Given a valence state,
         ! this stores the place within pdosAccum that this state should be
         ! stored.
   integer, allocatable, dimension (:) :: numAtomStates
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
   real (kind=double), allocatable, dimension (:)   :: pdosAccum
   real (kind=double), allocatable, dimension (:)   :: localizationIndex
   real (kind=double), allocatable, dimension (:)   :: totalSystemDos
   real (kind=double), allocatable, dimension (:)   :: totalTypeDos
   real (kind=double), allocatable, dimension (:)   :: totalAtomDos
   real (kind=double), allocatable, dimension (:)   :: energyValuesAvg
   real (kind=double), allocatable, dimension (:,:) :: pdosComplete
   real (kind=double), allocatable, dimension (:,:) :: electronNumber
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:) :: waveFnSqrd
#else
   real (kind=double), allocatable, dimension (:)    :: waveFnSqrdGamma
#endif


   ! Log the date and time we start.
   call timeStampStart (19)

   ! Normalization for Gaussian broadening.
   sigmaSqrtPi = sqrt(pi) * sigmaDOS

   ! Allocate arrays and matrices for this computation.
   allocate (pdosIndex       (valeDim))
   allocate (numAtomStates   (numAtomSites))
   ! Store DOS for each type's orbital sum.
   if     (detailCodePDOS == 0) then
      allocate (cumulNumDOS  (numAtomTypes + 1))
      allocate (pdosAccum    (valeDim))
   ! Store DOS for each atom's TDOS.
   elseif (detailCodePDOS == 1) then
      allocate (cumulNumDOS  (1)) ! Unused for this detailCodePDOS.
      allocate (pdosAccum    (numAtomSites))
   ! Store DOS for each QN_nl resolved atom.
   elseif (detailCodePDOS == 2) then
      allocate (cumulNumDOS  (numAtomSites + 1))
      allocate (pdosAccum    (valeDim))
   ! Store DOS for each QN_nlm resolved atom.
   elseif (detailCodePDOS == 3) then
      allocate (cumulNumDOS  (numAtomSites + 1))
      allocate (pdosAccum    (valeDim))
   endif
#ifndef GAMMA
   allocate      (waveFnSqrd (valeDim))
   allocate      (valeValeOL (valeDim,valeDim,1,1))   !Holds overlap.
   if (doDOS == 0) then  ! If true then we are not in the SCF phase.
      allocate   (valeVale(valeDim,numStates,1,1)) !Holds wave functions.
   endif
#else
   allocate      (waveFnSqrdGamma (valeDim))
   allocate      (valeValeOLGamma (valeDim,valeDim,1))   !Holds overlap.
   if (doDOS == 0) then  ! If true then we are not in the SCF phase.
      allocate   (valeValeGamma(valeDim,numStates,1)) !Holds wave fns
   endif
#endif

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
   if     (detailCodePDOS == 0) then ! Save DOS according to orbitals
      allocate (pdosComplete      (cumulDOSTotal,numEnergyPoints))
      allocate (totalTypeDos      (numEnergyPoints))
   elseif (detailCodePDOS == 1) then ! Save DOS according to atoms
      allocate (pdosComplete      (numAtomSites,numEnergyPoints))
   elseif (detailCodePDOS == 2) then ! Save DOS by to atoms & orbitals
      allocate (pdosComplete      (cumulDOSTotal,numEnergyPoints))
      allocate (totalAtomDos      (numEnergyPoints))
   elseif (detailCodePDOS == 3) then ! Save DOS by select atoms and lm QNs.
      allocate (pdosComplete      (cumulDOSTotal,numEnergyPoints)) !valeDim
      allocate (totalAtomDos      (numEnergyPoints))
   endif

   ! Assign values to the energy scale.
   do i = 1, numEnergyPoints
      energyScale(i) = eminDOS + (i-1) * deltaDOS
   enddo

   do h = 1, spin

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

      ! Begin accumulating the DOS values
      do i = 1, numKPoints


         ! Determine if we are doing the DOS in a post-SCF calculation, or
         !   within an SCF calculation.  (The doDOS flag is only set to true
         !   within an SCF calculation because that is the only place that this
         !   option makes sense to use.  In a post-SCF case, the DOS job was
         !   explicitly asked for by the user and so this request flag is not
         !   used.)
         if (doDOS == 1) then
            ! Read necessary data from SCF (setup,main) data structures.
            call readDataSCF(h,i,numStates)
         else
            ! Read necessary data from post SCF (intg,band) data structures.
            call readDataPSCF(h,i,numStates)
         endif


         do j = 1, numStates

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
                        & conjg(valeVale(valeDimIndex,j,1,1)) * &
                        & valeVale(:valeDim,j,1,1)

                  ! Compute the effects of overlap for the real part only
                  !   (real*real) + (imag*imag).
                  oneValeRealAccum = sum( &
                       & real(waveFnSqrd(:valeDim),double) * &
                       & real(valeValeOL(:valeDim,valeDimIndex,1,1),double) + &
                       & aimag(waveFnSqrd(:valeDim)) * &
                       & aimag(valeValeOL(:valeDim,valeDimIndex,1,1)))
#else
                  ! Compute the square of the wave function for each element.
                  waveFnSqrdGamma(:valeDim) = &
                        & valeValeGamma(valeDimIndex,j,1) * &
                        & valeValeGamma(:valeDim,j,1)

                  ! Compute the effects of overlap.
                  oneValeRealAccum = sum(waveFnSqrdGamma(:valeDim) * &
                        & valeValeOLGamma(:valeDim,valeDimIndex,1))
#endif

                  ! Store the current electron number assignment.
                  electronNumber(l,k) = electronNumber(l,k) + &
                        & oneValeRealAccum * kPointWeight(i) * occupancyNumber


                  pdosAccum(pdosIndex(valeDimIndex)) = &
                        & pdosAccum(pdosIndex(valeDimIndex)) + &
                        & oneValeRealAccum / real(spin,double)


                  ! Store the square of the accumulation as the localization
                  !   index for this state.
                  localizationIndex(j) = localizationIndex(j) + &
                        & oneValeRealAccum * oneValeRealAccum * &
                        & kPointWeight(i) / real(spin,double)
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
               if (expTerm < 50) then

                  ! Compute the exponential factor.  It is at this point that
                  !   the kpoint weighting factor is applied to make the energy
                  !   bucket we are considering be filled with the right number
                  !   of electrons.  The pdosAccum below already has the
                  !   1/spin factor included.  Note that the hartree conversion
                  !   must be included here because the sigmaDOS was in units
                  !   of hartree, but the expTerm was not (because the energy
                  !   values were also in hartree).  This left a hartree factor
                  !   that was uncompensated for in the sigmaSqrtPi.
                  expFactor = exp(-expTerm) / sigmaSqrtPi / hartree * &
                        & kPointWeight(i)

                  ! Broaden and store the complete pdos
                  pdosComplete(:,k) = pdosComplete(:,k) + pdosAccum(:) * &
                        & expFactor
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


      ! Adjust the electron number for each orbital and state by the ratio of
      !   the exact number of electrons in the system (systemElectrons) to the
      !   computed number of electrons.
      electronFactor = currentPopulation / totalElectronsComputed
      electronNumber(:,:) = electronNumber(:,:) * electronFactor


      ! Normalize the computed total DOS and the orbital DOS by the number of
      !   states computed.
      totalSystemDos(:) = totalSystemDos(:) * electronFactor

      pdosComplete(:,:) = pdosComplete(:,:) * electronFactor


      ! Find the total area under the computed TDOS.  This should be equal to
      !   the number of electron spin states available in the system
      !   (occupied + unoccupied) OVER THE REQUESTED ENERGY RANGE ONLY!
      integratedArea = sum(totalSystemDos(1:numEnergyPoints - 1) + &
            & totalSystemDos(2:numEnergyPoints)) * deltaDOS * 0.5_double

      ! Record the exact values and calculated values for electrons and states.
      write (20,*) 'Electrons Calculated:    ',totalElectronsComputed
      write (20,*) 'Electrons Expected:      ',currentPopulation
      write (20,*) 'Spin States Calculated:  ',integratedArea
      write (20,*) 'Spin States Expected:    ',numStatesInRange


      ! Begin recording the results to disk.

      ! Record the total system DOS, converting the scale to eV.
      do i = 1, numEnergyPoints
         write (59+h,fmt="(f14.4,f14.6)") energyScale(i) * hartree,&
               & totalSystemDos(i)
      enddo

      ! Loop over types for the types-based DOS
      if (detailCodePDOS == 0) then

         ! Print the key bits of information for the PDOS output.
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
            do j = 1, maxOrbitals  ! 1=s; 2=p; 3=d; 4=f
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
               write (69+h,ADVANCE="NO",fmt="(e12.6,1x)") &
                  & sum(pdosComplete(initIndex:finIndex,j))
               numCols = 1
               do k = initIndex,finIndex
                  numCols = numCols + 1
                  write (69+h,ADVANCE="NO",fmt="(e12.6,1x)") pdosComplete(k,j)
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
            do j = 1, maxOrbitals
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
               write (69+h,ADVANCE="NO",fmt="(e12.6,1x)") &
                     & sum(pdosComplete(initIndex:finIndex,j))
               numCols = 1
               do k = initIndex,finIndex
                  numCols = numCols + 1
                  write (69+h,ADVANCE="NO",fmt="(e12.6,1x)") pdosComplete(k,j)
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
            do j = 1, maxOrbitals
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
               write (69+h,ADVANCE="NO",fmt="(e12.6,1x)") &
                     & sum(pdosComplete(initIndex:finIndex,j))
               numCols = 1
               do k = initIndex,finIndex
                  numCols = numCols + 1
                  write (69+h,ADVANCE="NO",fmt="(e12.6,1x)") pdosComplete(k,j)
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
      !   kpoints.  This energy value is used as the energy for the
      !   localization index.

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
         write (79+h,fmt="(1x,f14.8,2x,f14.8)") energyValuesAvg(i) * hartree,&
               & localizationIndex(i)
      enddo
   enddo ! (h spin)

   ! Deallocate all the unnecessary matrices and arrays
   if (detailCodePDOS == 0) then
      deallocate (totalTypeDos)
   elseif (detailCodePDOS == 2) then
      deallocate (totalAtomDos)
   elseif (detailCodePDOS == 3) then
      deallocate (totalAtomDos)
   endif
   deallocate (cumulNumDOS)
   deallocate (pdosIndex)
   deallocate (numAtomStates)
   deallocate (pdosAccum)
   deallocate (energyScale)
   deallocate (localizationIndex)
   deallocate (totalSystemDos)
   deallocate (pdosComplete)
   deallocate (electronNumber)
#ifndef GAMMA
   deallocate (waveFnSqrd)
   if (doDOS == 0) then
      deallocate (valeVale)
      deallocate (valeValeOL)
   endif
#else
   deallocate (waveFnSqrdGamma)
   if (doDOS == 0) then
      deallocate (valeValeGamma)
      deallocate (valeValeOLGamma)
   endif
#endif

   ! Log the date and time we end.
   call timeStampEnd (19)

end subroutine computeDOS


end module O_DOS
