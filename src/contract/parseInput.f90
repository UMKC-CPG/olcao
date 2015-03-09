module ParseInputSubs

public

contains

subroutine parseContractInput

   ! Use the necessary program parameters
   use O_Kinds

   ! Import the necessary subroutine modules.
   use O_ReadDataSubs

   ! Import the necessary data modules.
   use ExecutionData
   use AtomData
   use GaussianBasisData
   use PotGaussianData

   ! Import the necessary modules.
   use O_RadialGrid
   use ChargeDensityMod

   implicit none

   integer    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Define local variables.
   integer :: i
   character*8 :: commandBuffer

   ! Read the command line parameter for plotting the basis numerically and
   !   creating a POVRAY scene file.
   call getarg(1,commandBuffer)
   read (commandBuffer,*) doVis

   ! Open the file needed to read the input.
   open(44,file='gauss.fit',status='old',form='formatted')
   open(5,file='contract.dat',status='old',form='formatted')
   readUnit = 5

   ! Open the file that will be written to as output for this program
   open(20,file='contract.out',status='new',form='formatted')
   writeUnit = 20

   ! Read the name of the element abbreviated from the periodic table.
   call readData(readUnit,writeUnit,len(elementName),elementName,&
         & len('ELEMENT_NAME'),'ELEMENT_NAME')

   ! Initialize the count of orbitals.
   numCoreOrbitals(:) = 0
   numValeOrbitalsPerBasis(:,:) = 0

   ! Read the list of core orbitals and valence orbitals for each basis set
   !   beyond the previous set for this atom.
   call readData(readUnit,writeUnit,4,numCoreOrbitals(:),&
         & len('NUM_CORE_ORBITALS'),'NUM_CORE_ORBITALS')
   call readData(readUnit,writeUnit,4,numValeOrbitalsPerBasis(:,1),&
         & len('NUM_VALE_ORBITALS_MB'),'NUM_VALE_ORBITALS_MB')
   call readData(readUnit,writeUnit,4,numValeOrbitalsPerBasis(:,2),&
         & len('NUM_VALE_ORBITALS_FB'),'NUM_VALE_ORBITALS_FB')
   call readData(readUnit,writeUnit,4,numValeOrbitalsPerBasis(:,3),&
         & len('NUM_VALE_ORBITALS_EB'),'NUM_VALE_ORBITALS_EB')

   ! Compute the total number of valence orbitals as the sum of the
   !   contributions of the individual basis additions.
   do i = 1,4
      numValeOrbitals(i) = sum(numValeOrbitalsPerBasis(i,:))
   enddo

   ! Compute the total orbitals for each lQN for this atom.
   numTotalOrbitals(:) = numCoreOrbitals(:) + numValeOrbitals(:)

   ! Read the max number of Gaussians to use for the orbitals (each orbital
   !   angular momentum can not use more than this although they can use less).
   call readData(readUnit,writeUnit,maxNumBasisGaussians,&
         & len('MAX_NUM_GAUSSIANS'),'MAX_NUM_GAUSSIANS')

   ! Read the minimum and maximum gaussian alphas.
   call readData(readUnit,writeUnit,minBasisAlpha,maxBasisAlpha,&
         & len('MIN_MAX_ALPHAS'),'MIN_MAX_ALPHAS')

   ! Read the atomic number and nuclear atomic alpha for this atom.
   call readData(readUnit,writeUnit,atomicNumber,len('ATOMIC_NUMBER'),&
         & 'ATOMIC_NUMBER')
   call readData(readUnit,writeUnit,nuclearAlpha,len('NUCLEAR_ALPHA'),&
         & 'NUCLEAR_ALPHA')

   ! Allocate space to hold the list of which gaussians should be selected for
   !    each of the s,p,d,f orbitals.
   allocate (selectGaussians(maxNumBasisGaussians,4))
   selectGaussians(:,:) = 0

   ! Read the list of gaussians to use for the s,p,d,f orbitals.
   call readOrbitalSelection(readUnit,writeUnit,1,'S_GAUSSIAN_LIST')
   call readOrbitalSelection(readUnit,writeUnit,2,'P_GAUSSIAN_LIST')
   call readOrbitalSelection(readUnit,writeUnit,3,'D_GAUSSIAN_LIST')
   call readOrbitalSelection(readUnit,writeUnit,4,'F_GAUSSIAN_LIST')


   ! Read the atomic potential function in its gaussian form.
   read (44,*) numPotGaussians
   numPotGaussians = numPotGaussians + 1 ! Need to include the nuclear pot.
   allocate (potAlphas(numPotGaussians))
   allocate (potCoeffs(numPotGaussians))
   do i = 2, numPotGaussians
      read (44,*) potCoeffs(i), potAlphas(i)
   enddo

   ! Insert the nuclear potential gaussian as the first term in the atomic
   !   potential function.
   potAlphas(1) = nuclearAlpha
   potCoeffs(1) = -1.0_double * atomicNumber



   ! Read data necessary for plotting the charge density and for comparison to
   !   the charge calculated with the atomSCF program.  (If requested.)
   call readData(readUnit,writeUnit,doCharge,len('CALCULATE_CHARGE'),&
         & 'CALCULATE_CHARGE')
   if (doCharge == 1) then

      ! Read the radial grid parameters.
      call readData(readUnit,writeUnit,radialMaxDist,aaWhatever,bbWhatever,&
            & len('RADIAL_GRID_PARAMETERS'),'RADIAL_GRID_PARAMETERS')

      ! Allocate space to hold charge present in each orbital and initialize it.
      allocate (orbitalCharge(4,numTotalOrbitals(1)))
      orbitalCharge (:,:) = 0.0_double

      ! Read the s,p,d,f orbital occupations.
      call readOrbitalOccupations (readUnit,writeUnit,1,'S_ORBITAL_OCCUPATIONS')
      call readOrbitalOccupations (readUnit,writeUnit,2,'P_ORBITAL_OCCUPATIONS')
      call readOrbitalOccupations (readUnit,writeUnit,3,'D_ORBITAL_OCCUPATIONS')
      call readOrbitalOccupations (readUnit,writeUnit,4,'F_ORBITAL_OCCUPATIONS')
   endif

   ! Close the input files.
   close (5)
   close (44)

   ! ALL DONE READING

end subroutine parseContractInput


subroutine readOrbitalSelection(readUnit,writeUnit,QN_l,tag)

   ! Use the necessary program parameters
   use O_Kinds

   ! Import the necessary subroutine modules.
   use O_ReadDataSubs

   ! Use necessary data modules.
   use GaussianBasisData
   use AtomData
   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Define local variables.
   integer :: QN_l
   character*15 :: tag

   if (numTotalOrbitals(QN_l) /= 0) then
      call readData(readUnit,writeUnit,maxNumBasisGaussians,&
            & selectGaussians(:,QN_l),15,tag)
   endif

end subroutine readOrbitalSelection


subroutine readOrbitalOccupations(readUnit, writeUnit,QN_l,tag)

   ! Use the necessary program parameters
   use O_Kinds

   ! Import the necessary subroutine modules.
   use O_ReadDataSubs

   ! Use necessary data modules.
   use ChargeDensityMod
   use AtomData

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Define local variables.
   integer :: QN_l
   character*21 :: tag

   if (numTotalOrbitals(QN_l) /= 0) then
      call readData(readUnit,writeUnit,numTotalOrbitals(QN_l), &
               & orbitalCharge(QN_l,1:numTotalOrbitals(QN_l)),21,tag)
   endif

end subroutine readOrbitalOccupations


subroutine implicitInput

   ! Use the necessary program parameters
   use O_Kinds

   ! Import the necessary data modules.
   use AtomData
   use PotGaussianData
   use GaussianBasisData
   use ExecutionData

   ! Import the necessary external modules.
   use O_RadialGrid

   implicit none

   ! Define local variables.
   integer :: i,j,k
   integer :: valeOrbitalCount
   real (kind=double) :: expTerm

   ! Set up the radial grid.
   if (doCharge == 1) then
      call setupRadialGrid(atomicNumber)
   endif

   ! Allocate space to hold the basis alphas and the basis ID numbers to say
   !   which basis each orbital belongs to.
   allocate (basisAlphas(maxNumBasisGaussians))
   allocate (basisID(maxval(numValeOrbitals(:)),4))

   ! Compute the number of s,p,d,f basis gaussians.
   numBasisGaussians(:)=0
   do i = 1, 4
      do j = 1, maxNumBasisGaussians
         if (selectGaussians(j,i) == 1) then
            numBasisGaussians(i) = numBasisGaussians(i) + 1
         endif
      enddo
   enddo

   ! Compute the number of s,p,d,f basis terms that must be printed.  This is
   !   the same as above + the number of terms with coeff=0 followed with a
   !   non-zero term at *any* higher index.  This is only needed for printing
   !   the final results.
   numBasisTerms(:)=0
   do i = 1,4
      do j = maxNumBasisGaussians,1,-1
         if (selectGaussians(j,i) == 1) then
            numBasisTerms(i) = j
            exit
         endif
      enddo
   enddo

   ! Compute the atomic gaussian alphas.
   expTerm = (maxBasisAlpha/minBasisAlpha) ** &
         & (1.0_double/real(maxNumBasisGaussians-1,double))
   do i = 1, maxNumBasisGaussians
      basisAlphas(i) = minBasisAlpha*expTerm**(i-1)
   enddo

   ! Compute the basis ID numbers for each orbital in the valence basis set.
   basisID(:,:) = 0
   do i = 1,4 ! s,p,d,f
      valeOrbitalCount = 0
      do j = 1,3 ! MB, FB, EB
         do k = 1, numValeOrbitalsPerBasis(i,j)
            valeOrbitalCount = valeOrbitalCount + 1
            basisID(valeOrbitalCount,i) = j
         enddo
      enddo
   enddo

end subroutine implicitInput

end module ParseInputSubs
