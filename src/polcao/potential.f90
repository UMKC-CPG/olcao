module O_Potential

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define variables for the solid state potential.
   integer :: numAlphas ! The total number of potential terms (alphas) in the
         !   system.
   integer :: potDim ! Total potential dimension determined as the sum number
         !   of alphas for all potential types.
   real (kind=double), allocatable, dimension (:)   :: potAlphas ! An ordered
         !   list of all the alphas of all types in the system.
   real (kind=double), allocatable, dimension (:)   :: intgConsts ! A list of
         !   constants of integration for the terms of each type in the system.
   real (kind=double), allocatable, dimension (:)   :: typeSpinSplit ! The
         !   initial splitting used for each type.  This will be applied to
         !   every term in potDim as the spinSplitFactor below.
   real (kind=double), allocatable, dimension (:)   :: spinSplitFactor ! A
         !   factor -1.0 to 1.0 for each potential alpha (term) that is used
         !   to create the initial spin splitting kick.  The number (x say) is
         !   used as follows:  (up-down) = (up+down)*x.
   real (kind=double), allocatable, dimension (:,:) :: potCoeffs ! An ordered
         !   list of the potential coefficients for each spin orientation.
         !   Index1=coeff; Index2=spin(1=up,2=dn)

   ! Define variables for controlling convergence of the potential.
   integer :: feedbackLevel
   integer :: lastIteration
   integer :: xcCode
   integer :: currIteration
   integer :: converged
   real (kind=double) :: relaxFactor
   real (kind=double) :: convgTest

   integer :: spin
   integer :: rel
   integer :: GGA

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


! Set up data structures and values associated with the total potential
!   function.
subroutine initPotStructures

   ! Use necessary modules.
   use O_Kinds
   use O_Constants, only: pi
   use O_PotTypes, only: numPotTypes, potTypes

   implicit none

   ! Define local variables.
   integer :: i,j
   real (kind=double) :: pi32  ! pi^(3/2)

   ! Initialize variables
   pi32            = pi ** 1.5_double
   potDim          = 0

   ! Compute the value for potDim.
   potDim = potTypes(numPotTypes)%cumulAlphaSum+potTypes(numPotTypes)%numAlphas

   ! Allocate space for the spin split factor, integration constants, and
   !   actual alphas in one long ordered list.
   allocate (spinSplitFactor (potDim))
   allocate (intgConsts      (potDim))
   allocate (potAlphas       (potDim))

   ! Compute the constants of integration, copy the potential alphas into the
   !   total list from the individual types, and copy the spin split factor
   !   from each type to all terms of that type.
   do i = 1, numPotTypes

      ! Calculate the integration constants where the constant (j) is the
      !   multiplicity * (pi^(3/2)) / alpha(j) * alpha(j)^1/2.
      ! Also store the alphas of all the types in one array.
      do j = 1, potTypes(i)%numAlphas

         spinSplitFactor(j + potTypes(i)%cumulAlphaSum) = typeSpinSplit(i)

         potAlphas(j + potTypes(i)%cumulAlphaSum) = potTypes(i)%alphas(j)

         intgConsts(j + potTypes(i)%cumulAlphaSum) = &
               & real(potTypes(i)%multiplicity,double) * &
               & pi32 / (potTypes(i)%alphas(j))**(3.0_double/2.0_double)
      enddo
   enddo

   ! Deallocate arrays that will not be used further.
   deallocate (typeSpinSplit)

end subroutine initPotStructures

subroutine setPotControlParameters (fbL,lastIt,corrCode,rlxFact,cTest,&
      & typeSpinSplitTemp)

   ! Use necessary modules.
   use O_Kinds

   implicit none

   ! Define dummy variables.
   integer :: fbL
   integer :: lastIt
   integer :: corrCode
   real (kind=double) :: rlxFact
   real (kind=double) :: cTest
   real (kind=double), dimension(:) :: typeSpinSplitTemp

   integer :: i
   integer :: info
   integer :: numCodes
   integer :: xcCodeParam
   integer :: spinParam
   integer :: relParam
   integer :: GGAParam
   character(len=25) :: functionalName
   character(len=100) :: dataDirectory
   character(len=100) :: xcCodeDataFile

   ! Set control parameters read from input.
   feedbacklevel   = fbL
   lastIteration   = lastIt
   xcCode          = corrCode
   relaxFactor     = rlxFact
   convgTest       = cTest

   ! Allocate space to hold the initially read in spin splitting for each type.
   allocate (typeSpinSplit(size(typeSpinSplitTemp)))

   typeSpinSplit(:) = typeSpinSplitTemp(:)

   ! Initialize control parameters.
   currIteration = 1
   converged     = 0

   ! open the xc_code.dat file
   call get_environment_variable("OLCAO_DATA", dataDirectory, info)
   xcCodeDataFile=trim(dataDirectory)//"/xc_code.dat"

   open (unit=9, file=xcCodeDataFile)

   ! Read past the header line and then get the total number of codes.
   read (9,*)
   read (9,*) numCodes

   ! Read each line(functional type). Once the correct one is found use
   ! it to set "spin", "rel", "GGA".
   do i = 1, numCodes
      read (9,*) functionalName, xcCodeParam, spinParam, &
                                     & relParam, GGAParam
      if (xcCodeParam.eq.xcCode) then
         spin = spinParam
         rel = relParam
         GGA = GGAParam
      endif
   enddo

end subroutine setPotControlParameters

subroutine initPotCoeffs

   ! Include the necessary modules
   use O_Kinds
   use O_PotTypes, only: numPotTypes, potTypes
   use MPI

   ! Define local variables
   integer :: i,j,k ! Loop index variables
   integer :: potTermCount
   real (kind=double) :: spaceHolder1 ! Gaussian exponential alphas
   real (kind=double) :: spaceHolder2 ! Total charge density
   real (kind=double) :: spaceHolder3 ! Valence charge density
   real (kind=double) :: spaceHolder4 ! Up - Down valence charge density.
   integer :: mpiRank,mpiErr

   ! Get MPI rank.
   call MPI_COMM_RANK (MPI_COMM_WORLD,mpiRank,mpierr)

   ! Allocate space for the potential coefficients.
   allocate (potCoeffs(potDim,spin))

   if (mpiRank == 0) then

      ! Open the potential file that will be read from in this program.
      open (unit=8,file='fort.8',status='old',form='formatted')
   
      ! Read the existing potential coefficients for each term, or for the
      !   spin up and then spin down terms separately if we are doing a spin
      !   polarized calculation.
      read (8,*)  ! File header that says number of types.
      do i = 1, spin
   
         ! Initialize the counter of potential terms.
         potTermCount = 0
   
         read (8,*)  ! Read tag indicating spin up or spin down.
         do j = 1, numPotTypes
            read (8,*)  ! Type header saying the number of terms for this type.
            do k = 1, potTypes(j)%numAlphas
   
               ! Increment the counter.
               potTermCount = potTermCount + 1
   
               read (8,*) potCoeffs(potTermCount,i),spaceHolder1,spaceHolder2,&
                  & spaceHolder3,spaceHolder4 
            enddo
         enddo
      enddo

      ! Close the potential file.
      close (8)
   endif

   ! Get everyone together
   MPI_BARRIER(MPI_COMM_WORLD,mpiErr)

   ! Broadcast the potential coefficients to the other processes.
   do i = 1, spin
      MPI_BCAST(potCoeffs(:,i),potDim,MPI_DOUBLE_PRECISION,0,&
            & MPI_COMM_WORLD,mpiErr)
   enddo

end subroutine initPotCoeffs


subroutine cleanUpPotential

   implicit none

   deallocate (spinSplitFactor)
   deallocate (intgConsts)
   deallocate (potAlphas)

   ! This is allocated for main and intg but not other programs.
   if (allocated(potCoeffs)) then
      deallocate (potCoeffs)
   endif

end subroutine cleanUpPotential


end module O_Potential
