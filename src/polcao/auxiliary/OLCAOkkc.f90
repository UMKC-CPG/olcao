program OLCAOkkc
!Written by Paul Rulis.  Adapted from EPS1 by Dr. Y.-N. Xu.
!Last modified on Apr. 11, 2007

!The purpose of the program is to calculate epsilon1 and the energy loss
! function from the given epsilon2 using the Kramers-Kronig conversion
! relation.

!The usual procedure for calculating the energy loss function (ELF) is as
! follows.  First form a new finer grained set of energies by linear
! interpolation.  Form a finer grained version of the epsilon2 data.  Perform
! Kramers Kronig conversion to determine epsilon1.  Use epsilon1 and epsilon2
! to find the ELF.

   use O_Kinds
   use O_Constants
   use MPI

   implicit none

   real (kind=double), allocatable, dimension(:) :: energy, fine_energy
   real (kind=double), allocatable, dimension(:) :: totalEps1,totalEps1i
   real (kind=double), allocatable, dimension(:) :: totalEps2
   real (kind=double), allocatable, dimension(:) :: totalElf,totalRefractIndex
   real (kind=double), allocatable, dimension(:,:) :: eps1,eps1i
   real (kind=double), allocatable, dimension(:,:) :: elf,refractIndex
   real (kind=double), allocatable, dimension(:,:) :: eps2,fine_eps2
   
   real (kind=double) :: coarseEnergyDiff,fineEnergyDiff
   integer :: numValues=0  ! This value will be set later to length*grain.
   integer :: grain=10     ! Grainularity factor for refining the integration.

   ! The variables that will be defined by a command line parameters
   integer :: length  ! Length of the optical conductivity file
   integer :: spin    ! 1 for spin up or total, 2 for spin down.

   character*40 :: buffer ! Character buffer for command line arguments

   ! Define global mpi parameters
   integer :: mpiRank
   integer :: mpiSize
   integer :: mpiErr

   ! Initialize the MPI interface
   call MPI_INIT (mpiErr)
   call MPI_COMM_RANK (MPI_COMM_WORLD,mpiRank,mpiErr)
   call MPI_COMM_SIZE (MPI_COMM_WORLD,mpiSize,mpiErr)

   ! Define the subroutine interfaces for proper passing of allocatable arrays.
   interface
      subroutine readOPTCData (length, energy, totalEps2, eps2, mpiRank)
         use O_Kinds
         use MPI
         integer, intent(in) :: length
         real (kind=double), dimension (:), intent(out) :: energy
         real (kind=double), dimension (:), intent(out) :: totalEps2
         real (kind=double), dimension (:,:), intent(out) :: eps2
         integer, intent(in) :: mpiRank
      end subroutine readOPTCData
      subroutine getCoarseEnergyDiff(length, energy, coarseEnergyDiff)
         use O_Kinds
         integer, intent(in) :: length
         real (kind=double), dimension (:), intent(in) :: energy
         real (kind=double), intent(out) :: coarseEnergyDiff
      end subroutine getCoarseEnergyDiff
      subroutine makeFineGrainEnergy(length,grain,energy,fine_energy)
         use O_Kinds
         integer, intent(in) :: length
         integer, intent(in) :: grain
         real (kind=double), dimension (:), intent(in)  :: energy
         real (kind=double), dimension (:), intent(out) :: fine_energy
      end subroutine makeFineGrainEnergy
      subroutine getFineEnergyDiff(numValues,grain,fine_energy,fineEnergyDiff)
         use O_Kinds
         integer, intent(in) :: numValues
         integer, intent(in) :: grain
         real (kind=double), dimension (:), intent(in) :: fine_energy
         real (kind=double), intent(out) :: fineEnergyDiff
      end subroutine getFineEnergyDiff
      subroutine makeFineEps2(length,grain,energy,fine_energy,eps2,&
                             &fine_eps2)
         use O_Kinds
         integer, intent(in) :: length
         integer, intent(in) :: grain
         real (kind=double), dimension (:), intent(in) :: energy
         real (kind=double), dimension (:), intent(in) :: fine_energy
         real (kind=double), dimension (:,:), intent(in) :: eps2
         real (kind=double), dimension (:,:), intent(out) :: fine_eps2
      end subroutine makeFineEps2
      subroutine kramersKronig(length,grain,numValues,mpiRank,mpiSize,&
                              &fineEnergyDiff,energy,fine_energy,eps1,eps1i,&
                              &fine_eps2)
         use O_Kinds
         use O_Constants
         use MPI
         integer, intent(in) :: length
         integer, intent(in) :: grain
         integer, intent(in) :: numValues
         integer, intent(in) :: mpiRank
         integer, intent(in) :: mpiSize
         real (kind=double), intent(in) :: fineEnergyDiff
         real (kind=double), dimension (:), intent(in) :: energy
         real (kind=double), dimension (:), intent(in) :: fine_energy
         real (kind=double), dimension (:,:), intent(out) :: eps1
         real (kind=double), dimension (:,:), intent(out) :: eps1i
         real (kind=double), dimension (:,:), intent(in) :: fine_eps2
      end subroutine kramersKronig
      subroutine getRefractIndex(length,eps1,eps2,refractIndex)
         use O_Kinds
         integer, intent(in) :: length
         real (kind=double), dimension (:,:), intent(in)  :: eps1
         real (kind=double), dimension (:,:), intent(in)  :: eps2
         real (kind=double), dimension (:,:), intent(out) :: refractIndex
      end subroutine getRefractIndex
      subroutine getELF(length,eps1,eps2,elf)
         use O_Kinds
         integer, intent(in) :: length
         real (kind=double), dimension (:,:), intent(in)  :: eps1
         real (kind=double), dimension (:,:), intent(in)  :: eps2
         real (kind=double), dimension (:,:), intent(out) :: elf
      end subroutine getELF
      subroutine averageFunctions(length,eps1,eps1i,refractIndex,&
                                 &totalEps1,totalEps1i,totalEps2,totalElf,&
                                 &totalRefractIndex)
         use O_Kinds
         integer, intent(in) :: length
         real (kind=double), dimension (:,:), intent(in)  :: eps1
         real (kind=double), dimension (:,:), intent(in)  :: eps1i
         real (kind=double), dimension (:,:), intent(in)  :: refractIndex
         real (kind=double), dimension (:),   intent(out) :: totalEps1
         real (kind=double), dimension (:),   intent(out) :: totalEps1i
         real (kind=double), dimension (:),   intent(inout) :: totalEps2
         real (kind=double), dimension (:),   intent(out) :: totalElf
         real (kind=double), dimension (:),   intent(out) :: totalRefractIndex
      end subroutine averageFunctions
      subroutine printData(length,energy,totalEps1,totalEps1i,totalEps2,&
                         & totalElf,totalRefractIndex,eps1,eps1i,elf,&
                         & refractIndex)
         use O_Kinds
         integer, intent(in) :: length
         real (kind=double), dimension (:), intent(in)   :: energy
         real (kind=double), dimension (:), intent(in)   :: totalEps1
         real (kind=double), dimension (:), intent(in)   :: totalEps1i
         real (kind=double), dimension (:), intent(in)   :: totalEps2
         real (kind=double), dimension (:), intent(in)   :: totalElf
         real (kind=double), dimension (:), intent(in)   :: totalRefractIndex
         real (kind=double), dimension (:,:), intent(in) :: eps1
         real (kind=double), dimension (:,:), intent(in) :: eps1i
         real (kind=double), dimension (:,:), intent(in) :: elf
         real (kind=double), dimension (:,:), intent(in) :: refractIndex
      end subroutine printData
   end interface

   ! Get command line parameter for the number of lines.
   call getarg(1,buffer)
   read (buffer,*) length

   ! Adjust the length to subtract the header
   length = length - 1

   ! Determine the number of values from the number of lines and the
   !   granularity we desire.
   numValues=length*grain

   ! Determine if this is a spin polarized calculation or not.
   call getarg(2,buffer)
   read (buffer,*) spin

   if (mpiRank == 0) then
      if (spin == 1) then
         open (unit=50, file='fort.50', status="old",form="formatted") ! Eps2
         open (unit=100,file='fort.100',status="new",form="formatted") ! Total
         open (unit=110,file='fort.110',status="new",form="formatted") ! Eps1
         open (unit=120,file='fort.120',status="new",form="formatted") ! ELF 
         open (unit=130,file='fort.130',status="new",form="formatted") ! Rfrct
         open (unit=140,file='fort.140',status="new",form="formatted") ! Eps1i
      else
         open (unit=50, file='fort.51', status="old",form="formatted") ! Eps2
         open (unit=100,file='fort.101',status="new",form="formatted") ! Total
         open (unit=110,file='fort.111',status="new",form="formatted") ! Eps1
         open (unit=120,file='fort.121',status="new",form="formatted") ! ELF 
         open (unit=130,file='fort.131',status="new",form="formatted") ! Rfrct
         open (unit=140,file='fort.141',status="new",form="formatted") ! Eps1i
      endif
   endif

   ! In the future, perhaps this could allocate only the space that is actually
   !   used by each process for eps1 and eps1i. For now though, we only
   !   allocate matrices according to which processes use them.
   allocate (energy(length))
   allocate (fine_energy(numValues))
   allocate (totalEps2(length))
   allocate (eps1(dim3,length))
   allocate (eps1i(dim3,length))
   allocate (eps2(dim3,length))
   allocate (fine_eps2(dim3,numValues))
   if (mpiRank == 0) then
      allocate (elf(dim3,length))
      allocate (refractIndex(dim3,length))
      allocate (totalElf(length))
      allocate (totalRefractIndex(length))
      allocate (totalEps1(length))
      allocate (totalEps1i(length))
   endif

   ! Read the energy and epsilon2 data (output of the optc program).
   call readOPTCData(length, energy, totalEps2, eps2, mpiRank)

   ! Determine the average energy difference for the coarse grain input.
   call getCoarseEnergyDiff(length,energy,coarseEnergyDiff)

   ! Generate a fine grained version of the energy data.
   call makeFineGrainEnergy(length,grain,energy,fine_energy)

   ! Determine the average energy difference for the fine grain energies.
   call getFineEnergyDiff(numValues,grain,fine_energy,fineEnergyDiff)

   ! Generate a fine grained version of the eps2 data.  This will be used to
   !  improve the accuracy of the integration as well as to avoid the poles in
   !  the integral.
   call makeFineEps2(length,grain,energy,fine_energy,eps2,fine_eps2)

   ! Do Kramers Kronig Conversion to get Epsilon1 from Epsilon2, and Epsilon1i
   !   from Epsiton2.
   call kramersKronig(length,grain,numValues,mpiRank,mpiSize,fineEnergyDiff,&
                      &energy,fine_energy,eps1,eps1i,fine_eps2)

   if (mpiRank == 0) then

      ! Calculate the refractive index from epsilon1.
      call getRefractIndex(length,eps1,eps2,refractIndex)
   
      ! Calculate the Energy Loss Function from Epsilon1 and Epsilon2
      call getELF(length,eps1,eps2,elf)
   
      ! Average Epsilon1, Epsilon1i, Epsilon2, and the ELF over 3 dimensions
      call averageFunctions(length,eps1,eps1i,refractIndex,&
                          & totalEps1,totalEps1i,totalEps2,totalElf,&
                          & totalRefractIndex)
   
      ! Print a data file of Epsilon1, Epsilon1i, Epsilon2 and the ELF
      call printData(length,energy,totalEps1,totalEps1i,totalEps2,totalElf,&
                    &totalRefractIndex,eps1,eps1i,elf,refractIndex)
   endif

   ! Deallocate.
   deallocate (energy)
   deallocate (fine_energy)
   deallocate (totalEps2)
   deallocate (eps1)
   deallocate (eps1i)
   deallocate (eps2)
   deallocate (fine_eps2)
   if (mpiRank == 0) then
      deallocate (elf)
      deallocate (refractIndex)
      deallocate (totalElf)
      deallocate (totalRefractIndex)
      deallocate (totalEps1)
      deallocate (totalEps1i)
   endif

   ! End the MPI interface
   call MPI_FINALIZE (mpiErr)

   if (mpiRank == 0) then
      stop 'Kramers-Kronig conversion relation applied.'
   endif

end program OLCAOkkc

subroutine readOPTCData (length, energy, totalEps2, eps2, mpiRank)

   use O_Kinds
   use MPI

   implicit none

   ! Define passed parameters
   integer, intent(in) :: length
   real (kind=double), dimension (:), intent(out) :: energy
   real (kind=double), dimension (:), intent(out) :: totalEps2
   real (kind=double), dimension (:,:), intent(out) :: eps2
   integer, intent(in) :: mpiRank

   ! Define local variables.
   integer :: i
   integer :: mpiErr

   if (mpiRank == 0) then

      ! Read past the header
      read (50,*)

      ! Read the energy and epsilon2 data (output of the optc program).
      do i=1,length
         read (50,*) energy(i),totalEps2(i),eps2(:,i)
      enddo
   endif

   ! Broadcast the data to all processes.
   call MPI_BCAST(energy(1:length),length,MPI_DOUBLE_PRECISION,0,&
         & MPI_COMM_WORLD,mpiErr)

end subroutine readOPTCData

subroutine getCoarseEnergyDiff(length, energy, coarseEnergyDiff)

   use O_Kinds

   implicit none

   integer, intent(in) :: length
   real (kind=double), dimension (:), intent(in) :: energy
   real (kind=double), intent(out) :: coarseEnergyDiff

   integer :: i

   !Determine the average energy difference for the coarse grain input.
   coarseEnergyDiff=0.0_double
   do i=1,length-1
      coarseEnergyDiff=coarseEnergyDiff+(energy(i+1)-energy(i))
   enddo
   coarseEnergyDiff=coarseEnergyDiff/(length-1)

end subroutine getCoarseEnergyDiff

subroutine makeFineGrainEnergy(length,grain,energy,fine_energy)

   use O_Kinds

   implicit none

   integer, intent(in) :: length
   integer, intent(in) :: grain
   real (kind=double), dimension (:), intent(in)  :: energy
   real (kind=double), dimension (:), intent(out) :: fine_energy

   integer :: i,j,k

   ! Generate a fine grained version of the energy data.
   do i=1,length-1
      do j=1,grain
         k=int((i-1)*grain+j)
         fine_energy(k) = energy(i)+(energy(i+1)-energy(i))*(j-1)/grain
      enddo
   enddo

end subroutine makeFineGrainEnergy

subroutine getFineEnergyDiff(numValues,grain,fine_energy,fineEnergyDiff)

   use O_Kinds

   implicit none

   integer, intent(in) :: numValues
   integer, intent(in) :: grain
   real (kind=double), dimension (:), intent(in) :: fine_energy
   real (kind=double), intent(out) :: fineEnergyDiff

   integer :: i

   fineEnergyDiff=0.0_double
   do i=1,numValues-grain-1
      fineEnergyDiff=fineEnergyDiff+fine_energy(i+1)-fine_energy(i)
   enddo
   fineEnergyDiff=fineEnergyDiff/(numValues-grain-1)

end subroutine getFineEnergyDiff

subroutine makeFineEps2(length,grain,energy,fine_energy,eps2,fine_eps2)

   use O_Kinds

   implicit none

   integer, intent(in) :: length
   integer, intent(in) :: grain
   real (kind=double), dimension (:), intent(in) :: energy
   real (kind=double), dimension (:), intent(in) :: fine_energy
   real (kind=double), dimension (:,:), intent(in) :: eps2
   real (kind=double), dimension (:,:), intent(out) :: fine_eps2

   integer :: i,j,k
   real (kind=double), allocatable, dimension(:,:) :: diffRatio

   allocate(diffRatio(dim3,length))

   do i=1,length-1
      diffRatio(:,i)=(eps2(:,i+1)-eps2(:,i))/(energy(i+1)-energy(i))
      do j=1,grain
         k=int((i-1)*grain+j)
         fine_eps2(:,k) = eps2(:,i)+diffRatio(:,i)*(fine_energy(k)-energy(i))
      enddo
   enddo

   deallocate(diffRatio)

end subroutine makeFineEps2

subroutine kramersKronig(length,grain,numValues,mpiRank,mpiSize,&
      & fineEnergyDiff,energy,fine_energy,eps1,eps1i,fine_eps2)

   use O_Kinds
   use O_Constants
   use MPI

   implicit none

   integer, intent(in) :: length
   integer, intent(in) :: grain
   integer, intent(in) :: numValues
   integer, intent(in) :: mpiRank
   integer, intent(in) :: mpiSize
   real (kind=double), intent(in) :: fineEnergyDiff
   real (kind=double), dimension (:), intent(in) :: energy
   real (kind=double), dimension (:), intent(in) :: fine_energy
   real (kind=double), dimension (:,:), intent(out) :: eps1
   real (kind=double), dimension (:,:), intent(out) :: eps1i
   real (kind=double), dimension (:,:), intent(in) :: fine_eps2

   ! Define local variables.
   integer :: i,j
   integer :: mpiErr
   integer :: rangeInit
   integer :: rangeFin
   integer :: workPerProcess
   integer :: remainder
   integer :: currRank
   integer, allocatable, dimension (:) :: rangeInitArray
   integer, allocatable, dimension (:) :: rangeFinArray
   real (kind=double) :: multFactor
   real (kind=double) :: fine_e,original_e
   real (kind=double), allocatable, dimension (:) :: totalSum,evenSum,oddSum
   real (kind=double), allocatable, dimension (:) :: totalSumi,evenSumi,oddSumi
   real (kind=double), allocatable, dimension (:,:) :: integrand,integrandi

   multFactor = 2.0_double/3.0_double/pi

   allocate (totalSum   (dim3))
   allocate (evenSum    (dim3))
   allocate (oddSum     (dim3))
   allocate (integrand  (dim3,numValues))
   allocate (totalSumi  (dim3))
   allocate (evenSumi   (dim3))
   allocate (oddSumi    (dim3))
   allocate (integrandi (dim3,numValues))
   allocate (rangeInitArray(mpiSize))
   allocate (rangeFinArray (mpiSize))

   ! Each process will compute a subset of the values for eps1 and eps1i. This
   !   bit that follows is the load balancing. Interestingly, processs zero
   !   needs to compute these values for all processes and store then in an
   !   array so that it can execute an MPI_GATHERV call to collect variable
   !   quantities of data from each process. (The quantities are variable
   !   because of load balancing).

   ! First we obtain the work per process in a crude way and record the
   !   default ranges and what the remainder of the division of work was.
   workPerProcess = int(length/mpiSize)
   remainder = mod(length,mpiSize)
   rangeInit = workPerProcess*mpiRank + 1
   rangeFin  = workPerProcess*(mpiRank+1)

   ! Now the idea is to spread one unit of extra work across a number of
   !   processes that is equal to the size of the remainder. (I.e. if the
   !   remainder is five, then five processors will have one unit of work
   !   extra to do.) The trick is that we want to make a direct assignment
   !   because each process decides its indices independently (without any
   !   communication). We can't take a serial approach and "shift" the ranges
   !   down over and over. (I.e. We can't find one new range and then adjust
   !   the next based on the previous because each process doesn't know the
   !   "previous" and we don't want to have to make them communicate that.)
   ! So, consider a job with 1000 cores (mpiSize=1000) and a remainder of 5.
   ! The first if statement just increments the final index for process 995 =
   !   (1000-5) by one while the initial index stays the same. Its new final
   !   index is (300*(995+1)) + 1 = 298801. (With a workPerProcess of 300.)
   if (mpiRank == (mpiSize - remainder)) then
      rangeFin = rangeFin + 1
   endif
   ! The next if statement will be entered for processes 996, 997, 998, and 999.
   !   (Note that there is no process with mpiRank=1000 because the first rank
   !   is 0.) Once inside, we recompute the initial and final array indices.
   !   We must take care that the initial index shift by an amount to
   ! For the example, let process 996 already have rangeInit=300*996+1=298801
   !   (with a workPerProcess of 300). This becomes 298801 + 5 - (1000-996) = 
   !   298802. The finalRange will increment by two, one due to the shift from
   !   the process before it (995) and one for the extra work this process will
   !   be doing.
   ! One more example. Consider process 998. Its original rangeInit is 299401
   !   and its original rangeFin is 299700. The new values are:
   !   rangeInit = 299401 + (5 - (1000-998))     = 299404
   !   rangeFin  = 299700 + (5 - (1000-(998+1))) = 299704
   ! As the rank increases, the amount of the shift in rangeInit increases.
   if (mpiRank > (mpiSize-remainder)) then
      rangeInit = rangeInit + (remainder - (mpiSize - mpiRank))
      rangeFin  = rangeFin  + (remainder - (mpiSize - (mpiRank+1)))
   endif

   ! Here we basically do the same thing over again, but just for the benefit
   !   of process zero.
   if (mpiRank == 0) then
      do i = 1, mpiSize
         currRank = i-1
         rangeInitArray(i) = workPerProcess*currRank + 1
         rangeFinArray(i)  = workPerProcess*currRank + 1
         if (currRank == (mpiSize - remainder)) then
            rangeFinArray(i) = rangeFinArray(i) + 1
         endif
         if (currRank > (mpiSize-remainder)) then
            rangeInitArray(i) = rangeInitArray(i) + &
                  & (remainder - (mpiSize - currRank))
            rangeFinArray(i)  = rangeFinArray(i)  + &
                  & (remainder - (mpiSize - (currRank+1)))
         endif
      enddo
   else
      rangeInitArray(:) = 0
      rangeFinArray(:)  = 0
   endif

   ! Start the loop over appropriate range for each process.
   do i = rangeInit, rangeFin
      totalSum(:)  = 0.0_double
      oddSum(:)    = 0.0_double
      evenSum(:)   = 0.0_double
      totalSumi(:) = 0.0_double
      oddSumi(:)   = 0.0_double
      evenSumi(:)  = 0.0_double

      ! Do the first point explicitly to save computation in the loop
      fine_e=fine_energy(1)
      original_e=energy(1)
      if (abs(fine_e - original_e) .gt. 0.00001_double) then
         integrand(:,1)=fine_eps2(:,1)*fine_e/ &
               & (original_e + fine_e)/(fine_e - original_e)
         totalSum(:)=totalSum(:)+integrand(:,1)
         integrandi(:,1)=fine_eps2(:,1)*fine_e/ &
               & (fine_e*fine_e + original_e*original_e)
         totalSumi(:)=totalSumi(:)+integrandi(:,1)
      endif

      ! Loop through the middle part
      do j=2,numValues-grain-1
         !Exclude all points from energy(i) in the summation to avoid
         !  poles in the integration.
         fine_e=fine_energy(j)
         original_e=energy(i)
         if (abs(fine_e - original_e) .gt. 0.00001_double) then
            integrand(:,j)=fine_eps2(:,j)*fine_e / &
                  & (original_e + fine_e)/(fine_e - original_e)
            integrandi(:,j)=fine_eps2(:,j)*fine_e / &
                  & (fine_e*fine_e + original_e*original_e)
            if ((j/2)*2 .eq. j) then
               evenSum(:)=evenSum(:)+integrand(:,j)
               evenSumi(:)=evenSumi(:)+integrandi(:,j)
            else
               oddSum(:)=oddSum(:)+integrand(:,j)
               oddSumi(:)=oddSumi(:)+integrandi(:,j)
            endif
            totalSum(:)=totalSum(:)+integrand(:,j)
            totalSumi(:)=totalSumi(:)+integrandi(:,j)
         endif
      enddo

      ! Do the last point explicitly to save computation in the loop
      fine_e=fine_energy(numValues-grain)
      original_e=energy(length)
      if (abs(fine_e - original_e) .gt. 0.00001_double) then
         integrand(:,numValues-grain)=fine_eps2(:,numValues-grain)*fine_e/ &
               & (original_e + fine_e)/(fine_e - original_e)
         totalSum(:)=totalSum(:)+integrand(:,numValues-grain)
         integrandi(:,numValues-grain)=fine_eps2(:,numValues-grain)*fine_e/ &
               & (fine_e*fine_e + original_e*original_e)
         totalSumi(:)=totalSumi(:)+integrandi(:,numValues-grain)
      endif

      eps1(:,i)=fineEnergyDiff*(integrand(:,int(numValues-grain)) + &
            & integrand(:,1) + 2.0_double * evenSum(:) + 4.0_double * &
            & oddSum(:)) * multFactor + 1.0_double
      eps1i(:,i)=fineEnergyDiff*(integrandi(:,int(numValues-grain)) + &
            & integrandi(:,1) + 2.0_double * evenSumi(:) + 4.0_double * &
            & oddSumi(:)) * multFactor + 1.0_double
   enddo

   ! Now we need to gather the results into process zero so that it can do
   !   a very few final computations (refractIndex and ELF) and then print
   !   out the total results.
   ! Presently, I'm going to do this in a rather inefficient and probably
   !   stupid way. My goal is (first) just to get it to work in parallel in the
   !   limited time that I have available. Basically, this will produce a
   !   complete message storm where each message is very small. Ideally this
   !   should be reworked to function with multidimensional matrices and
   !   probably a MPI_GATHERV type of call, but that is not going to happen
   !   at the moment. (A good resource is: http://stackoverflow.com/questions/17508647/sending-2d-arrays-in-fortran-with-mpi-gather.)

   if (mpiRank /= 0) then

      do i = rangeInit, rangeFin

         ! Send the eps1 results to process 0. (Argument 5 = tag = 1)
         call MPI_SEND(eps1(1:3,i),3,MPI_DOUBLE_PRECISION,0,1,&
               & MPI_COMM_WORLD,mpierr)

         ! Send the eps1i results to process 0. (Argument 5 = tag = 2)
         call MPI_SEND(eps1i(1:3,i),3,MPI_DOUBLE_PRECISION,0,2,&
               & MPI_COMM_WORLD,mpierr)
      enddo
   else

      ! Collect results from each process.
      do i = 1, mpiSize-1

         ! Each process has to send a set of results.
         do j = rangeInitArray(i), rangeFinArray(i)

            ! Receive the epsilon1 data. (Argument 5 = tag = 1)
            call MPI_RECV(eps1(1:3,j),3,MPI_DOUBLE_PRECISION,i,1,&
                  & MPI_WORLD_COMM,MPI_STATUS_IGNORE,mpiErr)

            ! Receive the imaginary epsilon1 data. (Argument 5 = tag = 2)
            call MPI_RECV(eps1i(1:3,j),3,MPI_DOUBLE_PRECISION,i,2,&
                  & MPI_WORLD_COMM,MPI_STATUS_IGNORE,mpiErr)
         enddo
      enddo
   endif

   deallocate (totalSum)
   deallocate (evenSum)
   deallocate (oddSum)
   deallocate (integrand)
   deallocate (totalSumi)
   deallocate (evenSumi)
   deallocate (oddSumi)
   deallocate (integrandi)
   deallocate (rangeInitArray)
   deallocate (rangeFinArray)

end subroutine kramersKronig

subroutine getRefractIndex(length,eps1,eps2,refractIndex)

   use O_Kinds

   implicit none

   integer, intent(in) :: length
   real (kind=double), dimension (:,:), intent(in)  :: eps1
   real (kind=double), dimension (:,:), intent(in)  :: eps2
   real (kind=double), dimension (:,:), intent(out) :: refractIndex

   integer :: i,j

   do i = 1,length
      do j = 1,3

         ! ~eps = eps1 + i eps2
         ! eps1 = n^2 - k^2  (n=refractive index, k=absorption coeff)
         ! eps2 = 2nk
         ! n = sqrt[(sqrt(eps1^2 + eps2^2)+eps1)/2]
         refractIndex(j,i) = sqrt((sqrt(eps1(j,i)*eps1(j,i) + &
               & eps2(j,i)*eps2(j,i))+eps1(j,i))/2.0_double)
      enddo
   enddo
end subroutine getRefractIndex

subroutine getELF(length,eps1,eps2,elf)

   use O_Kinds

   implicit none

   integer, intent(in) :: length
   real (kind=double), dimension (:,:), intent(in)  :: eps1
   real (kind=double), dimension (:,:), intent(in)  :: eps2
   real (kind=double), dimension (:,:), intent(out) :: elf

   integer :: i

   ! elf = eps2/(eps1^2 +eps2^2)
   do i = 1,length
      elf(:,i) = eps2(:,i)/(eps1(:,i)**2 + eps2(:,i)**2)
   enddo


end subroutine getELF

subroutine averageFunctions(length,eps1,eps1i,refractIndex,totalEps1,&
                          & totalEps1i,totalEps2,totalElf,totalRefractIndex)

   use O_Kinds

   implicit none

   integer, intent(in) :: length
   real (kind=double), dimension (:,:), intent(in)  :: eps1
   real (kind=double), dimension (:,:), intent(in)  :: eps1i
   real (kind=double), dimension (:,:), intent(in)  :: refractIndex
   real (kind=double), dimension (:),   intent(out) :: totalEps1
   real (kind=double), dimension (:),   intent(out) :: totalEps1i
   real (kind=double), dimension (:),   intent(inout) :: totalEps2
   real (kind=double), dimension (:),   intent(out) :: totalElf
   real (kind=double), dimension (:),   intent(out) :: totalRefractIndex

   integer :: i

   do i = 1,length
      totalEps1(i)         = sum(eps1(:,i))/3.0_double
      totalEps1i(i)        = sum(eps1i(:,i))/3.0_double
      totalElf(i)          = totalEps2(i)/(totalEps1(i)**2+totalEps2(i)**2)
      totalRefractIndex(i) = sum(refractIndex(:,i))/3.0_double
   enddo

end subroutine averageFunctions

subroutine printData(length,energy,totalEps1,totalEps1i,totalEps2,totalElf,&
      & totalRefractIndex,eps1,eps1i,elf,refractIndex)

   use O_Kinds

   implicit none

   integer, intent(in) :: length
   real (kind=double), dimension (:), intent(in)   :: energy
   real (kind=double), dimension (:), intent(in)   :: totalEps1
   real (kind=double), dimension (:), intent(in)   :: totalEps1i
   real (kind=double), dimension (:), intent(in)   :: totalEps2
   real (kind=double), dimension (:), intent(in)   :: totalElf
   real (kind=double), dimension (:), intent(in)   :: totalRefractIndex
   real (kind=double), dimension (:,:), intent(in) :: eps1
   real (kind=double), dimension (:,:), intent(in) :: eps1i
   real (kind=double), dimension (:,:), intent(in) :: elf
   real (kind=double), dimension (:,:), intent(in) :: refractIndex

   integer :: i

   write (100,*) "Energy   Epsilon1   Epsilon2   ELF"
   do i = 1,length
      write (100,701) energy(i),totalEps1(i),totalEps2(i),totalELF(i)
   enddo

   write (110,*) "Energy   totalEps1   xEps1   yEps1   zEps1"
   do i = 1,length
      write (110,700) energy(i),totalEps1(i),eps1(:,i)
   enddo

   write (120,*) "Energy   totalELF   xELF   yELF   zELF"
   do i = 1,length
      write (120,700) energy(i),totalELF(i),elf(:,i)
   enddo

   write (130,*) "Energy   totaln      xn      yn      zn"
   do i = 1,length
      write (130,700) energy(i),totalRefractIndex(i),refractIndex(:,i)
   enddo

   write (140,*) "Energy   totalEps1i   xEps1i   yEps1i   zEps1i"
   do i = 1,length
      write (140,700) energy(i),totalEps1i(i),eps1i(:,i)
   enddo

   700 format (5e15.7)
   701 format (4e15.7)

end subroutine printData
