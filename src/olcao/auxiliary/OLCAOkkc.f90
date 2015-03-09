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

   ! Define the subroutine interfaces for proper passing of allocatable arrays.
   interface
      subroutine readOPTCData (length, energy, totalEps2, eps2)
         use O_Kinds
         use O_Constants
         integer :: length
         real (kind=double), dimension (:) :: energy
         real (kind=double), dimension (:) :: totalEps2
         real (kind=double), dimension (:,:) :: eps2
      end subroutine readOPTCData
      subroutine getCoarseEnergyDiff(length, energy, coarseEnergyDiff)
         use O_Kinds
         use O_Constants
         integer :: length
         real (kind=double), dimension (:) :: energy
         real (kind=double) :: coarseEnergyDiff
      end subroutine getCoarseEnergyDiff
      subroutine makeFineGrainEnergy(length,grain,energy,fine_energy)
         use O_Kinds
         use O_Constants
         integer :: length, grain
         real (kind=double), dimension (:) :: energy, fine_energy
      end subroutine makeFineGrainEnergy
      subroutine getFineEnergyDiff(numValues,grain,fine_energy,fineEnergyDiff)
         use O_Kinds
         use O_Constants
         integer :: numValues, grain
         real (kind=double), dimension (:) :: fine_energy
         real (kind=double) :: fineEnergyDiff
      end subroutine getFineEnergyDiff
      subroutine makeFineEps2(length,grain,energy,fine_energy,eps2,&
                             &fine_eps2)
         use O_Kinds
         use O_Constants
         integer :: length,grain
         real (kind=double), dimension (:) :: energy,fine_energy
         real (kind=double), dimension (:,:) :: eps2,fine_eps2
      end subroutine makeFineEps2
      subroutine kramersKronig(length,grain,numValues,fineEnergyDiff,&
                              &energy,fine_energy,eps1,eps1i,fine_eps2)
         use O_Kinds
         use O_Constants
         integer :: length,grain,numValues
         real (kind=double) :: fineEnergyDiff
         real (kind=double), dimension (:) :: energy,fine_energy
         real (kind=double), dimension (:,:) :: eps1,eps1i,fine_eps2
      end subroutine kramersKronig
      subroutine getRefractIndex(length,eps1,eps2,refractIndex)
         use O_Kinds
         use O_Constants
         integer :: length
         real (kind=double), dimension (:,:) :: eps1,eps2,refractIndex
      end subroutine getRefractIndex
      subroutine getELF(length,eps1,eps2,elf)
         use O_Kinds
         use O_Constants
         integer :: length
         real (kind=double), dimension (:,:) :: eps1,eps2,elf
      end subroutine getELF
      subroutine averageFunctions(length,eps1,eps1i,refractIndex,&
                                 &totalEps1,totalEps1i,totalEps2,totalElf,&
                                 &totalRefractIndex)
         use O_Kinds
         use O_Constants
         integer :: length
         real (kind=double), dimension (:,:) :: eps1,eps1i
         real (kind=double), dimension (:,:) :: refractIndex
         real (kind=double), dimension (:)   :: totalEps1,totalEps1i
         real (kind=double), dimension (:)   :: totalEps2,totalElf
         real (kind=double), dimension (:)   :: totalRefractIndex
      end subroutine averageFunctions
      subroutine printData(length,energy,totalEps1,totalEps1i,totalEps2,&
                         & totalElf,totalRefractIndex,eps1,eps1i,elf,&
                         & refractIndex)
         use O_Kinds
         use O_Constants
         integer :: length
         real (kind=double), dimension (:)   :: energy
         real (kind=double), dimension (:)   :: totalEps1,totalEps1i,totalEps2
         real (kind=double), dimension (:)   :: totalElf,totalRefractIndex
         real (kind=double), dimension (:,:) :: eps1,eps1i,elf,refractIndex
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

   if (spin == 1) then
      open (unit=50, file='fort.50', status="old",form="formatted") ! Eps2
      open (unit=100,file='fort.100',status="new",form="formatted") ! Total
      open (unit=110,file='fort.110',status="new",form="formatted") ! Eps1
      open (unit=120,file='fort.120',status="new",form="formatted") ! ELF 
      open (unit=130,file='fort.130',status="new",form="formatted") ! Refract
      open (unit=140,file='fort.140',status="new",form="formatted") ! Eps1i
   else
      open (unit=50, file='fort.51', status="old",form="formatted") ! Eps2
      open (unit=100,file='fort.101',status="new",form="formatted") ! Total
      open (unit=110,file='fort.111',status="new",form="formatted") ! Eps1
      open (unit=120,file='fort.121',status="new",form="formatted") ! ELF 
      open (unit=130,file='fort.131',status="new",form="formatted") ! Refract
      open (unit=140,file='fort.141',status="new",form="formatted") ! Eps1i
   endif

   allocate (energy(length))
   allocate (fine_energy(numValues))
   allocate (totalEps1(length))
   allocate (totalEps1i(length))
   allocate (totalEps2(length))
   allocate (totalElf(length))
   allocate (totalRefractIndex(length))
   allocate (eps1(dim3,length))
   allocate (eps1i(dim3,length))
   allocate (eps2(dim3,length))
   allocate (fine_eps2(dim3,numValues))
   allocate (elf(dim3,length))
   allocate (refractIndex(dim3,length))

   ! Read the energy and epsilon2 data (output of the optc program).
   call readOPTCData(length, energy, totalEps2, eps2)

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
   call kramersKronig(length,grain,numValues,fineEnergyDiff,&
                      &energy,fine_energy,eps1,eps1i,fine_eps2)

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

   deallocate (energy)
   deallocate (fine_energy)
   deallocate (totalEps1)
   deallocate (totalEps1i)
   deallocate (totalEps2)
   deallocate (totalElf)
   deallocate (totalRefractIndex)
   deallocate (eps1)
   deallocate (eps1i)
   deallocate (eps2)
   deallocate (fine_eps2)
   deallocate (elf)
   deallocate (refractIndex)

   stop 'Kramers-Kronig conversion relation applied.'

end program OLCAOkkc

subroutine readOPTCData (length, energy, totalEps2, eps2)

   use O_Kinds
   use O_Constants

   implicit none

   integer :: length
   real (kind=double), dimension (:) :: energy
   real (kind=double), dimension (:) :: totalEps2
   real (kind=double), dimension (:,:) :: eps2

   integer :: i

   ! Read past the header
   read (50,*)

   ! Read the energy and epsilon2 data (output of the optc program).
   do i=1,length
      read (50,*) energy(i),totalEps2(i),eps2(:,i)
   enddo

end subroutine readOPTCData

subroutine getCoarseEnergyDiff(length, energy, coarseEnergyDiff)

   use O_Kinds
   use O_Constants

   implicit none

   integer :: length
   real (kind=double), dimension (:) :: energy
   real (kind=double) :: coarseEnergyDiff

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
   use O_Constants

   implicit none

   integer :: length, grain
   real (kind=double), dimension (:) :: energy, fine_energy

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
   use O_Constants

   implicit none

   integer :: numValues, grain
   real (kind=double), dimension (:) :: fine_energy
   real (kind=double) :: fineEnergyDiff

   integer :: i

   fineEnergyDiff=0.0_double
   do i=1,numValues-grain-1
      fineEnergyDiff=fineEnergyDiff+fine_energy(i+1)-fine_energy(i)
   enddo
   fineEnergyDiff=fineEnergyDiff/(numValues-grain-1)

end subroutine getFineEnergyDiff

subroutine makeFineEps2(length,grain,energy,fine_energy,eps2,fine_eps2)

   use O_Kinds
   use O_Constants

   implicit none

   integer :: length,grain
   real (kind=double), dimension (:) :: energy,fine_energy
   real (kind=double), dimension (:,:) :: eps2,fine_eps2

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

subroutine kramersKronig(length,grain,numValues,fineEnergyDiff,&
                        &energy,fine_energy,eps1,eps1i,fine_eps2)

   use O_Constants
   use O_Kinds

   implicit none

   integer :: length,grain,numValues
   real (kind=double) :: fineEnergyDiff
   real (kind=double), dimension (:) :: energy,fine_energy
   real (kind=double), dimension (:,:) :: eps1,eps1i,fine_eps2

   integer :: i,j
   real (kind=double) :: fine_e,original_e
   real (kind=double), allocatable, dimension (:) :: totalSum,evenSum,oddSum
   real (kind=double), allocatable, dimension (:) :: totalSumi,evenSumi,oddSumi
   real (kind=double), allocatable, dimension (:,:) :: integrand,integrandi

   real (kind=double) :: multFactor

   multFactor = 2.0_double/3.0_double/pi

   allocate (totalSum   (dim3))
   allocate (evenSum    (dim3))
   allocate (oddSum     (dim3))
   allocate (integrand  (dim3,numValues))
   allocate (totalSumi  (dim3))
   allocate (evenSumi   (dim3))
   allocate (oddSumi    (dim3))
   allocate (integrandi (dim3,numValues))

   do i=1,length
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
      if (abs(fine_e - original_e) .gt. 0.00001) then
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

   deallocate (totalSum)
   deallocate (evenSum)
   deallocate (oddSum)
   deallocate (integrand)
   deallocate (totalSumi)
   deallocate (evenSumi)
   deallocate (oddSumi)
   deallocate (integrandi)

end subroutine kramersKronig

subroutine getRefractIndex(length,eps1,eps2,refractIndex)

   use O_Kinds
   use O_Constants

   implicit none

   integer :: length
   real (kind=double), dimension (:,:) :: eps1,eps2,refractIndex

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
   use O_Constants

   implicit none

   integer :: length
   real (kind=double), dimension (:,:) :: eps1,eps2,elf

   integer :: i

   ! elf = eps2/(eps1^2 +eps2^2)
   do i = 1,length
      elf(:,i) = eps2(:,i)/(eps1(:,i)**2 + eps2(:,i)**2)
   enddo


end subroutine getELF

subroutine averageFunctions(length,eps1,eps1i,refractIndex,totalEps1,&
                          & totalEps1i,totalEps2,totalElf,totalRefractIndex)

   use O_Kinds
   use O_Constants

   implicit none

   integer :: length
   real (kind=double), dimension (:,:) :: eps1,eps1i,refractIndex
   real (kind=double), dimension (:) :: totalEps1,totalEps1i,totalEps2,&
         & totalElf,totalRefractIndex

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
   use O_Constants

   implicit none

   integer :: length
   real (kind=double), dimension (:)   :: energy,totalEps1,totalEps1i
   real (kind=double), dimension (:)   :: totalEps2,totalElf
   real (kind=double), dimension (:)   :: totalRefractIndex
   real (kind=double), dimension (:,:) :: eps1,eps1i,elf,refractIndex

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
