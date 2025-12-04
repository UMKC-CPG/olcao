module O_MTOP

   ! Import necessary modules.
   use O_Kinds
   use O_Constants
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ! Begin list of module subroutines.!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

#ifndef GAMMA
subroutine computeMTOPPolarization(inSCF,xyzP)

      ! Import necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Constants,   only: pi, eCharge, smallThresh, bohrRad
   use O_Potential,   only: spin
   use O_Populate,    only: electronPopulation
   use O_KPoints,     only: numKPoints, numAxialKPoints, mtopKPMap, kPoints
   use O_Lattice,     only: realVectors, realCellVolume
   use O_AtomicSites, only: valeDim, computeIonicMoment, xyzIonMoment, &
         & abcIonMoment
   use O_Input,       only: numStates
   use O_MatrixSubs,  only: fullMatrixElementMult
   use O_SecularEquation, only: valeVale, valeValeKO, readDataSCF, readDataPSCF
   use O_BLASZGERC

   ! Make sure that no funny variables are used.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: inSCF
   real(kind=double), dimension(3,2), intent(out) :: xyzP

   ! Define local variables.
   integer :: h,i,j,k,l
integer :: n,m
   integer :: maxOccupiedState
   integer :: currNumLines,numSteps
   integer :: energyLevelCounter
   integer :: kcur, knxt, axis
   integer :: kPointCount
   integer :: maxNumStrings
   integer, dimension(3) :: numStrings
   real(kind=double), dimension(3,2) :: xyzP_A
   real(kind=double), dimension(3,2) :: xyzP_B
   real(kind=double), dimension(3,2) :: xyzP_C
   real(kind=double), dimension(3,2) :: averagePhase_A
   real(kind=double), dimension(3,2) :: averagePhase_B
   real(kind=double), dimension(3,2) :: averagePhase_C
   real(kind=double), dimension (3,2) :: psi
   real(kind=double), allocatable, dimension (:) :: currentPopulation
   real (kind=double), allocatable, dimension (:,:,:) :: &
         & structuredElectronPopulation
   complex(kind=double), allocatable, dimension(:) :: phiString_B
   complex(kind=double), allocatable, dimension(:) :: phiString_C
   complex(kind=double), allocatable, dimension(:,:,:) :: valeValePsi
   complex(kind=double), allocatable, dimension(:,:,:) :: CKnxt
   complex(kind=double), allocatable, dimension(:,:,:) :: stateStateMat
   complex(kind=double), allocatable, dimension(:,:,:) :: prodM_A
   real(kind=double), allocatable, dimension(:,:) :: stringPhaseSet_A
   real(kind=double), allocatable, dimension(:,:) :: stringPhaseSet_B
   real(kind=double), allocatable, dimension(:,:) :: stringPhaseSet_C
complex(kind=double), allocatable, dimension(:,:) :: unitary

   complex(kind=double), dimension(2,2) :: ssm
   complex(kind=double), dimension(2,2) :: pma
   complex(kind=double) :: psb, psc
   real(kind=double) :: resultA
   real(kind=double) :: resultB
   real(kind=double) :: resultC

   ! Record the start of the calculation
   call timeStampStart(33)

   ! Compute the ionic contribution to the polarization.
   call computeIonicMoment(0)

   ! Demonstrate correctness of mathematical operations.
   write(20,*) "INITIALIZE"
   ssm = reshape((/(1.0,4.0), (2.0,3.0), (3.0,2.0), (4.0,1.0)/), shape(ssm))
   write(20,*) "ssm=",ssm(:,:)
   pma = reshape((/(1.0,0.0), (0.0,0.0), (0.0,0.0), (1.0,0.0)/), shape(pma))
   write(20,*) "pma=",pma(:,:)
   psb = (1.0, 0.0)
   write(20,*) "psb=",psb
   psc = (0.0, 0.0)
   write(20,*) "psc=",psc

   write(20,*) "STEP 1"
   pma = matmul(pma,ssm)
   psb = psb*matrixDet(ssm)
   psc = psc + log(matrixDet(ssm))
   write(20,*) "pma=",pma(:,:)
   write(20,*) "psb=",psb
   write(20,*) "psc=",psc
   write(20,*)

   write(20,*) "NEXT SSM"
   ssm = reshape((/(2.0,5.0), (1.0,1.0), (7.0,3.0), (8.0,4.0)/), shape(ssm))
   write(20,*) "ssm=",ssm(:,:)

   write(20,*) "STEP 2"
   pma = matmul(pma,ssm)
   psb = psb*matrixDet(ssm)
   psc = psc + log(matrixDet(ssm))
   write(20,*) "pma=",pma(:,:)
   write(20,*) "psb=",psb
   write(20,*) "psc=",psc

   write(20,*) "RESULT"
   resultA = -aimag(log(matrixDet(pma)))
   resultB = -aimag(log(psb))
   resultC = -aimag(psc)
   write(20,*) "A ", resultA
   write(20,*) "B ", resultB
   write(20,*) "C ", resultC

   write(20,*) "RESULT in -pi to pi"
   do while (resultA < -pi)
      resultA = resultA + 2.0d0*pi
   enddo
   do while (resultA > pi)
      resultA = resultA - 2.0d0*pi
   enddo
   do while (resultB < -pi)
      resultB = resultB + 2.0d0*pi
   enddo
   do while (resultB > pi)
      resultB = resultB - 2.0d0*pi
   enddo
   do while (resultC < -pi)
      resultC = resultC + 2.0d0*pi
   enddo
   do while (resultC > pi)
      resultC = resultC - 2.0d0*pi
   enddo
   write(20,*) "A ", resultA
   write(20,*) "B ", resultB
   write(20,*) "C ", resultC

   ! Compute the number of strings for each axis and also the maximum number
   !   of strings of all axes.
   numStrings(1) = numAxialKPoints(2)*numAxialKPoints(3)
   numStrings(2) = numAxialKPoints(1)*numAxialKPoints(3)
   numStrings(3) = numAxialKPoints(1)*numAxialKPoints(2)
   maxNumStrings = maxval(numStrings(:))


   ! allocate the Eigenvectors for sequences C
   ! CK and CnextK are eigenvectors at k and k+Δk
   allocate (CKnxt(valeDim,numStates,spin))
   allocate (valeValePsi(valeDim,valeDim,spin))
   allocate (valeValeKO(valeDim,valeDim,3))
write (20,*) "inSCF = ",inSCF
   if (inSCF == 0) then
      allocate (valeVale(valeDim,numStates,spin))
   endif
   allocate (structuredElectronPopulation (numStates,numKPoints,spin))
   allocate (phiString_B(spin))
   allocate (phiString_C(spin))
   allocate (currentPopulation(spin))
   allocate (stringPhaseSet_A(maxNumStrings,spin))
   allocate (stringPhaseSet_B(maxNumStrings,spin))
   allocate (stringPhaseSet_C(maxNumStrings,spin))

   valeVale(:,:,:) = cmplx(0.0_double,0.0_double,double)
   valeValeKO(:,:,:) = cmplx(0.0_double,0.0_double,double)
 
   psi(:,:) = 0.0d0
   xyzP(:,:) = 0.0d0
   xyzP_A(:,:) = 0.0d0
   xyzP_B(:,:) = 0.0d0
   xyzP_C(:,:) = 0.0d0

   ! Fill a matrix of electron populations from the electron population that
   !   was computed in populateLevels.  Note that electronPopulation is a one
   !   dimensional array that has some order, but is not sorted in the way
   !   that the energy eigen values were sorted.  Please read the comments in
   !   the populateLevels subroutine to understand the order.
   !   (You can also probably get it from the loop order here ;)
   energyLevelCounter=0
   do i = 1, numKPoints
      do j = 1, spin
         do k = 1, numStates
            energyLevelCounter = energyLevelCounter + 1
            structuredElectronPopulation (k,i,j) = &
                  & electronPopulation(energyLevelCounter)
         enddo
      enddo
   enddo

   ! Compute the index number of the highest occupied state for kpoint 1.
   maxOccupiedState = 0
   do i = 1, numStates
      if (sum(structuredElectronPopulation(i,1,:)) < smallThresh) then
         maxOccupiedState = i-1
         exit
      endif
   enddo

   ! Check that all other kPoints maintain the same highest occupied state.
   do i = 1, numKPoints
      if (sum(structuredElectronPopulation(maxOccupiedState,i,:)) < &
            & smallThresh) then
         stop "State index < maxOccupiedState is not occupied enough."
      endif
      if (sum(structuredElectronPopulation(maxOccupiedState+1,i,:)) > &
            & smallThresh) then
         stop "State index > maxOccupiedState occupied too much."
      endif
   enddo

   ! Now that we know the max number of occupied states, we can allocate
   !   space to hold the results as we compute them.
   allocate (stateStateMat(maxOccupiedState,maxOccupiedState,spin))
   allocate (prodM_A(maxOccupiedState,maxOccupiedState,spin))
allocate (unitary(maxOccupiedState,maxOccupiedState))
   write(20,*) "Expecting ", sum(numStrings(:)), " strings in three groups."
write(20,*) "maxOccState = ", maxOccupiedState

   do axis = 1, 3
      kPointCount = 0
      numSteps = numAxialKPoints(axis)
      currNumLines = numStrings(axis)
      stringPhaseSet_A(:,:) = 0.0_double
      stringPhaseSet_B(:,:) = 0.0_double
      stringPhaseSet_C(:,:) = 0.0_double

      do h = 1, currNumLines
         stateStateMat(:,:,:) = cmplx(0.0_double,0.0_double,double)
         phiString_B(:) = cmplx(1.0_double,0.0_double,double)
         phiString_C(:) = cmplx(0.0_double,0.0_double,double)
         prodM_A(:,:,:) = cmplx(0.0_double,0.0_double,double)
         do i = 1, maxOccupiedState ! Set to identity matrix.
            prodM_A(i,i,:) = cmplx(1.0_double,0.0_double,double)
         enddo

         do i = 1, numSteps
            kPointCount = kPointCount + 1
            kcur = mtopKPMap(axis,kPointCount)
            if (i < numSteps) then
               knxt = mtopKPMap(axis,kPointCount+1)
            else
               knxt = mtopKPMap(axis,kPointCount-numSteps+1)
            endif

write(20,*) "kcur knxt", kcur, knxt
write(20,*) "kcur", kPoints(:,kcur)
write(20,*) "knxt", kPoints(:,knxt)
write(20,*) "kdiff", kPoints(:,knxt) - kPoints(:,kcur)

            ! Read next kpoint waveFn coefficients and appropriate KOverlap.
            if (i < numSteps) then ! Code 3,4,5
               do j = 1, spin
                  if (inSCF == 1) then
                     call readDataSCF(j,knxt,numStates,axis+2)
                  else
                     call readDataPSCF(j,knxt,numStates,axis+2)
                  endif
               enddo
            else ! Code 6, 7, 8: PlusG version
               do j = 1, spin
                  if (inSCF == 1) then
                     call readDataSCF(j,knxt,numStates,axis+5)
                  else
                     call readDataPSCF(j,knxt,numStates,axis+5)
                  endif
               enddo
            endif

            ! Copy the valeVale (waveFn coeffs) to a temporary matrix because
            !   we will reuse valeVale for the other (current) kpoint.
            CKnxt(:valeDim,:numStates,:) = valeVale(:valeDim,:numStates,:)

            ! Read in only the current kpoint wave function coefficients. Code 0.
            do j = 1, spin
               if (inSCF == 1) then
                  call readDataSCF(j,kcur,numStates,0)
               else
                  call readDataPSCF(j,kcur,numStates,0)
               end if
            enddo

!write(20,*) "k,j,CKcur,CKnxt"
!do j = 1, numStates
!do k = 1, valeDim
!write(20,*) k,j,valeVale(k,j,1),CKnxt(k,j,1)
!enddo
!enddo
!
!write(20,*) "k,j,axis,KOverlap"
!do j = 1, valeDim
!do k = 1, valeDim
!write(20,*) k,j,axis,valeValeKO(k,j,axis)
!enddo
!enddo

            ! Initialize a matrix to hold the outer product of the waveFn
            !   coefficients from the current and next kpoint.
            valeValePsi(:,:,:) = cmplx(0.0_double,0.0_double,double)

!write(20,*) "numStates maxOcc=", numStates, maxOccupiedState

            do j = 1, maxOccupiedState
               do k = 1, maxOccupiedState
!do n = 1, valeDim
!write(20,*) "n,CKnxt,vV",n, CKnxt(n,j,1), valeVale(n,k,1)
!enddo
                  ! Get the average population of the curr
                  currentPopulation(:) = &
                        & (structuredElectronPopulation(j,kcur,:) + &
                        & structuredElectronPopulation(k,kcur,:)) / 2.0_double
!write(20,*) "currentPopulation = ", currentPopulation(:)

                  ! Create a valeValePsi for each pair of occupied states.
                  do l = 1, spin
valeValePsi(:,:,l) = cmplx(0.0_double,0.0_double,double)
                     ! Perform a rank1 update A := alpha*x*y**H + A where A is
                     !   valeValePsi, alpha is 1+i0, x is the wave function
                     !   coefficients of the current kpoint, y is the wave
                     !   function coefficients, and **H is the conjugate
                     !   transpose operation.
                     !call zgerc(valeDim,valeDim,&
                     !      & cmplx(1.0_double,0.0_double,double),valeVale(:,k,l),1,&
                     !      & CKnxt(:,j,l),1,valeValePsi(:,:,l),valeDim)
                     call zgerc(valeDim,valeDim,cmplx(1.0d0,0.0_double,double),&
                           & CKnxt(:,j,l),1,valeVale(:,k,l),1,valeValePsi(:,:,l),&
                           & valeDim)

!do m = 1, valeDim
!do n = 1, valeDim
!write(20,*) "n,m,vvP,vvKO",n,m, valeValePsi(n,m,l), valeValeKO(n,m,l)
!enddo
!enddo

                     call fullMatrixElementMult (stateStateMat(k,j,l), &
                           & valeValePsi(:,:,l), valeValeKO(:,:,axis), valeDim)
                  enddo ! spin
               enddo ! k: maxOccupiedStates
            enddo ! j: maxOccupiedStates
write(20,*) "UNITARY: i,h = ", i,h
unitary = matmul(conjg(transpose(stateStateMat(:,:,1))),stateStatemat(:,:,1))
do m = 1, maxOccupiedState
write(23+axis,*) unitary(m,m)
enddo
!do m = 1, maxOccupiedState
!do n = 1, maxOccupiedState
!write(20,*) "n,m,sSM,pMA = ", n, m, stateStateMat(n,m,1), prodM_A(n,m,1)
!enddo
!enddo


            ! Accumulate this stateStateMat into the final product.

            ! (A) Perform the product part of -im(ln(det(prod(M)))).
!            do j = 1, spin
!               prodM_A(:,:,j) = matmul(prodM_A(:,:,j),stateStateMat(:,:,j))
!            enddo

            ! (B) Perform the product (and determinant) parts of:
            !   -im(ln(prod(det(M))))
            do j = 1, spin
write(20,*) "B matrixDet"
               phiString_B(j) = phiString_B(j)*matrixDet(stateStateMat(:,:,j))
            enddo

            ! (C) Perform the sum (and natural log and determinant) parts of:
            !   -im(sum(ln(det(M))))
            do j = 1, spin
write(20,*) "C matrixDet"
               phiString_C(j) = phiString_C(j) + &
                     & log(matrixDet(stateStateMat(:,:,j)))
            enddo

         enddo ! i = 1,numSteps (completing a line)

         ! (A) Take the determinent of the product of matrices for this string
         !   followed by the log and extraction of the negative imaginary
         !   component.
         do i = 1, spin
write(20,*) "A matrixDet"
!            stringPhaseSet_A(h,i) = -aimag(log(matrixDet(prodM_A(:,:,i))))
         enddo

         ! (B) Take the log and extract the negative of the imaginary
         !   component.
         do i = 1, spin
            stringPhaseSet_B(h,i) = -aimag(log(phiString_B(i)))
         enddo

         ! (C) Extract the negative of the imaginary component.
         do i = 1, spin
            stringPhaseSet_C(h,i) = -aimag(phiString_C(i))
         enddo
write(20,*) "h, A,B,C Phase", h, stringPhaseSet_A(h,1), &
      & stringPhaseSet_B(h,1),stringPhaseSet_C(h,1)
call flush(20)

         ! Mark the completion of this string.
         !if (mod(h,10) .eq. 0) then
         !   write (20,ADVANCE="NO",FMT="(a1)") "|"
         !else
         !   write (20,ADVANCE="NO",FMT="(a1)") "."
         !endif
         !if (mod(h,50) .eq. 0) then
         !   write (20,*) " ",h
         !endif
         !call flush (20)

      enddo ! h (iterating over the set of lines)
      write(20,*)

!write(20,*) "Doing A"
!call flush(20)
!      call getAveragePhase(currNumLines,axis,stringPhaseSet_A,averagePhase_A)
write(20,*) "Doing B"
call flush(20)
      call getAveragePhase(currNumLines,axis,stringPhaseSet_B,averagePhase_B)
write(20,*) "Doing C"
call flush(20)
      call getAveragePhase(currNumLines,axis,stringPhaseSet_C,averagePhase_C)

!      do i = 1, spin
!!write(20,*) "axis psi currNumLines = ", axis, psi(axis,i), currNumLines
!         psi(axis,i)=psi(axis,i)/currNumLines
!
!!write(20,*) "axis psi = ", axis, psi(axis,i)
!         fn(axis,i) = psi(axis,i)/(2.0d0*pi)
!
!!write(20,*) "axis fn = ", axis, fn(axis,i)
!         fn(axis,i) = modulo(fn(axis,i) + 0.5d0, 1.0d0) - 0.5d0
!
!!write(20,*) "axis fn = ", axis, fn(axis,i)
!         xyzP(:,i) = xyzP(:,i) + (-eCharge/realCellVolume) &
!               & * fn(axis,i) * realVectors(:,axis)
!
!!write(20,*) "axis xyzP = ", axis, xyzP(:,i)
!      enddo

!         psi(axis) = psi(axis) - 2.0d0*pi*dnint( psi(axis)/(2.0d0*pi) )

      ! Convert to xyz and apply prefactors and conversions.
      do h = 1, spin
!         xyzP_A(:,h) = xyzP_A(:,h) + (-eCharge/realCellVolume) * &
!               & averagePhase_A(axis,h) / pi / spin * realVectors(:,axis)
         xyzP_B(:,h) = xyzP_B(:,h) + (-eCharge/realCellVolume) * &
               & averagePhase_B(axis,h) / pi / spin * realVectors(:,axis)
         xyzP_C(:,h) = xyzP_C(:,h) + (-eCharge/realCellVolume) * &
               & averagePhase_C(axis,h) / pi / spin * realVectors(:,axis)
      enddo
   enddo ! axis

   ! Apply prefactors and conversions to ionic contributions.
   xyzIonMoment(:) = xyzIonMoment(:) * eCharge/realCellVolume
   abcIonMoment(:) = abcIonMoment(:) * eCharge/realCellVolume

   do i = 1, spin
!      write(20,*) 'xyzP [C/m^2] = ', xyzP(:,i)
      !xyzP(:,i) = xyzP(:,i) * (10.0d0/(bohrRad**2))
!      write(20,*) 'xyzP_A [C/m^2] = ', xyzP_A(:,i)
      write(20,*) 'xyzP_B [C/m^2] = ', xyzP_B(:,i)
      write(20,*) 'xyzP_C [C/m^2] = ', xyzP_C(:,i)
      write(20,*) 'xyzIonMoment = ', xyzIonMoment(:)
      write(20,*) 'abcIonMoment = ', abcIonMoment(:)
!      write(20,*) 'Dipole Moment A = ', xyzIonMoment(:) - xyzP_A(:,i)
      write(20,*) 'Dipole Moment B+ = ', xyzIonMoment(:) + xyzP_B(:,i)
      write(20,*) 'Dipole Moment C+ = ', xyzIonMoment(:) + xyzP_C(:,i)
      write(20,*) 'Dipole Moment B- = ', xyzIonMoment(:) - xyzP_B(:,i)
      write(20,*) 'Dipole Moment C- = ', xyzIonMoment(:) - xyzP_C(:,i)
      write(20,*) 'Dipole Moment C/m^2 B+ = ', (xyzIonMoment(:) + xyzP_B(:,i)) * 10.0d0/(bohrRad**2)
      write(20,*) 'Dipole Moment C/m^2 C+ = ', (xyzIonMoment(:) + xyzP_C(:,i)) * 10.0d0/(bohrRad**2)
      write(20,*) 'Dipole Moment C/m^2 B- = ', (xyzIonMoment(:) - xyzP_B(:,i)) * 10.0d0/(bohrRad**2)
      write(20,*) 'Dipole Moment C/m^2 C- = ', (xyzIonMoment(:) - xyzP_C(:,i)) * 10.0d0/(bohrRad**2)
   enddo

   xyzP = xyzP_B

   ! Deallocate results. Note that valeVale is deallocated later.
   deallocate(CKnxt)
   deallocate(valeValePsi)
   deallocate(valeValeKO)
   deallocate(structuredElectronPopulation)
   deallocate(phiString_B)
   deallocate(phiString_C)
   deallocate(currentPopulation)
   deallocate(stringPhaseSet_A)
   deallocate(stringPhaseSet_B)
   deallocate(stringPhaseSet_C)
   deallocate(stateStateMat)
   deallocate(prodM_A)

   ! Record the end of the calculation
   call timeStampEnd(33)

end subroutine computeMTOPPolarization
#endif

subroutine getAveragePhase(currNumLines,axis,stringPhaseSet,averagePhase)

   ! Use necessary modules
   use O_Kinds
   use O_Constants, only: pi, eCharge, smallThresh
   use O_SortSubs, only: mergeSort
   use O_Potential, only: spin
   use O_AtomicSites, only: abcIonMoment
   use O_Lattice, only: realMag

   implicit none

   ! Define passed parameters.
   integer, intent(in) :: currNumLines
   integer, intent(in) :: axis
   real (kind=double), dimension(:,:), intent(inout) :: stringPhaseSet
   real (kind=double), dimension(3,spin), intent(inout) :: averagePhase

   ! Define local variables.
   integer :: h, i
   integer :: numSegments
   integer, allocatable, dimension(:) :: segmentBorders
   integer, allocatable, dimension(:) :: indexStringPhaseSet
   real(kind=double), allocatable, dimension(:) :: sortedStringPhaseSet


   ! Now that the phase for each line along this axis has been determined,
   !   we can take an average of the phase across all lines. However,
   !   before we can do that, we need to make sure that all the phases are
   !   part of the same bunch and that none of them are part of a separate
   !   bunch that was shifted by a quantum of polarization (2*pi).
   ! We do that through a series of steps. First, we sort the phases. Then,
   !   we shift all phases to be within 2pi of the first (lowest magnitude)
   !   phase by quanta of 2pi. Then, we check if one more shift by 2pi for
   !   each phase will bring it closer to the first phase or not. If it
   !   does, then that shifted phase may have a lower magnitude than the
   !   first, but we don't sort any thing. We just leave the phase in the
   !   same index it is currently in. After the per-phase shifting, to
   !   make all phases as close as possible to the first (originially
   !   lowest) phase, then we take the average of the phases. Finally,
   !   we shift that average to be as close as possible to 0.0.

   ! Use numSegments in place of currNumLines to preserve currNumLines.
   numSegments = currNumLines

   ! Allocate space for the sorting process.
   allocate(segmentBorders(currNumLines+1))
   allocate(indexStringPhaseSet(currNumLines))
   allocate(sortedStringPhaseSet(currNumLines))

   ! Initialize the segment borders.
   segmentBorders(1) = 0
   do i = 1, numSegments
      segmentBorders(i+1) = i
   enddo

   do h = 1, spin

write(20,*) "Got here 1a"
call flush(20)
do i = 1,currNumLines
write(20+axis,*) stringPhaseSet(i,h),1
enddo
      ! Sort the string phases and copy the sorted list over the original.
      !call mergeSort(stringPhaseSet(:,h),sortedStringPhaseSet,&
      !      & indexStringPhaseSet,segmentBorders,numSegments)
      !stringPhaseSet(:,h) = sortedStringPhaseSet(:)
write(20,*) "Got here 1b: axis abcIon", axis, abcIonMoment(axis)
call flush(20)
write(20,*) stringPhaseSet(:,h)
write(20,*) stringPhaseSet(:,h) / pi / spin * realMag(axis)

      ! Shift each phase by +/- 2pi * n with n = 0, 1, 2, 3, ...
      !   The goal is for each phase/pi/spin to be as close as possible to the
      !   abc ionic moment.
      do i = 1, currNumLines
write(20,*) "line i", i
         do while (stringPhaseSet(i,h)/pi/spin*realMag(axis) > abcIonMoment(axis))
            stringPhaseSet(i,h) = stringPhaseSet(i,h) - 2.0_double*pi
write(20,*) "Step1 ", stringPhaseSet(i,h) / pi / spin*realMag(axis)
         enddo
         do while (stringPhaseSet(i,h)/pi/spin*realMag(axis) < abcIonMoment(axis))
            stringPhaseSet(i,h) = stringPhaseSet(i,h) + 2.0_double*pi
write(20,*) "Step2 ", stringPhaseSet(i,h) / pi / spin*realMag(axis)
         enddo
         if (abs((stringPhaseSet(i,h)+2.0d0*pi)/pi/spin*realMag(axis) - abcIonMoment(axis))&
               & < abs(stringPhaseSet(i,h)/pi/spin*realMag(axis) - abcIonMoment(axis))) then
            stringPhaseSet(i,h) = stringPhaseSet(i,h) + 2.0_double*pi
write(20,*) "Step3 ", stringPhaseSet(i,h) / pi / spin*realMag(axis)
         endif
         if (abs((stringPhaseSet(i,h)-2.0d0*pi)/pi/spin*realMag(axis) - abcIonMoment(axis))&
               & < abs(stringPhaseSet(i,h)/pi/spin*realMag(axis) - abcIonMoment(axis))) then
            stringPhaseSet(i,h) = stringPhaseSet(i,h) - 2.0_double*pi
write(20,*) "Step4 ", stringPhaseSet(i,h) / pi / spin*realMag(axis)
         endif

!         do while ((stringPhaseSet(i,h) / pi / spin - abcIonMoment(axis)) > pi)
!            stringPhaseSet(i,h) = stringPhaseSet(i,h) - 2.0_double*pi
!write(20,*) "Step1 ", stringPhaseSet(i,h) / pi / spin
!         enddo
!         do while ((stringPhaseSet(i,h) / pi / spin - abcIonMoment(axis)) < -pi)
!            stringPhaseSet(i,h) = stringPhaseSet(i,h) + 2.0_double*pi
!write(20,*) "Step2 ",stringPhaseSet(i,h) / pi / spin
!         enddo
!         if ((stringPhaseSet(i,h) / pi / spin - abcIonMoment(axis)) > &
!               & (abcIonMoment(axis)-(stringPhaseSet(i,h)-pi))) then
!            stringPhaseSet(i,h) = stringPhaseSet(i,h) - 2.0_double*pi
!write(20,*) "Step3 ",stringPhaseSet(i,h) / pi / spin
!         endif
!         if ((stringPhaseSet(i,h) / pi / spin - abcIonMoment(axis)) > &
!               & (abcIonMoment(axis)-(stringPhaseSet(i,h)+pi))) then
!            stringPhaseSet(i,h) = stringPhaseSet(i,h) + 2.0_double*pi
!write(20,*) "Step4 ",stringPhaseSet(i,h) / pi / spin
!         endif
      enddo
!      ! Shift each phase to be as close as possible but still greater than the
!      !   first. Then, check if one more shift of 2pi will bring it closer,
!      !   even if it makes it less than the first phase in the sorted array.
!      do i = 2, currNumLines
!         do while ((stringPhaseSet(i,h) - stringPhaseSet(1,h)) > 2.0_double*pi)
!            stringPhaseSet(i,h) = stringPhaseSet(i,h) - 2.0_double*pi
!         enddo
!         if ((stringPhaseSet(i,h) - stringPhaseSet(1,h)) > &
!               & (stringPhaseSet(1,h)-(stringPhaseSet(i,h)-2.0d0*pi))) then
!            stringPhaseSet(i,h) = stringPhaseSet(i,h) - 2.0_double*pi
!         endif
!      enddo
write(20,*) "Got here 2b"
call flush(20)
do i = 1,currNumLines
write(30+axis,*) stringPhaseSet(i,h),1
enddo

      ! Compute the average phase and then ensure that it is >-pi, <+pi.
      averagePhase(axis,h) = sum(stringPhaseSet(1:currNumLines,h)) / &
            & currNumLines
write(20,*) "axis, Avg Phase = ", axis, averagePhase(axis,h)
      !do while (averagePhase(axis,h) < -pi)
      !   averagePhase(axis,h) = averagePhase(axis,h) + 2.0_double*pi
      !enddo
      !do while (averagePhase(axis,h) > pi)
      !   averagePhase(axis,h) = averagePhase(axis,h) - 2.0_double*pi
      !enddo
write(20,*) "axis, Avg Phase/pi/spin = ", axis, averagePhase(axis,h)/pi/spin*realMag(axis)
write(20,*) "Got here 3"
call flush(20)
   enddo

   ! Deallocate space for the sorting process.
   deallocate(segmentBorders)
   deallocate(indexStringPhaseSet)
   deallocate(sortedStringPhaseSet)

!   ! Now that the phase for each line along this axis has been determined,
!   !   we can take an average of the phase across all lines. However,
!   !   before we can do that, we need to make sure that all the phases are
!   !   part of the same bunch and that none of them are part of a separate
!   !   bunch that was shifted by a quantum of polarization (2*pi).
!   ! We do that by a series of steps. First, we shift all phases into the
!   !   range 0,+2pi. Then, we compare all phases to the first phase and
!   !   ensure that each is no more than +/- pi away from it. Finally, we
!   !   comute the average phase and then shift it to be within the range
!   !   -pi,+pi.
!   ! This will work regardless of the "spread" of phases in the inital
!   !   list. But, it will never let the phase "jump" by 2pi if we carry out
!   !   a series of calculations (e.g., compression, tension, etc.).
!   do h = 1, spin
!      ! Step one. Ensure all phases are in the 0 to +2pi range.
!      do i = 1, currNumLines
!write(20,*) "Stage1: i,h,sPS = ",i,h,stringPhaseSet(i,h)
!         ! Ensure that the phase is positive.
!         do while (stringPhaseSet(i,h) < 0.0_double)
!            stringPhaseSet(i,h) = stringPhaseSet(i,h) + 2.0_double*pi
!         enddo
!         ! Ensure that the phase is between 0 and 2pi.
!         stringPhaseSet(i,h) = modulo(stringPhaseSet(i,h),2.0_double*pi)
!      enddo
!
!      ! Step two. Ensure that the phases are no more than pi apart.
!write(20,*) "Stage2: i,h,sPS = ",1,h,stringPhaseSet(1,h)
!      do i = 2, currNumLines
!write(20,*) "Stage2: i,h,sPS = ",i,h,stringPhaseSet(i,h)
!         if (abs(stringPhaseSet(i,h) - stringPhaseSet(1,h)) > pi) then
!            if (stringPhaseSet(i,h) > stringPhaseSet(1,h)) then
!               stringPhaseSet(i,h) = stringPhaseSet(i,h) - 2.0_double*pi
!            else
!               stringPhaseSet(i,h) = stringPhaseSet(i,h) + 2.0_double*pi
!            endif
!         endif
!      enddo
!
!do i = 1, currNumLines
!write(20,*) "Stage3: i,h,sPS = ",i,h,stringPhaseSet(i,h)
!enddo
!
!      ! Step three. Compute the average and then ensure that it is >-pi, <+pi.
!      averagePhase(axis,h) = sum(stringPhaseSet(1:currNumLines,h)) / &
!            & currNumLines
!write(20,*) "axis, Avg Phase = ", axis, averagePhase(axis,h)
!      if (averagePhase(axis,h) < -pi) then
!         averagePhase(axis,h) = averagePhase(axis,h) + 2.0_double*pi
!      elseif (averagePhase(axis,h) > pi) then
!         averagePhase(axis,h) = averagePhase(axis,h) - 2.0_double*pi
!      endif
!write(20,*) "axis, Avg Phase = ", axis, averagePhase(axis,h)
!
!      ! Step four. Convert to xyz and apply prefactors and conversions from
!      !   atomic units (bohr radii) to ...
!      xyzP(:,h) = xyzP(:,h) + (-eCharge/realCellVolume) * &
!            & averagePhase(axis,h) / pi / spin * realVectors(:,axis)
!   enddo

end subroutine getAveragePhase


function matrixDet(A)

   ! Define modules to use.
   use O_Kinds

   ! Make sure no variables are accidentally declared.
   implicit none

   ! Declare return variable.
   complex(kind=double) :: matrixDet

   ! Declare passed parameters.
   complex(kind=double), intent(in) :: A(:,:)

   ! Define local variables.
   integer :: n,i,info,istat
   integer :: swaps,j,p
   complex(kind=double), allocatable :: Ac(:,:)
   integer,              allocatable :: ipiv(:)
   logical,              allocatable :: visited(:)
   real(kind=double), parameter :: tiny = 1.0d-14
   !real(kind=double), parameter :: tiny = 1.0d-300
   external :: zgetrf

   matrixDet = (1.0d0,0.0d0)

   if (size(A,1) /= size(A,2)) return
   n = size(A,1); if (n==0) return

   allocate(Ac(n,n), ipiv(n), visited(n), stat=istat)
      if (istat/=0) stop "Allocate failed."
   Ac = A
   call zgetrf(n,n,Ac,n,ipiv,info); if (info/=0) then
      deallocate(Ac,ipiv,visited)
      write (20,*) "info=",info
      write (20,*) "LU Decomposition failed."
      stop
   end if

   do i=1,n
      if (abs(Ac(i,i)) <= tiny) then
         deallocate(Ac,ipiv,visited)
         write(20,*) i, "Diagonal element is too small."
         stop
      end if
      matrixDet = matrixDet * Ac(i,i)
   end do

!write(20,*) "Squaring det.", matrixDet
!   matrixDet = matrixDet * matrixDet
!write(20,*) "Squared det.", matrixDet
!return

!write(20,*) "ipiv = ",ipiv(:)
   visited = .false.; swaps = 0
   do i=1,n
      if (.not. visited(i)) then
         j = i; p = 0
         do while (.not. visited(j))
            visited(j) = .true.; p = p + 1
            j = ipiv(j)
         end do
         if (p>0) swaps = swaps + (p-1)
      end if
   end do
   if (mod(swaps,2)==1) then
      write(20,*) "Odd swaps: ", swaps
      matrixDet = -matrixDet
   endif

   !matrixDet = matrixDet / &
   !      & cmplx(max(abs(matrixDet), tiny), 0.0d0, kind=double) ! unit phasor

   deallocate(Ac, ipiv, visited)
end function matrixDet

end module O_MTOP
