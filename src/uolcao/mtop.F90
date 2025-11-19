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
   use O_Lattice, only: recipVectors, realVectors, realCellVolume,&
         & invRealVectors
   use O_Potential,   only: spin
   use O_Populate,    only: electronPopulation
   use O_KPoints,     only: numKPoints, numAxialKPoints, mtopKPMap, kPoints
   use O_AtomicSites, only: valeDim, numAtomSites, atomSites, &
         & computeIonicMoment, xyzIonMoment
   use O_AtomicTypes, only: numAtomTypes, atomTypes, maxNumValeStates
   use O_Input,       only: numStates
   use O_MatrixSubs,  only: fullMatrixElementMult
   use O_SecularEquation, only: valeVale, valeValeKO, readDataSCF, readDataPSCF
   use O_BLASZGERU

   ! Make sure that no funny variables are used.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: inSCF
   real(kind=double), intent(out) :: xyzP(3,2)

   ! Define local variables.
   integer :: h,i,j,k,l
integer :: m,n
   integer :: maxOccupiedState
   real(kind=double), dimension(2) :: psiAxis
   real(kind=double), dimension (3,2) :: psi
   real(kind=double), allocatable, dimension (:) :: currentPopulation
   real (kind=double), allocatable, dimension (:,:,:) :: &
         & structuredElectronPopulation
   complex(kind=double) :: dotProduct 
   complex(kind=double) :: ld_step
   complex(kind=double), allocatable, dimension(:) :: phiString
   complex(kind=double), allocatable, dimension(:,:,:) :: valeValePsi
   complex(kind=double), allocatable, dimension(:,:,:) :: CKnxt
   complex(kind=double), allocatable, dimension(:,:,:) :: stateStateSum
   integer :: fixed2d,numStep
   integer :: energyLevelCounter,skipKP
   integer :: kcur, knxt, axis
   integer :: loop

   integer :: kPointCount

   real(kind=double) :: fn(3,2)


   ! Record the start of the calculation
   call timeStampStart(33)


   ! allocate the Eigenvectors for sequences C
   ! CK and CnextK are eigenvectors at k and k+Δk
   allocate (CKnxt(valeDim,numStates,spin))
   allocate (valeValePsi(valeDim,valeDim,spin))
   allocate (valeValeKO(valeDim,valeDim,3))
   if (inSCF == 0) then
      allocate (valeVale(valeDim,numStates,spin))
   endif
   allocate (structuredElectronPopulation (numStates,numKPoints,spin))
   allocate (phiString(spin))
   allocate (currentPopulation(spin))
   allocate (stateStateSum(numStates,numStates,spin))

   valeVale(:,:,:) = cmplx(0.0_double,0.0_double,double)
   valeValeKO(:,:,:) = cmplx(0.0_double,0.0_double,double)
   stateStateSum(:,:,:) = cmplx(0.0_double,0.0_double,double)
 
   psi(:,:) = 0.0d0
   xyzP(:,:) = 0.0d0

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

! For polarization along each axis the code fixes the other two k-coordinates
!   and generates a set of parallel closed loops.Along x (numStep), it holds Ky
!   and Kz constant (fixed2d), producing  Ky*Kz independent loops. Each loop
!   steps through Kx and wraps, so the number of steps is  Kx+1 wih the final
!   point equal to the first.Along each link between consecutive k-points the
!   code forms the occupied-subspace overlap matrix and takes its determinant,
!   which collapses the gauge freedom to a single complex link value.When the
!   loop closes, the imaginary part of the accumulated logarithms gives the
!   Berry phase for that loop.The same procedure applies along y and z . For
!   y, the loops fix Kx and Kz giving Kx*Kz loops, each with Ky+1 steps, and
!   the same determinant and logarithm accumulation. For z, the loops fix Kx
!   and Ky, giving Kx*Ky loops, each with Kz+1 steps, again accumulating log
!   det along the links and taking the imaginary part at closure.


   do axis = 1, 3
      if (axis == 1) then
         numStep = numAxialKPoints(1)
         fixed2d = numAxialKPoints(2) * numAxialKPoints(3)
      elseif (axis == 2) then
         numStep = numAxialKPoints(2)
         fixed2d = numAxialKPoints(1) * numAxialKPoints(3)
      else
         numStep = numAxialKPoints(3)
         fixed2d = numAxialKPoints(1) * numAxialKPoints(2)
      endif

      kPointCount = 0

      do loop = 1, fixed2d
         phiString(:) = (1.0d0, 0.0d0)
         psiAxis(:) = 0.0d0

         do i = 1, numStep
            kPointCount = kPointCount + 1
            kcur = mtopKPMap(axis,kPointCount)
            if (i < numStep) then
               knxt = mtopKPMap(axis,kPointCount+1)
            else
               knxt = mtopKPMap(axis,kPointCount-numStep+1)
            endif

!write(20,*) "kcur knxt", kcur, knxt
!write(20,*) "kcur", kPoints(:,kcur)
!write(20,*) "knxt", kPoints(:,knxt)

            ! Skip any kpoints with a negligable contribution for each state.
            skipKP = 1 ! Assume that we will skip these kpoints.
            do j = 1, numStates
               if ((sum(abs(structuredElectronPopulation(j,kcur,:))) > &
                     & smallThresh) .and. &
                     & (sum(abs(structuredElectronPopulation(j,knxt,:))) > &
                     & smallThresh)) then
                  skipKP = 0 ! If both kp contribute enough then don't skip
                  exit
               endif
            enddo
            if (skipKP == 1) then
               cycle
            endif

            ! Read next kpoint waveFn coefficients and KOverlap. Code 3,4,5
            do h = 1, spin
               if (inSCF == 1) then
                  call readDataSCF(h,knxt,numStates,axis+2)
               else
                  call readDataPSCF(h,knxt,numStates,axis+2)
               endif
            enddo

            ! Copy the valeVale (waveFn coeffs) to a temporary matrix because
            !   we will reuse valeVale for the other (current) kpoint.
            !if (i < numStep) then
               CKnxt(:valeDim,:numStates,:) = valeVale(:valeDim,:numStates,:)
            !else
            !   dotProduct = dot_product(recipVectors(:,axis),r)
            !   CKnxt(:valeDim,:numStates,:) = &
            !         & cmplx(cos(dotProduct),sin(dotProduct),double) * &
            !         & valeVale(:valeDim,:numStates,:)
            !endif

            ! Read in only the current kpoint wave function coefficients. Code 0.
            do h = 1, spin
               if (inSCF == 1) then
                  call readDataSCF(h,kcur,numStates,0)
               else
                  call readDataPSCF(h,kcur,numStates,0)
               end if
            enddo

            ! Initialize a matrix to hold the outer product of the waveFn
            !   coefficients from the current and next kpoint.
            valeValePsi(:,:,:) = cmplx(0.0_double,0.0_double,double)

            do j = 1, numStates

               currentPopulation(:) = structuredElectronPopulation(j,kcur,:)
               if (sum(abs(currentPopulation(:))) < smallThresh) then
                  maxOccupiedState = j-1
                  exit
               endif
            enddo

!write(20,*) "numStates maxOcc=", numStates, maxOccupiedState

            do j = 1, maxOccupiedState
               do k = 1, maxOccupiedState
!do n = 1, valeDim
!write(20,*) "n,CKnxt,vV",n, CKnxt(n,j,1), valeVale(n,k,1)
!enddo
                  currentPopulation(:) = &
                        & (structuredElectronPopulation(j,kcur,:) + &
                        & structuredElectronPopulation(k,kcur,:)) / 2.0_double

                  ! Create a valeValePsi for each pair of occupied states.
                  do l = 1, spin
                     call zgerc(valeDim,valeDim,currentPopulation(l),&
                           & CKnxt(:,j,l),1,valeVale(:,k,l),1,valeValePsi(:,:,l),&
                           & valeDim)

!do m = 1, valeDim
!do n = 1, valeDim
!write(20,*) "n,m,vvP,vvKO",n,m, valeValePsi(n,m,l), valeValeKO(n,m,l)
!enddo
!enddo

                     call fullMatrixElementMult (stateStateSum(k,j,l), &
                           & valeValePsi(:,:,l), valeValeKO(:,:,axis), valeDim)
                  enddo
               enddo
            enddo

!do n = 1, maxOccupiedState
!do m = 1, maxOccupiedState
!write(20,*) "m n sSS=",m,n,stateStateSum(m,n,1)
!enddo
!enddo

!write(20,*) "i axis phiString = ", i, axis, phiString(1)

            do j = 1, spin
               call phase_from_matrix(stateStateSum(:maxOccupiedState, &
                     & :maxOccupiedState,j), ld_step)
               phiString(j) = phiString(j) * ld_step
            enddo

!write(20,*) "ld_step phiString = ", ld_step, phiString(1)

         end do ! i = 1,numStep

         do i = 1, spin
            psi(axis,i) = psi(axis,i) + aimag(log(phiString(i)))
         enddo

!write(20,*) "axis psi = ", axis, psi(axis,1)

      end do ! loop

      do i = 1, spin
!write(20,*) "axis psi fixed2d = ", axis, psi(axis,i), fixed2d
         psi(axis,i)=psi(axis,i)/fixed2d

!write(20,*) "axis psi = ", axis, psi(axis,i)
         fn(axis,i) = psi(axis,i)/(2.0d0*pi)

!write(20,*) "axis fn = ", axis, fn(axis,i)
         fn(axis,i) = modulo(fn(axis,i) + 0.5d0, 1.0d0) - 0.5d0

!write(20,*) "axis fn = ", axis, fn(axis,i)
         xyzP(:,i) = xyzP(:,i) + (-eCharge/realCellVolume) &
               & * fn(axis,i) * realVectors(:,axis)

!write(20,*) "axis xyzP = ", axis, xyzP(:,i)
      enddo

!         psi(axis) = psi(axis) - 2.0d0*pi*dnint( psi(axis)/(2.0d0*pi) )
   enddo ! axis

   call computeIonicMoment

   do i = 1, spin
!      write(20,*) 'xyzP [C/m^2] = ', xyzP(:,i)
      xyzP(:,i) = xyzP(:,i) * (10.0d0*eCharge/(bohrRad**2))
      write(20,*) 'xyzP [C/m^2] = ', xyzP(:,i)
      write(20,*) 'xyzIonMoment = ', xyzIonMoment(:)
      write(20,*) 'Dipole Moment = ', xyzIonMoment(:) + xyzP(:,i)
   enddo

   deallocate(CKnxt)
   deallocate(valeValeKO)
   deallocate(structuredElectronPopulation)
   deallocate(phiString)
   deallocate(currentPopulation)

   ! Record the end of the calculation
   call timeStampEnd(33)

end subroutine computeMTOPPolarization
#endif

subroutine phase_from_matrix(A, detU)

   ! Define modules to use.
   use O_Kinds

   ! Make sure no variables are accidentally declared.
   implicit none

   ! Declare passed parameters.
   complex(kind=double), intent(in) :: A(:,:)
   complex(kind=double), intent(out) :: detU

   ! Define local variables.
   integer :: n,i,info,istat,swaps,j,p
   complex(kind=double), allocatable :: Ac(:,:)
   integer,              allocatable :: ipiv(:)
   logical,              allocatable :: visited(:)
   real(kind=double), parameter :: tiny = 1.0d-14
   external :: zgetrf

   detU = (1.0d0,0.0d0)

   if (size(A,1) /= size(A,2)) return
   n = size(A,1); if (n==0) return

!write(20,*) "n=",n
!do i=1,n
!   do j = 1, n
!      write(20,*) "j,i,A(j,i)=", j, i, A(j,i)
!   enddo
!enddo

   allocate(Ac(n,n), ipiv(n), visited(n), stat=istat)
      if (istat/=0) stop "Allocate failed."
   Ac = A
   call zgetrf(n,n,Ac,n,ipiv,info); if (info/=0) then
      deallocate(Ac,ipiv,visited)
      write (20,*) "info=",info
      stop "LU Decomposition failed."
   end if

!do i=1,n
!   do j = 1,n
!      write(20,*) "j,i,Ac(j,i)=", j, i, Ac(j,i)
!   enddo
!enddo
!
!do i=1,n
!   write(20,*) "i,Ac(i,i)=", i, Ac(i,i)
!enddo

   do i=1,n
      if (abs(Ac(i,i)) <= tiny) then
         deallocate(Ac,ipiv,visited)
         write(20,*) i, "Diagonal element is too small."
         return
      end if
      detU = detU * Ac(i,i)
   end do

!write(20,*) "Squaring det.", detU
!   detU = detU * detU

!   visited = .false.; swaps = 0
!   do i=1,n
!      if (.not. visited(i)) then
!         j = i; p = 0
!         do while (.not. visited(j))
!            visited(j) = .true.; p = p + 1
!            j = ipiv(j)
!         end do
!         if (p>0) swaps = swaps + (p-1)
!      end if
!   end do
!   if (mod(swaps,2)==1) detU = -detU
!
!   detU = detU / cmplx(max(abs(detU), tiny), 0.0d0, kind=double)  ! unit phasor

   deallocate(Ac, ipiv, visited)
end subroutine phase_from_matrix

end module O_MTOP
