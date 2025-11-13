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
subroutine computeMTOPPolarization(inSCF,P)

      ! Import necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Constants,   only: pi, eCharge, smallThresh, bohrRad
   use O_Lattice, only: recipVectors, realVectors, realCellVolume,&
         & invRealVectors
   use O_Potential,   only: spin
   use O_Populate,    only: electronPopulation
   use O_KPoints,     only: numKPoints, numAxialKPoints, MTOPIndexMap
   use O_AtomicSites, only: valeDim, numAtomSites, atomSites
   use O_AtomicTypes, only: numAtomTypes, atomTypes, maxNumValeStates
   use O_Input,       only: numStates
   use O_SecularEquation, only: valeVale, valeValeKO, readDataSCF, readDataPSCF
   use O_BLASZGERU

   ! Make sure that no funny variables are used.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: inSCF
   real(kind=double), intent(out) :: P(3,2)

   ! Define local variables.
   integer :: h,i,j,k
   real(kind=double), dimension(2) :: psiAxis
   real(kind=double), dimension (3,2) :: psi
   real(kind=double), allocatable, dimension (:) :: currentPopulation
   real (kind=double), allocatable, dimension (:,:,:) :: &
         & structuredElectronPopulation
   complex(kind=double) :: ld_step
   complex(kind=double), allocatable, dimension(:,:) :: phiString
   complex(kind=double), allocatable, dimension(:,:,:) :: valeValePsi

   ! Output
   complex(kind=double), allocatable, dimension(:,:,:) :: CK
   !complex(kind=double), allocatable, dimension(:,:) :: CnextK
   complex(kind=double), allocatable, dimension(:,:,:) :: Sdir
   !complex(kind=double), allocatable, dimension(:,:) :: Sstep

   !Define local variables .
   integer :: u,v, s0,s1 ! loop indices
!   integer :: n1,n2,n3,XY,YZ,XZ, linkYZ,linkXZ,linkXY
   integer, dimension(3) :: n
   integer :: Fixed2D,numstep,n1,n2,n3

   integer :: energyLevelCounter,skipKP,m,mp1
   integer :: s, s_start, s_end, kcur, knxt, axis,maxLoops, axis_eff
   integer :: loop, numLoops, numKPinLoop, s_start_block, TotLink
   integer :: dirA, dirB, Segment1, Segment2, e1_end, e2_end, e3_end
   integer :: s_in_axis, s_in_loop,l,numSeg1BySeg2
   integer :: step,s_loc,LoopCount,stepPerLoop

   real(kind=double) :: psi_axis, df(3), fn(3,2)
   real(kind=double) :: c12(3), c23(3), c31(3)

   complex(kind=double), allocatable, dimension(:,:) :: T,F,TG,FG

   real(kind=double), dimension(3,2) :: PCart
   complex(kind=double)    :: theta_step

!   write(20,*)'********************************************************************'


   ! allocate the Eigenvectors for sequences C
   ! CK and CnextK are eigenvectors at k and k+Δk
   allocate (CK(valeDim,numStates,spin))
   allocate (valeValePsi(valeDim,valeDim,spin))
   allocate (valeValeKO(valeDim,valeDim,3))
   if (inSCF == 0) then
      allocate (valeVale(valeDim,numStates,spin))
   endif
   allocate (structuredElectronPopulation (numStates,numKPoints,spin))
   allocate (phiString(3,spin))
   allocate (currentPopulation(spin))

   valeVale(:,:,:) = cmplx(0.0_double,0.0_double,double)
   valeValeKO(:,:,:) = cmplx(0.0_double,0.0_double,double)
 
   P(:,:) = 0.0d0
   psi(:,:) = 0.0d0

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

   n1=numAxialKPoints(1)
   n2=numAxialKPoints(2)
   n3=numAxialKPoints(3)


!For polarization along each axis the code fixes the other two k-coordinates and generates a set of 
!parallel closed loops.Along x (numstep), it holds Ky and Kz constant (Fixed2D), producing  Ky*Kz 
!independent loops. Each loop steps through Kx and wraps, so the number of steps is  Kx+1 wih the 
!final point equal to the first.Along each link between consecutive k-points the code forms the 
!occupied-subspace overlap matrix and takes its determinant, which collapses the gauge freedom to 
!a single complex link value.When the loop closes, the imaginary part of the accumulated logarithms 
!gives the Berry phase for that loop.The same procedure applies along y and z . For y, the loops fix 
!Kx and Kz giving Kx*Kz loops, each with Ky+1 steps, and the same determinant and logarithm accumulation. 
!For z, the loops fix Kx and Ky, giving Kx*Ky loops, each with Kz+1 steps, again accumulating log det 
!along the links and taking the imaginary part at closure.


   do axis = 1, 3
      if (axis == 1) then
         numstep = n1
         Fixed2D = n2*n3
      elseif (axis == 2) then
         numstep = n2
         Fixed2D = n1*n3
      else
         numstep = n3
         Fixed2D = n1*n2
      endif

      do loop = 1, Fixed2D
         phiString(axis,:) = (1.0d0, 0.0d0)
         psiAxis(:) = 0.0d0

         do i = 1, numstep 
            s_start = (loop-1)*(numstep+1) !+ 1

            kcur = MTOPIndexMap(s_start+i, axis)!; if (kcur<=0) exit
            knxt = MTOPIndexMap(s_start+i+1, axis)

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

            ! Read current kpoint waveFn coefficients and KOverlap. Code 3,4,5
            do h = 1, spin
               if (inSCF == 1) then
                  call readDataSCF(h,kcur,numStates,axis+2)
               else
                  call readDataPSCF(h,kcur,numStates,axis+2)
               endif
            enddo

            ! Copy the valeVale (waveFn coeffs) to a temporary matrix because
            !   we will reuse valeVale for the next kpoint.
            CK(:,1:numStates,:) = valeVale(1:valeDim,1:numStates,:)

            ! Read in only the next kpoint wave function coefficients. Code 0.
            do h = 1, spin
               if (inSCF == 1) then
                  call readDataSCF(h,knxt,numStates,0)
               else
                  call readDataPSCF(h,knxt,numStates,0)
               end if
            enddo

            ! Initialize a matrix to hold the outer product of the waveFn
            !   coefficients from the current and next kpoint.
            valeValePsi(:,:,:) = cmplx(0.0_double,0.0_double,double)

            ! Perform the rank-1 update: A = alpha * x * y**H + A.
            !   With A = valeValePsi, alpha = currentPop, and x,y being the
            !   wave function coefficients of the current and next kpoint
            !   respectively. (The **H means to take complex conjugate.)
            do j = 1, numStates

               ! Get the current average electron population between the two
               !   kpoint states.
               currentPopulation(:) = &
                     & (structuredElectronPopulation(j,kcur,:) + &
                     & structuredElectronPopulation(j,knxt,:)) / 2.0_double

               ! If this is not enough contribution, then skip.
               if (sum(abs(currentPopulation(:))) < smallThresh) cycle

               ! The the rank-1 update.
               do k = 1, spin
                  call zgerc(valeDim,valeDim,currentPopulation(k),&
                        & valeVale(:,j,k),1,CK(:,j,k),1,valeValePsi(:,:,k),&
                        & valeDim)
               enddo
            enddo

            do j = 1, spin
               valeValePsi(:,:,j) = matmul(valeValePsi(:,:,j),&
                     & valeValeKO(:,:,axis))
            enddo

            do j = 1, spin
               call phase_from_matrix(valeValePsi(:,:,j), ld_step)
               phiString(axis,j) = phiString(axis,j) * ld_step
            enddo

         end do  ! i = 1,numStep

         do i = 1, spin
            psiAxis(i) = aimag(log(phiString(axis,i)))
            psi(axis,i) = psi(axis,i) + psiAxis(i)
         enddo

      end do ! loop

      do i = 1, spin
         psi(axis,i)=psi(axis,i)/Fixed2D
         fn(axis,i) = psi(axis,i)/(2.0d0*pi)
         fn(axis,i) = modulo(fn(axis,i) + 0.5d0, 1.0d0) - 0.5d0
         P(:,i) = (-eCharge/realCellVolume)* fn(axis,i)*realVectors(:,axis) 
      enddo

!         psi(axis) = psi(axis) - 2.0d0*pi*dnint( psi(axis)/(2.0d0*pi) )
   enddo ! axis

   do i = 1, spin
      PCart(:,i) = P(:,i)* (10.0d0*eCharge/(bohrRad**2))
      write(20,*) 'P [C/m^2] = ', PCart(:,i)
   enddo
  
!   write(20,'(A,3ES24.16)') 'Polarization (a.u.): ', P(1,:)
!   write(20,'(A,3ES24.16)') 'Polarization (a.u.): ', P(2,:)
!   write(20,'(A,3ES24.16)') 'Polarization (a.u.): ', P(3,:)

   deallocate(CK)
   deallocate(valeValeKO)
   deallocate(structuredElectronPopulation)
   deallocate(phiString)
   deallocate(currentPopulation)

end subroutine computeMTOPPolarization
#endif

subroutine phase_from_matrix(A, detU)

   ! Define modules to use.
   use O_Kinds

   ! Make sure no variables are accidentally declared.
   implicit none

   ! Declare passed parameters.
   complex(kind=double), intent(in)    :: A(:,:)
   complex(kind=double), intent(out)  :: detU

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

   allocate(Ac(n,n), ipiv(n), visited(n), stat=istat); if (istat/=0) return
   Ac = A
   call zgetrf(n,n,Ac,n,ipiv,info); if (info/=0) then
   deallocate(Ac,ipiv,visited); return
   end if

   do i=1,n
      if (abs(Ac(i,i)) <= tiny) then
         deallocate(Ac,ipiv,visited); return
      end if
      detU = detU * Ac(i,i)
   end do

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
   if (mod(swaps,2)==1) detU = -detU

   detU = detU / cmplx(max(abs(detU), tiny), 0.0d0, kind=double)  ! unit phasor
!   theta_step = theta_step * detU

   deallocate(Ac, ipiv, visited)
end subroutine phase_from_matrix

end module O_MTOP
