module O_IntgSaving

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine multWithBasisFn1 (currentBasisFns,pairXBasisFn2,pairXBasisFn12,&
      & currentlmIndex,currentNumTotalStates,maxAlpha1Used)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules
   use O_AtomicTypes

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define variables passed to this subroutine
   real (kind=double), dimension (maxNumAtomAlphas,&
         & maxNumStates,2) :: currentBasisFns
   real (kind=double), dimension (16,maxNumAtomAlphas,&
         & maxNumStates) :: pairXBasisFn2
   real (kind=double), dimension (maxNumStates,maxNumStates) :: pairXBasisFn12
   integer, dimension (maxNumStates,2)  :: currentlmIndex
   integer, dimension(2) :: currentNumTotalStates
   integer :: maxAlpha1Used

   ! Define local variables
   integer :: l,m

   do l = 1, currentNumTotalStates(2)
      do m = 1, currentNumTotalStates(1)
         pairXBasisFn12(m,l) = &
               & sum(currentBasisFns(:maxAlpha1Used,m,1) * &
               & pairXBasisFn2(currentlmIndex(m,1),:maxAlpha1Used,l))
      enddo
   enddo

end subroutine multWithBasisFn1


#ifndef GAMMA


subroutine applyPhaseFactors (currentPair,pairXBasisFn12,statesDim1,statesDim2,&
      & k,runCode,currentKPoint)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules
   use O_KPoints
   use O_AtomicTypes

   ! Make sure no funny variables are defined accidentally.
   implicit none

   ! Define the variables passed to this subroutine
   integer :: statesDim1
   integer :: statesDim2
   complex (kind=double), dimension (maxNumStates,maxNumStates,&
        & numKPoints) :: currentPair
   real (kind=double), dimension (statesDim1,statesDim2) :: pairXBasisFn12
   integer :: k ! Cell loop index
   integer :: runCode
   integer :: currentKPoint  ! Only used for runCode>0

   ! Define the local variables
   integer :: l,m,n

   if (runCode == 0) then
      do l = 1, numKPoints
         do m = 1, statesDim2
            do n = 1, statesDim1
               currentPair(n,m,l) = currentPair(n,m,l) + &
                     & phaseFactor(l,k) * pairXBasisFn12(n,m)
            enddo
         enddo
      enddo
   elseif (runCode <= 2) then
      do m = 1, statesDim2
         do n = 1, statesDim1
            currentPair(n,m,1) = currentPair(n,m,1) + &
                  & phaseFactor(currentKPoint,k) * pairXBasisFn12(n,m)
         enddo
      enddo
   else  ! Include a -i factor.
      do m = 1, statesDim2
         do n = 1, statesDim1
            currentPair(n,m,1) = currentPair(n,m,1) + cmplx(0.0_double,&
                  & -1.0_double) * phaseFactor(currentKPoint,k) * &
                  & pairXBasisFn12(n,m)
         enddo
      enddo
   endif

end subroutine applyPhaseFactors


subroutine kPointLatticeOriginShift (currentNumTotalStates,currentPair,&
      & latticeVector,kPointCount,kPointIndex)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules
   use O_KPoints
   use O_AtomicTypes

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define variables passed to this subroutine
   integer, dimension(2) :: currentNumTotalStates
   integer :: KPointCount
   complex (kind=double), dimension (:,:,:) :: currentPair
   real (kind=double), dimension (3) :: latticeVector
   integer :: kPointIndex

   ! Define local variables for loop control.
   integer :: i,j,k

   ! Lattice origin shift kpoint effect variables
   complex (kind=double) :: dotProduct
   real (kind=double) :: kPointShiftDot

   do i = 1, kPointCount

      ! Find the kpoint, lattice vector dot product   
      kPointShiftDot = sum(kPoints(:,kPointIndex+i) * latticeVector(:))

      ! Check that the vectors are not perpendicular.
      if (kPointShiftDot .ne. 0.0_double) then

         ! Get the cosine and sine of the dot product
         dotProduct = cmplx(cos(kPointShiftDot),sin(kPointShiftDot),double)

         ! Loop over both atom wave functions.
         do j = 1, currentNumTotalStates(2)
            do k = 1, currentNumTotalStates(1)
               currentPair(k,j,i) = dotProduct * currentPair(k,j,i)
            enddo
         enddo
      endif
   enddo

end subroutine kPointLatticeOriginShift


subroutine saveCurrentPair (i,j,kPointCount,currentPair,&
      & ga_valeVale,ga_coreVale,ga_coreCore,currentNumTotalStates)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules
   use O_KPoints
   use O_AtomicTypes
   use O_AtomicSites

   Use MPI

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the variables passed to this subroutine.
   integer, intent(in) :: i,j  ! Atom1 and Atom2
   integer, intent(in) :: kPointCount
   complex (kind=double), intent(inout), dimension(:,:,:) :: currentPair
!   complex (kind=double), dimension (valeDim,valeDim,numKPoints) :: valeVale
!   complex (kind=double), dimension (coreDim,valeDim,numKPoints) :: coreVale
!   complex (kind=double), dimension (coreDim,coreDim,numKPoints) :: coreCore

   ! Define the matrix necessary for quickly accessing the complex conjugate
   !   of the currentPair.
   complex (kind=double), allocatable, dimension (:,:,:) :: currentPairDagger
   integer :: kPointCounter


   ! Index variables for tracking the locations in the large matrices where
   !   the accumulated results should be stored.
   integer, dimension (2) :: valeStateIndex
   integer, dimension (2) :: valeStateNum
   integer, dimension (2) :: coreStateIndex
   integer, dimension (2) :: coreStateNum

   ! Other local variables
   integer, dimension (2) :: currentAtomType
   integer,intent(inout),dimension(:) :: ga_valeVale,ga_coreVale,ga_coreCore
   integer, intent(in), dimension(2) :: currentNumTotalStates
   integer :: dat1,dat2

   integer :: mpiRank,mpierr
   ! Assign the current atom types
   currentAtomType(1) = atomSites(i)%atomTypeAssn
   currentAtomType(2) = atomSites(j)%atomTypeAssn


   ! Create the complex conjugate of the relevant currentPair section.
!   allocate (currentPairDagger(currentNumTotalStates(2), &
!    & currentNumTotalStates(1),kPointCount))
   allocate(currentpairDagger(maxNumStates,maxNumStates,kPointCount))
   do kPointCounter = 1, kPointCount
      currentPairDagger(:,:,kPointCounter) = &
            & transpose(conjg(currentPair(:,:,kPointCounter)))
!      currentPairDagger(1:valeStateNum(2)+coreStateNum(2),1:valeStateNum(1) + &
!            & coreStateNum(1),i) = &
!            & transpose(conjg(currentPair(1:valeStateNum(1)+coreStateNum(1),&
!            & 1:valeStateNum(2)+coreStateNum(2),i)))
   enddo

!   print *, "Current Pair"
!  do dat1=1,(size(currentPair,1))
!    do dat2=1,(size(currentPair,2))
!      print *, dat1,dat2, currentPair(dat1,dat2,1)
!    enddo
!  enddo
!   print *, "Cmplx Trans Curr Pair"
!  do dat1=1,(size(currentPairDagger,1))
!    do dat2=1,(size(currentPairDagger,2))
!      print *, dat1,dat2, currentPairDagger(dat1,dat2,1)
!    enddo
!  enddo

   ! Get the indices for where the states for these atoms should begin to
   !   be recorded in the core containing matrices.  Also get the number
   !   of states that these atoms occupy.
   valeStateIndex(1) = atomSites(i)%cumulValeStates
   valeStateIndex(2) = atomSites(j)%cumulValeStates
   valeStateNum(1)   = atomTypes(currentAtomType(1))%numValeStates
   valeStateNum(2)   = atomTypes(currentAtomType(2))%numValeStates

   ! Save the valence valence part.  The upper triangle contains the real
   !   part while the strict lower triangle contains the imaginary part.
   call valeValeSaving (i,j,kPointCount,valeStateIndex,valeStateNum,&
         & currentPair, currentPairDagger, ga_valeVale)

   ! If there is no core contribution then we don't have to save any
   !   core parts.
   if (coreDim == 0) return

   ! Get the indices for where the states for these atoms should begin to
   !   be recorded in the core containing matrices.  Also get the number
   !   of states that these atoms occupy.
   coreStateIndex(1) = atomSites(i)%cumulCoreStates
   coreStateIndex(2) = atomSites(j)%cumulCoreStates
   coreStateNum(1)   = atomTypes(currentAtomType(1))%numCoreStates
   coreStateNum(2)   = atomTypes(currentAtomType(2))%numCoreStates

   call MPI_Comm_rank(MPI_COMM_WORLD,mpiRank,mpierr)
   ! Save the overlap from the core of atom 1 with the valence of atom 2.
   !   This code is basically a replication of the old code.
   call coreValeSaving (kPointCount,valeStateIndex,valeStateNum,&
         & coreStateIndex,coreStateNum,currentPair,currentPairDagger,&
         & ga_coreVale)

   ! Now we do the Core Core part just the same as we did for the Vale
   !   Vale part.  The only difference is the shift of valeStateNum(:)
   !   in determining the start and end points of the matrix to copy.
   if ((coreStateNum(1) == 0) .or. (coreStateNum(2) == 0)) return

   call coreCoreSaving (i,j,kPointCount,valeStateNum,coreStateIndex,&
         & coreStateNum,currentPair,currentPairDagger, ga_coreCore)

   ! Deallocate the dagger of the current pair since it is no longer needed.
   deallocate (currentPairDagger)

end subroutine saveCurrentPair



subroutine valeValeSaving (atom1,atom2,kPointCount,valeStateIndex,valeStateNum,&
      & currentPair, currentPairDagger, ga_vv)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules.
   use O_KPoints
   use O_AtomicTypes
   use O_AtomicSites

   ! Make sure that there are not accidental variable declarations.
   implicit none
#include "mafdecls.fh"
#include "global.fh"

   ! Define variables passed to this subroutine
   integer, intent(inout),dimension(:) :: ga_vv
   integer, intent(in) :: atom1
   integer, intent(in) :: atom2
   integer, intent(in) :: kPointCount
   integer, intent(in), dimension (2) :: valeStateIndex
   integer, intent(in), dimension (2) :: valeStateNum
!   complex (kind=double), dimension (valeDim,valeDim,numKPoints) :: valeVale
   complex (kind=double), intent(in), dimension(:,:,:) :: currentPair
   complex (kind=double), intent(in), dimension(:,:,:) :: currentPairDagger
   integer, dimension(2) :: startIndices
   integer, dimension(2) :: endIndices
   integer, dimension(2) :: dagStart
   integer, dimension(2) :: dagEnd
   integer, dimension(2) :: ldReal
   integer, dimension(2) :: ldCmplx

   ! Define local variables for loop control and tracking.
   integer :: i,j,k

   startIndices(1) = valeStateIndex(1)+1
   startIndices(2) = valeStateIndex(2)+1
   endIndices(1) = valeStateIndex(1)+valeStateNum(1)
   endIndices(2) = valeStateIndex(2)+valeStateNum(2)
   dagStart(1) = startIndices(2)
   dagStart(2) = startIndices(1)
   dagEnd(1) = endIndices(2)
   dagEnd(2) = endIndices(1)
   ldReal(1) = valeStateNum(1)
   ldReal(2) = valeStateNum(2)
   ldCmplx(1) = valeStateNum(2)
   ldCmplx(2) = valeStateNum(1)
!   print *, "valeStatIndex:",valeStateIndex
!   print *, "valeStateNum: ",valeStateNum
!   print *, "valevale"
!   print *, "start:        ",startIndices
!   print *, "end:          ",endIndices
!   print *, "dagstart:     ",dagStart
!   print *, "dagend:       ",dagEnd
!   flush(6)
   if (atom1 == atom2) then
      do i = 1, kPointCount
!         valeVale(valeStateIndex(1)+1:valeStateIndex(1)+valeStateNum(1),&
!               & valeStateIndex(2)+1:valeStateIndex(2)+valeStateNum(2),i) = &
!               & currentPair(1:valeStateNum(1),1:valeStateNum(2),i)

!          do j=1, valeStateNum(1)
!            do k=1,valeStateNum(2)
!              print *, j,k, currentPair(j,k,1)
!            enddo
!          enddo
          ! Plop the full matrix on the diagonal
          call ga_put(ga_vv(i), startIndices(1),endIndices(1),&
            & startIndices(2), endIndices(2), &
            & currentPair(1:valeStateNum(1),1:valeStateNum(2),i),&
            & ldReal(1))

!         do j = valeStateIndex(1)+1,valeStateIndex(1)+valeStateNum(1)
!            do k = j+1,valeStateIndex(2)+valeStateNum(2)
!               valeVale(k,j,i) = (0.0_double,0.0_double)
!            enddo
!         enddo
      enddo
   else
      do i = 1, kPointCount
!         valeVale(valeStateIndex(1)+1:valeStateIndex(1)+valeStateNum(1),&
!               & valeStateIndex(2)+1:valeStateIndex(2)+valeStateNum(2),i) = &
!               & currentPair(1:valeStateNum(1),1:valeStateNum(2),i)
         ! First save to top half
         call ga_put(ga_vv(i),startIndices(1),endIndices(1), &
          & startIndices(2),endIndices(2), &
          & currentPair(1:valeStateNum(1),1:valeStateNum(2),i),&
          & ldReal(1))
         ! Put the transpose into the lower half
         call ga_put(ga_vv(i),dagStart(1),dagEnd(1), &
          & dagStart(2),dagEnd(2), &
          & currentPairDagger(1:valeStateNum(2),1:valeStatenum(1),i),&
          & ldCmplx(1))
      enddo
   endif
!   print *, "pre curentPair valeVale"
!   do i = 1,valeStateNum(1)
!    do j = 1,valeStateNum(2)
!      print *, currentPair(i,j,1)
!    enddo
!   enddo
!   print *, "pre currentPair Dagger valeVale"
!   do i = 1, valeStateNum(2)
!    do j = 1, valeStateNum(1)
!      print *, currentPairDagger(i,j,1)
!    enddo
!   enddo
   ! Print dat shit out negro.
!   call  ga_print(ga_vv)
end subroutine valeValeSaving



subroutine coreValeSaving (kPointCount,valeStateIndex,valeStateNum,&
      & coreStateIndex,coreStateNum,currentPair,currentPairDagger,&
      & ga_cv)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules.
   use O_KPoints
   use O_AtomicTypes
   use O_AtomicSites

   ! Make sure that there are not accidental variable declarations.
   implicit none
#include "mafdecls.fh"
#include "global.fh"

   ! Define variables passed to this subroutine
   integer, intent(inout),dimension(:) :: ga_cv
   integer :: kPointCount
   integer, dimension (2) :: valeStateIndex
   integer, dimension (2) :: valeStateNum
   integer, dimension (2) :: coreStateIndex
   integer, dimension (2) :: coreStateNum
!   complex (kind=double), dimension (coreDim,valeDim,numKPoints) :: coreVale
   complex (kind=double), intent(inout), dimension(:,:,:) :: currentPair
   complex (kind=double), intent(inout), dimension(:,:,:) :: currentPairDagger
   ! Define local variables for loop control.
   integer :: i,j

   integer, dimension(2) :: startIndices
   integer, dimension(2) :: endIndices
   integer, dimension(2) :: dagStart
   integer, dimension(2) :: dagEnd
   integer, dimension(2) :: ldReal, ldCmplx
   startIndices(1) = coreStateIndex(1)+1
   startIndices(2) = valeStateIndex(2)+1

   endIndices(1) = coreStateIndex(1) + coreStateNum(1)
   endIndices(2) = valeStateIndex(2) + valeStateNum(2)

   dagStart(1) = coreStateIndex(2)+1
   dagStart(2) = valeStateIndex(1)+1

   dagEnd(1) = coreStateIndex(2)+coreStateNum(2)
   dagEnd(2) = valeStateIndex(1)+valeStateNum(1)

   ldReal(1) = endIndices(1)-startIndices(1)+1
   ldReal(2) = endIndices(2)-startIndices(2)+1

   ldCmplx(1) = dagEnd(1)-dagStart(1)+1
   ldCmplx(2) = dagEnd(2)-dagStart(2)+1


   if (coreStateNum(1) /= 0) then
      do i = 1, kPointCount
!       coreVale(coreStateIndex(1)+1:coreStateIndex(1)+coreStateNum(1),&
!             & valeStateIndex(2)+1:valeStateIndex(2)+valeStateNum(2),i) =&
!             & currentPair(valeStateNum(1)+1:valeStateNum(1)+coreStateNum(1),&
!             & 1:valeStateNum(2),i)
         call ga_put(ga_cv(i),startIndices(1),endIndices(1),&
          & startIndices(2), endIndices(2),&
          & currentPair(valeStateNum(1)+1:valeStateNum(1)+coreStateNum(1),&
          & 1:valeStateNum(2),i),ldReal(1))
      enddo
   endif
   if (coreStateNum(2) /= 0) then
      do i = 1, kPointCount
!         coreVale(coreStateIndex(2)+1:coreStateIndex(2)+coreStateNum(2),&
!               & valeStateIndex(1)+1:valeStateIndex(1)+valeStateNum(1),i) = &
!               & currentPairDagger(valeStateNum(2)+1:valeStateNum(2)+&
!               & coreStateNum(2),1:valeStateNum(1),i)
         call ga_put(ga_cv(i),dagStart(1),dagEnd(1),&
          & dagStart(2), dagEnd(2), &
          & currentPairDagger(valeStateNum(2)+1:valeStateNum(2)+&
          & coreStateNum(2),1:valeStateNum(1),i), ldCmplx(1))
      enddo
   endif
!   print *, "pre curentPair coreVale"
!   do i = 1,valeStateNum(1)
!    do j = 1,valeStateNum(2)
!      print *, currentPair(i,j,1)
!    enddo
!   enddo
!   print *, "pre currentPair Dagger coreVale"
!   do i = valeStateNum(2)+1, valeStateNum(2)+coreStateNum(2)
!    do j = 1, valeStateNum(1)
!      print *, currentPairDagger(i,j,1)
!    enddo
!   enddo
!  call ga_print(ga_cv)
end subroutine coreValeSaving


subroutine coreCoreSaving (atom1,atom2,kPointCount,valeStateNum,&
      & coreStateIndex,coreStateNum,currentPair,currentPairDagger,&
      & ga_cc)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules.
   use O_KPoints
   use O_AtomicTypes
   use O_AtomicSites

   ! Make sure that there are not accidental variable declarations.
   implicit none
#include "mafdecls.fh"
#include "global.fh"

   ! Define variables passed to this subroutine
   integer, intent(inout),dimension(:) :: ga_cc
   integer :: atom1
   integer :: atom2
   integer :: kPointCount
   integer, dimension (2) :: valeStateNum
   integer, dimension (2) :: coreStateIndex
   integer, dimension (2) :: coreStateNum
!   complex (kind=double), dimension (coreDim,coreDim,numKPoints) :: coreCore
   complex (kind=double), intent(inout), dimension(:,:,:) :: currentPair
   complex (kind=double), intent(inout), dimension(:,:,:) :: currentPairDagger
   ! Define local variables for loop control.
   integer :: i
   
   integer, dimension(2) :: startIndices
   integer, dimension(2) :: endIndices
   integer, dimension(2) :: dagStart
   integer, dimension(2) :: dagEnd
   integer, dimension(2) :: ldReal, ldCmplx
   
   startIndices(1) = coreStateIndex(1)+1
   startIndices(2) = coreStateIndex(2)+1
   endIndices(1) = coreStateIndex(1)+coreStateNum(1)
   endIndices(2) = coreStateIndex(2)+coreStateNum(2)
   dagStart(1) = coreStateIndex(2)+1
   dagStart(2) = coreStateIndex(1)+1
   dagEnd(1) = coreStateIndex(2)+coreStateNum(2)
   dagEnd(2) = coreStateIndex(1)+coreStateNum(1)

   ldReal(1) = endIndices(1)-startIndices(1)+1
   ldReal(2) = endIndices(2)-startIndices(2)+1

   ldCmplx(1) = dagEnd(1)-dagStart(1)+1
   ldCmplx(2) = dagEnd(2)-dagStart(2)+1

   do i = 1, kPointCount

!      coreCore(coreStateIndex(1)+1:coreStateIndex(1)+coreStateNum(1),&
!            & coreStateIndex(2)+1:coreStateIndex(2)+coreStateNum(2),i) = &
!            & currentPair(valeStateNum(1)+1:valeStateNum(1)+coreStateNum(1),&
!            & valeStateNum(2)+1:valeStateNum(2)+coreStateNum(2),i)

      call ga_put(ga_cc(i), startIndices(1),endIndices(1),&
        & startIndices(2), endIndices(2), &
        & currentPair(valeStateNum(1)+1:valeStateNum(1)+coreStateNum(1),&
        & valeStateNum(2)+1:valeStateNum(2)+coreStateNum(2),i),&
        & ldReal(1))
      if (atom1 .ne. atom2) then
!         coreCore(coreStateIndex(2)+1:coreStateIndex(2)+coreStateNum(2),&
!               & coreStateIndex(1)+1:coreStateIndex(1)+coreStateNum(1),i)=&
!               & currentPairDagger(valeStateNum(2)+1:valeStateNum(2)+ &
!               & coreStateNum(2),valeStateNum(1)+1:valeStateNum(1)+ &
!               & coreStateNum(1),i)
        call ga_put(ga_cc(i),dagStart(1),dagEnd(1),&
          & dagStart(2),dagEnd(2), &
          & currentPairDagger(valeStateNum(2)+1:valeStateNum(2)+&
          & coreStateNum(2),valeStateNum(1)+1:valeStateNum(1)+&
          & coreStateNum(1),i),ldCmplx(1))
      endif
   enddo
!  call ga_print(ga_cc)
end subroutine coreCoreSaving

#else

subroutine applyPhaseFactorsGamma (currentPairGamma,pairXBasisFn12,statesDim1,&
      & statesDim2,k,runCode)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules
   use O_AtomicTypes

   ! Make sure no funny variables are defined
   implicit none

   ! Define the variables passed to this subroutine
   integer :: statesDim1
   integer :: statesDim2
   real (kind=double), dimension (maxNumStates,maxNumStates) :: currentPairGamma
   real (kind=double), dimension (statesDim1,statesDim2) :: pairXBasisFn12
   integer :: k ! Cell loop index
   integer :: runCode

   if (runCode <= 2) then
      currentPairGamma(:statesDim1,:statesDim2) = &
            & currentPairGamma(:statesDim1,:statesDim2) + &
            & pairXBasisFn12(:statesDim1,:statesDim2) 
   else  ! Inlude a -i factor and store the imaginary part (real=0 now).
      currentPairGamma(:statesDim1,:statesDim2) = &
            & currentPairGamma(:statesDim1,:statesDim2) - &
            & pairXBasisFn12(:statesDim1,:statesDim2)
   endif


end subroutine applyPhaseFactorsGamma

subroutine saveCurrentPairGamma (i,j,currentPairGamma,&
      & valeValeGamma,coreValeGamma,coreCoreGamma)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules
   use O_KPoints
   use O_AtomicTypes
   use O_AtomicSites

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the variables passed to this subroutine.
   integer :: i,j ! Atom1 and Atom2
   real (kind=double), dimension (maxNumStates,maxNumStates) :: currentPairGamma
   real (kind=double), dimension (valeDim,valeDim) :: valeValeGamma
   real (kind=double), dimension (coreDim,valeDim) :: coreValeGamma
   real (kind=double), dimension (coreDim,coreDim) :: coreCoreGamma

   ! Define the matrix for quickly accessing the transpose of the currentPair.
   real (kind=double), allocatable, dimension (:,:) :: currentPairGammaTranspose

   ! Index variables for tracking the locations in the large matrices where
   !   the accumulated results should be stored.
   integer, dimension (2) :: valeStateIndex
   integer, dimension (2) :: valeStateNum
   integer, dimension (2) :: coreStateIndex
   integer, dimension (2) :: coreStateNum

   ! Other local variables
   integer, dimension (2) :: currentAtomType

   ! Assign the current atom types
   currentAtomType(1) = atomSites(i)%atomTypeAssn
   currentAtomType(2) = atomSites(j)%atomTypeAssn

   ! Create the transpose of the relevant currentPairGamma section.
   allocate (currentPairGammaTranspose(maxNumStates,maxNumStates))
   currentPairGammaTranspose = transpose(currentPairGamma)
!   currentPairGammaTranspose(1:valeStateNum(2)+coreStateNum(2),&
!         & 1:valeStateNum(1)+coreStateNum(1)) = &
!         & transpose(currentPairGamma(1:valeStateNum(1)+coreStateNum(1),&
!         & 1:valeStateNum(2)+coreStateNum(2)))

   ! Get the indices for where the states for these atoms should begin to
   !   be recorded in the core containing matrices.  Also get the number
   !   of states that these atoms occupy.
   valeStateIndex(1) = atomSites(i)%cumulValeStates
   valeStateIndex(2) = atomSites(j)%cumulValeStates
   valeStateNum(1)   = atomTypes(currentAtomType(1))%numValeStates
   valeStateNum(2)   = atomTypes(currentAtomType(2))%numValeStates

   ! Save the valence valence part.  The upper triangle contains the real
   !   part while the strict lower triangle contains the imaginary part.
   call valeValeSavingGamma (i,j,valeStateIndex,valeStateNum,&
         & valeValeGamma,currentPairGamma)

   ! If there is no core contribution then we don't have to save any
   !   core parts.
   if (coreDim == 0) return

   ! Get the indices for where the states for these atoms should begin to
   !   be recorded in the core containing matrices.  Also get the number
   !   of states that these atoms occupy.
   coreStateIndex(1) = atomSites(i)%cumulCoreStates
   coreStateIndex(2) = atomSites(j)%cumulCoreStates
   coreStateNum(1)   = atomTypes(currentAtomType(1))%numCoreStates
   coreStateNum(2)   = atomTypes(currentAtomType(2))%numCoreStates

   ! Save the overlap from the core of atom 1 with the valence of atom 2.
   !   This code is basically a replication of the old code.
   call coreValeSavingGamma (valeStateIndex,valeStateNum,&
         & coreStateIndex,coreStateNum,coreValeGamma,currentPairGamma,&
         & currentPairGammaTranspose)

   ! Now we do the Core Core part just the same as we did for the Vale
   !   Vale part.  The only difference is the shift of valeStateNum(:)
   !   in determining the start and end points of the matrix to copy.

   if ((coreStateNum(1) == 0) .or. (coreStateNum(2) == 0)) return

   call coreCoreSavingGamma (i,j,valeStateNum,coreStateIndex,&
         & coreStateNum,coreCoreGamma,currentPairGamma,&
         & currentPairGammaTranspose)

   ! Deallocate the transpose of the current pair since it is no longer needed.
   deallocate (currentPairGammaTranspose)

end subroutine saveCurrentPairGamma



subroutine valeValeSavingGamma (atom1,atom2,valeStateIndex,valeStateNum,&
      & valeValeGamma,currentPairGamma)


   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules
   use O_AtomicTypes
   use O_AtomicSites

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define variables passed to this subroutine
   integer :: atom1
   integer :: atom2
   integer, dimension (2) :: valeStateIndex
   integer, dimension (2) :: valeStateNum
   real (kind=double), dimension (valeDim,valeDim) :: valeValeGamma
   real (kind=double), dimension (maxNumStates,&
         & maxNumStates) :: currentPairGamma

   ! Define local variables for loop control and tracking.
   integer :: i,j

   if (atom1 == atom2) then
      valeValeGamma(valeStateIndex(1)+1:valeStateIndex(1)+valeStateNum(1),&
            & valeStateIndex(2)+1:valeStateIndex(2)+valeStateNum(2)) = &
            & currentPairGamma(1:valeStateNum(1),1:valeStateNum(2))
      do i = valeStateIndex(1)+1,valeStateIndex(1)+valeStateNum(1)
         do j = i+1,valeStateIndex(2)+valeStateNum(2)
            valeValeGamma(j,i) = 0.0_double
         enddo
      enddo
   else
      valeValeGamma(valeStateIndex(1)+1:valeStateIndex(1)+valeStateNum(1),&
            & valeStateIndex(2)+1:valeStateIndex(2)+valeStateNum(2)) = &
            & currentPairGamma(1:valeStateNum(1),1:valeStateNum(2))
   endif
end subroutine valeValeSavingGamma



subroutine coreValeSavingGamma (valeStateIndex,valeStateNum,&
      & coreStateIndex,coreStateNum,coreValeGamma,currentPairGamma,&
      & currentPairGammaTranspose)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules
   use O_AtomicTypes
   use O_AtomicSites

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define variables passed to this subroutine
   integer, dimension (2) :: valeStateIndex
   integer, dimension (2) :: valeStateNum
   integer, dimension (2) :: coreStateIndex
   integer, dimension (2) :: coreStateNum
   real (kind=double), dimension (coreDim,valeDim) :: coreValeGamma
   real (kind=double), dimension (maxNumStates,&
         & maxNumStates) :: currentPairGamma
   real (kind=double), dimension (maxNumStates,&
         & maxNumStates) :: currentPairGammaTranspose

   ! Define local variables for loop control.
   integer :: i

   if (coreStateNum(1) .ne. 0) then
      coreValeGamma(coreStateIndex(1)+1:coreStateIndex(1)+coreStateNum(1),&
            & valeStateIndex(2)+1:valeStateIndex(2)+valeStateNum(2)) = &
            & currentPairGamma(valeStateNum(1)+1:&
            & valeStateNum(1)+coreStateNum(1),1:valeStateNum(2))
   endif
   if (coreStateNum(2) .ne. 0) then
      coreValeGamma(coreStateIndex(2)+1:coreStateIndex(2)+coreStateNum(2),&
            & valeStateIndex(1)+1:valeStateIndex(1)+valeStateNum(1)) = &
            & currentPairGammaTranspose(valeStateNum(2)+1:valeStateNum(2) + &
            & coreStateNum(2),1:valeStateNum(1))
   endif
end subroutine coreValeSavingGamma



subroutine coreCoreSavingGamma (atom1,atom2,valeStateNum,coreStateIndex,&
      & coreStateNum,coreCoreGamma,currentPairGamma,currentPairGammaTranspose)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules
   use O_AtomicTypes
   use O_AtomicSites

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define variables passed to this subroutine
   integer :: atom1
   integer :: atom2
   integer, dimension (2) :: valeStateNum
   integer, dimension (2) :: coreStateIndex
   integer, dimension (2) :: coreStateNum
   real (kind=double), dimension (coreDim,coreDim) :: coreCoreGamma
   real (kind=double), dimension (maxNumStates,&
         & maxNumStates) :: currentPairGamma
   real (kind=double), dimension (maxNumStates,&
         & maxNumStates) :: currentPairGammaTranspose


   coreCoreGamma(coreStateIndex(1)+1:coreStateIndex(1)+coreStateNum(1),&
         & coreStateIndex(2)+1:coreStateIndex(2)+coreStateNum(2)) = &
         & currentPairGamma(valeStateNum(1)+1:valeStateNum(1)+coreStateNum(1),&
         & valeStateNum(2)+1:valeStateNum(2)+coreStateNum(2))
   if (atom1 .ne. atom2) then
      coreCoreGamma(coreStateIndex(2)+1:coreStateIndex(2)+coreStateNum(2), &
            & coreStateIndex(1)+1:coreStateIndex(1)+coreStateNum(1)) = &
            & currentPairGammaTranspose(valeStateNum(2)+1:valeStateNum(2) + &
            & coreStateNum(2),valeStateNum(1)+1:valeStateNum(1)+coreStateNum(1))
   endif
end subroutine coreCoreSavingGamma

#endif

end module O_IntgSaving
