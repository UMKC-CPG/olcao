!! This entire module is a collection of subroutines which were written to
!! enable the parallelization of the OLCAO method (mostly OLCAOsetup).
!! It was written during a naive era for the author who was trying to 
!! blitzkreig the parallelization of a summer. In that light, code reuse and
!! best practices are lacking. However, everything was written with the
!! mind set of having a simple equation in which to load balance various
!! places of the code. It resulted in a couple different implenetations of 
!! the same thing, which could be cleaned up to one or two more simplistic
!! things. But we all know that research code is by and large produced on
!! a results basis, and this is where it eventaully ended up
!!
!! Author: James E. Currie
!! Email: jecyrd@mail.umkc.edu
module O_ParallelSubs
  use MPI
  implicit none
!  Public
!  integer :: mpiRank
!  integer :: mpiSize
!  integer :: mpierr

  contains

!  subroutine getMPIvars(mpiRank,mpiSize)
!    implicit none
!    integer, intent(inout) :: mpiRank,mpiSize
!    call MPI_COMM_RANK(MPI_COMM_WORLD,mpiRank,mpierr)
!    call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSize,mpierr)
!  end subroutine getMPIvars

  ! This subroutine is used to balance an loop for use with MPI.  The input
  ! (toBalance) is the number of things that needs to be split up. The
  ! output (initialVal, finalVal) are the start and stop of array indices.
  ! The subroutine currently works in the way that if a loop is not evenly
  ! distributable then the extra elements are added first to the n process,
  ! then to the n-1 process, then to the n-2 process, and so on.
subroutine loadBalMPI(toBalance, initialVal, finalVal,numProcs)
  implicit none
  
  integer, intent(in) :: toBalance,numProcs
  integer :: jobsPer, remainder
  integer, intent(out) :: initialVal, finalVal
  integer :: mpiRank, mpiSize, mpierr

  call MPI_COMM_RANK(MPI_COMM_WORLD, mpiRank, mpierr)
!  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpiSize, mpierr)
  mpiSize=numProcs

  jobsPer = int(toBalance / mpiSize)
  remainder = mod(toBalance,mpisize)

  initialVal = (jobsPer * mpiRank) + 1
  finalVal = (jobsPer * (mpiRank+1))

  if (mpiRank > (mpiSize - remainder)) then
    initialVal = initialVal + (remainder - (mpiSize - mpiRank))
    finalVal = finalVal + (remainder - (mpiSize - (mpiRank+1)))
  endif
  if (mpiRank == (mpiSize - remainder)) then
    finalVal = finalVal + 1
  endif

!    This is only needed for C type array indices. i.e. first index=0
!    initialVal = initialVal - 1
!    finalVal = finalVal - 1

end subroutine loadBalMPI

subroutine elecBalance(toBalance,initialVal,finalVal,numProcs,numPotSites)
  implicit none
  
  integer, intent(in) :: toBalance,numProcs,numPotSites
  integer :: jobsPer, remainder
  integer, intent(out) :: initialVal, finalVal
  integer :: actMpiRank, mpiSize, mpierr
  integer :: adjMpiRank

  call MPI_COMM_RANK(MPI_COMM_WORLD, actMpiRank, mpierr)
!  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpiSize, mpierr)
  adjMpiRank = (actMpiRank)-numPotSites
  mpiSize=numProcs

  jobsPer = int(toBalance / mpiSize)
  remainder = mod(toBalance,mpisize)

  initialVal = (jobsPer * adjMpiRank) + 1
  finalVal = (jobsPer * (adjMpiRank+1))

  if (adjMpiRank > (mpiSize - remainder)) then
    initialVal = initialVal + (remainder - (mpiSize - adjMpiRank))
    finalVal = finalVal + (remainder - (mpiSize - (adjMpiRank+1)))
  endif
  if (adjMpiRank == (mpiSize - remainder)) then
    finalVal = finalVal + 1
  endif

!    This is only needed for C type array indices. i.e. first index=0
!    initialVal = initialVal - 1
!    finalVal = finalVal - 1

end subroutine elecBalance

! This subroutine find the indices of a 2x2 matrix that correspond to a
! packed matrix (as defined by blas/blacs/lapack). This has only been
! tested for NxN matrices. The subroutine takes the index of a
! packed matrix, and outputs the upper triangle (x,y) coordinates for 
! a unpacked matrix.
subroutine findUnpackedIndices(packedIndex, x, y)
  implicit none

  integer, intent(in) :: packedIndex
  integer, intent(out) :: x, y
  integer :: pIndex 

  pIndex = packedIndex
  pIndex = pIndex - 1
  y = int((-1+ int(sqrt(8.0*real(pIndex)+1.0))) / 2)
  x = pIndex - y*(y+1)/2
  x = x + 1
  y = y + 1

  ! note for future, to go the reverse way
  ! packedIndex = x + (y*(y+1))/2
end subroutine findUnpackedIndices

! This subroutine does the reverse of the one above.
subroutine findPackedIndex(packedIndex, x, y)
  implicit none

  integer, intent(out) :: packedIndex
  integer, intent(in) :: x,y
  
  packedIndex = x + (y*(y+1))/2
  
end subroutine findPackedIndex

! Purpose: Create the global arrays necessary for the extremely naive
! implementation of the neutral and nuclear potential section of OLCAOsetup
subroutine gaSetupNeutAndNucQPot(ga_nlNeut,ga_locNeut,ga_nlNuc,ga_locNuc,ga_potAlphaOL,potDim)
  implicit none
#include "mafdecls.fh"
#include "global.fh"
  integer, intent(inout) :: ga_nlNeut,ga_locNeut,ga_nlNuc,ga_locNuc
  integer, intent(inout) :: ga_potAlphaOL
  integer, intent(in) :: potDim
  integer :: mpiSize,mpierr
  integer :: gastat

  call MPI_Comm_Size(MPI_COMM_WORLD,mpiSize,mpierr)

  ga_nlNeut = ga_create_handle()
  ga_locNeut = ga_create_handle()
  ga_nlNuc = ga_create_handle()
  ga_locNuc = ga_create_handle()
  ga_potAlphaOL = ga_create_handle()

  call ga_set_data(ga_nlNeut,2,(/potDim,potDim/),MT_F_DBL)
  call ga_set_data(ga_locNeut,2,(/potDim,potDim/),MT_F_DBL)
  call ga_set_data(ga_nlNuc,1,(/potDim/),MT_F_DBL)
  call ga_set_data(ga_locNuc,1,(/potDim/),MT_F_DBL)
  call ga_set_data(ga_potAlphaOL,2,(/potDim,potDim/),MT_F_DBL)

  call ga_set_array_name(ga_nlNeut,"nonLocalNeutQPot")
  call ga_set_array_name(ga_locNeut,"localNeutQPot")
  call ga_set_array_name(ga_nlNuc,"nonLocalNucQPot")
  call ga_set_array_name(ga_locNuc,"localNucQPot")
  call ga_set_array_name(ga_potAlphaOL,"potAlphaOverlap")

  gastat = ga_allocate(ga_nlNeut)
  gastat = ga_allocate(ga_locNeut)
  gastat = ga_allocate(ga_nlNuc)
  gastat = ga_allocate(ga_locNuc)
  gastat = ga_allocate(ga_potAlphaOL)

  call ga_fill(ga_nlNeut,(0.0d0))
  call ga_fill(ga_locNeut,(0.0d0))
  call ga_fill(ga_nlNuc,(0.0d0))
  call ga_fill(ga_locNuc,(0.0d0))
  call ga_fill(ga_potAlphaOL,(0.0d0))

end subroutine gaSetupneutAndNucQpot

! Purpose: Destroy (deallocate) the global arrays necessary for the 
! extremely naive implementation of the neutral and nuclear potential
! section of OLCAOsetup
subroutine destroyNeutAndNucQPot(ga_nlNeut,ga_locNeut,ga_nlNuc,&
  & ga_locNuc,ga_potAlphaOL)

  implicit none
#include "mafdecls.fh"
#include "global.fh"
  integer, intent(inout) :: ga_nlNeut,ga_locNeut,ga_nlNuc,ga_locNuc
  integer, intent(inout) :: ga_potAlphaOL
  integer :: gastat
 
  gastat = ga_destroy(ga_nlNeut)
  gastat = ga_destroy(ga_locNeut)
  gastat = ga_destroy(ga_nlNuc)
  gastat = ga_destroy(ga_locNuc)
  gastat = ga_destroy(ga_potAlphaOL)

end subroutine destroyNeutAndNucQPot

! Purpose: To setup (allocate) the global arrays necessary to naively
! parallelize the Residual Q section of OLCAOsetup
subroutine gaSetupResidQ(ga_residQ,numPotTypes,potDim)
  implicit none
#include "mafdecls.fh"
#include "global.fh"
  integer, intent(inout) :: ga_residQ
  integer, intent(in) :: numPotTypes, potDim
  integer :: mpiSize,mpierr
  integer :: gastat

  call MPI_Comm_Size(MPI_COMM_WORLD,mpiSize,mpierr)

  ga_residQ = ga_create_handle()

  call ga_set_data(ga_residQ,2,(/numPotTypes,potDim/),MT_F_DBL)

  call ga_set_array_name(ga_residQ,"nonLocalResidualQ")

  gastat = ga_allocate(ga_residQ)

  call ga_fill(ga_residQ,(0.0d0))

end subroutine gaSetupResidQ

! Purpose: To destroy (de-allocate) the global arrays necessary to naively
! parallelize the Residual Q section of OLCAOsetup
subroutine destroyResidQ(ga_residQ)
  implicit none
#include "mafdecls.fh"
#include "global.fh"
  integer, intent(inout) :: ga_residQ
  integer :: gastat
 
  gastat = ga_destroy(ga_residQ)

end subroutine destroyResidQ

! Purpose: To setup (allocate) the global arrays necessary to naively
! parallelize the exchange correlation section of OLCAOsetup
subroutine gaSetupExchCorr(ga_eCOL,ga_radWght,ga_exchRhoOP, & 
    & ga_numRayPoints, maxNumRayPoints,numPotSites,potDim,exchCorrOpCode)
  implicit none
#include "mafdecls.fh"
#include "global.fh"
  integer, intent(inout) :: ga_eCOL,ga_radWght,ga_exchRhoOP,ga_numRayPoints
  integer, intent(in) :: potDim,maxNumRayPoints,numPotSites,exchCorrOpCode
  integer :: gaStat
  integer :: mpiSize, mpierr

!  call MPI_Comm_size(MPI_COMM_WORLD,mpiSize,mpierr)

  ga_eCOL = ga_create_handle()
  ga_radWght = ga_create_handle()
  ga_exchRhoOP = ga_create_handle()
  ga_numRayPoints = ga_create_handle()

  call ga_set_data(ga_eCOL,2,(/potDim,potDim/),MT_F_DBL)
  call ga_set_data(ga_radWght,2,(/maxNumRayPoints,numPotSites/),MT_F_DBL)
  call ga_set_data(ga_exchRhoOP,4,(/potDim,maxNumRayPoints, &
      & exchCorrOpCode,numPotSites/),MT_F_DBL)
  call ga_set_data(ga_numRayPoints,1,(/numPotSites/),MT_F_DBL)

  call ga_set_array_name(ga_eCOL,"exchCorrOverlap")
  call ga_set_array_name(ga_radWght,"radialWeight")
  call ga_set_array_name(ga_exchRhoOP,"exchRhoOP")
  call ga_set_array_name(ga_numRayPoints,"numRayPoints")

  gastat = ga_allocate(ga_eCOL)
  gastat = ga_allocate(ga_radWght)
  gastat = ga_allocate(ga_exchRhoOP)
  gastat = ga_allocate(ga_numRayPoints)

  call ga_fill(ga_eCOL,(0.0d0))
  call ga_fill(ga_radWght,(0.0d0))
  call ga_fill(ga_exchRhoOP,(0.0d0))
  call ga_fill(ga_numRayPoints,(0.0d0))

end subroutine gaSetupExchCorr

! Purpose: To destroy (de-allocate) the global arrays necessary to naively
! parallelize the exchange correlation section of OLCAOsetup
subroutine destroyExchCorr(ga_eCOL,ga_radWght,ga_exchRhoOP,ga_numRayPoints)
  implicit none
#include "mafdecls.fh"
#include "global.fh"
  integer, intent(inout) :: ga_eCOL, ga_radWght, ga_exchRhoOP
  integer, intent(inout) :: ga_numRayPoints
  integer :: gastat

  gastat = ga_destroy(ga_eCOL)
  gastat = ga_destroy(ga_radWght)
  gastat = ga_destroy(ga_exchRhoOP)
  gastat = ga_destroy(ga_numRayPoints)

end subroutine destroyExchCorr

! The next couple of sections are to either setup (allocate) or destroy
! (de-allocate) the global arrays necessary for the integrals sections of 
! OLCAOsetup. The overlap matrix was done seperately from the rest because
! it is required to persist through the different sets of integral
! calculations.
subroutine gaSetupOL(cvOL,coreDim,valeDim,numKPoints)
  implicit none

! include GA and MA headers
#include "mafdecls.fh"
#include "global.fh"

  integer, intent(inout),dimension(:) :: cvOL
  integer, intent(in) :: coreDim,valeDim,numKPoints
  logical :: gastat

  integer :: mpiSize, mpierr
  integer :: valeBlockDim,coreBlockDim
  integer,dimension(2):: dims
  integer :: i
  call MPI_Comm_Size(MPI_COMM_WORLD,mpiSize,mpierr)

  do i=1,numKPoints
    cvOL(i) = ga_create_handle()
    call ga_set_data(cvOL(i),2,(/coreDim,valeDim/),MT_F_DCPL)
    
!    if (valeDim/mpiSize<1) then
!      valeBlockDim=valeDim
!    elseif(valeDim/mpiSize>1 .and. mod(valeDim,mpiSize)/=0) then
!      valeBlockDim=(valeDim/mpiSize) + 1
!    else
!      valeBlockDim=valeDim/mpiSize
!    endif
!    
!    if (coreDim/mpiSize<1) then
!      coreBlockDim=coreDim
!    elseif(coreDim/mpiSize>1 .and. mod(coreDim,mpiSize)/=0) then
!      coreBlockDim=(coreDim/mpiSize) + 1
!    else
!      coreBlockDim=coreDim/mpiSize
!    endif
!  
!    dims(1)=coreBlockDim
!    dims(2)=valeBlockDim
!    call ga_set_block_cyclic(cvOL(i),dims) 
  
    call ga_set_array_name(cvOL(i), "coreValeOL")
    gastat = ga_allocate(cvOL(i))
    call ga_fill(cvOL(i),(0.0d0,0.0d0))
  enddo
end subroutine gaSetupOL

! See above gaSetupOL subroutine
subroutine gaSetup(cc,vc,cv,vv,coreDim,valeDim,numKPoints)
  implicit none

! include GA and MA headers
#include "mafdecls.fh"
#include "global.fh"

  integer, intent(inout), dimension(:) :: cc,cv,vv!,cvOL
  integer, intent(inout) :: vc
  integer, intent(in) :: coreDim,valeDim,numKPoints
  logical :: gastat

  integer :: mpiSize, mpierr
  integer :: coreBlockDim,valeBlockDim
  integer, dimension(2) :: dims
  integer :: i

  call MPI_Comm_size(MPI_COMM_WORLD,mpiSize,mpierr)
!  if (valeDim/mpiSize<1) then
!    valeBlockDim=valeDim
!  elseif(valeDim/mpiSize>1 .and. mod(valeDim,mpiSize)/=0) then
!    valeBlockDim=(valeDim/mpiSize) + 1
!  else
!    valeBlockDim=valeDim/mpiSize
!  endif
!  
!  if (coreDim/mpiSize<1) then
!    coreBlockDim=coreDim
!  elseif(coreDim/mpiSize>1 .and. mod(coreDim,mpiSize)/=0) then
!    coreBlockDim=(coreDim/mpiSize) + 1
!  else
!    coreBlockDim=coreDim/mpiSize
!  endif

  vc = ga_create_handle()
  call ga_set_array_name(vc, "valeCore")
  call ga_set_data(vc, 2, (/valeDim,coreDim/),MT_F_DCPL)
!  dims(1)=valeBlockDim
!  dims(2)=coreBlockDim
!  call ga_set_block_cyclic(vc, dims)
  
  gastat = ga_allocate(vc)
  call ga_fill(vc,(0.0d0,0.0d0))

  do i=1,numKPoints
    cc(i) = ga_create_handle()
    cv(i) = ga_create_handle()
    vv(i) = ga_create_handle()
  
    call ga_set_array_name(cc(i), "coreCore")
    call ga_set_array_name(cv(i), "coreVale")
    call ga_set_array_name(vv(i), "valeVale")
  
    call ga_set_data(cc(i), 2, (/coreDim,coreDim/),MT_F_DCPL)
    call ga_set_data(cv(i), 2, (/coreDim,valeDim/),MT_F_DCPL)
    call ga_set_data(vv(i), 2, (/valeDim,valeDim/),MT_F_DCPL)
  
!    dims(1)=coreBlockDim
!    dims(2)=coreBlockDim
!    call ga_set_block_cyclic(cc(i), dims)
!   
!    dims(1)=coreBlockDim
!    dims(2)=valeBlockDim
!    call ga_set_block_cyclic(cv(i), dims)
!  
!    dims(1)=valeBlockDim
!    dims(2)=valeBlockDim
!    call ga_set_block_cyclic(vv(i), dims)
  

    gastat = ga_allocate(cc(i))
    gastat = ga_allocate(cv(i))
    gastat = ga_allocate(vv(i))
    
    call ga_fill(cc(i),(0.0d0,0.0d0))
    call ga_fill(cv(i),(0.0d0,0.0d0))
    call ga_fill(vv(i),(0.0d0,0.0d0))
  enddo

end subroutine gaSetup

! See above gaSetupOL subroutine
subroutine gaDestroy(cc,vc,cv,vv)
  implicit none

! include GA and MA headers
#include "mafdecls.fh"
#include "global.fh"

  integer :: gastat,i

  integer, intent(inout),dimension(:) :: cc,cv,vv
  integer, intent(inout) :: vc

  gastat = ga_destroy(vc)
  do i=1, size(cc,1)
    gastat = ga_destroy(cc(i))
    gastat = ga_destroy(cv(i))
    gastat = ga_destroy(vv(i))
  enddo
end subroutine gaDestroy

! See above gaSetupOL subroutine
subroutine gaDestroyOL(cvOL)
  implicit none

! include GA and MA headers
#include "mafdecls.fh"
#include "global.fh"

  integer :: gastat

  integer, intent(inout),dimension(:) :: cvOL
  integer :: i

  do i=1, size(cvOL,1)
    gastat = ga_destroy(cvOL(i))
  enddo

end subroutine gaDestroyOL

! This subroutine writes the interaction matrices to disk for the case
! that no orthogonaliation is done.
!!!!!!The subroutine needs to be edited so that it'll function after the
!!!!!! saveCurrentPair subroutine is written. !!!!!!
!subroutine writeValeValeNoCore(currentPair,hslabStart,hslabCount,writeOpCode,&
!  & opCode)
!  use HDF5
!  use MPI
!  use O_SetupIntegralsHDF5
!  use O_Kinds
!  use O_AtomicSites
!  use O_Potential
!  use O_PotTypes
!  use O_KPoints
!  use O_Orthogonalization
!  use O_ParallelSubs
!
!  implicit none
!#include "mafdecls.fh"
!#include "global.fh"
!
!  integer :: hdferr
!  integer :: x,y,z
!  integer(hsize_t), dimension(2) :: inverseStart
!  integer(hsize_t), dimension(2) :: inverseCount
!  integer(hsize_t), dimension(2) :: unionSlabCount
!  integer(hsize_t), dimension(2) :: unionSlabStart
!  integer(hid_t) :: memspace_dsid, memspaceComplex_dsid
!  integer(hid_t), dimension(numKPoints,potDim) :: datasetToWrite_did
!  
!  complex (kind=double), intent(in), dimension(:,:,:) :: currentPair
!  real    (kind=double), allocatable, dimension(:,:,:) :: diagonalPair
!!  real    (kind=double), dimension(size(currentPair,2),&
!!    & size(currentPair,1)) :: cmplxCurPair
!
!  integer, intent(in) :: writeOpCode
!  integer, intent(in) :: opCode
!  integer(hsize_t),intent(in), dimension(2) :: hslabStart
!  integer(hsize_t),intent(inout), dimension(2) :: hslabCount
! 
! 
!  ! Set the count buffer for the top half
!  hslabCount(1) = size(currentPair,1)
!  hslabCount(2) = size(currentPair,2)
!
!  ! We need a memory space ID to write this slab correctly
!  call h5screate_simple_f(2,hslabCount,memspace_dsid,hdferr)
!  
!  select case (opCode)
!  case (1)
!    datasetToWrite_did(:,1) = atomOverlap_did(:)
!  case (2)
!    datasetToWrite_did(:,1) = atomKEOverlap_did(:)
!  case (3)
!    datasetToWrite_did(:,1) = atomNucOverlap_did(:)
!  case (4)
!    datasetToWrite_did(:,:) = atomPotOverlap_did(:,:)
!  case default
!    print *, "wrong opCode passed to writeCurrentPair" 
!    stop
!  end select
!
!  if (writeOpCode == 1) then
!    allocate(diagonalPair (size(currentPair,1),size(currentPair,2),&
!      & numKPoints))
!
!    select case (opCode)
!
!    case (1,2,3)
!    ! need to place the complex parts of currentPair into the 
!    ! bottom half of the matrix that is to be written
!      do y=1,numKPoints 
!        do x=1,size(currentPair,2)
!          diagonalPair(x,x:,1) = real(currentPair(x,x:,1))
!          diagonalPair(x+1:,x,1) = aimag(currentPair(x+1:,x,1))
!        enddo
!      enddo
!
!    ! define the hyperslab to be written too
!      call h5sselect_hyperslab_f(valeVale_dsid,H5S_SELECT_SET_F,&
!        & hslabStart, hslabCount, hdferr)
!
!    ! Write the Data
!      do x=1, numKPoints
!        call h5dwrite_f(datasetToWrite_did(x,1), H5T_NATIVE_DOUBLE, &
!          & diagonalPair(:,:,x), hslabCount, hdferr, &
!            & file_space_id=valeVale_dsid, mem_space_id=memspace_dsid, &
!              & xfer_prp=valeVale_xferpid)
!      enddo
!    case(4)
!
!    ! need to place the complex parts of currentPair into the 
!    ! bottom half of the matrix that is to be written
!      do y=1,numKPoints 
!        do x=1,size(currentPair,2)
!          diagonalPair(x,x:,1) = real(currentPair(x,x:,1))
!          diagonalPair(x+1:,x,1) = aimag(currentPair(x+1:,x,1))
!        enddo
!      enddo
!
!    ! define the hyperslab to be written too
!      call h5sselect_hyperslab_f(valeVale_dsid,H5S_SELECT_SET_F,&
!        & hslabStart, hslabCount, hdferr)
!
!    ! Write the Data
!      do x=1, numKPoints
!        call h5dwrite_f(datasetToWrite_did(x, &
!          & potTypes(currPotTypeNumber)%cumulAlphaSum+currAlphaNumber)&
!            & , H5T_NATIVE_DOUBLE, diagonalPair(:,:,x), hslabCount,&
!              & hdferr,file_space_id=valeVale_dsid,&
!                & mem_space_id=memspace_dsid, xfer_prp=valeVale_xferpid)
!      enddo
!
!
!    case default
!      print *, "Incorrect opcode passed to writeCurrentPair i=j case"
!      stop
!    end select
!
!    deallocate(diagonalPair)
!  else
!
!    ! Set dimensions and starting point for the complex matrix
!    inverseStart(1) = hslabStart(2)
!    inverseStart(2) = hslabStart(1)
!    inverseCount(1) = hslabCount(2)
!    inverseCount(2) = hslabCount(1)
!   
!    ! need to define a dataspace for the complex matrix
!    call h5screate_simple_f(2,hslabCount,memspaceComplex_dsid,hdferr)
!
!    select case (opCode)
!    
!    case (1,2,3)
!      do x=1, numKPoints
!      ! Select slab and write the real values to the top half
!        call h5sselect_hyperslab_f(valeVale_dsid,H5S_SELECT_SET_F,&
!          & hslabStart, hslabCount, hdferr)
!        call h5dwrite_f(datasetToWrite_did(x,1),H5T_NATIVE_DOUBLE, &
!          & real(currentPair(:,:,x)),hslabCount,hdferr, &
!            & file_space_id=valeVale_dsid, mem_space_id=memspace_dsid, &
!              & xfer_prp=valeVale_xferpid)
!    
!      ! need to put the complex parts in a matrix of opposite dimensions
!      ! for writing into the bottom half of the matrix
!!        cmplxCurPair(:,:) = aimag(transpose(currentPair(:,:,x)))
! 
!      ! Select slab and write the complex values to the bottom half
!        call h5sselect_hyperslab_f(valeVale_dsid,H5S_SELECT_SET_F,&
!          & inverseStart, inverseCount, hdferr)
!        call h5dwrite_f(datasetToWrite_did(x,1),H5T_NATIVE_DOUBLE, &
!          & aimag(transpose(currentPair(:,:,x))),inverseCount,hdferr, &
!            & file_space_id=valeVale_dsid,&
!              & mem_space_id=memspaceComplex_dsid,xfer_prp=valeVale_xferpid)
!    enddo
!    case (4)
!      do x=1, numKPoints
!      ! Select slab and write the real values to the top half
!        call h5sselect_hyperslab_f(valeVale_dsid,H5S_SELECT_SET_F,&
!          & hslabStart, hslabCount, hdferr)
!        call h5dwrite_f(datasetToWrite_did(x,1),H5T_NATIVE_DOUBLE, &
!          & real(currentPair(:,:,x)),hslabCount,hdferr, &
!            & file_space_id=valeVale_dsid, mem_space_id=memspace_dsid, &
!              & xfer_prp=valeVale_xferpid)
!    
!      ! need to put the complex parts in a matrix of opposite dimensions
!      ! for writing into the bottom half of the matrix
!!        cmplxCurPair(:,:) = aimag(transpose(currentPair(:,:,x)))
! 
!      ! Select slab and write the complex values to the bottom half
!        call h5sselect_hyperslab_f(valeVale_dsid,H5S_SELECT_SET_F,&
!          & inverseStart, inverseCount, hdferr)
!        call h5dwrite_f(datasetToWrite_did(x,& 
!            & potTypes(currPotTypeNumber)%cumulAlphaSum+currAlphaNumber),& 
!              & H5T_NATIVE_DOUBLE, aimag(transpose(currentPair(:,:,x))),& 
!                & inverseCount,hdferr, file_space_id=valeVale_dsid, &
!                  & mem_space_id=memspaceComplex_dsid,&
!                    & xfer_prp=valeVale_xferpid)
!      enddo
!    case default
!      print *, "wrong opCode passed to writeCurrentPair i/=j case"
!      stop
!    end select
!    call h5sclose_f(memspaceComplex_dsid,hdferr)
!
!  endif
!
!  call h5sclose_f(memspace_dsid,hdferr)
!end subroutine writeValeValeNoCore

! Purpose: This subroutine undoes some of the load balancing that is done
! previously. It would be extremely easy to combine this with the two other
! load balancing routines earlier in this code. But due to time constraints
! for myself, (James) I'm taking the previouosly mentioned route of 
! "this works and it was easy to implement really quick"
subroutine writeResolve(toBalance, initialVal, finalVal,numProcs,tmpRank,valeDim)
  implicit none
  
  integer, intent(in) :: toBalance,numProcs,tmpRank
  integer, intent(in) :: valeDim
  integer :: jobsPer, remainder
  integer, intent(out) :: initialVal, finalVal
  integer :: mpiRank, mpiSize, mpierr

!  call MPI_COMM_RANK(MPI_COMM_WORLD, mpiRank, mpierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpiSize, mpierr)
!  mpiRank=tmpRank
!  print *, 'VALEDIM',valeDim
  jobsPer = int(toBalance / numProcs)
  remainder = mod(toBalance,mpisize)

  initialVal = (jobsPer * tmpRank) + 1
  finalVal = (jobsPer * (tmpRank+1))
  if (remainder>0) then
    if (tmpRank < (mpiSize-1)) then
      initialVal = initialVal + tmpRank
      finalVal = finalVal + tmpRank + 1
    endif
    if (tmpRank==(mpiSize-1)) then
      initialVal = initialVal + tmpRank
      finalVal = valeDim
    endif
  endif
!  if (tmpRank == (mpiSize-1) .and. remainder>0) then
!    initialVal= initialVal+1
!    finalVal = valeDim
!  endif

!    This is only needed for C type array indices. i.e. first index=0
!    initialVal = initialVal - 1
!    finalVal = finalVal - 1

end subroutine writeResolve

! This subroutine writes the valeVale matrix to disk, for the case that the
! orthogonaliztion is done. I realize it is kind of confusing and maybe
! not so straight forward so I will do my best to explain it.
! What happens is that for a small block of the matrix, we pull it down
! from a global array, to a local one on the single process execcuting
! the write to HDF5 sections from integralSCF. It then writes it to disk,
! then moves on to another small block and writes that to disk.
! The reason this was done instead of PHDF5 was because of
! NFS (Network File System) writes being inneficcient under the PHDF5
! paradigm as of HDF5 version 1.8.12.
subroutine writeValeVale(ga_valeVale,opCode,numKPoints,potDim, & 
    & currPotTypeNumber,currAlphaNumber,valeDim)
  use HDF5
  use O_Kinds
  use O_Constants
  use O_SetupIntegralsHDF5
  use O_PotTypes
  
  implicit none
#include "mafdecls.fh"
#include "global.fh"

  integer, intent(in),dimension(:) :: ga_valeVale
  integer, intent(in) :: opCode,valeDim
  integer, intent(in) :: currPotTypeNumber,currAlphaNumber
  integer, intent(in) :: numKPoints, potDim
  integer, dimension(2) :: hi, lo, blockDims, numBlocks
  integer :: mpiRank, mpiSize, mpierr, hdferr
  complex (kind=double), allocatable, dimension(:,:) :: valeValeGA
  real    (kind=double), allocatable, dimension(:,:) :: diskVV
  integer :: i,j,k,m,x,y, xDimCnt, yDimCnt
  integer(hid_t) :: memspace_dsid
  integer(hid_t), dimension(numKPoints,potDim) :: datasetToWrite_did
  integer(hsize_t), dimension(2) :: hslabCount,hslabStart,dims
  integer :: valeBlockDim
  integer :: minDim,maxDim,tmpRank
  ! Define small threshold for eliminating resultant values
  real (kind=double) :: smallThresh10

  smallThresh10 = real(1.0d-10,double)

  call MPI_Comm_rank(MPI_COMM_WORLD, mpiRank, mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD, mpiSize, mpierr)

  select case (opCode)
  case (1)
    datasetToWrite_did(:,1) = atomOverlap_did(:)
  case (2)
    datasetToWrite_did(:,1) = atomKEOverlap_did(:)
  case (3)
    datasetToWrite_did(:,1) = atomNucOverlap_did(:)
  case (4)
    datasetToWrite_did(:,:) = atomPotOverlap_did(:,:)
  case default
    print *, "wrong opCode passed to writeValeVale"
    stop
  end select
  
  if (valeDim/mpiSize<1) then
    valeBlockDim=valeDim
  elseif(valeDim/mpiSize>1 .and. mod(valeDim,mpiSize)/=0) then
    valeBlockDim=(valeDim/mpiSize) + 1
  else
    valeBlockDim=valeDim/mpiSize
  endif
  blockDims(1) = valeBlockDim
  blockDims(2) = valeBlockDim
  if (mod(valeDim,blockDims(1)) > 0) then
    numBlocks(1) = int(valeDim/blockDims(1)) + 1
  else
    numBlocks(1) = int(valeDim/blockDims(1))
  endif
  numBlocks(2) = numBlocks(1)
!  print *, "blockDims: ", blockDims
!  print *, "numBlocks: ", numBlocks(1)
!  print *, "Divide:    ", valeDim/blockDims(1)
!  print *, "mod:       ", mod(valeDim,blockDims(1))

  do i=1, numKPoints
    x=0
    y=0
    do m=1,numBlocks(1)**2
      if (x /=  numBlocks(1)-1 .and. y /= numBlocks(2)-1) then
        lo(1) = ((x)*blockDims(1))+1
        lo(2) = ((y)*blockDims(2))+1
        hi(1) = (x+1)*blockDims(1)
        hi(2) = (y+1)*blockDims(2)
      elseif (x == numBlocks(1)-1 .and. y /= numBlocks(2)-1) then
        lo(1) = valeDim-blockDims(1)+1
        lo(2) = ((y)*blockDims(2))+1
        hi(1) = valeDim
        hi(2) = (y+1)*blockDims(2)
      elseif (x /= numBlocks(1)-1 .and. y == numBlocks(2)-1) then
        lo(1) = ((x)*blockDims(1))+1
        lo(2) = valeDim-blockDims(2)+1
        hi(1) = (x+1)*blockDims(1)
        hi(2) = valeDim
      else
        lo(1) = valeDim-blockDims(1)+1
        lo(2) = valeDim-blockDims(2)+1
        hi(1) = valeDim
        hi(2) = valeDim
      endif
      
      allocate(valeValeGA(hi(1)-lo(1)+1, hi(2)-lo(2)+1))
      allocate(diskVV(hi(1)-lo(1)+1, hi(2)-lo(2)+1))
      call nga_get(ga_valeVale(i),lo,hi,valeValeGA, size(valeValeGA,1))
       
      yDimCnt=lo(2)
      do k=1,size(valeValeGA,2)
        xDimCnt=lo(1)
        do j=1,size(valeValeGA,1)
          if (yDimCnt<xDimCnt) then
            diskVV(j,k) = aimag(conjg(valeValeGA(j,k)))
          else
            diskVV(j,k) = real(valeValeGA(j,k))
          endif
          xDimCnt = xDimCnt + 1
        enddo
        yDimCnt = yDimCnt + 1
      enddo
      
      do k=1,size(diskVV,2)
        do j=1,size(diskVV,1)
         if(abs(diskVV(j,k))<smallThresh10) then
            diskVV(j,k) = 0.0_double
         endif
        enddo
      enddo

      hslabStart(1) = lo(1)-1
      hslabStart(2) = lo(2)-1
      
      hslabCount(1)=(hi(1)-lo(1))+1
      hslabCount(2)=(hi(2)-lo(2))+1
      call h5screate_simple_f(2,hslabCount,memspace_dsid,hdferr)
        
      ! define the hyperslab to be written to
      call h5sselect_hyperslab_f(valeVale_dsid,H5S_SELECT_SET_F, &
        & hslabStart,hslabCount,hdferr) 
      select case (opCode)
      case(1:3)
        ! Write slabs to disk
        call h5dwrite_f(datasetToWrite_did(i,1),H5T_NATIVE_DOUBLE, &
          & diskVV(:,:), hslabCount, hdferr, &
            & file_space_id=valeVale_dsid, mem_space_id=memspace_dsid)
      case(4)
        ! Write slabs to disk
        call h5dwrite_f(datasetToWrite_did(i, &
        & potTypes(currPotTypeNumber)%cumulAlphaSum+currAlphaNumber), &
        & H5T_NATIVE_DOUBLE, diskVV(:,:), hslabCount, &
        & hdferr, file_space_id=valeVale_dsid, &
        & mem_space_id=memspace_dsid)
      case default
        print *, "shits fucked up yo"
      end select

      x=x+1
      if (x==numBlocks(1)) then
        y = y +1
        x = 0
      endif

      deallocate(valeValeGA)
      deallocate(diskVV)
      call h5sclose_f(memspace_dsid,hdferr)
    enddo
  enddo

end subroutine writeValeVale

! This subroutine is for OLCAOMain to correctly read in the valeVale
! Matrix. Obviously it's not done yet.
subroutine readValeValeFromDisk()

end subroutine readValeValeFromDisk

end module O_ParallelSubs
