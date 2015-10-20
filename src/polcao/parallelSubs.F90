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
!  integer :: mpiErr

  contains

!  subroutine getMPIvars(mpiRank,mpiSize)
!    implicit none
!    integer, intent(inout) :: mpiRank,mpiSize
!    call MPI_COMM_RANK(MPI_COMM_WORLD,mpiRank,mpiErr)
!    call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSize,mpiErr)
!  end subroutine getMPIvars

  ! This subroutine is used to balance an loop for use with MPI.  The input
  ! (toBalance) is the number of things that needs to be split up. The
  ! output (initialVal, finalVal) are the start and stop of array indices.
  ! The subroutine currently works in the way that if a loop is not evenly
  ! distributable then the extra elements are added first to the n process,
  ! then to the n-1 process, then to the n-2 process, and so on.
subroutine loadBalMPI(toBalance,initialVal,finalVal,myProc,numProcs)
  implicit none

  ! Define passed parameters.
  integer, intent(in) :: toBalance,numProcs,myProc
  integer, intent(out) :: initialVal, finalVal

  ! Define local variables.
  integer :: jobsPer, remainder
  integer :: mpiRank, mpiSize, mpiErr

  mpiSize=numProcs
  mpiRank=myProc

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
  integer :: actMpiRank, mpiSize, mpiErr
  integer :: adjMpiRank

  call MPI_COMM_RANK(MPI_COMM_WORLD, actMpiRank, mpiErr)
!  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpiSize, mpiErr)
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
  integer :: mpiRank, mpiSize, mpiErr

!  call MPI_COMM_RANK(MPI_COMM_WORLD, mpiRank, mpiErr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpiSize, mpiErr)
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
  integer :: mpiRank, mpiSize, mpiErr, hdferr
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

  call MPI_Comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)
  call MPI_Comm_size(MPI_COMM_WORLD, mpiSize, mpiErr)

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
