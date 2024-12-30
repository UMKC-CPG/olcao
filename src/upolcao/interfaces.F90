module O_LAPACKParameters

   ! Make sure nothing funny is declared.
   implicit none

   ! The optimal blocksize; if this value is 1, an unblocked
   !   algorithm will give the best performance
   integer :: blockSize

   contains

subroutine setBlockSize(arrayDim)

   ! Make sure nothing funny is declared.
   implicit none

   ! Define passed parameters.
   integer :: arrayDim

   interface
      integer function ilaenv (ISPEC,NAME,OPTS,N1,N2,N3,N4)
         use O_Kinds
         integer       :: ISPEC
         character*(*) :: NAME
         character*(*) :: OPTS
         integer       :: N1
         integer       :: N2
         integer       :: N3
         integer       :: N4
      end function ilaenv
   end interface

#ifndef GAMMA
   blockSize = (ilaenv (1,'zhetrd','U',arrayDim,arrayDim,-1,-1)+1) * arrayDim
#else
   blockSize = (ilaenv (1,'dsytrd','U',arrayDim,arrayDim,-1,-1)+2) * arrayDim
#endif


end subroutine setBlockSize

end module O_LAPACKParameters




module O_LAPACKZHEGV

   ! Make sure nothing funny is declared by accident.
   implicit none

   interface
      subroutine zhegv (ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,RWORK,INFO)
         use O_Kinds      
         integer     :: LWORK
         integer     :: ITYPE   
         character*1 :: JOBZ
         character*1 :: UPLO
         integer     :: N
         integer     :: LDA
         integer     :: LDB
         complex     (kind=double), dimension (LDA,N) :: A
         complex     (kind=double), dimension (LDB,N) :: B
         real        (kind=double), dimension (N) :: W
         complex     (kind=double), dimension (LWORK) :: WORK
         real        (kind=double), dimension (3*N-2) :: RWORK
         integer     :: INFO
      end subroutine zhegv
   end interface

   contains

subroutine solveZHEGV(N,resultDim,A,B,eigenValues)

   ! Use necessary modules.
   use O_Kinds
   use O_LAPACKParameters, only: blockSize

   ! Make sure no funny variables are defined.
   implicit none

   ! Define passed parameters.
   integer, intent(IN) :: N
   integer, intent(IN) :: resultDim
   complex (kind=double), dimension (N,N), intent(INOUT) :: A
   complex (kind=double), dimension (N,N), intent(IN)    :: B
   real    (kind=double), dimension (resultDim) , intent(OUT)   :: eigenValues

   ! Define local variables.
   integer :: lwork !LAPACK
   integer :: info  !LAPACK
   real    (kind=double), allocatable, dimension (:)   :: rwork     !LAPACK
   complex (kind=double), allocatable, dimension (:)   :: work      !LAPACK
   real    (kind=double), allocatable, dimension (:)   :: tempLambda

   ! Get the optimal block size.
   lwork = blockSize

   ! Allocate space for temporary work
   allocate (work (lwork))
   allocate (rwork (3 * N - 2))
   allocate (tempLambda (N))


   call zhegv (1,'V','U',N,A(:,:),N,B(:,:),N,tempLambda,work,lwork,rwork,info)

   if (info /= 0) then
      write (20,*) 'ZHEGV error code=',info
      stop
   endif

   ! Copy the result from the temp array.
   eigenValues(:resultDim) = tempLambda(:resultDim)

   ! Deallocate unnecessary arrays and matrices.
   deallocate (tempLambda)
   deallocate (work)
   deallocate (rwork)

end subroutine solveZHEGV

end module O_LAPACKZHEGV




module O_LAPACKDSYGV

   ! Make sure nothing funny is declared by accident.
   implicit none

   interface
      subroutine dsygv (ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,INFO)
         use O_Kinds
         integer     :: LWORK
         integer     :: ITYPE
         character*1 :: JOBZ
         character*1 :: UPLO
         integer     :: N
         integer     :: LDA
         integer     :: LDB
         real        (kind=double), dimension (LDA,N) :: A
         real        (kind=double), dimension (LDB,N) :: B
         real        (kind=double), dimension (N) :: W
         real        (kind=double), dimension (LWORK) :: WORK
         integer     :: INFO
      end subroutine dsygv
   end interface

   contains

subroutine solveDSYGV(valeDim,numStates,valeValeGamma,valeValeOLGamma,&
      & eigenValues)

   ! Use necessary modules.
   use O_Kinds
   use O_LAPACKParameters, only: blockSize

   ! Make sure no funny variables are defined.
   implicit none

   ! Define passed parameters.
   integer :: valeDim
   integer :: numStates
   real (kind=double), dimension (valeDim,valeDim), intent(INOUT) :: &
         & valeValeGamma
   real (kind=double), dimension (valeDim,valeDim), intent(IN)    :: &
         & valeValeOLGamma
   real (kind=double), dimension (valeDim)        , intent(OUT)   :: &
         & eigenValues

   ! Define local variables.
   integer :: lwork !LAPACK
   integer :: info  !LAPACK
   real    (kind=double), allocatable, dimension (:)   :: workGamma !LAPACK
   real    (kind=double), allocatable, dimension (:)   :: tempEigenValues

   ! Get the optimal block size.
   lwork = blockSize

   ! Allocate space for temporary work.
   allocate (workGamma (lwork))
   allocate (tempEigenValues (valeDim))

   call dsygv (1,'V','U',valeDim,valeValeGamma(:,:),valeDim,&
         & valeValeOLGamma(:,:),valeDim,tempEigenValues,workGamma,lwork,info)

   if (info /= 0) then
      write (20,*) 'DSYGV error code=',info
      stop
   endif

   ! Copy the result from the temp array.
   eigenValues(:numStates) = tempEigenValues(:numStates)

   ! Deallocate unnecessary arrays and matrices.
   deallocate (workGamma)
   deallocate (tempEigenValues)

end subroutine solveDSYGV

end module O_LAPACKDSYGV



module O_LAPACKDPOSVX
   interface
      subroutine dposvx(FACT,UPLO,N,NRHS,A,LDA,AF,LDAF,EQUED,S,B,LDB,X,&
         & LDX,RCOND,FERR,BERR,WORK,IWORK,INFO)
         use O_Kinds
         character*1 :: FACT
         character*1 :: UPLO
         integer :: N
         integer :: NRHS
         integer :: LDA
         integer :: LDAF
         integer :: LDB
         integer :: LDX
         real (kind=double), dimension (LDA,N) :: A
         real (kind=double), dimension (LDAF,N) :: AF
         character*1 :: EQUED
         real (kind=double), dimension (N) :: S
         real (kind=double), dimension (LDB,NRHS) :: B
         real (kind=double), dimension (LDX,NRHS) :: X
         real (kind=double) :: RCOND
         real (kind=double), dimension (NRHS) :: FERR
         real (kind=double), dimension (NRHS) :: BERR
         real (kind=double), dimension (3*N) :: WORK
         integer, dimension (N) :: IWORK
         integer :: INFO
      end subroutine dposvx
   end interface

   contains

subroutine solveDPOSVX (N, NRHS, A, LD, B, INFO)

   ! Use necessary modules.
   use O_Kinds

   ! Make sure no funny variables are defined.
   implicit none

   ! Define passed parameters.
   integer, intent(IN) :: N
   integer, intent(IN) :: NRHS
   integer, intent(IN) :: LD
   integer, intent(OUT) :: INFO
   real (kind=double), dimension (LD,N), intent(IN)  :: A
   real (kind=double), dimension (LD,NRHS), intent(INOUT)  :: B

   ! Define local variables.
   character*1 :: equilibrated
   real (kind=double) :: conditionNumber
   real (kind=double), allocatable, dimension (:) :: scaleFactor
   real (kind=double), allocatable, dimension (:) :: forwardError
   real (kind=double), allocatable, dimension (:) :: backwardError
   real (kind=double), allocatable, dimension (:) :: realWorkSpace
   real (kind=double), allocatable, dimension (:,:) :: resultMatrix
   real (kind=double), allocatable, dimension (:,:) :: factor
   integer, allocatable, dimension (:) :: integerWorkSpace

   ! Allocate space to hold the working arrays.
   allocate (scaleFactor(N))
   allocate (forwardError(NRHS))
   allocate (backwardError(NRHS))
   allocate (realWorkSpace(3*N))
   allocate (resultMatrix(N,NRHS))
   allocate (factor (LD,LD))
   allocate (integerWorkSpace(N))

   call dposvx('E','U',N,NRHS,A,LD,factor,LD,equilibrated,scaleFactor,B,LD,&
         & resultMatrix,LD,conditionNumber,forwardError,backwardError,&
         & realWorkSpace,integerWorkSpace,info)

!   if (info /= 0) then
!      write (20, *) 'dposvx failed. INFO= ', info
!      stop
!   endif

   ! Recover the solution into the given B matrix.
   B(:,:) = resultMatrix(:,:)

   ! Deallocate space for the working arrays.
   deallocate (scaleFactor)
   deallocate (forwardError)
   deallocate (backwardError)
   deallocate (realWorkSpace)
   deallocate (resultMatrix)
   deallocate (factor)
   deallocate (integerWorkSpace)

end subroutine solveDPOSVX

end module O_LAPACKDPOSVX



module O_BLASZHER
   interface
      subroutine zher(UPLO,N,ALPHA,X,INCX,A,LDA)
         use O_Kinds
         character*1 :: UPLO
         integer :: N
         integer :: LDA
         real (kind=double) :: ALPHA
         complex (kind=double), dimension (N) :: X
         integer :: INCX
         complex (kind=double), dimension (LDA,N) :: A
      end subroutine zher
   end interface
end module O_BLASZHER

module O_BLASDSYR
   interface
      subroutine dsyr(UPLO,N,ALPHA,X,INCX,A,LDA)
         use O_Kinds
         character*1 :: UPLO
         integer :: N
         integer :: LDA
         real (kind=double) :: ALPHA
         real (kind=double), dimension (N) :: X
         integer :: INCX
         real (kind=double), dimension (LDA,N) :: A
      end subroutine dsyr
   end interface
end module O_BLASDSYR

module zherkInterface
   interface
      subroutine zherk (UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
         use O_Kinds
         character*1 :: UPLO
         character*1 :: TRANS
         integer :: N
         integer :: K
         integer :: LDA
         integer :: LDC
         real (kind=double) :: ALPHA
         complex (kind=double), dimension(LDA,N) :: A
         real (kind=double) :: BETA
         complex (kind=double), dimension(LDC,N) :: C
      end subroutine zherk
   end interface
end module zherkInterface

module zher2kInterface
   interface
      subroutine zher2k (UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
         use O_Kinds
         character*1 :: UPLO
         character*1 :: TRANS
         integer :: N
         integer :: K
         integer :: LDA
         integer :: LDB
         integer :: LDC
         complex (kind=double) :: ALPHA
         complex (kind=double), dimension(LDA,N) :: A
         complex (kind=double), dimension(LDB,N) :: B
         real (kind=double) :: BETA
         complex (kind=double), dimension(LDC,N) :: C
      end subroutine zher2k
   end interface
end module zher2kInterface

module dsyrkInterface
   interface
      subroutine dsyrk (UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
         use O_Kinds
         character*1 :: UPLO
         character*1 :: TRANS
         integer :: N
         integer :: K
         integer :: LDA
         integer :: LDC
         real (kind=double) :: ALPHA
         real (kind=double), dimension(LDA,N) :: A
         real (kind=double) :: BETA
         real (kind=double), dimension(LDC,N) :: C
      end subroutine dsyrk
   end interface
end module dsyrkInterface

module dsyr2kInterface
   interface
      subroutine dsyr2k (UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
         use O_Kinds
         character*1 :: UPLO
         character*1 :: TRANS
         integer :: N
         integer :: K
         integer :: LDA
         integer :: LDB
         integer :: LDC
         real (kind=double) :: ALPHA
         real (kind=double), dimension(LDA,N) :: A
         real (kind=double), dimension(LDB,N) :: B
         real (kind=double) :: BETA
         real (kind=double), dimension(LDC,N) :: C
      end subroutine dsyr2k
   end interface
end module dsyr2kInterface

