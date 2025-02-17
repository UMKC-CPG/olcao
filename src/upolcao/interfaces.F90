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

   ! Note that "zhetrd" is used per the documentation for zhegv. (Same for
   !   "dsytrd" and dsygv.)
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
   use O_MPI
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
      if (mpiRank == 0) then
         write (20,*) 'ZHEGV error code=',info
      endif
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


! Reference for this module: 
!    http://www.netlib.org/scalapack/explore-html/d7/dff/pzhegvx_8f_source.html
module O_SCALAPACKPZHEGVX

   ! Make sure nothing funny is declared by accident.
   implicit none

   interface
      subroutine pzhegvx(IBTYPE,JOBZ,RANGE,UPLO,N,A,IA,JA,DESCA,B,IB,JB,DESCB,&
            & VL,VU,IL,UI,ABSTOL,M,NZ,W,ORFAC,Z,IZ,JZ,DESCZ,WORK,LWORK,RWORK,&
            & LRWORK,IWORK,LIWORK,IFAIL,ICLUSTR,GAP,INFO)
         use O_KINDS
         integer   :: IBTYPE
         character :: JOBZ
         character :: RANGE
         character :: UPLO
         integer   :: N
         complex(kind=double), dimension(:,:) :: A
         integer   :: IA
         integer   :: JA
         integer, dimension(9) :: DESCA
         complex(kind=double), dimension(:,:) :: B
         integer   :: IB
         integer   :: JB
         integer, dimension(9) :: DESCB
         real(kind=double) :: VL
         real(kind=double) :: VU
         integer   :: IL
         integer   :: IU
         real(kind=double) :: ABSTOL
         integer   :: M
         integer   :: NZ
         real(kind=double), dimension(N) :: W
         real(kind=double) :: ORFAC
         complex(kind=double), dimension(:,:) :: Z
         integer   :: IZ
         integer   :: JZ
         integer, dimension(9) :: DESCZ
         complex(kind=double) :: WORK
         integer   :: LWORK
         real(kind=double), dimension(:) :: RWORK
         integer   :: LRWORK
         integer, dimension(:) :: IWORK
         integer   :: LIWORK
         integer, dimension(N) :: IFAIL
         integer, dimension(:) :: ICLUSTR
         real(kind=double), dimension(:) :: GAP
         integer   :: INFO
      end subroutine pzhegvx
   end interface

   contains

subroutine solvePZHEGVX(vvArr,vvOLArr,eVals,blcsinfo)

   use O_Kinds
   use O_MPI

   ! Make sure no funny variables are defined
   implicit none

   ! Define passed Parameters
   type(ArrayInfo), intent(inout) :: vvArr
   type(ArrayInfo), intent(inout) :: vvOLArr
   type(ArrayInfo), intent(inout) :: eVals
   type(BlacsInfo), intent(in) :: blcsinfo

   ! Define local input variables
   real(kind=double) :: VL,VU
   integer :: IL, IU
   real(kind=double) :: ABSTOL
   real(kind=double) :: ORFAC

   ! Define local output variables
   integer :: M, NZ
   real(kind=double), dimension(vvInfo%I) :: W
   complex(kind=double), dimension(vvInfo%I,vvInfo%I) :: Z

   ! Define local work variables
   complex(kind=double) :: WORK
   integer   :: LWORK
   real(kind=double), allocatable, dimension(:) :: RWORK
   integer   :: LRWORK
   integer, allocatable, dimension(:) :: IWORK
   integer   :: LIWORK

   ! Define other local variables
   integer, dimension(vvInfo%I) :: IFAIL
   integer, dimension(2*blcsinfo%prows*blcsinfo%pcols) :: ICLUSTR
   real(kind=double), dimension(blcsinfo%prows*blcsinfo%pcols) :: GAP
   integer   :: INFO

   ! Define external functions
   real(kind=double), external :: PDLAMCH
   integer, external :: numroc

   ! These are not referenced if RNGE='A'
   VL = 0
   VU = 0
   IL = 0
   IU = 0

   ! Eigvenvalues will be computed most accurately when ABSTOL is set to twice
   !   the underflow threshold 2*PDLAMCH('S') not zero.
   ABSTOL = 2*PDLAMCH(blcsinfo%context,'S')

   ! Global outputs, initializing to 0
   M = 0
   NZ = 0

   ! Allocate space for eigenvectors and intitialize
   allocate(W(vvInfo%I))
   W(:) = 0.0_double

   ! Specifies which eigenvectors should be reorthonalized. Eigenvectors
   ! that correspond to eigenvalues which are within tol=ORFAC*norm(A) of 
   ! each other are to be reorthogonalized.
   ! I think this should be set to 10^-3, but we'll use 0 for now, and no
   ! eigenvectors will be reorthogonalized.
   ORFAC = 0

   ! We will first call PZHEGVX with LWORK, LRWORK, and LIWORK specified as
   ! -1. this will cause scalapack to do a workspace query and output the
   ! optimal sizes of these parameters in the first index of the associated
   ! arrays. We'll then reallocate to the correct sizes.
   LWORK = -1
   allocate(WORK(1))
   WORK(:) = 0

   LRWORK = -1
   ALLOCATE(RWORK(1))
   RWORK(:) = 0

   LIWORK = 01
   allocate(IWORK(1))
   IWORK(:) = 0

   ! Other outputs initializing to 0
   IFAIL(:) = 0
   ICLUSTER(:) = 0
   GAP(:) = 0.0_double
   INFO = 0

   ! As stated above we call the PZHEGVX subroutine as a workspace query
   call pzhegvx(1,'V','A','U',vvInfo%I,vvArr%local,1,1,vvArr%desc, &
         & vvOLArr%local,1,1,vvOLArr%desc, VL,VU,IL,IU,ABSTOL,M,NZ,W,ORFAC, &
         & eVals%local,1,1,eVals%desc, WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK, &
         & IFAIL, ICLUSTR, GAP, INFO)

   ! Now we set the proper workspace parameters as needed, and resize.
   LWORK = WORK(1)
   LRWORK = RWORK(1)
   LIWORK = IWORK(1)

   deallocate(WORK)
   deallocate(RWORK)
   deallocate(IWORK)
   allocate(WORK(LWORK))
   allocate(RWORK(LRWORK))
   allocate(IWORK(LIWORK))
   WORK(:) = complex(0.0_double,0.0_double)
   RWORK(:) = 0.0_double
   IWORK(:) = 0

   ! Now we have sufficient workspace for scalapack and can now call PZHEGVX
   ! to actually solve our eigen problem.
   call pzhegvx(1,'V','A','U',vvInfo%I,vvArr%local,1,1,vvArr%desc, &
         & vvOLArr%local,1,1,vvOLArr%desc, VL,VU,IL,IU,ABSTOL,M,NZ,W,ORFAC, &
         & eVals%local,1,1,eVals%desc, WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK, &
         & IFAIL, ICLUSTR, GAP, INFO)

   deallocate(WORK)
   deallocate(RWORK)
   deallocate(IWORK)



end subroutine solvePZHEGVX

end module O_SCALAPACKPZHEGVX



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
   use O_MPI
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
      if (mpiRank == 0) then
         write (20,*) 'DSYGV error code=',info
      endif
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
   use O_MPI
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

   if (info /= 0) then
      if (mpiRank == 0) then
         write (20, *) 'dposvx failed. INFO= ', info
      endif
      stop
   endif

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

module pztrancInterface
   interface
      subroutine pztranc(M,N,ALPHA,A,IA,JA,DESCA,BETA,C,IC,JC,DESCC)
         use O_Kinds
         integer:: M
         integer:: N
         complex(kind=double) :: ALPHA
         complex(kind=double), dimension(N,M) :: A
         integer :: IA
         integer :: JA
         integer, dimension(9) :: DESCA
         complex(kind=double) :: BETA
         complex(kind=double), dimension(M,N) :: C
         integer :: IC
         integer :: JC
         integer, dimension(9) :: DESCC
      end subroutine pztranc
   end interface
end module pztrancInterface

module pdtranInterface
   interface
      subroutine pdtran(M,N,ALPHA,A,IA,JA,DESCA,BETA,C,IC,JC,DESCC)
         use O_Kinds
         integer:: M
         integer:: N
         real(kind=double) :: ALPHA
         real(kind=double), dimension(N,M) :: A
         integer :: IA
         integer :: JA
         integer, dimension(9) :: DESCA
         real(kind=double) :: BETA
         real(kind=double), dimension(M,N) :: C
         integer :: IC
         integer :: JC
         integer, dimension(9) :: DESCC
      end subroutine pdtran
   end interface
end module pdtranInterface

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

module pzherkInterface
   interface
      subroutine pzherk (UPLO,TRANS,N,K,ALPHA,A,IA,JA,DESCA,BETA,C,IC,JC,DESCC)
         use O_Kinds
         character*1 :: UPLO
         character*1 :: TRANS
         integer :: N
         integer :: K
         complex (kind=double) :: ALPHA
         complex (kind=double), dimension(K,N) :: A
         integer :: IA
         integer :: JA
         integer, dimension(9) :: DESCA
         complex (kind=double) :: BETA
         complex (kind=double), dimension(N,N) :: C
         integer :: IC
         integer :: JC
         integer, dimension(9) :: DESCC
      end subroutine pzherk
   end interface
end module pzherkInterface

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

module pzher2kInterface
   interface
      subroutine pzher2k (UPLO,TRANS,N,K,ALPHA,A,IA,JA,DESCA, &
            & B,IB,JB,DESCB,BETA,C,IC,JC,DESCC)
         use O_Kinds
         character*1 :: UPLO
         character*1 :: TRANS
         integer :: N
         integer :: K
         complex (kind=double) :: ALPHA
         complex (kind=double), dimension(K,N) :: A
         integer :: IA
         integer :: JA
         integer, dimension(9) :: DESCA
         complex (kind=double), dimension(K,N) :: B
         integer :: IB
         integer :: JB
         integer, dimension(9) :: DESCB
         complex (kind=double) :: BETA
         complex (kind=double), dimension(N,N) :: C
         integer :: IC
         integer :: JC
         integer, dimension(9) :: DESCC
      end subroutine pzher2k
   end interface
end module pzher2kInterface

module pzhemmInterface
   interface
      subroutine pzhemm (SIDE,UPLO,M,N,ALPHA,A,IA,JA,DESCA, &
            & B,IB,JB,DESCB,BETA,C,IC,JC,DESCC)
         use O_Kinds
         character*1 :: SIDE
         character*1 :: UPLO
         integer :: M
         integer :: N
         complex (kind=double) :: ALPHA
         complex (kind=double), dimension(N,N) :: A
         integer :: IA
         integer :: JA
         integer, dimension(9) :: DESCA
         complex (kind=double), dimension(M,N) :: B
         integer :: IB
         integer :: JB
         integer, dimension(9) :: DESCB
         complex (kind=double) :: BETA
         complex (kind=double), dimension(M,N) :: C
         integer :: IC
         integer :: JC
         integer, dimension(9) :: DESCC
      end subroutine pzhemm
   end interface
end module pzhemmInterface

module pzgemmInterface
   interface
      subroutine pzgemm (TRANSA,TRANSB,M,N,K,ALPHA,A,IA,JA,DESCA, &
            & B,IB,JB,DESCB,BETA,C,IC,JC,DESCC)
         use O_Kinds
         character*1 :: TRANSA
         character*1 :: TRANSB
         integer :: M
         integer :: N
         integer :: K
         complex (kind=double) :: ALPHA
         complex (kind=double), dimension(M,K) :: A
         integer :: IA
         integer :: JA
         integer, dimension(9) :: DESCA
         complex (kind=double), dimension(K,N) :: B
         integer :: IB
         integer :: JB
         integer, dimension(9) :: DESCB
         complex (kind=double) :: BETA
         complex (kind=double), dimension(M,N) :: C
         integer :: IC
         integer :: JC
         integer, dimension(9) :: DESCC
      end subroutine pzgemm
   end interface
end module pzgemmInterface

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

module pdsyr2kInterface
   interface
      subroutine pdsyr2k (UPLO,TRANS,N,K,ALPHA,A,IA,JA,DESCA,&
            & B,IB,JB,DESCB,BETA,C,IC,JC,DESCC)
         use O_Kinds
         character*1 :: UPLO
         character*1 :: TRANS
         integer :: N
         integer :: K
         real (kind=double) :: ALPHA
         real (kind=double), dimension(K,N) :: A
         integer :: IA
         integer :: JA
         integer, dimension(9) :: DESCA
         real (kind=double), dimension(K,N) :: B
         integer :: IB
         integer :: JB
         integer, dimension(9) :: DESCB
         real (kind=double) :: BETA
         real (kind=double), dimension(N,N) :: C
         integer :: IC
         integer :: JC
         integer, dimension(9) :: DESCC
      end subroutine pdsyr2k
   end interface
end module pdsyr2kInterface

module pdsymmInterface
   interface
      subroutine pdsymm (SIDE,UPLO,M,N,ALPHA,A,IA,JA,DESCA,&
            & B,IB,JB,DESCB,BETA,C,IC,JC,DESCC)
         use O_Kinds
         character*1 :: SIDE
         character*1 :: UPLO
         integer :: M
         integer :: N
         real (kind=double) :: ALPHA
         real (kind=double), dimension(N,N) :: A
         integer :: IA
         integer :: JA
         integer, dimension(9) :: DESCA
         real (kind=double), dimension(M,N) :: B
         integer :: IB
         integer :: JB
         integer, dimension(9) :: DESCB
         real (kind=double) :: BETA
         real (kind=double), dimension(M,N) :: C
         integer :: IC
         integer :: JC
         integer, dimension(9) :: DESCC
      end subroutine pdsymm
   end interface
end module pdsymmInterface

module pdgemmInterface
   interface
      subroutine pdgemm (TRANSA,TRANSB,M,N,K,ALPHA,A,IA,JA,DESCA,&
            & B,IB,JB,DESCB,BETA,C,IC,JC,DESCC)
         use O_Kinds
         character*1 :: TRANSA
         character*1 :: TRANSB
         integer :: M
         integer :: N
         integer :: K
         real (kind=double) :: ALPHA
         real (kind=double), dimension(M,K) :: A
         integer :: IA
         integer :: JA
         integer, dimension(9) :: DESCA
         real (kind=double), dimension(K,N) :: B
         integer :: IB
         integer :: JB
         integer, dimension(9) :: DESCB
         real (kind=double) :: BETA
         real (kind=double), dimension(M,N) :: C
         integer :: IC
         integer :: JC
         integer, dimension(9) :: DESCC
      end subroutine pdgemm
   end interface
end module pdgemmInterface
