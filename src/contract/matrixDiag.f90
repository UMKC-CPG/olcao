module DiagonalizationSubs

   private
   public :: diagonalizeSPDF

   contains

subroutine diagonalizeSPDF

   ! Import definition modules.
   use O_Kinds

   ! Import data modules.
   use AtomData
   use MatrixElementData
   use WaveFunctionData
   use GaussianBasisData

   implicit none

   ! Define local variables.
   integer :: i
   integer :: lwork
   integer :: info
   integer :: currDim
   real (kind=double), allocatable, dimension (:) :: work
   character*1, dimension (4) :: stateNames

   ! Define the lapack diagonalization interface.
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

   ! Initialize state names for printing.
   stateNames(1) = 's'
   stateNames(2) = 'p'
   stateNames(3) = 'd'
   stateNames(4) = 'f'

   ! Allocate space to hold the diagonalization results.
   allocate (eigenValues(maxNumBasisGaussians,4))

   ! Initialize all eigen values to zero.
   eigenValues(:,:) = 0.0_double

   ! Diagonalize each orbital type's block (s,p,d,f).
   do i = 1,4

      ! Skip this iteration if there are no orbitals of this type.
      if (numTotalOrbitals(i) == 0) cycle

      ! Make a short name for the current matrix dimension.
      currDim=numBasisGaussians(i)

      ! Initialize parameters for temporary work space.
      lwork = (ilaenv (1,'dsytrd','U',currDim,currDim,-1,-1)+2) * currDim

      ! Allocate space for temporary work space.
      allocate (work(lwork))

      ! Perform the diagonalization.  The eigenvectors are stored in the
      !   hamiltonianME matrix with the first index being the number of basis
      !   gaussians, the second index being the number of orbitals of the
      !   current angular momentum QN, and the third index being the angular
      !   momentum quantum number.
      call dsygv (1,'V','U',currDim,hamiltonianME(:currDim,:currDim,i),&
            & currDim,overlapME(:currDim,:currDim,i),currDim,&
            & eigenValues(:currDim,i),work,lwork,info)

      ! Abort if the diagonalization failed.
      if (info /= 0) then
         write (20,*) "Diagonalization failed with info=",info," i=",i
         stop
      endif

      ! Write out the eigenvalue results to standard output.
      write (20,*) ' Eigenvalues for ',stateNames(i),' states.'
      write (20,fmt="(10f12.6)") eigenValues(:numTotalOrbitals(i),i)

      ! Deallocate unnecessary work space.
      deallocate (work)
   enddo

   ! Correct the positions of the eigenvector values.  Recall that often the
   !   d or f orbitals will have fewer gaussians than the s and p orbitals.
   !   If (as can often be done) the first gaussian of the basis is not
   !   included in the d or f expansion we need to make sure that the
   !   indices of the hamiltonian match the indices of the alpha array.
   !   (e.g. The hamiltonian contains 3 terms corresponding to the last 3
   !   (of 4 say) gaussian alphas.  The hamiltonian needs to be shifted to
   !   have 4 terms with the first one being 0.  This way ham index 1 ->
   !   alpha index 1, and ham index 2 -> alpha index 2...)
   call correctHamiltonianIndex

end subroutine diagonalizeSPDF

subroutine correctHamiltonianIndex

   ! Import definition modules.
   use O_Kinds

   ! Import data modules.
   use GaussianBasisData
   use MatrixElementData

   implicit none

   ! Define passed parameters.
   

   ! Define local variables.
   integer :: i,j
   integer :: selectGaussianCount

   ! Consider each orbital type (spdf).
   do i = 1, 4

      ! Initialize the counter for the number of gaussian functions to be used
      !   by this orbital ang mom type.  (Will count down.)
      selectGaussianCount = numBasisGaussians(i)

      ! It is necessary to start at the back and work forward.  Recall that the
      !   hamiltonianME contains only the eigenvector coefficients for the
      !   non-zero selectGaussians and that they are all consecutive starting
      !   from the beginning of the first dimension index.  If
      !   numBasisGaussians(i) = maxNumBasisGaussians then this is pointless
      !   and redundent.  If numBasisGaussians(i) < maxNumBasisGaussians then
      !   the hamiltonianME gets filled from the end towards the front of the
      !   first dimension index using 0s when necessary and copying the coeffs
      !   when necessary.  Nothing gets overwritten.
      do j = maxNumBasisGaussians, 1, -1

         if (selectGaussians(j,i) == 0) then
            hamiltonianME(j,:,i) = 0.0_double
         else
            hamiltonianME(j,:,i) = hamiltonianME(selectGaussianCount,:,i)
            selectGaussianCount = selectGaussianCount - 1
         endif

      enddo
   enddo

end subroutine correctHamiltonianIndex

end module DiagonalizationSubs
