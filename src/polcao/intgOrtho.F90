module O_Orthogonalization

   ! Import necessary modules.
   use O_Kinds

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

! Reminder about how the whole orthogonalization thing works. Consult the book
!   by Ching and Rulis and the reference: Ching WY, Lin CC. Orthogonalized
!   linear combinations of atomic orbitals: Application to the calculation of
!   energy bands of Si III. Physical Review B, 12(12), 5536, (1975).

! We will look at the case of the overlap matrix first. Note first a bit of
!   unusual notation. Let |v'> be a single element of the set of new valence
!   basis functions (from all atoms), let |c> be a single element of the set
!   of core basis functions (from all atoms), and let |v> be a single element
!   of the set of current valence basis functions (from all atoms). We will let
!   <c|v'> represent a single matrix element of the overlap between a
!   particular core basis function and a particular valence basis function.
!   Don't think of these as state vectors. Instead, they are just a particular
!   member of the valeDim or coreDim set of basis functions for the system.

! The basic goal is to take the existing valence basis functions (that created
!   the existing VV CV, VC, and CC matrices of overlap matrix elements) and
!   create a new set of basis functions that are orthogonal to the core basis
!   functions.

!   That is <c|v'> = <v'|c> = 0 for all elements.

!   Essentially, each matrix element of interaction between a basis function
!   from the core and a basis function from the new valence is equal to zero.

! Looking at equation 3.45 from the Ching-Rulis book on OLCAO (2012) we have:
!   |v'> = |v> + SUM( C|c> ).

!   Here, the sum is over all core basis functions and there is a different
!   coefficient for each core:valence basis function pair. Thus, in words,
!   our new valence basis function elements are the same as the old ones plus a
!   modification that is based on some weighted contribution from *each*
!   core basis function. The weighting will be *just right* so that the
!   new valence basis function elements will be orthogonal to the core basis
!   functions.

! Our goal now is to define the value of those coefficients. We do that by the
!   following:

!   Consider the overlap of the new valence basis function with *a particular*
!   core basis function. (Multiply from the left by <c|.)
!   <c|v'> = <c|v> + SUM( <c|C|c> )

!   Here, the SUM is still a single summation over the elements of the *right*
!   set of core basis functions |c>, but we are now summing the overlap of each
!   of the core basis functions with *a particular* other core basis function.
!   The trick is that the core basis functions that are on the same atomic
!   site are orthogonal to each other, and the overlap of core basis functions
!   is *assumed* to be zero between basis functions on *different* atomic sites.
!   Thus the term in the SUM will act like a delta function and pull out only
!   one C for *the particular* core basis function that we multiplied from the
!   left with, <c|. (If the assumption is wrong, then the process is much more
!   difficult and we will have to do something different.)

!   In words, the overlap of a given core basis function with one of the new
!   valence basis functions is equal to the overlap of the given core basis
!   function and the old valence basis function plus a single coefficient that
!   is specific for that core:valence pair.

! Rewriting the equation in the light of that discussion we have:
!   <c|v'> = <c|v> + C

! Now, we *demand* orthogonality between |c> and |v'> so that <c|v'> = 0.
!   0 = <c|v> + C

!   Or, the orthogonalization coefficients between core and valence basis
!   functions are just the negative of the matrix elements of the overlap
!   between the core and valence basis functions.

! Expressly: C = -<c|v> and C* = -<v|c> recalling that these are all just
!   single elements and not matrices per se.

! Now, let us construct the new VV overlap matrix. We start with the original
!   definition of |v'> = |v> + SUM( C|c> ). and multiply from the left by <v'|
!   on the left-hand side and by |v> + SUM( C|c> ) from the left on the right
!   hand side. With C* = complex conjugate of C. That is:
!   <v'|v'> = (<v| + SUM( C*<c| )) * (|v> + SUM( C|c> ))

! Thus, distributing:
!   <v'|v'> = <v|v> + SUM( <v|C|c> ) + SUM( C*<c|v> ) + SUM(SUM( C*<c|C|c> ))

!   Note again that each SUM is over the terms of the core basis function
!   inside the SUM parentheses. The double sum will have two separate sets of
!   indices. We also remind the reader that the equation is to be understood
!   as expressing a single matrix element. Thus, the C inside the SUM can be
!   moved to any position within the SUM(). With that knowledge we rewrite the
!   equation as:

!   <v'|v'> = <v|v> + SUM( <v|c>C ) + SUM( C*<c|v> ) + SUM(SUM( C*<c|c>C ))

! We expand C = -<c|v> and C* = -<v|c> to get:
!   <v'|v'> = <v|v> - SUM( <v|c><c|v> ) - SUM( <v|c><c|v> ) + 
!             SUM(SUM( <v|c><c|c><c|v> ))

!   Note that the first two SUMs are now identical.

!   A very tricky part is that the double SUM is over two different indices,
!   and the symbols are all basically the same. So, to help distinguish I have
!   added some indices:

! Reducing and adding indices, we have:
!   <v'|v'> = <v|v> - 2*SUM( <v|c><c|v> ) + SUM(SUM( <v|cj><cj|ci><ci|v> ))

! We can now translate this equation into a sequence of matrix operations to
!   produce the entire V'V' matrix from the set of matrix elements <v'|v'>.

! We take the original VV matrix and first subtract 2*[VC*CV]. Note that VC is
!   the complex conjugate transpose (or dagger) of CV. Then we add VC*CC*CV.
!   It doesn't matter what order we do the multiplication in (meaning either
!   (VC*CC)*CV or VC*(CC*CV)) because the number of floating point operations
!   is the same in both cases.

! Some comments:
! 1. The VC and CV matrices are complex conjugate transposes of each other.
!    (a.k.a the dagger operation). 
! 2. If we only store CV during the integrals process, then the first operation
!    (2*[VC*CV]) will "be fine" because the zherk call will include the
!    conjugate transpose effect in its operation. However, the last operation
!    (VC*CC*CV) is a set of separate calls so we would need to compute and
!    store the conjugate transpose of CV to be able to do it. *Alternatively*
!    we can store both the VC and CV. This will prove to be the wiser option
!    because we will need to make repeated use of both VC and CV for later
!    orthogonalization operations.
! 3. An interesting thought. Because we *do* compute VC ahead of time we could
!    just use it directly to perform 2*[VC*CV] as a simple matrix
!    multiplication. This might be faster because it may eliminate the need to
!    compute the conjugate transpose *again*, *within* the zherk call. Need to
!    explore this.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Okay, now let's look at the case of orthogonalizing any VV matrix other than
!   the two-center overlap matrix.

! We have an similar situation, which is:
!   |v'> = |v> + SUM( C|c> ).

! The main difference is that now, we claim to already know the coefficients C.
!   The way that a given <v|v> kinetic energy matrix element needs to be
!   modified will be the same as the way that the overlap matrix needed to be
!   modified.

! Thus, we can immediately skip to:
!   <v'|v'> = (<v| + SUM( C*<c| )) * (|v> + SUM( C|c> ))

! Followed by:
!   <v'|v'> = <v|v> + SUM( <v|C|c> ) + SUM( C*<c|v> ) + SUM(SUM( C*<c|C|c> ))

!   Note again that each SUM is over the terms of the core basis function
!   inside the SUM parentheses. The double sum will have two separate sets of
!   indices. We also remind the reader that the equation is to be understood
!   as expressing a single matrix element (as opposed to a full matrix). Thus,
!   the C inside the SUM can be moved to any position within the SUM(). With
!   that knowledge we rewrite the equation as:

!   <v'|v'> = <v|v> + SUM( <v|c>C ) + SUM( C*<c|v> ) + SUM(SUM( C*<c|c>C ))

! Now, we can substitute in the expressions for C and C*, but we need to make
!   special note that those coefficients were obtained from the overlap
!   orthogonalization. That is:
!   C = -<c_OL|v_OL> and C* = -<v_OL|c_OL>

! Upon substitution:
!   <v'|v'> = <v|v> - SUM( <v|c><c_OL|v_OL> ) - SUM( <v_OL|c_OL><c|v> ) + 
!             SUM(SUM( <v_OL|c_OL><c|c><c_OL|v_OL> ))

!   A very tricky part is that the double SUM is over two different indices,
!   and the symbols are all basically the same. So, to help distinguish I have
!   added some indices:

! Adding indices, we have:
!   <v'|v'> = <v|v> - SUM( <v|c><c_OL|v_OL> ) - SUM( <v_OL|c_OL><c|v> ) + 
!             SUM(SUM( <v_OL|c_OLj><cj|ci><c_OLi|v_OL> ))

! We can now translate this equation into a sequence of matrix operations to
!   produce the entire V'V' matrix from the set of matrix elements <v'|v'>.

! We take the original VV matrix and first subtract [VC*CV_OL] followed by
!   subtracting [VC_OL*CV]. (That first step is what zher2k does.) Note that VC
!   is the complex conjugate transpose (or dagger) of CV. Then we add
!   VC_OL*CC*CV_OL.

! Some Comments:
! 1. Note that for the non-overlap matrices we do not need to store the
!    conjugate transpose of CV. We will make *use* of the VC_OL and CV_OL, but
!    we only need the CV, not the VC. The first operation -[VC*CV_OL]-[VC_OL*CV]
!    will "be fine" because it is all taken care of by the zher2k call. The
!    conjugate transpose of CV is included. However, the second operation does
!    not make use of VC, only VC_OL, CC, and CV_OL.

#ifndef GAMMA

subroutine valeCoreCoreValeOL (valeDim,coreDim,descriptCV_OL,descriptVV,&
      & localCV_OL,localVV)

   ! Import the necessary modules
   use O_Kinds
   use pzherkInterface
   use MPI

   ! Make sure no variables are accidently declared
   implicit none

   ! Define passed parameters
   integer, intent(in) :: valeDim ! N
   integer, intent(in) :: coreDim ! K
   integer, dimension(9), intent(in) :: descriptCV_OL ! DESCA
   integer, dimension(9), intent(in) :: descriptVV ! DESCC
   complex (kind=double), dimension (:,:), intent(inout) :: localCV_OL ! A
   complex (kind=double), dimension (:,:), intent(inout) :: localVV ! C

   ! Do a rank 1 update. That is: C := alpha*A**H*A + beta*C where H equals
   !   complex conjugate and transpose. Also, alpha = -2 and beta = 1. The
   !   matrix A is the coreVale and C is the valeVale.
   ! Note that the documentation on
   !   http://www.netlib.org/scalapack/pblas_qref.html#PvHERK seems to have an
   !   error. The second variable "TRANS" is said to have possible values of
   !   'N' or 'T' to indicate different operations. However, in the definition
   !   of the variable "K" it mentions that TRANS could also be 'C'. The
   !   text often says something like, "Do abc if TRANS='N' and do xyz
   !   otherwise." So, I suppose that it looks for 'N' or 'n' and the other
   !   character could be anything. But, be careful anyway. The documentation
   !   from IBM seems better: Search for Parallel ESSL.
   call pzherk('U','C',valeDim,coreDim,(-2.0_double,0.0_double),&
         & localCV_OL,0,0,descriptCV_OL,(1.0_double,0.0_double),localVV,0,0,&
         & descriptVV)

   ! For reference, the non-parallel approach would look like the following
   !   with (first) blas and (second) straight matmul:

   !call zherk('U','C',valeDim,coreDim,-2.0_double,coreValeOL,coreDim,&
   !      & 1.0_double,valeVale,valeDim)

   !valeVale = valeVale - 2.0_double * matmul(conjg(transpose(coreValeOL)),&
   !      & coreValeOL)

end subroutine valeCoreCoreValeOL

subroutine valeCoreCoreVale (valeDim,coreDim,descriptCV_OL,descriptCV,&
      & descriptVV,localCV_OL,localCV,localVV)

   ! Import the necessary modules
   use O_Kinds
   use pzher2kInterface
   use MPI

   ! Make sure no variables are accidently declared
   implicit none

   ! Define passed parameters
   integer, intent(in) :: valeDim ! N
   integer, intent(in) :: coreDim ! K
   integer, dimension(9), intent(in) :: descriptCV_OL ! DESCA
   integer, dimension(9), intent(in) :: descriptCV ! DESCB
   integer, dimension(9), intent(in) :: descriptVV ! DESCC
   complex (kind=double), dimension (:,:), intent(inout) :: localCV_OL ! A
   complex (kind=double), dimension (:,:), intent(inout) :: localCV ! B
   complex (kind=double), dimension (:,:), intent(inout) :: localVV ! C

   call pzher2k('U','C',valeDim,coreDim,(-1.0_double,0.0_double),&
         & localCV_OL,0,0,descriptCV_OL,localCV,0,0,descriptCV,&
         & (1.0_double,0.0_double),localVV,0,0,descriptVV)

   ! For reference, the non-parallel approach would look like the following
   !   with (first) blas and (second) straight matmul:

   !call zher2k('U','C',valeDim,coreDim,(-1.0_double,0.0_double),&
   !      & coreValeOL,coreDim,coreVale,coreDim,1.0_double,valeVale,valeDim)

   !valeVale = valeVale - matmul(conjg(transpose(coreValeOL)),coreVale) - &
   !      & matmul(conjg(transpose(coreVale)),coreValeOL)

end subroutine valeCoreCoreVale

subroutine valeCoreCoreCore (valeDim,coreDim,descriptVC,descriptCC,&
      & descriptVC_temp,localVC,localCC,localVC_temp)

   ! Import the necessary modules
   use O_Kinds
   use pzhemmInterface
   use MPI

   ! Make sure no variables are accidently declared
   implicit none

   ! Define passed parameters
   integer, intent(in) :: valeDim ! M
   integer, intent(in) :: coreDim ! N
   integer, dimension(9), intent(in) :: descriptVC ! DESCB
   integer, dimension(9), intent(in) :: descriptCC ! DESCA (Hermitian)
   integer, dimension(9), intent(in) :: descriptVC_temp ! DESCC (Resultant)
   complex (kind=double), dimension (:,:), intent(inout) :: localVC ! B
   complex (kind=double), dimension (:,:), intent(inout) :: localCC ! A (Herm)
   complex (kind=double), dimension (:,:), intent(inout) :: localVC_temp ! C

   ! Note that localVC_temp is the resultant matrix.

   call pzhemm ('R','U',valeDim,coreDim,(1.0_double,0.0_double),&
         & localCC,0,0,descriptCC,localVC,0,0,descriptVC,&
         & (0.0_double,0.0_double),localVC_temp,0,0,descriptVC_temp)

end subroutine valeCoreCoreCore

! This subroutine will execute a VC*CV matrix multiplication and store the
!   result in VV. (This should be Hermitian. This subroutine will also set
!   numbers less that 1e-10 to zero to improve compressed storage.
subroutine makeValeVale (valeDim,coreDim,descriptVC_temp,decriptCV,descriptVV,&
      & localVC_temp,localCV,localVV)

   ! Import the necessary modules
   use O_Kinds
   use pzgemmInterface
   use MPI

   ! Make sure no variables are accidently declared
   implicit none

   ! Define passed parameters
   integer, intent(in) :: valeDim
   integer, intent(in) :: coreDim
   integer, dimension(9), intent(in) :: descriptVC_temp ! DESCA
   integer, dimension(9), intent(in) :: descriptCV ! DESCB
   integer, dimension(9), intent(in) :: descriptVV ! DESCC (Resultant,Hermitian)
   complex (kind=double), dimension (:,:), intent(inout) :: localVC_temp ! A
   complex (kind=double), dimension (:,:), intent(inout) :: localCV ! B
   complex (kind=double), dimension (:,:), intent(inout) :: localVV ! C (Herm)

   ! Define loop control variables
   integer :: i,j
   integer :: currentIndex


   ! Define the small threshhold for eliminating resultant values
   real (kind=double) :: negligLimit
   negligLimit = real(1.0D-10,double)

   call pzgemm ('N','N',valeDim,valeDim,coreDim,(1.0_double,0.0_double),&
         & localVC_temp,0,0,descriptVC_temp,localCV,0,0,descriptCV,&
         & (1.0_double,0.0_double),localVV,0,0,descriptVV)

   ! Eliminate all values in localVV that are less than 1e-10.
   if (negligLimit /= 0.0_double) then
      do i = 1, size(localVV,2)
         do j = 1, size(localVV,1)
            if (abs(real(localVV(j,i)) <= negligLimit) then
               if (abs(aimag(localVV(j,i) <= negligLimit) then
                  localVV(j,i) == (0.0_double,0.0_double)
               else
                  localVV(j,i) == (0.0_double,aimag(localVV(j,i)),double)
               endif
            else
               if (abs(aimag(localVV(j,i) <= negligLimit) then
                  localVV(j,i) == (real(localVV(j,i)),0.0_double,double)
               endif
            endif
         enddo
      enddo
   endif

end subroutine makeValeVale

#else

subroutine valeCoreCoreValeOLGamma (valeDim,coreDim,descriptCV_OL,descriptVV,&
      & localCV_OL,localVV)

   ! Import the necessary modules
   use O_Kinds
   use pdsyrkInterface
   use MPI

   ! Make sure no variables are accidently declared
   implicit none

   ! Define passed parameters
   integer, intent(in) :: valeDim ! N
   integer, intent(in) :: coreDim ! K
   integer, dimension(9), intent(in) :: descriptCV_OL ! DESCA
   integer, dimension(9), intent(in) :: descriptVV ! DESCC
   real (kind=double), dimension (:,:), intent(inout) :: localCV_OL ! A
   real (kind=double), dimension (:,:), intent(inout) :: localVV ! C

   ! Do a rank 1 update. That is: C := alpha*A**T*A + beta*C where T equals
   !   transpose. Also, alpha = -2 and beta = 1. The matrix A is the coreVale
   !   and C is the valeVale.
   call pdsyrk('U','T',valeDim,coreDim,-2.0_double,localCV_OL,0,0,&
         & descriptCV_OL,1.0_double,localVV,0,0,descriptVV)

   ! For reference, the non-parallel approach would look like the following
   !   with (first) blas and (second) straight matmul:

   !call dsyrk('U','C',valeDim,coreDim,-2.0_double,coreValeOLGamma,coreDim,&
   !      & 1.0_double,valeValeOLGamma,valeDim)

   !valeValeOLGamma = valeValeOLGamma - 2.0_double * matmul(transpose( &
   !      & coreValeOLGamma),coreValeOLGamma)

end subroutine valeCoreCoreValeOLGamma


subroutine valeCoreCoreValeGamma (valeDim,coreDim,descriptCV_OL,descriptCV,&
      & descriptVV,localCV_OL,localCV,localVV)

   ! Import the necessary modules
   use O_Kinds
   use pdsyr2kInterface
   use MPI

   ! Make sure no variables are accidently declared
   implicit none

   ! Define passed parameters
   integer, intent(in) :: valeDim ! N
   integer, intent(in) :: coreDim ! K
   integer, dimension(9), intent(in) :: descriptCV_OL ! DESCA
   integer, dimension(9), intent(in) :: descriptCV ! DESCB
   integer, dimension(9), intent(in) :: descriptVV ! DESCC
   real (kind=double), dimension (:,:), intent(inout) :: localCV_OL ! A
   real (kind=double), dimension (:,:), intent(inout) :: localCV ! B
   real (kind=double), dimension (:,:), intent(inout) :: localVV ! C

   call pdsyr2k('U','T',valeDim,coreDim,-1.0_double,localCV_OL,0,0,&
         & descriptCV_OL,localCV,0,0,descriptCV,1.0_double,localVV,0,0,&
         & descriptVV)

   ! For reference, the non-parallel approach would look like the following
   !   with (first) blas and (second) straight matmul:

   !call dsyr2k('U','T',valeDim,coreDim,-1.0_double,coreValeOLGamma,coreDim,&
   !      & coreValeGamma,coreDim,1.0_double,valeValeGamma,valeDim)

   !valeValeGamma = valeValeGamma - matmul(transpose(coreValeOLGamma),&
   !      & coreValeGamma) - matmul(transpose(coreValeGamma),coreValeOLGamma)

end subroutine valeCoreCoreValeGamma



subroutine valeCoreCoreCoreGamma (valeDim,coreDim,descriptVC,descriptCC,&
      & descriptVC_temp,localVC,localCC,localVC_temp)

   ! Import the necessary modules
   use O_Kinds
   use pdsymmInterface
   use MPI

   ! Make sure no variables are accidently declared
   implicit none

   ! Define passed parameters
   integer, intent(in) :: valeDim ! M
   integer, intent(in) :: coreDim ! N
   integer, dimension(9), intent(in) :: descriptVC ! DESCB
   integer, dimension(9), intent(in) :: descriptCC ! DESCA (Symmetric)
   integer, dimension(9), intent(in) :: descriptVC_temp ! DESCC (Resultant)
   real (kind=double), dimension (:,:), intent(inout) :: localVC ! B
   real (kind=double), dimension (:,:), intent(inout) :: localCC ! A (Symm)
   real (kind=double), dimension (:,:), intent(inout) :: localVC_temp ! C

   ! Note that localVC_temp is the resultant matrix.

   call pdsymm ('R','U',valeDim,coreDim,1.0_double,localCC,0,0,descriptCC,&
         & localVC,0,0,descriptVC,0.0_double,localVC_temp,0,0,descriptVC_temp)

end subroutine valeCoreCoreCoreGamma


subroutine makeValeValeGamma (valeDim,coreDim,descriptVC_temp,decriptCV,&
      & descriptVV,localVC_temp,localCV,localVV)

   ! Import the necessary modules
   use O_Kinds
   use pdgemmInterface
   use MPI

   ! Make sure no variables are accidently declared
   implicit none

   ! Define passed parameters
   integer, intent(in) :: valeDim
   integer, intent(in) :: coreDim
   integer, dimension(9), intent(in) :: descriptVC_temp ! DESCA
   integer, dimension(9), intent(in) :: descriptCV ! DESCB
   integer, dimension(9), intent(in) :: descriptVV ! DESCC (Resultant,Symmetric)
   real (kind=double), dimension (:,:), intent(inout) :: localVC_temp ! A
   real (kind=double), dimension (:,:), intent(inout) :: localCV ! B
   real (kind=double), dimension (:,:), intent(inout) :: localVV ! C (Symm)

   ! Define loop control variables
   integer :: i,j
   integer :: currentIndex

   ! Define the small threshhold for eliminating resultant values
   real (kind=double) :: negligLimit
   negligLimit = real(1.0D-10,double)

   call pdgemm ('N','N',valeDim,valeDim,coreDim,1.0_double,localVC_temp,0,0,&
         & descriptVC_temp,localCV,0,0,descriptCV,1.0_double,localVV,0,0,&
         & descriptVV)

   ! Eliminate all values in localVV that are less than 1e-10.
   if (negligLimit /= 0.0_double) then
      do i = 1, size(localVV,2)
         do j = 1, size(localVV,1)
            if (abs(localVV(j,i) <= negligLimit) then
               localVV(j,i) == 0.0_double
            endif
         enddo
      enddo
   endif

end subroutine makeValeValeGamma

#endif

end module O_Orthogonalization
