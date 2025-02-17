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

#ifndef GAMMA

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
!   particular core basis function and a particular new valence basis function.
!   Don't think of these as state vectors. Instead, they are just a particular
!   member of the valeDim or coreDim set of basis functions for the system.

! The basic goal is to take the existing valence basis functions (that created
!   the existing VV, CV, VC, and CC matrices of overlap matrix elements) and
!   create a new set of basis functions that are orthogonal to the core basis
!   functions.

!   That is: we want <c|v'> = <v'|c> = 0 for all elements.

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
!   set of core basis functions C|c>, but we are now summing the overlap of each
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

! Rewriting the equation in light of that discussion we have:
!   <c|v'> = <c|v> + C

! Now, we *demand* orthogonality between |c> and |v'> so that <c|v'> = 0.
!   0 = <c|v> + C

!   Or, the orthogonalization coefficients between core and valence basis
!   functions are just the negative of the matrix elements of the overlap
!   between the core and valence basis functions. (I.e. Gramm-Schmidt
!   orthogonalization in which the projection of one vector on another is used
!   to determine the orthogonalization coefficient.)

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
!   produce the entire V'V' matrix from the set of matrix elements <v'|v'>,
!   which were produced from <v|v>, <c|c>, and <v|c>.

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

subroutine valeCoreCoreValeOL (valeDim,coreDim,valeVale,coreValeOL)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary interface.
   use zherkInterface

   ! Make sure no variables are accidently declared
   implicit none

   ! Define passed parameters
   integer :: valeDim
   integer :: coreDim
   complex (kind=double), dimension (valeDim,valeDim) :: valeVale
   complex (kind=double), dimension (coreDim,valeDim) :: coreValeOL

   call zherk('U','C',valeDim,coreDim,-2.0_double,coreValeOL,coreDim,&
         & 1.0_double,valeVale,valeDim)

!   valeVale = valeVale - 2.0_double * matmul(conjg(transpose(coreValeOL)),&
!         & coreValeOL)

end subroutine valeCoreCoreValeOL

subroutine valeCoreCoreVale (valeDim,coreDim,valeVale,coreVale,coreValeOL)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary interface.
   use zher2kInterface

   ! Make sure no variables are accidently declared
   implicit none

   ! Define passed parameters
   integer :: valeDim
   integer :: coreDim
   complex (kind=double), dimension (valeDim,valeDim) :: valeVale
   complex (kind=double), dimension (coreDim,valeDim) :: coreVale
   complex (kind=double), dimension (coreDim,valeDim) :: coreValeOL

   call zher2k('U','C',valeDim,coreDim,(-1.0_double,0.0_double),&
         & coreValeOL,coreDim,coreVale,coreDim,1.0_double,valeVale,valeDim)

!   valeVale = valeVale - matmul(conjg(transpose(coreValeOL)),coreVale) - &
!         & matmul(conjg(transpose(coreVale)),coreValeOL)

end subroutine valeCoreCoreVale

subroutine coreValeCoreCore (valeDim,coreDim,valeCore,coreVale,coreCore)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary blas interface
!   use zhemmInterface

   ! Make sure no variables are accidently declared
   implicit none

   ! Define passed parameters
   integer :: valeDim
   integer :: coreDim
   complex (kind=double), dimension (coreDim,valeDim) :: valeCore ! Transposed
   complex (kind=double), dimension (coreDim,ValeDim) :: coreVale
   complex (kind=double), dimension (coreDim,coreDim) :: coreCore

   integer :: i,j


!   call zhemm ('R','U',valeDim,coreDim,(1.0_double,0.0_double),&
!         & coreCore,coreDim,coreVale,coreDim,(1.0_double,0.0_double),&
!         & valeCore,valeDim)

   do i = 1,coreDim
      do j = 1, valeDim
         valeCore(i,j) = sum(conjg(coreVale(1:coreDim,j)) * &
              & coreCore(1:coreDim,i))
      enddo
   enddo

!   valeCore = matmul(conjg(transpose(coreVale)),coreCore)

end subroutine coreValeCoreCore

subroutine makeValeVale (valeDim,coreDim,packedValeDim,valeCore,coreVale,&
      & valeVale,packedValeVale,storeFlag,fullFlag)

   ! Import the necessary modules
   use O_Kinds

   ! Make sure no variables are accidently declared
   implicit none

   ! Define passed parameters
   integer :: valeDim
   integer :: coreDim
   integer :: packedValeDim
   complex (kind=double), dimension (coreDim,valeDim) :: valeCore ! Transposed
   complex (kind=double), dimension (coreDim,valeDim) :: coreVale
   complex (kind=double), dimension (valeDim,valeDim) :: valeVale
   real    (kind=double), dimension (2,packedValeDim * &
         & (packedValeDim+1)/2) :: packedValeVale
   integer :: storeFlag ! 1=Copy matrix for packed storage, 0=Do not.
   integer :: fullFlag  ! Only when storeFlag=0.  1=make full matrix, 0=do not.

   ! Define loop control variables
   integer :: i,j
   integer :: currentIndex


   ! Define the small threshhold for eliminating resultant values
   real (kind=double) :: smallThresh10
   smallThresh10 = real(1.0D-10,double)


   if (storeFlag == 0) then
      if (fullFlag == 0) then
         do i = 1, valeDim
            do j = 1, i
               valeVale(j,i) = valevale(j,i) + sum( &
                     & valeCore(1:coreDim,j) * coreVale(1:coreDim,i))
            enddo
         enddo
      else
         do i = 1, valeDim
            do j = 1, i
               valeVale(j,i) = valevale(j,i) + sum( &
                     & valeCore(1:coreDim,j) * coreVale(1:coreDim,i))
               valeVale(i,j) = conjg(valeVale(j,i))
            enddo
         enddo
      endif
   else

      ! Initialize the current index for packing the valeVale matrix.
      currentIndex = 0

      do i = 1, valeDim
         do j = 1, i

            ! Compute the result.
            valeVale(j,i) = valeVale(j,i) + sum( &
                  & valeCore(1:coreDim,j) * coreVale(1:coreDim,i))

            ! Increment the array index for packing the valeVale matrix.
            currentIndex = currentIndex + 1

            ! Pack it.
            packedValeVale(1,currentIndex) =  real(valeVale(j,i),double)
            packedValeVale(2,currentIndex) = aimag(valeVale(j,i))

            ! Check for negligable components and reduce them to 0.
            if (abs(packedValeVale(1,currentIndex)) < smallThresh10) then
               packedValeVale(1,currentIndex) = 0.0_double
            endif
            if (abs(packedValeVale(2,currentIndex)) < smallThresh10) then
               packedValeVale(2,currentIndex) = 0.0_double
            endif
         enddo
      enddo
   endif
end subroutine makeValeVale

#else

subroutine valeCoreCoreValeOLGamma (valeDim,coreDim,valeValeOLGamma,&
      & coreValeOLGamma)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary interface.
   use dsyrkInterface

   ! Make sure no variables are accidently declared
   implicit none

   ! Define passed parameters
   integer :: valeDim
   integer :: coreDim
   real (kind=double), dimension (valeDim,valeDim) :: valeValeOLGamma
   real (kind=double), dimension (coreDim,valeDim) :: coreValeOLGamma

   call dsyrk('U','C',valeDim,coreDim,-2.0_double,coreValeOLGamma,coreDim,&
         & 1.0_double,valeValeOLGamma,valeDim)

!   valeValeOLGamma = valeValeOLGamma - 2.0_double * matmul(transpose( &
!         & coreValeOLGamma),coreValeOLGamma)

end subroutine valeCoreCoreValeOLGamma



subroutine valeCoreCoreValeGamma (valeDim,coreDim,valeValeGamma,&
      & coreValeGamma,coreValeOLGamma)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary interface.
   use dsyr2kInterface

   ! Make sure no variables are accidently declared
   implicit none

   ! Define passed parameters
   integer :: valeDim
   integer :: coreDim
   real (kind=double), dimension (valeDim,valeDim) :: valeValeGamma
   real (kind=double), dimension (coreDim,valeDim) :: coreValeGamma
   real (kind=double), dimension (coreDim,valeDim) :: coreValeOLGamma

   call dsyr2k('U','C',valeDim,coreDim,-1.0_double,coreValeOLGamma,coreDim,&
         & coreValeGamma,coreDim,1.0_double,valeValeGamma,valeDim)

!   valeValeGamma = valeValeGamma - matmul(transpose(coreValeOLGamma),&
!         & coreValeGamma) - matmul(transpose(coreValeGamma),coreValeOLGamma)

end subroutine valeCoreCoreValeGamma



subroutine coreValeCoreCoreGamma (valeDim,coreDim,valeCoreGamma,&
      & coreValeGamma,coreCoreGamma)

   ! Import the necessary modules
   use O_Kinds

   ! Make sure no variables are accidently declared
   implicit none

   ! Define passed parameters
   integer :: valeDim
   integer :: coreDim
   real (kind=double), dimension (coreDim,valeDim) :: valeCoreGamma !Transposed
   real (kind=double), dimension (coreDim,valeDim) :: coreValeGamma
   real (kind=double), dimension (coreDim,coreDim) :: coreCoreGamma

   ! Define local loop variables
   integer :: i,j

!   valeCoreGamma = matmul(transpose(coreValeGamma),coreCoreGamma)

   do i = 1,coreDim
      do j = 1, valeDim
         valeCoreGamma(i,j) = sum(coreValeGamma(1:coreDim,j) * &
              & coreCoreGamma(1:coreDim,i))
      enddo
   enddo

end subroutine coreValeCoreCoreGamma



subroutine makeValeValeGamma (valeDim,coreDim,packedValeDim,valeCoreGamma,&
      & coreValeGamma,valeValeGamma,packedValeVale,storeFlag,fullFlag)

   ! Import the necessary modules
   use O_Kinds

   ! Make sure no variables are accidently declared
   implicit none

   ! Define passed parameters
   integer :: valeDim
   integer :: coreDim
   integer :: packedValeDim
   real (kind=double), dimension (coreDim,valeDim) :: valeCoreGamma !Transposed
   real (kind=double), dimension (coreDim,valeDim) :: coreValeGamma
   real (kind=double), dimension (valeDim,valeDim) :: valeValeGamma
   real (kind=double), dimension (1,packedValeDim * &
         & (packedValeDim+1)/2) :: packedValeVale
   integer :: storeFlag
   integer :: fullFlag

   ! Define local loop control
   integer :: i,j
   integer :: currentIndex

   ! Define the small threshhold for eliminating resultant values
   real (kind=double) :: smallThresh10
   smallThresh10 = real(1.0D-10,double)

   if (storeFlag == 0) then
      if (fullFlag == 0) then
         do i = 1, valeDim
            do j = 1, i
               valeValeGamma(j,i) = valeValeGamma(j,i) + sum( &
                     & valeCoreGamma(1:coreDim,j) * coreValeGamma(1:coreDim,i))
            enddo
         enddo
      else
         do i = 1, valeDim
            do j = 1, i
               valeValeGamma(j,i) = valeValeGamma(j,i) + sum( &
                     & valeCoreGamma(1:coreDim,j) * coreValeGamma(1:coreDim,i))
               valeValeGamma(i,j) = valeValeGamma(j,i)
            enddo
         enddo
      endif
   else

      ! Initialize the current index for packing the valeVale matrix.
      currentIndex = 0

      do i = 1, valeDim
         do j = 1, i

            ! Increment the array index for packing the valeVale matrix.
            currentIndex = currentIndex + 1

            ! Compute the result.
            valeValeGamma(j,i) = valeValeGamma(j,i) + sum( &
                  & valeCoreGamma(1:coreDim,j) * coreValeGamma(1:coreDim,i))

            ! Pack it.
            packedValeVale(1,currentIndex) = valeValeGamma(j,i)

            ! Check for negligable components and reduce them to 0.
            if (abs(packedValeVale(1,currentIndex)) < smallThresh10) then
               packedValeVale(1,currentIndex) = 0.0_double
            endif
         enddo
      enddo
   endif
end subroutine makeValeValeGamma

#endif

end module O_Orthogonalization
