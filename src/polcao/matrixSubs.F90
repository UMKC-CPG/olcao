module O_MatrixSubs

   ! Use necessary modules.
   use O_Kinds

   contains

#ifndef GAMMA

subroutine matrixElementMult(summation,matrix1,matrix2,dim1,dim2)

   ! Import the precision variables and the type definitions.
   use O_Kinds

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   real (kind=double), intent(inout) :: summation
   integer, intent(in) :: dim1
   integer, intent(in) :: dim2
   real (kind=double), dimension (dim1,dim2), intent(in) :: matrix1
   real (kind=double), dimension (dim1,dim2), intent(in) :: matrix2

   ! Define the local variables
   integer :: i,j
   real (kind=double) :: tempSummation

   do i=1,dim1
      tempSummation = 0.0_double
      do j=i+1,dim2
         tempSummation = tempSummation &
               & + matrix1(i,j) * matrix2(i,j) &
               & + matrix1(j,i) * matrix2(j,i)
     enddo
     summation = summation + tempSummation * 2.0_double &
           & + matrix1(i,i)*matrix2(i,i)
   enddo

end subroutine matrixElementMult

#else

subroutine matrixElementMultGamma(summation,matrix1,matrix2,dim1,dim2)

   ! Import the precision variables and the type definitions.
   use O_Kinds

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   real (kind=double), intent(inout) :: summation
   integer, intent(in) :: dim1
   integer, intent(in) :: dim2
   real (kind=double), dimension (dim1,dim2), intent(in) :: matrix1
   real (kind=double), dimension (dim1,dim2), intent(in) :: matrix2

   ! Define the local variables
   integer :: i
   integer :: diagonalIndex

   do i=1,dim1
      tempSummation = 0.0_double
      do j=i+1,dim2
         tempSummation = tempSummation &
               & + matrix1(j,i) * matrix2(j,i)
     enddo
     summation = summation + tempSummation * 2.0_double &
           & + matrix1(i,i)*matrix2(i,i)
   enddo
end subroutine matrixElementMultGamma

#endif

#ifndef GAMMA

subroutine readPartialWaveFns (datasetIDs,matrix,tempRealMatrix,&
      & tempImagMatrix,matrixDims,init,fin,dim1,dim2)

   ! Import the precision variables and the type definitions.
   use O_Kinds
   use HDF5

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters. (Note that the temp matrices could
   !   have been defined as local variables, but then we would need to
   !   repeatedly create and destroy them.  By using the current approach
   !   we just allocate once and then pass it to this subroutine a few times.
   integer, intent(in) :: dim1
   integer, intent(in) :: dim2
   integer (hid_t), dimension (2), intent(in) :: datasetIDs
   complex (kind=double), dimension (dim1,dim2), intent(inout) :: matrix
   real (kind=double), dimension (dim1,dim2), intent(out) :: tempRealMatrix
   real (kind=double), dimension (dim1,dim2), intent(out) :: tempImagMatrix
   integer (hsize_t), dimension (2), intent(in) :: matrixDims
   integer, intent(in) :: init
   integer, intent(in) :: fin

   ! Define local variables
   integer :: hdferr


   ! Read the real part.
   call h5dread_f (datasetIDs(1),H5T_NATIVE_DOUBLE,tempRealMatrix(:,:),&
         & matrixDims,hdferr)
   if (hdferr /= 0) stop 'Failed to read temp real matrix'

   ! Read the imaginary part.
   call h5dread_f (datasetIDs(2),H5T_NATIVE_DOUBLE,tempImagMatrix(:,:),&
         & matrixDims,hdferr)
   if (hdferr /= 0) stop 'Failed to read temp imag matrix'

   ! Copy the portion that was read into the wave function matrix.
   matrix(:,init:fin) = cmplx(tempRealMatrix(:,init:fin),&
         & tempImagMatrix(:,init:fin),double)

end subroutine readPartialWaveFns

#else

subroutine readPartialWaveFnsGamma (datasetID,matrix,tempRealMatrix,&
      & matrixDims,init,fin,dim1,dim2)

   ! Import the precision variables and the type definitions.
   use O_Kinds
   use HDF5

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters. (See note about temp matrix in
   !   readPartialWaveFns.)
   integer, intent(in) :: dim1
   integer, intent(in) :: dim2
   integer (hid_t), intent(in) :: datasetID
   real (kind=double), dimension (dim1,dim2), intent(inout) :: matrix
   real (kind=double), dimension (dim1,dim2), intent(out) :: tempRealMatrix
   integer (hsize_t), dimension (2), intent(in) :: matrixDims
   integer, intent(in) :: init
   integer, intent(in) :: fin

   ! Define local variables
   integer :: hdferr


   ! Read the real part.
   call h5dread_f (datasetID,H5T_NATIVE_DOUBLE,tempRealMatrix(:,:),&
         & matrixDims,hdferr)
   if (hdferr /= 0) stop 'Failed to read temp real matrix'

   ! Copy the portion that was read into the wave function matrix.
   matrix(:,init:fin) = tempRealMatrix(:,init:fin)

end subroutine readPartialWaveFnsGamma

#endif

#ifndef GAMMA

subroutine readMatrix (datasetIDs,matrix,tempRealMatrix,tempImagMatrix,&
         & matrixDims,dim1,dim2)

   ! Import the precision variables and the type definitions.
   use O_Kinds
   use HDF5

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters. (See note above about temp matrix.)
   integer, intent(in) :: dim1
   integer, intent(in) :: dim2
   integer (hid_t), dimension (2), intent(in) :: datasetIDs
   complex (kind=double), dimension (dim1,dim2), intent(inout) :: matrix
   real (kind=double), dimension (dim1,dim2), intent(out) :: tempRealMatrix
   real (kind=double), dimension (dim1,dim2), intent(out) :: tempImagMatrix
   integer (hsize_t), dimension (2), intent(in) :: matrixDims

   ! Define local variables
   integer :: hdferr

   ! Read the real part.
   call h5dread_f (datasetIDs(1),H5T_NATIVE_DOUBLE,tempRealMatrix(:,:),&
         & matrixDims,hdferr)
   if (hdferr /= 0) stop 'Failed to read temp real matrix'

   ! Read the imaginary part.
   call h5dread_f (datasetIDs(2),H5T_NATIVE_DOUBLE,tempImagMatrix(:,:),&
         & matrixDims,hdferr)
   if (hdferr /= 0) stop 'Failed to read temp imag matrix'

   ! Assign the wave function values.
   matrix(:,:) = cmplx(tempRealMatrix(:,:),tempImagMatrix(:,:),double)

end subroutine readMatrix

#else

subroutine readMatrixGamma (datasetID,matrix,matrixDims,dim1,dim2)

   ! Import the precision variables and the type definitions.
   use O_Kinds
   use HDF5

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer, intent(in) :: dim1
   integer, intent(in) :: dim2
   integer (hid_t), intent(in) :: datasetID
   real (kind=double), dimension (dim1,dim2), intent(inout) :: matrix
   integer (hsize_t), dimension (2), intent(in) :: matrixDims

   ! Define local variables
   integer :: hdferr

   ! Read the real part and store it.
   call h5dread_f (datasetID,H5T_NATIVE_DOUBLE,matrix(:,:),matrixDims,hdferr)
   if (hdferr /= 0) stop 'Failed to read real matrix'

end subroutine readMatrixGamma

#endif

! This subroutine does not need two versions for the GAMMA and non-GAMMA cases,
!   but it does need two versions for the accumulate and non-accumulate cases.
subroutine readPackedMatrixAccum (datasetID,packedMatrix,tempPackedMatrix,&
      & matrixDims,multFactor,dim1,dim2)

   ! Import the precision variables and the type definitions.
   use O_Kinds
   use HDF5

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer, intent(in) :: dim1
   integer, intent(in) :: dim2
   integer (hid_t), intent(in) :: datasetID
   real (kind=double), dimension (dim1,dim2), intent(inout) :: packedMatrix
   real (kind=double), dimension (dim1,dim2), intent(out) :: tempPackedMatrix
   integer (hsize_t), dimension (2), intent(in) :: matrixDims
   real (kind=double), intent(in) :: multFactor

   ! Define local variables
   integer :: hdferr

   ! Read the matrix into a temp matrix first so we can operate on it.
   call h5dread_f (datasetID,H5T_NATIVE_DOUBLE,tempPackedMatrix,matrixDims,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to read packed matrix into a temp matrix'

   ! When reading the hamiltonian terms for the scf iterations, it is
   !   necessary to multiply each matrix by the appropriate coefficient.
   !   If the multFactor coefficient is not zero, then that is done here.
   if (multFactor == 0.0_double) then
      packedMatrix = packedMatrix + tempPackedMatrix
   else
      packedMatrix = packedMatrix + tempPackedMatrix * multFactor
   endif

end subroutine readPackedMatrixAccum


subroutine readPackedMatrix (datasetID,packedMatrix,matrixDims,dim1,dim2)

   ! Import the precision variables and the type definitions.
   use O_Kinds
   use HDF5

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer, intent(in) :: dim1
   integer, intent(in) :: dim2
   integer (hid_t) :: datasetID
   real (kind=double), dimension (dim1,dim2), intent(out) :: packedMatrix
   integer (hsize_t), dimension (2), intent(in) :: matrixDims

   ! Define local variables
   integer :: hdferr

   ! Read the matrix directly.
   call h5dread_f (datasetID,H5T_NATIVE_DOUBLE,packedMatrix,matrixDims,hdferr)
   if (hdferr /= 0) stop 'Failed to read packed matrix'

end subroutine readPackedMatrix

#ifndef GAMMA

subroutine unpackMatrix (matrix,packedMatrix,dim2,fullFlag)

   ! Import the precision variables and the type definitions.
   use O_Kinds

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer, intent(in) :: dim2
   complex (kind=double), dimension (dim2,dim2), intent(out) :: matrix
   real (kind=double), dimension (dim2,dim2), intent(in) :: packedMatrix
   integer, intent(in) :: fullFlag ! The meaning of this flag is that when 0
         ! we only need to fill the upper triangle of matrix. When 1 we will
         ! fill the entire matrix.

   ! Recall that the packed matrix is stored as a single real matrix that
   !   represents a complex Hermitian matrix. The upper triangle of the packed
   !   matrix is the upper triangle of the real part of matrix(:,:). The lower
   !   triangle of the packed matrix is the upper triangle of the imaginary
   !   part of matrix(:,:).

   ! Define local variables.
   integer :: i,j

   ! Initialize the index counter for the packed matrix.
   currentIndex = 0

   if (fullFlag == 0) then ! We only fill the upper triangle of matrix(:,:).
      ! Iterate over the upper triangle of matrix(:,:).
      do i = 1, dim2
         do j = 1,i-1
            matrix(j,i) = cmplx(packedMatrix(j,i),packedMatrix(i,j),double)
         enddo
         matrix(i,i) = cmplx(packedMatrix(i,i),0.0_double,double)
      enddo
   else
      ! Iterate over the upper triangle of matrix(:,:).
      do i = 1, dim2
         do j = 1,i-1
            matrix(j,i) = cmplx(packedMatrix(j,i),packedMatrix(i,j),double)
            matrix(i,j) = conjg(cmplx(packedMatrix(j,i),packedMatrix(i,j),&
                  & double))
         enddo
         matrix(i,i) = cmplx(packedMatrix(i,i),0.0_double,double)
      enddo
   endif
end subroutine unpackMatrix

#else

subroutine unpackMatrixGamma (matrix,packedMatrix,dim2,fullFlag)

   ! Import the precision variables and the type definitions.
   use O_Kinds

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer, intent(in) :: dim2
   real (kind=double), dimension (dim2,dim2), intent(out) :: matrix
   real (kind=double), dimension (dim2,dim2), intent(in) :: packedMatrix
   integer, intent(in) :: fullFlag

   ! Define local variables.
   integer :: i,j

   if (fullFlag == 0) then
      do i = 1, dim2
         do j = 1,i
            matrix(j,i) = packedMatrix(j,i)
         enddo
      enddo
   else
      do i = 1, dim2
         do j = 1,dim2
            matrix(j,i) = packedMatrix(j,i)
         enddo
      enddo
   endif
end subroutine unpackMatrixGamma

#endif

#ifndef GAMMA

subroutine packMatrix (matrix,packedMatrix,dim2)

   ! Import the precision variables and the type definitions.
   use O_Kinds

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer, intent(in) :: dim2
   complex (kind=double), dimension (dim2,dim2), intent(in) :: matrix
   real (kind=double), dimension (dim2,dim2), intent(out) :: packedMatrix

   ! Note that the order in which the assignment is made makes a difference.
   !   The imaginary must be assigned first so that the real value will
   !   appear on the diagonal.  (The complex values on the diagonal are always
   !   equal to zero so they don't need to be stored explicitly.)

   ! As mentioned in unpackMatrix above, the packed matrix will store the real
   !   upper triangle of matrix(:,:) in its upper triangle (including the
   !   diagonal). The packed matrix will store the strict upper triangle of
   !   the imaginary part in its lower part.
   ! Loop over the upper triangle of matrix(:,:).
   do i = 1, dim2
      do j = 1,i
         packedMatrix(i,j) = aimag(matrix(j,i))
         packedMatrix(j,i) =  real(matrix(j,i))
      enddo
   enddo
end subroutine packMatrix

#else

subroutine packMatrixGamma (matrix,packedMatrix,dim2)

   ! Import the precision variables and the type definitions.
   use O_Kinds

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer, intent(in) :: dim2
   real (kind=double), dimension (dim2,dim2), intent(in) :: matrix
   real (kind=double), dimension (dim2,dim2), &
         & intent(out) :: packedMatrix

   ! Define local variables.
   integer :: i,j
   ! Initialize the lower triangle to zero. The upper triangle will be
   !   filled by the upper triangle of matrix(:,:). The reason to do
   !   this is to ensure that the matrix has a lot of (easily compressible)
   !   zeros in it before it is saved. I'm not sure if every compiler
   !   would ensure that, and thus if we tried to compress and save the
   !   matrix *without* doing this, we might be saving a lot of random
   !   noise which doesn't compress so well.
   do i = 1, dim2
      do j = i+1, dim2
         packedMatrix(j,i) = 0.0_double
      enddo
   enddo
   do i = 1, dim2
      do j = 1,i
         packedMatrix(j,i) =  matrix(j,i)
      enddo
   enddo
end subroutine packMatrixGamma

#endif

end module O_MatrixSubs
