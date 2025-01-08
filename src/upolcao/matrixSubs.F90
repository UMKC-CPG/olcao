module O_MatrixSubs

   ! Use necessary modules.
   use O_Kinds

   contains

#ifndef GAMMA

subroutine matrixElementMult(summation,matrix1,matrix2,dim1,dim2)

   ! Import the precision variables and the type definitions.
   use O_Kinds

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   real (kind=double), intent(inout) :: summation
   integer, intent(in) :: dim1
   integer, intent(in) :: dim2
   real (kind=double), intent(in), dimension (dim1,dim2*(dim2+1)/2) :: matrix1
   real (kind=double), intent(in), dimension (dim1,dim2*(dim2+1)/2) :: matrix2

   ! Define the local variables
   integer :: i,j
   integer :: diagonalIndex
   real (kind=double) :: tempSummation

   ! Initialize the tempSummation so that we do not interfere with the previous
   !   kpoint's result when we multiply by 2 at the end.
   tempSummation = 0.0_double


   ! Initialize the counter for when the loop is at the diagonal indices of the
   !   packed matrix.
   diagonalIndex = 0

   do i = 1,dim2
      do j = diagonalIndex+1, diagonalIndex+i-1
         tempSummation = tempSummation + &
               & matrix1(1,j)*matrix2(1,j) + matrix1(2,j)*matrix2(2,j)
      enddo
      diagonalIndex = diagonalIndex + i
      tempSummation = tempSummation + matrix1(1,diagonalIndex) * &
            & matrix2(1,diagonalIndex) / 2.0_double
   enddo


   ! Multiply by two for the lower half of the hermitian matrix.
   summation = summation + tempSummation * 2.0_double

end subroutine matrixElementMult

#else

subroutine matrixElementMultGamma(summation,matrix1,matrix2,dim1,dim2)

   ! Import the precision variables and the type definitions.
   use O_Kinds

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   real (kind=double), intent(inout) :: summation
   integer, intent(in) :: dim1
   integer, intent(in) :: dim2
   real (kind=double), intent(in), dimension (dim1,dim2*(dim2+1)/2) :: matrix1
   real (kind=double), intent(in), dimension (dim1,dim2*(dim2+1)/2) :: matrix2

   ! Define the local variables
   integer :: i
   integer :: diagonalIndex

   ! Here, we do not have to have a tempSummation or initialize the summation
   !   because there is only 1 kpoint and the summation was already
   !   initialized at the beginning of valeCharge.

   do i = 1, dim2*(dim2+1)/2
      summation = summation + matrix1(1,i)*matrix2(1,i)
   enddo

   ! Multiply by two for the lower half of the symmetric matrix.
   summation = summation * 2.0_double

   ! Correct for the fact that the diagonal does not need to be multiplied by 2.
   diagonalIndex = 0
   do i = 1, dim2
      diagonalIndex = diagonalIndex + i
      summation = summation - matrix1(1,diagonalIndex)*matrix2(1,diagonalIndex)
   enddo
end subroutine matrixElementMultGamma

#endif


#ifndef GAMMA

subroutine readMatrix (datasetIDs,matrix,tempRealMatrix,tempImagMatrix,&
         & matrixDims,dim1,dim2,immediateBroadcast)

   ! Import the precision variables and the type definitions.
   use O_Kinds
   use HDF5
   use MPI_F08
   use O_MPI

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer, intent(in) :: dim1
   integer, intent(in) :: dim2
   integer (hid_t), intent(in), dimension (2) :: datasetIDs
   complex (kind=double), intent(out), dimension (dim1,dim2) :: matrix
   real (kind=double), intent(out), dimension (dim1,dim2) :: tempRealMatrix
   real (kind=double), intent(out), dimension (dim1,dim2) :: tempImagMatrix
   integer (hsize_t), intent(in), dimension (2) :: matrixDims
   integer, intent(in) :: immediateBroadcast

   ! Define local variables
   integer :: hdferr

   ! Read the real part.
   if (mpiRank == 0) then
      call h5dread_f (datasetIDs(1),H5T_NATIVE_DOUBLE,tempRealMatrix(:,:),&
            & matrixDims,hdferr)
   endif
   if (immediateBroadcast == 1) then
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) then
         call stopMPI("Failed to read temp real matrix.")
      endif
   endif

   ! Read the imaginary part.
   if (mpiRank == 0) then
      call h5dread_f (datasetIDs(2),H5T_NATIVE_DOUBLE,tempImagMatrix(:,:),&
            & matrixDims,hdferr)
   endif
   if (immediateBroadcast == 1) then
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) then
         call stopMPI("Failed to read temp imag matrix.")
      endif
   endif

   ! Assign the wave function values.
   matrix(:,:) = cmplx(tempRealMatrix(:,:),tempImagMatrix(:,:),double)
   if (immediateBroadcast == 1) then
      call MPI_BCAST(matrix(:,:),dim1*dim2,MPI_DOUBLE_COMPLEX,0,&
            & MPI_COMM_WORLD,mpierr)
   endif

end subroutine readMatrix

#else

subroutine readMatrixGamma (datasetID,matrix,matrixDims,dim1,dim2,&
      & immediateBroadcast)

   ! Import the precision variables and the type definitions.
   use O_Kinds
   use HDF5
   use MPI_F08
   use O_MPI

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer, intent(in) :: dim1
   integer, intent(in) :: dim2
   integer (hid_t), intent(in) :: datasetID
   real (kind=double), intent(out), dimension (dim1,dim2) :: matrix
   integer (hsize_t), intent(in), dimension (2) :: matrixDims
   integer, intent(in) :: immediateBroadcast

   ! Define local variables
   integer :: hdferr

   ! Read the real part and store it.
   if (mpiRank == 0) then
      call h5dread_f (datasetID,H5T_NATIVE_DOUBLE,matrix(:,:),matrixDims,hdferr)
   endif
   if (immediateBroadcast == 1) then
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) then
         call stopMPI("Failed to read real matrix.")
      endif
      call MPI_BCAST(matrix(:,:),dim1*dim2,MPI_DOUBLE_PRECISION,0,&
            & MPI_COMM_WORLD,mpierr)
   endif

end subroutine readMatrixGamma

#endif

! This subroutine does not need two versions for the GAMMA and non-GAMMA cases,
!   but it does need two versions for the accumulate and non-accumulate cases.
subroutine readPackedMatrixAccum (datasetID,packedMatrix,tempPackedMatrix,&
      & matrixDims,multFactor,dim1,dim2,immediateBroadcast)

   ! Import the precision variables and the type definitions.
   use O_Kinds
   use HDF5
   use MPI_F08
   use O_MPI

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer, intent(in) :: dim1
   integer, intent(in) :: dim2
   integer (hid_t), intent(in) :: datasetID
   real (kind=double), intent(inout), dimension (dim1,dim2*(dim2+1)/2) &
         & :: packedMatrix
   real (kind=double), intent(inout), dimension (dim1,dim2*(dim2+1)/2) &
         & :: tempPackedMatrix
   integer (hsize_t), intent(in), dimension (2) :: matrixDims
   real (kind=double), intent(in) :: multFactor
   integer, intent(in) :: immediateBroadcast

   ! Define local variables
   integer :: hdferr

   ! Read the matrix into a temp matrix first so we can operate on it.
   if (mpiRank == 0) then
      call h5dread_f (datasetID,H5T_NATIVE_DOUBLE,tempPackedMatrix,matrixDims,&
            & hdferr)
   endif
   if (immediateBroadcast == 1) then
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) then
         call stopMPI("Failed to read packed matrix into a temp matrix.")
      endif
   endif

   ! When reading the hamiltonian terms for the scf iterations, it is
   !   necessary to multiply each matrix by the appropriate coefficient.
   !   If the multFactor coefficient is not zero, then that is done here.
   if (multFactor == 0.0_double) then
      packedMatrix = packedMatrix + tempPackedMatrix
   else
      packedMatrix = packedMatrix + tempPackedMatrix * multFactor
   endif
   if (immediateBroadcast == 1) then
      call MPI_BCAST(packedMatrix(:,:),dim1*dim2,MPI_DOUBLE_PRECISION,0,&
            & MPI_COMM_WORLD,mpierr)
   endif

end subroutine readPackedMatrixAccum


subroutine readPackedMatrix (datasetID,packedMatrix,matrixDims,dim1,dim2,&
      & immediateBroadcast)

   ! Import the precision variables and the type definitions.
   use O_Kinds
   use HDF5
   use MPI_F08
   use O_MPI

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer, intent(in) :: dim1
   integer, intent(in) :: dim2
   integer (hid_t), intent(in) :: datasetID
   real (kind=double), intent(out), dimension (dim1,dim2*(dim2+1)/2) &
         & :: packedMatrix
   integer (hsize_t), intent(in), dimension (2) :: matrixDims
   integer, intent(in) :: immediateBroadcast

   ! Define local variables
   integer :: hdferr

   ! Read the matrix directly.
   if (mpiRank == 0) then
      call h5dread_f (datasetID,H5T_NATIVE_DOUBLE,packedMatrix,matrixDims,hdferr)
   endif
   if (immediateBroadcast == 1) then
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) then
         call stopMPI("Failed to read packed matrix.")
      endif
      call MPI_BCAST(packedMatrix(:,:),dim1*dim2*(dim2+1)/2,&
            & MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
   endif

end subroutine readPackedMatrix

#ifndef GAMMA

subroutine unpackMatrix (matrix,packedMatrix,dim2,fullFlag)

   ! Import the precision variables and the type definitions.
   use O_Kinds

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer :: dim2
   complex (kind=double), dimension (dim2,dim2) :: matrix
   real (kind=double), dimension (2,dim2*(dim2+1)/2) :: packedMatrix
   integer :: fullFlag

   ! Define local variables.
   integer :: i,j,currentIndex

   ! Initialize the index counter for the packed matrix.
   currentIndex = 0

   if (fullFlag == 0) then
      do i = 1, dim2
         do j = 1,i
            currentIndex = currentIndex + 1
            matrix(j,i) = cmplx(packedMatrix(1,currentIndex),&
                              & packedMatrix(2,currentIndex),double)
         enddo
      enddo
   else
      do i = 1, dim2
         do j = 1,i
            currentIndex = currentIndex + 1
            matrix(j,i) = cmplx(packedMatrix(1,currentIndex),&
                             &  packedMatrix(2,currentIndex),double)
            matrix(i,j) = cmplx(packedMatrix(1,currentIndex),&
                             & -packedMatrix(2,currentIndex),double)
         enddo
!         matrix(i,i) = matrix(i,i)/2.0_double
      enddo
!      matrix(:,:) = matrix(:,:) + transpose(conjg(matrix(:,:)))
   endif
end subroutine unpackMatrix

#else

subroutine unpackMatrixGamma (matrix,packedMatrix,dim2,fullFlag)

   ! Import the precision variables and the type definitions.
   use O_Kinds

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer :: dim2
   real (kind=double), dimension (dim2,dim2) :: matrix
   real (kind=double), dimension (1,dim2*(dim2+1)/2) :: packedMatrix
   integer :: fullFlag

   ! Define local variables.
   integer :: i,j,currentIndex

   ! Initialize the index counter for the packed matrix.
   currentIndex = 0

   if (fullFlag == 0) then
      do i = 1, dim2
         do j = 1,i
            currentIndex = currentIndex + 1
            matrix(j,i) = packedMatrix(1,currentIndex)
         enddo
      enddo
   else
      do i = 1, dim2
         do j = 1,i
            currentIndex = currentIndex + 1
            matrix(j,i) = packedMatrix(1,currentIndex)
            matrix(i,j) = packedMatrix(1,currentIndex)
         enddo
!         matrix(i,i) = matrix(i,i)/2.0_double
      enddo
!      matrix(:,:) = matrix(:,:) + transpose(matrix(:,:))
   endif
end subroutine unpackMatrixGamma

#endif

#ifndef GAMMA

subroutine packMatrix (matrix,packedMatrix,dim2)

   ! Import the precision variables and the type definitions.
   use O_Kinds

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer :: dim2
   complex (kind=double), dimension (dim2,dim2) :: matrix
   real (kind=double), dimension (2,dim2*(dim2+1)/2) :: packedMatrix

   ! Define local variables.
   integer :: i,j,currentIndex

   ! Initialize the index counter for the packed matrix.
   currentIndex = 0

   do i = 1, dim2
      do j = 1,i
         currentIndex = currentIndex + 1
         packedMatrix(1,currentIndex) =  real(matrix(j,i),double)
         packedMatrix(2,currentIndex) = aimag(matrix(j,i))
      enddo
   enddo
end subroutine packMatrix

#else

subroutine packMatrixGamma (matrix,packedMatrix,dim2)

   ! Import the precision variables and the type definitions.
   use O_Kinds

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer :: dim2
   real (kind=double), dimension (dim2,dim2) :: matrix
   real (kind=double), dimension (1,dim2*(dim2+1)/2) :: packedMatrix

   ! Define local variables.
   integer :: i,j,currentIndex

   ! Initialize the index counter for the packed matrix.
   currentIndex = 0

   do i = 1, dim2
      do j = 1,i
         currentIndex = currentIndex + 1
         packedMatrix(1,currentIndex) =  matrix(j,i)
      enddo
   enddo
end subroutine packMatrixGamma

#endif


! This is currently not used by any program.
!subroutine unfold (inMatrix, outMatrix, matrixSize, addFlag)
!
!   ! Import the precision variables
!   use O_Kinds
!
!   ! Make sure that there are no accidental variable declarations.
!   implicit none
!
!   ! Define the dummy variables that are passed to this subroutine.
!   real    (kind=double), dimension (:,:) :: inMatrix
!   complex (kind=double), dimension (:,:) :: outMatrix
!   integer :: matrixSize
!   integer :: addFlag
!
!   ! Define local variables used for unfolding the matrix.
!   integer :: i,j ! Loop index variables
!
!   if (addFlag == 0) then
!      do i = 1, matrixSize - 1
!         do j = i + 1, matrixSize
!            outMatrix(i,j) = cmplx(inMatrix(i,j),-inMatrix(j,i),double)
!            outMatrix(j,i) = cmplx(inMatrix(i,j), inMatrix(j,i),double)
!         enddo
!         outMatrix(i,i) = cmplx(inMatrix(i,i),0.0_double,double)
!      enddo
!      outMatrix(matrixSize,matrixSize) = cmplx(inMatrix(matrixSize,matrixSize),&
!            & 0.0_double,double)
!   else
!      do i = 1, matrixSize - 1
!         do j = i + 1, matrixSize
!            outMatrix(i,j) = outMatrix(i,j) + &
!                  & cmplx(inMatrix(i,j),-inMatrix(j,i),double)
!            outMatrix(j,i) = outMatrix(j,i) + &
!                  & cmplx(inMatrix(i,j), inMatrix(j,i),double)
!         enddo
!         outMatrix(i,i) = outMatrix(i,i) + &
!               & cmplx(inMatrix(i,i),0.0_double,double)
!      enddo
!      outMatrix(matrixSize,matrixSize) = outMatrix(matrixSize,matrixSize) + &
!            & cmplx(inMatrix(matrixSize,matrixSize),0.0_double,double)
!   endif
!
!end subroutine unfold


end module O_MatrixSubs
