module matSubs
  use HDF5
  use dataStructs
  
  implicit none

  contains
subroutine unpackParMatrix (matrix,packedMatrix,dim2,fullFlag)
   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer :: dim2
   complex (kind=double), dimension (dim2,dim2) :: matrix
!   real (kind=double), dimension (2,dim2*(dim2+1)/2) :: packedMatrix
   real (kind=double), dimension (dim2,dim2) :: packedMatrix
   integer :: fullFlag

   ! Define local variables.
   integer :: i,j,currentIndex

   ! Initialize the index counter for the packed matrix.
   currentIndex = 0

!   if (fullFlag == 0) then
      do i = 1, dim2
         do j = 1,i
!            currentIndex = currentIndex + 1
            if (i==j) then
              matrix(j,i) = cmplx(packedMatrix(j,i),0.0_double,double)
            else
              matrix(j,i) = cmplx(packedMatrix(j,i),&
                              & packedMatrix(i,j),double)
            endif
         enddo
      enddo
  ! else
  !    do i = 1, dim2
  !       do j = 1,i
  !          currentIndex = currentIndex + 1
  !          matrix(j,i) = cmplx(packedMatrix(1,currentIndex),&
  !                           &  packedMatrix(2,currentIndex),double)
  !          matrix(i,j) = cmplx(packedMatrix(1,currentIndex),&
  !                           & -packedMatrix(2,currentIndex),double)
  !       enddo
! !        matrix(i,i) = matrix(i,i)/2.0_double
  !    enddo
! !     matrix(:,:) = matrix(:,:) + transpose(conjg(matrix(:,:)))
  ! endif
end subroutine unpackParMatrix

subroutine readPackedParMatrix (datasetID,packedMatrix,matrixDims,dim1,dim2)
   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer :: dim1
   integer :: dim2
   integer (hid_t) :: datasetID
!   real (kind=double), dimension (dim1,dim2*(dim2+1)/2) :: packedMatrix
   real (kind=double), dimension (dim1,dim2) :: packedMatrix
   integer (hsize_t), dimension (2) :: matrixDims

   ! Define local variables
   integer :: hdferr

   integer :: i,j

   ! Read the matrix directly.
   call h5dread_f (datasetID,H5T_NATIVE_DOUBLE,packedMatrix,matrixDims,hdferr)

   if (hdferr /= 0) stop 'Failed to read packed matrix'

end subroutine readPackedParMatrix

subroutine unpackSerMatrix (matrix,packedMatrix,dim2,fullFlag)
   ! Make sure that there are not accidental variable declarations.
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
end subroutine unpackSerMatrix

subroutine readPackedSerMatrix (datasetID,packedMatrix,matrixDims,dim1,dim2)
   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer :: dim1
   integer :: dim2
   integer (hid_t) :: datasetID
   real (kind=double), dimension (dim1,dim2*(dim2+1)/2) :: packedMatrix
   integer (hsize_t), dimension (2) :: matrixDims

   ! Define local variables
   integer :: hdferr

   ! Read the matrix directly.
   call h5dread_f (datasetID,H5T_NATIVE_DOUBLE,packedMatrix,matrixDims,hdferr)
   if (hdferr /= 0) stop 'Failed to read packed matrix'

end subroutine readPackedSerMatrix

end module matSubs
