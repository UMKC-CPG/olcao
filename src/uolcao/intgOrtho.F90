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
