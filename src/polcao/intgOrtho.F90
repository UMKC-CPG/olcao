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

! It would be nice if this subroutine didn't have to allocate a third
! global array in order to do it's matrix multiplication. Consider
! creating a global arrays patch which allows for the option of
! a complex conjugate transpose. Or better yet, get pblas working so that
! we can cut this down considerably.
subroutine valeCoreCoreValeOL (valeDim,coreDim,numKPoints,&
  & ga_valeVale,ga_coreValeOL,whichKP)

   ! Import the necessary modules
   use O_Kinds
   use MPI

   ! Import the necessary interface.
!   use zherkInterface

   ! Make sure no variables are accidently declared
   implicit none
#include "mafdecls.fh"
#include "global.fh"

   ! Define passed parameters
   integer, intent(in) :: whichKP
   integer, intent(inout) :: ga_valeVale, ga_coreValeOL
   integer, intent(in) :: numKPoints, coreDim, valeDim
   integer, dimension(3) :: numblocks, blockdims
   integer, dimension(3) :: alo,ahi,blo,bhi,clo,chi

!   alo=(/1,1,whichKP/)
!   ahi=(/valeDim,coreDim,whichKP/)
!   blo=(/1,1,whichKP/)
!   bhi=(/coreDim,valeDim,whichKP/)
!   clo=(/1,1,whichKP/)
!   chi=(/valeDim,valeDim,whichKP/)
!
!   call nga_matmul_patch('C','N',(-2.0_double,0.0_double),&
!    & (1.0_double,0.0_double),&
!    & ga_coreValeOL,alo,ahi, &
!    & ga_coreValeOL,blo,bhi, &
!    & ga_valeVale,  clo,chi)

   call ga_zgemm('C','N',valeDim,valeDim,coreDim,(-2.0_double,0.0_double),&
     & ga_coreValeOL,ga_coreValeOL,(1.0_double,0.0_double),ga_valeVale)

!    print *, "after ValeCoreCoreValeOL"
!    call ga_print(ga_valeVale)
!   call ga_zgemm('C','N', 64,64,8,&
!    & (-2.0_double,1.0_double),ga_conjCVOL,ga_tempCVOL, &
!    & (1.0_double,1.0_double),ga_tempVV)
end subroutine valeCoreCoreValeOL

subroutine valeCoreCoreVale (valeDim,coreDim,ga_valeVale,ga_coreVale,&
   & ga_coreValeOL, whichkp)

   ! Import the necessary modules
   use O_Kinds
   use MPI

   ! Import the necessary interface.
!   use zher2kInterface

   ! Make sure no variables are accidently declared
   implicit none
#include "mafdecls.fh"
#include "global.fh"

   ! Define passed parameters
   integer, intent(in) :: valeDim, coreDim
   integer, intent(in) :: whichKP
   integer, intent(inout) :: ga_valeVale
   integer, intent(inout) :: ga_coreVale
   integer, intent(inout) :: ga_coreValeOL
   integer, dimension (3) :: alo,ahi,blo,bhi,clo,chi

!   alo=(/1,1,whichKP/)
!   ahi=(/valeDim,coreDim,whichKP/)
!   blo=(/1,1,whichKP/)
!   bhi=(/coreDim,valeDim,whichKP/)
!   clo=(/1,1,whichKP/)
!   chi=(/valeDim,valeDim,whichKP/)
!   
!   call nga_matmul_patch('C','N',(-1.0_double,0.0_double),&
!    & (1.0_double,0.0_double), &
!    & ga_coreValeOL,alo,ahi, &
!    & ga_coreVale,blo,bhi, &
!    & ga_valeVale,clo,chi)
   
   call ga_zgemm('C','N',valeDim,valeDim,coreDim,(-1.0_double,0.0_double),&
     & ga_coreValeOL, ga_coreVale,(1.0_double,0.0_double),ga_valeVale)
   
!   call nga_matmul_patch('C','N',(-1.0_double,0.0_double),&
!    & (1.0_double,0.0_double), &
!    & ga_coreVale,alo,ahi, &
!    & ga_coreValeOL,blo,bhi, &
!    & ga_valeVale,clo,chi)

   call ga_zgemm('C','N',valeDim,valeDim,coreDim,(-1.0_double,0.0_double),&
     & ga_coreVale,ga_coreValeOL,(1.0_double,0.0_double),ga_valeVale)

!   call ga_zgemm('C','N', valeDim,valeDim,coreDim,&
!    & cmplx(-1.0_double,0.0_double),ga_conjCVOL,ga_tempCVOL, &
!    & 1.0_double,ga_tempVV)
!
!   call ga_zgemm('C','N', valeDim,valeDim,coreDim,&
!    & cmplx(-1.0_double,-0.0_double),ga_conjCV,ga_tempCVOL, &
!    & 1.0_double,ga_tempVV)

!   call zher2k('U','C',valeDim,coreDim,(-1.0_double,0.0_double),&
!         & coreValeOL,coreDim,coreVale,coreDim,1.0_double,valeVale,valeDim)

!   valeVale = valeVale - matmul(conjg(transpose(coreValeOL)),coreVale) - &
!         & matmul(conjg(transpose(coreVale)),coreValeOL)

!   call zherk('U','C',valeDim,coreDim,-2.0_double,coreValeOL,coreDim,&
!         & 1.0_double,valeVale,valeDim)
   ! Create temp arrays to do matrix routines
end subroutine valeCoreCoreVale

subroutine coreValeCoreCore (valeDim,coreDim,ga_valeCore,ga_coreVale,&
  & ga_coreCore, whichKP)

   ! Import the necessary modules
   use O_Kinds
   use MPI

   ! Import the necessary blas interface
!   use zhemmInterface

   ! Make sure no variables are accidently declared
   implicit none
#include "mafdecls.fh"
#include "global.fh"

   ! Define passed parameters
   integer, intent(in) :: valeDim, coreDim, whichKP
   integer, intent(inout) :: ga_valeCore
   integer, intent(inout) :: ga_coreVale,ga_coreCore
   integer, dimension (3) :: alo,ahi,blo,bhi,clo,chi

!   alo=(/1,1,whichKP/)
!   ahi=(/valeDim,coreDim,whichKP/)
!   blo=(/1,1,whichKP/)
!   bhi=(/coreDim,coreDim,whichKP/)
!   clo=(/1,1,whichKP/)
!   chi=(/valeDim,coreDim,whichKP/)
!   
!   call nga_matmul_patch('C','N',(1.0_double,0.0_double),&
!    & (1.0_double,0.0_double), &
!    & ga_coreVale,alo,ahi, &
!    & ga_coreCore,blo,bhi, &
!    & ga_valeCore,clo,chi)

   call ga_zgemm('C','N',valeDim,coreDim,coreDim,(1.0_double,0.0_double),&
     & ga_coreVale, ga_coreCore,(1.0_double,0.0_double),ga_valeCore)

!    print *, "after coreValeCoreCore"
!    call ga_print(ga_valeCore)
!   call ga_zgemm('C','N', valeDim,coreDim,coreDim, &
!    & 1.0_double,ga_tempCV,ga_tempCC, 1.0_double, ga_valeCore)

!   call zhemm ('R','U',valeDim,coreDim,(1.0_double,0.0_double),&
!         & coreCore,coreDim,coreVale,coreDim,(1.0_double,0.0_double),&
!         & valeCore,valeDim)

  ! do i = 1,coreDim
  !    do j = 1, valeDim
  !       valeCore(i,j) = sum(conjg(coreVale(1:coreDim,j)) * &
  !            & coreCore(1:coreDim,i))
  !    enddo
  ! enddo

!   valeCore = matmul(conjg(transpose(coreVale)),coreCore)
end subroutine coreValeCoreCore


subroutine makeValeVale (valeDim,coreDim,ga_valeCore, &
      & ga_coreVale, ga_valeVale,whichKP)

   ! Import the necessary modules
   use O_Kinds
   use MPI

   ! Make sure no variables are accidently declared
   implicit none
#include "mafdecls.fh"
#include "global.fh"

   ! Define passed parameters
   integer, intent(in) :: valeDim, coreDim, whichKP
   ! Define loop control variables
   integer, intent(inout) :: ga_valeCore
   integer, intent(inout) :: ga_coreVale, ga_valeVale

   integer, dimension (3) :: alo,ahi,blo,bhi,clo,chi

   ! Define the small threshhold for eliminating resultant values
   real (kind=double) :: smallThresh10
   smallThresh10 = real(1.0D-10,double)
   
!   alo=(/1,1,whichKP/)
!   ahi=(/valeDim,coreDim,whichKP/)
!   blo=(/1,1,whichKP/)
!   bhi=(/coreDim,valeDim,whichKP/)
!   clo=(/1,1,whichKP/)
!   chi=(/valeDim,valeDim,whichKP/)
!   
!   call nga_matmul_patch('N','N',(1.0_double,0.0_double),&
!    & (1.0_double,0.0_double), &
!    & ga_valeCore,alo,ahi, &
!    & ga_coreVale,blo,bhi, &
!    & ga_valeVale,clo,chi)
    call ga_zgemm('N','N',valeDim,valeDim,coreDim,(1.0_double,0.0_double),&
      & ga_valeCore,ga_coreVale,(1.0_double,0.0_double),ga_valeVale)

!   print *, "after makeValeVale"
!   call ga_print(ga_valeVale)
   
!   call ga_zgemm('N','N',valeDim,valeDim,coreDim, &
!    & 1.0_double,ga_valeCore,ga_tempCV,1.0_double, ga_tempVV)

!   if (storeFlag == 0) then
!      if (fullFlag == 0) then
!         do i = 1, valeDim
!            do j = 1, i
!               valeVale(j,i) = valevale(j,i) + sum( &
!                     & valeCore(1:coreDim,j) * coreVale(1:coreDim,i))
!            enddo
!         enddo
!      else
!         do i = 1, valeDim
!            do j = 1, i
!               valeVale(j,i) = valevale(j,i) + sum( &
!                     & valeCore(1:coreDim,j) * coreVale(1:coreDim,i))
!               valeVale(i,j) = conjg(valeVale(j,i))
!            enddo
!         enddo
!      endif
!   else
!
!      ! Initialize the current index for packing the valeVale matrix.
!      currentIndex = 0
!
!      do i = 1, valeDim
!         do j = 1, i
!
!            ! Compute the result.
!            valeVale(j,i) = valeVale(j,i) + sum( &
!                  & valeCore(1:coreDim,j) * coreVale(1:coreDim,i))
!
!            ! Increment the array index for packing the valeVale matrix.
!            currentIndex = currentIndex + 1
!
!            ! Pack it.
!            packedValeVale(1,currentIndex) =  real(valeVale(j,i),double)
!            packedValeVale(2,currentIndex) = aimag(valeVale(j,i))
!
!            ! Check for negligable components and reduce them to 0.
!            if (abs(packedValeVale(1,currentIndex)) < smallThresh10) then
!               packedValeVale(1,currentIndex) = 0.0_double
!            endif
!            if (abs(packedValeVale(2,currentIndex)) < smallThresh10) then
!               packedValeVale(2,currentIndex) = 0.0_double
!            endif
!         enddo
!      enddo
!   endif
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
