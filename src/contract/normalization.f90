module NormalizationSubs

public

contains

subroutine computeAngNorm

   ! Use definition modules.
   use O_Kinds
   use O_Constants

   ! Use data modules.
   use GaussianBasisData

   implicit none

   ! Define local variables.
   integer :: i

   ! Allocate space to hold the angular normalization coefficients where the
   !   second is simply the same as the first with a factor of 1/sqrt(4pi) *
   !   sqrt(1,3,15,105) for s,p,d,f (*1/2 for f only).
   allocate (angNorm1(4,maxNumBasisGaussians))
   allocate (angNorm2(4,maxNumBasisGaussians))

   ! Assign the angular normalization coefficients for each s,p,d,f orbitals
   !   for each basis gaussian.  The angular normalization coefficient for a
   !   gaussian function depends on the exponential alpha of the gaussian and
   !   the angular momentum quantum number.  angNorm1 consists of the first
   !   two terms while angNorm2 factors in the third term.
   ! S = 2                      * (8*alpha^3/pi)^0.25   *  sqrt(1/4pi)
   ! P = sqrt(8/3)              * (32*alpha^5/pi)^0.25  *  sqrt(3/4pi)
   ! D = (8*alpha/sqrt(15))     * (8*alpha^3/pi)^0.25   *  sqrt(15/4pi)
   ! F = (16*alpha^2/sqrt(105)) * (8*alpha/pi)^0.25     *  sqrt(105/4pi)

   ! S = (16 * alpha^3 / pi)^0.25 * sqrt(1/4pi)
   ! P = (2048/9 * alpha^5 / pi)^0.25 * sqrt(3/4pi)
   ! D = (32768/225 * alpha^8 / pi)^0.25 * sqrt(15/4pi)
   ! F = (524288/11025 * alpha^9 / pi)^0.25 * sqrt(105/4pi)
   do i = 1, maxNumBasisGaussians
      angNorm1(1,i) = (2.0_double * &
            & (8.0_double * basisAlphas(i)**3 / pi)**0.25_double)
      angNorm1(2,i) = (sqrt(8.0_double/3.0_double) * &
            & (32.0_double * basisAlphas(i)**5 / pi)**0.25_double)
      angNorm1(3,i) = (8.0_double * basisAlphas(i) / sqrt(15.0_double) * &
            & (8.0_double * basisAlphas(i)**3 / pi)**0.25_double)
      angNorm1(4,i) = (16.0_double * basisAlphas(i)**2 / sqrt(105.0_double)*&
            & (8.0_double * basisAlphas(i) / pi)**0.25_double)

      angNorm2(1,i) = angNorm1(1,i) * sqrt(1.0_double / pi / 4.0_double)
      angNorm2(2,i) = angNorm1(2,i) * sqrt(3.0_double / pi / 4.0_double)
      angNorm2(3,i) = angNorm1(3,i) * sqrt(15.0_double / pi / 4.0_double)
      angNorm2(4,i) = angNorm1(4,i) * sqrt(105.0_double / pi / 4.0_double)
   enddo


end subroutine computeAngNorm

subroutine renormWaveFns

   ! Use definition modules.
   use O_Kinds
   use O_Constants

   ! Use data modules.
   use AtomData
   use GaussianBasisData
   use MatrixElementData  ! At this point the hamiltonianME hold waveFns.

   implicit none

   ! Define local variables.
   integer :: i,j,k,l
   real (kind=double) :: tempCoeff
   real (kind=double) :: coeff
   real (kind=double) :: factor
   real (kind=double) :: alphaSum

   ! Allocate space to hold the normalization constants for each gaussian term
   !   of each orbital.
!   allocate (normCoeff(max(numTotalOrbitals(:)),4))

   do i = 1, 4  ! The four angular momentum QNs.
      do j = 1, numTotalOrbitals(i)

         ! Initialize the normalization coefficient for this orbital.
         coeff = 0.0_double

         do k = 1, maxNumBasisGaussians
            do l = 1, k

               ! Since results are symmetric we will need to modify the coeff
               !   for non-double counted terms.
               if (k == l) then
                  factor = 1.0_double
               else
                  factor = 2.0_double
               endif

               ! Precompute the sum of basis gaussian exponential alphas.
               alphaSum = basisAlphas(k) + basisAlphas(l)

               ! Compute the integration coefficient.
               select case (i)
                  case (1)
                     tempCoeff = sqrt(pi)/2.0_double/(2.0_double*alphaSum*&
                           & sqrt(alphaSum))*angNorm1(1,k)*angNorm1(1,l)
                  case (2)
                     tempCoeff = sqrt(pi)*3.0_double/4.0_double/(2.0_double*&
                           & alphaSum*alphaSum*sqrt(alphaSum))*&
                           & angNorm1(2,k)*angNorm1(2,l)
                  case (3)
                     tempCoeff = sqrt(pi)*15.0_double/8.0_double/(2.0_double*&
                           & alphaSum*alphaSum*alphaSum*sqrt(alphaSum))*&
                           & angNorm1(3,k)*angNorm1(3,l)
                  case (4)
                     tempCoeff = sqrt(pi)*105.0_double/16.0_double/(2.0_double*&
                           & alphaSum*alphaSum*alphaSum*alphaSum*&
                           & sqrt(alphaSum))*angNorm1(4,k)*angNorm1(4,l)
                  case default
                     write(6,*) "normalization.f90 subroutine: renormWaveFns"
                     write(6,*) "Only cases 1,2,3,4 for s,p,d,f are implemented"
               end select

               ! Accumulate the final renormalization coefficient.
               coeff = coeff + tempCoeff * factor * &
                     & hamiltonianME(k,j,i)*hamiltonianME(l,j,i)
            enddo
         enddo

         ! Write the normalization coefficient for this orbital.
         write (20,*) "QN_l=",i-1," orbIndex=",j," Norm=",coeff

         ! Apply the normalization coefficient.
         do k = 1, numBasisGaussians(i)
            hamiltonianME(k,j,i) = angNorm1(i,k)*hamiltonianME(k,j,i) / &
                  & sqrt(coeff)
         enddo
      enddo
   enddo
end subroutine renormWaveFns

end module NormalizationSubs
