module matrixElementSubs

   private
   public :: prepMatrixElements

contains

subroutine prepMatrixElements

   ! Import necessary definitions.
   use O_Kinds

   ! Import necessary data modules.
   use MatrixElementData
   use GaussianBasisData
   use PotGaussianData
   use AtomData

   implicit none

   ! Define local variables.
   integer :: i,j,k,l  ! Loop indexing variables.
   integer, dimension (4) :: basisCount1, basisCount2, doAngMom
   real (kind=double) :: alpha1, alpha2
   real (kind=double), dimension (4)   :: overlap !(spdf)
   real (kind=double), dimension (4,3) :: hamiltonian !(spdf,ke_nucpot_epot)



   ! Allocate space to hold the matrix elements.
   allocate (overlapME(maxNumBasisGaussians,maxNumBasisGaussians,4))
   allocate (hamiltonianME(maxNumBasisGaussians,maxNumBasisGaussians,4))

   ! Initialize both matrices to 0.
   overlapME(:,:,:) = 0.0_double
   hamiltonianME(:,:,:) = 0.0_double

   ! Initialize counters for the number of gaussians used for s,p,d,f orbitals
   !   in the outer basis gaussian loop.  Take care to include all terms with
   !   an implicit 0 coefficient.
   basisCount1(:) = 0

   do i = 1, maxNumBasisGaussians

      ! Increment the count of s,p,d,f basis gaussians if necessary.
      do j = 1,4
         if (selectGaussians(i,j) == 1) basisCount1(j) = basisCount1(j) + 1
      enddo

      ! Initialize counters for the number of gaussians used for s,p,d,f
      !   orbitals in the inner basis gaussian loop.
      basisCount2(:) = 0

      do j = 1, maxNumBasisGaussians

         ! Increment the count of s,p,d,f basis gaussians if necessary.
         do k = 1,4
            if (selectGaussians(j,k) == 1) basisCount2(k) = basisCount2(k) + 1
         enddo

         ! Set flags to do (1) or not do (0) s,p,d,f ang mom orbitals.
         doAngMom(:) = 0
         do l = 1,4
            if ((selectGaussians(i,l) == 1) .and. &
                  & (selectGaussians(j,l) == 1)) doAngMom(l) = 1
         enddo

         ! Initialize the hamiltonian terms for each orbital.
         hamiltonian(:,:) = 0.0_double

         ! Save short names for the current gaussian basis alphas.
         alpha1 = basisAlphas(i)
         alpha2 = basisAlphas(j)

         ! Compute overlap between gaussian basis fns of the same orbital type.
         call computeOverlap (alpha1,alpha2,overlap(:))

         ! Compute KE between gaussian basis fns of the same orbital type.
         call computeKineticEnergy (alpha1,alpha2,hamiltonian(:,1))

         ! Compute nuclear potential gaussian with basis gaussians integral.
         !   This is of the form 1/r * exp(-alpha * r^2).
         call computeNucPot (alpha1,alpha2,potAlphas(1),potCoeffs(1),&
               & hamiltonian(:,2))

         ! Compute electronic potential gaussians with basis gaussians integral.
         !   This is of the form exp(-alpha * r^2).
         do k = 2, numPotGaussians
            call computeElecPot (alpha1,alpha2,potAlphas(k),potCoeffs(k),&
                  & hamiltonian(:,3))
         enddo

         ! Collect the kinetic energy, nuclear potential, and electronic
         !   potential terms of the hamiltonian together.
         hamiltonian(:,1) = hamiltonian(:,1) + hamiltonian(:,2) + &
               & hamiltonian(:,3)

         ! Save the overlap and hamiltonian matrix element values for the
         !   current gaussian pair for each orbital angular momentum QN.
         do k = 1,4
            if (doAngMom(k) == 1) then
               overlapME(basisCount1(k),basisCount2(k),k) = &
                     & angNorm2(k,i)*angNorm2(k,j)*overlap(k)
               hamiltonianME(basisCount1(k),basisCount2(k),k) = &
                     & angNorm2(k,i)*angNorm2(k,j)*hamiltonian(k,1)
            endif
         enddo
      enddo
   enddo
end subroutine prepMatrixElements

! This computes the integrals of two Gaussian functions of the same angular
!   momentum character centered at the same site.  The angular momentum
!   character is indicated by the prefixed r^somepower, the a is the combined
!   exponential alpha of the two Gaussian functions.:
!(ss)=int(int(int((r^0exp(-a*r^2))*r^2*sin(theta),
!      r=0..inf),theta=0..pi),phi=0..2pi)
!(pp)=int(int(int((r^2*exp(-a*r^2))*r^2*sin(theta),
!      r=0..inf),theta=0..pi),phi=0..2pi)
!(dd)=int(int(int((r^4*exp(-a*r^2))*r^2*sin(theta),
!      r=0..inf),theta=0..pi),phi=0..2pi)
!(ff)=int(int(int((r^6*exp(-a*r^2))*r^2*sin(theta),
!      r=0..inf),theta=0..pi),phi=0..2pi)
!The resultant coefficients are NOT direction included here.  E.g. The answer
!   for dd is 15/4 * pi^(3/4) / a^(7/2) when done in full.  Here however, the
!   factor of 15 is not included.
subroutine computeOverlap (alpha1,alpha2,overlap)

   ! Import necessary definitions.
   use O_Kinds
   use O_Constants

   implicit none

   ! Define passed parameters.
   real (kind=double) :: alpha1
   real (kind=double) :: alpha2
   real (kind=double), dimension(4) :: overlap

   ! Define local variables.
   real (kind=double) :: alphaSum
   real (kind=double) :: alphaFactor

   alphaSum    = alpha1+alpha2
   alphaFactor = (sqrt(pi/alphaSum))**3

   overlap(1) = alphaFactor
   overlap(2) = alphaFactor/(2.0_double*alphaSum)
   overlap(3) = alphaFactor/(4.0_double*alphaSum*alphaSum)
   overlap(4) = alphaFactor/(8.0_double*alphaSum*alphaSum*alphaSum)

end subroutine computeOverlap


! This computes the integrals of two Gaussian functions of the same angular
!   momentum character centered at the same site with a -1/2 Laplacian
!   operator between them.  The angular momentum
!   character is indicated by the prefixed r^somepower, the a is the combined
!   exponential alpha of the two Gaussian functions.:
subroutine computeKineticEnergy (alpha1,alpha2,kineticEnergy)

   ! Import necessary definitions.
   use O_Kinds
   use O_Constants

   implicit none

   ! Define passed parameters.
   real (kind=double) :: alpha1
   real (kind=double) :: alpha2
   real (kind=double), dimension(4) :: kineticEnergy

   ! Define local variables.
   real (kind=double) :: alphaSum
   real (kind=double) :: alphaFactor1
   real (kind=double) :: alphaFactor2

   alphaSum     = alpha1+alpha2
   alphaFactor1 = alpha1*alpha2/alphaSum
   alphaFactor2 = (sqrt(pi/alphaSum))**3

   kineticEnergy(1) = 3.0_double * alphaFactor1 * alphaFactor2
   kineticEnergy(2) = 5.0_double * alphaFactor1 * alphaFactor2 / &
         & (2.0_double*alphaSum)
   kineticEnergy(3) = 7.0_double * alphaFactor1 * alphaFactor2 / &
         & (4.0_double*alphaSum*alphaSum)
   kineticEnergy(4) = 9.0_double * alphaFactor1 * alphaFactor2 / &
         & (8.0_double*alphaSum*alphaSum*alphaSum)

end subroutine computeKineticEnergy


subroutine computeNucPot (alpha1,alpha2,nucPotAlpha,nucPotCoeff,nucPot)

   ! Import necessary definitions.
   use O_Kinds
   use O_Constants

   implicit none

   ! Define passed parameters.
   real (kind=double) :: alpha1
   real (kind=double) :: alpha2
   real (kind=double) :: nucPotAlpha
   real (kind=double) :: nucPotCoeff
   real (kind=double), dimension(4) :: nucPot

   ! Define local variables.
   real (kind=double) :: alphaSum
   real (kind=double) :: alphaFactor
   
   alphaSum    = alpha1+alpha2+nucPotAlpha
   alphaFactor = 2.0_double * pi/alphaSum

   nucPot(1) = nucPotCoeff*alphaFactor
   nucPot(2) = nucPotCoeff*alphaFactor/(3.0_double*alphaSum)
   nucPot(3) = nucPotCoeff*alphaFactor/(7.5_double*alphaSum*alphaSum)
   nucPot(4) = nucPotCoeff*alphaFactor/(17.5_double*alphaSum*alphaSum*alphaSum)

end subroutine computeNucPot


subroutine computeElecPot (alpha1,alpha2,elecPotAlpha,elecPotCoeff,elecPot)

   ! Import necessary definitions.
   use O_Kinds
   use O_Constants

   implicit none

   ! Define passed parameters.
   real (kind=double) :: alpha1
   real (kind=double) :: alpha2
   real (kind=double) :: elecPotAlpha
   real (kind=double) :: elecPotCoeff
   real (kind=double), dimension(4) :: elecPot

   ! Define local variables.
   real (kind=double) :: alphaSum
   real (kind=double) :: alphaFactor

   alphaSum    = alpha1+alpha2+elecPotAlpha
   alphaFactor = (sqrt(pi/alphaSum))**3

   elecPot(1) = elecPot(1) + elecPotCoeff*alphaFactor
   elecPot(2) = elecPot(2) + elecPotCoeff*alphaFactor/(2.0_double*alphaSum)
   elecPot(3) = elecPot(3) + elecPotCoeff*alphaFactor/(4.0_double*alphaSum*&
         & alphaSum)
   elecPot(4) = elecPot(4) + elecPotCoeff*alphaFactor/(8.0_double*alphaSum*&
         & alphaSum*alphaSum)

end subroutine computeElecPot


end module MatrixElementSubs
