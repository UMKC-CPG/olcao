module ChargeDensityMod

use O_Kinds

public

real (kind=double) :: negligLimit ! Negligability limit between two gaussian
      ! functions in the basis for determining the contribution to the charge.

real (kind=double), dimension (2) :: accumCharge ! Accumulate the core and
      ! valence contributions for each radial point for each orbital of each
      ! angular momentum quantum number.

real (kind=double), allocatable, dimension (:,:) :: radialCharge ! Index 1
      ! is 1:4 for radial value, core, valence, and total charge.  Index 2 is
      ! 1:numRadialPoints.
real (kind=double), allocatable, dimension (:,:) :: orbitalCharge ! Index 1
      ! is 1:4 (max num l_QN in program), Index 2 is from 1:numTotalOrbitals(1)
      ! (Number of 1s orbitals because 1s is the most numerous orbital type.)
      ! This holds the charge per orbital.

contains

subroutine makeRho

   ! Use encessary definitions.
   use O_Kinds
   use O_Constants

   ! Use necessary data.
   use AtomData
   use GaussianBasisData
   use MatrixElementData

   ! Use necessary modules.
   use O_RadialGrid

   implicit none

   ! Define local variables.
   integer :: i,j,k,l,m
   real (kind=double) :: expTerm
   real (kind=double), allocatable, dimension (:,:) :: factor ! Used to double
         ! count cross interactions and singly count self interactions between
         ! Gaussian terms.

   ! Allocate space to hold the results of the computation and the factor used
   !   to help comput it.
   allocate (radialCharge(4,numRadialPoints))
   allocate (factor(maxNumBasisGaussians,maxNumBasisGaussians))

   ! Initialize all the radial charge values to zero.
   radialCharge(:,:) = 0.0_double

   ! Initialize all cross terms to 2 and all self terms to 1.
   factor(:,:) = 2.0_double
   do i = 1, maxNumBasisGaussians
      factor(i,i) = 1.0_double
   enddo

   ! Initialize the value of the negligability limit.
   negligLimit = radialMaxDist

   do i = 1, numRadialPoints

      ! Save the radial position for easy printing later.
      radialCharge(1,i) = radialPoints(i)

      ! Loop over the max number of angular momentum QNs.
      do j = 1, 4 ! s,p,d,f

         ! Loop over the max of the core or valence orbitals for this angmom
         !   QN.  In this way we treat both orbital sets as efficiently as
         !   possible.  (I think).
         do k = 1,max(numCoreOrbitals(j),numValeOrbitals(j))

            ! Initialize an accumulator for the charge contributed by the
            !   basis gaussians of this orbital.  Index 1=core, Index 2=vale.
            accumCharge(:) = 0.0_double

            do l = 1, numBasisGaussians(j)
               do m = 1, l

                  ! Compute the exponential alpha of the product of two
                  !   basis gaussian functions.
                  expTerm = (basisAlphas(l)+basisAlphas(m)) * &
                        & radialPoints(i)**2

                  if (expTerm < negligLimit) then

                     if (k <= numCoreOrbitals(j)) then
                        accumCharge(1) = accumCharge(1) + &
                           & exp(-expTerm) * factor(m,l) * &
                           & hamiltonianME(l,k,j)*hamiltonianME(m,k,j)
                     endif

                     if (k <= numValeOrbitals(j)) then
                        accumCharge(2) = accumCharge(2) + &
                           & exp(-expTerm) * factor(m,l) * &
                           & hamiltonianME(l,k+numCoreOrbitals(j),j) * &
                           & hamiltonianME(m,k+numCoreOrbitals(j),j)
                     endif
                  endif
               enddo
            enddo

            ! Accumulate the charge contributed by this orbital to the total
            !   core and valence charge.  Note that we scale the accumCharge by
            !   the number of electrons present in the orbital.
            radialCharge(2,i) = radialCharge(2,i) + accumCharge(1) * &
                  & orbitalCharge(j,k)
            radialCharge(3,i) = radialCharge(3,i) + accumCharge(2) * &
                  & orbitalCharge(j,k+numCoreOrbitals(j))
         enddo
      enddo

      ! Multiply in a factor of 4*pi to the core and valence charge.
      radialCharge(2,i) = radialCharge(2,i) / 4.0_double / pi
      radialCharge(3,i) = radialCharge(3,i) / 4.0_double / pi

      ! Save the total charge for this radial point.
      radialCharge(4,i) = radialCharge(2,i) + radialCharge(3,i)
   enddo

end subroutine makeRho

subroutine printRho

   ! Use encessary definitions.
   use O_Kinds

   ! Use necessary subroutine modules.
   use O_WriteDataSubs

   ! Use necessary modules.
   use O_RadialGrid

   implicit none

   ! Define local variables.
   integer :: i
   integer :: fileUnit

   ! Define the file unit to use for output.
   fileUnit = 50

   ! Open the wave function file for printing.
   open(unit=fileUnit,form="formatted",file="charge.plot",status="new")

   do i = 1, numRadialPoints
      call writeData(4,radialCharge(:,i),fileUnit)
   enddo

   close (fileUnit)
end subroutine printRho

end module ChargeDensityMod
