module ExecutionData

   ! Import definition modules.
   use O_Kinds

   integer :: doCharge
   integer :: doVis
end module ExecutionData



module AtomData

   ! Import definition modules.
   use O_Kinds

   character*2 :: elementName
   integer, dimension(4)   :: numCoreOrbitals
   integer, dimension(4)   :: numValeOrbitals
   integer, dimension(4)   :: numTotalOrbitals
   integer, dimension(4,3) :: numValeOrbitalsPerBasis
   real (kind=double) :: atomicNumber
   real (kind=double) :: nuclearAlpha
end module AtomData



module GaussianBasisData

   ! Import definition modules.
   use O_Kinds

   integer :: maxNumBasisGaussians
   integer, dimension (4) :: numBasisGaussians ! Num gaussians for spdf type.
   integer, dimension (4) :: numBasisTerms ! Includes imbedded coeff=0 terms.
   integer, allocatable, dimension(:,:) :: selectGaussians
         ! Dim1 = maxNumBasisGaussian; Dim2 = 4 for spdf.
   integer, allocatable, dimension(:,:) :: basisID ! Dim1 = numValeOrbitals for
         ! a given s,p,d,f.  Dim2 = s,p,d,f
   real (kind=double) :: minBasisAlpha
   real (kind=double) :: maxBasisAlpha
   real (kind=double), allocatable, dimension (:) :: basisAlphas
   real (kind=double), allocatable, dimension (:,:) :: angNorm1
   real (kind=double), allocatable, dimension (:,:) :: angNorm2
end module GaussianBasisData



module PotGaussianData

   ! Import definition modules.
   use O_Kinds

   integer :: numPotGaussians
   real (kind=double), allocatable, dimension (:) :: potAlphas
   real (kind=double), allocatable, dimension (:) :: potCoeffs
end module PotGaussianData




module MatrixElementData

   ! Import definition modules.
   use O_Kinds

   real (kind=double), allocatable, dimension (:,:,:) :: overlapME ! Indicies 1
         ! and 2 are 1:numBasisGaussians, and index 3 is 1:4 for the s,p,d,f
         ! orbitals.
   real (kind=double), allocatable, dimension (:,:,:) :: hamiltonianME !
         ! Indicies 1 and 2 are 1:numBasisGaussians, and index 3 is 1:4 for the
         ! s,p,d,f orbitals.  The first index is used for the number of basis
         ! gaussians and the second is used for the number of orbitals of the
         ! current angular momentum QN when the wave function coefficients
         ! are saved here after diagonalization.
end module MatrixElementData



module WaveFunctionData

   ! Import definition modules.
   use O_Kinds

   real (kind=double), allocatable, dimension (:,:) :: eigenValues ! Index 1
         ! is the max number of basis functions and index 2 is the angMomQN.

   ! Annoying, but note that the eigen vectors are stored in hamiltonianME and
   !   I did not bother to copy them to another named array, although I should
   !   probably just make a pointer here that points to it.  It would be slower
   !   to evaluate but more clear, and the program only takes < 2 seconds to
   !   run anyway.  If I have time and the reason I will do it later.

end module WaveFunctionData
