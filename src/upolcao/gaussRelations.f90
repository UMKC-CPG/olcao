module O_GaussianRelations

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Atomic elemental data.
   integer, allocatable, dimension (:) :: numElementAtomicAlphas
   real (kind=double), allocatable, dimension (:,:) :: elementAtomicAlphas

   ! Electronic potential elemental data.
   integer, allocatable, dimension (:) :: numElementPotAlphas
   real (kind=double), allocatable, dimension (:,:) :: elementPotAlphas

   ! Nuclear potential elemental data
   real (kind=double), allocatable, dimension (:) :: elementNucAlpha


   ! Large dimension matrices holding interaction results.
   real (kind=double), allocatable, dimension (:,:,:,:)     :: alphaDist
   real (kind=double), allocatable, dimension (:,:,:,:)     :: alphaCenter
   real (kind=double), allocatable, dimension (:,:,:,:,:)   :: alphaNucDist
   real (kind=double), allocatable, dimension (:,:,:,:,:,:) :: alphaPotDist

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

! Recall that the product of two gaussians is also a gaussian.  Our goal
!   is to solve the problem of determining how far apart two gaussians need
!   to be to make the peak height of the product of the gaussians less
!   than some threshold.  At that distance we consider the two gaussians
!   in the product to be non-interacting.  We will then use this distance
!   to determine if two atoms will interact or not.  We calculate the
!   magnitude of the distance between the two atoms squared and compare it
!   to the non-interacting distance squared.  The reason that we compare
!   the squared values is so that we can avoid taking a square root in
!   both cases because it is computational costly and if the answer is
!   independent of whether we square both sides or not.
! Consider the product of two gaussians:
!   exp(-a*(r-A) * exp(-b*(r-B)^2) = Kexp(-c*(r-C)^2)
!   = exp((-ab(A-B)^2) / (a+b)) * exp(-(a+b)(x-((aA+bB)/(a+b)))^2)
!   So we have:  K = exp((-ab(A-B)^2) / (a+b))
!                c = a+b
!                C = ((aA+bB) / (a+b))
! The derivation is pretty easy and can be easily searched for online too.
!   This process is then repeated for the three center integrals where we take
!   the already combined gaussian and multiply it with another.  Again, the
!   goal is to find the cut off ranges etc.  To speed up the calculation we
!   pre-compute and store the results for the various combinations of alphas.
!   We can do this because the alphas vary only by the different elements in
!   the system and often there are only 5 or so different elements.

subroutine makeAlphaDist

   ! Import the necessary modules
   use O_Kinds
   use O_Lattice, only: logBasisFnThresh
   use O_AtomicTypes, only: atomTypes, numAtomTypes, numElements, &
         & maxNumAtomAlphas

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables
   integer :: i,j,k,l
   integer :: currElementID


   ! Allocate space to hold the number of alphas for each element, and the
   !   actual alphas for each element.
   allocate (numElementAtomicAlphas(numElements))
   allocate (elementAtomicAlphas(maxNumAtomAlphas,numElements))

   ! Initialize the above arrays.
   numElementAtomicAlphas (:)   = 0
   elementAtomicAlphas    (:,:) = 0.0_double

   ! Determine the number of atomic alphas for each element and record the
   !   alphas for that element.
   do i = 1, numAtomTypes

      currElementID = atomTypes(i)%elementID

      ! Check if this element has been assigned yet.
      if (numElementAtomicAlphas(currElementID) == 0) then

         ! Record the number of alphas for this element.
         numElementAtomicAlphas(currElementID) = (atomTypes(i)%numOrbAlphas(1))

         ! Record the alphas for this element.
         elementAtomicAlphas(:numElementAtomicAlphas(currElementID), &
               & currElementID) = atomTypes(i)%alphas(:)
     endif
   enddo

   ! Allocate space to hold the alpha distance matrix and alpha center matrix.
   allocate (alphaDist(maxNumAtomAlphas,maxNumAtomAlphas,numElements,&
         & numElements))
   allocate (alphaCenter(maxNumAtomAlphas,maxNumAtomAlphas,numElements,&
         & numElements))

   ! Initialize the alphaDist and alphaCenter matrices to 0.
   alphaDist   (:,:,:,:) = 0.0_double
   alphaCenter (:,:,:,:) = 0.0_double

   ! Loop over the element combinations to make the alpha distance matrix.
   do i = 1, numElements
   do j = 1, numElements
      do k = 1, numElementAtomicAlphas(i)
      do l = 1, numElementAtomicAlphas(j)

         alphaCenter(l,k,j,i) = elementAtomicAlphas(l,j) / &
               & (elementAtomicAlphas(l,j) + elementAtomicAlphas(k,i))

         ! Produce a quantity like: a*(alpha1 + alpha2) / (alpha1*alpha2) which
         !   is a Gaussian product theorem prefactor.
         alphaDist(l,k,j,i) = logBasisFnThresh / &
               & elementAtomicAlphas(k,i) / alphaCenter(l,k,j,i)
      enddo
      enddo
   enddo
   enddo

end subroutine makeAlphaDist

subroutine makeAlphaPotDist

   ! Import the necessary modules
   use O_Kinds
   use O_Lattice, only: logBasisFnThresh
   use O_PotTypes, only: maxNumPotAlphas, numPotTypes, potTypes
   use O_AtomicTypes, only: numElements, maxNumAtomAlphas

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables
   integer :: i,j,k,l,m,n
   integer :: currElementID


   ! Allocate space to hold the number of alphas for each element, and the
   !   actual alphas for each element.
   allocate (numElementPotAlphas(numElements))
   allocate (elementPotAlphas(maxNumPotAlphas,numElements))

   ! Initialize the above arrays.
   numElementPotAlphas (:)   = 0
   elementPotAlphas    (:,:) = 0.0_double

   ! Determine the number of potential alphas for each element and record the
   !   alphas for that element.
   do i = 1, numPotTypes

      currElementID = potTypes(i)%elementID

      ! Check if this element has been assigned yet.
      if (numElementPotAlphas(currElementID) == 0) then

         ! Record the number of alphas for this element.
         numElementPotAlphas(currElementID) = (potTypes(i)%numAlphas)

         ! Record the alphas for this element.
         elementPotAlphas(:numElementPotAlphas(currElementID), &
               & currElementID) = potTypes(i)%alphas(:)
     endif
   enddo

   ! Allocate space to hold the alpha distance matrix and alpha center matrix.
   allocate (alphaPotDist(maxNumAtomAlphas,maxNumAtomAlphas,maxNumPotAlphas,&
         & numElements,numElements,numElements))

   ! Initialize the alphaPotDist matrix to 0.
   alphaPotDist   (:,:,:,:,:,:) = 0.0_double

   ! Loop over the element combinations to make the alpha distance matrix.
   do i = 1, numElements
   do j = 1, numElements
   do k = 1, numElements
      do l = 1, numElementPotAlphas(i)
      do m = 1, numElementAtomicAlphas(j)
      do n = 1, numElementAtomicAlphas(k)
         alphaPotDist(n,m,l,k,j,i) = logBasisFnThresh * &
               & (elementAtomicAlphas(n,k) + elementAtomicAlphas(m,j) + &
               & elementPotAlphas(l,i)) / (elementAtomicAlphas(n,k) + &
               & elementAtomicAlphas(m,j)) / elementPotAlphas(l,i)
      enddo
      enddo
      enddo
   enddo
   enddo
   enddo

   ! Deallocate arrays that are no longer necessary
   deallocate (numElementPotAlphas)
   deallocate (elementPotAlphas)
   deallocate (numElementAtomicAlphas)
   deallocate (elementAtomicAlphas)

end subroutine makeAlphaPotDist


subroutine makeAlphaNucDist

   ! Import the necessary modules
   use O_Kinds
   use O_Lattice, only: logBasisFnThresh
   use O_PotTypes, only: numPotTypes, potTypes
   use O_AtomicTypes, only: numElements, maxNumAtomAlphas

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables
   integer :: i,j,k,l,m
   integer :: currElementID


   ! Allocate space to hold the number of alphas for each element, and the
   !   actual alphas for each element.
   allocate (elementNucAlpha(numElements))

   ! Initialize the above arrays.
   elementNucAlpha (:) = 0.0_double

   ! Determine the number of potential alphas for each element and record the
   !   alphas for that element.
   do i = 1, numPotTypes

      currElementID = potTypes(i)%elementID

      ! Check if this element has been assigned yet.
      if (elementNucAlpha(currElementID) == 0) then

         ! Record the nuclear alpha for this element.
         elementNucAlpha(currElementID) = (potTypes(i)%nucAlpha)
     endif
   enddo

   ! Allocate space to hold the alpha distance matrix and alpha center matrix.
   allocate (alphaNucDist(maxNumAtomAlphas,maxNumAtomAlphas,numElements,&
         & numElements,numElements))

   ! Initialize the alphaPotDist matrix to 0.
   alphaNucDist   (:,:,:,:,:) = 0.0_double

   ! Loop over the element combinations to make the alpha distance matrix.
   do i = 1, numElements
   do j = 1, numElements
   do k = 1, numElements
      do l = 1, numElementAtomicAlphas(i)
      do m = 1, numElementAtomicAlphas(j)

         alphaNucDist(m,l,k,j,i) = logBasisFnThresh * &
               & (elementAtomicAlphas(m,j) + elementAtomicAlphas(l,i) + &
               & elementNucAlpha(k)) / (elementAtomicAlphas(m,j) + &
               & elementAtomicAlphas(l,i)) / elementNucAlpha(k)
      enddo
      enddo
   enddo
   enddo
   enddo

   ! Deallocate arrays that are no longer necessary
   deallocate (elementNucAlpha)

end subroutine makeAlphaNucDist


subroutine cleanUpGaussRelations

   implicit none

   deallocate (alphaDist)
   deallocate (alphaNucDist)
   deallocate (alphaPotDist)
   deallocate (alphaCenter)

end subroutine cleanUpGaussRelations


end module O_GaussianRelations
