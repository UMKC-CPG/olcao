module O_ValeCharge

   ! Import necessary modules.
   use O_Kinds
   use O_Constants

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Charge fitting structures.  The second index for potRho and the only
   !   index for the others is reserved for spin.
   real (kind=double), allocatable, dimension (:,:) :: potRho
   real (kind=double), allocatable, dimension (:)   :: chargeDensityTrace
   real (kind=double), allocatable, dimension (:)   :: nucPotTrace
   real (kind=double), allocatable, dimension (:)   :: kineticEnergyTrace


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine makeValenceRho

   ! Import the necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_AtomicSites,     only: valeDim
   use O_Input,           only: numStates
   use O_KPoints,         only: numKPoints
   use O_Constants,       only: smallThresh
   use O_Potential,       only: spin, potDim, potCoeffs, numPlusUJAtoms
   use O_MainEVecHDF5,    only: valeStates, eigenVectors_did
   use O_Populate,        only: electronPopulation, cleanUpPopulation
   use O_SetupIntegralsHDF5, only: atomOverlap_did, atomKEOverlap_did, &
         & atomNucOverlap_did, atomPotOverlap_did, atomDims
#ifndef GAMMA
   use O_BLASZHER
   use O_SecularEquation, only: valeVale, cleanUpSecularEqn, energyEigenValues,&
         & update1UJ
   use O_MatrixSubs, only: readMatrix, readPackedMatrix, matrixElementMult, &
         & packMatrix
#else
   use O_BLASDSYR
   use O_SecularEquation, only: valeValeGamma, cleanUpSecularEqn, &
         & energyEigenValues, update1UJ
   use O_MatrixSubs, only: readMatrixGamma, readPackedMatrix, &
         & matrixElementMultGamma, packMatrixGamma
#endif

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the local variables used in this subroutine.
   integer :: i,j,k ! Loop index variables
   integer :: skipKP
   integer :: dim1
   integer :: energyLevelCounter
   real (kind=double) :: sumElecEnergy
   real (kind=double), allocatable, dimension (:)     :: tempDensity
   real (kind=double), allocatable, dimension (:)     :: electronEnergy
   real (kind=double), allocatable, dimension (:)     :: currentPopulation
   real (kind=double), allocatable, dimension (:,:)   :: tempRealValeVale
   real (kind=double), allocatable, dimension (:,:)   :: tempImagValeVale
   real (kind=double), allocatable, dimension (:,:)   :: packedValeVale
   real (kind=double), allocatable, dimension (:,:,:) :: packedValeValeSpin
   real (kind=double), allocatable, dimension (:,:,:) :: &
         & structuredElectronPopulation
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: valeValeRho
#else
   real    (kind=double), allocatable, dimension (:,:,:) :: valeValeRhoGamma
#endif

   ! Log the date and time that we start.
   call timeStampStart (17)

   ! Define whether the packed arrays have two rows (complex) or one (real).
   !   This also defines the size of some other small arrays.
#ifndef GAMMA
   dim1 = 2
#else
   dim1 = 1
#endif

   ! Allocate the main valence valence density matrix.
#ifndef GAMMA
   allocate (valeValeRho(valeDim,valeDim,spin)) ! Complex
#else
   allocate (valeValeRhoGamma(valeDim,valeDim,spin)) ! Real
#endif

   ! Allocate a temporary packed valeVale matrix with a spin component.
   allocate (packedValeValeSpin(dim1,valeDim*(valeDim+1)/2,spin))

   ! Allocate a temporary holder for the charge density for the case that the
   !   calculation is spin polarized.  This is needed when rewriting the charge
   !   in a up+down, up-down form.
   allocate (tempDensity(dim1))

   ! Allocate space for the valence charge density (as represented by a
   !   summation of atom centered Gaussian functions in the same way as the
   !   potential function except with different coefficients).
   allocate (potRho (potDim,spin)) ! This will be deallocated in the makeSCFPot
         ! subroutine since after its values are copied to a local array there
         ! it will not be needed again until here.

   ! Allocate space to hold the trace of the charge density, nuclear potential,
   !   and kinetic energy.
   allocate (chargeDensityTrace(spin))
   allocate (nucPotTrace(spin))
   allocate (kineticEnergyTrace(spin))

   ! Allocate space to hold the currentPopulation based on spin
   allocate (currentPopulation (spin))
   allocate (electronEnergy    (spin))
   allocate (structuredElectronPopulation (numStates,numKPoints,spin))

   ! Initialize variables from data modules.
   potRho(:,:)           = 0.0_double
   nucPotTrace(:)        = 0.0_double
   kineticEnergyTrace(:) = 0.0_double
   chargeDensityTrace(:) = 0.0_double

   ! Initialize local variables
   electronEnergy(:)     = 0.0_double


   ! Fill a matrix of electron populations from the electron population that
   !   was computed in populateLevels.  Note that electronPopulation is a one
   !   dimensional array that has some order, but is not sorted in the way
   !   that the energy eigen values were sorted.  Please read the comments in
   !   the populateLevels subroutine to understand the order.
   !   (You can also probably get it from the loop order here ;)
   energyLevelCounter=0
   do i = 1, numKPoints
      do j = 1, spin
         do k = 1, numStates
            energyLevelCounter = energyLevelCounter + 1
            structuredElectronPopulation (k,i,j) = &
                  & electronPopulation(energyLevelCounter)
         enddo
      enddo
   enddo

   do i = 1, numKPoints

      ! Initialize space to read the wave functions.  Note that this matrix is
      !   a double complex matric because the wave function data contains both
      !   real and imaginary parts and has no symmetry.
      ! Note that it is only necessary to go through the initialization action
      !   if there are more than 1 kpoint.  For the 1 kpoint case, the valeVale
      !   matrix was not changed so it can still be used here.  Also note
      !   that the gammaKPoint option will never enter the "if" block.
#ifndef GAMMA
      if (numKPoints > 1) then

         ! Skip any kpoints with a negligable contribution for each state.
         skipKP = 0
         do j = 1, numStates
            if (sum(abs(structuredElectronPopulation(j,i,:)))>smallThresh) then
               skipKP = 1
               exit
            endif
         enddo
         if (skipKP == 0) then
            cycle
         endif

         ! Read the computed wave function into the valeVale matrix.  I know
         !   that the name is misleading, but so far this seems to work out
         !   the best in terms of not having too many names or too many
         !   allocate deallocate calls.  (And also we need to do this for the
         !   case where there is only 1 kpoint because we never write the
         !   wave function to disk and so the valeVale *still* holds it from
         !   the diagonalization call.  Therefore we are stuck using the same
         !   name even in the >1 kpoint case.)
         allocate (tempRealValeVale(valeDim,numStates))
         allocate (tempImagValeVale(valeDim,numStates))
         do j = 1, spin
            call readMatrix (eigenVectors_did(:,i,j),&
                  & valeVale(:,:numStates,1,j),&
                  & tempRealValeVale(:,:numStates),&
                  & tempImagValeVale(:,:numStates),&
                  & valeStates,valeDim,numStates)
         enddo
         deallocate (tempRealValeVale)
         deallocate (tempImagValeVale)
      endif

      ! Initialize matrix to receive the valeVale density matrix (square of the
      !   wave function).
      valeValeRho(:,:,:) = cmplx(0.0_double,0.0_double,double)
#else
      ! All of the information we need is already available in system memory.
      !   The only thing we need to do is initialize this matrix to zero.
      valeValeRhoGamma(:,:,:) = 0.0_double
#endif


      ! Accumulate the valeValeRho matrix upper triangle and electron energy.
#ifndef GAMMA
      do j = 1, numStates
         currentPopulation(:) = structuredElectronPopulation(j,i,:)

         if (sum(abs(currentPopulation(:))) < smallThresh) cycle

         electronEnergy(:) = electronEnergy(:) + currentPopulation(:) * &
               & energyEigenValues(j,i,:)

         do k = 1, spin
            call zher('U',valeDim,currentPopulation(k),valeVale(:,j,1,k),1,&
                  & valeValeRho(:,:,k),valeDim)
         enddo
      enddo

      ! In the event that the calculation includes plusUJ terms, then we need
      !   to compute the plusUJ potential from each atom with such a
      !   contribution. Of course, we also want to do it as efficiently as
      !   possible. The computation follows equation #9 from Anisimov VI,
      !   Zaanen J, Andersen OK. Band theory and Mott insulators: Hubbard U
      !   instead of Stoner I. Physical Review B, 1991;44(3):943. Available
      !   from: http://dx.doi.org/10.1103/PhysRevB.44.943.
      ! We will perform the update in two phases so that we are only doing work
      !   when the appropriate data structures are most available (to avoid
      !   having to do any extra data-reading or data-transfer just for the
      !   purpose of updating the plusUJ terms). We will do phase 1 here
      !   because at this point in the program we have access to the charge
      !   density matrix with all kpoint and electron energy level population
      !   effects accounted for already. The charge density matrix is not
      !   currently packed so it is fairly easy to reference the matrix
      !   elements. Also, at this point the charge denstiy matrix still has a
      !   up and down spin representation instead of a total and up minus down
      !   representation.

      ! Start the computation of the plusUJ terms if there are any.
      if (numPlusUJAtoms > 0) then
         call update1UJ(i,valeValeRho)
      endif

      ! Pack the matrix for easy comparison with the hamiltonian terms
      !   to be read in next.  Store the result in the appropriate packed
      !   valeVale spin array.
      do j = 1, spin
         call packMatrix(valeValeRho(:,:,j),packedValeValeSpin(:,:,j),&
               & valeDim)
      enddo
#else
      do j = 1, numStates
         currentPopulation(:) = structuredElectronPopulation(j,i,:)
         if (sum(abs(currentPopulation(:))) < smallThresh) cycle

         electronEnergy(:) = electronEnergy(:) + currentPopulation(:) * &
               & energyEigenValues(j,i,:)

         do k = 1, spin
            call dsyr('U',valeDim,currentPopulation(k),&
                  & valeValeGamma(:,j,k),1,valeValeRhoGamma(:,:,k),valeDim)
         enddo
      enddo

      ! Note: the documentation written above for the plusUJ applies here too.
      if (numPlusUJAtoms > 0) then
         call update1UJ(i,valeValeRhoGamma)
      endif

      ! Pack the matrix for easy comparison with the hamiltonian terms
      !   to be read in next.  Store the result in the appropriate packed
      !   valeVale spin array.
      do j = 1, spin
         call packMatrixGamma(valeValeRhoGamma(:,:,j),&
            & packedValeValeSpin(:,:,j),valeDim)
      enddo
#endif


      ! Allocate space to hold the overlap matrix. Later, this will also be used
      !   to read in the Hamiltonian matrix terms (KE, nuclear, and electronic
      !   potential).
      allocate (packedValeVale(dim1,valeDim*(valeDim+1)/2))

      ! Read the overlap matrix into the packedValeVale representation.
      call readPackedMatrix (atomOverlap_did(i),packedValeVale,&
            & atomDims,dim1,valeDim)

      ! In the case that the calculation is spin polarized (spin=2) then we
      !   need to convert the values in the packedValeValeSpin density matrix
      !   from being spin up and spin down to being spin up + spin down and
      !   spin up - spin down.  This is accomplished with a temporary variable.
      if (spin == 2) then
         do j = 1, valeDim*(valeDim+1)/2
            tempDensity(:) = packedValeValeSpin(:,j,1)
            packedValeValeSpin(:,j,1)=tempDensity(:)+packedValeValeSpin(:,j,2)
            packedValeValeSpin(:,j,2)=tempDensity(:)-packedValeValeSpin(:,j,2)
         enddo
      endif

      ! In the next section, we will read in the hamiltonian terms (and
      !   overlap) and perform an element by element multiplication with the
      !   density matrix for each term.

      ! Compute the integration of the charge density to show that <Psi|Psi> is
      !   actually equal to the expected number of eletrons. Note that in the
      !   spin polarized case the first index will refer to the total number of
      !   electrons and the second index will refer to the spin difference.

      do j = 1, spin
#ifndef GAMMA
         call matrixElementMult (chargeDensityTrace(j),packedValeVale,&
               & packedValeValeSpin(:,:,j),dim1,valeDim)
#else
         call matrixElementMultGamma (chargeDensityTrace(j),packedValeVale,&
               & packedValeValeSpin(:,:,j),dim1,valeDim)
#endif
      enddo


      ! Compute the nuclear contribution to the fitted potential first.
      call readPackedMatrix (atomNucOverlap_did(i),packedValeVale,&
            & atomDims,dim1,valeDim)
      do j = 1, spin
#ifndef GAMMA
         call matrixElementMult (nucPotTrace(j),packedValeVale,&
               & packedValeValeSpin(:,:,j),dim1,valeDim)
#else
         call matrixElementMultGamma (nucPotTrace(j),packedValeVale,&
               & packedValeValeSpin(:,:,j),dim1,valeDim)
#endif
      enddo

      ! Now compute the kinetic energy.
      call readPackedMatrix (atomKEOverlap_did(i),packedValeVale,&
            & atomDims,dim1,valeDim)
      do j = 1, spin
#ifndef GAMMA
         call matrixElementMult (kineticEnergyTrace(j),packedValeVale,&
               & packedValeValeSpin(:,:,j),dim1,valeDim)
#else
         call matrixElementMultGamma (kineticEnergyTrace(j),packedValeVale,&
               & packedValeValeSpin(:,:,j),dim1,valeDim)
#endif
      enddo

      ! Loop over atomic potential terms next.
      do j = 1, potDim
         call readPackedMatrix (atomPotOverlap_did(i,j),packedValeVale,&
               & atomDims,dim1,valeDim)
         do k = 1, spin
#ifndef GAMMA
            call matrixElementMult (potRho(j,k),packedValeVale,&
                  & packedValeValeSpin(:,:,k),dim1,valeDim)
#else
            call matrixElementMultGamma (potRho(j,k),packedValeVale,&
                  & packedValeValeSpin(:,:,k),dim1,valeDim)
#endif
         enddo
      enddo

      ! Deallocate the space no longer needed.
      deallocate (packedValeVale)

   enddo ! i   numKPoints

   ! Now that all k-points have been accumulated, record the total and spin
   !   difference number of electrons in the system.
   write (20,*) "Total number of electrons from <psi|psi> = ",&
         & chargeDensityTrace(1)
   if (spin == 2) then
      write (20,*) "Electron spin difference from <psi|psi> = ",&
            & chargeDensityTrace(2)
   endif

   if (spin == 1) then

      ! Obtain a summation of the electrons
      sumElecEnergy = sum(potCoeffs(:,1) * potRho(:,1)) + &
            & kineticEnergyTrace(1) + nucPotTrace(1)


      ! Write the two electron numbers (energy, and ???)
      write (20,fmt='(a17,f18.8)') 'Electron Energy = ',electronEnergy(1)
      write (20,fmt='(a17,f18.8)') 'Electron Sum    = ',sumElecEnergy
   else

      ! Obtain a summation of the electrons for each spin.

      ! Do total plus the spin difference to obtain the spin up.
      sumElecEnergy = 0.5_double * (sum(potCoeffs(:,1) * &
            & (potRho(:,1) + potRho(:,2))) + kineticEnergyTrace(1) + &
            & kineticEnergyTrace(2) + nucPotTrace(1) + nucPotTrace(2))
      ! Write the two electron numbers (energy, and ???)
      write (20,fmt='(a24,f18.8)') '(UP)   Electron Energy = ',electronEnergy(1)
      write (20,fmt='(a24,f18.8)') '(UP)   Electron Sum    = ',sumElecEnergy


      ! Do total minus the spin difference to obtain the spin down.
      sumElecEnergy = 0.5_double * (sum(potCoeffs(:,2) * &
            & (potRho(:,1) - potRho(:,2))) + kineticEnergyTrace(1) - &
            & kineticEnergyTrace(2) + nucPotTrace(1) - nucPotTrace(2))
      ! Write the spin down electronEnergy???
      write (20,fmt='(a24,f18.8)') '(DOWN) Electron Energy = ',electronEnergy(2)
      write (20,fmt='(a24,f18.8)') '(DOWN) Electron Sum    = ',sumElecEnergy

   endif

   ! Deallocate arrays associated with the electron population.
   call cleanUpPopulation

   ! Deallocate arrays associated with the secular equation and its solution.
   call cleanUpSecularEqn

   ! Deallocate other arrays and matrices defined in this subroutine.
   deallocate (structuredElectronPopulation)
   deallocate (packedValeValeSpin)
   deallocate (currentPopulation)
   deallocate (electronEnergy)
   deallocate (tempDensity)
#ifndef GAMMA
   deallocate (valeValeRho)
#else
   deallocate (valeValeRhoGamma)
#endif

   ! Log the date and time we end.
   call timeStampEnd (17)

end subroutine makeValenceRho


end module O_ValeCharge
