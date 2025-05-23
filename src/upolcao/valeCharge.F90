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
   real (kind=double), allocatable, dimension (:)   :: massVelocityTrace
   real (kind=double), allocatable, dimension (:,:) :: dipoleMomentTrace

   real (kind=double), allocatable, dimension (:,:)   :: packedValeVale
   real (kind=double), allocatable, dimension (:,:,:) :: packedValeValeRho

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine makeValenceRho(inSCF)

   ! Import the necessary modules.
   use O_Kinds
   use O_MPI
   use O_TimeStamps
   use O_CommandLine, only: doDIMO_SCF, doDIMO_PSCF, doForce_SCF, doForce_PSCF
   use O_AtomicSites, only: valeDim
   use O_Input, only: numStates
   use O_KPoints, only: numKPoints
   use O_Constants, only: smallThresh
   use O_Potential, only: rel,spin,potDim,potCoeffs,&
         & numPlusUJAtoms, converged
   use O_Populate, only: electronPopulation,cleanUpPopulation
   use O_SCFIntegralsHDF5, only: atomOverlap_did,atomKEOverlap_did, &
         & atomMVOverlap_did,atomNPOverlap_did,atomDMOverlap_did, &
         & atomPotOverlap_did,packedVVDims
   use O_PSCFIntegralsHDF5, only: atomDMOverlapPSCF_did
#ifndef GAMMA
   use O_BLASZHER
   use O_SecularEquation, only: valeVale,cleanUpSecularEqn,energyEigenValues,&
         & update1UJ, readDataSCF, readDataPSCF
   use O_MatrixSubs, only: readPackedMatrix,matrixElementMult,packMatrix
   use O_Force, only: computeForce
#else
   use O_BLASDSYR
   use O_SecularEquation, only: valeValeGamma, cleanUpSecularEqn, &
         & energyEigenValues, update1UJ, readDataSCF, readDataPSCF
   use O_MatrixSubs, only: readPackedMatrix,matrixElementMultGamma,&
         & packMatrixGamma
   use O_Force, only: computeForceGamma
#endif

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: inSCF

   ! Define the local variables used in this subroutine.
   integer :: i,j,k ! Loop index variables
   integer :: skipKP
   integer :: dim1
   integer :: energyLevelCounter
   real (kind=double) :: sumElecEnergy
   real (kind=double), allocatable, dimension (:)     :: tempDensity
   real (kind=double), allocatable, dimension (:)     :: electronEnergy
   real (kind=double), allocatable, dimension (:)     :: currentPopulation
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
   allocate (packedValeValeRho(dim1,valeDim*(valeDim+1)/2,spin))

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
   !   kinetic energy, mass velocity, and dipole moment.
   if (inSCF == 1) then
      allocate (chargeDensityTrace(spin))
      allocate (nucPotTrace(spin))
      allocate (kineticEnergyTrace(spin))
      if (rel == 1) then
         allocate (massVelocityTrace(spin))
      endif
   endif
   if (((doDIMO_SCF == 1) .and. (converged == 1)) .or. &
         & ((doDIMO_PSCF == 1) .and. (inSCF == 0))) then
      allocate (dipoleMomentTrace(3, spin))
   endif

   ! Allocate space to hold the currentPopulation based on spin
   allocate (currentPopulation (spin))
   allocate (electronEnergy (spin))
   allocate (structuredElectronPopulation (numStates,numKPoints,spin))

   ! Initialize variables from data modules.
   potRho(:,:) = 0.0_double
   if (inSCF == 1) then
      nucPotTrace(:) = 0.0_double
      kineticEnergyTrace(:) = 0.0_double
      chargeDensityTrace(:) = 0.0_double
      if (rel == 1) then
         massVelocityTrace(:) = 0.0_double
      endif
   endif
   if (((doDIMO_SCF == 1) .and. (converged == 1)) .or. &
         & ((doDIMO_PSCF == 1) .and. (inSCF == 0))) then
      dipoleMomentTrace(:,:) = 0.0_double
   endif

   ! Initialize local variables
   electronEnergy(:) = 0.0_double


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
      !   a double complex matrix because the wave function data contains both
      !   real and imaginary parts and has no symmetry that might permit
      !   packing or use of triangular form (if it were a Hermitian matrix
      !   instead).
      ! Note that it is only necessary to go through the initialization action
      !   if there are more than 1 kpoint.  For the 1 kpoint case, the valeVale
      !   matrix was not changed so it can still be used here.  Also note
      !   that the gammaKPoint option will never enter the "if" block.
#ifndef GAMMA
!      if (numKPoints > 1) then

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

      ! Determine if we are doing the valeCharge in a post-SCF calculation
      !   or within an SCF calculation.
      do j = 1, spin
         if (inSCF == 1) then
            call readDataSCF(j,i,numStates,0) ! Read wave functions only.
         else
            call readDataPSCF(j,i,numStates,0) ! Read wave functions only.
         endif
      enddo
!      endif

      ! Initialize matrix to receive the valeVale density matrix (square of the
      !   wave function).
      valeValeRho(:,:,:) = cmplx(0.0_double,0.0_double,double)
#else
      skipKP = 0 ! Avoid compiler warning about unused variable.
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
            call zher('U',valeDim,currentPopulation(k),valeVale(:,j,k),1,&
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
         call packMatrix(valeValeRho(:,:,j),packedValeValeRho(:,:,j),&
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
               & packedValeValeRho(:,:,j),valeDim)
      enddo
#endif


      ! Allocate space to hold the overlap matrix. Later, this will also be
      !   used to read in the Hamiltonian matrix terms (KE, nuclear, electronic
      !   potential, and (if needed) the mass velocity).
      allocate (packedValeVale(dim1,valeDim*(valeDim+1)/2))

      ! Read the overlap matrix into the packedValeVale representation.
      call readPackedMatrix (atomOverlap_did(i),packedValeVale,&
            & packedVVDims,dim1,valeDim,1)

      ! In the case that the calculation is spin polarized (spin=2) then we
      !   need to convert the values in the packedValeValeRho density matrix
      !   from being spin up and spin down to being spin up + spin down and
      !   spin up - spin down.  This is accomplished with a temporary variable.
      if (spin == 2) then
         do j = 1, valeDim*(valeDim+1)/2
            tempDensity(:) = packedValeValeRho(:,j,1)
            packedValeValeRho(:,j,1)=tempDensity(:)+packedValeValeRho(:,j,2)
            packedValeValeRho(:,j,2)=tempDensity(:)-packedValeValeRho(:,j,2)
         enddo
      endif

      ! In the next section, we will read in the hamiltonian terms (and
      !   overlap) and perform an element by element multiplication with the
      !   density matrix for each term.

      ! Compute the integration of the charge density to show that <Psi|Psi> is
      !   actually equal to the expected number of eletrons. Note that in the
      !   spin polarized case the first index will refer to the total number of
      !   electrons and the second index will refer to the spin difference.

      if (inSCF == 1) then

         do j = 1, spin ! j=1 -> Total; j=2 -> Difference
#ifndef GAMMA
            call matrixElementMult (chargeDensityTrace(j),packedValeVale,&
                  & packedValeValeRho(:,:,j),dim1,valeDim)
#else
            call matrixElementMultGamma (chargeDensityTrace(j),&
                  & packedValeVale,packedValeValeRho(:,:,j),dim1,valeDim)
#endif
         enddo

         ! Compute the nuclear contribution to the fitted potential first.
         call readPackedMatrix (atomNPOverlap_did(i),packedValeVale,&
               & packedVVDims,dim1,valeDim,1)
         do j = 1, spin ! j=1 -> Total; j=2 -> Difference
#ifndef GAMMA
            call matrixElementMult (nucPotTrace(j),packedValeVale,&
                  & packedValeValeRho(:,:,j),dim1,valeDim)
#else
            call matrixElementMultGamma (nucPotTrace(j),packedValeVale,&
                  & packedValeValeRho(:,:,j),dim1,valeDim)
#endif
         enddo

         ! Now compute the kinetic energy.
         call readPackedMatrix (atomKEOverlap_did(i),packedValeVale,&
               & packedVVDims,dim1,valeDim,1)
         do j = 1, spin ! j=1 -> Total; j=2 -> Difference
#ifndef GAMMA
            call matrixElementMult (kineticEnergyTrace(j),packedValeVale,&
                  & packedValeValeRho(:,:,j),dim1,valeDim)
#else
            call matrixElementMultGamma (kineticEnergyTrace(j),packedValeVale,&
                  & packedValeValeRho(:,:,j),dim1,valeDim)
#endif
         enddo

         ! If needed, compute the mass velocity.
         if (rel == 1) then
            call readPackedMatrix (atomMVOverlap_did(i),packedValeVale,&
                  & packedVVDims,dim1,valeDim,1)
            do j = 1, spin ! j=1 -> Total; j=2 -> Difference
#ifndef GAMMA
               call matrixElementMult (massVelocityTrace(j),packedValeVale,&
                     & packedValeValeRho(:,:,j),dim1,valeDim)
#else
               call matrixElementMultGamma (massVelocityTrace(j),packedValeVale,&
                     & packedValeValeRho(:,:,j),dim1,valeDim)
#endif
            enddo
         endif

         ! Loop over atomic potential terms next.
         do j = 1, potDim
            call readPackedMatrix (atomPotOverlap_did(i,j),packedValeVale,&
                  & packedVVDims,dim1,valeDim,1)
            do k = 1, spin ! j=1 -> Total; j=2 -> Difference
#ifndef GAMMA
               call matrixElementMult (potRho(j,k),packedValeVale,&
                     & packedValeValeRho(:,:,k),dim1,valeDim)
#else
               call matrixElementMultGamma (potRho(j,k),packedValeVale,&
                     & packedValeValeRho(:,:,k),dim1,valeDim)
#endif
            enddo
         enddo
      endif ! inSCF == 1

      ! If needed, compute the dipole moment.
      if (((doDIMO_SCF == 1) .and. (converged == 1)) .or. &
            & ((doDIMO_PSCF == 1) .and. (inSCF == 0))) then
         do j = 1, 3 ! xyz directions

            if (inSCF == 0) then
               call readPackedMatrix (atomDMOverlap_did(i,j),packedValeVale,&
                     & packedVVDims,dim1,valeDim,1)
            else
               call readPackedMatrix (atomDMOverlapPSCF_did(i,j),&
                     & packedValeVale,packedVVDims,dim1,valeDim,1)
            endif
            do k = 1, spin
#ifndef GAMMA
               call matrixElementMult (dipoleMomentTrace(j,k),packedValeVale,&
                     & packedValeValeRho(:,:,k),dim1,valeDim)
#else
               call matrixElementMultGamma (dipoleMomentTrace(j,k),&
                     & packedValeVale,packedValeValeRho(:,:,k),dim1,valeDim)
#endif
            enddo
         enddo
      endif

      ! If needed, compute the forces
      if (((doForce_SCF == 1) .and. (converged == 1)) .or. &
            & ((doForce_PSCF == 1) .and. (inSCF == 0))) then
         do j = 1, 3 ! xyz directions
            do k = 1, spin
#ifndef GAMMA
               call computeForce(valeValeRho,i,k,j)
#else
               call computeForceGamma(valeValeRhoGamma,k,j)
#endif
            enddo
         enddo
      endif

      ! Deallocate the space no longer needed.
      deallocate (packedValeVale)

   enddo ! i   numKPoints

   ! Now that all k-points have been accumulated, record the total and spin
   !   difference number of electrons in the system.
   if (inSCF == 1) then
      if (mpiRank == 0) then
         write (20,*) "Total number of electrons from <psi|psi> = ",&
               & chargeDensityTrace(1)
         if (spin == 2) then
            write (20,*) "Electron spin difference from <psi|psi> = ",&
                  & chargeDensityTrace(2)
         endif
      endif

      if (spin == 1) then

         ! Obtain a summation of the electrons.
         if (rel == 0) then
            sumElecEnergy = sum(potCoeffs(:,1) * potRho(:,1)) + &
                  & kineticEnergyTrace(1) + nucPotTrace(1)
         else
            sumElecEnergy = sum(potCoeffs(:,1) * potRho(:,1)) + &
                  & kineticEnergyTrace(1) + massVelocityTrace(1) + &
                  & nucPotTrace(1)
         endif


         ! Write the two electron numbers (both are energy values, but they
         !   are computed through different means.)
         if (mpiRank == 0) then
            write (20,fmt='(a17,f18.8)') 'Electron Energy = ',electronEnergy(1)
            write (20,fmt='(a17,f18.8)') 'Electron Sum    = ',sumElecEnergy
         endif
      else

         ! Obtain a summation of the electrons for each spin.

         ! Do total plus the spin difference to obtain the spin up.
         if (rel == 0) then
            sumElecEnergy = 0.5_double * (sum(potCoeffs(:,1) * &
                  & (potRho(:,1) + potRho(:,2))) + kineticEnergyTrace(1) + &
                  & kineticEnergyTrace(2) + nucPotTrace(1) + nucPotTrace(2))
         else
            sumElecEnergy = 0.5_double * (sum(potCoeffs(:,1) * &
                  & (potRho(:,1) + potRho(:,2))) + kineticEnergyTrace(1) + &
                  & kineticEnergyTrace(2) + massVelocityTrace(1) + &
                  & massVelocityTrace(2) + nucPotTrace(1) + nucPotTrace(2))
         endif
         ! Write the two electron numbers. (Both are energy values, but they
         !   are computed thorugh different means.
         if (mpiRank == 0) then
            write (20,fmt='(a24,f18.8)') '(UP)   Electron Energy = ', &
                  & electronEnergy(1)
            write (20,fmt='(a24,f18.8)') '(UP)   Electron Sum    = ', &
                  & sumElecEnergy
         endif


         ! Do total minus the spin difference to obtain the spin down.
         if (rel == 0) then
            sumElecEnergy = 0.5_double * (sum(potCoeffs(:,2) * &
                  & (potRho(:,1) - potRho(:,2))) + kineticEnergyTrace(1) - &
                  & kineticEnergyTrace(2) + nucPotTrace(1) - nucPotTrace(2))
         else
            sumElecEnergy = 0.5_double * (sum(potCoeffs(:,2) * &
                  & (potRho(:,1) - potRho(:,2))) + kineticEnergyTrace(1) - &
                  & kineticEnergyTrace(2) + massVelocityTrace(1) - &
                  & massVelocityTrace(2) + nucPotTrace(1) - nucPotTrace(2))
         endif

         ! Write the spin down electronEnergy???
         if (mpiRank == 0) then
            write (20,fmt='(a24,f18.8)') '(DOWN) Electron Energy = ', &
                  & electronEnergy(2)
            write (20,fmt='(a24,f18.8)') '(DOWN) Electron Sum    = ', &
                  & sumElecEnergy
         endif
      endif
   endif


   ! If the dipole moment calculation was requested, then print it. Note that
   !   the units of the dipole moment are e * a_0 where e is the electronic
   !   charge and a_0 is the bohr radius.
   ! The units of the dipole moment calculation are ea0 where e is the number
   !    of valence electrons and a0 is the bohr radius. 
   if (mpiRank == 0) then
      if (((doDIMO_SCF == 1) .and. (converged == 1)) .or. &
            & ((doDIMO_PSCF == 1) .and. (inSCF == 0))) then
         write (74,*) "Dipole Moment x-direction: ", dipoleMomentTrace(1,1)
         write (74,*) "Dipole Moment y-direction: ", dipoleMomentTrace(2,1)
         write (74,*) "Dipole Moment z-direction: ", dipoleMomentTrace(3,1)
         if (spin == 2) then
            write (75,*) "Dipole Moment x-direction: ", dipoleMomentTrace(1,2) 
            write (75,*) "Dipole Moment y-direction: ", dipoleMomentTrace(2,2)
            write (75,*) "Dipole Moment z-direction: ", dipoleMomentTrace(3,2)
         endif
      endif
   endif

   ! Deallocate arrays associated with the electron population.
   call cleanUpPopulation

   ! Deallocate arrays associated with the secular equation and its solution.
   call cleanUpSecularEqn

   ! Deallocate other arrays and matrices defined in this subroutine.
   deallocate (structuredElectronPopulation)
   deallocate (packedValeValeRho)
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
