module O_SecularEquation

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define module data.
   real    (kind=double), allocatable, dimension (:,:,:) :: energyEigenValues

#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:,:) :: valeVale
   complex (kind=double), allocatable, dimension (:,:,:,:) :: valeValeOL
#else
   real    (kind=double), allocatable, dimension (:,:,:) :: valeValeGamma
   real    (kind=double), allocatable, dimension (:,:,:) :: valeValeOLGamma
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


subroutine secularEqnAllKP(spinDirection, numStates)

   ! Import necessary modules.
   use HDF5
   use O_Kinds
   use O_TimeStamps
   use O_KPoints, only: numKPoints
   use O_Potential, only: spin, potDim, potCoeffs, numPlusUJAtoms, currIteration
   use O_AtomicSites, only: valeDim
   use O_MainEValHDF5, only: eigenValues_did, states
   use O_MainEVecHDF5, only: eigenVectors_did, valeStates
   use O_SetupIntegralsHDF5, only: atomDims, atomOverlap_did, &
         & atomKEOverlap_did, atomNucOverlap_did, atomPotOverlap_did
#ifndef GAMMA
   use O_LAPACKZHEGV
   use O_MatrixSubs, only: readPackedMatrix, readPackedMatrixAccum, unpackMatrix
#else
   use O_LAPACKDSYGV
   use O_MatrixSubs, only: readPackedMatrix, readPackedMatrixAccum, &
         & unpackMatrixGamma
#endif

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the passed parameters.
   integer :: spinDirection
   integer :: numStates

   ! Define the local variables used in this subroutine.
   integer :: i,j ! Loop index variables
   integer :: hdferr
   integer :: dim1
   real    (kind=double), allocatable, dimension (:,:)   :: packedValeVale
   real    (kind=double), allocatable, dimension (:,:)   :: tempPackedValeVale

   ! Record the date and time that we start.
   call timeStampStart (15)

   ! Initialize the dimension of the packed matrices to include two components
   !   (real,imaginary) or just a real component.
#ifndef GAMMA
   dim1 = 2
#else
   dim1 = 1
#endif

   ! Only allocate for first spin to prevent double allocation.  These arrays
   !   will be deallocated in the makeValenceRho subroutine to accomodate the
   !   common case of one kpoint SCF calculations where it is not necessary
   !   (or efficient) to write the energy eigen values and wave function to
   !   disk with each iteration.  These matrices can simply be kept in memory
   !   and used directly when needed in makeValenceRho.  Note that only 1
   !   kpoint is needed at a time, that is the meaning of the "1" in the
   !   valeVale and valeValeOL allocation statements.
   if (spinDirection == 1) then
      allocate(energyEigenValues (numStates,numKPoints,spin))
#ifndef GAMMA
      allocate(valeVale          (valeDim,valeDim,1,spin)) ! Complex
#else
      allocate(valeValeGamma     (valeDim,valeDim,spin)) ! Real
#endif
   endif

   ! These matrices are deallocated at the end of this subroutine.
#ifndef GAMMA
   allocate(valeValeOL         (valeDim,valeDim,1,spin)) ! Complex
#else
   allocate(valeValeOLGamma    (valeDim,valeDim,spin)) ! Real
#endif
   allocate(packedValeVale     (dim1,valeDim*(valeDim+1)/2))
   allocate(tempPackedValeVale (dim1,valeDim*(valeDim+1)/2))


   ! Begin loop over all kpoints.
   do i = 1,numKPoints

      ! Prepare the matrices.
      packedValeVale(:,:) = 0.0_double
      tempPackedValeVale(:,:) = 0.0_double

#ifndef GAMMA
      valeValeOL(:,:,1,spinDirection)=cmplx(0.0_double,0.0_double,double)

      valeVale(:,:,1,spinDirection) = cmplx(0.0_double,0.0_double,double)
#else
      valeValeGamma(:,:,spinDirection) = 0.0_double
      valeValeOLGamma(:,:,spinDirection) = 0.0_double
#endif

      ! Read the nuclear potential term into packed hamiltonian.
      call readPackedMatrix(atomNucOverlap_did(i),packedValeVale,&
            & atomDims,dim1,valeDim)

      ! Read the kinetic energy term into the still packed hamiltonian.
      call readPackedMatrixAccum(atomKEOverlap_did(i),packedValeVale,&
            & tempPackedValeVale,atomDims,0.0_double,dim1,valeDim)

      ! Read the atomic potential terms into the still packed hamiltonian.
      do j = 1, potDim
         call readPackedMatrixAccum(atomPotOverlap_did(i,j),packedValeVale,&
               & tempPackedValeVale,atomDims,potCoeffs(j,spinDirection),&
               & dim1,valeDim)
      enddo

      ! Unpack the hamiltonian matrix.
#ifndef GAMMA
      call unpackMatrix(valeVale(:,:,1,spinDirection),packedValeVale,valeDim,0)
#else
      call unpackMatrixGamma(valeValeGamma(:,:,spinDirection),packedValeVale,&
            & valeDim,0)
#endif

      ! Read the atomic overlap matrix. 
      call readPackedMatrix(atomOverlap_did(i),packedValeVale,&
            & atomDims,dim1,valeDim)

      ! Unpack the overlap matrix.
#ifndef GAMMA
      call unpackMatrix(valeValeOL(:,:,1,1),packedValeVale,valeDim,0)
#else
      call unpackMatrixGamma(valeValeOLGamma(:,:,1),packedValeVale,valeDim,0)
#endif

      ! For each atom with a Hubbard U and Hund J term, we need to apply its
      !   effect on relevant matrix elements of the Hamiltonian. This is only
      !   needed in the event that we actually have atoms with these terms.
      ! Prior to application of the UJ effect on the Hamiltonian, we need to
      !   complete the update process that was started in the makeValenceRho
      !   subroutine of the valeCharge.F90 file. At that point we had the charge
      !   density matrix in a convenient structure (i.e. an unpacked matrix as
      !   opposed to a packed matrix), but we didn't have the overlap matrix in
      !   a convenient form (i.e. it would have been packed). Therefore, to
      !   multiply the charge density matrix terms by the appropriate overlap
      !   matrix terms would have required us to do a bit of annoying triangle
      !   math which at the moment I don't have time to develop. (Although it is
      !   probably fairly easy-ish.) At any rate, we have access to the overlap
      !   matrix elements in an unpacked matrix now. So, we will multiply them
      !   against the stored charge density matrix elements from the previous
      !   iteration before we apply the result to the Hamiltonian.
      ! On the first iteration of an SCF calculation there are no stored results
      !   from the previous iteration. Therefore, we will just have zeros for
      !   the plusUJ elements and there will be no need to modify the
      !   Hamiltonian.
      ! As a reminder, the reason that we need two phases for the update process
      !   is that although the charge density matrix is obtained through the
      !   Psi* * Psi operation 
      if ((numPlusUJAtoms > 0) .and. (currIteration > 1)) then
         call update2AndApplyUJ(i,spinDirection)
      endif

      ! Solve the eigen problem with a LAPACK routine.
#ifndef GAMMA
      call solveZHEGV(valeDim,numStates,valeVale(:,:,1,spinDirection),&
            & valeValeOL(:,:,1,1),energyEigenValues(:,i,spinDirection))
#else
      call solveDSYGV(valeDim,numStates,valeValeGamma(:,:,spinDirection),&
            & valeValeOLGamma(:,:,1),energyEigenValues(:,i,spinDirection))
#endif

      ! Write the energy eigenValues onto disk in HDF5 format in a.u.
      call h5dwrite_f (eigenValues_did(i,spinDirection),H5T_NATIVE_DOUBLE,&
            & energyEigenValues(:,i,spinDirection),states,hdferr)
      if (hdferr /= 0) stop 'Cannot write energy eigen values.'

      ! In the event that we have some atoms with plusUJ terms, we need to
      !   update the terms. The update depends on the charge in each of the
      !   highest d or f orbitals of the affected atoms. Therefore, we will
      !   need to read in the overlap matrix again because it was just
      !   destroyed in the ZHEGV solution. An alternative approach might be
      !   to make a copy of the overlap matrix and hold it in reserve until we
      !   need it for this calculation.
      if (numPlusUJAtoms > 0) then

         ! Read the atomic overlap matrix. 
         call readPackedMatrix(atomOverlap_did(i),packedValeVale,&
               & atomDims,dim1,valeDim)

         ! Unpack the overlap matrix.
#ifndef GAMMA
         call unpackMatrix(valeValeOL(:,:,1,1),packedValeVale,valeDim,0)
#else
         call unpackMatrixGamma(valeValeOLGamma(:,:,1),packedValeVale,valeDim,0)
#endif
      endif

#ifndef GAMMA
      ! Write the wave function only if there is more than 1 kpoint.
      !   Clearly this will exclude gamma kpoint calculations too.
      if (numKPoints >= 1) then
         ! Write the eigenVectors onto disk in HDF5 format for this
         !   kpoint and spin direction.
         call h5dwrite_f(eigenVectors_did(1,i,spinDirection),&
               & H5T_NATIVE_DOUBLE,real(valeVale(:,:numStates,1,&
               & spinDirection),double),valeStates,hdferr)
         if (hdferr /= 0) stop 'Cannot write real energy eigen vectors.'
         call h5dwrite_f(eigenVectors_did(2,i,spinDirection),&
               & H5T_NATIVE_DOUBLE,aimag(valeVale(:,:numStates,1,&
               & spinDirection)),valeStates,hdferr)
         if (hdferr /= 0) stop 'Cannot write imag energy eigen vectors.'
      endif
#endif
   enddo ! Loop i over kpoints.

   ! Deallocate unnecessary arrays and matrices.
   deallocate (packedValeVale)
   deallocate (tempPackedValeVale)
#ifndef GAMMA
   deallocate(valeValeOL)
#else
   deallocate(valeValeOLGamma)
#endif

   ! Record the date and time that we finish.
   call timeStampEnd (15)

end subroutine secularEqnAllKP


subroutine secularEqnOneKP (spinDirection,currKPoint,numStates,doSYBD)

   ! Import necessary modules.
   use HDF5
   use O_Kinds
   use O_TimeStamps
   use O_Potential, only: numPlusUJAtoms
   use O_AtomicSites, only: valeDim
   use O_PSCFBandHDF5, only: valeStatesBand, statesBand, eigenVectorsBand_did, &
         & eigenValuesBand_did
#ifndef GAMMA
   use O_LAPACKZHEGV
#else
   use O_LAPACKDSYGV
#endif

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the passed parameters.
   integer :: spinDirection
   integer :: currKPoint
   integer :: numStates
   integer :: doSYBD

   ! Define the local variables used in this subroutine.
   integer :: hdferr

   ! In the event that there are some atoms that have a plusUJ term, perform
   !   the final update to the plusUJ term and then apply the plusUJ
   !   modification to the Hamiltonian.
   if (numPlusUJAtoms > 0) then
      call update2AndApplyUJ(currKPoint,spinDirection)
   endif

   ! Solve the eigen problem with a LAPACK routine.
#ifndef GAMMA
   call solveZHEGV(valeDim,numStates,valeVale(:,:,1,spinDirection),&
         & valeValeOL(:,:,1,1),energyEigenValues(:,currKPoint,spinDirection))
#else
   call solveDSYGV(valeDim,numStates,valeValeGamma(:,:,spinDirection),&
         & valeValeOLGamma(:,:,1),energyEigenValues(:,currKPoint,spinDirection))
#endif

   if (doSYBD == 0) then
#ifndef GAMMA
      ! Write the eigenVector results to disk.
      call h5dwrite_f (eigenVectorsBand_did(1,currKPoint,spinDirection),&
            & H5T_NATIVE_DOUBLE,real(valeVale(:,1:numStates,1,spinDirection),&
            & double),valeStatesBand,hdferr)
      if (hdferr /= 0) stop 'Cannot write real energy eigen vectors.'
      call h5dwrite_f (eigenVectorsBand_did(2,currKPoint,spinDirection),&
            & H5T_NATIVE_DOUBLE,aimag(valeVale(:,1:numStates,1,spinDirection)),&
            & valeStatesBand,hdferr)
      if (hdferr /= 0) stop 'Cannot write imag energy eigen vectors.'
#else
      ! Write the eigenVector results to disk.
      call h5dwrite_f (eigenVectorsBand_did(1,currKPoint,spinDirection),&
            & H5T_NATIVE_DOUBLE,valeValeGamma(:,1:numStates,spinDirection),&
            & valeStatesBand,hdferr)
      if (hdferr /= 0) stop 'Cannot write real energy eigen vectors.'
#endif


      ! Write the eigenValue results to disk, in a.u.
      call h5dwrite_f (eigenValuesBand_did(currKPoint,spinDirection),&
            & H5T_NATIVE_DOUBLE,energyEigenValues(:,currKPoint,&
            & spinDirection),statesBand,hdferr)
      if (hdferr /= 0) stop 'Cannot write energy eigen values.'

   endif

end subroutine secularEqnOneKP

! This subroutine will update the plusUJ term values on the basis of the charge
!   density matrix (obtained from the product of wave function coefficients
!   (psi* * psi) but before the overlap matrix elements have been multiplied
!   against the charge density matrix elements. That multiplication will be done
!   later (in update2AndApplyUJ) when the overlap matrix is available in an
!   unpacked form.
! The action of this algorithm will be to multiply each appropriate 5x5 or 7x7
!   block of the charge density matrix by the U coefficient
subroutine update1UJ (currKPoint, valeValeRho)

   ! Use necessary modules.
   use O_Kinds
   use O_Potential, only: numPlusUJAtoms, plusUJAtomSize, plusUJAtomValeIndex, &
         & plusUJAtomGSElectrons, plusUJAtomValue, plusUJ, spin

   ! Make sure that no unnecessary variables are declared.
   implicit none

   ! Declare the passed parameters.
   integer, intent(in) :: currKPoint
#ifndef GAMMA
   complex (kind=double), intent(inout), dimension (:,:,:) :: valeValeRho
#else
   real (kind=double),intent(inout), dimension (:,:,:) :: valeValeRho
#endif

   ! Declare local variables.
   integer :: i,j,k ! Loop index variables.
   integer :: currUJIndex ! The valence dimension number that is one before the
         ! first d or f orbital of whichever atom is being treated at the time.
   integer :: currUJSize ! The number of orbitals (5 for d, 7 for f) for
         ! whichever atom is being treated at the time.


   ! The essential physical goal that we need to perform for each atom with a
   !   plusUJ contribution follows equation (9) from Anisimov which is: V_ms =
   !   U*SUM_m'(n_m'(-s) - n^0) + (U-J)*SUM_m'/=m(n_m's - n^0) + V^LDA. Note
   !   that m = d orbital magnetic quantum number (1-5) or f orbital magnetic
   !   quantum number (1-7); m' is a dummy index number that runs over the
   !   magnetic quantum numbers; n_m's and n_m'(-s) are the occupation numbers
   !   of the indexed magnetic quantum number and spin states; n^0 is a
   !   reference point defined to assume even distribution of all d (f)
   !   electrons across d (f) orbitals; V^LDA is the already applied LDA
   !   potential; and s = one spin direction and (-s) is the other.
   ! We will compute a spin-up V_ms and a spin-down V_m(-s) that is constructed
   !   on the basis of the actual occupation of the spin-up and spin-down
   !   orbitals.
   ! See: Anisimov VI, Zaanen J, Andersen OK. Band theory and Mott insulators:
   !   Hubbard U instead of Stoner I. Physical Review B, 1991;44(3):943.
   !   Available from: http://dx.doi.org/10.1103/PhysRevB.44.943

   ! Start collecting the update to plusUJ from each atom with a UJ term.
   do i = 1, spin
      do j = 1, numPlusUJAtoms

         currUJIndex = plusUJAtomValeIndex(j)
         currUJSize = plusUJAtomSize(j)

         do k = 1, currUJSize ! Either 5 or 7
#ifndef GAMMA
            plusUJ(k,j,i,currKPoint) = &
                  & valeValeRho(currUJIndex+k,currUJIndex+k,i)
#else
            plusUJ(k,j,i,currKPoint) = &
                  & valeValeRho(currUJIndex+k,currUJIndex+k,i)
#endif
         enddo
      enddo
   enddo

end subroutine update1UJ

! The big picture for this subroutine is that we need to modify specifc matrix
!   elements of the Hamiltonian matrix.
! Before we can do that though, we need to complete the construction of the
!   plusUJ terms. If this is the first iteration of an SCF cycle, then this
!   subroutine should not be called because the plusUJ terms are all zero. Only
!   on the second iteration will we have any knowledge about the plusUJ terms.
!   (Basically this is because I am assuming that the scfV.dat input file will
!   not have any initial values for the plusUJ terms and that the plusUJ terms
!   can only be created *after* we have the charge density matrix. The charge
!   density matrix can only be created from the single particle wave functions.)
! However, if this is part of a non-SCF calculation, then there are no
!   "iterations" and we need to complete the plusUJ right away.
subroutine update2AndApplyUJ(currKPoint,spinDirection)

   ! Use necessary modules.
   use O_Kinds
   use O_KPoints, only: numKPoints, kPointWeight
   use O_PotTypes, only: numPotTypes
   use O_AtomicSites, only: numAtomSites
   use O_Potential, only: plusUJ, plusUJAtomSize, numPlusUJAtoms, &
         & plusUJAtomValeIndex, plusUJAtomValue, plusUJAtomGSElectrons, spin

   ! Make sure that no funny variables are accidentally defined.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: currKPoint
   integer, intent(in) :: spinDirection

   ! Define local variables.
   integer :: i,j,k ! Loop index variables.
   integer :: currUJIndex ! The valence dimension number that is one before the
         ! first d or f orbital of whichever atom is being treated at the time.
   integer :: currUJSize ! The number of orbitals (5 for d, 7 for f) for
         ! whichever atom is being treated at the time.
   integer :: oppositeSpin ! If spinDirection==1 then this is 2 and vice-versa.
   real (kind=double) :: sum1, sum2 ! Two intermediate summations used in the
         ! computation of the plusUJ term for localized d- or f-orbitals in
         ! select atoms.

   ! An important thing to point out is that the plusUJ terms depend on *both*
   !   spin-up and spin-down charge densities. At this point we have access to
   !   both spin-up and spin-down charge densities in the plusUJ matrix from the
   !   update1 subroutine called in the previous SCF iteration (in the
   !   makeValenceRho subroutine). So, we *could* finalize the plusUJ term for
   !   both spin directions, but we will only do one (the current spinDirection)
   !   at a time. The reason is that at this point in the program, the valeVale
   !   matrix holds the wave function coefficients for only one spin. So, we
   !   can only apply the contributions to the up or down Hamiltonian.

   ! The essential physical goal that we need to perform for each atom with a
   !   plusUJ contribution follows equation (9) from Anisimov which is: V_ms =
   !   U*SUM_m'(n_m'(-s) - n^0) + (U-J)*SUM_m'/=m(n_m's - n^0) + V^LDA. Note
   !   that m = d orbital magnetic quantum number (1-5) or f orbital magnetic
   !   quantum number (1-7); m' is a dummy index number that runs over the
   !   magnetic quantum numbers; n_m's and n_m'(-s) are the occupation numbers
   !   of the indexed magnetic quantum number and spin states; n^0 is a
   !   reference point defined to assume even distribution of all d (f)
   !   electrons across d (f) orbitals; V^LDA is the already applied LDA
   !   potential; and s = one spin direction and (-s) is the other;
   ! We will compute a spin-up V_ms and a spin-down V_m(-s) that is constructed
   !   on the basis of the actual occupation of the spin-up and spin-down
   !   orbitals.
   ! See: Anisimov VI, Zaanen J, Andersen OK. Band theory and Mott insulators:
   !   Hubbard U instead of Stoner I. Physical Review B, 1991;44(3):943.
   !   Available from: http://dx.doi.org/10.1103/PhysRevB.44.943

   ! As another note, the division of plusUJAtomGSElectrons by (2.0*currUJSize)
   !   is designed to be a division by either 10 or 14 depending on whether
   !   this is a d or f orbital so that the value of n^0 will be (Actual number
   !   of electrons) / 10 or (Actual number of electrons) / 14. Of course, this
   !   quantity is also divided by the current kPoint weight to scale the number
   !   of electrons to the current kPoint.

   ! Establisht he value of the opposite spin direction.
   oppositeSpin = mod(spinDirection,2) + 1

   ! Update the plusUJ term and apply it to the Hamiltonian for each atom with
   !   some plusUJ contribution.
   do i = 1, numPlusUJAtoms

      ! Get the valence dimension index and size of the plusUJ term for the
      !   current atom.
      currUJIndex = plusUJAtomValeIndex(i)
      currUJSize = plusUJAtomSize(i)

      ! Update the plusUJ term by multiplying the charge density matrix by the
      !   overlap matrix elements to get the actual charge density. (Previously,
      !   the terms only held the products of wave function coefficients. This
      !   is close to the actual charge, but it lacks the overlap.)
      do j = 1, currUJSize

         ! Compute the first temporary sum for this orbital of this atom.
         sum1 = 0.0_double
         do k = 1, currUJSize
#ifndef GAMMA
            sum1 = sum1 + plusUJ(k,i,oppositeSpin,currKPoint) * &
                  & valeValeOL(currUJIndex+k,currUJIndex+k,1,1) - &
                  & plusUJAtomGSElectrons(i) / (2.0_double * &
                  & real(currUJSize,double)) * kPointWeight(currKPoint) / &
                  & real(spin,double)
#else
            sum1 = sum1 + plusUJ(k,i,oppositeSpin,currKPoint) * &
                  & valeValeOLGamma(currUJIndex+k,currUJIndex+k,1) - &
                  & plusUJAtomGSElectrons(i) / (2.0_double * &
                  & real(currUJSize,double)) * kPointWeight(currKPoint) / &
                  & real(spin,double)
#endif
         enddo
         sum1 = sum1 * plusUJAtomValue(1,i)

         ! Compute the second temporary sum for this orbital of this atom.
         sum2 = 0.0_double
         do k = 1, currUJSize

            ! Don't include m=m' terms in the summation.
            if (k == j) then
               cycle
            endif

            ! Accumulate summation terms.
#ifndef GAMMA
            sum2 = sum2 + plusUJ(k,i,spinDirection,currKPoint) * &
                  & valeValeOL(currUJIndex+k,currUJIndex+k,1,1) - &
                  & plusUJAtomGSElectrons(i) / (2.0_double * &
                  & real(currUJSize,double)) * kPointWeight(currKPoint) / &
                  & real(spin,double)
#else
            sum2 = sum2 + plusUJ(k,i,spinDirection,currKPoint) * &
                  & valeValeOLGamma(currUJIndex+k,currUJIndex+k,1) - &
                  & plusUJAtomGSElectrons(i) / (2.0_double * &
                  & real(currUJSize,double)) * kPointWeight(currKPoint) / &
                  & real(spin,double)
#endif
         enddo
         sum2 = sum2 * (plusUJAtomValue(1,i) - plusUJAtomValue(2,i))

      ! Now that the plusUJ term has been fully updated (for the current
      !   spinDirection) we will apply it to the Hamiltonian with the current
      !   spinDirection.
#ifndef GAMMA
         valeVale(currUJIndex+j,currUJIndex+j,1,spinDirection) = &
               & valeVale(currUJIndex+j,currUJIndex+j,1,spinDirection) + &
               & sum1 + sum2
#else
         valeValeGamma(currUJIndex+j,currUJIndex+j,spinDirection) = &
               & valeValeGamma(currUJIndex+j,currUJIndex+j,spinDirection) + &
               & sum1 + sum2
#endif
      enddo
   enddo

end subroutine update2AndApplyUJ


subroutine shiftEnergyEigenValues(energyShift,numStates)

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define dummy variables passed to this subroutine.
   real (kind=double) :: energyShift
   integer, intent (in) :: numStates

   ! Shift the energyEigenValues down by the requested about.
   energyEigenValues(:,:,:) = energyEigenValues(:,:,:) - energyShift

end subroutine shiftEnergyEigenValues


subroutine readEnergyEigenValuesBand(numStates)

   ! Use necessary modules
   use HDF5
   use O_Kinds
   use O_KPoints, only: numKPoints
   use O_Potential, only: spin
   use O_PSCFBandHDF5, only: statesBand, eigenValuesBand_did

   ! Make sure that no funny variables are defined.
   implicit none

   ! define variables passed to this subroutine.
   integer :: numStates

   ! Define local variables used in this subroutine.
   integer :: i,j
   integer :: hdferr
   real (kind=double), allocatable, dimension (:)   :: energyValuesTemp

   ! Allocate space for a reading buffer.
   allocate (energyValuesTemp (numStates))

   ! Read the ground state energy eigen values.

   ! Loop over each kpoint and spin direction to read the energy values.
   do i = 1, numKPoints
      do j = 1, spin
         call h5dread_f (eigenValuesBand_did(i,j),H5T_NATIVE_DOUBLE,&
               & energyValuesTemp(:),statesBand,hdferr)
         if (hdferr /= 0) stop 'Failed to read energy eigen values'

         ! Copy the necessary values into the final array.
         energyEigenValues(:,i,j) = energyValuesTemp(:)
      enddo
   enddo

   ! Deallocate the reading buffer
   deallocate (energyValuesTemp)

end subroutine readEnergyEigenValuesBand


subroutine appendExcitedEValsBand (firstStateIndex,numStates)

   ! Use necessary modules
   use HDF5
   use O_Kinds
   use O_KPoints, only: numKPoints
   use O_Potential, only: spin
   use O_PSCFBandHDF5, only: statesBand, eigenValuesBand2_did

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define variables passed to this subroutine.
   integer :: firstStateIndex
   integer :: numStates

   ! Define local variables used in this subroutine.
   integer :: i,j
   integer :: hdferr
   real (kind=double), allocatable, dimension (:)   :: energyValuesTemp

   ! Allocate space for a reading buffer.
   allocate (energyValuesTemp (numStates))

   ! It is assumed that we know the highest occupied state for each kpoint.
   !   We now read the excited state energy eigen values and merge them into
   !   the ground state ones with the merge points being the highest occupied
   !   states.  We must be careful when dealing with degenerate states.
   do i = 1, numKPoints
      do j = 1, spin

         ! Get the excited state energy values.
         call h5dread_f (eigenValuesBand2_did(i,j),H5T_NATIVE_DOUBLE,&
               & energyValuesTemp(:),statesBand,hdferr)
         if (hdferr /= 0) stop 'Failed to read excited energy eigen values'

         energyEigenValues(firstStateIndex:numStates,i,j) = &
               & energyValuesTemp(firstStateIndex:numStates)

      enddo
   enddo

   ! Deallocate the reading buffer
   deallocate (energyValuesTemp)

end subroutine appendExcitedEValsBand


subroutine preserveValeValeOL

   ! Use necessary modules.
   use O_Potential, only: spin, numPlusUJAtoms
   use O_AtomicSites, only: valeDim

   ! Make sure that no funny variables are defined.
   implicit none

   ! In some cases (e.g. during a spin polarized calculation) it is necessary
   !   to preserve the valeValeOL during the diagonalization.  This is
   !   because the LAPACK routine will destroy the overlap matrix because it
   !   is the same for spin up or spin down (only one copy is made).  We will
   !   copy the valeValeOL (spin index 1) into the spin index 2 part of the
   !   same matrix (valeValeOL).
   if ((spin == 2) .or. (numPlusUJAtoms > 0)) then
#ifndef GAMMA
      valeValeOL(:valeDim,:valeDim,1,2) = valeValeOL(:valeDim,:valeDim,1,1)
#else
      valeValeOLGamma(:valeDim,:valeDim,2) = &
            & valeValeOLGamma(:valeDim,:valeDim,1)
#endif
   endif

   ! IMPORTANT NOTE:  This whole thing could be "improved" in terms of memory
   !   usage if necessary.  Instead of copying to the spin index 2 of the
   !   valeValeOL we could copy to spin index 2 of the valeVale because it
   !   has not yet been used.  Along with this, the valeValeOL would only need
   !   to be allocated 1 spin index worth of memory.  This could be important,
   !   but it will not be pursued now because it is an ugly and potentially
   !   more confusing thing to do.  (The next subroutine would also need to be
   !   changed accordingly.)

end subroutine preserveValeValeOL


subroutine restoreValeValeOL

   ! Use necessary modules.
   use O_Potential, only: spin
   use O_AtomicSites, only: valeDim

   ! Make sure that no funny variables are defined.
   implicit none

   ! As with the above subroutine, the valeValeOL matrix must sometimes be
   !   saved from destruction during the diagonalization.  There is no spin
   !   aspect to the overlap matrix and so the matrix was copied from the
   !   spin index 1 to the spin index 2 for safe keeping.  This subroutine
   !   will restore the matrix from spin index 2 to spin index 1.
   if (spin == 2) then
#ifndef GAMMA
      valeValeOL(:valeDim,:valeDim,1,1) = valeValeOL(:valeDim,:valeDim,1,2)
#else
      valeValeOLGamma(:valeDim,:valeDim,1) = &
            & valeValeOLGamma(:valeDim,:valeDim,2)
#endif
   endif
end subroutine restoreValeValeOL



subroutine readDataSCF(h,i,numStates)

   ! Import necessary data modules.
   use O_KPoints, only: numKPoints
   use O_AtomicSites, only: valeDim
   use O_SetupIntegralsHDF5, only: atomDims, atomOverlap_did
   use O_MainEVecHDF5, only: valeStates, eigenVectors_did
#ifndef GAMMA
   use O_MatrixSubs, only: readMatrix, readPackedMatrix, unpackMatrix
#else
   use O_MatrixSubs, only: readPackedMatrix, unpackMatrixGamma
#endif

   ! Define passed parameters.
   integer, intent(in) :: h ! Spin variable.
   integer, intent(in) :: i ! KPoint variable
   integer, intent(in) :: numStates

   ! Define local variables.
   integer :: dim1
   real (kind=double), allocatable, dimension (:,:) :: packedValeVale
#ifndef GAMMA
   real (kind=double), allocatable, dimension (:,:) :: tempRealValeVale
   real (kind=double), allocatable, dimension (:,:) :: tempImagValeVale
#endif

   ! Initialize the dimension of the packed matrices to include two components
   !   (real,imaginary) or just a real component.
#ifndef GAMMA
   dim1 = 2
#else
   dim1 = 1
#endif

   ! Allocate space to read the packed overlap matrix.
   allocate (packedValeVale(dim1,valeDim*(valeDim+1)/2))

   ! Read the overlap matrix.  The tempPackedMatrix is not used.
   call readPackedMatrix(atomOverlap_did(i),packedValeVale,&
         & atomDims,dim1,valeDim)

   ! Unpack the matrix.
#ifndef GAMMA
   call unpackMatrix(valeValeOL(:,:,1,1),packedValeVale,valeDim,1)
#else
   call unpackMatrixGamma(valeValeOLGamma(:,:,1),packedValeVale,valeDim,1)
#endif
   ! Deallocate the packed matrix used to read and unpack the data.
   deallocate (packedValeVale)

#ifndef GAMMA
   ! Read the wave functions for this kpoint from the datasets into
   !   the valeVale matrix.  If numKPoints==1, the wave functions are
   !   already in the valeVale(1:valeDim,1:numStates) matrix.
   if (numKPoints > 1) then

      ! Allocate space to read the wave functions.
      allocate (tempRealValeVale (valeDim,numStates))
      allocate (tempImagValeVale (valeDim,numStates))

      call readMatrix(eigenVectors_did(:,i,h),valeVale(:,:numStates,1,1),&
            & tempRealValeVale(:,:),tempImagValeVale(:,:),&
            & valeStates,valeDim,numStates)

      ! Deallocate the space to read the wave functions.
      deallocate (tempRealValeVale)
      deallocate (tempImagValeVale)
   endif
#endif

end subroutine readDataSCF



subroutine readDataPSCF(h,i,numStates)

   ! Use necessary modules.
   use O_AtomicSites, only: valeDim
   use O_PSCFBandHDF5, only: valeStatesBand, valeValeBand, &
         & eigenVectorsBand_did, valeValeBand_did
#ifndef GAMMA
   use O_MatrixSubs, only: readMatrix, readPackedMatrix, unpackMatrix
#else
   use O_MatrixSubs, only: readMatrixGamma, readPackedMatrix, &
         & unpackMatrixGamma
#endif

   ! Define passed parameters.
   integer :: h ! Spin variable.
   integer :: i ! KPoint variable
   integer :: numStates

   ! Define local variables.
   integer :: dim1
   real (kind=double), allocatable, dimension (:,:) :: packedValeVale
#ifndef GAMMA
   real (kind=double), allocatable, dimension (:,:) :: tempRealValeVale
   real (kind=double), allocatable, dimension (:,:) :: tempImagValeVale
#endif

   ! Initialize the dimension of the packed matrices to include two components
   !   (real,imaginary) or just a real component.
#ifndef GAMMA
   dim1 = 2
#else
   dim1 = 1
#endif

#ifndef GAMMA
   ! Allocate space to read the complex wave function.
   allocate (tempRealValeVale (valeDim,numStates))
   allocate (tempImagValeVale (valeDim,numStates))

   ! Read the complex wave function from the datasets.
   call readMatrix(eigenVectorsBand_did(:,i,h),valeVale(:,:,1,1),&
         & tempRealValeVale(:,:),tempImagValeVale(:,:),&
         & valeStatesBand,valeDim,numStates)

   ! Deallocate the space to read the complex wave function.
   deallocate (tempRealValeVale)
   deallocate (tempImagValeVale)
#else
   ! Read the real wave function from the datasets.
   call readMatrixGamma(eigenVectorsBand_did(1,i,h),&
         & valeValeGamma(:,:,1),valeStatesBand,valeDim,numStates)
#endif

   ! Allocate space to read the packed overlap matrix.
   allocate (packedValeVale(dim1,valeDim*(valeDim+1)/2))

   ! Read the overlap matrix.
   call readPackedMatrix(valeValeBand_did(i),packedValeVale,&
         & valeValeBand,dim1,valeDim)
#ifndef GAMMA
   call unpackMatrix(valeValeOL(:,:,1,1),packedValeVale,valeDim,1)
#else
   call unpackMatrixGamma(valeValeOLGamma(:,:,1),packedValeVale,valeDim,1)
#endif
   ! Deallocate the space used to read the matrix.
   deallocate (packedValeVale)


end subroutine readDataPSCF



subroutine cleanUpSecularEqn

   ! Make sure that no funny variables are defined.
   implicit none

   deallocate (energyEigenValues)
#ifndef GAMMA
   deallocate (valeVale)
#else
   deallocate (valeValeGamma)
#endif

end subroutine cleanUpSecularEqn

end module O_SecularEquation
