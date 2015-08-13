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
   use O_Potential, only: spin, potDim, potCoeffs
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
      valeVale(:,:,1,spinDirection)=cmplx(0.0_double,0.0_double,double)
      valeValeOL(:,:,1,spinDirection)=cmplx(0.0_double,0.0_double,double)
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
            & energyEigenValues(:numStates,i,spinDirection),states,hdferr)
      if (hdferr /= 0) stop 'Can not write energy eigen values.'

#ifndef GAMMA
      ! Write the wave function only if there is more than 1 kpoint.
      !   Clearly this will exclude gamma kpoint calculations too.
      if (numKPoints >= 1) then
         ! Write the eigenVectors onto disk in HDF5 format for this
         !   kpoint and spin direction.
         call h5dwrite_f(eigenVectors_did(1,i,spinDirection),&
               & H5T_NATIVE_DOUBLE,real(valeVale(:,:numStates,1,&
               & spinDirection),double),valeStates,hdferr)
         if (hdferr /= 0) stop 'Can not write real energy eigen vectors.'
         call h5dwrite_f(eigenVectors_did(2,i,spinDirection),&
               & H5T_NATIVE_DOUBLE,aimag(valeVale(:,:numStates,1,&
               & spinDirection)),valeStates,hdferr)
         if (hdferr /= 0) stop 'Can not write imag energy eigen vectors.'
      endif
#endif
   enddo

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


subroutine secularEqnOneKP (spinDirection,choiceKP,numStates,doSYBD)

   ! Import necessary modules.
   use HDF5
   use O_Kinds
   use O_TimeStamps
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
   integer :: choiceKP
   integer :: numStates
   integer :: doSYBD

   ! Define the local variables used in this subroutine.
   integer :: hdferr

   ! Solve the eigen problem with a LAPACK routine.
#ifndef GAMMA
   call solveZHEGV(valeDim,numStates,valeVale(:,:,1,spinDirection),&
         & valeValeOL(:,:,1,1),energyEigenValues(:,choiceKP,spinDirection))
#else
   call solveDSYGV(valeDim,numStates,valeValeGamma(:,:,spinDirection),&
         & valeValeOLGamma(:,:,1),energyEigenValues(:,choiceKP,spinDirection))
#endif

   if (doSYBD == 0) then
#ifndef GAMMA
      ! Write the eigenVector results to disk.
      call h5dwrite_f (eigenVectorsBand_did(1,choiceKP,spinDirection),&
            & H5T_NATIVE_DOUBLE,real(valeVale(:,1:numStates,1,spinDirection),&
            & double),valeStatesBand,hdferr)
      if (hdferr /= 0) stop 'Can not write real energy eigen vectors.'
      call h5dwrite_f (eigenVectorsBand_did(2,choiceKP,spinDirection),&
            & H5T_NATIVE_DOUBLE,aimag(valeVale(:,1:numStates,1,spinDirection)),&
            & valeStatesBand,hdferr)
      if (hdferr /= 0) stop 'Can not write imag energy eigen vectors.'
#else
      ! Write the eigenVector results to disk.
      call h5dwrite_f (eigenVectorsBand_did(1,choiceKP,spinDirection),&
            & H5T_NATIVE_DOUBLE,valeValeGamma(:,1:numStates,spinDirection),&
            & valeStatesBand,hdferr)
      if (hdferr /= 0) stop 'Can not write real energy eigen vectors.'
#endif


      ! Write the eigenValue results to disk, in a.u.
      call h5dwrite_f (eigenValuesBand_did(choiceKP,spinDirection),&
            & H5T_NATIVE_DOUBLE,energyEigenValues(:numStates,choiceKP,&
            & spinDirection),statesBand,hdferr)
      if (hdferr /= 0) stop 'Can not write energy eigen values.'

   endif

end subroutine secularEqnOneKP


subroutine shiftEnergyEigenValues(energyShift,numStates)

   ! Make sure that no funny variables are defined.
   implicit none

   ! define vriables passed to this subroutine
   integer, intent (in) :: numStates

   ! Define passed dummy parameters.
   real (kind=double) :: energyShift

   ! Shift the energyEigenValues down by the requested about.
   energyEigenValues(:numStates,:,:) = &
         & energyEigenValues(:numStates,:,:) - energyShift

end subroutine shiftEnergyEigenValues


!subroutine convertEnergyEigenValuesToeV
!
!   ! Use necessary modules.
!   use O_Constants
!   use O_Input
!
!   ! Make sure that no funny variables are defined.
!   implicit none
!
!   energyEigenValues(:numStates,:,:) = energyEigenValues(:numStates,:,:) * &
!         & hartree
!
!end subroutine convertEnergyEigenValuesToeV


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
               & energyValuesTemp(:numStates),statesBand,hdferr)
         if (hdferr /= 0) stop 'Failed to read energy eigen values'

         ! Copy the necessary values into the final array.
         energyEigenValues(:numStates,i,j) = energyValuesTemp(:numStates)
      enddo
   enddo

   ! Deallocate the reading buffer
   deallocate (energyValuesTemp)

end subroutine readEnergyEigenValuesBand


subroutine appendExcitedEnergyEigenValuesBand (firstStateIndex,numStates)

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
   integer :: i,j,k
   integer :: hdferr
   integer :: firstUnoccupiedState
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
               & energyValuesTemp(:numStates),statesBand,hdferr)
         if (hdferr /= 0) stop 'Failed to read excited energy eigen values'

         energyEigenValues(firstStateIndex:numStates,i,j) = &
               & energyValuesTemp(firstStateIndex:numStates)

      enddo
   enddo

   ! Deallocate the reading buffer
   deallocate (energyValuesTemp)

end subroutine appendExcitedEnergyEigenValuesBand


subroutine preserveValeValeOL

   ! Use necessary modules.
   use O_Potential, only: spin
   use O_AtomicSites, only: valeDim

   ! Make sure that no funny variables are defined.
   implicit none

   ! In some cases during a spin polarized calculations it is necessary
   !   to preserve the valeValeOL during the diagonalization.  This is
   !   because the LAPACK routine will destroy the overlap matrix because it
   !   is the same for spin up or spin down (only one copy is made).  We will
   !   copy the valeValeOL (spin index 1) into the spin index 2 part of the
   !   same matrix (valeValeOL).
   if (spin == 2) then
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
   integer :: h ! Spin variable.
   integer :: i ! KPoint variable
   integer :: numStates

   ! Define local variables.
   integer :: dim1
   real (kind=double), allocatable, dimension (:,:) :: packedValeVale
   real (kind=double), allocatable, dimension (:,:) :: tempRealValeVale
   real (kind=double), allocatable, dimension (:,:) :: tempImagValeVale

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
   real (kind=double), allocatable, dimension (:,:) :: tempRealValeVale
   real (kind=double), allocatable, dimension (:,:) :: tempImagValeVale

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
