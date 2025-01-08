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
   real (kind=double), allocatable, dimension (:,:,:) :: energyEigenValues
         ! States, KPoints, Spin

#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: valeVale
   complex (kind=double), allocatable, dimension (:,:)   :: valeValeOL
   complex (kind=double), allocatable, dimension (:,:,:) :: valeValeMM
#else
   real (kind=double), allocatable, dimension (:,:,:) :: valeValeGamma
   real (kind=double), allocatable, dimension (:,:)   :: valeValeOLGamma
   real (kind=double), allocatable, dimension (:,:,:) :: valeValeMMGamma
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


subroutine secularEqnSCF(spinDirection, numStates)

   ! Import necessary modules.
   use HDF5
   use MPI_F08
   use O_MPI
   use O_Kinds
   use O_TimeStamps
   use O_KPoints, only: numKPoints
   use O_Potential, only: rel, spin, potDim, potCoeffs, numPlusUJAtoms, &
         & currIteration
   use O_AtomicSites, only: valeDim
   use O_SCFEigValHDF5, only: eigenValues_did, states
   use O_SCFEigVecHDF5, only: eigenVectors_did, eigenVectors_aid, valeStates
   use O_SCFIntegralsHDF5, only: packedVVDims, atomOverlap_did, &
         & atomKEOverlap_did, atomMVOverlap_did, atomNPOverlap_did, &
         & atomPotOverlap_did
#ifndef GAMMA
   use O_LAPACKZHEGV
   use O_MatrixSubs, only: readPackedMatrix, readPackedMatrixAccum,&
         & unpackMatrix
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
!integer :: k, l, m, n, o, p
   integer :: hdferr
   integer :: hdf5Status
   integer(hsize_t), dimension (1) :: attribIntDims ! Attribute dataspace dim
   integer :: dim1
   real    (kind=double), allocatable, dimension (:,:)   :: packedValeVale
   real    (kind=double), allocatable, dimension (:,:)   :: tempPackedValeVale
!complex(kind=double), allocatable, dimension(:,:) :: identity

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
      allocate(valeVale          (valeDim,valeDim,spin)) ! Complex
#else
      allocate(valeValeGamma     (valeDim,valeDim,spin)) ! Real
#endif
   endif

   ! These matrices are deallocated at the end of this subroutine.
#ifndef GAMMA
   allocate(valeValeOL         (valeDim,valeDim)) ! Complex
#else
   allocate(valeValeOLGamma    (valeDim,valeDim)) ! Real
#endif
   allocate(packedValeVale     (dim1,valeDim*(valeDim+1)/2))
   allocate(tempPackedValeVale (dim1,valeDim*(valeDim+1)/2))


   ! Begin loop over all kpoints.
   do i = 1,numKPoints

      ! If the eigen vectors for this kpoint and spin direction are already
      !   computed, then skip.
      call checkAttributeHDF5(eigenVectors_aid(i,spinDirection),&
            & "Eigen vector SCF",hdf5Status)
      if (hdf5Status == 1) cycle

      ! Prepare the matrices.
      packedValeVale(:,:) = 0.0_double
      tempPackedValeVale(:,:) = 0.0_double

#ifndef GAMMA
      valeVale(:,:,spinDirection) = cmplx(0.0_double,0.0_double,double)
      valeValeOL(:,:) = cmplx(0.0_double,0.0_double,double)
#else
      valeValeGamma(:,:,spinDirection) = 0.0_double
      valeValeOLGamma(:,:) = 0.0_double
#endif

      if (mpiRank == 0) then

         ! Read the nuclear potential term into packed hamiltonian.
         call readPackedMatrix(atomNPOverlap_did(i),packedValeVale,&
               & packedVVDims,dim1,valeDim,0)

         ! Read the kinetic energy term into the still packed hamiltonian.
         call readPackedMatrixAccum(atomKEOverlap_did(i),packedValeVale,&
               & tempPackedValeVale,packedVVDims,0.0_double,dim1,valeDim,0)

         ! Read the mass velocity term into the still packed hamiltonian if
         !   we are doing a scalar relativistic calculation.
         if (rel == 1) then
            ! Note that the -1.0 introduce a negative sign to the term. In the
            !   future, the sign should be incorporated into the matrix
            !   calculation itself to avoid the extra work here.
            call readPackedMatrixAccum(atomMVOverlap_did(i),packedValeVale,&
                  & tempPackedValeVale,packedVVDims,0.0_double,dim1,valeDim,0)
         endif

         ! Read the atomic potential terms into the still packed hamiltonian.
         do j = 1, potDim
            call readPackedMatrixAccum(atomPotOverlap_did(i,j),packedValeVale,&
                  & tempPackedValeVale,packedVVDims,potCoeffs(j,spinDirection),&
                  & dim1,valeDim,0)
         enddo

         ! Unpack the hamiltonian matrix.
#ifndef GAMMA
         call unpackMatrix(valeVale(:,:,spinDirection),packedValeVale,valeDim,0)
      endif
      call MPI_BCAST(valeVale(:,:,spinDirection),valeDim*valeDim,&
            & MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpierr)
#else
         call unpackMatrixGamma(valeValeGamma(:,:,spinDirection),&
               & packedValeVale,valeDim,0)
      endif
      call MPI_BCAST(valeValeGamma(:,:,spinDirection),valeDim*valeDim,&
            & MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
#endif

      if (mpiRank == 0) then
         ! Read the atomic overlap matrix. 
         call readPackedMatrix(atomOverlap_did(i),packedValeVale,&
               & packedVVDims,dim1,valeDim,0)

         ! Unpack the overlap matrix.
#ifndef GAMMA
         call unpackMatrix(valeValeOL(:,:),packedValeVale,valeDim,0)
      endif
      call MPI_BCAST(valeValeOL(:,:),valeDim*valeDim,&
            & MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpierr)
#else
         call unpackMatrixGamma(valeValeOLGamma(:,:),packedValeVale,valeDim,0)
      endif
      call MPI_BCAST(valeValeOLGamma(:,:),valeDim*valeDim,&
            & MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
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
!write(20,*) "valeValeSCF"
!do j = 1,valeDim
!   do k = 1,valeDim
!      write(20,fmt="(2e12.4)",advance="NO") valeVale(k,j,spinDirection)
!   enddo
!   write(20,*)
!enddo
      call solveZHEGV(valeDim,numStates,valeVale(:,:,spinDirection),&
            & valeValeOL(:,:),energyEigenValues(:,i,spinDirection))
#else
      call solveDSYGV(valeDim,numStates,valeValeGamma(:,:,spinDirection),&
            & valeValeOLGamma(:,:),energyEigenValues(:,i,spinDirection))
#endif

!! Test normalization.
!allocate (identity(numStates,numStates))
!! Read the atomic overlap matrix. 
!call readPackedMatrix(atomOverlap_did(i),packedValeVale,&
!      & packedVVDims,dim1,valeDim)
!
!      ! Unpack the overlap matrix.
!#ifndef GAMMA
!call unpackMatrix(valeValeOL(:,:,1,1),packedValeVale,valeDim,0)
!#else
!call unpackMatrixGamma(valeValeOLGamma(:,:,1),packedValeVale,valeDim,0)
!#endif
!
!identity(:numStates,:numStates) = &
!      & matmul(&
!      & matmul(transpose(conjg(valeVale(:valeDim,:numStates,1,spinDirection))),&
!             & valeValeOL(:,:,1,1)),valeVale(:,:numStates,1,spinDirection))
!write(24,*) "Normalization i = ",i
!do k = 1, num_states
!write(24,*) "k,I(k,k) = ",k,identity(k,k)
!enddo
!
!deallocate (identity)

      ! Write the energy eigenValues onto disk in HDF5 format in a.u.
      if (mpiRank == 0) then
         call h5dwrite_f (eigenValues_did(i,spinDirection),H5T_NATIVE_DOUBLE,&
               & energyEigenValues(:,i,spinDirection),states,hdferr)
      endif
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) then
         call stopMPI("Cannot write energy eigen values.")
      endif

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
               & packedVVDims,dim1,valeDim,1)

         ! Unpack the overlap matrix.
#ifndef GAMMA
         call unpackMatrix(valeValeOL(:,:),packedValeVale,valeDim,0)
#else
         call unpackMatrixGamma(valeValeOLGamma(:,:),packedValeVale,valeDim,0)
#endif
      endif

#ifndef GAMMA
      ! Write the eigenVectors onto disk in HDF5 format for this
      !   kpoint and spin direction.
      if (mpiRank == 0) then
         call h5dwrite_f(eigenVectors_did(1,i,spinDirection),&
               & H5T_NATIVE_DOUBLE,real(valeVale(:,:numStates,spinDirection),&
               & double),valeStates,hdferr)
      endif
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) then
         call stopMPI("Cannot write real energy eigen vectors.")
      endif

      if (mpiRank == 0) then
         call h5dwrite_f(eigenVectors_did(2,i,spinDirection),&
               & H5T_NATIVE_DOUBLE,aimag(valeVale(:,:numStates,spinDirection)),&
               & valeStates,hdferr)
      endif
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) then
         call stopMPI("Cannot write imag energy eigen vectors.")
      endif
#else
      if (mpiRank == 0) then
         call h5dwrite_f(eigenVectors_did(1,i,spinDirection),&
               & H5T_NATIVE_DOUBLE,valeValeGamma(:,:numStates,spinDirection),&
               & valeStates,hdferr)
      endif
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) then
         call stopMPI("Cannot write real energy eigen vectors.")
      endif
#endif

      ! Record that this calculation is complete.
      if (mpiRank == 0) then
         call h5awrite_f(eigenVectors_aid(i,spinDirection),H5T_NATIVE_INTEGER,&
               & 1,attribIntDims,hdferr)
      endif
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) then
         call stopMPI("Failed to record eigenvector success.")
      endif
   enddo ! Loop i over kpoints.

   ! Once all kpoints are done, we need to reset the attributes so that the
   !   next iteration doesn't think that all the kpoints are done already.
   do i = 1, numKPoints
      if (mpiRank == 0) then
         call h5awrite_f(eigenVectors_aid(i,spinDirection),H5T_NATIVE_INTEGER,&
               & 0,attribIntDims,hdferr)
      endif
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) then
         call stopMPI("Failed to reset eigenvector success to 0.")
      endif
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

end subroutine secularEqnSCF


subroutine secularEqnPSCF(spinDirection,numStates,numComponents,ol_did,&
      & ham_did,eVal_did,eVec_did,eVec_aid)

   ! Import necessary modules.
   use HDF5
   use MPI_F08
   use O_MPI
   use O_Kinds
   use O_TimeStamps
   use O_KPoints, only: numKPoints
   use O_Potential, only: spin, numPlusUJAtoms, currIteration
   use O_AtomicSites, only: valeDim
   use O_PSCFEigValHDF5, only: states
   use O_PSCFEigVecHDF5, only: valeStatesPSCF
   use O_PSCFIntegralsHDF5, only: packedVVDimsPSCF
#ifndef GAMMA
   use O_LAPACKZHEGV
   use O_MatrixSubs, only: readPackedMatrix, readPackedMatrixAccum,&
         & unpackMatrix
#else
   use O_LAPACKDSYGV
   use O_MatrixSubs, only: readPackedMatrix, readPackedMatrixAccum, &
         & unpackMatrixGamma
#endif

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the passed parameters.
   integer, intent(in) :: spinDirection
   integer, intent(in) :: numStates
   integer, intent(in) :: numComponents
   integer(hid_t), dimension(numKPoints), intent(in) :: ol_did
   integer(hid_t), dimension(numKPoints,spin), intent(in) :: ham_did
   integer(hid_t), dimension(numKPoints,spin), intent(in) :: eVal_did
   integer(hid_t), dimension(numComponents,numKPoints,spin), &
         & intent(in) :: eVec_did
   integer(hid_t), dimension(numKPoints,spin), intent(in) :: eVec_aid

   ! Define the local variables used in this subroutine.
   integer :: i ! Loop index variables
!integer :: k, l, m, n, o, p
   integer :: hdferr
   integer :: hdf5Status
   integer(hsize_t), dimension (1) :: attribIntDims ! Attribute dataspace dim
   integer :: dim1
   real    (kind=double), allocatable, dimension (:,:)   :: packedValeVale
   real    (kind=double), allocatable, dimension (:,:)   :: tempPackedValeVale
!complex(kind=double), allocatable, dimension(:,:) :: identity

   ! Record the date and time that we start.
   call timeStampStart (15)

   ! Initialize the dimension of the packed matrices to include two components
   !   (real,imaginary) or just a real component.
#ifndef GAMMA
   dim1 = 2
#else
   dim1 = 1
#endif

   ! Only allocate for first spin to prevent double allocation because we do
   !   spin sequentially. All these arrays are deallocated !  . These arrays
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
      allocate(valeVale          (valeDim,valeDim,spin)) ! Complex
#else
      allocate(valeValeGamma     (valeDim,valeDim,spin)) ! Real
#endif
   endif

   ! These matrices are deallocated at the end of this subroutine.
#ifndef GAMMA
   allocate(valeValeOL         (valeDim,valeDim)) ! Complex
#else
   allocate(valeValeOLGamma    (valeDim,valeDim)) ! Real
#endif
   allocate(packedValeVale     (dim1,valeDim*(valeDim+1)/2))
   allocate(tempPackedValeVale (dim1,valeDim*(valeDim+1)/2))


   ! Begin loop over all kpoints.
   do i = 1,numKPoints

      ! If the eigenvectors for this kpoint and spin direction are already
      !   computed, then skip.
      call checkAttributeHDF5(eVec_aid(i,spinDirection),&
            & "Eigen vector PSCF",hdf5Status)
      if (hdf5Status == 1) cycle

      ! Prepare the matrices.
      packedValeVale(:,:) = 0.0_double
      tempPackedValeVale(:,:) = 0.0_double

#ifndef GAMMA
      valeVale(:,:,spinDirection) = cmplx(0.0_double,0.0_double,double)
      valeValeOL(:,:) = cmplx(0.0_double,0.0_double,double)
#else
      valeValeGamma(:,:,spinDirection) = 0.0_double
      valeValeOLGamma(:,:) = 0.0_double
#endif

      ! Read the combined Hamiltonian term into the packed valeVale.
      call readPackedMatrix(ham_did(i,spinDirection),&
            & packedValeVale,packedVVDimsPSCF,dim1,valeDim,1)

      ! Unpack the hamiltonian matrix.
#ifndef GAMMA
      call unpackMatrix(valeVale(:,:,spinDirection),packedValeVale,valeDim,0)
#else
      call unpackMatrixGamma(valeValeGamma(:,:,spinDirection),&
            & packedValeVale,valeDim,0)
#endif

      ! Read the atomic overlap matrix. 
      call readPackedMatrix(ol_did(i),packedValeVale,packedVVDimsPSCF,dim1,&
            & valeDim,1)

      ! Unpack the overlap matrix.
#ifndef GAMMA
      call unpackMatrix(valeValeOL(:,:),packedValeVale,valeDim,0)
#else
      call unpackMatrixGamma(valeValeOLGamma(:,:),packedValeVale,valeDim,0)
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
!write(20,*) "valeValePSCF"
!do j = 1,valeDim
!   do k = 1,valeDim
!      write(20,fmt="(2e12.4)",advance="NO") valeVale(k,j,spinDirection)
!   enddo
!   write(20,*)
!enddo
      call solveZHEGV(valeDim,numStates,valeVale(:,:,spinDirection),&
            & valeValeOL(:,:),energyEigenValues(:,i,spinDirection))
#else
      call solveDSYGV(valeDim,numStates,valeValeGamma(:,:,spinDirection),&
            & valeValeOLGamma(:,:),energyEigenValues(:,i,spinDirection))
#endif

!! Test normalization.
!allocate (identity(numStates,numStates))
!! Read the atomic overlap matrix. 
!call readPackedMatrix(atomOverlap_did(i),packedValeVale,&
!      & packedVVDims,dim1,valeDim)
!
!      ! Unpack the overlap matrix.
!#ifndef GAMMA
!call unpackMatrix(valeValeOL(:,:,1,1),packedValeVale,valeDim,0)
!#else
!call unpackMatrixGamma(valeValeOLGamma(:,:,1),packedValeVale,valeDim,0)
!#endif
!
!identity(:numStates,:numStates) = &
!      & matmul(&
!      & matmul(transpose(conjg(valeVale(:valeDim,:numStates,1,spinDirection))),&
!             & valeValeOL(:,:,1,1)),valeVale(:,:numStates,1,spinDirection))
!write(24,*) "Normalization i = ",i
!do k = 1, num_states
!write(24,*) "k,I(k,k) = ",k,identity(k,k)
!enddo
!
!deallocate (identity)

      ! Write the energy eigenValues onto disk in HDF5 format in a.u.
      if (mpiRank == 0) then
         call h5dwrite_f (eVal_did(i,spinDirection),H5T_NATIVE_DOUBLE,&
               & energyEigenValues(:,i,spinDirection),states,hdferr)
      endif
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) then
         call stopMPI("Cannot write energy eigen values.")
      endif

      ! If we have some atoms with plusUJ terms, we need to update the just
      !   obtained eVectors. The update depends on the charge in each of the
      !   highest d or f orbitals of the affected atoms. Therefore, we will
      !   need to reread in the overlap matrix again because it was just
      !   destroyed in the ZHEGV solution. An alternative approach might be
      !   to make a copy of the overlap matrix and hold it in reserve until we
      !   need it for this calculation.
      if (numPlusUJAtoms > 0) then

         ! Read the atomic overlap matrix. 
         call readPackedMatrix(ol_did(i),packedValeVale,packedVVDimsPSCF,dim1,&
               & valeDim,1)

         ! Unpack the overlap matrix.
#ifndef GAMMA
         call unpackMatrix(valeValeOL(:,:),packedValeVale,valeDim,0)
#else
         call unpackMatrixGamma(valeValeOLGamma(:,:),packedValeVale,valeDim,0)
#endif
      endif

#ifndef GAMMA
      ! Write the eigenVectors onto disk in HDF5 format for this
      !   kpoint and spin direction.
      if (mpiRank == 0) then
         call h5dwrite_f(eVec_did(1,i,spinDirection),&
               & H5T_NATIVE_DOUBLE,real(valeVale(:,:numStates,spinDirection),&
               & double), valeStatesPSCF,hdferr)
      endif
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) then
         call stopMPI("Cannot write real energy eigen vectors.")
      endif

      if (mpiRank == 0) then
         call h5dwrite_f(eVec_did(2,i,spinDirection),&
               & H5T_NATIVE_DOUBLE,aimag(valeVale(:,:numStates,spinDirection)),&
               & valeStatesPSCF,hdferr)
      endif
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) then
         call stopMPI("Cannot write imag energy eigen vectors.")
      endif
#else
      if (mpiRank == 0) then
         call h5dwrite_f(eVec_did(1,i,spinDirection),&
               & H5T_NATIVE_DOUBLE,real(valeValeGamma(:,:numStates,&
               & spinDirection),double),valeStatesPSCF,hdferr)
      endif
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) then
         call stopMPI("Cannot write real energy eigen vectors.")
      endif
#endif

      ! Record that this kpoint has been finished.
      if (mpiRank == 0) then
         if (mod(i,10) .eq. 0) then
            write (20,ADVANCE="NO",FMT="(a1)") "|"
         else
            write (20,ADVANCE="NO",FMT="(a1)") "."
         endif
         if (mod(i,50) .eq. 0) then
            write (20,*) " ",i
         endif
         call flush (20)
      endif

      ! Record that this calculation is complete.
      if (mpiRank == 0) then
         call h5awrite_f(eVec_aid(i,spinDirection),H5T_NATIVE_INTEGER,&
               & 1,attribIntDims,hdferr)
      endif
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) then
         call stopMPI("Failed to record eigenvector success.")
      endif

      if (mpiRank == 0) then
         call h5aclose_f(eVec_aid(i,spinDirection),hdferr)
      endif
      call MPI_BCAST(hdferr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (hdferr /= 0) then
         call stopMPI("Failed to close eigenvector attribute.")
      endif
   enddo ! Loop i over kpoints.

   ! Deallocate unnecessary arrays and matrices.
   deallocate (packedValeVale)
   deallocate (tempPackedValeVale)
#ifndef GAMMA
   deallocate(valeVale)
   deallocate(valeValeOL)
#else
   deallocate(valeValeGamma)
   deallocate(valeValeOLGamma)
#endif

   ! Record the date and time that we finish.
   call timeStampEnd (15)

end subroutine secularEqnPSCF


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
         & plusUJ, spin

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
   use O_KPoints, only: kPointWeight
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
   complex (kind=double) :: sum1, sum2 ! Two intermediate summations used in
         ! thencomputation of the plusUJ term for localized d- or f-orbitals
         ! in select atoms.

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
         sum1 = cmplx(0.0_double,0.0_double,double)
         do k = 1, currUJSize
#ifndef GAMMA
            sum1 = sum1 + plusUJ(k,i,oppositeSpin,currKPoint) * &
                  & valeValeOL(currUJIndex+k,currUJIndex+k) - &
                  & plusUJAtomGSElectrons(i) / (2.0_double * &
                  & real(currUJSize,double)) * kPointWeight(currKPoint) / &
                  & real(spin,double)
#else
            sum1 = sum1 + plusUJ(k,i,oppositeSpin,currKPoint) * &
                  & valeValeOLGamma(currUJIndex+k,currUJIndex+k) - &
                  & plusUJAtomGSElectrons(i) / (2.0_double * &
                  & real(currUJSize,double)) * kPointWeight(currKPoint) / &
                  & real(spin,double)
#endif
         enddo
         sum1 = sum1 * plusUJAtomValue(1,i)

         ! Compute the second temporary sum for this orbital of this atom.
         sum2 = cmplx(0.0_double,0.0_double,double)
         do k = 1, currUJSize

            ! Don't include m=m' terms in the summation.
            if (k == j) then
               cycle
            endif

            ! Accumulate summation terms.
#ifndef GAMMA
            sum2 = sum2 + plusUJ(k,i,spinDirection,currKPoint) * &
                  & valeValeOL(currUJIndex+k,currUJIndex+k) - &
                  & plusUJAtomGSElectrons(i) / (2.0_double * &
                  & real(currUJSize,double)) * kPointWeight(currKPoint) / &
                  & real(spin,double)
#else
            sum2 = sum2 + plusUJ(k,i,spinDirection,currKPoint) * &
                  & valeValeOLGamma(currUJIndex+k,currUJIndex+k) - &
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
         valeVale(currUJIndex+j,currUJIndex+j,spinDirection) = &
               & valeVale(currUJIndex+j,currUJIndex+j,spinDirection) + &
               & sum1 + sum2
#else
         valeValeGamma(currUJIndex+j,currUJIndex+j,spinDirection) = &
               & valeValeGamma(currUJIndex+j,currUJIndex+j,spinDirection) + &
               & real(sum1,double) + real(sum2,double)
#endif
      enddo
   enddo

end subroutine update2AndApplyUJ


subroutine shiftEnergyEigenValues(energyShift)

!use O_KPoints
   ! Make sure that no funny variables are defined.
   implicit none

   ! Define dummy variables passed to this subroutine.
   real (kind=double) :: energyShift
!integer :: i,j,k

   ! Shift the energyEigenValues down by the requested about.
   energyEigenValues(:,:,:) = energyEigenValues(:,:,:) - energyShift

!write(20,*) "energyShift:  ", energyShift
!do i = 1, 1
!do j = 1, numKPoints
!do k = 1, numStates
!write(20,*) i,j,k, energyEigenValues(k,j,i)
!enddo
!enddo
!enddo

end subroutine shiftEnergyEigenValues


subroutine readDataSCF(h,i,numStates,matrixCode)

   ! Import necessary data modules.
   use MPI_F08
   use O_MPI
   !use O_KPoints, only: numKPoints
   use O_AtomicSites, only: valeDim
   use O_SCFIntegralsHDF5, only: packedVVDims,atomOverlap_did,atomMMOverlap_did
#ifndef GAMMA
   use O_SCFEigVecHDF5, only: valeStates,eigenVectors_did
   use O_MatrixSubs, only: readMatrix,readPackedMatrix,unpackMatrix
#else
   use O_SCFEigVecHDF5, only: valeStates
   use O_MatrixSubs, only: readPackedMatrix,unpackMatrixGamma
#endif

   ! Define passed parameters.
   integer, intent(in) :: h ! Spin variable.
   integer, intent(in) :: i ! KPoint variable
   integer, intent(in) :: numStates
   integer, intent(in) :: matrixCode

   ! Define local variables.
   integer :: dim1
   integer :: j ! Loop index (usually xyz).
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

   ! Avoid compiler warning about unused variables.
#ifdef GAMMA
   if (h > 2) stop
   if (numStates > valeDim) stop
   valeStates = valeStates
#endif

   if (matrixCode > 0) then

      ! Allocate space to read a packed matrix.
      allocate (packedValeVale(dim1,valeDim*(valeDim+1)/2))

      if (matrixCode == 1) then
         if (mpiRank == 0) then

            ! Read the overlap matrix.  The tempPackedMatrix is not used.
            call readPackedMatrix(atomOverlap_did(i),packedValeVale,&
                  & packedVVDims,dim1,valeDim,0)

            ! Unpack the matrix.
#ifndef GAMMA
            call unpackMatrix(valeValeOL(:,:),packedValeVale,valeDim,1)
         endif
         call MPI_BCAST(valeValeOL(:,:),valeDim*valeDim,&
               & MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpierr)
#else
            call unpackMatrixGamma(valeValeOLGamma(:,:),packedValeVale,&
                  & valeDim,1)
         endif
         call MPI_BCAST(valeValeOLGamma(:,:),valeDim*valeDim,&
               & MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
#endif
      elseif (matrixCode == 2) then
         if (mpiRank == 0) then
            do j = 1, 3
               ! Read the xyz momentum matrix elements.
               call readPackedMatrix(atomMMOverlap_did(i,j),packedValeVale,&
                     & packedVVDims,dim1,valeDim,0)

               ! Unpack the matrix.
#ifndef GAMMA
               call unpackMatrix(valeValeMM(:,:,j),packedValeVale,valeDim,1)
#else
               call unpackMatrixGamma(valeValeMMGamma(:,:,j),packedValeVale,&
                   & valeDim,1)
#endif
            enddo
         endif
#ifndef GAMMA
         call MPI_BCAST(valeValeMM(:,:,:),valeDim*valeDim*3,&
               & MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpierr)
#else
         call MPI_BCAST(valeValeMMGamma(:,:,:),valeDim*valeDim*3,&
               & MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
#endif
      endif

      ! Deallocate the packed matrix used to read and unpack the data.
      deallocate (packedValeVale)

   endif

#ifndef GAMMA
   ! Read the wave functions for this kpoint from the datasets into
   !   the valeVale matrix.  If numKPoints==1, the wave functions should
   !   already be in the valeVale(1:valeDim,1:numStates) matrix.
!   if (numKPoints > 1) then
   if (mpiRank == 0) then

      ! Allocate space to read the wave functions.
      allocate (tempRealValeVale (valeDim,numStates))
      allocate (tempImagValeVale (valeDim,numStates))

      call readMatrix(eigenVectors_did(:,i,h),valeVale(:,:numStates,h),&
            & tempRealValeVale(:,:),tempImagValeVale(:,:),&
            & valeStates,valeDim,numStates,0) ! No internal broadcast.

      ! Deallocate the space to read the wave functions.
      deallocate (tempRealValeVale)
      deallocate (tempImagValeVale)
   endif
   call MPI_BCAST(valeVale(:,:numStates,h),valeDim*numStates,&
         & MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpierr)
!   endif

#endif

end subroutine readDataSCF



subroutine readDataPSCF(h,i,numStates,matrixCode)

   ! Use necessary modules.
   use MPI_F08
   use O_MPI
   !use O_KPoints, only: numKPoints
   use O_AtomicSites, only: valeDim
   use O_PSCFIntegralsHDF5, only: packedVVDimsPSCF,atomOverlapPSCF_did,&
         & atomMMOverlapPSCF_did
   use O_PSCFEigVecHDF5, only: valeStatesPSCF,eigenVectorsPSCF_did
#ifndef GAMMA
   use O_MatrixSubs, only: readMatrix, readPackedMatrix, unpackMatrix
#else
   use O_MatrixSubs, only: readMatrixGamma, readPackedMatrix, &
         & unpackMatrixGamma
#endif

   ! Define passed parameters.
   integer, intent(in) :: h ! Spin variable.
   integer, intent(in) :: i ! KPoint variable
   integer, intent(in) :: numStates
   integer, intent(in) :: matrixCode

   ! Define local variables.
   integer :: dim1
   integer :: j ! Loop index (usually xyz).
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

   if (matrixCode > 0) then

      ! Allocate space to read a packed matrix.
      allocate (packedValeVale(dim1,valeDim*(valeDim+1)/2))

      if (matrixCode == 1) then
         if (mpiRank == 0) then
            ! Read the overlap matrix.  The tempPackedMatrix is not used.
            call readPackedMatrix(atomOverlapPSCF_did(i),packedValeVale,&
                  & packedVVDimsPSCF,dim1,valeDim,0)

            ! Unpack the matrix.
#ifndef GAMMA
            call unpackMatrix(valeValeOL(:,:),packedValeVale,valeDim,1)
         endif
         call MPI_BCAST(valeValeOL(:,:),valeDim*valeDim,&
               & MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpierr)
#else
            call unpackMatrixGamma(valeValeOLGamma(:,:),packedValeVale,valeDim,1)
         endif
         call MPI_BCAST(valeValeOLGamma(:,:),valeDim*valeDim,&
               & MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
#endif

      elseif (matrixCode == 2) then
         if (mpiRank == 0) then
            do j = 1, 3
               ! Read the xyz momentum matrix elements.
!write(20,*) i,j,shape(packedValeVale),packedVVDimsPSCF,dim1,valeDim
               call readPackedMatrix(atomMMOverlapPSCF_did(i,j),packedValeVale,&
                     & packedVVDimsPSCF,dim1,valeDim,0)
#ifndef GAMMA
!write(20,*) i,packedValeVale(:,:)
#endif

               ! Unpack the matrix.
#ifndef GAMMA
               call unpackMatrix(valeValeMM(:,:,j),packedValeVale,valeDim,1)
#else
               call unpackMatrixGamma(valeValeMMGamma(:,:,j),packedValeVale,&
                     & valeDim,1)
#endif
            enddo
         endif
#ifndef GAMMA
         call MPI_BCAST(valeValeMM(:,:,:),valeDim*valeDim*3,&
               & MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpierr)
#else
         call MPI_BCAST(valeValeMMGamma(:,:,:),valeDim*valeDim*3,&
               & MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
#endif
      endif

      ! Deallocate the packed matrix used to read and unpack the data.
      deallocate (packedValeVale)
   endif

#ifndef GAMMA
   ! Read the wave functions for this kpoint from the datasets into
   !   the valeVale matrix.  If numKPoints==1, the wave functions should
   !   already be in the valeVale(1:valeDim,1:numStates) matrix, unless it
   !   was already computed once in which case we just have not yet read it
   !   in (same for the eigen values).
!   if (numKPoints > 1) then

   if (mpiRank == 0) then
      ! Allocate space to read the complex wave function.
      allocate (tempRealValeVale (valeDim,numStates))
      allocate (tempImagValeVale (valeDim,numStates))

      ! Read the complex wave function from the datasets.
      call readMatrix(eigenVectorsPSCF_did(:,i,h),valeVale(:,:numStates,h),&
         & tempRealValeVale(:,:),tempImagValeVale(:,:),&
         & valeStatesPSCF,valeDim,numStates,0) ! No internal broadcast

      ! Deallocate the space to read the complex wave function.
      deallocate (tempRealValeVale)
      deallocate (tempImagValeVale)
   endif
   call MPI_BCAST(valeVale(:,:numStates,h),valeDim*numStates,&
         & MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpierr)
!   endif
#else
   if (mpiRank == 0) then
      ! Read the real wave function from the datasets.
      call readMatrixGamma(eigenVectorsPSCF_did(1,i,h),&
            & valeValeGamma(:,:numStates,h),valeStatesPSCF,valeDim,&
            & numStates,0) ! No internal broadcast
   endif
   call MPI_BCAST(valeValeGamma(:,:numStates,h),valeDim*numStates,&
         & MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
#endif

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
