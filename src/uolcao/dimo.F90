module O_DIMO

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no variables are declared accidentally.
   implicit none

   contains

subroutine computeDIMO (inSCF)

   ! Import necessary modules.
   use O_TimeStamps
   use O_AtomicSites, only: valeDim, computeIonicMoment, xyzIonMoment
   use O_Input, only: numStates, numElectrons
   use O_KPoints, only: numKPoints
   use O_Constants, only: smallThresh, eCharge, bohrRad
   use O_Potential, only: spin
   use O_Populate, only: electronPopulation,cleanUpPopulation
   use O_SCFIntegralsHDF5, only: atomOverlap_did,atomDMOverlap_did, &
         & packedVVDims
   use O_PSCFIntegralsHDF5, only: atomDMOverlapPSCF_did,atomOverlapPSCF_did,&
         & packedVVDimsPSCF
   use O_Lattice, only: realCellVolume
#ifndef GAMMA
   use O_BLASZHER
   use O_MatrixSubs, only: readPackedMatrix,matrixElementMult,packMatrix
   use O_SecularEquation, only: valeVale,energyEigenValues,&
         & readDataSCF, readDataPSCF
#else
   use O_BLASDSYR
   use O_MatrixSubs, only: readPackedMatrix, &
         & matrixElementMultGamma,packMatrixGamma
   use O_SecularEquation, only: valeValeGamma,energyEigenValues,&
         & readDataSCF, readDataPSCF
#endif

   ! Define passed parameters.
   integer :: inSCF

   ! Define local variables.
   integer :: h,i,j,k
   integer :: dim1
   integer :: skipKP
   integer :: energyLevelCounter
   real (kind=double), allocatable, dimension (:) :: tempDensity
   real (kind=double), allocatable, dimension (:,:) :: dipoleMomentTrace
   real (kind=double), allocatable, dimension (:) :: currentPopulation
   real (kind=double), allocatable, dimension (:,:) :: packedValeVale
   real (kind=double), allocatable, dimension (:,:,:) :: packedValeValeRho
   real (kind=double), allocatable, dimension (:,:,:) :: &
         & structuredElectronPopulation
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: valeValeRho
#else
   real    (kind=double), allocatable, dimension (:,:,:) :: valeValeRhoGamma
#endif

   ! Log the date and time that we start.
   call timeStampStart (34)

   ! Set working parameters.
#ifndef GAMMA
   dim1 = 2
#else
   dim1 = 1
#endif

   ! Allocate working space
#ifndef GAMMA
   if (inSCF == 0) then
      allocate (valeVale(valeDim,numStates,spin))
   endif
   allocate (valeValeRho(valeDim,valeDim,spin))
#else
   if (inSCF == 0) then
      allocate (valeValeGamma(valeDim,numStates,spin))
   endif
   allocate (valeValeRhoGamma(valeDim,valeDim,spin))
#endif
   allocate (tempDensity(dim1))
   allocate (currentPopulation (spin))
   allocate (dipoleMomentTrace(3,spin))
   allocate (structuredElectronPopulation (numStates,numKPoints,spin))
   allocate (packedValeVale(dim1,valeDim*(valeDim+1)/2))
   allocate (packedValeValeRho(dim1,valeDim*(valeDim+1)/2,spin))

   ! Initialize arrays.
#ifndef GAMMA
   valeVale(:,:,:) = cmplx(0.0_double,0.0_double,double)
#else
   valeValeGamma(:,:,:) = 0.0_double
#endif
   dipoleMomentTrace(:,:) = 0.0_double

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
         skipKP = 1 ! Assume that we will skip this kpoint.
         do j = 1, numStates
            if (sum(abs(structuredElectronPopulation(j,i,:)))>smallThresh) then
               skipKP = 0 ! Enough contribution to not skip.
               exit
            endif
         enddo
         if (skipKP == 1) then
            cycle
         endif

         ! Determine if we are doing the valeCharge in a post-SCF calculation
         !   or within an SCF calculation.
         do h = 1, spin
            if (inSCF == 1) then
               call readDataSCF(h,i,numStates,0) ! Read wave functions only.
            else
               call readDataPSCF(h,i,numStates,0) ! Read wave functions only.
            endif
         enddo
!      endif

      ! Initialize matrix to receive the valeVale density matrix (square of the
      !   wave function).
      valeValeRho(:,:,:) = cmplx(0.0_double,0.0_double,double)
#else
      ! All of the information we need is already available in system memory.
      !   The only thing we need to do is initialize this matrix to zero.
      valeValeRhoGamma(:,:,:) = 0.0_double
#endif


      ! Accumulate the valeValeRho matrix upper triangle and electron energy.
      do j = 1, numStates
         currentPopulation(:) = structuredElectronPopulation(j,i,:)

         if (sum(abs(currentPopulation(:))) < smallThresh) cycle

         do k = 1, spin
#ifndef GAMMA
            call zher('U',valeDim,currentPopulation(k),valeVale(:,j,k),1,&
                  & valeValeRho(:,:,k),valeDim)
#else
            call dsyr('U',valeDim,currentPopulation(k),&
                  & valeValeGamma(:,j,k),1,valeValeRhoGamma(:,:,k),valeDim)
#endif
         enddo
      enddo

      ! Pack the matrix for easy comparison with the hamiltonian terms
      !   to be read in next.  Store the result in the appropriate packed
      !   valeVale spin array.
      do j = 1, spin
#ifndef GAMMA
         call packMatrix(valeValeRho(:,:,j),packedValeValeRho(:,:,j),&
               & valeDim)
#else
         call packMatrixGamma(valeValeRhoGamma(:,:,j),&
               & packedValeValeRho(:,:,j),valeDim)
#endif
      enddo

      ! Read the overlap matrix into the packedValeVale representation.
      if (inSCF == 1) then
         call readPackedMatrix (atomOverlap_did(i),packedValeVale,&
               & packedVVDims,dim1,valeDim)
      else
         call readPackedMatrix (atomOverlapPSCF_did(i),packedValeVale,&
               & packedVVDimsPSCF,dim1,valeDim)
      endif

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


      do j = 1, 3 ! xyz directions

         if (inSCF == 1) then
            call readPackedMatrix (atomDMOverlap_did(i,j),packedValeVale,&
                  & packedVVDims,dim1,valeDim)
         else
            call readPackedMatrix (atomDMOverlapPSCF_did(i,j),&
                  & packedValeVale,packedVVDimsPSCF,dim1,valeDim)
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
   enddo ! kpoint

   ! Compute the nuclear contribution to the dipole
   call computeIonicMoment(0)

   ! The units of the dipole moment are e * a_0 where e is the electronic
   !   charge and a_0 is the bohr radius.
   ! The units of the dipole moment calculation are ea0 where e is the number
   !    of valence electrons and a0 is the bohr radius.
   do i = 1, spin
      if ((i == 1) .and. (spin == 1)) then
         write(74,fmt="(a)") "Total:"
      elseif (i == 1) then
         write(74,fmt="(a)") "Spin Up:"
      else
         write(75,fmt="(a)") "Spin Dn:"
      endif

      write (73+i,*) "Elec Mom (x,y,z): ", dipoleMomentTrace(:,i)
      write (73+i,*) "Elec. Polarization (x,y,z) [C/m^2]: ", eCharge * &
            & dipoleMomentTrace(:,i)/realCellVolume * 10.0_double / &
            & bohrRad**2
      write (74,*) "Nuclear Moment (x,y,z):    ", xyzIonMoment(:)
      write (74,*) "Dipole Moment (a.u.):    ", (xyzIonMoment(:) - &
            & dipoleMomentTrace(:,i))
      write (74,*) "Dipole Moment (Debye):    ", (xyzIonMoment(:) - &
            & dipoleMomentTrace(:,i)) * 2.541746473_double
      write (74,*) "Total Dipole Moment (Debye):    ", &
            & sqrt(sum((xyzIonMoment(:) - dipoleMomentTrace(:,i))**2)) &
            * 2.541746473_double
      write (74,*) "Polarization (x,y,z) [C/m^2]:    ", eCharge * &
            & (xyzIonMoment(:) - dipoleMomentTrace(:,i)) / &
            & realCellVolume * 10.0_double / bohrRad**2
      write (74,*) "Total Polarization [C/m^2]:    ", eCharge * &
            & sqrt(sum((xyzIonMoment(:) - dipoleMomentTrace(:,i))**2)) / &
            & realCellVolume * 10.0_double / bohrRad**2
   enddo

   ! Deallocate arrays associated with the electron population.
   call cleanUpPopulation

   deallocate (structuredElectronPopulation)
   deallocate (packedValeValeRho)
   deallocate (currentPopulation)
   deallocate (tempDensity)
#ifndef GAMMA
   deallocate (valeValeRho)
#else
   deallocate (valeValeRhoGamma)
#endif

   ! Log the date and time we end.
   call timeStampEnd (34)

end subroutine computeDIMO

end module O_DIMO
