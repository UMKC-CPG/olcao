module O_SCFEigVecHDF5

   ! Use the HDF5 module (for hsize_t and hid_t).
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Declare array that holds the dimensions of the eigenvector dataset.
   integer(hsize_t), dimension (2) :: valeStates

   ! Declare array that holds the dimensions of the chunk.
   integer(hsize_t), dimension (2) :: valeStatesChunk

   ! Declare the eigenvector subgroup of the scf_fid.
   integer(hid_t) :: eigenVectors_gid

   ! Declare the eigenvector dataspace.
   integer(hid_t) :: valeStates_dsid

   ! Declare the property list for the eigenvector data space.
   integer(hid_t) :: valeStates_plid

   ! Declare the eigenvector dataset array.  The number of datasets will depend
   !   on the number of kpoints and so it will vary.  Those datasets will be
   !   given ID numbers dynamically.  (Index1=real,imag; Index2=1..kpoints;
   !   Index3=1..spin)
   integer(hid_t), allocatable, dimension (:,:,:) :: eigenVectors_did

   ! Declare the attribute IDs that will be used to track completion of the
   !   eigenvector calculations. We need one for each kpoint and spin (because
   !   the real and imaginary will always be completed at the same time).
   integer(hid_t), allocatable, dimension (:,:) :: eigenVectors_aid

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


subroutine initSCFEigVecHDF5 (scf_fid,attribInt_dsid,attribIntDims,numStates)

   ! Import any necessary definition modules.
   use HDF5
   use O_Kinds
   use O_MPI

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints
   use O_AtomicSites, only: valeDim
   use O_Potential, only: spin

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Define the passed parameters.
   integer(hid_t) :: scf_fid
   integer(hid_t) :: attribInt_dsid
   integer(hsize_t), dimension (1) :: attribIntDims
   integer, intent(in) :: numStates

   ! Define local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize data structure dimensions.
   valeStates(1)    = valeDim
   valeStates(2)    = numStates

   ! Check that the chunk size is not too large (the assumption here is that
   !   the number being stored are 8 byte reals and that we should not go over
   !   2 billion bytes. Note that if x*y = >250M and we want a*b = 250M then
   !   the additional requirement x/y = a/b leads to b = sqrt(250M/>250M)*y.
   !   Thus a = 250M/b.
   if (valeStates(1) * valeStates(2) > 250000000) then
      valeStatesChunk(2) = int(sqrt(real(250000000,double) / &
            & real(valeStates(1) * valeStates(2),double)) * valeStates(2))
      valeStatesChunk(1) = int(250000000 / valeStatesChunk(2))
   else
      valeStatesChunk(1) = valeStates(1)
      valeStatesChunk(2) = valeStates(2)
   endif

   ! Allocate space to hold IDs for the datasets in the eigenvectors group.
#ifndef GAMMA 
   allocate (eigenVectors_did(2,numKPoints,spin)) ! Complex
#else
   allocate (eigenVectors_did(1,numKPoints,spin)) ! Real
#endif

   ! Allocate space to hold the attributes for tracking completion.
   allocate (eigenVectors_aid(numKPoints,spin))

   ! Only process 0 opens the HDF5 structure.
   if (mpiRank /= 0) return

   ! Create the eigenVectors group in the scf_fid.
   call h5gcreate_f (scf_fid,"eigenVectors",eigenVectors_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create eigenvectors group SCF'

   ! Create the dataspace that will be used for the energy eigen vectors.
   call h5screate_simple_f(2,valeStates,valeStates_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create valeStates dataspace SCF'

   ! Create the property list first.  Then set the properties one at a time.
   call h5pcreate_f(H5P_DATASET_CREATE_F,valeStates_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create valeStates plid SCF'
   call h5pset_layout_f(valeStates_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set valeStates plid layout as chunked SCF'
   call h5pset_chunk_f(valeStates_plid,2,valeStatesChunk,hdferr)
   if (hdferr /= 0) stop 'Failed to set valeStates plid chunk size SCF'
!   call h5pset_shuffle_f(valeStates_plid,hdferr)
   call h5pset_deflate_f   (valeStates_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set valeStates for deflation SCF'

   ! Create the datasets for the eigenvectors.
   do i = 1, spin
      do j = 1, numKPoints
#ifndef GAMMA
         write (currentName,fmt="(a4,i7.7,i7.7)") "real",j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectors_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeStates_dsid,eigenVectors_did(1,j,i),hdferr,&
               & valeStates_plid)
         if (hdferr /= 0) stop 'Failed to create eigenVectors real did SCF'
         write (currentName,fmt="(a4,i7.7,i7.7)") "imag",j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectors_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeStates_dsid,eigenVectors_did(2,j,i),hdferr,&
               & valeStates_plid)
         if (hdferr /= 0) stop 'Failed to create eigenVectors imaginary did SCF'
#else
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectors_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeStates_dsid,eigenVectors_did(1,j,i),hdferr,&
               & valeStates_plid)
         if (hdferr /= 0) stop 'Failed to create eigenVectors real did SCF'
#endif
      enddo
   enddo

   ! Create the eigenvector tracking attributes.
   do i = 1, spin
      do j = 1, numKPoints
         call h5acreate_f (eigenVectors_did(1,j,i),"status",&
               & H5T_NATIVE_INTEGER,attribInt_dsid,eigenVectors_aid(j,i),&
               & hdferr)
            if (hdferr /= 0) stop 'Failed to create eigenvectos_aid SCF.'
      enddo
   enddo

   ! Init eigenvector tracking attributes to uncomputed (status=0) state.
   do i = 1, spin
      do j = 1, numKPoints
         call h5awrite_f (eigenVectors_aid(j,i),H5T_NATIVE_INTEGER,0,&
               & attribIntDims,hdferr)
            if (hdferr /= 0) stop 'Failed to init eigenvectos_aid SCF.'
      enddo
   enddo

end subroutine initSCFEigVecHDF5 


subroutine accessSCFEigVecHDF5 (scf_fid,numStates)

   ! Import any necessary definition modules.
   use HDF5
   use O_Kinds
   use O_MPI

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints
   use O_AtomicSites, only: valeDim
   use O_Potential, only: spin

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Define the passed parameters.
   integer(hid_t) :: scf_fid
   integer, intent(in) :: numStates

   ! Define local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize data structure dimensions.
   valeStates(1)    = valeDim
   valeStates(2)    = numStates

   ! Check that the chunk size is not too large (the assumption here is that
   !   the numbers being stored are 8 byte reals and that we should not go
   !   over 2 billion bytes. Note that if x*y = >250M and we want a*b = 250M
   !   then the additional requirement x/y = a/b leads to:
   !   b = sqrt(250M/>250M)*y. Thus a = 250M/b.
   if (valeStates(1) * valeStates(2) > 250000000) then
      valeStatesChunk(2) = int(sqrt(real(250000000,double) / &
            & real(valeStates(1) * valeStates(2),double)) * valeStates(2))
      valeStatesChunk(1) = int(250000000 / valeStatesChunk(2))
   else
      valeStatesChunk(1) = valeStates(1)
      valeStatesChunk(2) = valeStates(2)
   endif

   ! Allocate space to hold IDs for the datasets in the eigenvectors group.
#ifndef GAMMA 
   allocate (eigenVectors_did(2,numKPoints,spin)) ! Complex
#else
   allocate (eigenVectors_did(1,numKPoints,spin)) ! Real
#endif

   ! Allocate space to hold the attributes for tracking completion.
   allocate (eigenVectors_aid(numKPoints,spin))

   ! Only process 0 opens the HDF5 structure.
   if (mpiRank /= 0) return

   ! Open the eigenVectors group in the scf_fid.
   call h5gopen_f (scf_fid,"/eigenVectors",eigenVectors_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open eigenvectors group SCF.'

   ! Open the datasets for the eigenvectors.
   do i = 1, spin
      do j = 1, numKPoints
#ifndef GAMMA
         write (currentName,fmt="(a4,i7.7,i7.7)") "real",j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenVectors_gid,currentName,&
               & eigenVectors_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVectors real did SCF.'
         write (currentName,fmt="(a4,i7.7,i7.7)") "imag",j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenVectors_gid,currentName,&
               & eigenVectors_did(2,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVectors imaginary did SCF.'
#else
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenVectors_gid,currentName,&
               & eigenVectors_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVectors real did SCF.'
#endif
      enddo
   enddo

   ! Obtain the property list for the energy eigen vectors.
   call h5dget_create_plist_f(eigenVectors_did(1,1,1),valeStates_plid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to obtain valeStates property list SCF.'

   ! Obtain the dataspace for the energy eigen vectors.
   call h5dget_space_f(eigenVectors_did(1,1,1),valeStates_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain valeStates dataspace SCF.'

   ! Open the eigenvector tracking attributes.
   do i = 1, spin
      do j = 1, numKPoints
         call h5aopen_f (eigenVectors_did(1,j,i),"status",&
               & eigenVectors_aid(j,i),hdferr)
            if (hdferr /= 0) stop 'Failed to open eigenvectos_aid SCF.'
      enddo
   enddo

end subroutine accessSCFEigVecHDF5 


subroutine closeSCFEigVecHDF5

   ! Import any necessary definition modules.
   use HDF5
   use O_MPI

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints
   use O_Potential, only: spin

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Declare local loop variables and error control.
   integer :: i,j
   integer :: hdferr

   if (mpiRank == 0) then

      ! Close the eigenvector dataspace.
      call h5sclose_f (valeStates_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to close valeStates_dsid SCF.'

      ! Close the eigenvector datasets next.
      do i = 1, spin
         do j = 1, numKPoints
#ifndef GAMMA
            call h5dclose_f (eigenVectors_did(1,j,i),hdferr)
            if (hdferr /= 0) stop 'Failed to close eigenVectors_did real SCF'
            call h5dclose_f (eigenVectors_did(2,j,i),hdferr)
            if (hdferr /= 0) stop &
                  & 'Failed to close eigenVectors_did imaginary SCF'
#else
            call h5dclose_f (eigenVectors_did(1,j,i),hdferr)
            if (hdferr /= 0) stop 'Failed to close eigenVectors_did real SCF'
#endif
         enddo
      enddo

      ! Closed the eigen vector attributes.
      do i = 1, spin
         do j = 1, numKPoints
            call h5aclose_f(eigenVectors_aid(j,i),hdferr)
            if (hdferr /= 0) stop 'Failed to close eigenVectors_aid SCF'
         enddo
      enddo

      ! Close the eigenvector property list.
      call h5pclose_f (valeStates_plid,hdferr)
      if (hdferr /= 0) stop 'Failed to close valeStates_plid SCF.'

      ! Close the eigenvector group.
      call h5gclose_f (eigenVectors_gid,hdferr)
   endif

   ! Deallocate dataset and attribute arrays.
   deallocate (eigenVectors_did)
   deallocate (eigenVectors_aid)

end subroutine closeSCFEigVecHDF5


end module O_SCFEigVecHDF5
