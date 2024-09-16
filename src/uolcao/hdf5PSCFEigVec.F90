module O_PSCFEigVecHDF5

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

   ! Declare the eigenvector subgroup of the pscf_fid.
   integer(hid_t) :: eigenVectors_gid
   integer(hid_t) :: eigenVectorsSYBD_gid

   ! Declare the eigenvector dataspace.
   integer(hid_t) :: valeStates_dsid

   ! Declare the property list for the eigenvector data space.
   integer(hid_t) :: valeStates_plid

   ! Declare the eigenvector dataset array.  The number of datasets will depend
   !   on the number of kpoints and so it will vary.  Those datasets will be
   !   given ID numbers dynamically.  (Index1=real,imag; Index2=1..kpoints;
   !   Index3=1..spin)
   integer(hid_t), allocatable, dimension (:,:,:) :: eigenVectors_did
   integer(hid_t), allocatable, dimension (:,:,:) :: eigenVectorsSYBD_did

   ! Declare the attribute IDs that will be used to track completion of the
   !   eigenvector calculations. We need one for each kpoint and spin (because
   !   the real and imaginary will always be completed at the same time).
   integer(hid_t), allocatable, dimension (:,:) :: eigenVectors_aid
   integer(hid_t), allocatable, dimension (:,:) :: eigenVectorsSYBD_aid

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initPSCFEigVecHDF5 (pscf_fid,attribInt_dsid,attribIntDims,numStates)

   ! Import any necessary definition modules.
   use HDF5
   use O_Kinds

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints, numPathKP
   use O_AtomicSites, only: valeDim
   use O_Potential, only: spin

   ! Define the passed parameters.
   integer(hid_t) :: pscf_fid
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

   ! Create the eigenVectors group in the pscf_fid.
   call h5gcreate_f (pscf_fid,"eigenVectors",eigenVectors_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create eigenvectors group'
   call h5gcreate_f (pscf_fid,"eigenVectorsSYBD",eigenVectorsSYBD_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create eigenvectors SYBD group'

   ! Create the dataspace that will be used for the energy eigen vectors.
   call h5screate_simple_f(2,valeStates,valeStates_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create valeStates dataspace'

   ! Create the property list first.  Then set the properties one at a time.
   call h5pcreate_f(H5P_DATASET_CREATE_F,valeStates_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create valeStates plid'
   call h5pset_layout_f(valeStates_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set valeStates plid layout as chunked'
   call h5pset_chunk_f(valeStates_plid,2,valeStatesChunk,hdferr)
   if (hdferr /= 0) stop 'Failed to set valeStates plid chunk size'
!   call h5pset_shuffle_f(valeStates_plid,hdferr)
   call h5pset_deflate_f   (valeStates_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set valeStates for deflation'

   ! Allocate space to hold the IDs for the datasets in the eigenvectors group.
#ifndef GAMMA 
   allocate (eigenVectors_did(2,numKPoints,spin)) ! Complex
   allocate (eigenVectors_did(2,numPathKP,spin)) ! Complex
#else
   allocate (eigenVectors_did(1,numKPoints,spin)) ! Real
   allocate (eigenVectors_did(1,numPathKP,spin)) ! Real (unlikely to be used)
#endif

   ! Create the datasets for the eigenvectors.
   do i = 1, spin
      do j = 1, numKPoints
#ifndef GAMMA
         write (currentName,fmt="(a4,i7.7,i7.7)") "real",j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectors_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeStates_dsid,eigenVectors_did(1,j,i),hdferr,&
               & valeStates_plid)
         if (hdferr /= 0) stop 'Failed to create eigenVectors real did'
         write (currentName,fmt="(a4,i7.7,i7.7)") "imag",j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectors_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeStates_dsid,eigenVectors_did(2,j,i),hdferr,&
               & valeStates_plid)
         if (hdferr /= 0) stop 'Failed to create eigenVectors imaginary did'
#else
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectors_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeStates_dsid,eigenVectors_did(1,j,i),hdferr,&
               & valeStates_plid)
         if (hdferr /= 0) stop 'Failed to create eigenVectors real did'
#endif
      enddo
      do j = 1, numPathKP
#ifndef GAMMA
         write (currentName,fmt="(a4,i7.7,i7.7)") "real",j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectorsSYBD_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeStates_dsid,eigenVectorsSYBD_did(1,j,i),hdferr,&
               & valeStates_plid)
         if (hdferr /= 0) stop 'Failed to create eigenVectors SYBD real did'
         write (currentName,fmt="(a4,i7.7,i7.7)") "imag",j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectorsSYBD_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeStates_dsid,eigenVectorsSYBD_did(2,j,i),hdferr,&
               & valeStates_plid)
         if (hdferr /= 0) stop 'Failed to create eigenVectors SYBD imag. did'
#else
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectorsSYBD_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeStates_dsid,eigenVectorsSYBD_did(1,j,i),hdferr,&
               & valeStates_plid)
         if (hdferr /= 0) stop 'Failed to create eigenVectors SYBD real did'
#endif
      enddo
   enddo

   ! Allocate space to hold the attributes for tracking completion.
   allocate (eigenVectors_aid(numKPoints,spin))
   allocate (eigenVectors_aid(numPathKP,spin))

   ! Create the eigenvector tracking attributes.
   do i = 1, spin
      do j = 1, numKPoints
         call h5acreate_f (eigenVectors_did(1,j,i),"status",&
               & H5T_NATIVE_INTEGER,attribInt_dsid,eigenVectors_aid(j,i),&
               & hdferr)
            if (hdferr /= 0) stop 'Failed to create eigenVectors_aid.'
      enddo
      do j = 1, numPathKP
         call h5acreate_f (eigenVectorsSYBD_did(1,j,i),"status",&
               & H5T_NATIVE_INTEGER,attribInt_dsid,eigenVectorsSYBD_aid(j,i),&
               & hdferr)
            if (hdferr /= 0) stop 'Failed to create eigenVectorsSYBD_aid.'
      enddo
   enddo

   ! Init eigenvector tracking attributes to uncomputed (status=0) state.
   do i = 1, spin
      do j = 1, numKPoints
         call h5awrite_f (eigenVectors_aid(j,i),H5T_NATIVE_INTEGER,0,&
               & attribIntDims,hdferr)
            if (hdferr /= 0) stop 'Failed to init eigenVectors_aid.'
      enddo
      do j = 1, numPathKP
         call h5awrite_f (eigenVectorsSYBD_aid(j,i),H5T_NATIVE_INTEGER,0,&
               & attribIntDims,hdferr)
            if (hdferr /= 0) stop 'Failed to init eigenVectorsSYBD_aid.'
      enddo
   enddo

end subroutine initPSCFEigVecHDF5 


subroutine accessPSCFEigVecHDF5 (pscf_fid,attribInt_dsid,attribIntDims,numStates)

   ! Import any necessary definition modules.
   use HDF5
   use O_Kinds

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints, numPathKP
   use O_AtomicSites, only: valeDim
   use O_Potential, only: spin

   ! Define the passed parameters.
   integer(hid_t) :: pscf_fid
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

   ! Open the eigenVectors group in the pscf_fid.
   call h5gopen_f (pscf_fid,"/eigenVectors",eigenVectors_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open eigenvectors group.'
   call h5gopen_f (pscf_fid,"/eigenVectorsSYBD",eigenVectorsSYBD_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open eigenvectorsSYBD group.'

   ! Open the datasets for the eigenvectors.
   do i = 1, spin
      do j = 1, numKPoints
#ifndef GAMMA
         write (currentName,fmt="(a4,i7.7,i7.7)") "real",j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenVectors_gid,currentName,&
               & eigenVectors_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVectors real did.'
         write (currentName,fmt="(a4,i7.7,i7.7)") "imag",j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenVectors_gid,currentName,&
               & eigenVectors_did(2,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVectors imaginary did.'
#else
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenVectors_gid,currentName,&
               & eigenVectors_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVectors real did.'
#endif
      enddo
      do j = 1, numPathKP
#ifndef GAMMA
         write (currentName,fmt="(a4,i7.7,i7.7)") "real",j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenVectorsSYBD_gid,currentName,&
               & eigenVectorsSYBD_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVectorsSYBD real did.'
         write (currentName,fmt="(a4,i7.7,i7.7)") "imag",j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenVectorsSYBD_gid,currentName,&
               & eigenVectorsSYBD_did(2,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVectorsSYBD imaginary did.'
#else
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenVectorsSYBD_gid,currentName,&
               & eigenVectorsSYBD_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVectorsSYBD real did.'
#endif
      enddo
   enddo

   ! Obtain the property list for the energy eigen vectors. (Used for both
   !   SYBD and regular KPoint meshes.)
   call h5dget_create_plist_f(eigenVectors_did(1,1,1),valeStates_plid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to obtain valeStates property list.'

   ! Obtain the dataspace for the energy eigen vectors. (Used for both SYBD
   !   and regular KPoint meshes.)
   call h5dget_space_f(eigenVectors_did(1,1,1),valeStates_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain valeStates dataspace.'

   ! Allocate space to hold the attributes for tracking completion.
   allocate (eigenVectors_aid(numKPoints,spin))
   allocate (eigenVectors_aid(numPathKP,spin))

   ! Open the eigenvector tracking attributes.
   do i = 1, spin
      do j = 1, numKPoints
         call h5aopen_f (eigenVectors_did(1,j,i),"status",&
               & eigenVectors_aid(j,i),hdferr)
            if (hdferr /= 0) stop 'Failed to open eigenVectors_aid.'
      enddo
      do j = 1, numPathKP
         call h5aopen_f (eigenVectorsSYBD_did(1,j,i),"status",&
               & eigenVectorsSYBD_aid(j,i),hdferr)
            if (hdferr /= 0) stop 'Failed to open eigenVectorsSYBD_aid.'
      enddo
   enddo

end subroutine accessPSCFEigVecHDF5


subroutine closePSCFEigVecHDF5

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints, numPathKP
   use O_Potential, only: spin

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Declare local loop variables and error control.
   integer :: i,j
   integer :: hdferr

   ! Close the eigenvector dataspace.
   call h5sclose_f (valeStates_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close valeStates_dsid.'

   ! Close the eigenvector datasets next.
   do i = 1, spin
      do j = 1, numKPoints
#ifndef GAMMA
         call h5dclose_f (eigenVectors_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to close eigenVectors_did real'
         call h5dclose_f (eigenVectors_did(2,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to close eigenVectors_did imaginary'
#else
         call h5dclose_f (eigenVectors_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to close eigenVectors_did real'
#endif
      enddo
      do j = 1, numPathKP
#ifndef GAMMA
         call h5dclose_f (eigenVectorsSYBD_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to close eigenVectors_did SYBD real'
         call h5dclose_f (eigenVectorsSYBD_did(2,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to close eigenVectors_did SYBD imag.'
#else
         call h5dclose_f (eigenVectorsSYBD_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to close eigenVectorsSYBD_did real'
#endif
      enddo
   enddo

   ! Attributes are closed when the data is written.

   ! Close the eigenvector property list.
   call h5pclose_f (valeStates_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close valeStates_plid.'

   ! Close the eigenvector group.
   call h5gclose_f (eigenVectors_gid,hdferr)
   call h5gclose_f (eigenVectorsSYBD_gid,hdferr)

   ! Deallocate unnecessary array.
   deallocate (eigenVectors_did)
   deallocate (eigenVectorsSYBD_did)

end subroutine closePSCFEigVecHDF5

end module O_PSCFEigVecHDF5
