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
   integer(hsize_t), dimension (2) :: valeStatesPSCF

   ! Declare array that holds the dimensions of the chunk.
   integer(hsize_t), dimension (2) :: valeStatesPSCFChunk

   ! Declare the eigenvector subgroup of the pscf_fid.
   integer(hid_t) :: eigenVectorsPSCF_gid
   integer(hid_t) :: eigenVectorsSYBD_PSCF_gid

   ! Declare the eigenvector dataspace.
   integer(hid_t) :: valeStatesPSCF_dsid

   ! Declare the property list for the eigenvector data space.
   integer(hid_t) :: valeStatesPSCF_plid

   ! Declare the eigenvector dataset array.  The number of datasets will depend
   !   on the number of kpoints and so it will vary.  Those datasets will be
   !   given ID numbers dynamically.  (Index1=real,imag; Index2=1..kpoints;
   !   Index3=1..spin)
   integer(hid_t), allocatable, dimension (:,:,:) :: eigenVectorsPSCF_did
   integer(hid_t), allocatable, dimension (:,:,:) :: eigenVectorsSYBD_PSCF_did

   ! Declare the attribute IDs that will be used to track completion of the
   !   eigenvector calculations. We need one for each kpoint and spin (because
   !   the real and imaginary will always be completed at the same time).
   integer(hid_t), allocatable, dimension (:,:) :: eigenVectorsPSCF_aid
   integer(hid_t), allocatable, dimension (:,:) :: eigenVectorsSYBD_PSCF_aid

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initPSCFEigVecHDF5 (pscf_fid,attribIntPSCF_dsid,&
      & attribIntDimsPSCF,numStates)

   ! Import any necessary definition modules.
   use HDF5
   use O_Kinds

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints_HDF5, numPathKP_HDF5
   use O_AtomicSites, only: valeDim
   use O_Potential, only: spin

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Define the passed parameters.
   integer(hid_t) :: pscf_fid
   integer(hid_t) :: attribIntPSCF_dsid
   integer(hsize_t), dimension (1) :: attribIntDimsPSCF
   integer, intent(in) :: numStates

   ! Define local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize data structure dimensions.
   valeStatesPSCF(1)    = valeDim
   valeStatesPSCF(2)    = numStates

   ! Check that the chunk size is not too large (the assumption here is that
   !   the number being stored are 8 byte reals and that we should not go over
   !   2 billion bytes. Note that if x*y = >250M and we want a*b = 250M then
   !   the additional requirement x/y = a/b leads to b = sqrt(250M/>250M)*y.
   !   Thus a = 250M/b.
   if (valeStatesPSCF(1) * valeStatesPSCF(2) > 250000000) then
      valeStatesPSCFChunk(2) = int(sqrt(real(250000000,double) / &
            & real(valeStatesPSCF(1) * valeStatesPSCF(2),double)) * &
            & valeStatesPSCF(2))
      valeStatesPSCFChunk(1) = int(250000000 / valeStatesPSCFChunk(2))
   else
      valeStatesPSCFChunk(1) = valeStatesPSCF(1)
      valeStatesPSCFChunk(2) = valeStatesPSCF(2)
   endif

   ! Create the eigenVectors group in the pscf_fid.
   call h5gcreate_f (pscf_fid,"eigenVectors",eigenVectorsPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create eigenvectorsPSCF group'
   call h5gcreate_f (pscf_fid,"eigenVectorsSYBD",eigenVectorsSYBD_PSCF_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to create eigenvectorsSYBD_PSCF group'

   ! Create the dataspace that will be used for the energy eigen vectors.
   call h5screate_simple_f(2,valeStatesPSCF,valeStatesPSCF_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create valeStatesPSCF dataspace'

   ! Create the property list first.  Then set the properties one at a time.
   call h5pcreate_f(H5P_DATASET_CREATE_F,valeStatesPSCF_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create valeStatesPSCF plid'
   call h5pset_layout_f(valeStatesPSCF_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set valeStatesPSCF plid layout as chunked'
   call h5pset_chunk_f(valeStatesPSCF_plid,2,valeStatesPSCFChunk,hdferr)
   if (hdferr /= 0) stop 'Failed to set valeStatesPSCF plid chunk size'
!   call h5pset_shuffle_f(valeStatesPSCF_plid,hdferr)
   call h5pset_deflate_f   (valeStatesPSCF_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set valeStatesPSCF for deflation'

   ! Allocate space to hold IDs for the datasets in the eigenvectors group.
#ifndef GAMMA 
   allocate (eigenVectorsPSCF_did(2,numKPoints_HDF5,spin)) ! Complex
   allocate (eigenVectorsSYBD_PSCF_did(2,numPathKP_HDF5,spin)) ! Complex
#else
   allocate (eigenVectorsPSCF_did(1,numKPoints_HDF5,spin)) ! Real
   allocate (eigenVectorsSYBD_PSCF_did(1,numPathKP_HDF5,spin)) ! Real, useless?
#endif

   ! Create the datasets for the eigenvectors.
   do i = 1, spin
      do j = 1, numKPoints_HDF5
#ifndef GAMMA
         write (currentName,fmt="(a4,i7.7,i7.7)") "real",j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectorsPSCF_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeStatesPSCF_dsid,eigenVectorsPSCF_did(1,j,i),hdferr,&
               & valeStatesPSCF_plid)
         if (hdferr /= 0) stop 'Failed to create eigenVectors real did PSCF'
         write (currentName,fmt="(a4,i7.7,i7.7)") "imag",j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectorsPSCF_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeStatesPSCF_dsid,eigenVectorsPSCF_did(2,j,i),hdferr,&
               & valeStatesPSCF_plid)
         if (hdferr /= 0) stop 'Failed to create eigenVectors imag did PSCF'
#else
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectorsPSCF_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeStatesPSCF_dsid,eigenVectorsPSCF_did(1,j,i),hdferr,&
               & valeStatesPSCF_plid)
         if (hdferr /= 0) stop 'Failed to create eigenVectors real did PSCF'
#endif
      enddo

      ! Repeat for the SYBD eigenvectors.
      do j = 1, numPathKP_HDF5
#ifndef GAMMA
         write (currentName,fmt="(a4,i7.7,i7.7)") "real",j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectorsSYBD_PSCF_gid,currentName,&
               & H5T_NATIVE_DOUBLE,valeStatesPSCF_dsid,&
               & eigenVectorsSYBD_PSCF_did(1,j,i),hdferr,valeStatesPSCF_plid)
         if (hdferr /= 0) stop 'Failed create eigenVectors SYBD real did PSCF'
         write (currentName,fmt="(a4,i7.7,i7.7)") "imag",j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectorsSYBD_PSCF_gid,currentName,&
               & H5T_NATIVE_DOUBLE,valeStatesPSCF_dsid,&
               & eigenVectorsSYBD_PSCF_did(2,j,i),hdferr,valeStatesPSCF_plid)
         if (hdferr /= 0) stop 'Failed create eigenVectors SYBD imag. did PSCF'
#else
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectorsSYBD_PSCF_gid,currentName,&
               & H5T_NATIVE_DOUBLE,valeStatesPSCF_dsid,&
               & eigenVectorsSYBD_PSCF_did(1,j,i),hdferr,valeStatesPSCF_plid)
         if (hdferr /= 0) stop 'Failed to create eigenVectors SYBD real did PSCF'
#endif
      enddo
   enddo

   ! Allocate space to hold the attributes for tracking completion.
   allocate (eigenVectorsPSCF_aid(numKPoints_HDF5,spin))
   allocate (eigenVectorsSYBD_PSCF_aid(numPathKP_HDF5,spin))

   ! Create the eigenvector tracking attributes.
   do i = 1, spin
      do j = 1, numKPoints_HDF5
         call h5acreate_f (eigenVectorsPSCF_did(1,j,i),"status",&
               & H5T_NATIVE_INTEGER,attribIntPSCF_dsid,&
               & eigenVectorsPSCF_aid(j,i),hdferr)
            if (hdferr /= 0) stop 'Failed to create eigenVectorsPSCF_aid.'
      enddo
      do j = 1, numPathKP_HDF5
         call h5acreate_f (eigenVectorsSYBD_PSCF_did(1,j,i),"status",&
               & H5T_NATIVE_INTEGER,attribIntPSCF_dsid,&
               & eigenVectorsSYBD_PSCF_aid(j,i),hdferr)
            if (hdferr /= 0) stop 'Failed to create eigenVectorsSYBD_PSCF_aid.'
      enddo
   enddo

   ! Init eigenvector tracking attributes to uncomputed (status=0) state.
   do i = 1, spin
      do j = 1, numKPoints_HDF5
         call h5awrite_f (eigenVectorsPSCF_aid(j,i),H5T_NATIVE_INTEGER,0,&
               & attribIntDimsPSCF,hdferr)
            if (hdferr /= 0) stop 'Failed to init eigenVectorsPSCF_aid.'
      enddo
      do j = 1, numPathKP_HDF5
         call h5awrite_f (eigenVectorsSYBD_PSCF_aid(j,i),H5T_NATIVE_INTEGER,0,&
               & attribIntDimsPSCF,hdferr)
            if (hdferr /= 0) stop 'Failed to init eigenVectorsSYBD_PSCF_aid.'
      enddo
   enddo

end subroutine initPSCFEigVecHDF5 


subroutine accessPSCFEigVecHDF5 (pscf_fid,attribIntPSCF_dsid,&
      & attribIntDimsPSCF,numStates)

   ! Import any necessary definition modules.
   use HDF5
   use O_Kinds

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints_HDF5, numPathKP_HDF5
   use O_AtomicSites, only: valeDim
   use O_Potential, only: spin

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Define the passed parameters.
   integer(hid_t) :: pscf_fid
   integer(hid_t) :: attribIntPSCF_dsid
   integer(hsize_t), dimension (1) :: attribIntDimsPSCF
   integer, intent(in) :: numStates

   ! Define local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize data structure dimensions.
   valeStatesPSCF(1)    = valeDim
   valeStatesPSCF(2)    = numStates

   ! Check that the chunk size is not too large (the assumption here is that
   !   the number being stored are 8 byte reals and that we should not go over
   !   2 billion bytes. Note that if x*y = >250M and we want a*b = 250M then
   !   the additional requirement x/y = a/b leads to b = sqrt(250M/>250M)*y.
   !   Thus a = 250M/b.
   if (valeStatesPSCF(1) * valeStatesPSCF(2) > 250000000) then
      valeStatesPSCFChunk(2) = int(sqrt(real(250000000,double) / &
            & real(valeStatesPSCF(1) * valeStatesPSCF(2),double)) * &
            & valeStatesPSCF(2))
      valeStatesPSCFChunk(1) = int(250000000 / valeStatesPSCFChunk(2))
   else
      valeStatesPSCFChunk(1) = valeStatesPSCF(1)
      valeStatesPSCFChunk(2) = valeStatesPSCF(2)
   endif

   ! Open the eigenVectors group in the pscf_fid.
   call h5gopen_f (pscf_fid,"/eigenVectors",eigenVectorsPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open eigenvectorsPSCF group.'
   call h5gopen_f (pscf_fid,"/eigenVectorsSYBD",eigenVectorsSYBD_PSCF_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to open eigenvectorsSYBD_PSCF group.'

   ! Allocate space to hold IDs for the datasets in the eigenvectors group.
#ifndef GAMMA 
   allocate (eigenVectorsPSCF_did(2,numKPoints_HDF5,spin)) ! Complex
   allocate (eigenVectorsSYBD_PSCF_did(2,numPathKP_HDF5,spin)) ! Complex
#else
   allocate (eigenVectorsPSCF_did(1,numKPoints_HDF5,spin)) ! Real
   allocate (eigenVectorsSYBD_PSCF_did(1,numPathKP_HDF5,spin)) ! Real, useless?
#endif

   ! Open the datasets for the eigenvectors.
   do i = 1, spin
      do j = 1, numKPoints_HDF5
#ifndef GAMMA
         write (currentName,fmt="(a4,i7.7,i7.7)") "real",j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenVectorsPSCF_gid,currentName,&
               & eigenVectorsPSCF_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVectorsPSCF real did.'
         write (currentName,fmt="(a4,i7.7,i7.7)") "imag",j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenVectorsPSCF_gid,currentName,&
               & eigenVectorsPSCF_did(2,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVectorsPSCF imag did.'
#else
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenVectorsPSCF_gid,currentName,&
               & eigenVectorsPSCF_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVectorsPSCF real did.'
#endif
      enddo

      ! Repeat for the SYBD eigenvectors.
      do j = 1, numPathKP_HDF5
#ifndef GAMMA
         write (currentName,fmt="(a4,i7.7,i7.7)") "real",j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenVectorsSYBD_PSCF_gid,currentName,&
               & eigenVectorsSYBD_PSCF_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVectorsSYBD_PSCF real did.'
         write (currentName,fmt="(a4,i7.7,i7.7)") "imag",j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenVectorsSYBD_PSCF_gid,currentName,&
               & eigenVectorsSYBD_PSCF_did(2,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVectorsSYBD_PSCF imag did.'
#else
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenVectorsSYBD_PSCF_gid,currentName,&
               & eigenVectorsSYBD_PSCF_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVectorsSYBD_PSCF real did.'
#endif
      enddo
   enddo

   ! Obtain the property list for the energy eigen vectors. (Used for both
   !   SYBD and regular KPoint meshes.)
   call h5dget_create_plist_f(eigenVectorsPSCF_did(1,1,1),valeStatesPSCF_plid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to obtain valeStatesPSCF property list.'

   ! Obtain the dataspace for the energy eigen vectors. (Used for both SYBD
   !   and regular KPoint meshes.)
   call h5dget_space_f(eigenVectorsPSCF_did(1,1,1),valeStatesPSCF_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain valeStatesPSCF dataspace.'

   ! Allocate space to hold the attributes for tracking completion.
   allocate (eigenVectorsPSCF_aid(numKPoints_HDF5,spin))
   allocate (eigenVectorsSYBD_PSCF_aid(numPathKP_HDF5,spin))

   ! Open the eigenvector tracking attributes.
   do i = 1, spin
      do j = 1, numKPoints_HDF5
         call h5aopen_f (eigenVectorsPSCF_did(1,j,i),"status",&
               & eigenVectorsPSCF_aid(j,i),hdferr)
            if (hdferr /= 0) stop 'Failed to open eigenVectorsPSCF_aid.'
      enddo
      do j = 1, numPathKP_HDF5
         call h5aopen_f (eigenVectorsSYBD_PSCF_did(1,j,i),"status",&
               & eigenVectorsSYBD_PSCF_aid(j,i),hdferr)
            if (hdferr /= 0) stop 'Failed to open eigenVectorsSYBD_PSCF_aid.'
      enddo
   enddo

end subroutine accessPSCFEigVecHDF5


subroutine closePSCFEigVecHDF5

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints_HDF5, numPathKP_HDF5
   use O_Potential, only: spin

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Declare local loop variables and error control.
   integer :: i,j
   integer :: hdferr

   ! Close the eigenvector dataspace.
   call h5sclose_f (valeStatesPSCF_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close valeStatesPSCF_dsid.'

   ! Close the eigenvector datasets next.
   do i = 1, spin
      do j = 1, numKPoints_HDF5
#ifndef GAMMA
         call h5dclose_f (eigenVectorsPSCF_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to close eigenVectorsPSCF_did real'
         call h5dclose_f (eigenVectorsPSCF_did(2,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to close eigenVectorsPSCF_did imaginary'
#else
         call h5dclose_f (eigenVectorsPSCF_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to close eigenVectorsPSCF_did real'
#endif
      enddo
      do j = 1, numPathKP_HDF5
#ifndef GAMMA
         call h5dclose_f (eigenVectorsSYBD_PSCF_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to close eigenVectorsPSCF_did SYBD real'
         call h5dclose_f (eigenVectorsSYBD_PSCF_did(2,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to close eigenVectorsPSCF_did SYBD imag.'
#else
         call h5dclose_f (eigenVectorsSYBD_PSCF_did(1,j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to close eigenVectorsSYBD_PSCF_did real'
#endif
      enddo
   enddo

   ! Attributes are closed when the data is written.

   ! Close the eigenvector property list.
   call h5pclose_f (valeStatesPSCF_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close valeStatesPSCF_plid.'

   ! Close the eigenvector group.
   call h5gclose_f (eigenVectorsPSCF_gid,hdferr)
   call h5gclose_f (eigenVectorsSYBD_PSCF_gid,hdferr)

   ! Deallocate unnecessary array.
   deallocate (eigenVectorsPSCF_did)
   deallocate (eigenVectorsSYBD_PSCF_did)

end subroutine closePSCFEigVecHDF5

end module O_PSCFEigVecHDF5
