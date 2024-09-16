module O_PSCFEigValHDF5

   ! Use the HDF5 module (for hsize_t and hid_t).
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Declare arrays that hold the dimensions of the eigenvalue dataset.
   integer(hsize_t), dimension (1) :: states

   ! Declare the eigenvalue subgroup of the pscf_fid.
   integer(hid_t) :: eigenValues_gid
   integer(hid_t) :: eigenValuesSYBD_gid

   ! Declare the eigenvalue dataspace.
   integer(hid_t) :: states_dsid

   ! Declare the property list for the eigenvalue data space.
   integer(hid_t) :: states_plid

   ! Declare the eigenvalue dataset array.  The number of datasets will depend
   !   on the number of kpoints and so it will vary.  Those datasets will be
   !   given ID numbers dynamically.  (Index1=1..kpoints; Index2=1..spin)
   integer(hid_t), allocatable, dimension (:,:)   :: eigenValues_did
   integer(hid_t), allocatable, dimension (:,:)   :: eigenValuesSYBD_did

   ! No status attributes are used for the eigenvalues. If the eigenvectors
   !   are completed, then it is assumed that the eigenvalues are also
   !   complete.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initPSCFEigValHDF5 (pscf_fid,numStates)

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints, numPathKP
   use O_Potential, only: spin

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer(hid_t) :: pscf_fid
   integer, intent(in) :: numStates

   ! Define local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize data structure dimensions.
   states(1) = numStates

   ! Create the eigenValues group in the pscf_fid.
   call h5gcreate_f (pscf_fid,"eigenValues",eigenValues_gid,hdferr)
   call h5gcreate_f (pscf_fid,"eigenValuesSYBD",eigenValuesSYBD_gid,hdferr)

   ! Create the dataspace that will be used for the energy eigen values.
   call h5screate_simple_f(1,states,states_dsid,hdferr)

   ! Create the property list first.  Then set the properties one at a time.
   call h5pcreate_f(H5P_DATASET_CREATE_F,states_plid,hdferr)
   call h5pset_layout_f(states_plid,H5D_CHUNKED_F,hdferr)
   call h5pset_chunk_f(states_plid,1,states,hdferr)
!   call h5pset_shuffle_f(states_plid,hdferr)
   call h5pset_deflate_f   (states_plid,1,hdferr)

   ! Allocate space to hold the IDs for the datasets in the eigenvalues group.
   allocate (eigenValues_did(numKPoints,spin))
   allocate (eigenValuesSYBD_did(numPathKP,spin))

   ! Create the datasets for the eigen values.
   do i = 1, spin
      do j = 1, numKPoints
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenValues_gid,currentName,H5T_NATIVE_DOUBLE,&
               & states_dsid,eigenValues_did(j,i),hdferr,states_plid)
      enddo
      do j = 1, numPathKP
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenValuesSYBD_gid,currentName,H5T_NATIVE_DOUBLE,&
               & states_dsid,eigenValuesSYBD_did(j,i),hdferr,states_plid)
      enddo
   enddo

end subroutine initPSCFEigValHDF5


subroutine accessPSCFEigValHDF5 (pscf_fid, numStates)

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints, numPathKP
   use O_Potential, only: spin

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   integer(hid_t) :: pscf_fid
   integer, intent(in) :: numStates

   ! Define local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize data structure dimensions.
   states(1) = numStates

   ! Open the eigenValues group in the scf_fid.
   call h5gopen_f (pscf_fid,"/eigenValues",eigenValues_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open eigenValues group.'
   call h5gopen_f (pscf_fid,"/eigenValuesSYBD",eigenValuesSYBD_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open eigenValuesSYBD group.'

   ! Allocate space to hold the IDs for the datasets in the eigenvalues group.
   allocate (eigenValues_did(numKPoints,spin))
   allocate (eigenValuesSYBD_did(numPathKP,spin))

   ! Open the datasets for the eigen vectors.
   do i = 1, spin
      do j = 1, numKPoints
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenValues_gid,currentName,eigenValues_did(j,i),&
               & hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenValues dataset.'
      enddo
      do j = 1, numPathKP
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenValuesSYBD_gid,currentName,&
               & eigenValuesSYBD_did(j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenValuesSYBD dataset.'
      enddo
   enddo

   ! Checkpointing attributes are only used for the eigen vectors.

   ! Obtain the property list for the eigenvalues. (Used for both SYBD
   !   and regular KPoint meshes.)
   call h5dget_create_plist_f(eigenValues_did(1,1),states_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain eigenvalues property list.'

   ! Obtain the dataspace for the eigenvalues. (Used for both SYBD and
   !   regular KPoint meshes.)
   call h5dget_space_f(eigenValues_did(1,1),states_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain eigenvalues dataspace.'

end subroutine accessPSCFEigValHDF5


subroutine closePSCFEigValHDF5

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

   ! Close the eigenvalue dataspace.
   call h5sclose_f (states_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close states_dsid.'

   ! Close the eigenvalue datasets next.
   do i = 1, spin
      do j = 1, numKPoints
         call h5dclose_f (eigenValues_did(j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to close eigenValues_did'
      enddo
      do j = 1, numPathKP
         call h5dclose_f (eigenValuesSYBD_did(j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to close eigenValuesSYBD_did'
      enddo
   enddo

   ! Dellocate space used to hold IDs for datasets in the eigenvalues group.
   deallocate (eigenValues_did)
   deallocate (eigenValuesSYBD_did)

   ! Close the eigenvalue property list.
   call h5pclose_f (states_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close states_plid.'

   ! Close the eigenvalues group.
   call h5gclose_f (eigenValues_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close eigen values gid.'
   call h5gclose_f (eigenValuesSYBD_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close eigen values SYBD gid.'

end subroutine closePSCFEigValHDF5

end module O_PSCFEigValHDF5
