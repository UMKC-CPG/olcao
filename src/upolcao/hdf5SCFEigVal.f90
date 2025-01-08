module O_SCFEigValHDF5

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

   ! Declare the eigenvalue subgroup of the scf_fid.
   integer(hid_t) :: eigenValues_gid

   ! Declare the eigenvalue dataspace.
   integer(hid_t) :: states_dsid

   ! Declare the property list for the eigenvalue data space.
   integer(hid_t) :: states_plid

   ! Declare the eigenvalue dataset array.  The number of datasets will depend
   !   on the number of kpoints and so it will vary.  Those datasets will be
   !   given ID numbers dynamically.  (Index1=1..kpoints; Index2=1..spin)
   integer(hid_t), allocatable, dimension (:,:)   :: eigenValues_did

   ! No status attributes are used for the eigenvalues. If the eigenvectors
   !   are completed, then it is assumed that the eigenvalues are also
   !   complete.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

      
subroutine initSCFEigValHDF5 (scf_fid,numStates)

   ! Import any necessary definition modules.
   use HDF5
   use O_MPI

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints
   use O_Potential, only: spin

   ! Define the passed parameters.
   integer(hid_t) :: scf_fid
   integer, intent(in) :: numStates

   ! Define local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize data structure dimensions.
   states(1) = numStates

   ! Allocate space to hold the IDs for the datasets in the eigenvalues group.
   allocate (eigenValues_did(numKPoints,spin))

   ! Only process 0 opens the HDF5 structure.
   if (mpiRank /= 0) return

   ! Create the eigenValues group in the scf_fid.
   call h5gcreate_f (scf_fid,"eigenValues",eigenValues_gid,hdferr)

   ! Create the dataspace that will be used for the energy eigen values.
   call h5screate_simple_f(1,states,states_dsid,hdferr)

   ! Create the property list first.  Then set the properties one at a time.
   call h5pcreate_f(H5P_DATASET_CREATE_F,states_plid,hdferr)
   call h5pset_layout_f(states_plid,H5D_CHUNKED_F,hdferr)
   call h5pset_chunk_f(states_plid,1,states,hdferr)
!   call h5pset_shuffle_f(states_plid,hdferr)
   call h5pset_deflate_f   (states_plid,1,hdferr)

   ! Create the datasets for the eigen values.
   do i = 1, spin
      do j = 1, numKPoints
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenValues_gid,currentName,H5T_NATIVE_DOUBLE,&
               & states_dsid,eigenValues_did(j,i),hdferr,states_plid)
      enddo
   enddo


end subroutine initSCFEigValHDF5


subroutine accessSCFEigValHDF5 (scf_fid,numStates)

   ! Import any necessary definition modules.
   use HDF5
   use O_MPI

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints
   use O_Potential, only: spin

   ! Define the passed parameters.
   integer(hid_t) :: scf_fid
   integer, intent(in) :: numStates

   ! Define local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize data structure dimensions.
   states(1) = numStates

   ! Allocate space to hold the IDs for the datasets in the eigenvalues group.
   allocate (eigenValues_did(numKPoints,spin))

   ! Only process 0 opens the HDF5 structure.
   if (mpiRank /= 0) return

   ! Open the eigenValues group in the scf_fid.
   call h5gopen_f (scf_fid,"/eigenValues",eigenValues_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open eigenValues group.'

   ! Open the datasets for the eigen vectors.
   do i = 1, spin
      do j = 1, numKPoints
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenValues_gid,currentName,eigenValues_did(j,i),&
               & hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenValues dataset.'
      enddo
   enddo

   ! Checkpointing attributes are only used for the eigen vectors.

   ! Obtain the property list for the eigenvalues.
   call h5dget_create_plist_f(eigenValues_did(1,1),states_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain eigenvalues property list.'

   ! Obtain the dataspace for the eigenvalues.
   call h5dget_space_f(eigenValues_did(1,1),states_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain eigenvalues dataspace.'

end subroutine accessSCFEigValHDF5


subroutine closeSCFEigValHDF5

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

      ! Close the eigenvalue dataspace.
      call h5sclose_f (states_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to close states_dsid.'

      ! Close the eigenvalue datasets next.
      do i = 1, spin
         do j = 1, numKPoints
            call h5dclose_f (eigenValues_did(j,i),hdferr)
            if (hdferr /= 0) stop 'Failed to close eigenValues_did'
         enddo
      enddo

      ! Close the eigenvalue property list.
      call h5pclose_f (states_plid,hdferr)
      if (hdferr /= 0) stop 'Failed to close states_plid.'

      ! Close the eigenvalues group.
      call h5gclose_f (eigenValues_gid,hdferr)
   endif

   ! Dellocate space used to hold IDs for datasets in the eigenvalues group.
   deallocate (eigenValues_did)

end subroutine closeSCFEigValHDF5


end module O_SCFEigValHDF5
