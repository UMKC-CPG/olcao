module O_MainEValHDF5

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

   ! Declare the eigenvalue subgroup of the main_fid.
   integer(hid_t) :: eigenValues_gid

   ! Declare the eigenvalue dataspace.
   integer(hid_t) :: states_dsid

   ! Declare the property list for the eigenvalue data space.
   integer(hid_t) :: states_plid

   ! Declare the eigenvalue dataset array.  The number of datasets will depend
   !   on the number of kpoints and so it will vary.  Those datasets will be
   !   given ID numbers dynamically.  (Index1=1..kpoints; Index2=1..spin)
   integer(hid_t), allocatable, dimension (:,:)   :: eigenValues_did

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initMainEValHDF5 (main_fid,numStates)

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints
   use O_Potential, only: spin

   ! Define the passed parameters.
   integer(hid_t) :: main_fid
   integer, intent(in) :: numStates

   ! Define local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize data structure dimensions.
   states(1) = numStates

   ! Create the eigenValues group in the main_fid.
   call h5gcreate_f (main_fid,"eigenValues",eigenValues_gid,hdferr)

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

   ! Create the datasets for the eigen values.
   do i = 1, numKPoints
      do j = 1, spin
         write (currentName,fmt="(i7.7,i7.7)") i,j
         currentName = trim (currentName)
         call h5dcreate_f(eigenValues_gid,currentName,H5T_NATIVE_DOUBLE,&
               & states_dsid,eigenValues_did(i,j),hdferr,states_plid)
      enddo
   enddo


end subroutine initMainEValHDF5


subroutine closeMainEValHDF5

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints
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
   enddo

   ! Close the eigenvalue property list.
   call h5pclose_f (states_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close states_plid.'

   ! Close the eigenvalues group.
   call h5gclose_f (eigenValues_gid,hdferr)

end subroutine closeMainEValHDF5


end module O_MainEValHDF5
