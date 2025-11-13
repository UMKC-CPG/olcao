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
   integer(hid_t) :: eigenValuesPSCF_gid

   ! Declare the eigenvalue dataspace.
   integer(hid_t) :: statesPSCF_dsid

   ! Declare the property list for the eigenvalue data space.
   integer(hid_t) :: statesPSCF_plid

   ! Declare the eigenvalue dataset array.  The number of datasets will depend
   !   on the number of kpoints and so it will vary.  Those datasets will be
   !   given ID numbers dynamically.  (Index1=1..kpoints; Index2=1..spin)
   integer(hid_t), allocatable, dimension (:,:)   :: eigenValuesPSCF_did

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
   use O_KPoints, only: numKPoints
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
   call h5gcreate_f (pscf_fid,"eigenValues",eigenValuesPSCF_gid,hdferr)

   ! Create the dataspace that will be used for the energy eigen values.
   call h5screate_simple_f(1,states,statesPSCF_dsid,hdferr)

   ! Create the property list first.  Then set the properties one at a time.
   call h5pcreate_f(H5P_DATASET_CREATE_F,statesPSCF_plid,hdferr)
   call h5pset_layout_f(statesPSCF_plid,H5D_CHUNKED_F,hdferr)
   call h5pset_chunk_f(statesPSCF_plid,1,states,hdferr)
!   call h5pset_shuffle_f(statesPSCF_plid,hdferr)
   call h5pset_deflate_f   (statesPSCF_plid,1,hdferr)

   ! Allocate space to hold the IDs for the datasets in the eigenvalues group.
   allocate (eigenValuesPSCF_did(numKPoints,spin))

   ! Create the datasets for the eigen values.
   do i = 1, spin
      do j = 1, numKPoints
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dcreate_f(eigenValuesPSCF_gid,currentName,H5T_NATIVE_DOUBLE,&
               & statesPSCF_dsid,eigenValuesPSCF_did(j,i),hdferr,&
               & statesPSCF_plid)
      enddo
   enddo

end subroutine initPSCFEigValHDF5


subroutine accessPSCFEigValHDF5 (pscf_fid, numStates)

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints
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
   call h5gopen_f (pscf_fid,"eigenValues",eigenValuesPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open eigenValues group PSCF'

   ! Allocate space to hold the IDs for the datasets in the eigenvalues group.
   allocate (eigenValuesPSCF_did(numKPoints,spin))

   ! Open the datasets for the eigen vectors.
   do i = 1, spin
      do j = 1, numKPoints
         write (currentName,fmt="(i7.7,i7.7)") j,i
         currentName = trim (currentName)
         call h5dopen_f(eigenValuesPSCF_gid,currentName,&
               & eigenValuesPSCF_did(j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenValues dataset PSCF'
      enddo
   enddo

   ! Checkpointing attributes are only used for the eigen vectors.

   ! Obtain the property list for the eigenvalues.
   call h5dget_create_plist_f(eigenValuesPSCF_did(1,1),statesPSCF_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain eigenvalues property list PSCF'

   ! Obtain the dataspace for the eigenvalues.
   call h5dget_space_f(eigenValuesPSCF_did(1,1),statesPSCF_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain eigenvalues dataspace PSCF'

end subroutine accessPSCFEigValHDF5


subroutine closePSCFEigValHDF5

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
   call h5sclose_f (statesPSCF_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close statesPSCF_dsid PSCF'

   ! Close the eigenvalue datasets next.
   do i = 1, spin
      do j = 1, numKPoints
         call h5dclose_f (eigenValuesPSCF_did(j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to close eigenValuesPSCF_did PSCF'
      enddo
   enddo

   ! Dellocate space used to hold IDs for datasets in the eigenvalues group.
   deallocate (eigenValuesPSCF_did)

   ! Close the eigenvalue property list.
   call h5pclose_f (statesPSCF_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close statesPSCF_plid PSCF'

   ! Close the eigenvalues group.
   call h5gclose_f (eigenValuesPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close eigen values gid PSCF'

end subroutine closePSCFEigValHDF5

end module O_PSCFEigValHDF5
