module O_MainEVecHDF5

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

   ! Declare the eigenvector subgroup of the main_fid.
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

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initMainEVecHDF5 (main_fid,numStates)

   ! Import any necessary definition modules.
   use HDF5
   use O_Kinds

   ! Import necessary object modules.
   use O_KPoints, only: numKPoints
   use O_AtomicSites, only: valeDim
   use O_Potential, only: spin

   ! Define the passed parameters.
   integer(hid_t) :: main_fid
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
      valeStatesChunk(2) = int(sqrt(real(250000000/(valeStates(1) * &
            & valeStates(2)),double)) * valeStates(2))
      valeStatesChunk(1) = int(250000000 / valeStatesChunk(2))
   else
      valeStatesChunk(1) = valeStates(1)
      valeStatesChunk(2) = valeStates(2)
   endif

   ! Create the eigenVectors group in the main_fid.
   call h5gcreate_f (main_fid,"eigenVectors",eigenVectors_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create eigenvectors group'

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
   allocate (eigenVectors_did(2,numKPoints,spin))
#else
   allocate (eigenVectors_did(1,numKPoints,spin))
#endif

   ! Create the datasets for the eigenvectors.
   do i = 1, numKPoints
      do j = 1, spin
#ifndef GAMMA
         write (currentName,fmt="(a4,i7.7,i7.7)") "real",i,j
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectors_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeStates_dsid,eigenVectors_did(1,i,j),hdferr,&
               & valeStates_plid)
         if (hdferr /= 0) stop 'Failed to create eigenVectors real did'
         write (currentName,fmt="(a4,i7.7,i7.7)") "imag",i,j
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectors_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeStates_dsid,eigenVectors_did(2,i,j),hdferr,&
               & valeStates_plid)
         if (hdferr /= 0) stop 'Failed to create eigenVectors imaginary did'
#else
         write (currentName,fmt="(i7.7,i7.7)") i,j
         currentName = trim (currentName)
         call h5dcreate_f(eigenVectors_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeStates_dsid,eigenVectors_did(1,i,j),hdferr,&
               & valeStates_plid)
         if (hdferr /= 0) stop 'Failed to create eigenVectors real did'
#endif
      enddo
   enddo

end subroutine initMainEVecHDF5 


subroutine closeMainEVecHDF5

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
   enddo

   ! Close the eigenvector property list.
   call h5pclose_f (valeStates_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close valeStates_plid.'

   ! Close the eigenvector group.
   call h5gclose_f (eigenVectors_gid,hdferr)

   ! Deallocate unnecessary array.
   deallocate (eigenVectors_did)

end subroutine closeMainEVecHDF5


end module O_MainEVecHDF5
