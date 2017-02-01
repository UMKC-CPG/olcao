! A module for storing the potential function coefficients and charge density
!   coefficients.
module O_MainPotRhoHDF5

   ! Use the HDF5 module (for hsize_t and hid_t).
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Declare array to hold the dimension of a given coefficient dataset.
   integer(hsize_t), dimension (1) :: terms

   ! Declare the potential and various fitted charge density coeffs subgroups
   !   of the main_fid. Also make a group for the alphas.
   integer(hid_t) :: potCoeffs_gid
   integer(hid_t) :: totalRhoCoeffs_gid
   integer(hid_t) :: valeRhoCoeffs_gid
   integer(hid_t) :: spinDiffRhoCoeffs_gid
   integer(hid_t) :: alphas_gid

   ! Declare the potRhoCoeffs dataspace. The alphas will use the same dataspace.
   !   (They all have the same dataspace.)
   integer(hid_t) :: potRhoCoeffs_dsid

   ! Declare the property list for the potential and charge density
   !   coefficients data space. The alphas will share this property list too.
   !   (They all have the same property list.)
   integer(hid_t) :: potRhoCoeffs_plid

   ! Declare dataset arrays for the potential coefficients, the total fitted
   !   charge density, the valence fitted charge density, and the fitted spin
   !   difference (only if needed).  The number of datasets of each will depend
   !   on the number of iterations and so it will vary.  Those datasets will be
   !   given ID numbers dynamically.  Note that only the potCoeffs will have a
   !   spin up and spin down set. The charge density will have a total (up +
   !   down), a total valence (up + down), and a spin difference (up - down).
   !   There is only one alpha dataset because it is a constant for all
   !   coefficient sets and it does not change with the iterations.
   !   (Index1=1..iteration; Index2=1..spin)
   integer(hid_t), allocatable, dimension (:,:) :: potCoeffs_did
   integer(hid_t), allocatable, dimension (:)   :: totalRhoCoeffs_did
   integer(hid_t), allocatable, dimension (:)   :: valeRhoCoeffs_did
   integer(hid_t), allocatable, dimension (:)   :: spinDiffRhoCoeffs_did
   integer(hid_t) :: alphas_did

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initMainPotRhoHDF5 (main_fid)

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_Potential, only: spin, potDim, lastIteration

   ! Define the passed parameters.
   integer(hid_t) :: main_fid

   ! Define local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize data structure dimensions.
   terms(1) = potDim

   ! Create the potential coefficient and various fitted charge density
   !   coefficient groups in the main_fid.
   call h5gcreate_f (main_fid,"potCoeffs",potCoeffs_gid,hdferr)
   call h5gcreate_f (main_fid,"totalRhoCoeffs",totalRhoCoeffs_gid,hdferr)
   call h5gcreate_f (main_fid,"valeRhoCoeffs",valeRhoCoeffs_gid,hdferr)
   call h5gcreate_f (main_fid,"spinDiffRhoCoeffs",spinDiffRhoCoeffs_gid,hdferr)
   call h5gcreate_f (main_fid,"alphas",alphas_gid,hdferr)

   ! Create the dataspace that will be used for the potential and charge density
   !   coefficients.
   call h5screate_simple_f(1,terms,potRhoCoeffs_dsid,hdferr)

   ! Create the property list first.  Then set the properties one at a time.
   call h5pcreate_f(H5P_DATASET_CREATE_F,potRhoCoeffs_plid,hdferr)
   call h5pset_layout_f(potRhoCoeffs_plid,H5D_CHUNKED_F,hdferr)
   call h5pset_chunk_f(potRhoCoeffs_plid,1,terms,hdferr)
!   call h5pset_shuffle_f(states_plid,hdferr)
   call h5pset_deflate_f   (potRhoCoeffs_plid,1,hdferr)

   ! Allocate space to hold the IDs for the datasets in the potRhoCoeffs group.
   allocate (potCoeffs_did(lastIteration,spin))
   allocate (totalRhoCoeffs_did(lastIteration))
   allocate (valeRhoCoeffs_did(lastIteration))
   allocate (spinDiffRhoCoeffs_did(lastIteration))

   ! Create the datasets for the potential coefficients.
   do i = 1, spin
      do j = 1, lastIteration
         write (currentName,fmt="(i7.7,i7.7)") i,j
         currentName = trim (currentName)
         call h5dcreate_f(potCoeffs_gid,currentName,H5T_NATIVE_DOUBLE,&
               & potRhoCoeffs_dsid,potCoeffs_did(j,i),hdferr,&
               & potRhoCoeffs_plid)
      enddo
   enddo

   ! Create the datasets for the various charge density coeffs.
   do i = 1, lastIteration
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)
      call h5dcreate_f(totalRhoCoeffs_gid,currentName,H5T_NATIVE_DOUBLE,&
            & potRhoCoeffs_dsid,totalRhoCoeffs_did(i),hdferr,&
            & potRhoCoeffs_plid)
      call h5dcreate_f(valeRhoCoeffs_gid,currentName,H5T_NATIVE_DOUBLE,&
            & potRhoCoeffs_dsid,valeRhoCoeffs_did(i),hdferr,&
            & potRhoCoeffs_plid)
      call h5dcreate_f(spinDiffRhoCoeffs_gid,currentName,H5T_NATIVE_DOUBLE,&
            & potRhoCoeffs_dsid,spinDiffRhoCoeffs_did(i),hdferr,&
            & potRhoCoeffs_plid)
   enddo

   ! Create the dataset for the alpahs. (Recall that we use the potRhoCoeffs
   !   dataspace and property list.)
   write (currentName,fmt="(i7.7)") 1
   currentName = trim (currentName)
   call h5dcreate_f(alphas_gid,currentName,H5T_NATIVE_DOUBLE,&
         & potRhoCoeffs_dsid,alphas_did,hdferr,potRhoCoeffs_plid)

end subroutine initMainPotRhoHDF5


subroutine closeMainPotRhoHDF5

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_Potential, only: spin, lastIteration

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Declare local loop variables and error control.
   integer :: i,j,k
   integer :: hdferr

   ! Close the potRhoCoeffs dataspace.
   call h5sclose_f (potRhoCoeffs_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close potRhoCoeffs_dsid.'

   ! Close the potential coefficient datasets.
   do i = 1, spin
      do j = 1, lastIteration
         call h5dclose_f (potCoeffs_did(j,i),hdferr)
         if (hdferr /= 0) stop 'Failed to close potCoeffs_did'
      enddo
   enddo

   ! Close the charge density coefficient datasets.
   do i = 1, lastIteration
      call h5dclose_f (totalRhoCoeffs_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close totalRhoCoeffs_did'
      call h5dclose_f (valeRhoCoeffs_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close valeRhoCoeffs_did'
      call h5dclose_f (spinDiffRhoCoeffs_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close spinDiffRhoCoeffs_did'
   enddo

   ! Close the alphas dataset.
   call h5dclose_f (alphas_did,hdferr)
   if (hdferr /= 0) stop 'Failed to close alphas_did'


   ! Deallocate space used to hold dataset IDs for the potential and charge
   !   density coefficients.
   deallocate (potCoeffs_did)
   deallocate (totalRhoCoeffs_did)
   deallocate (valeRhoCoeffs_did)
   deallocate (spinDiffRhoCoeffs_did)


   ! Close the potRhoCoeffs property list.
   call h5pclose_f (potRhoCoeffs_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close potRhoCoeffs_plid.'

   ! Close the potential coefficients group.
   call h5gclose_f (potCoeffs_gid,hdferr)
   call h5gclose_f (totalRhoCoeffs_gid,hdferr)
   call h5gclose_f (valeRhoCoeffs_gid,hdferr)
   call h5gclose_f (spinDiffRhoCoeffs_gid,hdferr)
   call h5gclose_f (alphas_gid,hdferr)

end subroutine closeMainPotRhoHDF5


end module O_MainPotRhoHDF5
