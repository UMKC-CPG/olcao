module O_SetupIntegralsHDF5

   ! Import any necessary definition modules.
   use HDF5
   use O_Potential, only: rel
   use O_CommandLine, only: doDIMO

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define the main subgroup from setup_fid that holds the multicenter
   !   Gaussian integral results.
   integer(hid_t) :: atomIntgGroup_gid ! Atom and potential overlap group.

   ! Begin group, dataspace, dataset definitions for atomIntgGroup_gid.

   ! Define the group IDs under atomIntgGroup_gid
   integer(hid_t) :: atomOverlap_gid
   integer(hid_t) :: atomKEOverlap_gid
   integer(hid_t) :: atomMVOverlap_gid
   integer(hid_t) :: atomNucOverlap_gid
   integer(hid_t) :: atomDMOverlap_gid
   integer(hid_t) :: atomPotOverlap_gid

   ! Define the group IDs of the dynamically numbered subgroups of
   !   atomPotOverlap_gid.  (Number of kpoints)
   integer(hid_t), allocatable, dimension (:) :: atomPotKPointOL_gid
   integer(hid_t), allocatable, dimension (:) :: atomDMKPointOL_gid

   ! The dataspaces of each dataset in atomIntgGroup_gid are the same in
   !   all characteristics (type, dimension, etc.) and therefore can be
   !   given a static ID definition now.
   integer(hid_t) :: valeVale_dsid

   ! Each of the below datasets for atomIntgGroup will use the same property
   !   list.
   integer(hid_t) :: valeVale_plid

   ! Define array that holds the dimensions of the dataset.
   integer(hsize_t), dimension (2) :: atomDims

   ! Define array that holds the dimensions of the chunk.
   integer (hsize_t), dimension (2) :: atomDimsChunk

   ! The number of datasets under each of the atomIntgGroup groups will
   !   vary depending on the problem.  (# of kpoints, potential dimension).
   !   Those datasets will be given IDs dynamically.
   integer(hid_t), allocatable, dimension (:)   :: atomOverlap_did
   integer(hid_t), allocatable, dimension (:)   :: atomKEOverlap_did
   integer(hid_t), allocatable, dimension (:)   :: atomMVOverlap_did
   integer(hid_t), allocatable, dimension (:)   :: atomNucOverlap_did
   integer(hid_t), allocatable, dimension (:,:) :: atomDMOverlap_did
   integer(hid_t), allocatable, dimension (:,:) :: atomPotOverlap_did

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initSetupIntegralHDF5 (setup_fid)

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_KPoints     ! For numKPoints
   use O_Potential   ! For potDim
   use O_AtomicSites ! For valeDim

   ! Define the passed parameters.
   integer(hid_t) :: setup_fid

   ! Define local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize data structure dimensions.
#ifndef GAMMA
   atomDims(1) = 2 ! Real and imaginary are needed.
#else
   atomDims(1) = 1 ! Real needed only
#endif
   atomDims(2) = valeDim*(valeDim+1)/2 ! Linear storage of 1/2 matrix.

   ! Check that the chunk size is not too large. The assumption here is that
   !   the number being stored are 8 byte reals and that we should not go over
   !   2 billion bytes. Note that if x*y = >250M and we want a*b = 250M then
   !   the additional requirement x/y = a/b leads to b = sqrt(250M/>250M)*y.
   !   Thus a = 250M/b.
   if (atomDims(1) * atomDims(2) > 250000000) then
      atomDimsChunk(2) = int(250000000/atomDims(1))
      atomDimsChunk(1) = atomDims(1)
   else
      atomDimsChunk(1) = atomDims(1)
      atomDimsChunk(2) = atomDims(2)
   endif

   ! Create the Integral group within the setup HDF5 file.
   call h5gcreate_f (setup_fid,"/atomIntgGroup",atomIntgGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create atom intg group'

   ! Create the subgroups within the atomIntgGroup.
   call h5gcreate_f (atomIntgGroup_gid,"atomOverlap",atomOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create atom overlap group'
   call h5gcreate_f (atomIntgGroup_gid,"atomKEOverlap",atomKEOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create kinetic energy overlap group'
   if (rel == 1) then
      call h5gcreate_f (atomIntgGroup_gid,"atomMVOverlap",atomMVOverlap_gid,&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to create mass velocity overlap group'
   endif
   call h5gcreate_f (atomIntgGroup_gid,"atomNucOverlap",atomNucOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to create atom nuclear overlap group'
   if (doDIMO == 1) then
      call h5gcreate_f (atomIntgGroup_gid,"atomDMOverlap",atomDMOverlap_gid,&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to create dipole moment overlap group'
   endif
   call h5gcreate_f (atomIntgGroup_gid,"atomPotOverlap",atomPotOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to create atom potential overlap group'

   ! Create the subgroups in the atomPotOverlap subgroup.  There will be one
   !   subgroup for each kpoint in the atomPotOverlap subgroup.  The other
   !   atom*Overlap groups will not have subgroups.  Instead, each kpoint will
   !   be a dataset.

   ! First, sufficient space must be allocated to hold the gid values for the
   !   subgroups that exist for each kpoint for the atomicPotential
   !   hamiltonian terms.  Also allocate space to hold the dataset IDs.
   allocate (atomPotKPointOL_gid (numKPoints))
   allocate (atomOverlap_did     (numKPoints))
   allocate (atomKEOverlap_did   (numKPoints))
   if (rel == 1) then
      allocate (atomMVOverlap_did   (numKPoints))
   endif
   allocate (atomNucOverlap_did  (numKPoints))
   if (doDIMO == 1) then
      allocate (atomDMKPointOL_gid  (numKPoints))
      allocate (atomDMOverlap_did   (numKPoints,3))
   endif
   allocate (atomPotOverlap_did  (numKPoints,potDim))

   ! Then loop over the kpoints and assign a gid name equal to the kpoint # for
   !   the atomic-potential hamiltonian terms.
   do i = 1, numKPoints
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)

      call h5gcreate_f (atomPotOverlap_gid,currentName,atomPotKPointOL_gid(i),&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to create potential kpoint group'

      if (doDIMO == 1) then
         call h5gcreate_f (atomDMOverlap_gid,currentName,&
               & atomDMKPointOL_gid(i),hdferr)
         if (hdferr /= 0) stop 'Failed to create dipole kpoint group'
      endif
   enddo


   ! Create the dataspace that will be used for each dataset in atomIntgGroup
   !   and all of its subgroups.  The same dataspace definition works for all
   !   of the datasets.
   call h5screate_simple_f(2,atomDims,valeVale_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create vale vale dsid'

   ! Define the properties of the datasets to be made.

   ! Create the property list first.  Then set the properties one at a time.
   call h5pcreate_f      (H5P_DATASET_CREATE_F,valeVale_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create vale vale plid'
   call h5pset_layout_f  (valeVale_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set vale vale plid layout'
   call h5pset_chunk_f   (valeVale_plid,2,atomDimsChunk,hdferr)
   if (hdferr /= 0) stop 'Failed to set vale vale plid chunk size'
!   call h5pset_shuffle_f (valeVale_plid,hdferr)
   call h5pset_deflate_f (valeVale_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set vale vale for deflation'

   ! Create the datasets that will be used for the subgroups of atomIntgGroup.
   do i = 1, numKPoints
      write (currentName,fmt="(i7.7)") i
      call h5dcreate_f (atomOverlap_gid,currentName,H5T_NATIVE_DOUBLE,&
            & valeVale_dsid,atomOverlap_did(i),hdferr,valeVale_plid)
      if (hdferr /= 0) stop 'Failed to create atom overlap did'
      call h5dcreate_f (atomKEOverlap_gid,currentName,H5T_NATIVE_DOUBLE,&
            & valeVale_dsid,atomKEOverlap_did(i),hdferr,valeVale_plid)
      if (hdferr /= 0) stop 'Failed to create KE overlap did'
      if (rel == 1) then
         call h5dcreate_f (atomMVOverlap_gid,currentName,H5T_NATIVE_DOUBLE,&
            & valeVale_dsid,atomMVOverlap_did(i),hdferr,valeVale_plid)
         if (hdferr /= 0) stop 'Failed to create MV overlap did'
      endif
      call h5dcreate_f (atomNucOverlap_gid,currentName,H5T_NATIVE_DOUBLE,&
            & valeVale_dsid,atomNucOverlap_did(i),hdferr,valeVale_plid)
      if (hdferr /= 0) stop 'Failed to create nuclear overlap did'
      if (doDIMO == 1) then
         do j = 1, 3
            write (currentName,fmt="(i7.7)") j
            call h5dcreate_f (atomDMKPointOL_gid(i),currentName,&
               & H5T_NATIVE_DOUBLE,valeVale_dsid,atomDMOverlap_did(i,j),&
               & hdferr,valeVale_plid)
            if (hdferr /= 0) stop 'Failed to create DM overlap did'
         enddo
      endif
      do j = 1, potDim
         write (currentName,fmt="(i7.7)") j
         call h5dcreate_f(atomPotKPointOL_gid(i),currentName,H5T_NATIVE_DOUBLE,&
               & valeVale_dsid,atomPotOverlap_did(i,j),hdferr,valeVale_plid)
         if (hdferr /= 0) stop 'Failed to create potential overlap did'
      enddo
   enddo

end subroutine initSetupIntegralHDF5


subroutine accessSetupIntegralHDF5 (setup_fid)

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_KPoints,     only: numKPoints
   use O_Potential,   only: rel, potDim
   use O_AtomicSites, only: valeDim

   ! Define the passed parameters.
   integer(hid_t) :: setup_fid

   ! Define local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize data structure dimensions.
#ifndef GAMMA
   atomDims(1) = 2 ! Real and imaginary are needed.
#else
   atomDims(1) = 1 ! Real needed only
#endif
   atomDims(2) = valeDim*(valeDim+1)/2 ! Linear storage of 1/2 matrix.

   ! Open the Integral group within the setup HDF5 file.
   call h5gopen_f (setup_fid,"/atomIntgGroup",atomIntgGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom intg group'

   ! Open the subgroups within the atomIntgGroup.
   call h5gopen_f (atomIntgGroup_gid,"atomOverlap",atomOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom overlap group'
   call h5gopen_f (atomIntgGroup_gid,"atomKEOverlap",atomKEOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open kinetic energy overlap group'
   if (rel == 1) then
      call h5gopen_f (atomIntgGroup_gid,"atomMVOverlap",atomMVOverlap_gid,&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to open mass velocity overlap group'
   endif
   call h5gopen_f (atomIntgGroup_gid,"atomNucOverlap",atomNucOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to open atom nuclear overlap group'
   if (doDIMO == 1) then
      call h5gopen_f (atomIntgGroup_gid,"atomDMOverlap",atomDMOverlap_gid,&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to open dipole moment overlap group'
   endif
   call h5gopen_f (atomIntgGroup_gid,"atomPotOverlap",atomPotOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to open atom potential overlap group'

   ! Open the subgroups in the atomPotOverlap subgroup.  There will be one
   !   subgroup for each kpoint in the atomPotOverlap subgroup.  The other
   !   atom*Overlap groups will not have subgroups.  Instead, each kpoint will
   !   be a dataset.

   ! First, sufficient space must be allocated to hold the gid values for the
   !   subgroups that exist for each kpoint for the atomicPotential
   !   hamiltonian terms.  Also allocate space to hold the dataset IDs.
   allocate (atomPotKPointOL_gid (numKPoints))
   allocate (atomDMKPointOL_gid (numKPoints))
   allocate (atomOverlap_did (numKPoints))
   allocate (atomKEOverlap_did (numKPoints))
   if (rel == 1) then
      allocate (atomMVOverlap_did (numKPoints))
   endif
   allocate (atomNucOverlap_did (numKPoints))
   if (doDIMO == 1) then
      allocate (atomDMOverlap_did (numKPoints,3))
   endif
   allocate (atomPotOverlap_did (numKPoints,potDim))

   ! Then loop over the kpoints and assign a gid name equal to the kpoint # for
   !   the atomic-potentian hamiltonian terms.
   do i = 1, numKPoints
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)

      call h5gopen_f (atomPotOverlap_gid,currentName,atomPotKPointOL_gid(i),&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to open potential kpoint group'

      if (doDIMO == 1) then
         call h5gopen_f (atomDMOverlap_gid,currentName,atomDMKPointOL_gid(i),&
               & hdferr)
         if (hdferr /= 0) stop 'Failed to open dipole kpoint group'
      endif
   enddo


   ! Open the datasets that will be used for the subgroups of atomIntgGroup.
   do i = 1, numKPoints
      write (currentName,fmt="(i7.7)") i
      call h5dopen_f (atomOverlap_gid,currentName,atomOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open atom overlap did'
      call h5dopen_f (atomKEOverlap_gid,currentName,atomKEOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open KE overlap did'
      if (rel == 1) then
         call h5dopen_f (atomMVOverlap_gid,currentName,atomMVOverlap_did(i),&
               & hdferr)
        if (hdferr /= 0) stop 'Failed to open MV overlap did'
      endif
      call h5dopen_f (atomNucOverlap_gid,currentName,atomNucOverlap_did(i),&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to open nuclear overlap did'
      if (doDIMO == 1) then
         do j = 1, 3
            write (currentName,fmt="(i7.7)") j
            call h5dopen_f (atomDMKPointOL_gid(i),currentName,&
                  & atomDMOverlap_did(i,j),hdferr)
            if (hdferr /= 0) stop 'Failed to open DM overlap did'
         enddo
      endif
      do j = 1, potDim
         write (currentName,fmt="(i7.7)") j
         call h5dopen_f(atomPotKPointOL_gid(i),currentName,&
               & atomPotOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open potential overlap did'
      enddo
   enddo

   ! Obtain the properties of the datasets that were just opened.  They are all
   !   the same and so only one copy is necessary.  (Actually, this value is
   !   not really used, but in the "close" subroutine we close this id so we
   !   should make sure to have it for both the setup and main calls to the
   !   close routine.)
   call h5dget_create_plist_f (atomOverlap_did(1),valeVale_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain vale vale plid'

   ! Obtain the dataspace that is used for each dataset in atomIntgGroup
   !   and all of its subgroups.  The same dataspace definition works for all
   !   of the datasets.  (Same as for the plist above.)
   call h5dget_space_f(atomOverlap_did(1),valeVale_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain vale vale dsid'

end subroutine accessSetupIntegralHDF5


subroutine closeSetupIntegralHDF5

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_KPoints,   only: numKPoints
   use O_Potential, only: rel, potDim
   use O_CommandLine, only: doDIMO

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Declare local loop variables and error control.
   integer :: i,j
   integer :: hdferr

   ! Close the property list first.
   call h5pclose_f (valeVale_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close valeVale_plid.'

   ! Close the datasets next.
   do i = 1, numKPoints
      call h5dclose_f (atomOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close atomOverlap_did.'

      call h5dclose_f (atomKEOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close atomKEOverlap_did.'

      if (rel == 1) then
         call h5dclose_f (atomMVOverlap_did(i),hdferr)
        if (hdferr /= 0) stop 'Failed to close atomMVOverlap_did.'
      endif

      call h5dclose_f (atomNucOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close atomNucOverlap_did.'

      if (doDIMO == 1) then
         do j = 1, 3
            call h5dclose_f (atomDMOverlap_did(i,j),hdferr)
            if (hdferr /= 0) stop 'Failed to close atomDMOverlap_did.'
         enddo
      endif

      do j = 1, potDim
         call h5dclose_f (atomPotOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to close atomPotOverlap_did.'
      enddo
   enddo

   ! Close the data spaces next.
   call h5sclose_f (valeVale_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close valeVale_dsid.'

   ! Close the groups.
   do i = 1, numKPoints
      call h5gclose_f (atomPotKPointOL_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close atomPotKPointOL_gid.'
   enddo
   call h5gclose_f (atomOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomOverlap_gid.'
   call h5gclose_f (atomKEOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomKEOverlap_gid.'
   if (rel == 1) then
      call h5gclose_f (atomMVOverlap_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atomMVOverlap_gid.'
   endif
   call h5gclose_f (atomNucOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomNucOverlap_gid.'
   if (doDIMO == 1) then
      do i = 1, numKPoints
         call h5gclose_f (atomDMKPointOL_gid(i),hdferr)
         if (hdferr /= 0) stop 'Failed to close atomDMKPointOL_gid.'
      enddo
      call h5gclose_f (atomDMOverlap_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to close atomDMOverlap_gid.'
   endif
   call h5gclose_f (atomPotOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomPotOverlap_gid.'
   call h5gclose_f (atomIntgGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomIntgGroup_gid.'

   ! Deallocate the arrays that hold the id numbers.
   deallocate (atomPotKPointOL_gid)
   deallocate (atomOverlap_did)
   deallocate (atomKEOverlap_did)
   if (rel == 1) then
      deallocate (atomMVOverlap_did)
   endif
   deallocate (atomNucOverlap_did)
   if (doDIMO == 1) then
      deallocate (atomDMKPointOL_gid)
      deallocate (atomDMOverlap_did)
   endif
   deallocate (atomPotOverlap_did)
   
end subroutine closeSetupIntegralHDF5

end module O_SetupIntegralsHDF5
