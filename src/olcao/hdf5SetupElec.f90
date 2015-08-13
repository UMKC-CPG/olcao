module O_SetupElecStatHDF5

   ! Import necessary modules.
   use HDF5
   use O_Potential

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define the main subgroup from setup_fid that holds the mesh information
   !   for the electrostatic calculation.
   integer(hid_t) :: elecStatGroup_gid ! Electrostatic group.

   ! There are no subgroups of the elecStatGroup, only datasets.

   ! The dataspace of each dataset in the elecStatGroup will be one of
   !   three possible kinds, a 2D matrix or a 1D array of potDim size, or a
   !   matrix of potDim and numPotTypes dimension.
   integer(hid_t) :: potTypesPot_dsid
   integer(hid_t) :: potPot_dsid
   integer(hid_t) :: pot_dsid

   ! Each of the below datasets will use one of these property lists.  They
   !   correspond to the above dataspaces.
   integer(hid_t) :: potTypesPot_plid
   integer(hid_t) :: potPot_plid
   integer(hid_t) :: pot_plid

   ! Define arrays that hold the dimensions of the datasets.
   integer(hsize_t), dimension (2) :: potTypesPot
   integer(hsize_t), dimension (2) :: potPot
   integer(hsize_t), dimension (1) :: pot

   ! Define the dataset IDs under elecStatGroup_gid.  There are no
   !   dynamically given dataset IDs.
   integer(hid_t) :: potAlphaOverlap_did
   integer(hid_t) :: nonLocalNeutQPot_did
   integer(hid_t) :: nonLocalNucQPot_did
   integer(hid_t) :: localNeutQPot_did
   integer(hid_t) :: localNucQPot_did
   integer(hid_t) :: nonLocalResidualQ_did
   integer(hid_t) :: coreChargeDensity_did

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initSetupElecStatHDF5 (setup_fid)

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_PotTypes    ! For potDim, numPotTypes

   ! Define the passed parameters.
   integer(hid_t) :: setup_fid

   ! Define local variables.
   integer :: hdferr

   ! Initialize data structure dimensions.
   potPot(1)    = potDim
   potPot(2)    = potDim
   potTypesPot(1) = numPotTypes
   potTypesPot(2) = potDim
   pot(1)    = potDim

   ! Create the electrostatic group within the HDF5 file.
   call h5gcreate_f (setup_fid,"/elecStatGroup",elecStatGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create electrostatics group'

   ! Create the dataspaces that will be used for each dataset in the
   !   elecStatGroup.
   call h5screate_simple_f (2,potTypesPot,potTypesPot_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create potTypesPot dsid'
   call h5screate_simple_f (2,potPot,potPot_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create potPot dsid'
   call h5screate_simple_f (1,pot,pot_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create pot dsid'

   ! Create the property lists first.  Then set the properties for each list
   !   one at a time.
   call h5pcreate_f      (H5P_DATASET_CREATE_F,potTypesPot_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create potTypesPot plid'
   call h5pcreate_f      (H5P_DATASET_CREATE_F,potPot_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create potPot plid'
   call h5pcreate_f      (H5P_DATASET_CREATE_F,pot_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create pot plid'
   call h5pset_layout_f  (potTypesPot_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set potTypesPot chunked layout'
   call h5pset_layout_f  (potPot_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set potPot chunked layout'
   call h5pset_layout_f  (pot_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set pot chunked layout'
   call h5pset_chunk_f   (potTypesPot_plid,2,potTypesPot,hdferr)
   if (hdferr /= 0) stop 'Failed to set potTypesPot chunked property'
   call h5pset_chunk_f   (potPot_plid,2,potPot,hdferr)
   if (hdferr /= 0) stop 'Failed to set potPot chunked property'
   call h5pset_chunk_f   (pot_plid,1,pot,hdferr)
   if (hdferr /= 0) stop 'Failed to set pot chunked property'
!   call h5pset_shuffle_f (potTypesPot_plid,hdferr)
!   call h5pset_shuffle_f (potPot_plid,hdferr)
!   call h5pset_shuffle_f (pot_plid,hdferr)
   call h5pset_deflate_f   (potTypesPot_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set potTypesPot deflate property'
   call h5pset_deflate_f   (potPot_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set potPot deflate property'
   call h5pset_deflate_f   (pot_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set pot deflate property'

   ! Create the datasets that will be used for the elecStatGroup.
   call h5dcreate_f(elecStatGroup_gid,"potAlphaOverlap",H5T_NATIVE_DOUBLE,&
         & potPot_dsid,potAlphaOverlap_did,hdferr,potPot_plid)
   if (hdferr /= 0) stop 'Failed to create potAlphaOverlap did.'
   call h5dcreate_f(elecStatGroup_gid,"nonLocalNeutQPot",H5T_NATIVE_DOUBLE,&
         & potPot_dsid,nonLocalNeutQPot_did,hdferr,potPot_plid)
   if (hdferr /= 0) stop 'Failed to create nonLocalNeutQPot did.'
   call h5dcreate_f(elecStatGroup_gid,"nonLocalNucQPot",H5T_NATIVE_DOUBLE,&
         & pot_dsid,nonLocalNucQPot_did,hdferr,pot_plid)
   if (hdferr /= 0) stop 'Failed to create nonLocalNucQPot did.'
   call h5dcreate_f(elecStatGroup_gid,"localNeutQPot",H5T_NATIVE_DOUBLE,&
         & potPot_dsid,localNeutQPot_did,hdferr,potPot_plid)
   if (hdferr /= 0) stop 'Failed to create localNeutQPot did.'
   call h5dcreate_f(elecStatGroup_gid,"localNucQPot",H5T_NATIVE_DOUBLE,&
         & pot_dsid,localNucQPot_did,hdferr,pot_plid)
   if (hdferr /= 0) stop 'Failed to create localNucQPot did.'
   call h5dcreate_f(elecStatGroup_gid,"nonLocalResidualQ",H5T_NATIVE_DOUBLE,&
         & potTypesPot_dsid,nonLocalResidualQ_did,hdferr,potTypesPot_plid)
   if (hdferr /= 0) stop 'Failed to create nonLocalResidualQ did.'
   call h5dcreate_f(elecStatGroup_gid,"coreChargeDensity",H5T_NATIVE_DOUBLE,&
         & pot_dsid,coreChargeDensity_did,hdferr,pot_plid)
   if (hdferr /= 0) stop 'Failed to create coreChargeDensity did.'

end subroutine initSetupElecStatHDF5



subroutine accessSetupElecStatHDF5 (setup_fid)

   ! Import any necessary definition modules.
   use HDF5

   ! Import necessary object modules.
   use O_PotTypes    ! For potDim, numPotTypes

   ! Define the passed parameters.
   integer(hid_t) :: setup_fid

   ! Define local variables.
   integer :: hdferr

   ! Initialize data structure dimensions.
   potPot(1)      = potDim
   potPot(2)      = potDim
   potTypesPot(1) = numPotTypes
   potTypesPot(2) = potDim
   pot(1)         = potDim

   ! Open the electrostatic group within the HDF5 file.
   call h5gopen_f (setup_fid,"/elecStatGroup",elecStatGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open electrostatics group'

   ! Open the datasets that will be used for the elecStatGroup.
   call h5dopen_f(elecStatGroup_gid,"potAlphaOverlap",potAlphaOverlap_did,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to open potAlphaOverlap did.'
   call h5dopen_f(elecStatGroup_gid,"nonLocalNeutQPot",nonLocalNeutQPot_did,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to open nonLocalNeutQPot did.'
   call h5dopen_f(elecStatGroup_gid,"nonLocalNucQPot",nonLocalNucQPot_did,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to open nonLocalNucQPot did.'
   call h5dopen_f(elecStatGroup_gid,"localNeutQPot",localNeutQPot_did,hdferr)
   if (hdferr /= 0) stop 'Failed to open localNeutQPot did.'
   call h5dopen_f(elecStatGroup_gid,"localNucQPot",localNucQPot_did,hdferr)
   if (hdferr /= 0) stop 'Failed to open localNucQPot did.'
   call h5dopen_f(elecStatGroup_gid,"nonLocalResidualQ",nonLocalResidualQ_did,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to open nonLocalResidualQ did.'
   call h5dopen_f(elecStatGroup_gid,"coreChargeDensity",coreChargeDensity_did,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to open coreChargeDensity did.'

   ! Obtain the property lists.  Read the note in accessSetupIntegralHDF5.
   call h5dget_create_plist_f (nonLocalResidualQ_did,potTypesPot_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain potTypesPot plid'
   call h5dget_create_plist_f (potAlphaOverlap_did,potPot_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain potPot plid'
   call h5dget_create_plist_f (coreChargeDensity_did,pot_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain pot plid'

   ! Obtain the dataspaces that are used for the elecStatGroup datasets.  Read
   !   the note in accessSetupINtegralHDF5.
   call h5dget_space_f (nonLocalResidualQ_did,potTypesPot_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain potTypesPot dsid'
   call h5dget_space_f (potAlphaOverlap_did,potPot_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain potPot dsid'
   call h5dget_space_f (coreChargeDensity_did,pot_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain pot dsid'

end subroutine accessSetupElecStatHDF5


subroutine closeSetupElecStatHDF5

   ! Import any necessary definition modules.
   use HDF5

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Declare local loop variables and error control.
   integer :: hdferr

   ! Close all access to the electrostatic calculation parts of the HDF file
   !   for setup.

   ! Close the property lists first.
   call h5pclose_f (potTypesPot_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close potTypesPot_plid.'
   call h5pclose_f (potPot_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close potPot_plid.'
   call h5pclose_f (pot_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pot_plid.'

   ! Close the datasets next.
   call h5dclose_f (potAlphaOverlap_did,hdferr)
   if (hdferr /= 0) stop 'Failed to close potAlphaOverlap_did.'
   call h5dclose_f (nonLocalNeutQPot_did,hdferr)
   if (hdferr /= 0) stop 'Failed to close nonLocalNeutQPot_did.'
   call h5dclose_f (nonLocalNucQPot_did,hdferr)
   if (hdferr /= 0) stop 'Failed to close nonLocalNucQPot_did.'
   call h5dclose_f (localNeutQPot_did,hdferr)
   if (hdferr /= 0) stop 'Failed to close localNeutQPot_did.'
   call h5dclose_f (localNucQPot_did,hdferr)
   if (hdferr /= 0) stop 'Failed to close localNucQPot_did.'
   call h5dclose_f (nonLocalResidualQ_did,hdferr)
   if (hdferr /= 0) stop 'Failed to close nonLocalResidualQ_did.'

   ! Close the data spaces next.
   call h5sclose_f (potTypesPot_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close potTypesPot_dsid.'
   call h5sclose_f (potPot_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close potPot_dsid.'
   call h5sclose_f (pot_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pot_dsid.'

   ! Close the group.
   call h5gclose_f (elecStatGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close elecStatGroup_gid.'

end subroutine closeSetupElecStatHDF5


end module O_SetupElecStatHDF5
