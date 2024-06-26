module O_SCFHDF5

   ! Use the HDF5 module for HDF5 defined types (e.g. size_t and hid_t).
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Declare the file ID.
   integer(hid_t) :: scf_fid

   ! Declare the property list for the scf file and its associated parameters.
   integer(hid_t)  :: scf_plid
   integer         :: mdc_nelmts  ! Meta-data cache num elements.
   integer(size_t) :: rdcc_nelmts ! Raw-data chunk cache num elements.
   integer(size_t) :: rdcc_nbytes ! Raw-data chunk cache size in bytes.
   real            :: rdcc_w0     ! Raw-data chunk cache weighting parameter.

   ! Declare shared variables that are used for creating attributes that will
   !   track the completion of all datasets.
   integer(hid_t) :: attribInt_dsid ! Attribute dataspace.
   integer(hsize_t), dimension (1) :: attribIntDims ! Dataspace dimensionality

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initSCFDF5 (maxNumRayPoints)

   ! Use necessary modules.
   use O_TimeStamps

   ! Use the HDF5 module.
   use HDF5

   ! Use the subsection object modules for scf.
   use O_SCFIntegralsHDF5
   use O_SCFElecStatHDF5
   use O_SCFExchCorrHDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed dummy variables.
   integer :: maxNumRayPoints

   ! Declare local variables.
   integer :: hdferr
   logical :: file_exists

   ! Log the time we start to setup the SCF HDF5 files.
   call timeStampStart(31)

   ! Initialize the Fortran 90 HDF5 interface.
   call h5open_f(hdferr)
   if (hdferr < 0) stop 'Failed to open HDF library'

   ! Create the property list for the scf hdf5 file and turn off
   !   chunk caching.
   call h5pcreate_f (H5P_FILE_ACCESS_F,scf_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create scf plid.'
   call h5pget_cache_f (scf_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to get scf plid cache settings.'
   call h5pset_cache_f (scf_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to set scf plid cache settings.'

   ! Determine if an HDF5 file already exists for this calculation.
   inquire (file="scf-temp.hdf5", exist=file_exists)

   ! If it does, then access the existing file. If not, then create one.
   if (file_exists .eqv. .true.) then
      ! We are continuing a previous calculation.

      ! Open the HDF5 file for reading / writing.
      call h5fopen_f ("scf-temp.hdf5",H5F_ACC_RDWR_F,scf_fid,hdferr,&
            & scf_plid)
      if (hdferr /= 0) stop 'Failed to open scf-temp.hdf5 file.'

      ! Access the groups of the HDF5 file.
      call accessSCFIntegralHDF5 (scf_fid)
      call accessSCFElecStatHDF5 (scf_fid)
      call accessSCFExchCorrHDF5 (scf_fid,maxNumRayPoints)

   else
      ! We are starting a new calculation.

      ! Create the HDF5 file that will hold all the computed results. This
      !   uses the default file creation and file access properties.
      call h5fcreate_f ("scf-temp.hdf5",H5F_ACC_EXCL_F,scf_fid,hdferr,&
            & H5P_DEFAULT_F,scf_plid)
      if (hdferr /= 0) stop 'Failed to create scf-temp.hdf5 file.'

      ! All datasets will have an attached attribute logging that the
      !   calculation has successfully completed. (Checkpointing.) Thus,
      !   we need to create the shared attribute dataspace.
      attribIntDims(1) = 1
      call h5screate_simple_f (1,attribIntDims(1),attribInt_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to create the attribInt_dsid'

      ! Create the subgroups of the scf hdf5 file. This must be done in this
      !   order due to dependencies on potPot_dsid and others.
      call initSCFIntegralHDF5 (scf_fid,attribInt_dsid,attribIntDims)
      call initSCFElecStatHDF5 (scf_fid,attribInt_dsid,attribIntDims)
      call initSCFExchCorrHDF5 (scf_fid,attribInt_dsid,attribIntDims,&
            & maxNumRayPoints)
   endif


   ! Log the time we finish setting up the SCF HDF5 files.
   call timeStampEnd(31)

end subroutine initSCFHDF5


! I think we may not need this.
!subroutine accessSetupHDF5
!
!   ! Use the HDF5 module.
!   use HDF5
!
!   ! Use the subsection object modules for setup.
!   use O_SetupIntegralsHDF5, only: accessSetupIntegralHDF5
!   use O_SetupExchCorrHDF5, only: accessSetupExchCorrHDF5
!   use O_SetupElecStatHDF5, only: accessSetupElecStatHDF5
!
!   ! Make sure that no funny variables are defined.
!   implicit none
!
!   ! Declare local variables.
!   integer :: hdferr
!
!   ! Create the property list for the setup hdf5 file and turn off
!   !   chunk caching.
!   call h5pcreate_f    (H5P_FILE_ACCESS_F,setup_plid,hdferr)
!   if (hdferr /= 0) stop 'Failed to create setup plid in accessSetupHDF5.'
!   call h5pget_cache_f (setup_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
!         & hdferr)
!   if (hdferr /= 0) stop 'Failed to get setup plid cache settings.'
!   call h5pset_cache_f (setup_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
!   if (hdferr /= 0) stop 'Failed to set setup plid cache settings.'
!
!
!   ! Open the HDF5 file that will hold all the computed results.  This will
!   !   die if the file already exists.  It opens read-only.
!   call h5fopen_f ("setup-temp.hdf5",H5F_ACC_RDONLY_F,setup_fid,hdferr,&
!         & setup_plid)
!   if (hdferr /= 0) stop 'Failed to open setup-temp.hdf5 file.'
!
!   ! Access the subgroups of the setup hdf5 file.
!   call accessSetupIntegralHDF5 (setup_fid)
!   call accessSetupExchCorrHDF5 (setup_fid)
!   call accessSetupElecStatHDF5 (setup_fid)
!
!end subroutine accessSetupHDF5


subroutine closeSCFHDF5

   ! Use the HDF5 module.
   use HDF5

   ! Use the subsection object modules for scf.
   use O_SCFIntegralsHDF5, only: closeSCFIntegralHDF5
   use O_SCFExchCorrHDF5,  only: closeSCFExchCorrHDF5
   use O_SCFElecStatHDF5,  only: closeSCFElecStatHDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare local variables.
   integer :: hdferr

   ! Close access to all the subgroups and their parts.
   call closeSCFIntegralHDF5
   call closeSCFExchCorrHDF5
   call closeSCFElecStatHDF5

   ! Close the property list.
   call h5pclose_f (scf_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close scf_plid.'

   ! Close the file.
   call h5fclose_f (scf_fid,hdferr)
   if (hdferr /= 0) stop 'Failed to close scf_fid.'

end subroutine closeSCFHDF5


end module O_SCFHDF5
