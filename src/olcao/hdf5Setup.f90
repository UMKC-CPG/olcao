module O_SetupHDF5

   ! Use the HDF5 module for HDF5 defined types (e.g. size_t and hid_t).
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define the file ID
   integer(hid_t) :: setup_fid

   ! Define the property list for the setup file and its associated parameters.
   integer(hid_t)  :: setup_plid
   integer         :: mdc_nelmts  ! Meta-data cache num elements.
   integer(size_t) :: rdcc_nelmts ! Raw-data chunk cache num elements.
   integer(size_t) :: rdcc_nbytes ! Raw-data chunk cache size in bytes.
   real            :: rdcc_w0     ! Raw-data chunk cache weighting parameter.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initSetupHDF5 (maxNumRayPoints)

   ! Use necessary modules.
   use O_TimeStamps

   ! Use the HDF5 module.
   use HDF5

   ! Use the subsection object modules for setup.
   use O_SetupIntegralsHDF5, only: initSetupIntegralHDF5
   use O_SetupElecStatHDF5,  only: initSetupElecStatHDF5
   use O_SetupExchCorrHDF5,  only: initSetupExchCorrHDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed dummy variables.
   integer :: maxNumRayPoints

   ! Declare local variables.
   integer :: hdferr

   ! Log the time we start to setup the HDF5 files.
   call timeStampStart(6)

   ! Initialize the Fortran 90 HDF5 interface.
   call h5open_f(hdferr)
   if (hdferr < 0) stop 'Failed to open HDF library'

   ! Create the property list for the setup hdf5 file and turn off
   !   chunk caching.
   call h5pcreate_f    (H5P_FILE_ACCESS_F,setup_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create setup plid.'
   call h5pget_cache_f (setup_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to get setup plid cache settings.'
   call h5pset_cache_f (setup_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to set setup plid cache settings.'


   ! Create the HDF5 file that will hold all the computed results.  This will
   !   die if the file already exists.  It also uses the default file creation
   !   and file access properties.
   call h5fcreate_f ("setup-temp.hdf5",H5F_ACC_EXCL_F,setup_fid,hdferr,&
         & H5P_DEFAULT_F,setup_plid)
   if (hdferr /= 0) stop 'Failed to create setup-temp.hdf5 file.'

   ! Create the subgroups of the setup hdf5 file and define all their parts.
   !   This must be done in this order due to dependencies on potPot_dsid and
   !   others.
   call initSetupIntegralHDF5 (setup_fid)
   call initSetupElecStatHDF5 (setup_fid)
   call initSetupExchCorrHDF5 (setup_fid,maxNumRayPoints)

   ! Log the time we finish setting up the HDF5 files.
   call timeStampEnd(6)

end subroutine initSetupHDF5


subroutine accessSetupHDF5

   ! Use the HDF5 module.
   use HDF5

   ! Use the subsection object modules for setup.
   use O_SetupIntegralsHDF5, only: accessSetupIntegralHDF5
   use O_SetupExchCorrHDF5, only: accessSetupExchCorrHDF5
   use O_SetupElecStatHDF5, only: accessSetupElecStatHDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare local variables.
   integer :: hdferr

   ! Create the property list for the setup hdf5 file and turn off
   !   chunk caching.
   call h5pcreate_f    (H5P_FILE_ACCESS_F,setup_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create setup plid in accessSetupHDF5.'
   call h5pget_cache_f (setup_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to get setup plid cache settings.'
   call h5pset_cache_f (setup_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to set setup plid cache settings.'


   ! Open the HDF5 file that will hold all the computed results.  This will
   !   die if the file already exists.  It opens read-only.
   call h5fopen_f ("setup-temp.hdf5",H5F_ACC_RDONLY_F,setup_fid,hdferr,&
         & setup_plid)
   if (hdferr /= 0) stop 'Failed to open setup-temp.hdf5 file.'

   ! Access the subgroups of the setup hdf5 file.
   call accessSetupIntegralHDF5 (setup_fid)
   call accessSetupExchCorrHDF5 (setup_fid)
   call accessSetupElecStatHDF5 (setup_fid)

end subroutine accessSetupHDF5


subroutine closeSetupHDF5

   ! Use the HDF5 module.
   use HDF5

   ! Use the subsection object modules for setup.
   use O_SetupIntegralsHDF5, only: closeSetupIntegralHDF5
   use O_SetupExchCorrHDF5,  only: closeSetupExchCorrHDF5
   use O_SetupElecStatHDF5,  only: closeSetupElecStatHDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare local variables.
   integer :: hdferr

   ! Close access to all the subgroups and their parts.
   call closeSetupIntegralHDF5
   call closeSetupExchCorrHDF5
   call closeSetupElecStatHDF5

   ! Close the property list.
   call h5pclose_f (setup_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close setup_plid.'

   ! Close the file.
   call h5fclose_f (setup_fid,hdferr)
   if (hdferr /= 0) stop 'Failed to close setup_fid.'

end subroutine closeSetupHDF5


end module O_SetupHDF5
