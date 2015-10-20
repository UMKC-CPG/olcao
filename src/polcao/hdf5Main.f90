module O_MainHDF5

   ! Use the HDF5 module for HDF5 defined types (e.g. hsize_t and hid_t).
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define the file ID
   integer(hid_t) :: main_fid

   ! Define the property list for main_fid and its associated parameters.
   integer(hid_t)  :: main_plid
   integer         :: mdc_nelmts  ! Meta-data cache num elements.
   integer(size_t) :: rdcc_nelmts ! Raw-data chunk cache num elements.
   integer(size_t) :: rdcc_nbytes ! Raw-data chunk cache size in bytes.
   real            :: rdcc_w0     ! Raw-data chunk cache weighting parameter.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initMainHDF5(numStates)

   ! Use the HDF5 module.
   use HDF5
   use O_MainEValHDF5, only: initMainEValHDF5
   use O_MainEVecHDF5, only: initMainEVecHDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed variable.
   integer, intent(inout) :: numStates

   ! Define local variables.
   integer :: hdferr

   ! Initialize the Fortran 90 HDF5 interface.
   call h5open_f(hdferr)
   if (hdferr < 0) stop 'Fail to open HDF library'

   ! Begin creating files, groups, dataspaces, property lists, and datasets for
   !   the eigen values/vectors.

   ! Create the property list for the main hdf5 file and turn off
   !   chunk caching.
   call h5pcreate_f    (H5P_FILE_ACCESS_F,main_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create main plid.'
   call h5pget_cache_f (main_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to get main plid cache settings.'
   call h5pset_cache_f (main_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to set main plid cache settings.'


   ! Create the HDF5 file that will hold the computed main results.  This will
   !   die if the file already exists.
   call h5fcreate_f ("main-temp.hdf5",H5F_ACC_EXCL_F,main_fid,hdferr,&
         & H5P_DEFAULT_F,main_plid)
   if (hdferr /= 0) stop 'Failed to create main hdf5 file.  Already exists?'


   ! Then create the necessary subgroups of the main hdf5 file.
   call initMainEVecHDF5 (main_fid,numStates)
   call initMainEValHDF5  (main_fid,numStates)

end subroutine initMainHDF5


subroutine closeMainHDF5

   ! Use the HDF5 module.
   use HDF5
   use O_MainEValHDF5, only: closeMainEValHDF5
   use O_MainEVecHDF5, only: closeMainEVecHDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define local variables.
   integer :: hdferr

   ! Close the main_fid property list.
   call h5pclose_f (main_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the main_plid'

   ! Close the main file.
   call h5fclose_f (main_fid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the main_fid'

   ! Close the subgroups of the main hdf5 file.
   call closeMainEVecHDF5
   call closeMainEValHDF5

end subroutine closeMainHDF5


end module O_MainHDF5
