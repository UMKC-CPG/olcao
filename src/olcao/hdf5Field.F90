module O_FieldHDF5

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
   integer(hid_t) :: field_fid

   ! Define the property list for the field file and its associated parameters.
   integer(hid_t)  :: field_plid
   integer         :: mdc_nelmts  ! Meta-data cache num elements.
   integer(size_t) :: rdcc_nelmts ! Raw-data chunk cache num elements.
   integer(size_t) :: rdcc_nbytes ! Raw-data chunk cache size in bytes.
   real            :: rdcc_w0     ! Raw-data chunk cache weighting parameter.

   ! Define the group IDs of the subgroups for the wave function, charge
   !   density, and potential function data.
   integer(hid_t) :: wave_gid
   integer(hid_t) :: rho_gid
   integer(hid_t) :: pot_gid

   ! The dataspaces of the wave, rho, and pot datasets are all the same in all
   !   key respects (type, dimension, etc.) and therefore can be given a static
   !   shared dsid definition now for abc dimensional axes.
   integer(hid_t) :: abc_dsid

   ! Similarly, each of the datasets uses the same property list.
   integer(hid_t) :: abc_plid

   ! Define the array that holds the dimensions of the dataset.
   integer(hsize_t), dimension (3) :: abcDims

   ! Define the array that holds the dimensions of the data chunk.
   integer(hsize_t), dimension (3) :: abcDimsChunk

   ! Presently, the number of datasets under each group is fixed, but as always
   !   it might be a good idea to make it flexible for potential unforseen
   !   future uses.
   integer(hid_t), dimension (4) :: wave_did ! (up,down,real,imag)
   integer(hid_t), dimension (4) :: rho_did !  (up,down,live,diff)
   integer(hid_t), dimension (4) :: pot_did !  (up,down,live,diff)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initFieldHDF5

   ! Use necessary modules.
   use O_TimeStamps
   use O_Lattice, only: numMeshPoints

   ! Use the HDF5 module.
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare local variables.
   integer :: numDataPoints
   integer :: hdferr

   ! Log the time we start to field the HDF5 files.
   call timeStampStart(27)

   ! Initialize the Fortran 90 HDF5 interface.
   call h5open_f(hdferr)
   if (hdferr < 0) stop 'Failed to open HDF library'

   ! Create the property list for the field hdf5 file and turn off
   !   chunk caching.
   call h5pcreate_f    (H5P_FILE_ACCESS_F,field_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field plid.'
   call h5pget_cache_f (field_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to get field plid cache settings.'
   call h5pset_cache_f (field_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to set field plid cache settings.'

   ! Create the HDF5 file that will hold all the computed results.  This will
   !   die if the file already exists.  It also uses the default file creation
   !   and file access properties.
   call h5fcreate_f ("field-temp.hdf5",H5F_ACC_EXCL_F,field_fid,hdferr,&
         & H5P_DEFAULT_F,field_plid)
   if (hdferr /= 0) stop 'Failed to create field-temp.hdf5 file.'

   ! Create the groups of the field hdf5 file.
   call h5gcreate_f (field_fid,"/waveGroup",wave_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create wave group'
   call h5gcreate_f (field_fid,"/rhoGroup",rho_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create charge density (rho) group'
   call h5gcreate_f (field_fid,"/potGroup",pot_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create potential group'

   ! Initialize data structure dimensions.
   abcDims(:) = numMeshPoints(:)

   ! Establish the upper limit for data chunk size. We start with an assumption
   !   that the numbers being stored are 8 byte reals and that we should not go
   !   over 2 billion bytes in a given chunk. Thus, for the three dimensional
   !   data array we require a*b*c < 250M. (From 2e9 / 8.) Without thinking
   !   anything through we will just demand dimensions that meet that
   !   requirement through a simple scale factor. Note that this may not work
   !   correctly for sufficiently large abcDims(:) values because of integer
   !   overflow. Need to write this is a way that is immune to integer overflow
   !   problems.
   numDataPoints = abcDims(1)*abcDims(2)
   if (numDataPoints > 250000000) then
      abcDimsChunk(1:2) = int(250000000/numDataPoints * abcDims(1:2))
   else
      abcDimsChunk(1:2) = abcDims(1:2)
   endif
   abcDimsChunk(3) = 1

   ! Create the dataspace that will be used for each dataset. The same
   !   dataspace definition works for all of the datasets.
   call h5screate_simple_f(3,abcDims,abc_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create abc_dsid'

   ! Define the properties of the datasets to be made.

   ! Create the property list, then set the properties one at a time.
   call h5pcreate_f      (H5P_DATASET_CREATE_F,abc_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create abc_plid'
   call h5pset_layout_f  (abc_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set abc_plid layout'
   call h5pset_chunk_f   (abc_plid,3,abcDimsChunk,hdferr)
   if (hdferr /= 0) stop 'Failed to set abc chunk size'
   call h5pset_deflate_f (abc_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set abc plid for deflation'

   ! Create the datasets that will be used.

   ! First we make all the wave function datasets. This includes the real and
   !   imaginary spin up and spin down data. Note that when a spin non-
   !   polarized calculation is done then the total (up+down) will be stored in
   !   the "up" dataset. Note further that when a gamma k-point calculation is
   !   done only the real_up will be used.
   call h5dcreate_f (wave_gid,"real_up",H5T_NATIVE_DOUBLE,abc_dsid,&
         & wave_did(1),hdferr,abc_plid)
   if (hdferr /= 0) stop 'Failed to create real spin up wave did'
   call h5dcreate_f (wave_gid,"real_dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & wave_did(2),hdferr,abc_plid)
   if (hdferr /= 0) stop 'Failed to create real spin down wave did'
   call h5dcreate_f (wave_gid,"imag_up",H5T_NATIVE_DOUBLE,abc_dsid,&
         & wave_did(3),hdferr,abc_plid)
   if (hdferr /= 0) stop 'Failed to create imaginary spin up wave did'
   call h5dcreate_f (wave_gid,"imag_dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & wave_did(4),hdferr,abc_plid)
   if (hdferr /= 0) stop 'Failed to create imaginary spin down wave did'

   ! Now we will make the charge density (rho) datasets. This includes the
   !   up and down spin charge density obtained from the interacting (live)
   !   material and also the up and down spin constructed from the difference
   !   between the interacting material and that which would be created from
   !   isolated (i.e. non-interacting) atoms. Further, as with the wave
   !   function above, if the calculation is non spin-polarized then the total
   !   (up+down) will be stored in the "up" portion.
   call h5dcreate_f (rho_gid,"live_up",H5T_NATIVE_DOUBLE,abc_dsid,&
         & rho_did(1),hdferr,abc_plid)
   if (hdferr /= 0) stop 'Failed to create interacting (live) spin up rho did'
   call h5dcreate_f (rho_gid,"live_dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & rho_did(2),hdferr,abc_plid)
   if (hdferr /= 0) stop 'Failed to create interacting (live) spin down rho did'
   call h5dcreate_f (rho_gid,"diff_up",H5T_NATIVE_DOUBLE,abc_dsid,&
         & rho_did(3),hdferr,abc_plid)
   if (hdferr /= 0) stop &
         & 'Failed to create interacting-isolated spin up rho difference did'
   call h5dcreate_f (rho_gid,"diff_dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & rho_did(4),hdferr,abc_plid)
   if (hdferr /= 0) stop &
         & 'Failed to create interacting-isolated spin down rho difference did'

   ! Finally, we make the potential datasets. As with the charge density, the
   !   up, down, interacting (live), and interacting-isolated difference atomic
   !   potential function is obtained.
   call h5dcreate_f (pot_gid,"live_up",H5T_NATIVE_DOUBLE,abc_dsid,&
         & pot_did(1),hdferr,abc_plid)
   if (hdferr /= 0) stop 'Failed to create interacting (live) spin up pot did'
   call h5dcreate_f (pot_gid,"live_dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & pot_did(2),hdferr,abc_plid)
   if (hdferr /= 0) stop 'Failed to create interacting (live) spin down pot did'
   call h5dcreate_f (pot_gid,"diff_up",H5T_NATIVE_DOUBLE,abc_dsid,&
         & pot_did(3),hdferr,abc_plid)
   if (hdferr /= 0) stop &
         & 'Failed to create interacting-isolated spin up pot difference did'
   call h5dcreate_f (pot_gid,"diff_dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & pot_did(4),hdferr,abc_plid)
   if (hdferr /= 0) stop &
         & 'Failed to create interacting-isolated spin down pot difference did'

   ! Log the time we finish setting up the HDF5 files.
   call timeStampEnd(27)

end subroutine initFieldHDF5


subroutine accessFieldHDF5

   ! Use necessary modules.
   use O_Lattice, only: numMeshPoints

   ! Use the HDF5 module.
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare local variables.
   integer :: hdferr

   ! Create the property list for the field hdf5 file and turn off
   !   chunk caching.
   call h5pcreate_f    (H5P_FILE_ACCESS_F,field_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field plid in accessFieldHDF5.'
   call h5pget_cache_f (field_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to get field plid cache settings.'
   call h5pset_cache_f (field_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to set field plid cache settings.'


   ! Open the HDF5 file that will hold all the computed results.  This will
   !   die if the file already exists.  It opens read-only.
   call h5fopen_f ("field-temp.hdf5",H5F_ACC_RDONLY_F,field_fid,hdferr,&
         & field_plid)
   if (hdferr /= 0) stop 'Failed to open field-temp.hdf5 file.'

   ! Initialize data structure dimensions.
   abcDims(:) = numMeshPoints(:)

   ! Open the groups of the hdf5 field file.
   call h5gopen_f (field_fid,"/waveGroup",wave_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open waveGroup'
   call h5gopen_f (field_fid,"/rhoGroup",rho_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open rhoGroup'
   call h5gopen_f (field_fid,"/potGroup",pot_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open potGroup'

   ! Open the datasets.

   ! Wave function datasets first.
   call h5dopen_f (wave_gid,"real_up",wave_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to open real up wave did'
   call h5dopen_f (wave_gid,"real_dn",wave_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to open real down wave did'
   call h5dopen_f (wave_gid,"imag_up",wave_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to open imaginary up wave did'
   call h5dopen_f (wave_gid,"imag_dn",wave_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to open imaginary down wave did'

   ! Charge density (rho) datasets second.
   call h5dopen_f (rho_gid,"live_up",rho_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to open live up rho did'
   call h5dopen_f (rho_gid,"live_dn",rho_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to open live down rho did'
   call h5dopen_f (rho_gid,"diff_up",rho_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to open live-isolated up rho difference did'
   call h5dopen_f (rho_gid,"diff_dn",rho_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to open live-isolated down rho difference did'

   ! Potential function datasets last.
   call h5dopen_f (pot_gid,"live_up",pot_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to open live up pot did'
   call h5dopen_f (pot_gid,"live_dn",pot_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to open live down pot did'
   call h5dopen_f (pot_gid,"diff_up",pot_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to open live-isolated up pot difference did'
   call h5dopen_f (pot_gid,"diff_dn",pot_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to open live-isolated down pot difference did'

   ! Obtain the properties of the datasets that were just opened. They are all
   !   the same and so only one copy is necessary. (Actually, this value is
   !   not really used, but in the "close" subroutine we close this id so we
   !   should make sure to have it open before attempting to close it.
   call h5dget_create_plist_f (wave_did(1),abc_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain abc_plid'

   ! Obtain the dataspace that is used for each dataset. The same dataspace
   !   definition works for all of the datasets.
   call h5dget_space_f (wave_did(1),abc_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain abc_dsid'

end subroutine accessFieldHDF5


subroutine closeFieldHDF5

   ! Use the HDF5 module.
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare local variables.
   integer :: hdferr

   ! Close the property list.
   call h5pclose_f (abc_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close abc_plid'

   ! Close the datasets next.

   ! Close the wave function datasets first.
   call h5dclose_f (wave_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to close wave_did(1)'
   call h5dclose_f (wave_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to close wave_did(2)'
   call h5dclose_f (wave_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to close wave_did(3)'
   call h5dclose_f (wave_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to close wave_did(4)'

   ! Close the charge density (rho) datasets second.
   call h5dclose_f (rho_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to close rho_did(1)'
   call h5dclose_f (rho_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to close rho_did(2)'
   call h5dclose_f (rho_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to close rho_did(3)'
   call h5dclose_f (rho_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to close rho_did(4)'

   ! Close the potential function datasets first.
   call h5dclose_f (pot_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to close pot_did(1)'
   call h5dclose_f (pot_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to close pot_did(2)'
   call h5dclose_f (pot_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to close pot_did(3)'
   call h5dclose_f (pot_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to close pot_did(4)'

   ! Close the dataspace next.
   call h5sclose_f (abc_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close abc_dsid'

   ! Close the groups.
   call h5gclose_f (wave_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close wave_gid'
   call h5gclose_f (rho_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close rho_gid'
   call h5gclose_f (pot_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pot_gid'

   ! Close the field property list.
   call h5pclose_f (field_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close field_plid.'

   ! Close the field file.
   call h5fclose_f (field_fid,hdferr)
   if (hdferr /= 0) stop 'Failed to close field_fid.'

end subroutine closeFieldHDF5


end module O_FieldHDF5
