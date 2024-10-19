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

   ! Declare shared variables that are used for creating attributes that will
   !   track the completion of all datasets.
   integer(hid_t) :: attribInt_dsid ! Attribute dataspace.
   integer(hsize_t), dimension (1) :: attribIntDims ! Dataspace dimensionality

   ! Define the group IDs of the subgroups for the wave function, charge
   !   density, and potential function data.
   integer(hid_t) :: wav_gid
   integer(hid_t) :: rho_gid
   integer(hid_t) :: pot_gid

   ! The dataspaces of the wav, rho, and pot datasets are all the same in all
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
   !   future uses. Note: "sum" = spin sum (up+dn) and "diff" = spin
   !   difference (up-dn). Also: "live" = from interacting system and "N" =
   !   from neutral atoms.
   integer(hid_t), dimension (4) :: wav_did ! (sum,diff,live-N sum,live-N diff)
   integer(hid_t), dimension (4) :: rho_did ! (sum,diff,live-N sum,live-N diff)
   integer(hid_t), dimension (4) :: pot_did ! (sum,diff,live-N sum,live-N diff)

   ! Integer to identify which axis is the one that needs to be tracked for the
   !   determination of when to write data chunks.
   integer :: triggerAxis

   ! Because the data may be large we don't store it all in memory at the same
   !   time. Thus, we need to write only portions to disk at a time. Therefore,
   !   we need to define a chunk dataspace for the file.
   integer (hid_t) :: fileFieldChunk_dsid

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


subroutine prepFieldHDF5 (inSCF)
   
   ! Use necessary modules.
   use HDF5
   use O_CommandLine, only: excitedQN_n, excitedQN_l, basisCode_SCF, &
         & basisCode_PSCF

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: inSCF

   ! Declare local variables.
   integer :: hdferr
   logical :: file_exists
   character*14 :: fileName
   character*2 :: edge

   ! The HDF5 interface should already be open.

   ! Identify the file name of the Field HDF5 file to create/access.
   if (excitedQN_n == 0) then
      write(edge,fmt="(a)") "gs"
   else
      if (excitedQN_l == 0) then
         write(edge,fmt="(i1,a1)") excitedQN_n, "s"
      elseif (excitedQN_l == 1) then
         write(edge,fmt="(i1,a1)") excitedQN_n, "p"
      elseif (excitedQN_l == 2) then
         write(edge,fmt="(i1,a1)") excitedQN_n, "d"
      elseif (excitedQN_l == 3) then
         write(edge,fmt="(i1,a1)") excitedQN_n, "f"
      elseif (excitedQN_l == 4) then
         write(edge,fmt="(i1,a1)") excitedQN_n, "g"
      endif
   endif
   if (inSCF == 1) then
      if (basisCode_SCF == 1) then
         write(fileName,fmt="(a2,a12)") edge,"_scf+field-mb.hdf5"
      elseif (basisCode_SCF == 2) then
         write(fileName,fmt="(a2,a12)") edge,"_scf+field-fb.hdf5"
      elseif (basisCode_SCF == 3) then
         write(fileName,fmt="(a2,a12)") edge,"_scf+field-eb.hdf5"
      endif
   else
      if (basisCode_PSCF == 1) then
         write(fileName,fmt="(a2,a12)") edge,"_pscf+field-mb.hdf5"
      elseif (basisCode_PSCF == 2) then
         write(fileName,fmt="(a2,a12)") edge,"_pscf+field-fb.hdf5"
      elseif (basisCode_PSCF == 3) then
         write(fileName,fmt="(a2,a12)") edge,"_pscf+field-eb.hdf5"
      endif
   endif

   ! Create the property list for the field hdf5 file and turn off
   !   chunk caching.
   call h5pcreate_f (H5P_FILE_ACCESS_F,field_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field plid.'
   call h5pget_cache_f (field_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to get field plid cache settings.'
   call h5pset_cache_f (field_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to set field plid cache settings.'

   ! Determine if an HDF5 file already exists for this calculation.
   inquire (file=fileName, exist=file_exists)

   ! If it does, then access the existing file. If not, then create one.
   if (file_exists .eqv. .true.) then
      ! We are continuing a previous calculation.

      ! Open the HDF5 file for reading / writing.
      call h5fopen_f (fileName,H5F_ACC_RDWR_F,field_fid,hdferr,field_plid)
      if (hdferr /= 0) stop 'Failed to open field hdf5 file.'

      ! Access the groups of the HDF5 file.
      call accessFieldHDF5 (field_fid)
   else
      ! We are starting a new calculation.

      ! Create the HDF5 file that will hold all the computed results. This
      !   uses the default file creation and file access properties.
      call h5fcreate_f (fileName,H5F_ACC_EXCL_F,field_fid,hdferr,&
            & H5P_DEFAULT_F,field_plid)
      if (hdferr /= 0) stop 'Failed to create field hdf5 file.'

      ! All datasets will have an attached attribute logging that the
      !   calculation has successfully completed. (Checkpointing.) Thus,
      !   we need to create the shared attribute dataspace.
      attribIntDims(1) = 1
      call h5screate_simple_f (1,attribIntDims(1),attribInt_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to create the attribInt_dsid'

      ! Create the subgroups of the field hdf5 file.
      call initFieldHDF5 (field_fid)
   endif

end subroutine prepFieldHDF5


subroutine initFieldHDF5 (field_fid)

   ! Use necessary modules.
   use O_Lattice, only: numMeshPoints

   ! Use the HDF5 module.
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed parameters.
   integer(hid_t) :: field_fid

   ! Declare local variables.
   integer :: maxNumDataPoints
   integer :: hdferr

!   ! Initialize the Fortran 90 HDF5 interface.
!   call h5open_f(hdferr)
!   if (hdferr < 0) stop 'Failed to open HDF library'
!
!   ! Create the property list for the field hdf5 file and turn off
!   !   chunk caching.
!   call h5pcreate_f    (H5P_FILE_ACCESS_F,field_plid,hdferr)
!   if (hdferr /= 0) stop 'Failed to create field plid.'
!   call h5pget_cache_f (field_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
!         & hdferr)
!   if (hdferr /= 0) stop 'Failed to get field plid cache settings.'
!   call h5pset_cache_f (field_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
!   if (hdferr /= 0) stop 'Failed to set field plid cache settings.'

!   ! Create the HDF5 file that will hold all the computed results.  This will
!   !   die if the file already exists.  It also uses the default file creation
!   !   and file access properties.
!   call h5fcreate_f ("field-temp.hdf5",H5F_ACC_EXCL_F,field_fid,hdferr,&
!         & H5P_DEFAULT_F,field_plid)
!   if (hdferr /= 0) stop 'Failed to create field-temp.hdf5 file.'

   ! Create the groups of the field hdf5 file.
   call h5gcreate_f (field_fid,"/wavGroup",wav_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field wav group'
   call h5gcreate_f (field_fid,"/rhoGroup",rho_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field charge density (rho) group'
   call h5gcreate_f (field_fid,"/potGroup",pot_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field potential group'

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
   ! Compute the largest chunk size that contains fewer than the hard-coded
   !   maximum number of datapoints yet that still retains the largest extent
   !   along the a, b, and c axes as possible in that priority order.
   maxNumDataPoints = 250000000

   ! Check that the minimum space requirement is met to store data easily
   if (abcDims(1) > maxNumDataPoints) then  ! Note the greater-than sign.
      stop 'Increase maxNumDataPoints'
   else
      ! Assume that the chunk dimensions will be limited by the a-axis.
      abcDimsChunk(1) = abcDims(1)
      abcDimsChunk(2:3) = 1
      triggerAxis = 1

      ! Check if we can store two dimensions at a time.
      if (abcDims(1)*abcDims(2) < maxNumDataPoints) then ! Note less-than.
         abcDimsChunk(1:2) = abcDims(1:2)
         abcDimsChunk(3) = 1
         triggerAxis = 2

         ! Check if we can store all three dimensions.
         if (abcDims(1)*abcDims(2)*abcDims(3) < maxNumDataPoints) then ! Note <
            abcDimsChunk(1:3) = abcDims(1:3)
            triggerAxis = 3
         endif
      endif
   endif

   ! Create the dataspace that will be used for each dataset within the file.
   !   The same dataspace definition works for all of the datasets.
   call h5screate_simple_f(3,abcDims,abc_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field abc_dsid'

   ! Create the dataspace that will be used to describe the data in memory
   !   before being written to a file.
   call h5screate_simple_f(3,abcDimsChunk,fileFieldChunk_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field fileFieldChunk_dsid'

   ! Define the properties of the datasets to be made.

   ! Create the property list, then set the properties one at a time.
   call h5pcreate_f      (H5P_DATASET_CREATE_F,abc_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field abc_plid'
   call h5pset_layout_f  (abc_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set field abc_plid layout'
   call h5pset_chunk_f   (abc_plid,3,abcDimsChunk,hdferr)
   if (hdferr /= 0) stop 'Failed to set field abc chunk size'
   call h5pset_deflate_f (abc_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set field abc plid for deflation'

   ! Create the datasets that will be used.

   ! Possible Future: Not currently implemented this way.
   ! First we make all the wave function datasets. This includes the real and
   !   imaginary spin sum and spin diff data. Note that when a spin non-
   !   polarized calculation is done then the total (up+down) will be stored in
   !   the "sum" dataset. Note further that when a gamma k-point calculation is
   !   done only the real_up+dn and (if needed for gamma spin-polarized
   !   calculations) real_up-dn will be used.
   !
   ! Current implementation:
   ! First we make all the wave function squared datasets. This includes the
   !   sum and difference spin psi^2 obtained from the interacting (live)
   !   material and also sum and diff spin psi^2 constructed from the difference
   !   between the interacting material and that which would be created from
   !   isolated (i.e. non-interacting) atoms. Note that when a spin non-
   !   polarized calculation is done then the total (up+down) will be stored in
   !   the "sum" dataset.
   call h5dcreate_f (wav_gid,"wav_live_up+dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & wav_did(1),hdferr,abc_plid)
   if (hdferr /= 0) stop 'Failed to create field live spin sum wav did'
   call h5dcreate_f (wav_gid,"wav_live_up-dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & wav_did(2),hdferr,abc_plid)
   if (hdferr /= 0) stop 'Failed to create field live spin diff wav did'
   call h5dcreate_f (wav_gid,"wav_diff_up+dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & wav_did(3),hdferr,abc_plid)
   if (hdferr /= 0) stop &
         & 'Failed to create field psi^2 spin sum wav difference did'
   call h5dcreate_f (wav_gid,"wav_diff_up-dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & wav_did(4),hdferr,abc_plid)
   if (hdferr /= 0) stop &
         & 'Failed to create field psi^2 spin diff wav difference did'

   ! Now we will make the charge density (rho) datasets. This includes the
   !   sum and diff spin charge density obtained from the interacting (live)
   !   material and also sum and diff spin data constructed from the difference
   !   between the interacting material and that which would be created from
   !   isolated (i.e. non-interacting) atoms. Further, as with the wave
   !   function above, if the calculation is non spin-polarized then the total
   !   (up+down) will be stored in the "up+dn" portion.
   call h5dcreate_f (rho_gid,"rho_live_up+dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & rho_did(1),hdferr,abc_plid)
   if (hdferr /= 0) stop &
         & 'Failed to create field interacting (live) spin sum rho did'
   call h5dcreate_f (rho_gid,"rho_live_up-dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & rho_did(2),hdferr,abc_plid)
   if (hdferr /= 0) stop &
         & 'Failed to create field interacting (live) spin diff rho did'
   call h5dcreate_f (rho_gid,"rho_diff_up+dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & rho_did(3),hdferr,abc_plid)
   if (hdferr /= 0) stop 'Failed to create field interacting-isolated ' // &
         & 'spin sum rho difference did'
   call h5dcreate_f (rho_gid,"rho_diff_up-dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & rho_did(4),hdferr,abc_plid)
   if (hdferr /= 0) stop 'Failed to create field interacting-isolated ' // &
         & 'spin diff rho difference did'

   ! Finally, we make the potential datasets. As with the charge density, the
   !   sum, diff, interacting (live), and interacting-isolated difference
   !   atomic potential function is obtained.
   call h5dcreate_f (pot_gid,"pot_live_up+dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & pot_did(1),hdferr,abc_plid)
   if (hdferr /= 0) stop 'Failed to create field interacting (live) spin ' // &
         & 'sum pot did'
   call h5dcreate_f (pot_gid,"pot_live_up-dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & pot_did(2),hdferr,abc_plid)
   if (hdferr /= 0) stop 'Failed to create field interacting (live) spin ' // &
         & 'diff pot did'
   call h5dcreate_f (pot_gid,"pot_diff_up+dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & pot_did(3),hdferr,abc_plid)
   if (hdferr /= 0) stop 'Failed to create field interacting-isolated ' // &
         & 'spin sum pot difference did'
   call h5dcreate_f (pot_gid,"pot_diff_up-dn",H5T_NATIVE_DOUBLE,abc_dsid,&
         & pot_did(4),hdferr,abc_plid)
   if (hdferr /= 0) stop 'Failed to create field interacting-isolated ' // &
         & 'spin diff pot difference did'

end subroutine initFieldHDF5


subroutine accessFieldHDF5(field_fid)

   ! Use necessary modules.
   use O_Lattice, only: numMeshPoints

   ! Use the HDF5 module.
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed parameters.
   integer(hid_t) :: field_fid

   ! Declare local variables.
   integer :: hdferr

!   ! Create the property list for the field hdf5 file and turn off
!   !   chunk caching.
!   call h5pcreate_f    (H5P_FILE_ACCESS_F,field_plid,hdferr)
!   if (hdferr /= 0) stop 'Failed to create field plid in accessFieldHDF5.'
!   call h5pget_cache_f (field_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
!         & hdferr)
!   if (hdferr /= 0) stop 'Failed to get field plid cache settings.'
!   call h5pset_cache_f (field_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
!   if (hdferr /= 0) stop 'Failed to set field plid cache settings.'
!
!
!   ! Open the HDF5 file that will hold all the computed results.  This will
!   !   die if the file already exists.  It opens read-only.
!   call h5fopen_f ("field-temp.hdf5",H5F_ACC_RDONLY_F,field_fid,hdferr,&
!         & field_plid)
!   if (hdferr /= 0) stop 'Failed to open field-temp.hdf5 file.'

   ! Initialize data structure dimensions.
   abcDims(:) = numMeshPoints(:)

   ! Open the groups of the hdf5 field file.
   call h5gopen_f (field_fid,"/wavGroup",wav_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open field wavGroup'
   call h5gopen_f (field_fid,"/rhoGroup",rho_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open field rhoGroup'
   call h5gopen_f (field_fid,"/potGroup",pot_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open field potGroup'

   ! Open the datasets.

   ! Wave function |psi|^2 datasets first.
   call h5dopen_f (wav_gid,"wav_live_up+dn",wav_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to open  live spin sum wav did'
   call h5dopen_f (wav_gid,"wav_live_up-dn",wav_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to open  live spin diff wav did'
   call h5dopen_f (wav_gid,"wav_diff_up+dn",wav_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to open  live-iso spin sum wav did'
   call h5dopen_f (wav_gid,"wav_diff_up-dn",wav_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to open  live-iso spin diff wav did'

   ! Charge density (rho) datasets second.
   call h5dopen_f (rho_gid,"rho_live_up+dn",rho_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to open  live spin sum rho did'
   call h5dopen_f (rho_gid,"rho_live_up-dn",rho_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to open  live spin diff rho did'
   call h5dopen_f (rho_gid,"rho_diff_up+dn",rho_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to open  live-iso spin sum rho did'
   call h5dopen_f (rho_gid,"rho_diff_up-dn",rho_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to open  live-iso spin diff rho did'

   ! Potential function datasets last.
   call h5dopen_f (pot_gid,"pot_live_up+dn",pot_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to open  live spin sum pot did'
   call h5dopen_f (pot_gid,"pot_live_up-dn",pot_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to open  live spin diff pot did'
   call h5dopen_f (pot_gid,"pot_diff_up+dn",pot_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to open  live-iso spin sum pot did'
   call h5dopen_f (pot_gid,"pot_diff_up-dn",pot_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to open  live-iso spin diff pot did'

   ! Obtain the properties of the datasets that were just opened. They are all
   !   the same and so only one copy is necessary. (Actually, this value is
   !   not really used, but in the "close" subroutine we close this id so we
   !   should make sure to have it open before attempting to close it.
   call h5dget_create_plist_f (wav_did(1),abc_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain field abc_plid'

   ! Obtain the dataspace that is used for each dataset. The same dataspace
   !   definition works for all of the datasets.
   call h5dget_space_f (wav_did(1),abc_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain field abc_dsid'

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
   if (hdferr /= 0) stop 'Failed to close field abc_plid'

   ! Close the datasets next.

   ! Close the wave function datasets first.
   call h5dclose_f (wav_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to close field wav_did(1)'
   call h5dclose_f (wav_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to close field wav_did(2)'
   call h5dclose_f (wav_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to close field wav_did(3)'
   call h5dclose_f (wav_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to close field wav_did(4)'

   ! Close the charge density (rho) datasets second.
   call h5dclose_f (rho_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to close field rho_did(1)'
   call h5dclose_f (rho_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to close field rho_did(2)'
   call h5dclose_f (rho_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to close field rho_did(3)'
   call h5dclose_f (rho_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to close field rho_did(4)'

   ! Close the potential function datasets first.
   call h5dclose_f (pot_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to close field pot_did(1)'
   call h5dclose_f (pot_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to close field pot_did(2)'
   call h5dclose_f (pot_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to close field pot_did(3)'
   call h5dclose_f (pot_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to close field pot_did(4)'

   ! Close the dataspace next.
   call h5sclose_f (abc_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close field abc_dsid'

   ! Close the groups.
   call h5gclose_f (wav_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close field wav_gid'
   call h5gclose_f (rho_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close field rho_gid'
   call h5gclose_f (pot_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close field pot_gid'

   ! Close the field property list.
   call h5pclose_f (field_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close field field_plid.'

   ! Close the field file.
   call h5fclose_f (field_fid,hdferr)
   if (hdferr /= 0) stop 'Failed to close field_fid.'

end subroutine closeFieldHDF5


end module O_FieldHDF5
