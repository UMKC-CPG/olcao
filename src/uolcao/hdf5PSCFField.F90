module O_PSCFFieldHDF5

   ! Use the HDF5 module for HDF5 defined types (e.g. size_t and hid_t).
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   ! Define the file ID
!   integer(hid_t) :: pscf_fid
!
!   ! Define the property list for the field file and its associated parameters.
!   integer(hid_t)  :: fieldPSCF_plid
!   integer         :: mdc_nelmts  ! Meta-data cache num elements.
!   integer(size_t) :: rdcc_nelmts ! Raw-data chunk cache num elements.
!   integer(size_t) :: rdcc_nbytes ! Raw-data chunk cache size in bytes.
!   real            :: rdcc_w0     ! Raw-data chunk cache weighting parameter.

   ! Define the group IDs of the subgroups for the wave function, charge
   !   density, and potential function data. Then make the group IDs for the
   !   mesh.
   integer(hid_t) :: wavPSCF_gid
   integer(hid_t) :: rhoPSCF_gid
   integer(hid_t) :: potPSCF_gid
   integer(hid_t) :: meshPSCF_gid

   ! The dataspaces of the wav, rho, and pot datasets are all the same in all
   !   key respects (type, dimension, etc.) and therefore can be given a static
   !   shared dsid definition now for abc dimensional axes.
   integer(hid_t) :: abcPSCF_dsid

   ! Similarly, each of the datasets uses the same property list.
   integer(hid_t) :: abcPSCF_plid

   ! Now, make the dsid and plid for the mesh itself.
   integer(hid_t), dimension(3) :: abcMeshPSCF_dsid
   integer(hid_t) :: abcMeshPSCF_plid

   ! Define the array that holds the dimensions of the dataset.
   integer(hsize_t), dimension (3) :: abcDims

   ! Define the array that holds the dimensions of the data chunk.
   integer(hsize_t), dimension (3) :: abcDimsChunkPSCF
   integer(hsize_t), dimension (1) :: abcMaxDimChunkPSCF

   ! Presently, the number of datasets under each group is fixed, but as always
   !   it might be a good idea to make it flexible for potential unforseen
   !   future uses. Note: "sum" = spin sum (up+dn) and "diff" = spin
   !   difference (up-dn). Also: "live" = from interacting system and "N" =
   !   from neutral atoms. Each did dimension is for:
   !   (sum,diff,live-N sum,live-N diff). Then make the mesh did.
   integer(hid_t), dimension (4) :: wavPSCF_did
   integer(hid_t), dimension (4) :: rhoPSCF_did
   integer(hid_t), dimension (4) :: potPSCF_did
   integer(hid_t), dimension (3) :: meshPSCF_did

   ! Integer to identify which axis is the one that needs to be tracked for the
   !   determination of when to write data chunks.
   integer :: triggerAxisPSCF

   ! Because the data may be large we don't store it all in memory at the same
   !   time. Thus, we need to write only portions to disk at a time. Therefore,
   !   we need to define a chunk dataspace for the file.
   integer (hid_t) :: fileFieldChunkPSCF_dsid

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initPSCFFieldHDF5 (pscf_fid)

   ! Use necessary modules.
   use O_Lattice, only: numMeshPoints

   ! Use the HDF5 module.
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed parameters.
   integer(hid_t) :: pscf_fid

   ! Declare local variables.
   integer :: maxNumDataPoints
   integer :: hdferr

!   ! Initialize the Fortran 90 HDF5 interface.
!   call h5open_f(hdferr)
!   if (hdferr < 0) stop 'Failed to open HDF library'
!
!   ! Create the property list for the field hdf5 file and turn off
!   !   chunk caching.
!   call h5pcreate_f    (H5P_FILE_ACCESS_F,fieldPSCF_plid,hdferr)
!   if (hdferr /= 0) stop 'Failed to create field plid.'
!   call h5pget_cache_f (fieldPSCF_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
!         & hdferr)
!   if (hdferr /= 0) stop 'Failed to get field plid cache settings.'
!   call h5pset_cache_f (fieldPSCF_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
!   if (hdferr /= 0) stop 'Failed to set field plid cache settings.'
!
!   ! Create the HDF5 file that will hold all the computed results.  This will
!   !   die if the file already exists.  It also uses the default file creation
!   !   and file access properties.
!   call h5fcreate_f ("field-temp.hdf5",H5F_ACC_EXCL_F,pscf_fid,hdferr,&
!         & H5P_DEFAULT_F,fieldPSCF_plid)
!   if (hdferr /= 0) stop 'Failed to create field-temp.hdf5 file.'

   ! Create the groups of the field hdf5 file.
   call h5gcreate_f (pscf_fid,"/wavGroup",wavPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create pscf wav group'
   call h5gcreate_f (pscf_fid,"/rhoGroup",rhoPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create pscf charge density (rho) group'
   call h5gcreate_f (pscf_fid,"/potGroup",potPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create pscf potential group'

   ! Create the group that define the mesh.
   call h5gcreate_f (pscf_fid,"/mesh",meshPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create pscf mesh group'


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
      abcDimsChunkPSCF(1) = abcDims(1)
      abcDimsChunkPSCF(2:3) = 1
      triggerAxisPSCF = 1

      ! Check if we can store two dimensions at a time.
      if (abcDims(1)*abcDims(2) < maxNumDataPoints) then ! Note less-than.
         abcDimsChunkPSCF(1:2) = abcDims(1:2)
         abcDimsChunkPSCF(3) = 1
         triggerAxisPSCF = 2

         ! Check if we can store all three dimensions.
         if (abcDims(1)*abcDims(2)*abcDims(3) < maxNumDataPoints) then ! Note <
            abcDimsChunkPSCF(1:3) = abcDims(1:3)
            triggerAxisPSCF = 3
         endif
      endif
   endif

   ! Get the largest dimension as use that for the chunk size of the mesh.
   abcMaxDimChunkPSCF(1) = maxval(abcDimsChunkPSCF(:))

   ! Create the dataspace that will be used for each mesh dataset within the
   !   file. The same dataspace definition works for all of the datasets.
   call h5screate_simple_f(3,abcDims,abcPSCF_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create pscf abcPSCF_dsid'

   ! Create the dataspace that will be used to describe the data in memory
   !   before being written to a file.
   call h5screate_simple_f(3,abcDimsChunkPSCF,fileFieldChunkPSCF_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create pscf fileFieldChunkPSCF_dsid'

   ! Similarly, create the dataspaces for the mesh in the file.
   call h5screate_simple_f(1,abcDims(1),abcMeshPSCF_dsid(1),hdferr)
   if (hdferr /= 0) stop 'Failed to create pscf abcMeshXPSCF_dsid'
   call h5screate_simple_f(1,abcDims(2),abcMeshPSCF_dsid(2),hdferr)
   if (hdferr /= 0) stop 'Failed to create pscf abcMeshYPSCF_dsid'
   call h5screate_simple_f(1,abcDims(3),abcMeshPSCF_dsid(3),hdferr)
   if (hdferr /= 0) stop 'Failed to create pscf abcMeshZPSCF_dsid'

   ! Define the properties of the datasets to be made.

   ! Create the property list, then set the properties one at a time.
   call h5pcreate_f      (H5P_DATASET_CREATE_F,abcPSCF_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create pscf abcPSCF_plid'
   call h5pset_layout_f  (abcPSCF_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set pscf abcPSCF_plid layout'
   call h5pset_chunk_f   (abcPSCF_plid,3,abcDimsChunkPSCF,hdferr)
   if (hdferr /= 0) stop 'Failed to set pscf abc chunk size'
   call h5pset_deflate_f (abcPSCF_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set pscf abc plid for deflation'

   ! Create a similar property list for the mesh axes.
   call h5pcreate_f      (H5P_DATASET_CREATE_F,abcMeshPSCF_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create pscf abcMeshPSCF_plid'
   call h5pset_layout_f  (abcMeshPSCF_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set pscf abcMeshPSCF_plid layout'
   call h5pset_chunk_f   (abcMeshPSCF_plid,1,abcMaxDimChunkPSCF,hdferr)
   if (hdferr /= 0) stop 'Failed to set pscf abc chunk size'
   call h5pset_deflate_f (abcMeshPSCF_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set pscf abcMesh plid for deflation'


   ! Create the datasets that will be used.

   ! Make the mesh data sets.
   call h5dcreate_f (meshPSCF_gid,"meshX",H5T_NATIVE_DOUBLE,&
         & abcMeshPSCF_dsid(1),meshPSCF_did(1),hdferr,abcMeshPSCF_plid)
   if (hdferr /= 0) stop 'Failed to create pscf meshX did'
   call h5dcreate_f (meshPSCF_gid,"meshY",H5T_NATIVE_DOUBLE,&
         & abcMeshPSCF_dsid(2),meshPSCF_did(2),hdferr,abcMeshPSCF_plid)
   if (hdferr /= 0) stop 'Failed to create pscf meshY did'
   call h5dcreate_f (meshPSCF_gid,"meshZ",H5T_NATIVE_DOUBLE,&
         & abcMeshPSCF_dsid(3),meshPSCF_did(3),hdferr,abcMeshPSCF_plid)
   if (hdferr /= 0) stop 'Failed to create pscf meshZ did'

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
   call h5dcreate_f (wavPSCF_gid,"wav_live_up+dn",H5T_NATIVE_DOUBLE,&
         & abcPSCF_dsid,wavPSCF_did(1),hdferr,abcPSCF_plid)
   if (hdferr /= 0) stop 'Failed to create pscf live spin sum wav did'
   call h5dcreate_f (wavPSCF_gid,"wav_live_up-dn",H5T_NATIVE_DOUBLE,&
         & abcPSCF_dsid,wavPSCF_did(2),hdferr,abcPSCF_plid)
   if (hdferr /= 0) stop 'Failed to create pscf live spin diff wav did'
   call h5dcreate_f (wavPSCF_gid,"wav_diff_up+dn",H5T_NATIVE_DOUBLE,&
         & abcPSCF_dsid,wavPSCF_did(3),hdferr,abcPSCF_plid)
   if (hdferr /= 0) stop &
         & 'Failed to create pscf psi^2 spin sum wav difference did'
   call h5dcreate_f (wavPSCF_gid,"wav_diff_up-dn",H5T_NATIVE_DOUBLE,&
         & abcPSCF_dsid,wavPSCF_did(4),hdferr,abcPSCF_plid)
   if (hdferr /= 0) stop &
         & 'Failed to create pscf psi^2 spin diff wav difference did'

   ! Now we will make the charge density (rho) datasets. This includes the
   !   sum and diff spin charge density obtained from the interacting (live)
   !   material and also sum and diff spin data constructed from the difference
   !   between the interacting material and that which would be created from
   !   isolated (i.e. non-interacting) atoms. Further, as with the wave
   !   function above, if the calculation is non spin-polarized then the total
   !   (up+down) will be stored in the "up+dn" portion.
   call h5dcreate_f (rhoPSCF_gid,"rho_live_up+dn",H5T_NATIVE_DOUBLE,&
         & abcPSCF_dsid,rhoPSCF_did(1),hdferr,abcPSCF_plid)
   if (hdferr /= 0) stop &
         & 'Failed to create pscf interacting (live) spin sum rho did'
   call h5dcreate_f (rhoPSCF_gid,"rho_live_up-dn",H5T_NATIVE_DOUBLE,&
         & abcPSCF_dsid,rhoPSCF_did(2),hdferr,abcPSCF_plid)
   if (hdferr /= 0) stop &
         & 'Failed to create pscf interacting (live) spin diff rho did'
   call h5dcreate_f (rhoPSCF_gid,"rho_diff_up+dn",H5T_NATIVE_DOUBLE,&
         & abcPSCF_dsid,rhoPSCF_did(3),hdferr,abcPSCF_plid)
   if (hdferr /= 0) stop 'Failed to create pscf interacting-isolated ' // &
         & 'spin sum rho difference did'
   call h5dcreate_f (rhoPSCF_gid,"rho_diff_up-dn",H5T_NATIVE_DOUBLE,&
         & abcPSCF_dsid,rhoPSCF_did(4),hdferr,abcPSCF_plid)
   if (hdferr /= 0) stop 'Failed to create pscf interacting-isolated ' // &
         & 'spin diff rho difference did'

   ! Finally, we make the potential datasets. As with the charge density, the
   !   sum, diff, interacting (live), and interacting-isolated difference
   !   atomic potential function is obtained.
   call h5dcreate_f (potPSCF_gid,"pot_live_up+dn",H5T_NATIVE_DOUBLE,&
         & abcPSCF_dsid,potPSCF_did(1),hdferr,abcPSCF_plid)
   if (hdferr /= 0) stop 'Failed to create pscf interacting (live) spin ' // &
         & 'sum pot did'
   call h5dcreate_f (potPSCF_gid,"pot_live_up-dn",H5T_NATIVE_DOUBLE,&
         & abcPSCF_dsid,potPSCF_did(2),hdferr,abcPSCF_plid)
   if (hdferr /= 0) stop 'Failed to create pscf interacting (live) spin ' // &
         & 'diff pot did'
   call h5dcreate_f (potPSCF_gid,"pot_diff_up+dn",H5T_NATIVE_DOUBLE,&
         & abcPSCF_dsid,potPSCF_did(3),hdferr,abcPSCF_plid)
   if (hdferr /= 0) stop 'Failed to create pscf interacting-isolated ' // &
         & 'spin sum pot difference did'
   call h5dcreate_f (potPSCF_gid,"pot_diff_up-dn",H5T_NATIVE_DOUBLE,&
         & abcPSCF_dsid,potPSCF_did(4),hdferr,abcPSCF_plid)
   if (hdferr /= 0) stop 'Failed to create pscf interacting-isolated ' // &
         & 'spin diff pot difference did'

end subroutine initPSCFFieldHDF5


subroutine accessPSCFFieldHDF5(pscf_fid)

   ! Use necessary modules.
   use O_Lattice, only: numMeshPoints

   ! Use the HDF5 module.
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed parameters.
   integer(hid_t) :: pscf_fid

   ! Declare local variables.
   integer :: hdferr

!   ! Create the property list for the field hdf5 file and turn off
!   !   chunk caching.
!   call h5pcreate_f    (H5P_FILE_ACCESS_F,fieldPSCF_plid,hdferr)
!   if (hdferr /= 0) stop 'Failed to create field plid in accessFieldHDF5.'
!   call h5pget_cache_f (fieldPSCF_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
!         & hdferr)
!   if (hdferr /= 0) stop 'Failed to get field plid cache settings.'
!   call h5pset_cache_f (fieldPSCF_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
!   if (hdferr /= 0) stop 'Failed to set field plid cache settings.'
!
!
!   ! Open the HDF5 file that will hold all the computed results.  This will
!   !   die if the file already exists.  It opens read-only.
!   call h5fopen_f ("field-temp.hdf5",H5F_ACC_RDONLY_F,pscf_fid,hdferr,&
!         & fieldPSCF_plid)
!   if (hdferr /= 0) stop 'Failed to open field-temp.hdf5 file.'

   ! Initialize data structure dimensions.
   abcDims(:) = numMeshPoints(:)

   ! Open the groups of the hdf5 field file.
   call h5gopen_f (pscf_fid,"/wavGroup",wavPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open pscf wavGroup'
   call h5gopen_f (pscf_fid,"/rhoGroup",rhoPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open pscf rhoGroup'
   call h5gopen_f (pscf_fid,"/potGroup",potPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open pscf potGroup'

   ! Open the mesh group.
   call h5gopen_f (pscf_fid,"/mesh",meshPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open the mesh Group'

   ! Open the datasets.

   ! Wave function |psi|^2 datasets first.
   call h5dopen_f (wavPSCF_gid,"wav_live_up+dn",wavPSCF_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to open pscf live spin sum wav did'
   call h5dopen_f (wavPSCF_gid,"wav_live_up-dn",wavPSCF_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to open pscf live spin diff wav did'
   call h5dopen_f (wavPSCF_gid,"wav_diff_up+dn",wavPSCF_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to open pscf live-iso spin sum wav did'
   call h5dopen_f (wavPSCF_gid,"wav_diff_up-dn",wavPSCF_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to open pscf live-iso spin diff wav did'

   ! Charge density (rho) datasets second.
   call h5dopen_f (rhoPSCF_gid,"rho_live_up+dn",rhoPSCF_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to open pscf live spin sum rho did'
   call h5dopen_f (rhoPSCF_gid,"rho_live_up-dn",rhoPSCF_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to open pscf live spin diff rho did'
   call h5dopen_f (rhoPSCF_gid,"rho_diff_up+dn",rhoPSCF_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to open pscf live-iso spin sum rho did'
   call h5dopen_f (rhoPSCF_gid,"rho_diff_up-dn",rhoPSCF_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to open pscf live-iso spin diff rho did'

   ! Potential function datasets last.
   call h5dopen_f (potPSCF_gid,"pot_live_up+dn",potPSCF_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to open pscf live spin sum pot did'
   call h5dopen_f (potPSCF_gid,"pot_live_up-dn",potPSCF_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to open pscf live spin diff pot did'
   call h5dopen_f (potPSCF_gid,"pot_diff_up+dn",potPSCF_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to open pscf live-iso spin sum pot did'
   call h5dopen_f (potPSCF_gid,"pot_diff_up-dn",potPSCF_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to open pscf live-iso spin diff pot did'

   ! Mesh datasets.
   call h5dopen_f (meshPSCF_gid,"meshX",meshPSCF_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to open meshX did'
   call h5dopen_f (meshPSCF_gid,"meshY",meshPSCF_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to open meshY did'
   call h5dopen_f (meshPSCF_gid,"meshZ",meshPSCF_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to open meshZ did'

   ! Obtain the properties of the datasets that were just opened. They are all
   !   the same and so only one copy is necessary. (Actually, this value is
   !   not really used, but in the "close" subroutine we close this id so we
   !   should make sure to have it open before attempting to close it.
   call h5dget_create_plist_f (wavPSCF_did(1),abcPSCF_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain pscf abcPSCF_plid'

   ! Obtain the dataspace that is used for each dataset. The same dataspace
   !   definition works for all of the datasets.
   call h5dget_space_f (wavPSCF_did(1),abcPSCF_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain pscf abcPSCF_dsid'

   ! Similarly, obtain the plid and dsid for the mesh.
   call h5dget_create_plist_f (meshPSCF_did(1),abcMeshPSCF_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain pscf abcMeshPSCF_plid'
   call h5dget_space_f (meshPSCF_did(1),abcMeshPSCF_dsid(1),hdferr)
   if (hdferr /= 0) stop 'Failed to obtain pscf abcMeshPSCF_dsid X'
   call h5dget_space_f (meshPSCF_did(2),abcMeshPSCF_dsid(2),hdferr)
   if (hdferr /= 0) stop 'Failed to obtain pscf abcMeshPSCF_dsid Y'
   call h5dget_space_f (meshPSCF_did(3),abcMeshPSCF_dsid(3),hdferr)
   if (hdferr /= 0) stop 'Failed to obtain pscf abcMeshPSCF_dsid Z'

   ! Get the size of the chunks from the property list.
   call h5pget_chunk_f (abcPSCF_plid, 3, abcDimsChunkPSCF, hdferr)

end subroutine accessPSCFFieldHDF5


subroutine closePSCFFieldHDF5

   ! Use the HDF5 module.
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare local variables.
   integer :: hdferr

   ! Close the property list.
   call h5pclose_f (abcPSCF_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf abcPSCF_plid'

   ! Close the datasets next.

   ! Close the wave function datasets first.
   call h5dclose_f (wavPSCF_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf wavPSCF_did(1)'
   call h5dclose_f (wavPSCF_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf wavPSCF_did(2)'
   call h5dclose_f (wavPSCF_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf wavPSCF_did(3)'
   call h5dclose_f (wavPSCF_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf wavPSCF_did(4)'

   ! Close the charge density (rho) datasets second.
   call h5dclose_f (rhoPSCF_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf rhoPSCF_did(1)'
   call h5dclose_f (rhoPSCF_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf rhoPSCF_did(2)'
   call h5dclose_f (rhoPSCF_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf rhoPSCF_did(3)'
   call h5dclose_f (rhoPSCF_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf rhoPSCF_did(4)'

   ! Close the potential function datasets first.
   call h5dclose_f (potPSCF_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf potPSCF_did(1)'
   call h5dclose_f (potPSCF_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf potPSCF_did(2)'
   call h5dclose_f (potPSCF_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf potPSCF_did(3)'
   call h5dclose_f (potPSCF_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf potPSCF_did(4)'

   ! Close the mesh datasets.
   call h5dclose_f (meshPSCF_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to close meshPSCF X'
   call h5dclose_f (meshPSCF_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to close meshPSCF Y'
   call h5dclose_f (meshPSCF_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to close meshPSCF Z'


   ! Close the dataspaces next.
   call h5sclose_f (abcPSCF_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf abcPSCF_dsid'
   call h5sclose_f (abcMeshPSCF_dsid(1),hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf abcMeshPSCF_dsid X'
   call h5sclose_f (abcMeshPSCF_dsid(2),hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf abcMeshPSCF_dsid Y'
   call h5sclose_f (abcMeshPSCF_dsid(3),hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf abcMeshPSCF_dsid Z'

   ! Close the groups.
   call h5gclose_f (wavPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf wavPSCF_gid'
   call h5gclose_f (rhoPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf rhoPSCF_gid'
   call h5gclose_f (potPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf potPSCF_gid'
   call h5gclose_f (meshPSCF_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf meshPSCF_gid'

!   ! Close the field property list.
!   call h5pclose_f (fieldPSCF_plid,hdferr)
!   if (hdferr /= 0) stop 'Failed to close pscf fieldPSCF_plid.'
!
!   ! Close the field file.
!   call h5fclose_f (pscf_fid,hdferr)
!   if (hdferr /= 0) stop 'Failed to close pscf_fid.'

end subroutine closePSCFFieldHDF5


end module O_PSCFFieldHDF5
