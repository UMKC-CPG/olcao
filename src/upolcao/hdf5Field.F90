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

   ! Define the file ID and file name.
   integer(hid_t) :: file_fid
   character*21 :: fileName

   ! Define the property list for the field file and its associated parameters.
   integer(hid_t)  :: file_plid
   integer         :: mdc_nelmts  ! Meta-data cache num elements.
   integer(size_t) :: rdcc_nelmts ! Raw-data chunk cache num elements.
   integer(size_t) :: rdcc_nbytes ! Raw-data chunk cache size in bytes.
   real            :: rdcc_w0     ! Raw-data chunk cache weighting parameter.

   ! Declare shared variables that are used for creating attributes that will
   !   track the completion of all datasets.
   integer(hid_t) :: attribInt_dsid ! Attribute dataspace.
   integer(hsize_t), dimension (1) :: attribIntDims ! Dataspace dimensionality

   ! Define the group IDs of the subgroups for the complex wave function,
   !   square of the wave function, charge density, potential function data,
   !   and the mesh.
   integer(hid_t) :: psiR_gid
   integer(hid_t) :: psiI_gid
   integer(hid_t) :: wav_gid
   integer(hid_t) :: rho_gid
   integer(hid_t) :: pot_gid
   integer(hid_t) :: mesh_gid

   ! The dataspaces of the wav, rho, and pot, datasets are all the same
   !   in all key respects (type, dimension, etc.) and therefore can be given
   !   a static shared dsid definition now for abc dimensional axes. The
   !   mesh is similar.
   integer(hid_t) :: field_dsid
   integer(hid_t) :: mesh_dsid

   ! Similarly, each of the wav, rho, and pot, datasets uses the same
   !   property list. The mesh uses a similar (but different) property list.
   integer(hid_t) :: field_plid
   integer(hid_t) :: mesh_plid

   ! Because the data may be large we don't store it all in memory at the same
   !   time. Thus, we need to write only portions to disk at a time. Therefore,
   !   we need to define a chunk dataspace for the file (Field and Mesh).
   integer (hid_t) :: memFieldChunk_dsid
   integer (hid_t) :: memMeshChunk_dsid

   ! Define the array that holds the dimensions of the wav, rho, and pot
   !   datasets, which are all the same. The mesh dataset has slightly
   !   different dimensions
   integer(hsize_t), dimension (3) :: fieldDims
   integer(hsize_t), dimension (4) :: meshDims

   ! Define the array that holds the dimensions of the wav, rho, and pot
   !   data chunk, which are all the same. The mesh array is similar.
   integer(hsize_t), dimension (3) :: fieldDimsChunk
   integer(hsize_t), dimension (4) :: meshDimsChunk

   ! Presently, the number of datasets under each group is fixed, but as always
   !   it might be a good idea to make it flexible for potential unforseen
   !   future uses. Note: "sum" = spin sum (up+dn) and "diff" = spin
   !   difference (up-dn). Also: "live" = from interacting system and "N" =
   !   from neutral atoms. There is only one dataset for the mesh.
   ! For a spin non-polarized (spin degenerate) calculation, each did array
   !   will hold the following:
   !   (1) The total from the interacting (live) system.
   !   (2) The interacting (live) system total minus the non-interacting
   !       (neutral, N) system total.
   !   (3) The non-interacting (neutral, N) system total.
   ! For a spin polarized (spin non-degenerate) calculation, each did array
   !   will hold the following:
   !   (1) The spin sum (up + down) from the interacting (live) system.
   !   (2) The spin difference (up - down) from the interacting (live) system.
   !   (3) The interacting (live) system spin sum (up + down) minus the non-
   !       interacting (neutral, N) system spin sum (up + down).
   !   (4) The non-interacting (neutral, N) system spin sum (up + down).
   ! Note that we do not yet include the live spin difference minus the
   !   neutral spin difference because presently the neutral system will have
   !   no spin difference. Therefore, subtracting it from the live system will
   !   yield any new result. In the future, if the basis functions and the
   !   atomic calculations lead to clear spin differences, then we may want to
   !   FIX this.
   integer(hid_t), dimension (4) :: psiR_did ! (sum,diff,live-N sum,live-N diff)
   integer(hid_t), dimension (4) :: psiI_did ! (sum,diff,live-N sum,live-N diff)
   integer(hid_t), dimension (4) :: wav_did ! (sum,diff,live-N sum,live-N diff)
   integer(hid_t), dimension (4) :: rho_did ! (sum,diff,live-N sum,live-N diff)
   integer(hid_t), dimension (4) :: pot_did ! (sum,diff,live-N sum,live-N diff)
   integer(hid_t) :: mesh_did

   ! Integer to identify which axis is the one that needs to be tracked for the
   !   determination of when to write data chunks.
   integer :: triggerAxis
   integer :: meshTriggerAxis

   ! Arrays of names and order of datasets (profile columns) and groups.
   integer :: numDataSets
   character*16, dimension(20) :: dataSetNames
   character*9, dimension(5) :: groupNames

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
         write(fileName,fmt="(a2,a18)") edge,"_scf+field-mb.hdf5"
      elseif (basisCode_SCF == 2) then
         write(fileName,fmt="(a2,a18)") edge,"_scf+field-fb.hdf5"
      elseif (basisCode_SCF == 3) then
         write(fileName,fmt="(a2,a18)") edge,"_scf+field-eb.hdf5"
      endif
   else
      if (basisCode_PSCF == 1) then
         write(fileName,fmt="(a2,a19)") edge,"_pscf+field-mb.hdf5"
      elseif (basisCode_PSCF == 2) then
         write(fileName,fmt="(a2,a19)") edge,"_pscf+field-fb.hdf5"
      elseif (basisCode_PSCF == 3) then
         write(fileName,fmt="(a2,a19)") edge,"_pscf+field-eb.hdf5"
      endif
   endif

   ! Define the names of the datasets (profile columns) and groups.
   numDataSets = 20

   dataSetNames(1) = "psi_r_live_up+dn" ! Real part of the wave fn.
   dataSetNames(2) = "psi_r_live_up-dn"
   dataSetNames(3) = "psi_r_diff_up+dn"
   dataSetNames(4) = "psi_r_neutral"
   dataSetNames(5) = "psi_i_live_up+dn" ! Imaginary part of the wave fn.
   dataSetNames(6) = "psi_i_live_up-dn"
   dataSetNames(7) = "psi_i_diff_up+dn"
   dataSetNames(8) = "psi_i_neutral"
   dataSetNames(9) = "wav_live_up+dn"
   dataSetNames(10) = "wav_live_up-dn"
   dataSetNames(11) = "wav_diff_up+dn"
   dataSetNames(12) = "wav_neutral"
   dataSetNames(13) = "rho_live_up+dn"
   dataSetNames(14) = "rho_live_up-dn"
   dataSetNames(15) = "rho_diff_up+dn"
   dataSetNames(16) = "rho_neutral"
   dataSetNames(17) = "pot_live_up+dn"
   dataSetNames(18) = "pot_live_up-dn"
   dataSetNames(19) = "pot_diff_up+dn"
   dataSetNames(20) = "pot_neutral"
!   dataSetNames(1) = "psi_r_live_up+dn" ! Real part of the wave fn.
!   dataSetNames(2) = "psi_r_live_up-dn"
!   dataSetNames(3) = "psi_r_diff_up+dn"
!   dataSetNames(4) = "psi_r_diff_up-dn"
!   dataSetNames(5) = "psi_i_live_up+dn" ! Imaginary part of the wave fn.
!   dataSetNames(6) = "psi_i_live_up-dn"
!   dataSetNames(7) = "psi_i_diff_up+dn"
!   dataSetNames(8) = "psi_i_diff_up-dn"
!   dataSetNames(9) = "wav_live_up+dn"
!   dataSetNames(10) = "wav_live_up-dn"
!   dataSetNames(11) = "wav_diff_up+dn"
!   dataSetNames(12) = "wav_diff_up-dn"
!   dataSetNames(13) = "rho_live_up+dn"
!   dataSetNames(14) = "rho_live_up-dn"
!   dataSetNames(15) = "rho_diff_up+dn"
!   dataSetNames(16) = "rho_diff_up-dn"
!   dataSetNames(17) = "pot_live_up+dn"
!   dataSetNames(18) = "pot_live_up-dn"
!   dataSetNames(19) = "pot_diff_up+dn"
!   dataSetNames(20) = "pot_diff_up-dn"
   groupNames(1) = "psiRGroup"
   groupNames(2) = "psiIGroup"
   groupNames(3) = "wavGroup"
   groupNames(4) = "rhoGroup"
   groupNames(5) = "potGroup"

   ! Create the property list for the field hdf5 file and turn off
   !   chunk caching.
   call h5pcreate_f (H5P_FILE_ACCESS_F,file_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field plid.'
   call h5pget_cache_f (file_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to get field plid cache settings.'
   call h5pset_cache_f (file_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to set field plid cache settings.'

   ! Determine if an HDF5 file already exists for this calculation.
   inquire (file=trim(fileName), exist=file_exists)

   ! If it does, then access the existing file. If not, then create one.
   if (file_exists .eqv. .true.) then
      ! We are continuing a previous calculation.

      ! Open the HDF5 file for reading / writing.
      call h5fopen_f (trim(fileName),H5F_ACC_RDWR_F,file_fid,hdferr,file_plid)
      if (hdferr /= 0) stop 'Failed to open field hdf5 file.'

      ! Access the groups of the HDF5 file.
      call accessFieldHDF5 (file_fid)
   else
      ! We are starting a new calculation.

      ! Create the HDF5 file that will hold all the computed results. This
      !   uses the default file creation and file access properties.
      call h5fcreate_f (trim(fileName),H5F_ACC_EXCL_F,file_fid,hdferr,&
            & H5P_DEFAULT_F,file_plid)
      if (hdferr /= 0) stop 'Failed to create field hdf5 file.'

      ! All datasets will have an attached attribute logging that the
      !   calculation has successfully completed. (Checkpointing.) Thus,
      !   we need to create the shared attribute dataspace.
      attribIntDims(1) = 1
      call h5screate_simple_f (1,attribIntDims(1),attribInt_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to create the attribInt_dsid'

      ! Create the subgroups of the field hdf5 file.
      call initFieldHDF5 (file_fid)
   endif

end subroutine prepFieldHDF5


subroutine initFieldHDF5 (file_fid)

   ! Use necessary modules.
   use O_Lattice, only: numMeshPoints

   ! Use the HDF5 module.
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed parameters.
   integer(hid_t) :: file_fid

   ! Declare local variables.
   integer :: maxNumDataPoints
   integer :: hdferr

!   ! Initialize the Fortran 90 HDF5 interface.
!   call h5open_f(hdferr)
!   if (hdferr < 0) stop 'Failed to open HDF library'
!
!   ! Create the property list for the field hdf5 file and turn off
!   !   chunk caching.
!   call h5pcreate_f    (H5P_FILE_ACCESS_F,file_plid,hdferr)
!   if (hdferr /= 0) stop 'Failed to create field plid.'
!   call h5pget_cache_f (file_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
!         & hdferr)
!   if (hdferr /= 0) stop 'Failed to get field plid cache settings.'
!   call h5pset_cache_f (file_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
!   if (hdferr /= 0) stop 'Failed to set field plid cache settings.'

!   ! Create the HDF5 file that will hold all the computed results.  This will
!   !   die if the file already exists.  It also uses the default file creation
!   !   and file access properties.
!   call h5fcreate_f ("field-temp.hdf5",H5F_ACC_EXCL_F,file_fid,hdferr,&
!         & H5P_DEFAULT_F,file_plid)
!   if (hdferr /= 0) stop 'Failed to create field-temp.hdf5 file.'

   ! Create the groups of the field hdf5 file.
   call h5gcreate_f (file_fid,groupNames(1),psiR_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field complex wave fn (psiR) group'
   call h5gcreate_f (file_fid,groupNames(2),psiI_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field complex wave fn (psiI) group'
   call h5gcreate_f (file_fid,groupNames(3),wav_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field psi^2 (wav) group'
   call h5gcreate_f (file_fid,groupNames(4),rho_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field charge density (rho) group'
   call h5gcreate_f (file_fid,groupNames(5),pot_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field potential (pot) group'
   call h5gcreate_f (file_fid,"meshGroup",mesh_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field mesh group'

   ! Initialize data structure dimensions for on-disk representation and the
   !   one in-memory (chunk) representation that we know for sure so far.
   fieldDims(:) = numMeshPoints(:)
   meshDims(1) = 3
   meshDims(2:4) = fieldDims(1:3)
   meshDimsChunk(1) = 3

   ! Establish the upper limit for data chunk size. We start with an assumption
   !   that the numbers being stored are 8 byte reals and that we should not go
   !   over 4 billion bytes in a given chunk. Thus, for the three dimensional
   !   data array we require a*b*c < 500M. (From 4e9 / 8.) Without thinking
   !   anything through we will just demand dimensions that meet that
   !   requirement through a simple scale factor. Note that this may not work
   !   correctly for sufficiently large fieldDims(:) values because of integer
   !   overflow. Need to write this in a way that is immune to integer overflow
   !   problems.
   ! Compute the largest chunk size that contains fewer than the hard-coded
   !   maximum number of datapoints yet that still retains the largest extent
   !   along the a, b, and c axes as possible in that priority order. The
   !   number of data points depends on the number of fields that will be
   !   simultaneously computed.
   maxNumDataPoints = 500000000

   ! Check that the minimum space requirement is met to store data easily
   if (fieldDims(1) > maxNumDataPoints) then  ! Note the greater-than sign.
      stop 'Increase maxNumDataPoints'
   else
      ! Assume that the chunk dimensions will be limited by the a-axis.
      fieldDimsChunk(1) = fieldDims(1)
      fieldDimsChunk(2:3) = 1
      triggerAxis = 1

      ! Check if we can store two dimensions at a time.
      if (fieldDims(1)*fieldDims(2) < maxNumDataPoints) then ! Note less-than.
         fieldDimsChunk(1:2) = fieldDims(1:2)
         fieldDimsChunk(3) = 1
         triggerAxis = 2

         ! Check if we can store all three dimensions.
         if (fieldDims(1)*fieldDims(2)*fieldDims(3) < &
               & maxNumDataPoints) then ! Note less-than.
            fieldDimsChunk(1:3) = fieldDims(1:3)
            triggerAxis = 3
         endif
      endif
   endif

   ! Repeat the process for the mesh with a factor of 3 included to keep
   !   the chunk size limited.
   if (meshDims(2)*3 > maxNumDataPoints) then
      stop 'Increase maxNumDataPoints'
   else
      meshDimsChunk(2) = meshDims(2)
      meshDimsChunk(3:4) = 1
      meshTriggerAxis = 1

      if (meshDims(2)*meshDims(3)*3 < maxNumDataPoints) then
         meshDimsChunk(2:3) = meshDims(2:3)
         meshDimsChunk(4) = 1
         meshTriggerAxis = 2

         if (meshDims(1)*meshDims(2)*meshDims(3)*3 < maxNumDataPoints) then
            meshDimsChunk(2:4) = meshDims(2:4)
            meshTriggerAxis = 3
         endif
      endif
   endif

   ! Create the dataspace that will be used for each dataset within the file.
   !   The same dataspace definition works for all of the datasets.
   call h5screate_simple_f(3,fieldDims,field_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field_dsid'
   call h5screate_simple_f(4,meshDims,mesh_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create mesh_dsid'

   ! Create the dataspace that will be used to describe the data in memory
   !   before being written to a file.
   call h5screate_simple_f(3,fieldDimsChunk,memFieldChunk_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field memFieldChunk_dsid'
   call h5screate_simple_f(4,meshDimsChunk,memMeshChunk_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field memMeshChunk_dsid'

   ! Define the properties of the datasets to be made.

   ! Create the field property list, then set the properties one at a time.
   call h5pcreate_f      (H5P_DATASET_CREATE_F,field_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field field_plid'
   call h5pset_layout_f  (field_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set field field_plid layout'
   call h5pset_chunk_f   (field_plid,3,fieldDimsChunk,hdferr)
   if (hdferr /= 0) stop 'Failed to set field abc chunk size'
   call h5pset_deflate_f (field_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set field abc plid for deflation'

   ! Create the mesh property list, then set the properties one at a time.
   call h5pcreate_f      (H5P_DATASET_CREATE_F,mesh_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create field mesh_plid'
   call h5pset_layout_f  (mesh_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set field mesh_plid layout'
   call h5pset_chunk_f   (mesh_plid,4,meshDimsChunk,hdferr)
   if (hdferr /= 0) stop 'Failed to set field mesh chunk size'
   call h5pset_deflate_f (mesh_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set field mesh plid for deflation'

   ! Create the datasets that will be used.

   ! First we make all the wave function datasets. This includes the real and
   !   imaginary spin sum and spin diff data. Note that when a spin non-
   !   polarized calculation is done then the total (up+down) will be stored in
   !   the "sum" dataset.

   ! Real part first. (1-4)
   call h5dcreate_f (psiR_gid,dataSetNames(1),H5T_NATIVE_DOUBLE,field_dsid,&
         & psiR_did(1),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field real live spin sum psi did'
   call h5dcreate_f (psiR_gid,dataSetNames(2),H5T_NATIVE_DOUBLE,field_dsid,&
         & psiR_did(2),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field real live spin diff psi did'
   call h5dcreate_f (psiR_gid,dataSetNames(3),H5T_NATIVE_DOUBLE,field_dsid,&
         & psiR_did(3),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field real interacting' // &
         & '-isolated spin sum psi difference did'
   call h5dcreate_f (psiR_gid,dataSetNames(4),H5T_NATIVE_DOUBLE,field_dsid,&
         & psiR_did(4),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field real neutral psi did.'

   ! Imaginary part second. (5-8)
   call h5dcreate_f (psiI_gid,dataSetNames(5),H5T_NATIVE_DOUBLE,field_dsid,&
         & psiI_did(1),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field imag live spin sum psi did'
   call h5dcreate_f (psiI_gid,dataSetNames(6),H5T_NATIVE_DOUBLE,field_dsid,&
         & psiI_did(2),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field imag live spin diff psi did'
   call h5dcreate_f (psiI_gid,dataSetNames(7),H5T_NATIVE_DOUBLE,field_dsid,&
         & psiI_did(3),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field imag interacting' // &
         & '-isolated spin sum psi difference did'
   call h5dcreate_f (psiI_gid,dataSetNames(8),H5T_NATIVE_DOUBLE,field_dsid,&
         & psiI_did(4),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field imag neutral psi did.'

   ! Now we make all the wave function squared datasets. This includes the
   !   sum and difference spin psi^2 obtained from the interacting (live)
   !   material and also sum and diff spin psi^2 constructed from the difference
   !   between the interacting material and that which would be created from
   !   isolated (i.e. non-interacting) atoms. Note that when a spin non-
   !   polarized calculation is done then the total (up+down) will be stored in
   !   the "sum" dataset.
   call h5dcreate_f (wav_gid,dataSetNames(9),H5T_NATIVE_DOUBLE,field_dsid,&
         & wav_did(1),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field live spin sum wav did'
   call h5dcreate_f (wav_gid,dataSetNames(10),H5T_NATIVE_DOUBLE,field_dsid,&
         & wav_did(2),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field live spin diff wav did'
   call h5dcreate_f (wav_gid,dataSetNames(11),H5T_NATIVE_DOUBLE,field_dsid,&
         & wav_did(3),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field interacting-isolated ' // &
         & 'spin sum wav difference did'
   call h5dcreate_f (wav_gid,dataSetNames(12),H5T_NATIVE_DOUBLE,field_dsid,&
         & wav_did(4),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field wav neutral did.'

   ! Now we will make the charge density (rho) datasets. This includes the
   !   sum and diff spin charge density obtained from the interacting (live)
   !   material and also sum and diff spin data constructed from the difference
   !   between the interacting material and that which would be created from
   !   isolated (i.e. non-interacting) atoms. Further, as with the wave
   !   function above, if the calculation is non spin-polarized then the total
   !   (up+down) will be stored in the "up+dn" portion.
   call h5dcreate_f (rho_gid,dataSetNames(13),H5T_NATIVE_DOUBLE,field_dsid,&
         & rho_did(1),hdferr,field_plid)
   if (hdferr /= 0) stop &
         & 'Failed to create field interacting (live) spin sum rho did'
   call h5dcreate_f (rho_gid,dataSetNames(14),H5T_NATIVE_DOUBLE,field_dsid,&
         & rho_did(2),hdferr,field_plid)
   if (hdferr /= 0) stop &
         & 'Failed to create field interacting (live) spin diff rho did'
   call h5dcreate_f (rho_gid,dataSetNames(15),H5T_NATIVE_DOUBLE,field_dsid,&
         & rho_did(3),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field interacting-isolated ' // &
         & 'spin sum rho difference did'
   call h5dcreate_f (rho_gid,dataSetNames(16),H5T_NATIVE_DOUBLE,field_dsid,&
         & rho_did(4),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field rho neutral did.'

   ! Then, we make the potential datasets. As with the charge density, the
   !   sum, diff, interacting (live), and interacting-isolated difference
   !   atomic potential function is obtained.
   call h5dcreate_f (pot_gid,dataSetNames(17),H5T_NATIVE_DOUBLE,field_dsid,&
         & pot_did(1),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field interacting (live) spin ' // &
         & 'sum pot did'
   call h5dcreate_f (pot_gid,dataSetNames(18),H5T_NATIVE_DOUBLE,field_dsid,&
         & pot_did(2),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field interacting (live) spin ' // &
         & 'diff pot did'
   call h5dcreate_f (pot_gid,dataSetNames(19),H5T_NATIVE_DOUBLE,field_dsid,&
         & pot_did(3),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field interacting-isolated ' // &
         & 'spin sum pot difference did'
   call h5dcreate_f (pot_gid,dataSetNames(20),H5T_NATIVE_DOUBLE,field_dsid,&
         & pot_did(4),hdferr,field_plid)
   if (hdferr /= 0) stop 'Failed to create field pot neutral did.'

   ! Finally, we make the mesh dataset.
   call h5dcreate_f (mesh_gid,"mesh",H5T_NATIVE_DOUBLE,mesh_dsid,mesh_did,&
         & hdferr,mesh_plid)
   if (hdferr /= 0) stop 'Failed to create the mesh did'

end subroutine initFieldHDF5


subroutine accessFieldHDF5(file_fid)

   ! Use necessary modules.
   use O_Lattice, only: numMeshPoints

   ! Use the HDF5 module.
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed parameters.
   integer(hid_t) :: file_fid

   ! Declare local variables.
   integer :: hdferr
   integer :: maxNumDataPoints

!   ! Create the property list for the field hdf5 file and turn off
!   !   chunk caching.
!   call h5pcreate_f    (H5P_FILE_ACCESS_F,file_plid,hdferr)
!   if (hdferr /= 0) stop 'Failed to create field plid in accessFieldHDF5.'
!   call h5pget_cache_f (file_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
!         & hdferr)
!   if (hdferr /= 0) stop 'Failed to get field plid cache settings.'
!   call h5pset_cache_f (file_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
!   if (hdferr /= 0) stop 'Failed to set field plid cache settings.'
!
!
!   ! Open the HDF5 file that will hold all the computed results.  This will
!   !   die if the file already exists.  It opens read-only.
!   call h5fopen_f ("field-temp.hdf5",H5F_ACC_RDONLY_F,file_fid,hdferr,&
!         & file_plid)
!   if (hdferr /= 0) stop 'Failed to open field-temp.hdf5 file.'

   ! Initialize data structure dimensions for on-disk representation and the
   !   one in-memory (chunk) representation that we know for sure so far.
   fieldDims(:) = numMeshPoints(:)
   meshDims(1) = 3
   meshDims(2:4) = fieldDims(1:3)
   meshDimsChunk(1) = 3

   ! Establish the upper limit for data chunk size. We start with an assumption
   !   that the numbers being stored are 8 byte reals and that we should not go
   !   over 2 billion bytes in a given chunk. Thus, for the three dimensional
   !   data array we require a*b*c < 250M. (From 2e9 / 8.) Without thinking
   !   anything through we will just demand dimensions that meet that
   !   requirement through a simple scale factor. Note that this may not work
   !   correctly for sufficiently large fieldDims(:) values because of integer
   !   overflow. Need to write this is a way that is immune to integer overflow
   !   problems.
   ! Compute the largest chunk size that contains fewer than the hard-coded
   !   maximum number of datapoints yet that still retains the largest extent
   !   along the a, b, and c axes as possible in that priority order.
   maxNumDataPoints = 250000000

   ! Check that the minimum space requirement is met to store data easily
   if (fieldDims(1) > maxNumDataPoints) then  ! Note the greater-than sign.
      stop 'Increase maxNumDataPoints'
   else
      ! Assume that the chunk dimensions will be limited by the a-axis.
      fieldDimsChunk(1) = fieldDims(1)
      fieldDimsChunk(2:3) = 1
      triggerAxis = 1

      ! Check if we can store two dimensions at a time.
      if (fieldDims(1)*fieldDims(2) < maxNumDataPoints) then ! Note less-than.
         fieldDimsChunk(1:2) = fieldDims(1:2)
         fieldDimsChunk(3) = 1
         triggerAxis = 2

         ! Check if we can store all three dimensions.
         if (fieldDims(1)*fieldDims(2)*fieldDims(3) < &
               & maxNumDataPoints) then ! Note <
            fieldDimsChunk(1:3) = fieldDims(1:3)
            triggerAxis = 3
         endif
      endif
   endif

   if (meshDims(2)*3 > maxNumDataPoints) then
      stop 'Increase maxNumDataPoints'
   else
      meshDimsChunk(2) = meshDims(2)
      meshDimsChunk(3:4) = 1
      meshTriggerAxis = 1

      if (meshDims(2)*meshDims(3)*3 < maxNumDataPoints) then
         meshDimsChunk(2:3) = meshDims(2:3)
         meshDimsChunk(4) = 1
         meshTriggerAxis = 2

         if (meshDims(1)*meshDims(2)*meshDims(3)*3 < maxNumDataPoints) then
            meshDimsChunk(2:4) = meshDims(2:4)
            meshTriggerAxis = 3
         endif
      endif
   endif

   ! Open the groups of the hdf5 field file.
   call h5gopen_f (file_fid,groupNames(1),psiR_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open field psiRGroup'
   call h5gopen_f (file_fid,groupNames(2),psiI_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open field psiIGroup'
   call h5gopen_f (file_fid,groupNames(3),wav_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open field wavGroup'
   call h5gopen_f (file_fid,groupNames(4),rho_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open field rhoGroup'
   call h5gopen_f (file_fid,groupNames(5),pot_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open field potGroup'
   call h5gopen_f (file_fid,"meshGroup",mesh_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open field mesh group'

   ! Open the datasets.

   ! Wave function (complex psi real part) datasets.
   call h5dopen_f (psiR_gid,dataSetNames(1),psiR_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to open real live spin sum psi did'
   call h5dopen_f (psiR_gid,dataSetNames(2),psiR_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to open real live spin diff psi did'
   call h5dopen_f (psiR_gid,dataSetNames(3),psiR_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to open real live-iso spin sum psi did'
   call h5dopen_f (psiR_gid,dataSetNames(4),psiR_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to open real psi neutral did'

   ! Wave function (complex psi imaginary part) datasets.
   call h5dopen_f (psiI_gid,dataSetNames(5),psiI_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to open imag live spin sum psi did'
   call h5dopen_f (psiI_gid,dataSetNames(6),psiI_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to open imag live spin diff psi did'
   call h5dopen_f (psiI_gid,dataSetNames(7),psiI_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to open imag live-iso spin sum psi did'
   call h5dopen_f (psiI_gid,dataSetNames(8),psiI_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to open imag psi neutral did'

   ! Wave function |psi|^2 datasets.
   call h5dopen_f (wav_gid,dataSetNames(9),wav_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to open live spin sum wav did'
   call h5dopen_f (wav_gid,dataSetNames(10),wav_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to open live spin diff wav did'
   call h5dopen_f (wav_gid,dataSetNames(11),wav_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to open live-iso spin sum wav did'
   call h5dopen_f (wav_gid,dataSetNames(12),wav_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to open wav neutral did.'

   ! Charge density (rho) datasets.
   call h5dopen_f (rho_gid,dataSetNames(13),rho_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to open live spin sum rho did'
   call h5dopen_f (rho_gid,dataSetNames(14),rho_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to open live spin diff rho did'
   call h5dopen_f (rho_gid,dataSetNames(15),rho_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to open live-iso spin sum rho did'
   call h5dopen_f (rho_gid,dataSetNames(16),rho_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to open rho neutral did.'

   ! Potential function datasets.
   call h5dopen_f (pot_gid,dataSetNames(17),pot_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to open live spin sum pot did'
   call h5dopen_f (pot_gid,dataSetNames(18),pot_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to open live spin diff pot did'
   call h5dopen_f (pot_gid,dataSetNames(19),pot_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to open live-iso spin sum pot did'
   call h5dopen_f (pot_gid,dataSetNames(20),pot_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to open pot neutral did.'

   ! Mesh dataset last.
   call h5dopen_f (mesh_gid,"mesh",mesh_did,hdferr)
   if (hdferr /= 0) stop 'Failed to open mesh did'

   ! Obtain the properties of the datasets that were just opened. They are all
   !   the same and so only one copy is necessary. (Actually, this value is
   !   not really used, but in the "close" subroutine we close this id so we
   !   should make sure to have it open before attempting to close it.
   call h5dget_create_plist_f (wav_did(1),field_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain field_plid'
   call h5dget_create_plist_f (mesh_did,mesh_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain mesh_plid'

   ! Obtain the dataspace that is used for the all field datasets.
   call h5dget_space_f (wav_did(1),field_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain field_dsid'

   ! Obtain the dataspace that is used for the mesh dataset.
   call h5dget_space_f (mesh_did,mesh_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain mesh_dsid'

end subroutine accessFieldHDF5


subroutine closeFieldHDF5

   ! Use the HDF5 module.
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare local variables.
   integer :: hdferr

   ! Close the property list.
   call h5pclose_f (field_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close field_plid'
   call h5pclose_f (mesh_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close mesh_plid'

   ! Close the datasets next.

   ! Close the complex wave function real part datasets.
   call h5dclose_f (psiR_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to close field psiR_did(1)'
   call h5dclose_f (psiR_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to close field psiR_did(2)'
   call h5dclose_f (psiR_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to close field psiR_did(3)'
   call h5dclose_f (psiR_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to close field psiR_did(4)'

   ! Close the complex wave function imaginary part datasets.
   call h5dclose_f (psiI_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to close field psiI_did(1)'
   call h5dclose_f (psiI_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to close field psiI_did(2)'
   call h5dclose_f (psiI_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to close field psiI_did(3)'
   call h5dclose_f (psiI_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to close field psiI_did(4)'

   ! Close the wave function squared datasets.
   call h5dclose_f (wav_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to close field wav_did(1)'
   call h5dclose_f (wav_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to close field wav_did(2)'
   call h5dclose_f (wav_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to close field wav_did(3)'
   call h5dclose_f (wav_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to close field wav_did(4)'

   ! Close the charge density (rho) datasets.
   call h5dclose_f (rho_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to close field rho_did(1)'
   call h5dclose_f (rho_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to close field rho_did(2)'
   call h5dclose_f (rho_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to close field rho_did(3)'
   call h5dclose_f (rho_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to close field rho_did(4)'

   ! Close the potential function datasets.
   call h5dclose_f (pot_did(1),hdferr)
   if (hdferr /= 0) stop 'Failed to close field pot_did(1)'
   call h5dclose_f (pot_did(2),hdferr)
   if (hdferr /= 0) stop 'Failed to close field pot_did(2)'
   call h5dclose_f (pot_did(3),hdferr)
   if (hdferr /= 0) stop 'Failed to close field pot_did(3)'
   call h5dclose_f (pot_did(4),hdferr)
   if (hdferr /= 0) stop 'Failed to close field pot_did(4)'

   ! Close the mesh dataset last.
   call h5dclose_f (mesh_did,hdferr)
   if (hdferr /= 0) stop 'Failed to close mesh_did'

   ! Close the dataspace next.
   call h5sclose_f (field_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close field_dsid'
   call h5sclose_f (mesh_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close mesh_dsid'

   ! Close the groups.
   call h5gclose_f (wav_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close field wav_gid'
   call h5gclose_f (rho_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close field rho_gid'
   call h5gclose_f (pot_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close field pot_gid'
   call h5gclose_f (mesh_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close mesh_gid'

   ! Close the file property list.
   call h5pclose_f (file_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close file_plid.'

   ! Close the file HDF ID.
   call h5fclose_f (file_fid,hdferr)
   if (hdferr /= 0) stop 'Failed to close file_fid.'

end subroutine closeFieldHDF5


end module O_FieldHDF5
