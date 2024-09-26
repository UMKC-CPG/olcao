module O_PSCFHDF5

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
   integer(hid_t) :: pscf_fid

   ! Declare the property list for the pscf file and its associated parameters.
   integer(hid_t)  :: pscf_plid
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

subroutine initHDF5_PSCF (numStates)

   ! Use necessary modules.
   use O_TimeStamps

   ! Use the HDF5 module.
   use HDF5

   ! Use the subsection object modules for pscf.
   use O_PSCFIntegralsHDF5
   use O_PSCFEigValHDF5
   use O_PSCFEigVecHDF5
   use O_CommandLine, only: excitedQN_n, excitedQN_l, basisCode_PSCF

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed dummy variables.
   integer, intent(in) :: numStates

   ! Declare local variables.
   integer :: hdferr
   logical :: file_exists
   character*15 :: fileName
   character*2 :: edge

   ! Log the time we start to setup the PSCF HDF5 files.
   call timeStampStart(31)

   ! Initialize the Fortran 90 HDF5 interface.
   call h5open_f(hdferr)
   if (hdferr < 0) stop 'Failed to open HDF library'

   ! Identify the file name of the hdf5 file we need to create/open.
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
   if (basisCode_PSCF == 1) then
      write(fileName,fmt="(a2,a13)") edge,"_pscf-mb.hdf5"
   elseif (basisCode_PSCF == 2) then
      write(fileName,fmt="(a2,a13)") edge,"_pscf-fb.hdf5"
   elseif (basisCode_PSCF == 3) then
      write(fileName,fmt="(a2,a13)") edge,"_pscf-eb.hdf5"
   endif

   ! Create the property list for the pscf hdf5 file and turn off
   !   chunk caching.
   call h5pcreate_f (H5P_FILE_ACCESS_F,pscf_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create pscf plid.'
   call h5pget_cache_f (pscf_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to get pscf plid cache settings.'
   call h5pset_cache_f (pscf_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to set pscf plid cache settings.'

   ! Determine if an HDF5 file already exists for this calculation.
   inquire (file=fileName, exist=file_exists)

   ! If it does, then access the existing file. If not, then create one.
   if (file_exists .eqv. .true.) then
      ! We are continuing a previous calculation.

      ! Open the HDF5 file for reading / writing.
      call h5fopen_f (fileName,H5F_ACC_RDWR_F,pscf_fid,hdferr,pscf_plid)
      if (hdferr /= 0) stop 'Failed to open pscf hdf5 file.'

      ! Access the groups of the HDF5 file.
      call accessPSCFIntegralHDF5 (pscf_fid)
      call accessPSCFEigVecHDF5 (pscf_fid,attribInt_dsid,attribIntDims,numStates)
      call accessPSCFEigValHDF5 (pscf_fid,numStates)

   else
      ! We are starting a new calculation.

      ! Create the HDF5 file that will hold all the computed results. This
      !   uses the default file creation and file access properties.
      call h5fcreate_f (fileName,H5F_ACC_EXCL_F,pscf_fid,hdferr,&
            & H5P_DEFAULT_F,pscf_plid)
      if (hdferr /= 0) stop 'Failed to create pscf hdf5 file.'

      ! All datasets will have an attached attribute logging that the
      !   calculation has successfully completed. (Checkpointing.) Thus,
      !   we need to create the shared attribute dataspace.
      attribIntDims(1) = 1
      call h5screate_simple_f (1,attribIntDims(1),attribInt_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to create the attribInt_dsid'

      ! Create the subgroups of the pscf hdf5 file. This must be done in this
      !   order due to dependencies on potPot_dsid and others.
      call initPSCFIntegralHDF5 (pscf_fid,attribInt_dsid,attribIntDims)
      call initPSCFEigVecHDF5 (pscf_fid,attribInt_dsid,attribIntDims,numStates)
      call initPSCFEigValHDF5 (pscf_fid,numStates)
   endif


   ! Log the time we finish setting up the PSCF HDF5 files.
   call timeStampEnd(31)

end subroutine initHDF5_PSCF


subroutine closeHDF5_PSCF

   ! Use the HDF5 module.
   use HDF5

   ! Use the subsection object modules for pscf.
   use O_PSCFIntegralsHDF5, only: closePSCFIntegralHDF5
   use O_PSCFEigVecHDF5,  only: closePSCFEigVecHDF5
   use O_PSCFEigValHDF5,  only: closePSCFEigValHDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare local variables.
   integer :: hdferr

   ! Close access to all the subgroups and their parts.
   call closePSCFIntegralHDF5
   call closePSCFEigVecHDF5
   call closePSCFEigValHDF5

   ! Close the property list.
   call h5pclose_f (pscf_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf_plid.'

   ! Close the file.
   call h5fclose_f (pscf_fid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf_fid.'

   ! Close access to the HDF5 interface.
   call h5close_f (hdferr)
   if (hdferr /= 0) stop 'Failed to close the HDF5 interface.'

end subroutine closeHDF5_PSCF


end module O_PSCFHDF5
