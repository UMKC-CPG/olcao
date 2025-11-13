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

   ! Declare the kPoint group ID.
   integer(hid_t) :: kPoint_gid

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
   use O_PSCFFieldHDF5
   use O_CommandLine, only: excitedQN_n, excitedQN_l, basisCode_PSCF, &
         & doSYBD_PSCF, doMTOP_PSCF
   use O_KPoints, only: numAxialKPoints, numPathKP, numPaths, &
         & numTotalHighSymKP, numHighSymKP, highSymKP, highSymKPChar

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed dummy variables.
   integer, intent(in) :: numStates

   ! Declare local variables.
   integer :: i,j
   integer :: hdferr
   integer :: charCount
   logical :: fileExists
   logical :: groupExists
   character*15 :: fileName
   character*17 :: kPointName
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

   ! Create the name for the kPoint dependent data group.
   if (doSYBD_PSCF == 1) then
      write(kPointName,fmt="(i4.4,a1)") numPathKP, "_"
      charCount = 5
      pathLoop: do i = 1, numPaths
         do j = 1, numHighSymKP(i)
            charCount = charCount + 1
            kPointName(charCount:charCount) = highSymKPChar(j,i)
            if (charCount == 17) then
               exit pathLoop
            endif
         enddo
      enddo pathLoop
   elseif (doMTOP_PSCF == 1) then
      write(kPointName,fmt="(i4.4,a1,i4.4,a1,i4.4,a3)") numAxialKPoints(1), &
            & "_", numAxialKPoints(2), "_", numAxialKPoints(3), "_MP"
   else
      write(kPointName,fmt="(i5.5,a1,i5.5,a1,i5.5)") numAxialKPoints(1), &
            & "_", numAxialKPoints(2), "_", numAxialKPoints(3)
   endif

   ! Determine if an HDF5 file already exists for this calculation.
   inquire (file=fileName, exist=fileExists)

   ! If it does, then access the existing file. If not, then create one.
   if (fileExists .eqv. .true.) then
      ! We are continuing a previous calculation, but we might not be
      !   continuing with the same set of kPoints.

      ! Open the HDF5 file for reading / writing.
      call h5fopen_f (fileName,H5F_ACC_RDWR_F,pscf_fid,hdferr,pscf_plid)
      if (hdferr /= 0) stop 'Failed to open pscf hdf5 file.'

      ! Check if a top-level group for the current kPoint set exists. If so,
      !   then we are continuing that kPoint set. Otherwise, we are starting
      !   a new kPoint set and will need to initialize it.
      call h5lexists_f(pscf_fid,kPointName,groupExists,hdferr)

      ! If the group exists, then access the kPoint dependent data.
      if (groupExists .eqv. .true.) then

         ! Open the group.
         call h5gopen_f(pscf_fid,kPointName,kPoint_gid,hdferr)

         ! Access the groups of the HDF5 file.
         call accessPSCFIntegralHDF5 (kPoint_gid)
         call accessPSCFEigVecHDF5 (kPoint_gid,attribInt_dsid,attribIntDims,&
               & numStates)
         call accessPSCFEigValHDF5 (kPoint_gid,numStates)
      else

         ! All datasets will have an attached attribute logging that the
         !   calculation has successfully completed. (Checkpointing.) Thus,
         !   we need to create the shared attribute dataspace.
         attribIntDims(1) = 1
         call h5screate_simple_f (1,attribIntDims(1),attribInt_dsid,hdferr)
         if (hdferr /= 0) stop 'Failed to create the attribInt_dsid PSCF'

         ! Create the kPoint group that will hold all kPoint dependent results.
         call h5gcreate_f(pscf_fid,kPointName,kPoint_gid,hdferr)

         ! Create the subgroups of the pscf hdf5 file. This must be done in this
         !   order due to dependencies on potPot_dsid and others.
         call initPSCFIntegralHDF5 (kPoint_gid,attribInt_dsid,attribIntDims)
         call initPSCFEigVecHDF5 (kPoint_gid,attribInt_dsid,attribIntDims,&
               & numStates)
         call initPSCFEigValHDF5 (kPoint_gid,numStates)
      endif

      ! Access the kPoint independent data
      call accessPSCFFieldHDF5 (pscf_fid)

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
      if (hdferr /= 0) stop 'Failed to create the attribInt_dsid PSCF'

      ! Create the kPoint group that will hold all kPoint dependent results.
      call h5gcreate_f(pscf_fid,kPointName,kPoint_gid,hdferr)

      ! Create the subgroups of the pscf hdf5 file. This must be done in this
      !   order due to dependencies on potPot_dsid and others.
      call initPSCFIntegralHDF5 (kPoint_gid,attribInt_dsid,attribIntDims)
      call initPSCFEigVecHDF5 (kPoint_gid,attribInt_dsid,attribIntDims,&
            & numStates)
      call initPSCFEigValHDF5 (kPoint_gid,numStates)
      call initPSCFFieldHDF5 (pscf_fid)
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
   use O_PSCFFieldHDF5, only: closePSCFFieldHDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare local variables.
   integer :: hdferr

   ! Close access to all the subgroups and their parts.
   call closePSCFIntegralHDF5
   call closePSCFEigVecHDF5
   call closePSCFEigValHDF5
   call closePSCFFieldHDF5

   ! Close the property list.
   call h5pclose_f (pscf_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf_plid.'

   ! Close the file.
   call h5fclose_f (pscf_fid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf_fid.'

   ! Close access to the HDF5 interface.
   call h5close_f (hdferr)
   if (hdferr /= 0) stop 'Failed to close the HDF5 interface PSCF.'

end subroutine closeHDF5_PSCF


end module O_PSCFHDF5
