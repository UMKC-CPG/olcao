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

subroutine initHDF5_SCF (maxNumRayPoints, numStates)

   ! Use necessary modules.
   use O_TimeStamps

   ! Use the HDF5 module.
   use HDF5

   ! Use the subsection object modules for scf.
   use O_SCFIntegralsHDF5
   use O_SCFElecStatHDF5
   use O_SCFExchCorrHDF5
   use O_SCFEigValHDF5
   use O_SCFEigVecHDF5
   use O_SCFPotRhoHDF5
   use O_CommandLine, only: excitedQN_n, excitedQN_l, basisCode_SCF, &
         & doSYBD_SCF, doMTOP_SCF
   use O_KPoints, only: numAxialKPoints, numPathKP, numPaths, &
         & numTotalHighSymKP, numHighSymKP, highSymKP, highSymKPChar

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed dummy variables.
   integer, intent(in) :: maxNumRayPoints
   integer, intent(in) :: numStates

   ! Declare local variables.
   integer :: i,j
   integer :: hdferr
   integer :: charCount
   logical :: fileExists
   logical :: groupExists
   character*14 :: fileName
   character*17 :: kPointName
   character*2 :: edge

   ! Log the time we start to setup the SCF HDF5 files.
   call timeStampStart(6)

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
   if (basisCode_SCF == 1) then
      write(fileName,fmt="(a2,a12)") edge,"_scf-mb.hdf5"
   elseif (basisCode_SCF == 2) then
      write(fileName,fmt="(a2,a12)") edge,"_scf-fb.hdf5"
   elseif (basisCode_SCF == 3) then
      write(fileName,fmt="(a2,a12)") edge,"_scf-eb.hdf5"
   endif


   ! Create the property list for the scf hdf5 file and turn off
   !   chunk caching.
   call h5pcreate_f (H5P_FILE_ACCESS_F,scf_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create scf plid.'
   call h5pget_cache_f (scf_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to get scf plid cache settings.'
   call h5pset_cache_f (scf_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to set scf plid cache settings.'

   ! Create the name for the kPoint dependent data group.
   if (doSYBD_SCF == 1) then
      write(kPointName,fmt="(i4.4)") numPathKP
      charCount = 4
      pathLoop: do i = 1, numPaths
         do j = 1, numHighSymKP(i)
            write(kPointName,fmt="(a,1a)") kPointName,highSymKPChar(j,i)
            charCount = charCount + 1
            if (charCount == 17) then
               exit pathLoop
            endif
         enddo
      enddo pathLoop
   elseif (doMTOP_SCF == 1) then
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
      ! We are continuing a previous calculation.

      ! Open the HDF5 file for reading / writing.
      call h5fopen_f (fileName,H5F_ACC_RDWR_F,scf_fid,hdferr,&
            & scf_plid)
      if (hdferr /= 0) stop 'Failed to open scf hdf5 file.'

      ! Check if a top-level group for the current kPoint set exists. If so,
      !   then we are continuing that kPoint set. Otherwise, we are starting
      !   a new kPoint set and will need to initialize it.
      call h5lexists_f(scf_fid,kPointName,groupExists,hdferr)

      ! If the group exists, then access the kPoint dependent data.
      if (groupExists .eqv. .true.) then

         ! Open the group.
         call h5gopen_f(scf_fid,kPointName,kPoint_gid,hdferr)

         ! Access the groups of the HDF5 file.
         call accessSCFIntegralHDF5 (kPoint_gid)
         call accessSCFEigVecHDF5 (kPoint_gid,attribInt_dsid,attribIntDims,&
               & numStates)
         call accessSCFEigValHDF5 (kPoint_gid,numStates)
      else

         ! All datasets will have an attached attribute logging that the
         !   calculation has successfully completed. (Checkpointing.) Thus,
         !   we need to create the shared attribute dataspace.
         attribIntDims(1) = 1
         call h5screate_simple_f (1,attribIntDims(1),attribInt_dsid,hdferr)
         if (hdferr /= 0) stop 'Failed to create the attribInt_dsid'

         ! Create the kPoint group that will hold all kPoint dependent results.
         call h5gcreate_f(scf_fid,kPointName,kPoint_gid,hdferr)

         ! Create the subgroups of the pscf hdf5 file. This must be done in this
         !   order due to dependencies on potPot_dsid and others.
         call initSCFIntegralHDF5 (kPoint_gid,attribInt_dsid,attribIntDims)
         call initSCFEigVecHDF5 (kPoint_gid,attribInt_dsid,attribIntDims,&
               & numStates)
         call initSCFEigValHDF5 (kPoint_gid,numStates)
      endif

      ! Access the kPoint independent groups of the HDF5 file.
      call accessSCFElecStatHDF5 (scf_fid)
      call accessSCFExchCorrHDF5 (scf_fid)
      call accessSCFPotRhoHDF5 (scf_fid)

   else
      ! We are starting a new calculation.

      ! Create the HDF5 file that will hold all the computed results. This
      !   uses the default file creation and file access properties.
      call h5fcreate_f (fileName,H5F_ACC_EXCL_F,scf_fid,hdferr,&
            & H5P_DEFAULT_F,scf_plid)
      if (hdferr /= 0) stop 'Failed to create scf hdf5 file.'

      ! All datasets will have an attached attribute logging that the
      !   calculation has successfully completed. (Checkpointing.) Thus,
      !   we need to create the shared attribute dataspace.
      attribIntDims(1) = 1
      call h5screate_simple_f (1,attribIntDims(1),attribInt_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to create the attribInt_dsid'

      ! Create the kPoint group that will hold all kPoint dependent results.
      call h5gcreate_f(scf_fid,kPointName,kPoint_gid,hdferr)

      ! Create the subgroups of the scf hdf5 file. This must be done in this
      !   order due to dependencies on potPot_dsid and others.
      call initSCFIntegralHDF5 (kPoint_gid,attribInt_dsid,attribIntDims)
      call initSCFEigVecHDF5 (kPoint_gid,attribInt_dsid,attribIntDims,numStates)
      call initSCFEigValHDF5 (kPoint_gid,numStates)
      call initSCFElecStatHDF5 (scf_fid,attribInt_dsid,attribIntDims)
      call initSCFExchCorrHDF5 (scf_fid,attribInt_dsid,attribIntDims,&
            & maxNumRayPoints)
      call initSCFPotRhoHDF5 (scf_fid)
   endif


   ! Log the time we finish setting up the SCF HDF5 files.
   call timeStampEnd(6)

end subroutine initHDF5_SCF


subroutine closeHDF5_SCF

   ! Use the HDF5 module.
   use HDF5

   ! Use the subsection object modules for scf.
   use O_SCFIntegralsHDF5, only: closeSCFIntegralHDF5
   use O_SCFExchCorrHDF5,  only: closeSCFExchCorrHDF5
   use O_SCFElecStatHDF5,  only: closeSCFElecStatHDF5
   use O_SCFEigVecHDF5,  only: closeSCFEigVecHDF5
   use O_SCFEigValHDF5,  only: closeSCFEigValHDF5
   use O_SCFPotRhoHDF5,  only: closeSCFPotRhoHDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare local variables.
   integer :: hdferr

   ! Close access to all the subgroups and their parts.
   call closeSCFIntegralHDF5
   call closeSCFExchCorrHDF5
   call closeSCFElecStatHDF5
   call closeSCFEigVecHDF5
   call closeSCFEigValHDF5
   call closeSCFPotRhoHDF5

   ! Close the property list.
   call h5pclose_f (scf_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close scf_plid.'

   ! Close the file.
   call h5fclose_f (scf_fid,hdferr)
   if (hdferr /= 0) stop 'Failed to close scf_fid.'

   ! Close access to the HDF5 interface.
   call h5close_f (hdferr)
   if (hdferr /= 0) stop 'Failed to close the HDF5 interface.'

end subroutine closeHDF5_SCF


end module O_SCFHDF5
