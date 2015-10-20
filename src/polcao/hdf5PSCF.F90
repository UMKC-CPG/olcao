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

   ! Define the file ID
   integer(hid_t) :: pscf_fid

   ! Define the property list for the pscf file and its associated parameters.
   integer(hid_t)  :: pscf_plid
   integer         :: mdc_nelmts  ! Meta-data cache num elements.
   integer(size_t) :: rdcc_nelmts ! Raw-data chunk cache num elements.
   integer(size_t) :: rdcc_nbytes ! Raw-data chunk cache size in bytes.
   real            :: rdcc_w0     ! Raw-data chunk cache weighting parameter.

   ! Define the group IDs. (Each group will have a dataset for each kpoint.)
   integer(hid_t) :: overlapVV_gid
   integer(hid_t) :: overlapCV_gid
   integer(hid_t) :: momentumX_gid
   integer(hid_t) :: momentumY_gid
   integer(hid_t) :: momentumZ_gid
   integer(hid_t) :: eigenVec_gid
   integer(hid_t) :: eigenVal_gid

   ! The dataspaces of each dataset in the overlapVV and momentum groups
   !   are the same in all characteristics (type, dimension, etc.) and
   !   therefore can be given a shared ID definition. The dataspaces for the
   !   overlapCV, eigenVectors, and eigenValues are different.
   integer(hid_t) :: valeVale_dsid ! overlapVV and Momentum
   integer(hid_t) :: coreVale_dsid ! overlapCV
   integer(hid_t) :: valeStates_dsid ! Eigenvectors
   integer(hid_t) :: states_dsid ! Eigenvalues

   ! This dataspace will be used for attributes that have a single integer
   !   value. (Specifically, the flag for whether or not the momentum matrices
   !   have been computed.)
   integer(hid_t) :: attribInt_dsid

   ! Define an attribute ID that will be used in the momentumX_gid object. It
   !   will say whether it has (1) or has not (0) been computed.
   integer(hid_t) :: momentumX_aid

   ! Each of the overlapVV and momentum datasets for each kpoint will use the
   !   same property list. Each of the overlapCV, eigenVector, and eigenValue
   !   datasets will use another property list.
   integer(hid_t) :: valeVale_plid ! overlapVV and Momentum
   integer(hid_t) :: coreVale_plid ! overlapCV
   integer(hid_t) :: valeStates_plid ! Eigenvector
   integer(hid_t) :: states_plid ! Eigenvalues

   ! Define arrays that hold the dimensions of the datasets.
   integer(hsize_t), dimension (2) :: valeValeDims ! OverlapVV and Momentum
   integer(hsize_t), dimension (2) :: coreValeDims ! OverlapCV
   integer(hsize_t), dimension (2) :: valeStatesDims ! Eigenvectors
   integer(hsize_t), dimension (1) :: statesDims ! Eigenvalues
   integer(hsize_t), dimension (1) :: attribIntDims

   ! Define array that holds the dimensions of the chunk.
   integer (hsize_t), dimension (2) :: valeValeDimsChunk ! OverlapVV & Momentum
   integer (hsize_t), dimension (2) :: coreValeDimsChunk ! OverlapCV
   integer (hsize_t), dimension (2) :: valeStatesDimsChunk ! Eigenvectors

   ! The number of datasets under each of the groups will depend on a few
   !   things. For the overlapCV, the first array index is either 1 or 2 to
   !   allow for storage of real and imaginary values. (Gamma kpoints
   !   calculatoins will only nee 1.) The second index is for the number of
   !   kpoints. For the overlapVV and momentum, the array index is for the
   !   number of kpoints. For the eigenVectors, the first index is either 1
   !   or 2 to allow for storage of real and imaginary values. (Gamma kpoint
   !   calculations will only need 1.) The second index is for the number of
   !   kpoints. The third index is for the spin (either 1 or 2). For the
   !   eigenValues the first index is for the number of kpoints and the second
   !   index is for spin (either 1 or 2). Note that the overlap and
   !   momentum matrices are square and hermitian matrices so that we will be
   !   storing the data in a compacted form (real in the upper triangle and
   !   imaginary in the strict lower triangle).
   integer(hid_t), allocatable, dimension (:,:)   :: overlapCV_did
   integer(hid_t), allocatable, dimension (:)     :: overlapVV_did
   integer(hid_t), allocatable, dimension (:)     :: momentumX_did
   integer(hid_t), allocatable, dimension (:)     :: momentumY_did
   integer(hid_t), allocatable, dimension (:)     :: momentumZ_did
   integer(hid_t), allocatable, dimension (:,:,:) :: eigenVec_did
   integer(hid_t), allocatable, dimension (:,:)   :: eigenVal_did

   ! In the case of reading the band data it is possible that two band data
   !   files will need to be opened at the same time for PACS calculations.
   !   In this case a second set of group, dataset, property list, and file IDs
   !   will be needed. These datasets can share the same dataspace IDs and
   !   matrix dimension definitions.

   ! Second file ID.
   integer(hid_t) :: pscf2_fid

   ! Second property list ID.
   integer(hid_t) :: pscf2_plid

   ! Second group IDs.
   integer(hid_t) :: eigenVec2_gid
   integer(hid_t) :: eigenVal2_gid
   integer(hid_t) :: overlapVV2_gid

   ! Second dataset IDs.
   integer(hid_t), allocatable, dimension (:)     :: overlapVV2_did
   integer(hid_t), allocatable, dimension (:,:,:) :: eigenVec2_did
   integer(hid_t), allocatable, dimension (:,:)   :: eigenVal2_did

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

! This is only called by the pscf program. Other programs that use the data
!   in the pscf hdf5 file will use the access and closeAccess subroutines.
subroutine initPSCFHDF5

   ! Use necessary modules.
   use O_Kinds
   use O_Input,       only: numStates
   use O_Potential,   only: spin
   use O_KPoints,     only: numKPoints
   use O_AtomicSites, only: coreDim,valeDim
   use O_TimeStamps

   ! Use the HDF5 module.
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Log the time we start to pscf the HDF5 files.
   call timeStampStart(27)

   ! Initialize the Fortran 90 HDF5 interface.
   call h5open_f(hdferr)
   if (hdferr < 0) stop 'Failed to open HDF interface'

   ! Initialize the values of the dimension variables.
   valeValeDims(1)   = valeDim
   valeValeDims(2)   = valeDim
   coreValeDims(1)   = coreDim
   coreValeDims(2)   = valeDim
   valeStatesDims(1) = valeDim
   valeStatesDims(2) = numStates
   statesDims(1)     = numStates
   attribIntDims(1)  = 1

   ! We have a couple of complications.
   ! First, in the current parallel algorithm there is a distributed global
   !   array that each process will contribute to. Initially each process
   !   computes a subset of the matrix elements and throws that subset up
   !   into the GA. Then each process grabs a different portion of the GA
   !   for the purpose of orthogonalization, participates in the
   !   orthogonalization calculation, and stores its result back into the
   !   GA. Finally, *one* process will pull down pieces of the GA and
   !   write them to disk serially. The only time that HDF5 enters is for
   !   the time that the one process writes to disk.
   ! Second complication: the one process cannot hold the entire matrix
   !   in memory at once because it will often be too large. Therefore,
   !   it needs to only pull down chunks of the matrix from the GA and
   !   write those. Obviously, it would like to pull down chunk sizes
   !   that are as large as possible to make efficient use of the cache
   !   and the HDF5 compression routines (that operate only across a
   !   chunk).
   ! Third complication: there is a maximum chunk size (as of the current
   !   version of HDF5 that is on the order of 2GB). Therefore, I don't want
   !   to have a number of chunk elements that is greater than 250M (assuming
   !   that 8-byte reals are being stored).
   ! Resolution: check that the chunk size is not too large (the assumption
   !   here is that the number being stored are 8 byte reals and that we
   !   should not go over 2 billion bytes. E.g. Let x=y=valeDim & a=b=chunkDim
   !   for simple notation here. So if x*y > 250M and we want a*b < 250M but
   !   of a size such that n*(a*b) is only slightly > x*y where n is a perfect
   !   square (4,9,16, etc.) then we need to think carefully. We want to
   !   maximize the chunk size, but not go overboard or else we will have
   !   some huge chunks just to cover some thin border regions.
   call getChunkDims(valeValeDims(1),valeValeDims(2),valeValeDimsChunk, &
         & 250000000)
   call getChunkDims(coreValeDims(1),coreValeDims(2),coreValeDimsChunk, &
         & 250000000)
   call getChunkDims(valeStatesDims(1),valeStatesDims(2),valeStatesDimsChunk, &
         & 250000000)

   ! Create the property list for the pscf hdf5 file and turn off chunk caching.
   call h5pcreate_f    (H5P_FILE_ACCESS_F,pscf_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create pscf plid.'
   call h5pget_cache_f (pscf_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to get pscf plid cache settings.'
   call h5pset_cache_f (pscf_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to set pscf plid cache settings.'

   ! Turn on a strong close so that when a close is requested all open objects
   !   in the file are also closed.
   call h5pset_fclose_degree_f(pscf_plid, H5F_CLOSE_STRONG_F, hdferr)
   if (hdferr /= 0) stop 'Failed to set fclose degree in pscf_plid'

   ! Create the HDF5 file that will hold all the computed results.  This will
   !   die if the file already exists.  It also uses the default file creation
   !   and file access properties.
   call h5fcreate_f ("pscf-temp.hdf5",H5F_ACC_EXCL_F,pscf_fid,hdferr,&
         & H5P_DEFAULT_F,pscf_plid)
   if (hdferr /= 0) stop 'Failed to create pscf-temp.hdf5 file.'

   ! Create the top level groups in the pscf output file.
   call h5gcreate_f (pscf_fid,"eigenVectors",eigenVec_gid,hdferr)
   if (hdferr == -1) stop 'Cannot create eigenVec_gid'
   call h5gcreate_f (pscf_fid,"eigenValues",eigenVal_gid,hdferr)
   if (hdferr == -1) stop 'Cannot create eigenVal_gid'
   call h5gcreate_f (pscf_fid,"overlap",overlapVV_gid,hdferr)
   if (hdferr == -1) stop 'Cannot create overlapVV_gid'
   if (coreDim /= 0) then
      call h5gcreate_f (pscf_fid,"orthoCoeff",overlapCV_gid,hdferr)
      if (hdferr == -1) stop 'Cannot create overlapCV_gid'
   endif
   call h5gcreate_f (pscf_fid,"momentumX",momentumX_gid,hdferr)
   if (hdferr == -1) stop 'Cannot create momentumX_gid'
   call h5gcreate_f (pscf_fid,"momentumX",momentumY_gid,hdferr)
   if (hdferr == -1) stop 'Cannot create momentumY_gid'
   call h5gcreate_f (pscf_fid,"momentumX",momentumZ_gid,hdferr)
   if (hdferr == -1) stop 'Cannot create momentumZ_gid'

   ! Create the dataspaces that will be used for each group.
   call h5screate_simple_f (2,valeStatesDims,valeStates_dsid,hdferr) !eigenVec
   if (hdferr == -1) stop 'Cannot create valeStates_dsid'
   call h5screate_simple_f (1,statesDims,states_dsid,hdferr) !eigenVal
   if (hdferr == -1) stop 'Cannot create states_dsid'
   call h5screate_simple_f (2,valeValeDims,valeVale_dsid,hdferr) !overlap+mom
   if (hdferr == -1) stop 'Cannot create valeVale_dsid'
   if (coreDim /= 0) then
      callh5screate_simple_f (2,coreValeDims,coreVale_dsid,hdferr) !orthoCoeffs
      if (hdferr == -1) stop 'Cannot create coreVale_dsid'
   endif
   call h5screate_simple_f(1,attribIntDims,attribInt_dsid,hdferr) !mom comptd
   if (hdferr /= 0) stop 'Failed to create attribInt_dsid.'

   ! Create the property list for each dataspace (except the attribute, it
   !   just uses a default property list for no reason other than lazyness).
   call h5pcreate_f (H5P_DATASET_CREATE_F,valeStates_plid,hdferr) !eigenVec
   if (hdferr == -1) stop 'Cannot create valeStates_plid'
   call h5pcreate_f (H5P_DATASET_CREATE_F,states_plid,hdferr) !eigenVal
   if (hdferr == -1) stop 'Cannot create states_plid'
   call h5pcreate_f (H5P_DATASET_CREATE_F,valeVale_plid,hdferr) !overlap+mom
   if (hdferr == -1) stop 'Cannot create valeVale_plid'
   if (coreDim /= 0) then
      call h5pcreate_f (H5P_DATASET_CREATE_F,coreVale_plid,hdferr) !mom comptd
      if (hdferr == -1) stop 'Cannot create coreVale_plid'
   endif
   call h5pset_layout_f (valeStates_plid,H5D_CHUNKED,hdferr)
   call h5pset_layout_f (states_plid,H5D_CHUNKED,hdferr)
   call h5pset_layout_f (valeVale_plid,H5D_CHUNKED,hdferr)
   if (coreDim /= 0) then
      call h5pset_layout_f (coreVale_plid,H5D_CHUNKED,hdferr)
   endif
   call h5pset_chunk_f (valeStates_plid,2,valeStatesDimsChunk,hdferr)
   call h5pset_chunk_f (states_plid,2,statesDims,hdferr) ! No sub-chunk.
   call h5pset_chunk_f (valeVale_plid,2,valeValeDimsChunk,hdferr)
   if (coreDim /= 0) then
      call h5pset_chunk_f (coreVale_plid,2,coreValeDimsChunk,hdferr)
   endif
   call h5pset_deflate_f (valeStates_plid,1,hdferr)
   call h5pset_deflate_f (states_plid,1,hdferr)
   call h5pset_deflate_f (valeVale_plid,1,hdferr)
   if (coreDim /= 0) then
      call h5pset_deflate_f (coreVale_plid,1,hdferr)
   endif

   ! Create the computed status attribute for the momentumX_gid.
   call h5acreate_f (momentumX_gid,"computed",H5T_NATIVE_INTEGER,&
         & attribInt_dsid,momentumX_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create momentumX_aid.'

   ! Allocate space to hold the IDs for the datasets in each group.
#ifndef GAMMA
   allocate (eigenVec_did(2,numKPoints,spin))
#else
   allocate (eigenVec_did(1,numKPoints,spin))
#endif
   allocate (eigenVal_did(numKPoints,spin))
   allocate (overlapVV_did(numKPoints))
   if (coreDim /= 0) then
#ifndef GAMMA
      allocate (overlapCV_did(2,numKPoints))
#else
      allocate (overlapCV_did(1,numKPoints))
#endif
   endif
   ! Even if we don't need the momentum matrices we will make the datasets.
   allocate (momentumX_did(numKPoints))
   allocate (momentumY_did(numKPoints))
   allocate (momentumZ_did(numKPoints))

   ! Create the datasets.
   do i = 1, numKPoints

      ! Make the eigenVector datasets (real and imaginary (if needed)) and the
      !   eigenValue datasets.
      do j = 1, spin
         write (currentName,fmt="(a4,i7.7,i7.7)") "real",i,j
         currentName = trim(currentName)
         call h5dcreate_f(eigenVec_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeStates_dsid,eigenVec_did(1,i,j),hdferr,valeStates_plid)
         if (hdferr == -1) stop 'Cannot create eigenVec_did 1'
#ifndef GAMMA
         write (currentName,fmt="(a4,i7.7,i7.7)") "imag",i,j
         currentName = trim(currentName)
         call h5dcreate_f(eigenVec_gid,currentName,H5T_NATIVE_DOUBLE,&
               & valeStates_dsid,eigenVec_did(2,i,j),hdferr,valeStates_plid)
         if (hdferr == -1) stop 'Cannot create eigenVec_did 2'
#endif
         write (currentName,fmt="(i7.7,i7.7)") i,j
         currentName = trim(currentName)
         call h5dcreate_f(eigenVal_gid,currentName,H5T_NATIVE_DOUBLE,&
               & states_dsid,eigenVal_did(i,j),hdferr,states_plid)
         if (hdferr == -1) stop 'Cannot create eigenVal_did'
      enddo

      ! Make overlapVV, overlapCV, and momentum[XYZ] datasets (even if we don't
      !   need the momentum datasets at this point in time because it doesn't
      !   hurt and it makes life later much easier when/if we eventually need
      !   to fill them). 
      write (currentName,fmt="(i7.7)") i
      currentName = trim(currentName)
      call h5dcreate_f(overlapVV_gid,currentName,H5T_NATIVE_DOUBLE,&
            & valeVale_dsid,overlapVV_did(i),hdferr,valeVale_plid)
      if (hdferr == -1) stop 'Cannot create overlapVV_did'
      if (coreDim /= 0) then
         write (currentName,fmt="(a4,i7.7)") "real",i
         currentName = trim(currentName)
         call h5dcreate_f(overlapCV_gid,currentName,H5T_NATIVE_DOUBLE,&
               & coreVale_dsid,overlapCV_did(1,i),hdferr,coreVale_plid)
         if (hdferr == -1) stop 'Cannot create overlapCV_did 1'
#ifndef GAMMA
         write (currentName,fmt="(a4,i7.7)") "imag",i
         currentName = trim(currentName)
         call h5dcreate_f(overlapCV_gid,currentName,H5T_NATIVE_DOUBLE,&
               & coreVale_dsid,overlapCV_did(2,i),hdferr,coreVale_plid)
         if (hdferr == -1) stop 'Cannot create overlapCV_did 2'
#endif
      endif
      write (currentName,fmt="(i7.7)") i
      currentName = trim(currentName)
      call h5dcreate_f(momentumX_gid,currentName,H5T_NATIVE_DOUBLE,&
            & valeVale_dsid,momentumX_did(i),hdferr,valeVale_plid)
      if (hdferr == -1) stop 'Cannot create momentumX_did'
      call h5dcreate_f(momentumY_gid,currentName,H5T_NATIVE_DOUBLE,&
            & valeVale_dsid,momentumY_did(i),hdferr,valeVale_plid)
      if (hdferr == -1) stop 'Cannot create momentumY_did'
      call h5dcreate_f(momentumZ_gid,currentName,H5T_NATIVE_DOUBLE,&
            & valeVale_dsid,momentumZ_did(i),hdferr,valeVale_plid)
      if (hdferr == -1) stop 'Cannot create momentumZ_did'
   enddo

   ! Log the time we finish setting up the HDF5 files.
   call timeStampEnd(27)

end subroutine initPSCFHDF5

! This is only called by the pscf program. Other programs that use the data
!   in the pscf hdf5 file will use the access and closeAccess subroutines.
subroutine closePSCFHDF5

   ! Use the HDF5 module.
   use HDF5

   ! Import other necessary data modules.
   use O_Potential,   only: spin
   use O_KPoints,     only: numKPoints
   use O_AtomicSites, only: coreDim

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare local variables.
   integer :: i,j
   integer :: hdferr

   ! Close the dataspaces first.
   call h5close_f (valeStates_dsid,hdferr)
   if (hdferr == -1) stop 'Cannot close valeStates_dsid'
   call h5close_f (states_dsid,hdferr)
   if (hdferr == -1) stop 'Cannot close states_dsid'
   call h5close_f (valeVale_dsid,hdferr)
   if (hdferr == -1) stop 'Cannot close valeVale_dsid'
   if (coreDim /= 0) then
      call h5close_f (coreVale_dsid,hdferr)
      if (hdferr == -1) stop 'Cannot close coreVale_dsid'
   endif

   ! Close the datasets.
   do i = 1, numKPoints

      ! Close the eigenVector and eigenValue datasets.
      do i = 1, spin
         call h5dclose_f (eigenVec_did(1,i,j),hdferr)
         if (hdferr == -1) stop 'Cannot close eigenVec_did 1'
#ifndef GAMMA
         call h5dclose_f (eigenVec_did(1,i,j),hdferr)
         if (hdferr == -1) stop 'Cannot close eigenVec_did 1'
#endif
         call h5dclose_f (eigenVal_did(i,j),hdferr)
         if (hdferr == -1) stop 'Cannot close eigenVal_did'
      enddo

      ! Close the overlapVV, overlapCV, and momentum[XYZ] datasets.
      call h5dclose_f (overlapVV_did(i),hdferr)
      if (hdferr == -1) stop 'Cannot close overlapVV_did'
      if (coreDim /= 0) then
         call h5dclose_f (overlapCV_did(1,i),hdferr)
         if (hdferr == -1) stop 'Cannot close overlapCV_did 1'
#ifndef GAMMA
         call h5dclose_f (overlapCV_did(2,i),hdferr)
         if (hdferr == -1) stop 'Cannot close overlapCV_did 2'
#endif
      endif
      call h5dclose_f (momentumX_did(i),hdferr)
      if (hdferr == -1) stop 'Cannot close momentumX_did'
      call h5dclose_f (momentumY_did(i),hdferr)
      if (hdferr == -1) stop 'Cannot close momentumY_did'
      call h5dclose_f (momentumZ_did(i),hdferr)
      if (hdferr == -1) stop 'Cannot close momentumZ_did'
   enddo

   ! Close the property lists.
   call h5pclose_f (valeStates_plid,hdferr)
   if (hdferr == -1) stop 'Cannot close valeStates_plid'
   call h5pclose_f (states_plid,hdferr)
   if (hdferr == -1) stop 'Cannot close states_plid'
   call h5pclose_f (valeVale_plid,hdferr)
   if (hdferr == -1) stop 'Cannot close valeVale_plid'
   if (coreDim /= 0) then
      call h5pclose_f (coreVale_plid,hdferr)
      if (hdferr == -1) stop 'Cannot close coreVale_plid'
   endif

   ! Close the top level groups.
   call h5gclose_f (eigenVec_gid,hdferr)
   if (hdferr == -1) stop 'Cannot close eigenVec_gid'
   call h5gclose_f (eigenVal_gid,hdferr)
   if (hdferr == -1) stop 'Cannot close eigenVal_gid'
   call h5gclose_f (overlapVV_gid,hdferr)
   if (hdferr == -1) stop 'Cannot close overlapVV_gid'
   if (coreDim /= 0) then
      call h5gclose_f (overlapCV_gid,hdferr)
      if (hdferr == -1) stop 'Cannot close overlapCV_gid'
   endif
   call h5gclose_f (momentumX,hdferr)
   if (hdferr == -1) stop 'Cannot close momentumX_gid'
   call h5gclose_f (momentumY,hdferr)
   if (hdferr == -1) stop 'Cannot close momentumY_gid'
   call h5gclose_f (momentumZ,hdferr)
   if (hdferr == -1) stop 'Cannot close momentumZ_gid'

   ! Close the file property list.
   call h5pclose_f (pscf_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf_plid.'

   ! Close the file.
   call h5fclose_f (pscf_fid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf_fid.'

   ! Closer the HDF5 interface.
   call h5close_f(hdferr)
   if (hdferr /=0) stop 'Failed to close the HDF5 interace.'

   ! Deallocate arrays.
   deallocate (eigenVec_did)
   deallocate (eigenVal_did)
   deallocate (overlapVV_did)
   if (coreDim /= 0) then
      deallocate (overlapCV_did)
   endif
   deallocate (momentumX_did)
   deallocate (momentumY_did)
   deallocate (momentumZ_did)

end subroutine closePSCFHDF5

subroutine accessPSCFHDF5

   ! Use the HDF5 module.
   use HDF5

   ! Import necessary modules.
   use O_Potential,   only: spin
   use O_CommandLine, only: stateSet
   use O_KPoints,     only: numKPoints
   use O_AtomicSites, only: coreDim, valeDim

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize the Fortran 90 HDF5 interface.
   call h5open_f (hdferr)
   if (hdferr /= 0) stop 'Failed to open the HDF5 interface'

   ! Initialize the values of the dimension variables.
   valeValeDims(1)   = valeDim
   valeValeDims(2)   = valeDim
   coreValeDims(1)   = coreDim
   coreValeDims(2)   = valeDim
   valeStatesDims(1) = valeDim
   valeStatesDims(2) = numStates
   statesDims(1)     = numStates
   attribIntDims(1)  = 1

   ! Create the property list for the pscf hdf5 file and turn off chunk caching.
   call h5pcreate_f (H5P_FILE_ACCESS_F,pscf_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create pscf_plid.'
   call h5pget_cache_f (pscf_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to get pscf_plid cache settings.'
   call h5pset_cache_f (pscf_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to set pscf_plid cache settings.'


   ! Open the HDF5 file that will hold all the computed results.  Read only.
   call h5fopen_f ("pscf-temp.hdf5",H5F_ACC_RDONLY_F,pscf_fid,hdferr,&
         & pscf_plid)
   if (hdferr /= 0) stop 'Failed to open pscf-temp.hdf5 file.'


   ! Open the top level groups in the pscf output file.
   call h5gopen_f (pscf_fid,"eigenVectors",eigenVec_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open eigenVec_gid.'
   call h5gopen_f (pscf_fid,"eigenValues",eigenVal_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open eigenVal_gid.'
   call h5gopen_f (pscf_fid,"overlap",overlapVV_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open overlapVV_gid.'
   if (coreDim /= 0) then
      call h5gopen_f (pscf_fid,"orthoCoeff",overlapCV_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to open overlapCV_gid.'
   endif
   call h5gopen_f (pscf_fid,"momentumX",momentumX_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open momentumX_gid.'
   call h5gopen_f (pscf_fid,"momentumX",momentumY_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open momentumY_gid.'
   call h5gopen_f (pscf_fid,"momentumX",momentumZ_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open momentumZ_gid.'


   ! Allocate space to hold the IDs for the datasets in each group.
#ifndef GAMMA
   allocate (eigenVec_did(2,numKPoints,spin))
#else
   allocate (eigenVec_did(1,numKPoints,spin))
#endif
   allocate (eigenVal_did(numKPoints,spin))
   allocate (overlapVV_did(numKPoints))
   if (coreDim /= 0) then
#ifndef GAMMA
      allocate (overlapCV_did(2,numKPoints))
#else
      allocate (overlapCV_did(1,numKPoints))
#endif
   endif
   allocate (momentumX_did(numKPoints))
   allocate (momentumY_did(numKPoints))
   allocate (momentumZ_did(numKPoints))


   ! Open the datasets (in the case of the momentum we may have to create).
   do i = 1, numKPoints

      ! Open the datasets for the overlapCV, overlapVV, and momentum matrices.
      if (coreDim /= 0) then
         write (currentName,fmt="(a4,i7.7)") "real",i
         currentName = trim (currentName)
         call h5dopen_f (overlapCV_gid,currentName,overlapCV_did(1,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open overlapCV_did 1'
#ifndef GAMMA
         write (currentName,fmt="(a4,i7.7)") "imag",i
         currentName = trim (currentName)
         call h5dopen_f (overlapCV_gid,currentName,overlapCV_did(2,i),hdferr)
         if (hdferr /= 0) stop 'Failed to open overlapCV_did 2'
#endif
      endif
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)
      call h5dopen_f (overlapVV_gid,currentName,overlapVV_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open overlapVV_did'
      call h5dopen_f (momentumX_gid,currentName,momentumX_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open momentumX_did'
      call h5dopen_f (momentumX_gid,currentName,momentumY_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open momentumY_did'
      call h5dopen_f (momentumX_gid,currentName,momentumZ_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open momentumZ_did'

      ! Open the datasets for the eigenVectors and eigenValues.
      do j = 1, spin
         write (currentName,fmt="(a4,i7.7,i7.7)") "real",i,j
         currentName = trim (currentName)
         call h5dopen_f (eigenVec_gid,currentName,eigenVec_did(1,i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVec_did 1'
#ifndef GAMMA
         write (currentName,fmt="(a4,i7.7,i7.7)") "imag",i,j
         currentName = trim (currentName)
         call h5dopen_f (eigenVec_gid,currentName,eigenVec_did(2,i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVec_did 2'
#endif
         write (currentName,fmt="(i7.7,i7.7)") i,j
         currentName = trim (currentName)
         call h5dopen_f (eigenVal_gid,currentName,eigenVal_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open eigenVal_did'
      enddo ! j = spin

   enddo ! i = numKPoints


   ! In the case that we are doing a PACS type optc calculation, then we need
   !   to also access data in a second pscf hdf5 file.
   if (stateSet == 1) then

      ! Create the property list for the band hdf5 file and turn off
      !   chunk caching.
      call h5pcreate_f    (H5P_FILE_ACCESS_F,pscf2_plid,hdferr)
      if (hdferr /= 0) stop 'Failed to create pscf2 plid.'
      call h5pget_cache_f (pscf2_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,&
            & rdcc_w0,hdferr)
      if (hdferr /= 0) stop 'Failed to get pscf2 plid cache settings.'
      call h5pset_cache_f (pscf2_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to set pscf2 plid cache settings.'

      ! Turn on a strong close so that when a close is requested all open
      !   objects in the file are also closed.
      call h5pset_fclose_degree_f(pscf_plid, H5F_CLOSE_STRONG_F, hdferr)
      if (hdferr /= 0) stop 'Failed to set fclose degree in pscf2 plid'

      ! Open the second band calculation file.
      call h5fopen_f ("band-temp.2.hdf5",H5F_ACC_RDONLY_F,pscf2_fid,hdferr,&
         & pscf2_plid)
      if (hdferr /= 0) stop 'Failed to open the pscf2 file'

      ! Open the top level groups of the second pscf calculation file.
      call h5gopen_f (pscf2_fid,"eigenVectors",eigenVec2_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to open eigenVec2_gid.'
      call h5gopen_f (pscf2_fid,"eigenValues",eigenVal2_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to open eigenVal2_gid.'
      call h5gopen_f (pscf2_fid,"overlap",overlapVV2_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to open overlapVV2_gid'

      ! Allocate space for the dataset ID numbers in the pscf output file.
      allocate (eigenVec2_did  (2,numKPoints,spin))
      allocate (eigenVal2_did  (numKPoints,spin))
      allocate (overlapVV2_did (numKPoints))

      ! Open the datasets in the second pscf output file for the overlap.
      do i = 1, numKPoints
         write (currentName,fmt="(i7.7)") i
         currentName = trim (currentName)
         call h5dopen_f (overlapVV2_gid,currentName,overlapVV2_did(i),hdferr)
         if (hdferr /= 0) stop 'Failed to open overlapVV2_did'
      enddo


      ! Open the datasets in the second pscf output file for the wave functions.
      do i = 1, numKPoints
         do j = 1, spin
            write (currentName,fmt="(a4,i7.7,i7.7)") "real",i,j
            currentName = trim (currentName)
            call h5dopen_f (eigenVec2_gid,currentName,eigenVec2_did(1,i,j),&
                  & hdferr)
            if (hdferr /= 0) stop 'Failed to open eigenVec2_did 1'
#ifndef GAMMA
            write (currentName,fmt="(a4,i7.7,i7.7)") "imag",i,j
            currentName = trim (currentName)
            call h5dopen_f (eigenVec2_gid,currentName,eigenVec2_did(2,i,j),&
                  & hdferr)
            if (hdferr /= 0) stop 'Failed to open eigenVec2_did 2'
#endif
            write (currentName,fmt="(i7.7,i7.7)") i,j
            currentName = trim (currentName)
            call h5dopen_f (eigenVal2_gid,currentName,eigenVal2_did(i,j),hdferr)
            if (hdferr /= 0) stop 'Failed to open eigenVal2_did'
         enddo
      enddo
   endif

end subroutine accessPSCFHDF5

subroutine closeAccessPSCFHDF5

   ! Import the necessary HDF data modules
   use HDF5

   ! Import the necessary data modules.
   use O_Potential,   only: spin
   use O_CommandLine, only: stateSet
   use O_KPoints,     only: numKPoints
   use O_AtomicSites, only: coreDim

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Define local loop control and error control variables
   integer :: i,j
   integer :: hdferr

   ! Close pscf datasets for the overlapCV, overlapVV, and momentum matrices.
   do i = 1, numKPoints
      if (coreDim /= 0) then
         call h5dclose_f (overlapCV_did(1,i),hdferr)
         if (hdferr == -1) stop 'Cannot close ovarlapCV_did 1'
#ifndef GAMMA
         call h5dclose_f (overlapCV_did(2,i),hdferr)
         if (hdferr == -1) stop 'Cannot close ovarlapCV_did 1'
#endif
      endif
      call h5dclose_f (overlapVV_did(i),hdferr)
      if (hdferr == -1) stop 'Cannot close overlapVV_did'
      call h5dclose_f (momentumX_did(i),hdferr)
      if (hdferr == -1) stop 'Cannot close momentumX_did'
      call h5dclose_f (momentumY_did(i),hdferr)
      if (hdferr == -1) stop 'Cannot close momentumY_did'
      call h5dclose_f (momentumZ_did(i),hdferr)
      if (hdferr == -1) stop 'Cannot close momentumZ_did'
   enddo

   ! Close the pscf datasets for the eigenvalues and eigenvectors next.
   do i = 1, numKPoints
      do j = 1, spin
         call h5dclose_f (eigenVec_did(1,i,j),hdferr)
         if (hdferr == -1) stop 'Cannot close eigenVec_did 1'
#ifndef GAMMA
         call h5dclose_f (eigenVec_did(2,i,j),hdferr)
         if (hdferr == -1) stop 'Cannot close eigenVec_did 2'
#endif
         call h5dclose_f (eigenVal_did(i,j),hdferr)
         if (hdferr == -1) stop 'Cannot close eigenVal_did'
      enddo
   enddo

   ! Close the pscf property lists.
   call h5pclose_f (valeStates_plid,hdferr)
   if (hdferr == -1) stop 'Cannot close valeStates_plid'
   call h5pclose_f (states_plid,hdferr)
   if (hdferr == -1) stop 'Cannot close states_plid'
   call h5pclose_f (valeVale_plid,hdferr)
   if (hdferr == -1) stop 'Cannot close valeVale_plid'
   if (coreDim /= 0) then
      call h5pclose_f (coreVale_plid,hdferr)
      if (hdferr == -1) stop 'Cannot close coreVale_plid'
   endif

   ! Close the pscf top level groups.
   call h5gclose_f (eigenVec_gid,hdferr)
   if (hdferr == -1) stop 'Cannot close eigenVec_gid'
   call h5gclose_f (eigenVal_gid,hdferr)
   if (hdferr == -1) stop 'Cannot close eigenVal_gid'
   if (coreValeBand(1) /= 0) then
      call h5gclose_f (overlapCV_gid,hdferr)
      if (hdferr == -1) stop 'Cannot close overlapCV_gid'
   endif
   call h5gclose_f (overlapVV_gid,hdferr)
   if (hdferr == -1) stop 'Cannot close overlapVV_gid'
   call h5gclose_f (momentumX_gid,hdferr)
   if (hdferr == -1) stop 'Cannot close momentumX'
   call h5gclose_f (momentumY_gid,hdferr)
   if (hdferr == -1) stop 'Cannot close momentumY'
   call h5gclose_f (momentumZ_gid,hdferr)
   if (hdferr == -1) stop 'Cannot close momentumZ'

   ! Close the file property list.
   call h5pclose_f (pscf_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pscf_plid.'

   ! Close the pscf file.
   call h5fclose_f (pscf_fid,hdferr)
   if (hdferr == -1) stop 'Cannot close pscf_fid'

   ! Deallocate arrays.
   deallocate (eigenVec_did)
   deallocate (eigenVal_did)
   deallocate (overlapVV_did)
   if (coreDim /= 0) then
      deallocate (overlapCV_did)
   endif
   deallocate (momentumX_did)
   deallocate (momentumY_did)
   deallocate (momentumZ_did)

   ! In the case that we are closing access to the pscf hdf5 file from a
   !   PACS type calculation we need to also close the second file that was
   !   opened for the excited state.
   if (stateSet == 1) then

      ! Close the overlap datasets in the second pscf file.
      do i = 1, numKPoints
         call h5dclose_f (overlapVV2_did(i),hdferr)
         if (hdferr /= 0) stop 'Failed to close overlapVV2_did'
      enddo

      ! Close the wavefunction datasets in the second band file.
      do i = 1, numKPoints
         do j = 1, spin
            call h5dclose_f (eigenVec2_did(1,i,j),hdferr)
            if (hdferr /= 0) stop 'Failed to close eigenVec2_did 1'
#ifndef GAMMA
            call h5dclose_f (eigenVec2_did(2,i,j),hdferr)
            if (hdferr /= 0) stop 'Failed to close eigenVec2_did 2'
#endif
            call h5dclose_f (eigenVal2_did(i,j),hdferr)
            if (hdferr /= 0) stop 'Failed to close eigenVal2_did'
         enddo
      enddo

      ! Close the top level groups.
      call h5gclose_f (eigenVec2_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to close eigenVec2_gid'
      call h5gclose_f (eigenVal2_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to close eigenVal2_gid'
      call h5gclose_f (overlapVV2_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to close overlapVV2_gid'

      ! Close the HDF5 file access property list
      call h5pclose_f (pscf2_plid,hdferr)
      if (hdferr /= 0) stop 'Failed to close pscf2_plid'

      ! Close the HDF5 file that holds the results of the band calculation.
      call h5fclose_f (pscf2_fid,hdferr)
      if (hdferr /= 0) stop 'Failed to close pscf2_fid'

      ! Deallocate arrays.
      deallocate (eigenVec2_did)
      deallocate (eigenVal2_did)
      deallocate (overlapVV2_did)
   endif


end subroutine closeAccessPSCFHDF5

subroutine getChunkDims (dim1, dim2, chunkDims, maxChunkBytes)

   ! Use necessary modules
   use HDF5

   ! Define passed parameters
   integer(hsize_t), intent(in) :: dim1
   integer(hsize_t), intent(in) :: dim2
   integer(hsize_t), dimension(2), intent(out) :: chunkDims
   integer, intent(in) :: maxChunkBytes

   ! Define local variables
   integer :: i

   ! We need to deal with the general situation that the dimensions of the
   !   matrix that will be stored on disk may not be equal. We want the largest
   !   chunk that is less than maxChunkBytes (nominally 250 million as of HDF5
   !   version 1.8.6). So, if we can't use one chunk to cover the whole matrix
   !   then we divide the full matrix dimensions by increasingly large
   !   integers. We err on the side of making a chunk just a bit too large so
   !   that only a small "sliver" will extend over the full stored matrix.
   !   This will attempt to minimize waste on the disk.
   if ((dim1 * dim2) > maxChunkBytes) then
      i = 1
      do while (1)

         ! The current divisor. (Start with chunks of size 1/2 the full matrix,
         !   then 1/3, then 1/4,... until we find a chunk size that is < 250M.)
         i = i + 1

         ! Get the first dimension and then the second.
         if (mod(dims1,i) == 0) then
            chunkDims(1) = (dim1 / i)
         else
            chunkDims(1) = (dim1 / i) + mod(dim1,i)
         endif
         if (mod(dims2,i) == 0) then
            chunkDims(2) = (dim2 / i)
         else
            chunkDims(2) = (dim2 / i) + mod(dim2,i)
         endif

         ! Compare to the max number of matrix element (bytes).
         if ((chunkDims(1)*chunkDims(2)) < maxChunkBytes) then
             exit
         endif
      enddo
   else ! One chunk will cover the whole matrix.
      chunkDims(1) = dim1
      chunkDims(2) = dim2
   endif

end subroutine getChunkDims

end module O_PSCFHDF5
