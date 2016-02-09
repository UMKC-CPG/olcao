module O_PSCFBandHDF5

   ! Use the HDF5 module for HDF5 defined types (e.g. size_t and hid_t).
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define arrays that hold the dimensions of the datasets.
   integer(hsize_t), dimension (2) :: valeStatesBand
   integer(hsize_t), dimension (2) :: valeValeBand
   integer(hsize_t), dimension (2) :: coreValeBand
   integer(hsize_t), dimension (1) :: statesBand

   ! Define arrays that hold the dimensions of the chunks (if they might lead
   !   to chunks that are too big).
   integer(hsize_t), dimension (2) :: valeStatesBandChunk
   integer(hsize_t), dimension (2) :: valeValeBandChunk
   integer(hsize_t), dimension (2) :: coreValeBandChunk

   ! Define the file ID.
   integer(hid_t) :: band_fid

   ! Define the input file access property list ID. 
   integer(hid_t)   :: band_plid
   integer          :: mdc_nelmts  ! Meta-data cache num elements.
   integer (size_t) :: rdcc_nelmts ! Raw-data chunk cache num elements.
   integer (size_t) :: rdcc_nbytes ! Raw-data chunk cache size in bytes.
   real             :: rdcc_w0     ! Raw-data chunk cache weighting parameter.

   ! For the band output there will be four groups.
   integer(hid_t) :: eigenVectorsBand_gid
   integer(hid_t) :: eigenValuesBand_gid
   integer(hid_t) :: valeValeBand_gid
   integer(hid_t) :: coreValeBand_gid

   ! For the band output there will be four types of dataspaces for the four 
   !   types of data in the above groups.
   integer(hid_t) :: valeStatesBand_dsid
   integer(hid_t) :: statesBand_dsid
   integer(hid_t) :: valeValeBand_dsid
   integer(hid_t) :: coreValeBand_dsid

   ! For the band output dataspaces there will also be four associated
   !   property lists.
   integer(hid_t) :: valeStatesBand_plid
   integer(hid_t) :: statesBand_plid
   integer(hid_t) :: valeValeBand_plid
   integer(hid_t) :: coreValeBand_plid

   ! There will also be four groups of data sets for the band output.
   integer(hid_t), allocatable, dimension (:,:,:) :: eigenVectorsBand_did
   integer(hid_t), allocatable, dimension (:,:)   :: eigenValuesBand_did
   integer(hid_t), allocatable, dimension (:)     :: valeValeBand_did
   integer(hid_t), allocatable, dimension (:,:)   :: coreValeBand_did

   ! In the case of reading the band data it is possible that two band data
   !   files will need to be opened at the same time for PACS calculations.
   !   In this case a second set of group, property list, dataset, and file IDs
   !   will be needed.

   ! Second file ID.
   integer(hid_t) :: band2_fid

   ! Second property list ID.
   integer(hid_t)   :: band2_plid

   ! Second group IDs.
   integer(hid_t) :: eigenVectorsBand2_gid
   integer(hid_t) :: eigenValuesBand2_gid
   integer(hid_t) :: valeValeBand2_gid

   ! Second dataset IDs.
   integer(hid_t), allocatable, dimension (:,:,:) :: eigenVectorsBand2_did
   integer(hid_t), allocatable, dimension (:,:)   :: eigenValuesBand2_did
   integer(hid_t), allocatable, dimension (:)     :: valeValeBand2_did

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subs.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

subroutine initPSCFBandHDF5(numStates)

   ! Import the necessary HDF modules
   use HDF5

   ! Import the necessary data modules.
   use O_Kinds
   use O_Potential,   only: spin
   use O_CommandLine, only: doSYBD
   use O_KPoints,     only: numKPoints
   use O_AtomicSites, only: coreDim, valeDim

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: numStates

   ! Define local variables that will be used to create the dataspaces etc.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize the Fortran 90 HDF5 interface.
   call h5open_f(hdferr)

   ! Initialize the values of the dimension variables.
   valeStatesBand(1)    = valeDim
   valeStatesBand(2)    = numStates
#ifndef GAMMA
   valeValeBand(1)      = 2
#else
   valeValeBand(1)      = 1
#endif
   valeValeBand(2)      = valeDim*(valeDim+1)/2
   coreValeBand(1)      = coreDim
   coreValeBand(2)      = valeDim
   statesBand(1)        = numStates

   ! Check that the chunk size is not too large (the assumption here is that
   !   the number being stored are 8 byte reals and that we should not go over
   !   2 billion bytes. Note that if x*y = >250M and we want a*b = 250M then
   !   the additional requirement x/y = a/b leads to b = sqrt(250M/>250M)*y.
   !   Thus a = 250M/b.
   if (valeStatesBand(1) * valeStatesBand(2) > 250000000) then
      valeStatesBandChunk(2) = int(sqrt(real(250000000,double) / &
            & real(valeStatesBand(1) * valeStatesBand(2),double)) * &
            & valeStatesBand(2))
      valeStatesBandChunk(1) = int(250000000 / valeStatesBandChunk(2))
   else
      valeStatesBandChunk(1) = valeStatesBand(1)
      valeStatesBandChunk(2) = valeStatesBand(2)
   endif

   ! We will do a similar procedure for the coreValeBandChunk.
   if (coreValeBand(1) * coreValeBand(2) > 250000000) then
      coreValeBandChunk(2) = int(sqrt(real(250000000,double) / &
            & real(coreValeBand(1) * coreValeBand(2),double)) * &
            & coreValeBand(2))
      coreValeBandChunk(1) = int(250000000 / coreValeBandChunk(2))
   else
      coreValeBandChunk(1) = coreValeBand(1)
      coreValeBandChunk(2) = coreValeBand(2)
   endif

   ! The procedure for valeValeBand is a bit different because the
   !   valeValeBand(1) will always be 1 or 2.
   if (valeValeBand(1) * valeValeBand(2) > 250000000) then
      valeValeBandChunk(2) = int(250000000/valeValeBand(1))
      valeValeBandChunk(1) = valeValeBand(1)
   else
      valeValeBandChunk(1) = valeValeBand(1)
      valeValeBandChunk(2) = valeValeBand(2)
   endif

   ! The symmetric band calculation will not use these HDF5 settings for output.
   if (doSybd .eq. 0) then

      ! Create the property list for the band hdf5 file and turn off
      !   chunk caching.
      call h5pcreate_f    (H5P_FILE_ACCESS_F,band_plid,hdferr)
      if (hdferr == -1) stop 'Can not create band_plid'
      call h5pget_cache_f (band_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,&
            & rdcc_w0,hdferr)
      if (hdferr == -1) stop 'Can not get h5p_cache'
      call h5pset_cache_f (band_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,&
            & hdferr)
      if (hdferr == -1) stop 'Can not set h5p_cache'
   

      ! Create the HDF5 file that will hold the results of the band calculation
      !   to be used by other programs (bond,optc,dos).
      call h5fcreate_f ("band-temp.hdf5",H5F_ACC_EXCL_F,band_fid,hdferr,&
            H5P_DEFAULT_F,band_plid)
      if (hdferr == -1) stop 'Can not create band_fid'


      ! Create the top level groups in the band output file.
      call h5gcreate_f (band_fid,"eigenVectors",eigenVectorsBand_gid,hdferr)
      if (hdferr == -1) stop 'Can not create eigenVectorsBand_gid'
      call h5gcreate_f (band_fid,"eigenValues",eigenValuesBand_gid,hdferr)
      if (hdferr == -1) stop 'Can not create eigenValuesBand_gid'
      call h5gcreate_f (band_fid,"overlap",valeValeBand_gid,hdferr)
      if (hdferr == -1) stop 'Can not create valeValeBand_gid'
      if (coreDim /= 0) then
         call h5gcreate_f (band_fid,"orthoCoeff",coreValeBand_gid,hdferr)
         if (hdferr == -1) stop 'Can not create coreValeBand_gid'
      endif

      ! Create the dataspaces that will be used for each group.
      call h5screate_simple_f (2,valeStatesBand,valeStatesBand_dsid,hdferr)
      if (hdferr == -1) stop 'Can not create valeStatesBand_dsid'
      call h5screate_simple_f (1,statesBand,statesBand_dsid,hdferr)
      if (hdferr == -1) stop 'Can not create statesBand_dsid'
      call h5screate_simple_f (2,valeValeBand,valeValeBand_dsid,hdferr)
      if (hdferr == -1) stop 'Can not create valeValeBand_dsid'
      if (coreDim /= 0) then
         call h5screate_simple_f (2,coreValeBand,coreValeBand_dsid,hdferr)
         if (hdferr == -1) stop 'Can not create coreValeBand_dsid'
      endif

      ! Create the property lists first, then set them properties one at a time.
      call h5pcreate_f      (H5P_DATASET_CREATE_F,valeStatesBand_plid,hdferr)
      if (hdferr == -1) stop 'Can not create valeStatesBand_plid'
      call h5pcreate_f      (H5P_DATASET_CREATE_F,statesBand_plid,hdferr)
      if (hdferr == -1) stop 'Can not create statesBand_plid'
      call h5pcreate_f      (H5P_DATASET_CREATE_F,valeValeBand_plid,hdferr)
      if (hdferr == -1) stop 'Can not create valeValeBand_plid'
      if (coreDim /= 0) then
         call h5pcreate_f      (H5P_DATASET_CREATE_F,coreValeBand_plid,hdferr)
         if (hdferr == -1) stop 'Can not create coreValeBand_plid'
      endif
      call h5pset_layout_f  (valeStatesBand_plid,H5D_CHUNKED_F,hdferr)
      call h5pset_layout_f  (statesBand_plid,H5D_CHUNKED_F,hdferr)
      call h5pset_layout_f  (valeValeBand_plid,H5D_CHUNKED_F,hdferr)
      if (coreDim /= 0) then
         call h5pset_layout_f  (coreValeBand_plid,  H5D_CHUNKED_F,hdferr)
      endif
      call h5pset_chunk_f   (valeStatesBand_plid,2,valeStatesBandChunk,hdferr)
      call h5pset_chunk_f   (statesBand_plid,1,statesBand,hdferr)
      call h5pset_chunk_f   (valeValeBand_plid,2,valeValeBandChunk,hdferr)
      if (coreDim /= 0) then
         call h5pset_chunk_f   (coreValeBand_plid,2,coreValeBandChunk,hdferr)
      endif
!      call h5pset_shuffle_f (valeStatesBand_plid,hdferr)
!      call h5pset_shuffle_f (statesBand_plid,    hdferr)
!      call h5pset_shuffle_f (valeValeBand_plid,  hdferr)
!      if (coreDim(1) /= 0) then
!         call h5pset_shuffle_f (coreValeBand_plid,  hdferr)
!      endif
      call h5pset_deflate_f (valeStatesBand_plid,1,hdferr)
      call h5pset_deflate_f (statesBand_plid,    1,hdferr)
      call h5pset_deflate_f (valeValeBand_plid,  1,hdferr)
      if (coreDim /= 0) then
         call h5pset_deflate_f (coreValeBand_plid,  1,hdferr)
      endif

      ! Allocate space to hold the IDs for the datasets in the eigen groups and
      !   overlap.
      allocate (eigenVectorsBand_did(2,numKPoints,spin))
      allocate (eigenValuesBand_did(numKPoints,spin))
      allocate (valeValeBand_did(numKPoints))
      if (coreDim /= 0) then
         allocate (coreValeBand_did(2,numKPoints))
      endif


      ! Create the datasets for the overlap coreVale and valeVale matrices.
      do i = 1, numKPoints
         if (coreDim /= 0) then
            write (currentName,fmt="(a4,i7.7)") "real",i
            currentName = trim (currentName)
            call h5dcreate_f(coreValeBand_gid,currentName,&
                  & H5T_NATIVE_DOUBLE,coreValeBand_dsid,coreValeBand_did(1,i),&
                  & hdferr,coreValeBand_plid)
            if (hdferr == -1) stop 'Can not create coreValeBand_did1'
#ifndef GAMMA
            write (currentName,fmt="(a4,i7.7)") "imag",i
            currentName = trim (currentName)
            call h5dcreate_f(coreValeBand_gid,currentName,&
                  & H5T_NATIVE_DOUBLE,coreValeBand_dsid,&
                  & coreValeBand_did(2,i),hdferr,coreValeBand_plid)
            if (hdferr == -1) stop 'Can not create coreValeBand_did2'
#endif
         endif
         write (currentName,fmt="(i7.7)") i
         currentName = trim (currentName)
         call h5dcreate_f(valeValeBand_gid,currentName,&
               & H5T_NATIVE_DOUBLE,valeValeBand_dsid,valeValeBand_did(i),&
               & hdferr,valeValeBand_plid)
            if (hdferr == -1) stop 'Can not create valeValeBand_did'
      enddo


      ! Create the datasets for the eigen values and vectors.
      do i = 1, numKPoints
         do j = 1, spin
            write (currentName,fmt="(a4,i7.7,i7.7)") "real",i,j
            currentName = trim (currentName)
            call h5dcreate_f(eigenVectorsBand_gid,currentName,&
                  & H5T_NATIVE_DOUBLE,valeStatesBand_dsid,&
                  & eigenVectorsBand_did(1,i,j),hdferr,valeStatesBand_plid)
            if (hdferr == -1) stop 'Can not create eigenVectorsBand_did1'
#ifndef GAMMA
            write (currentName,fmt="(a4,i7.7,i7.7)") "imag",i,j
            currentName = trim (currentName)
            call h5dcreate_f(eigenVectorsBand_gid,currentName,&
                  & H5T_NATIVE_DOUBLE,valeStatesBand_dsid,&
                  & eigenVectorsBand_did(2,i,j),hdferr,valeStatesBand_plid)
            if (hdferr == -1) stop 'Can not create eigenVectorsBand_did2'
#endif
            write (currentName,fmt="(i7.7,i7.7)") i,j
            currentName = trim (currentName)
            call h5dcreate_f(eigenValuesBand_gid,currentName,H5T_NATIVE_DOUBLE,&
                  & statesBand_dsid,eigenValuesBand_did(i,j),hdferr,&
                  & statesBand_plid)
            if (hdferr == -1) stop 'Can not create eigenValuesBand_did'
         enddo
      enddo
   endif

end subroutine initPSCFBandHDF5

subroutine closePSCFBandHDF5

   ! Import the necessary HDF data modules
   use HDF5

   ! Import the necessary data modules.
   use O_Potential,   only: spin
   use O_CommandLine, only: doSYBD, stateSet
   use O_KPoints,     only: numKPoints

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Define local loop control and error control variables
   integer :: i,j
   integer :: hdferr

   ! Close all the groups, datasets, dataspaces, and property lists for band
   !   only if the symmetric band calculation was not done since it does not
   !   use any of these.

   if (doSybd .eq. 0) then

      ! Close the band dataspaces first (only if this was a band calculation
      !   where data was written and not a dos, optc, bond, wave, etc type of
      !   calculation where data was only read).
      call h5sclose_f (valeStatesBand_dsid,hdferr)
      if (hdferr == -1) stop 'Can not close valeStatesBand_dsid'
      call h5sclose_f (statesBand_dsid,hdferr)
      if (hdferr == -1) stop 'Can not close statesBand_dsid'
      call h5sclose_f (valeValeBand_dsid,hdferr)
      if (hdferr == -1) stop 'Can not close valeValeBand_dsid'
      if (coreValeBand(1) /= 0) then
         call h5sclose_f (coreValeBand_dsid,hdferr)
         if (hdferr == -1) stop 'Can not close coreValeBand_dsid'
      endif

      ! Close the band datasets for the overlap coreVale and valeVale matrices.
      do i = 1, numKPoints
         if (coreValeBand(1) /= 0) then
            call h5dclose_f (coreValeBand_did(1,i),hdferr)
            if (hdferr == -1) stop 'Can not close coreValeBand_did1'
#ifndef GAMMA
            call h5dclose_f (coreValeBand_did(2,i),hdferr)
            if (hdferr == -1) stop 'Can not close coreValeBand_did2'
#endif
         endif
         call h5dclose_f (valeValeBand_did(i),hdferr)
         if (hdferr == -1) stop 'Can not close valeValeBand_did'
      enddo

      ! Close the band datasets for the eigenvalues and eigenvectors next.
      do i = 1, numKPoints
         do j = 1, spin
            call h5dclose_f (eigenVectorsBand_did(1,i,j),hdferr)
            if (hdferr == -1) stop 'Can not close eigenVectorsBand_did1'
#ifndef GAMMA
            call h5dclose_f (eigenVectorsBand_did(2,i,j),hdferr)
            if (hdferr == -1) stop 'Can not close eigenVectorsBand_did2'
#endif
            call h5dclose_f (eigenValuesBand_did(i,j),hdferr)
            if (hdferr == -1) stop 'Can not close eigenValuesBand_did'
         enddo
      enddo

      ! Close the band property lists.
      call h5pclose_f (valeStatesBand_plid,hdferr)
      call h5pclose_f (statesBand_plid,hdferr)
      call h5pclose_f (valeValeBand_plid,hdferr)
      if (coreValeBand(1) /= 0) then
         call h5pclose_f (coreValeBand_plid,hdferr)
      endif
      call h5pclose_f (band_plid,hdferr)

      ! Close the band top level groups.
      call h5gclose_f (eigenVectorsBand_gid,hdferr)
      call h5gclose_f (eigenValuesBand_gid,hdferr)
      call h5gclose_f (valeValeBand_gid,hdferr)
      if (coreValeBand(1) /= 0) then
         call h5gclose_f (coreValeBand_gid,hdferr)
      endif

      ! Close the band file.
      call h5fclose_f (band_fid,hdferr)


      ! In the case that we are closing access to the band hdf5 file from a
      !   PACS type calculation we need to also close the second file that was
      !   opened for the excited state.
      if (stateSet == 1) then

         ! Close the overlap datasets in the second band file.
         do i = 1, numKPoints
            call h5dclose_f (valeValeBand2_did(i),hdferr)
         enddo

         ! Close the wavefunction datasets in the second band file.
         do i = 1, numKPoints
            do j = 1, spin
               call h5dclose_f (eigenVectorsBand2_did(1,i,j),hdferr)
#ifndef GAMMA
               call h5dclose_f (eigenVectorsBand2_did(2,i,j),hdferr)
#endif
               call h5dclose_f (eigenValuesBand2_did(i,j),hdferr)
            enddo
         enddo

         ! Close the band top level groups.
         call h5gclose_f (eigenVectorsBand2_gid,hdferr)
         call h5gclose_f (eigenValuesBand2_gid,hdferr)
         call h5gclose_f (valeValeBand2_gid,hdferr)

         ! Close the HDF5 file access property list
         call h5pclose_f (band2_plid,hdferr)

         ! Close the HDF5 file that holds the results of the band calculation.
         call h5fclose_f (band2_fid,hdferr)
      endif
   endif

   ! Close access to the HDF interface.
   call h5close_f (hdferr)

end subroutine closePSCFBandHDF5


subroutine accessPSCFBandHDF5(numStates)

   ! Import necessary HDF5 modules.
   use HDF5

   ! Import the necessary data modules.
   use O_Potential,   only: spin
   use O_CommandLine, only: stateSet
   use O_KPoints,     only: numKPoints
   use O_AtomicSites, only: coreDim, valeDim

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Define passed parameters
   integer, intent(in) :: numStates

   ! Define the local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize the Fortran 90 HDF5 interface.
   call h5open_f(hdferr)

   ! Initialize the values of the dimension variables.
   valeStatesBand(1)    = valeDim
   valeStatesBand(2)    = numStates
#ifndef GAMMA
   valeValeBand(1)      = 2
#else
   valeValeBand(1)      = 1
#endif
   valeValeBand(2)      = valeDim*(valeDim+1)/2
   coreValeBand(1)      = coreDim
   coreValeBand(2)      = valeDim
   statesBand(1)        = numStates

   ! Create the property list for the band hdf5 file and turn off
   !   chunk caching.
   call h5pcreate_f    (H5P_FILE_ACCESS_F,band_plid,hdferr)
   call h5pget_cache_f (band_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   call h5pset_cache_f (band_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
   

   ! Open the HDF5 file that holds the results of the band calculation.
   call h5fopen_f ("band-temp.hdf5",H5F_ACC_RDONLY_F,band_fid,hdferr,&
         & band_plid)

   ! Open the top level groups in the band output file.
   call h5gopen_f (band_fid,"eigenVectors",eigenVectorsBand_gid,hdferr)
   call h5gopen_f (band_fid,"eigenValues",eigenValuesBand_gid,hdferr)
   call h5gopen_f (band_fid,"overlap",valeValeBand_gid,hdferr)
   if (coreDim /= 0) then
      call h5gopen_f (band_fid,"orthoCoeff",coreValeBand_gid,hdferr)
   endif

   ! Allocate space for the dataset ID numbers in the band output file.
   allocate (eigenValuesBand_did  (numKPoints,spin))
   allocate (eigenVectorsBand_did (2,numKPoints,spin))
   allocate (valeValeBand_did     (numKPoints))
   allocate (coreValeBand_did     (2,numKPoints))

   ! Open the datasets in the band output file for the overlap matrices.
   do i = 1, numKPoints
      if (coreDim /= 0) then
         write (currentName,fmt="(a4,i7.7)") "real",i
         currentName = trim (currentName)
         call h5dopen_f (coreValeBand_gid,currentName,coreValeBand_did(1,i),&
               & hdferr)
#ifndef GAMMA
         write (currentName,fmt="(a4,i7.7)") "imag",i
         currentName = trim (currentName)
         call h5dopen_f (coreValeBand_gid,currentName,coreValeBand_did(2,i),&
               & hdferr)
#endif
      endif
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)
      call h5dopen_f (valeValeBand_gid,currentName,valeValeBand_did(i),hdferr)
   enddo


   ! Open the datasets in the band output file for the wave functions.
   do i = 1, numKPoints
      do j = 1, spin
         write (currentName,fmt="(a4,i7.7,i7.7)") "real",i,j
         currentName = trim (currentName)
         call h5dopen_f (eigenVectorsBand_gid,currentName,&
               & eigenVectorsBand_did(1,i,j),hdferr)
#ifndef GAMMA
         write (currentName,fmt="(a4,i7.7,i7.7)") "imag",i,j
         currentName = trim (currentName)
         call h5dopen_f (eigenVectorsBand_gid,currentName,&
               & eigenVectorsBand_did(2,i,j),hdferr)
#endif
         write (currentName,fmt="(i7.7,i7.7)") i,j
         currentName = trim (currentName)
         call h5dopen_f (eigenValuesBand_gid,currentName,&
               & eigenValuesBand_did(i,j),hdferr)
      enddo
   enddo

   if (stateSet == 1) then  ! Doing a PACS type optc calculation.

      ! Create the property list for the band hdf5 file and turn off
      !   chunk caching.
      call h5pcreate_f    (H5P_FILE_ACCESS_F,band2_plid,hdferr)
      call h5pget_cache_f (band2_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,&
            & rdcc_w0,hdferr)
      call h5pset_cache_f (band2_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,&
            & hdferr)
   
      ! Open the second band calculation file.
      call h5fopen_f ("band-temp.2.hdf5",H5F_ACC_RDONLY_F,band2_fid,hdferr,&
         & band2_plid)

      ! Open the top level groups of the second band calculation file.
      call h5gopen_f (band2_fid,"eigenVectors",eigenVectorsBand2_gid,hdferr)
      call h5gopen_f (band2_fid,"eigenValues",eigenValuesBand2_gid,hdferr)
      call h5gopen_f (band2_fid,"overlap",valeValeBand2_gid,hdferr)

      ! Allocate space for the dataset ID numbers in the band output file.
      allocate (eigenVectorsBand2_did (2,numKPoints,spin))
      allocate (eigenValuesBand2_did  (numKPoints,spin))
      allocate (valeValeBand2_did     (numKPoints))


      ! Open the datasets in the second band output file for the overlap.
      do i = 1, numKPoints
         write (currentName,fmt="(i7.7)") i
         currentName = trim (currentName)
         call h5dopen_f (valeValeBand2_gid,currentName,&
               & valeValeBand2_did(i),hdferr)
      enddo


      ! Open the datasets in the second band output file for the wave functions.
      do i = 1, numKPoints
         do j = 1, spin
            write (currentName,fmt="(a4,i7.7,i7.7)") "real",i,j
            currentName = trim (currentName)
            call h5dopen_f (eigenVectorsBand2_gid,currentName,&
                  & eigenVectorsBand2_did(1,i,j),hdferr)
#ifndef GAMMA
            write (currentName,fmt="(a4,i7.7,i7.7)") "imag",i,j
            currentName = trim (currentName)
            call h5dopen_f (eigenVectorsBand2_gid,currentName,&
                  & eigenVectorsBand2_did(2,i,j),hdferr)
#endif
            write (currentName,fmt="(i7.7,i7.7)") i,j
            currentName = trim (currentName)
            call h5dopen_f (eigenValuesBand2_gid,currentName,&
                  & eigenValuesBand2_did(i,j),hdferr)
         enddo
      enddo
   endif
end subroutine accessPSCFBandHDF5


subroutine closeAccessPSCFBandHDF5

   ! Import the necessary HDF data modules
   use HDF5

   ! Import the necessary data modules.
   use O_Potential,   only: spin
   use O_CommandLine, only: stateSet
   use O_KPoints,     only: numKPoints

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Define local loop control and error control variables
   integer :: i,j
   integer :: hdferr

   ! Close all the groups, datasets, dataspaces, and property lists for band.

   ! Close the band datasets for the overlap coreVale and valeVale matrices.
   do i = 1, numKPoints
      if (coreValeBand(1) /= 0) then
         call h5dclose_f (coreValeBand_did(1,i),hdferr)
         if (hdferr == -1) stop 'Can not close coreValeBand_did1'
#ifndef GAMMA
         call h5dclose_f (coreValeBand_did(2,i),hdferr)
         if (hdferr == -1) stop 'Can not close coreValeBand_did2'
#endif
      endif
      call h5dclose_f (valeValeBand_did(i),hdferr)
      if (hdferr == -1) stop 'Can not close valeValeBand_did'
   enddo

   ! Close the band datasets for the eigenvalues and eigenvectors next.
   do i = 1, numKPoints
      do j = 1, spin
         call h5dclose_f (eigenVectorsBand_did(1,i,j),hdferr)
         if (hdferr == -1) stop 'Can not close eigenVectorsBand_did1'
#ifndef GAMMA
         call h5dclose_f (eigenVectorsBand_did(2,i,j),hdferr)
         if (hdferr == -1) stop 'Can not close eigenVectorsBand_did2'
#endif
         call h5dclose_f (eigenValuesBand_did(i,j),hdferr)
         if (hdferr == -1) stop 'Can not close eigenValuesBand_did'
      enddo
   enddo

   ! Close the band property lists.
   call h5pclose_f (valeStatesBand_plid,hdferr)
   if (hdferr == -1) stop 'Can not close valeStatesBand_plid'
   call h5pclose_f (statesBand_plid,hdferr)
   if (hdferr == -1) stop 'Can not close statesBand_plid'
   call h5pclose_f (valeValeBand_plid,hdferr)
   if (hdferr == -1) stop 'Can not close valeValeBand_plid'
   if (coreValeBand(1) /= 0) then
      call h5pclose_f (coreValeBand_plid,hdferr)
      if (hdferr == -1) stop 'Can not close coreValeBand_plid'
   endif
   call h5pclose_f (band_plid,hdferr)
   if (hdferr == -1) stop 'Can not close band_plid'

   ! Close the band top level groups.
   call h5gclose_f (eigenVectorsBand_gid,hdferr)
   if (hdferr == -1) stop 'Can not close eigenVectorsBand_gid'
   call h5gclose_f (eigenValuesBand_gid,hdferr)
   if (hdferr == -1) stop 'Can not close eigenValuesBand_gid'
   call h5gclose_f (valeValeBand_gid,hdferr)
   if (hdferr == -1) stop 'Can not close valeValeBand_gid'
   if (coreValeBand(1) /= 0) then
      call h5gclose_f (coreValeBand_gid,hdferr)
      if (hdferr == -1) stop 'Can not close coreValeBand_gid'
   endif

   ! Close the band file.
   call h5fclose_f (band_fid,hdferr)
   if (hdferr == -1) stop 'Can not close band_fid'


   ! In the case that we are closing access to the band hdf5 file from a
   !   PACS type calculation we need to also close the second file that was
   !   opened for the excited state.
   if (stateSet == 1) then

      ! Close the overlap datasets in the second band file.
      do i = 1, numKPoints
         call h5dclose_f (valeValeBand2_did(i),hdferr)
         if (hdferr == -1) stop 'Can not close valeValeBand2_did'
      enddo

      ! Close the wavefunction datasets in the second band file.
      do i = 1, numKPoints
         do j = 1, spin
            call h5dclose_f (eigenVectorsBand2_did(1,i,j),hdferr)
            if (hdferr == -1) stop 'Can not close eigenVectorsBand2_did1'
#ifndef GAMMA
            call h5dclose_f (eigenVectorsBand2_did(2,i,j),hdferr)
            if (hdferr == -1) stop 'Can not close eigenVectorsBand2_did2'
#endif
            call h5dclose_f (eigenValuesBand2_did(i,j),hdferr)
            if (hdferr == -1) stop 'Can not close eigenValuesBand2_did'
         enddo
      enddo

      ! Close the band top level groups.
      call h5gclose_f (eigenVectorsBand2_gid,hdferr)
      if (hdferr == -1) stop 'Can not close eigenVectorsBand2_gid'
      call h5gclose_f (eigenValuesBand2_gid,hdferr)
      if (hdferr == -1) stop 'Can not close eigenValuesBand2_gid'
      call h5gclose_f (valeValeBand2_gid,hdferr)
      if (hdferr == -1) stop 'Can not close valeValeBand2_gid'

      ! Close the HDF5 file access property list
      call h5pclose_f (band2_plid,hdferr)
      if (hdferr == -1) stop 'Can not close band2_plid'

      ! Close the HDF5 file that holds the results of the band calculation.
      call h5fclose_f (band2_fid,hdferr)
      if (hdferr == -1) stop 'Can not close band2_fid'
   endif

   ! Close access to the HDF interface.
   call h5close_f (hdferr)
   if (hdferr == -1) stop 'Can not close HDF5 interface'

end subroutine closeAccessPSCFBandHDF5


#ifndef GAMMA
subroutine saveCoreValeOL (coreValeOL,i)
#else
subroutine saveCoreValeOL (coreValeOLGamma,i)
#endif

   ! Use necessary modules.
   use HDF5
   use O_Kinds
   use O_CommandLine, only: doSYBD
   use O_AtomicSites, only: coreDim, valeDim

   ! Define passed dummy parameters.
#ifndef GAMMA
   complex (kind=double), dimension (:,:,:) :: coreValeOL
#else
   real (kind=double), dimension (:,:)      :: coreValeOLGamma
#endif
   integer :: i ! Current KPoint index number.
   
   ! Define local variables.
   integer :: hdferr

   if ((doSYBD == 0) .and. (coreDim /= 0)) then
#ifndef GAMMA
      call h5dwrite_f (coreValeBand_did(1,i),H5T_NATIVE_DOUBLE,&
            & real(coreValeOL(:coreDim,:valeDim,1),double),coreValeBand,hdferr)
      if (hdferr /= 0) stop 'Can not write real core vale band.'
      call h5dwrite_f (coreValeBand_did(2,i),H5T_NATIVE_DOUBLE,&
            & aimag(coreValeOL(:coreDim,:valeDim,1)),coreValeBand,hdferr)
      if (hdferr /= 0) stop 'Can not write imag core vale band.'
#else
      call h5dwrite_f (coreValeBand_did(1,i),H5T_NATIVE_DOUBLE,&
            & coreValeOLGamma(:coreDim,:ValeDim),coreValeBand,hdferr)
      if (hdferr /= 0) stop 'Can not write real core vale band.'
#endif
   endif

end subroutine saveCoreValeOL



end module O_PSCFBandHDF5
