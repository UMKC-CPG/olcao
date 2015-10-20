module O_PSCFIntgHDF5

   ! Use the HDF5 module for HDF5 defined types (e.g. size_t and hid_t).
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define variables that hold information used to construct the HDF5 file.

   ! A target value for the size of a chunk to record.  This is a lower limit
   !   on the data to hold in memory.  Once this limit is surpassed, then the
   !   data is recorded to disk.
   integer(hsize_t) :: targetChunkSize

   ! A number that limits the memory footprint of the accumulated integral
   !   data.  It will be the maximum number of floating point numbers
   !   to be stored.
   integer(hsize_t) :: maxMemory


   ! Define arrays that hold the dimensions of the datasets.
   integer(hsize_t), dimension (1) :: attribIntDims
   integer(hsize_t), dimension (2) :: totalDataDims
   integer(hsize_t), dimension (2) :: totalLoopIndexDims
   integer(hsize_t), dimension (2) :: loopIndexDims
   integer(hsize_t), dimension (2) :: chunkDims
   integer(hsize_t), dimension (2) :: dataChunkSize
   integer(hsize_t), dimension (2) :: loopIndexChunkSize

   ! Define arrays that hold the starting positions for reading the file data.
   integer(hsize_t), dimension (2) :: currentStart
   integer(hsize_t), dimension (2) :: currentIndexStart


   ! Define the HDF5 file's structure and components.


   ! Define the file ID
   integer(hid_t) :: intg_fid

   ! Define the input file access property list ID. 
   integer(hid_t)   :: intg_plid
   integer          :: mdc_nelmts  ! Meta-data cache num elements.
   integer (size_t) :: rdcc_nelmts ! Raw-data chunk cache num elements.
   integer (size_t) :: rdcc_nbytes ! Raw-data chunk cache size in bytes.
   real             :: rdcc_w0     ! Raw-data chunk cache weighting parameter.

   ! Define the group IDs under intg_fid
   integer(hid_t) :: overlap_gid
   integer(hid_t) :: hamiltonian_gid
   integer(hid_t) :: momentumX_gid
   integer(hid_t) :: momentumY_gid
   integer(hid_t) :: momentumZ_gid
   integer(hid_t) :: loopIndices_gid
   integer(hid_t) :: intg_gid ! Generic group ID used for reading.

   ! This dataspace will be used for single value integer attributes.
   integer(hid_t) :: attribInt_dsid

   ! Define an attribute ID for the momentumX_gid object that will say whether
   !   it has been computed (1) or not (0).
   integer(hid_t) :: momentumX_aid

   ! This dataspace will be used to reference the memory data when writing.  It
   !   will be the size of the data chunk for the current i loop atom.
   integer(hid_t) :: dataChunk_dsid

   ! This dataspace will be used to track the loop indices in memory just like
   !   the above dataChunk_dsid.
   integer(hid_t) :: loopIndexChunk_dsid

   ! These dataspaces will be used to write to datasets on file.  (These are
   !   the dataspaces on the file, not it memory.)
   integer(hid_t), dimension(2) :: fileHamiltonian_dsid
   integer(hid_t) :: fileOverlap_dsid
   integer(hid_t) :: fileMomentumX_dsid
   integer(hid_t) :: fileMomentumY_dsid
   integer(hid_t) :: fileMomentumZ_dsid
   integer(hid_t) :: fileLoopIndices_dsid
   integer(hid_t) :: fileIntg_dsid  !  Generic_dsid for reading data.
   integer(hid_t) :: loopIndices_dsid  ! Dataspace for loops in memory.

   ! The property lists for the file representations of the loop indices and
   !   the data will be made and destroyed with each i loop iteration.  These
   !   are used for writing the data.
   integer(hid_t) :: fileLoopIndices_plid
   integer(hid_t) :: fileData_plid
   integer(hid_t) :: fileIntg_plid  !  Generic _plid for reading data.

   ! For each iteration of the intgAndMom loop there will need to be a dataset
   !   id number for the overlap, hamiltonian, and momentum matrices.  The
   !   number will be released after each write so we only need one at a time.
   integer(hid_t), dimension(2) :: hamiltonian_did
   integer(hid_t) :: overlap_did
   integer(hid_t) :: momentumX_did
   integer(hid_t) :: momentumY_did
   integer(hid_t) :: momentumZ_did
   integer(hid_t) :: intg_did      ! Generic _did for reading the data.

   ! Also, for each iteration of the atom 1 loop there will need to be a
   !   dataset that keeps track of the loop indices for later reading.
   integer(hid_t) :: loopIndices_did

contains

subroutine initPSCFIntgHDF5

   ! Import necessary modules.
   use O_Kinds
   use HDF5

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Define local variables that will be used to create the dataspaces etc.
   integer :: hdferr

   ! Initialize the Fortran 90 HDF5 interface.
   call h5open_f(hdferr)
   if (hdferr /= 0) stop 'Failed to open HDF5 interface'

   ! Define the memory foot print parameter.  This controls the maximum amount
   !   of memory that is allocated in the post SCF integral calculation.
   maxMemory = 10000000

   
   ! Create the property list for the intg hdf5 file and turn off
   !   chunk caching.
   call h5pcreate_f    (H5P_FILE_ACCESS_F,intg_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create property list for intg file.'

   call h5pget_cache_f (intg_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to get default cache settings for intg file.'

   call h5pset_cache_f (intg_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to set cache settings for intg file.'


   ! Create the HDF5 file that will hold the computed intg results.  This will
   !   die if the file already exists.  It also uses the default file creation
   !   and file access properties.
   call h5fcreate_f ("intg-temp.hdf5",H5F_ACC_EXCL_F,intg_fid,hdferr,&
         & H5P_DEFAULT_F,intg_plid)
   if (hdferr /= 0) stop 'Failed to create PSCF intg hdf5 file.'


   ! Then create all the groups.  (The momentum matrix elements are only
   !   computed if it is requested, but the group is always made so that the
   !   attribute tag can mark its status as computed or not.)
   call h5gcreate_f (intg_fid,"overlap",overlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create overlap_gid'
   call h5gcreate_f (intg_fid,"hamiltonian",hamiltonian_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create hamiltonian_gid'
   call h5gcreate_f (intg_fid,"momentumX",momentumX_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create momentumX_gid'
   call h5gcreate_f (intg_fid,"momentumY",momentumY_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create momentumY_gid'
   call h5gcreate_f (intg_fid,"momentumZ",momentumZ_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create momentumZ_gid'
   call h5gcreate_f (intg_fid,"loopIndices",loopIndices_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create loopIndices_gid'

   ! Set the dimension of the attribute dataspace.
   attribIntDims(1) = 1

   ! Create the dataspace for the momentumX_gid computed status attribute.
   call h5screate_simple_f(1,attribIntDims,attribInt_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create attribInt_dsid.'

   ! Create the computed status attribute for the momentumX_gid.
   call h5acreate_f (momentumX_gid,"computed",H5T_NATIVE_INTEGER,&
         & attribInt_dsid,momentumX_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create momentumX_aid.'


end subroutine initPSCFIntgHDF5


subroutine openMOMEPSCFIntgHDF5

   ! Import the necessary HDF modules
   use HDF5

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Define the error control variable.
   integer :: hdferr

   ! Define the memory foot print parameter.  This controls the maximum amount
   !   of memory that is allocated in the post SCF integral calculation.
   maxMemory = 10000000

   ! Create the property list for the intg hdf5 file and turn off
   !   chunk caching.
   call h5pcreate_f    (H5P_FILE_ACCESS_F,intg_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create property list for intg file.'
   call h5pget_cache_f (intg_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to get default cache settings for intg file.'
   call h5pset_cache_f (intg_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to set cache settings for intg file.'

   ! Open the HDF5 file that holds the computed intg results.
   call h5fopen_f ("intg-temp.hdf5",H5F_ACC_RDWR_F,intg_fid,hdferr,intg_plid)
   if (hdferr /= 0) stop 'Failed to open intg-temp.hdf5'

   ! Open the momentumX,Y,Z groups.
   call h5gopen_f (intg_fid,'momentumX',momentumX_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open momentumX_gid'
   call h5gopen_f (intg_fid,'momentumY',momentumY_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open momentumY_gid'
   call h5gopen_f (intg_fid,'momentumZ',momentumZ_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open momentumZ_gid'

   ! Open the loop indices group.
   call h5gopen_f (intg_fid,"/loopIndices",loopIndices_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open loopIndices_gid.'

   ! Open the 'computed' attribute by name in case it needs to be set.
   call h5aopen_f (momentumX_gid,'computed',momentumX_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open momentumX_aid'

end subroutine openMOMEPSCFIntgHDF5


subroutine closeMOMEPSCFIntgHDF5

   ! Import the necessary HDF modules
   use HDF5

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Define the error control variable.
   integer :: hdferr

   ! Close the momentumX,Y,Z groups.
   call h5gclose_f (momentumX_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the momentumX_gid.'
   call h5gclose_f (momentumY_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the momentumY_gid.'
   call h5gclose_f (momentumZ_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the momentumZ_gid.'

   ! Close the loopIndices group.
   call h5gclose_f (loopIndices_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the loopIndices_gid.'

   ! Close the property list.
   call h5pclose_f (intg_plid,hdferr)

   ! Close the HDF5 MOME attribute.
   call h5aclose_f (momentumX_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to close momentumX_aid'

   ! Close the access to the integral file.
   call h5fclose_f (intg_fid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the intg_fid.'

   ! Close access to the HDF5 interface.
   call h5close_f (hdferr)
   if (hdferr /= 0) stop 'Failed to close the HDF5 interface.'

end subroutine closeMOMEPSCFIntgHDF5


subroutine closePSCFIntgHDF5

   ! Import the necessary HDF modules
   use HDF5

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Define the error control variable.
   integer :: hdferr

   ! Close all the groups, datasets, dataspaces, and property lists for intg.

   ! The dataspaces are mostly closed after the end of each outer atom loop.
   !   The exception is the attribute dataspace for the momentumX_gid that says
   !   whether or not the momentum matrix elements were calculated.  This is
   !   closed now.
   call h5sclose_f(attribInt_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the attribInt_dsid.'

   ! The datasets are closed at the end of each outer atom loop.

   ! The property lists are close at the end of each outer atom loop.

   ! Close the attributes of the groups.
   call h5aclose_f (momentumX_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the momentumX_aid.'

   ! Close the groups.
   call h5gclose_f (overlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the overlap_gid.'
   call h5gclose_f (hamiltonian_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the hamiltonian_gid.'
   call h5gclose_f (momentumX_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the momentumX_gid.'
   call h5gclose_f (momentumY_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the momentumY_gid.'
   call h5gclose_f (momentumZ_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the momentumZ_gid.'
   call h5gclose_f (loopIndices_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the loopIndices_gid.'

   ! Close the file access property list
   call h5pclose_f (intg_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the intg_plid.'

   ! Close the file.
   call h5fclose_f (intg_fid,hdferr)
   if (hdferr /= 0) stop 'Failed to close the intg_fid.'

   ! Close access to the HDF interface.
   call h5close_f (hdferr)
   if (hdferr /= 0) stop 'Failed to close the HDF5 interface.'

end subroutine closePSCFIntgHDF5


subroutine setMOMEStatus(doMOME)

   ! Use necessary modules.
   use HDF5

   ! Make sure no variables are accidentally defined.
   implicit none

   ! Define passed parameter.
   integer, intent (IN) :: doMOME

   ! Define local variables.
   integer :: hdferr

   ! Write the given value for whether or not the momentum matrix elements
   !   have been calculated.
   call h5awrite_f(momentumX_aid,H5T_NATIVE_INTEGER,doMOME,attribIntDims,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to write momentumX_aid'

end subroutine setMOMEStatus


subroutine getMOMEStatus(doneMOME)

   ! Use necessary modules.
   use HDF5

   ! Make sure no variables are accidentally defined.
   implicit none

   ! Define passed parameter.
   integer, intent (OUT) :: doneMOME

   ! Define local variables.
   integer :: hdferr

   ! Set the dimension of the attribute dataspace.
   attribIntDims(1) = 1

   ! Gain access to the HDF5 interface.
   call h5open_f (hdferr)
   if (hdferr /= 0) stop 'Failed to open HDF5 interface'

   ! Open the HDF5 file that holds the computed intg results.
   call h5fopen_f ("intg-temp.hdf5",H5F_ACC_RDWR_F,intg_fid,hdferr,intg_plid)
   if (hdferr /= 0) stop 'Failed to open intg-temp.hdf5'

   ! Open the momentumX group.
   call h5gopen_f (intg_fid,'momentumX',momentumX_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open momentumX_gid'

   ! Open the 'computed' attribute by name to learn if the momentum matrix
   !   elements have been computed.
   call h5aopen_f (momentumX_gid,'computed',momentumX_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to open momentumX_aid'

   ! Read the variable that indicates whether or not the momentum matrix
   !   elements have been calculated.
   call h5aread_f(momentumX_aid,H5T_NATIVE_INTEGER,doneMOME,attribIntDims,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to read momentumX_aid'

   ! Close the HDF5 MOME attribute.
   call h5aclose_f (momentumX_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to close momentumX_aid'

   ! Close the momentumX group.
   call h5gclose_f (momentumX_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close momentumX_gid'

   ! Close the integral file.
   call h5fclose_f (intg_fid,hdferr)
   if (hdferr /= 0) stop 'Failed to close intg_fid'

end subroutine getMOMEStatus


subroutine initReadIntgHDF5 (runCode)

   ! Import the necessary HDF data modules
   use HDF5

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Define the dummy variables.
   integer :: runCode

   ! Define the local error control variable.
   integer :: hdferr

   ! Assign the values for the intg input.
   loopIndexDims(1) = 2
   loopIndexDims(2) = 0

   ! Create the property list for the intg hdf5 file and turn off
   !   chunk caching.
   call h5pcreate_f    (H5P_FILE_ACCESS_F,intg_plid,hdferr)
   call h5pget_cache_f (intg_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   call h5pset_cache_f (intg_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)

   ! Open the HDF5 file that holds the results of the intg calculation.
   call h5fopen_f ("intg-temp.hdf5",H5F_ACC_RDONLY_F,intg_fid,hdferr,intg_plid)


   ! Open the top level group in the HDF5 file that hold the individual
   !   calculated components from intg.
   select case (runCode)
   case (1)
      call h5gopen_f (intg_fid,"/overlap",intg_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to open overlap group.'
   case (2)
      call h5gopen_f (intg_fid,"/hamiltonian",intg_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to open hamiltonian group.'
   case (3)
      call h5gopen_f (intg_fid,"/momentumX",intg_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to open momentumX group.'
   case (4)
      call h5gopen_f (intg_fid,"/momentumY",intg_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to open momentumY group.'
   case (5)
      call h5gopen_f (intg_fid,"/momentumZ",intg_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to open momentumZ group.'
   end select


   ! There is one dataset in the above groups for each atom.  They will be
   !   opened and closed in the appropriate loops.  They hold the data from
   !   the computed integrals (overlap, hamiltonian).

   ! Open the loop indices group.
   call h5gopen_f (intg_fid,"/loopIndices",loopIndices_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open loopIndices group.'

   ! There is one dataset in the loopIndices group for each atom.  They
   !   will be opened and closed in the appropriate loops.  This collection of
   !   datasets tracks the loop indices.

   ! The data spaces for the data and loop indices will be defined for each
   !   i loop atom iteration.

end subroutine initReadIntgHDF5


end module O_PSCFIntgHDF5
