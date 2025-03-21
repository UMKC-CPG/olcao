module O_SCFExchCorrHDF5

   ! Import necessary modules.
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define the main subgroup from scf_fid that holds the mesh information
   !   for the exchange correlation calculation.
   integer(hid_t) :: exchCorrGroup_gid ! Exchange Correlation group.

   ! Begin group, dataspace, dataset definitions for exchCorrGroup_gid.

   ! There are no subgroups of the exchCorrGroup, only datasets.

   ! The dataspace of each dataset in the exchCorrGroup will be one of
   !   three possible kinds.  Only two are listed because the third is the
   !   same as the one in the elecStatGroup, potPot_dsid.  The "points"
   !   here is the max number of ray points.
   integer(hid_t) :: numPoints_dsid
   integer(hid_t) :: potPoints_dsid
   integer(hid_t) :: potPot_dsid
   integer(hid_t) :: points_dsid

   ! Each of the below datasets will use one of these property lists.  They
   !   correspond to the above dataspaces.  The potPot_plid is not listed here,
   !   it is given in the elecStat module and imported when needed.
   integer(hid_t) :: potPoints_plid
   integer(hid_t) :: potPot_plid
   integer(hid_t) :: points_plid

   ! Define arrays that hold the dimensions of the datasets.  The potDims2
   !   dimension set is not given here because it is predominantly used in the
   !   elecStat module.  Import it from there when needed.
   integer(hsize_t), dimension (3) :: potPoints
   integer(hsize_t), dimension (2) :: potPot
   integer(hsize_t), dimension (1) :: points
   integer(hsize_t), dimension (1) :: numPoints

   ! The number of datasets under exchCorrGroup_gid will also vary depending
   !   on the nature of the problem (potential dimension).  Those datasets
   !   will be given IDs dynamically.  There is also one more dataset for the
   !   exchange correclation overlap.
   integer(hid_t), allocatable, dimension (:) :: numPoints_did
   integer(hid_t), allocatable, dimension (:) :: exchRhoOp_did
   integer(hid_t), allocatable, dimension (:) :: radialWeight_did
   integer(hid_t) :: exchCorrOverlap_did

   ! Define the attribute ID that will be used for tracking completion.
   integer(hid_t) :: exchCorrGroup_aid

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initSCFExchCorrHDF5 (scf_fid,attribInt_dsid,attribIntDims,&
      & maxNumRayPoints)

   ! Import any necessary definition modules.
   use HDF5
   use O_MPI
   use O_SCFElecStatHDF5, only: potPot_plid, potPot_dsid

   ! Import necessary object modules.
   use O_PotSites, only: numPotSites
   use O_Potential, only: potDim, GGA

   ! Define the passed parameters.
   integer(hid_t), intent(in) :: scf_fid
   integer(hid_t) :: attribInt_dsid
   integer(hsize_t), dimension (1) :: attribIntDims
   integer, intent(in) :: maxNumRayPoints

   ! Define local variables.
   integer :: i
   integer :: hdferr
   character*30 :: currentName

   ! Initialize data structure dimensions.
   potPoints(1)    = potDim
   potPoints(2)    = maxNumRayPoints
   ! For nonGGA only the rho op value is stored. For GGA the first and second
   !   derivatives, as well as the rho op value are stored. Thus, for GGA 10
   !   indices are required instead of 1.
   if (GGA == 0) then
      potPoints(3) = 1
   else
      potPoints(3) = 10
   endif
   potPot(1)       = potDim
   potPot(2)       = potDim
   points(1)       = maxNumRayPoints
   numPoints(1)    = 1

   ! Allocate space to hold the IDs for the datasets in the exchCorrGroup.
   allocate (numPoints_did(numPotSites))
   allocate (exchRhoOp_did(numPotSites))
   allocate (radialWeight_did(numPotSites))

   ! Only process 0 opens the HDF5 structure.
   if (mpiRank /= 0) return

   ! Create the exchange correlation group within the scf_fid.
   call h5gcreate_f (scf_fid,"/exchCorrGroup",exchCorrGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to create exchange correlation group'

   ! Begin group, dataspace, dataset definitions for exchCorrGroup_gid.
   
   ! Create the dataspaces that will be used for each dataset in the
   !   exchCorrGroup.
   call h5screate_simple_f(1,numPoints,numPoints_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create numPoints dsid.'
   call h5screate_simple_f(3,potPoints,potPoints_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create potPoints dsid.'
   call h5screate_simple_f(1,points,points_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to create points dsid.'

   ! Create the property list first.  Then set the properties one at a time.
   call h5pcreate_f(H5P_DATASET_CREATE_F,potPoints_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create potPoints plid.'
   call h5pcreate_f(H5P_DATASET_CREATE_F,points_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create points plid.'
   call h5pset_layout_f(potPoints_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set potPoints chunked layout.'
   call h5pset_layout_f(points_plid,H5D_CHUNKED_F,hdferr)
   if (hdferr /= 0) stop 'Failed to set points chunked layout.'
   call h5pset_chunk_f(potPoints_plid,3,potPoints,hdferr)
   if (hdferr /= 0) stop 'Failed to set potPoints chunked property.'
   call h5pset_chunk_f(points_plid,1,points,hdferr)
   if (hdferr /= 0) stop 'Failed to set points chunked property.'
!   call h5pset_shuffle_f(potPoints_plid,hdferr)
!   call h5pset_shuffle_f(points_plid,hdferr)
   call h5pset_deflate_f   (potPoints_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set potPoints deflate property.'
   call h5pset_deflate_f   (points_plid,1,hdferr)
   if (hdferr /= 0) stop 'Failed to set points deflate property.'

   ! Create the datasets that will be used for the exchCorrGroup.
   call h5dcreate_f(exchCorrGroup_gid,"exchCorrOverlap",H5T_NATIVE_DOUBLE,&
         & potPot_dsid,exchCorrOverlap_did,hdferr,potPot_plid)
   if (hdferr /= 0) stop 'Failed to create the exchCorrOverlap did.'
   do i = 1, numPotSites
      write (currentName,*) i,"nP"
      currentName = trim (currentName)
      call h5dcreate_f(exchCorrGroup_gid,currentName,H5T_NATIVE_INTEGER,&
            & numPoints_dsid,numPoints_did(i),hdferr) ! No plid necessary.
            !  The value stored is a single number.
      if (hdferr /= 0) stop 'Failed to create the numPoints did.'
      write (currentName,*) i,"rW"
      currentName = trim (currentName)
      call h5dcreate_f(exchCorrGroup_gid,currentName,H5T_NATIVE_DOUBLE,&
            & points_dsid,radialWeight_did(i),hdferr,points_plid)
      if (hdferr /= 0) stop 'Failed to create the radialWeight did.'
      write (currentName,*) i,"eRO"
      currentName = trim (currentName)
      call h5dcreate_f(exchCorrGroup_gid,currentName,H5T_NATIVE_DOUBLE,&
            & potPoints_dsid,exchRhoOp_did(i),hdferr,potPoints_plid)
      if (hdferr /= 0) stop 'Failed to create the exchRhoOp did.'
   enddo

   ! Create the attribute that will indicate the completion.
   call h5acreate_f (exchCorrGroup_gid,"status",H5T_NATIVE_INTEGER,&
         & attribInt_dsid,exchCorrGroup_aid,hdferr)
   if (hdferr /= 0) stop 'Failed to create exchCorrGroup aid'

   ! Initialize the status to incomplete (zero).
   call h5awrite_f (exchCorrGroup_aid,H5T_NATIVE_INTEGER,0,attribIntDims,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to initialize exchCorrGroup aid'


end subroutine initSCFExchCorrHDF5



subroutine accessSCFExchCorrHDF5 (scf_fid)

   ! Import any necessary definition modules.
   use HDF5
   use O_MPI

   ! Import necessary object modules.
   use O_PotSites, only: numPotSites
   use O_Potential, only: potDim, GGA

   ! Define the passed parameters.
   integer(hid_t), intent(in) :: scf_fid

   ! Define local variables.
   integer :: i
   integer :: hdferr
   integer :: tempNumRayPoints
   integer :: tempMaxNumRayPoints
   character*30 :: currentName

   ! Initialize the tempMaxNumRayPoints to a small number.
   tempMaxNumRayPoints = 0

   ! Initialize data structure dimensions.  Note that on access, we do not
   !   already have "maxNumRayPoints" computed.  Here we will simply read the
   !   dataset values to determine the array dimensions.  Therefore, the array
   !   dimension variables are defined at the end of this subroutine instead of
   !   at the beginning as is done for other similar subroutines.  (Except for
   !   numPoints because we need it to find out the others.  There's always
   !   something isn't there.)
   numPoints(1) = 1

   ! Allocate space to hold the IDs for the datasets in the exchCorrGroup.
   allocate (numPoints_did(numPotSites))
   allocate (exchRhoOp_did(numPotSites))
   allocate (radialWeight_did(numPotSites))

   if (mpiRank == 0) then

      ! Open the exchange correlation group within the scf_fid.
      call h5gopen_f (scf_fid,"/exchCorrGroup",exchCorrGroup_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to open exchange correlation group'

      ! Open the datasets that will be used for the exchCorrGroup.
      call h5dopen_f(exchCorrGroup_gid,"exchCorrOverlap",exchCorrOverlap_did,&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to open the exchCorrOverlap did.'
      do i = 1, numPotSites
         write (currentName,*) i,"nP"
         currentName = trim (currentName)
         call h5dopen_f(exchCorrGroup_gid,currentName,numPoints_did(i),hdferr)
         if (hdferr /= 0) stop 'Failed to open the numPoints did.'
         write (currentName,*) i,"rW"
         currentName = trim (currentName)
         call h5dopen_f(exchCorrGroup_gid,currentName,radialWeight_did(i),&
               & hdferr)
         if (hdferr /= 0) stop 'Failed to open the radialWeight did.'
         write (currentName,*) i,"eRO"
         currentName = trim (currentName)
         call h5dopen_f(exchCorrGroup_gid,currentName,exchRhoOp_did(i),hdferr)
         if (hdferr /= 0) stop 'Failed to open the exchRhoOp did.'
      enddo

      ! Obtain the property lists.  Read note for integrals and plid numbers in
      !   accessSetupIntegralHDF5.
      call h5dget_create_plist_f (exchRhoOp_did(1),potPoints_plid,hdferr)
      if (hdferr /= 0) stop 'Failed to obtain pot points plid'
      call h5dget_create_plist_f (radialWeight_did(1),points_plid,hdferr)
      if (hdferr /= 0) stop 'Failed to obtain points plid'

      ! Obtain the dataspaces that will be used for each dataset in the
      !   exchCorrGroup.  Read note for integrals and dsid numbers in
      !   accessSetupIntegralHDF5.
      call h5dget_space_f(numPoints_did(1),numPoints_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to obtain numPoints dsid.'
      call h5dget_space_f(exchRhoOp_did(1),potPoints_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to obtain potPoints dsid.'
      call h5dget_space_f(radialWeight_did(1),points_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to obtain points dsid.'

      ! Obtain the maximum number of ray points of all potential sites in the
      !   system in order to define the size of the stored matrices.
      do i = 1, numPotSites

         ! Num ray points for this potSite.
         call h5dread_f (numPoints_did(i),H5T_NATIVE_INTEGER,&
               & tempNumRayPoints,numPoints,hdferr)

         ! Max of all read in.
         tempMaxNumRayPoints = max(tempMaxNumRayPoints,tempNumRayPoints)
      enddo
   endif
   call MPI_BCAST(tempMaxNumRayPoints,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

   ! Define the dimensions of the remaining HDF5 stored matrices.
   potPoints(1) = potDim
   potPoints(2) = tempMaxNumRayPoints
   if (GGA == 0) then
      potPoints(3) = 1
   else
      potPoints(3) = 10
   endif
   points(1) = tempMaxNumRayPoints

   ! Open the attribute that records the completion status.
   if (mpiRank == 0) then
      call h5aopen_f(exchCorrGroup_gid,"status",exchCorrGroup_aid,hdferr)
      if (hdferr /= 0) stop 'Failed to open exchCorrGroup_aid'
   endif

end subroutine accessSCFExchCorrHDF5



subroutine closeSCFExchCorrHDF5

   ! Import any necessary definition modules.
   use HDF5
   use O_MPI

   ! Import necessary object modules.
   use O_PotSites, only: numPotSites

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Declare local loop variables and error control.
   integer :: i
   integer :: hdferr

   ! Close all access to the exchange correlation parts of the HDF file
   !   for setup.
   if (mpiRank == 0) then

      ! Close the property lists first.
      call h5pclose_f (potPoints_plid,hdferr)
      if (hdferr /= 0) stop 'Failed to close potPoints_plid.'
      call h5pclose_f (points_plid,hdferr)
      if (hdferr /= 0) stop 'Failed to close points_plid.'

      ! Close the datasets next.
      do i = 1, numPotSites
         call h5dclose_f (numPoints_did(i),hdferr)
         if (hdferr /= 0) stop 'Failed to close numPoints_did.'
         call h5dclose_f (exchRhoOp_did(i),hdferr)
         if (hdferr /= 0) stop 'Failed to close exchRhoOp_did.'
         call h5dclose_f (radialWeight_did(i),hdferr)
         if (hdferr /= 0) stop 'Failed to close radialWeight_did.'
      enddo

      ! Attributes are closed when checked or written.

      ! Close the data spaces next.
      call h5sclose_f (potPoints_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to close potPoints_dsid.'
      call h5sclose_f (points_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to close points_dsid.'

      ! Close the groups.
      call h5gclose_f (exchCorrGroup_gid,hdferr)
      if (hdferr /= 0) stop 'Failed to close exchCorrGroup_gid.'
   endif

   ! Deallocate space that holds the IDs for the datasets in the exchCorrGroup.
   deallocate (numPoints_did)
   deallocate (exchRhoOp_did)
   deallocate (radialWeight_did)
 
end subroutine closeSCFExchCorrHDF5

end module O_SCFExchCorrHDF5
