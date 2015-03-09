module hdf5Routines
  use HDF5
  use dataStructs

  implicit none

  public

  integer :: valeDim
  integer :: numPotSites
  integer :: numKPoints
  integer :: potDim
  integer :: numPotTypes

  integer(hid_t) :: setup_fid
  integer(hid_t)  :: setup_plid
  integer         :: mdc_nelmts  ! Meta-data cache num elements.
  integer(size_t) :: rdcc_nelmts ! Raw-data chunk cache num elements.
  integer(size_t) :: rdcc_nbytes ! Raw-data chunk cache size in bytes.
  real            :: rdcc_w0     ! Raw-data chunk cache weighting parameter.
  integer(hid_t) :: atomIntgGroup_gid ! Atom and potential overlap group.
  integer(hid_t) :: atomOverlap_gid
  integer(hid_t) :: atomKEOverlap_gid
  integer(hid_t) :: atomNucOverlap_gid
  integer(hid_t) :: atomPotOverlap_gid
  integer(hid_t), allocatable, dimension (:) :: atomPotKPointOL_gid
  integer(hid_t) :: valeVale_dsid
  integer(hid_t) :: valeVale_plid
  integer (hid_t) :: valeVale_xferpid
  integer(hsize_t), dimension (2) :: atomDims
  integer(hid_t), allocatable, dimension (:)   :: atomOverlap_did
  integer(hid_t), allocatable, dimension (:)   :: atomKEOverlap_did
  integer(hid_t), allocatable, dimension (:)   :: atomNucOverlap_did
  integer(hid_t), allocatable, dimension (:,:) :: atomPotOverlap_did
  integer(hid_t) :: elecStatGroup_gid ! Electrostatic group.
  integer(hid_t) :: potTypesPot_dsid
  integer(hid_t) :: potPot_dsid
  integer(hid_t) :: pot_dsid
  integer(hid_t) :: potTypesPot_plid
  integer(hid_t) :: potPot_plid
  integer(hid_t) :: pot_plid
  integer(hid_t) :: potTypesPot_xferpid
  integer(hid_t) :: potPot_xferpid
  integer(hid_t) :: pot_xferpid
  integer(hsize_t), dimension (2) :: potDims2
  integer(hsize_t), dimension (2) :: potTypesPot
  integer(hsize_t), dimension (1) :: potDims1
  integer(hid_t) :: potAlphaOverlap_did
  integer(hid_t) :: nonLocalNeutQPot_did
  integer(hid_t) :: nonLocalNucQPot_did
  integer(hid_t) :: localNeutQPot_did
  integer(hid_t) :: localNucQPot_did
  integer(hid_t) :: nonLocalResidualQ_did
  integer(hid_t) :: coreChargeDensity_did
  integer(hid_t) :: exchCorrGroup_gid ! Exchange Correlation group.
  integer(hid_t) :: numPoints_dsid
  integer(hid_t) :: potPoints_dsid
  integer(hid_t) :: points_dsid
  integer(hid_t) :: potPoints_plid
  integer(hid_t) :: points_plid
  integer(hid_t):: potPOints_xferpid, points_xferpid
  integer(hsize_t), dimension(3) :: potPoints
  integer(hsize_t), dimension(1) :: points
  integer(hsize_t), dimension(1) :: numPoints
  integer(hid_t), allocatable, dimension(:) :: numPoints_did, exchRhoOp_did
  integer(hid_t), allocatable, dimension(:) :: radialWeight_did
  integer(hid_t) :: exchCorrOverlap_did

  contains

subroutine accessSetupHDF5(fileName)

   ! Use the HDF5 module.
   use HDF5
   
   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare local variables.
   integer :: hdferr
   character(len=*) :: filename

   valeDim     = systemVars%valeDim
   numPotSites = systemVars%numPotSites
   numKPoints  = systemVars%numKPoints
   potDim      = systemVars%potDim
   numPotTypes = systemVars%numPotTypes

   ! init hdf5 interface
   call h5open_f(hdferr)
   if (hdferr < 0) stop 'Failed to open HDF library'

   ! Create the property list for the setup hdf5 file and turn off
   !   chunk caching.
   call h5pcreate_f    (H5P_FILE_ACCESS_F,setup_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create setup plid in accessSetupHDF5.'
   call h5pget_cache_f (setup_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to get setup plid cache settings.'
   call h5pset_cache_f (setup_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to set setup plid cache settings.'


   ! Open the HDF5 file that will hold all the computed results.  This will
   !   die if the file already exists.  It opens read-only.
   call h5fopen_f (fileName,H5F_ACC_RDONLY_F,setup_fid,hdferr,&
         & setup_plid)
   if (hdferr /= 0) stop 'Failed to open setup-temp.hdf5 file.'

   ! Access the subgroups of the setup hdf5 file.
   call accessSetupIntegralHDF5 (setup_fid)
   call accessSetupExchCorrHDF5 (setup_fid)
   call accessSetupElecStatHDF5 (setup_fid)

end subroutine accessSetupHDF5

subroutine closeSetupHDF5

   ! Use the HDF5 module.
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   integer :: hdferr

   ! Close access to all the subgroups and their parts.
   call closeSetupIntegralHDF5
   call closeSetupExchCorrHDF5
   call closeSetupElecStatHDF5

   ! Close the property list.
   call h5pclose_f (setup_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close setup_plid.'

   ! Close the file.
   call h5fclose_f (setup_fid,hdferr)
   if (hdferr /= 0) stop 'Failed to close setup_fid.'

   ! Closer the HDF5 interface.
   call h5close_f(hdferr)
   if (hdferr /=0) stop 'Failed to close the HDF5 interace.'
end subroutine closeSetupHDF5

subroutine accessSetupIntegralHDF5 (setup_fid)

   ! Import any necessary definition modules.
   use HDF5

   ! Define the passed parameters.
   integer(hid_t) :: setup_fid

   ! Define local variables.
   integer :: i,j
   integer :: hdferr
   character*30 :: currentName

   ! Initialize data structure dimensions.
   atomDims(1) = 2 ! Real and imaginary are needed.
   atomDims(2) = valeDim*(valeDim+1)/2 ! Linear storage of 1/2 matrix.

   ! Open the Integral group within the setup HDF5 file.
   call h5gopen_f (setup_fid,"/atomIntgGroup",atomIntgGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom intg group'

   ! Open the subgroups within the atomIntgGroup.
   call h5gopen_f (atomIntgGroup_gid,"atomOverlap",atomOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open atom overlap group'
   call h5gopen_f (atomIntgGroup_gid,"atomKEOverlap",atomKEOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open kinetic energy overlap group'
   call h5gopen_f (atomIntgGroup_gid,"atomNucOverlap",atomNucOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to open atom nuclear overlap group'
   call h5gopen_f (atomIntgGroup_gid,"atomPotOverlap",atomPotOverlap_gid,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to open atom potential overlap group'

   ! Open the subgroups in the atomPotOverlap subgroup.  There will be one
   !   subgroup for each kpoint in the atomPotOverlap subgroup.  The other
   !   atom*Overlap groups will not have subgroups.  Instead, each kpoint will
   !   be a dataset.

   ! First, sufficient space must be allocated to hold the gid values for the
   !   subgroups that exist for each kpoint for the atomicPotential
   !   hamiltonian terms.  Also allocate space to hold the dataset IDs.
   allocate (atomPotKPointOL_gid (numKPoints))
   allocate (atomOverlap_did     (numKPoints))
   allocate (atomKEOverlap_did   (numKPoints))
   allocate (atomNucOverlap_did  (numKPoints))
   allocate (atomPotOverlap_did  (numKPoints,potDim))

   ! Then loop over the kpoints and assign a gid name equal to the kpoint # for
   !   the atomic-potentian hamiltonian terms.
   do i = 1, numKPoints
      write (currentName,fmt="(i7.7)") i
      currentName = trim (currentName)
      call h5gopen_f (atomPotOverlap_gid,currentName,atomPotKPointOL_gid(i),&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to open potential kpoint group'
   enddo


   ! Open the datasets that will be used for the subgroups of atomIntgGroup.
   do i = 1, numKPoints
      write (currentName,fmt="(i7.7)") i
      call h5dopen_f (atomOverlap_gid,currentName,atomOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open atom overlap did'
      call h5dopen_f (atomKEOverlap_gid,currentName,atomKEOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to open KE overlap did'
      call h5dopen_f (atomNucOverlap_gid,currentName,atomNucOverlap_did(i),&
            & hdferr)
      if (hdferr /= 0) stop 'Failed to open nuclear overlap did'
      do j = 1, potDim
         write (currentName,fmt="(i7.7)") j
         call h5dopen_f(atomPotKPointOL_gid(i),currentName,&
               & atomPotOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to open potential overlap did'
      enddo
   enddo

   ! Obtain the properties of the datasets that were just opened.  They are all
   !   the same and so only one copy is necessary.  (Actually, this value is
   !   not really used, but in the "close" subroutine we close this id so we
   !   should make sure to have it for both the setup and main calls to the
   !   close routine.)
   call h5dget_create_plist_f (atomOverlap_did(1),valeVale_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain vale vale plid'

   ! Obtain the dataspace that is used for each dataset in atomIntgGroup
   !   and all of its subgroups.  The same dataspace definition works for all
   !   of the datasets.  (Same as for the plist above.)
   call h5dget_space_f(atomOverlap_did(1),valeVale_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain vale vale dsid'

end subroutine accessSetupIntegralHDF5


subroutine closeSetupIntegralHDF5

   ! Import any necessary definition modules.
   use HDF5

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Declare local loop variables and error control.
   integer :: i,j
   integer :: hdferr

   ! Close the property list first.
   call h5pclose_f (valeVale_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close valeVale_plid.'

   ! Close the datasets next.
   do i = 1, numKPoints
      call h5dclose_f (atomOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close atomOverlap_did.'
      call h5dclose_f (atomKEOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close atomKEOverlap_did.'
      call h5dclose_f (atomNucOverlap_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close atomNucOverlap_did.'
      do j = 1, potDim
         call h5dclose_f (atomPotOverlap_did(i,j),hdferr)
         if (hdferr /= 0) stop 'Failed to close atomPotOverlap_did.'
      enddo
   enddo

   ! Close the data spaces next.
   call h5sclose_f (valeVale_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close valeVale_dsid.'

   ! Close the groups.
   do i = 1, numKPoints
      call h5gclose_f (atomPotKPointOL_gid(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close atomPotKPointOL_gid.'
   enddo
   call h5gclose_f (atomOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomOverlap_gid.'
   call h5gclose_f (atomKEOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomKEOverlap_gid.'
   call h5gclose_f (atomNucOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomNucOverlap_gid.'
   call h5gclose_f (atomPotOverlap_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomPotOverlap_gid.'
   call h5gclose_f (atomIntgGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close atomIntgGroup_gid.'

   ! Deallocate the arrays that hold the id numbers.
   deallocate (atomPotKPointOL_gid)
   deallocate (atomOverlap_did)
   deallocate (atomKEOverlap_did)
   deallocate (atomNucOverlap_did)
   deallocate (atomPotOverlap_did)
   
end subroutine closeSetupIntegralHDF5

subroutine accessSetupElecStatHDF5 (setup_fid)

   ! Import any necessary definition modules.
   use HDF5

   ! Define the passed parameters.
   integer(hid_t) :: setup_fid

   ! Define local variables.
   integer :: hdferr

   ! Initialize data structure dimensions.
   potDims2(1)    = potDim
   potDims2(2)    = potDim
   potTypesPot(1) = numPotTypes
   potTypesPot(2) = potDim
   potDims1(1)    = potDim

   ! Open the electrostatic group within the HDF5 file.
   call h5gopen_f (setup_fid,"/elecStatGroup",elecStatGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open electrostatics group'

   ! Open the datasets that will be used for the elecStatGroup.
   call h5dopen_f(elecStatGroup_gid,"potAlphaOverlap",potAlphaOverlap_did,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to open potAlphaOverlap did.'
   call h5dopen_f(elecStatGroup_gid,"nonLocalNeutQPot",nonLocalNeutQPot_did,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to open nonLocalNeutQPot did.'
   call h5dopen_f(elecStatGroup_gid,"nonLocalNucQPot",nonLocalNucQPot_did,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to open nonLocalNucQPot did.'
   call h5dopen_f(elecStatGroup_gid,"localNeutQPot",localNeutQPot_did,hdferr)
   if (hdferr /= 0) stop 'Failed to open localNeutQPot did.'
   call h5dopen_f(elecStatGroup_gid,"localNucQPot",localNucQPot_did,hdferr)
   if (hdferr /= 0) stop 'Failed to open localNucQPot did.'
   call h5dopen_f(elecStatGroup_gid,"nonLocalResidualQ",nonLocalResidualQ_did,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to open nonLocalResidualQ did.'
   call h5dopen_f(elecStatGroup_gid,"coreChargeDensity",coreChargeDensity_did,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to open coreChargeDensity did.'

   ! Obtain the property lists.  Read the note in accessSetupIntegralHDF5.
   call h5dget_create_plist_f (nonLocalResidualQ_did,potTypesPot_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain potTypesPot plid'
   call h5dget_create_plist_f (potAlphaOverlap_did,potPot_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain potPot plid'
   call h5dget_create_plist_f (coreChargeDensity_did,pot_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain pot plid'

   ! Obtain the dataspaces that are used for the elecStatGroup datasets.  Read
   !   the note in accessSetupINtegralHDF5.
   call h5dget_space_f (nonLocalResidualQ_did,potTypesPot_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain potTypesPot dsid'
   call h5dget_space_f (potAlphaOverlap_did,potPot_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain potPot dsid'
   call h5dget_space_f (coreChargeDensity_did,pot_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to obtain pot dsid'

end subroutine accessSetupElecStatHDF5


subroutine closeSetupElecStatHDF5

   ! Import any necessary definition modules.
   use HDF5

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Declare local loop variables and error control.
   integer :: hdferr

   ! Close all access to the electrostatic calculation parts of the HDF file
   !   for setup.
   ! Close MPI property lists
   call h5pclose_f (potTypesPot_xferpid,hdferr)
   if (hdferr /= 0) stop 'Failed to close potTypesPot_xferpid.'
   call h5pclose_f (potPot_xferpid,hdferr)
   if (hdferr /= 0) stop 'Failed to close potPot_xferpid.'
   call h5pclose_f (pot_xferpid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pot_xferpid.'

   ! Close the property lists first.
   call h5pclose_f (potTypesPot_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close potTypesPot_plid.'
   call h5pclose_f (potPot_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close potPot_plid.'
   call h5pclose_f (pot_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pot_plid.'


   ! Close the datasets next.
   call h5dclose_f (potAlphaOverlap_did,hdferr)
   if (hdferr /= 0) stop 'Failed to close potAlphaOverlap_did.'
   call h5dclose_f (nonLocalNeutQPot_did,hdferr)
   if (hdferr /= 0) stop 'Failed to close nonLocalNeutQPot_did.'
   call h5dclose_f (nonLocalNucQPot_did,hdferr)
   if (hdferr /= 0) stop 'Failed to close nonLocalNucQPot_did.'
   call h5dclose_f (localNeutQPot_did,hdferr)
   if (hdferr /= 0) stop 'Failed to close localNeutQPot_did.'
   call h5dclose_f (localNucQPot_did,hdferr)
   if (hdferr /= 0) stop 'Failed to close localNucQPot_did.'
   call h5dclose_f (nonLocalResidualQ_did,hdferr)
   if (hdferr /= 0) stop 'Failed to close nonLocalResidualQ_did.'

   ! Close the data spaces next.
   call h5sclose_f (potTypesPot_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close potTypesPot_dsid.'
   call h5sclose_f (potPot_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close potPot_dsid.'
   call h5sclose_f (pot_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close pot_dsid.'

   ! Close the group.
   call h5gclose_f (elecStatGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close elecStatGroup_gid.'

end subroutine closeSetupElecStatHDF5

subroutine accessSetupExchCorrHDF5 (setup_fid)

   ! Import any necessary definition modules.
   use HDF5

   ! Define the passed parameters.
   integer(hid_t) :: setup_fid

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

   ! Open the exchange correlation group within the setup_fid.
   call h5gopen_f (setup_fid,"/exchCorrGroup",exchCorrGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to open exchange correlation group'


   ! Allocate space to hold the IDs for the datasets in the exchCorrGroup.
   allocate (numPoints_did(numPotSites))
   allocate (exchRhoOp_did(numPotSites))
   allocate (radialWeight_did(numPotSites))

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
      call h5dopen_f(exchCorrGroup_gid,currentName,radialWeight_did(i),hdferr)
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

   ! Define the dimensions of the remaining HDF5 stored matrices.
   potPoints(1) = potDim
   potPoints(2) = tempMaxNumRayPoints
   ! For nonGGA only the rho op value is stored.
   ! For GGA the first and second derivatives, as well as
   ! the rho op value are stored.
   points(1)    = tempMaxNumRayPoints

end subroutine accessSetupExchCorrHDF5



subroutine closeSetupExchCorrHDF5

   ! Import any necessary definition modules.
   use HDF5

   ! Make sure that no variables are implicitly declared.
   implicit none

   ! Declare local loop variables and error control.
   integer :: i
   integer :: hdferr

   ! Close all access to the exchange correlation parts of the HDF file
   !   for setup.

   ! Close the property lists first.
   call h5pclose_f (potPoints_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close potPoints_plid.'
   call h5pclose_f (points_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close points_plid.'

   ! Close MPI property lists
   call h5pclose_f (potPoints_xferpid, hdferr)
   if (hdferr /= 0) stop 'Failed to close potPoints_xferpid'
   call h5pclose_f (points_xferpid, hdferr)
   if (hdferr /= 0) stop 'Failed to close points_xferpid'

   ! Close the datasets next.
   do i = 1, numPotSites
      call h5dclose_f (numPoints_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close numPoints_did.'
      call h5dclose_f (exchRhoOp_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close exchRhoOp_did.'
      call h5dclose_f (radialWeight_did(i),hdferr)
      if (hdferr /= 0) stop 'Failed to close radialWeight_did.'
   enddo

   ! Close the data spaces next.
   call h5sclose_f (potPoints_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close potPoints_dsid.'
   call h5sclose_f (points_dsid,hdferr)
   if (hdferr /= 0) stop 'Failed to close points_dsid.'

   ! Close the groups.
   call h5gclose_f (exchCorrGroup_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to close exchCorrGroup_gid.'

   ! Deallocate space that holds the IDs for the datasets in the exchCorrGroup.
   deallocate (numPoints_did)
   deallocate (exchRhoOp_did)
   deallocate (radialWeight_did)
 
end subroutine closeSetupExchCorrHDF5

end module hdf5Routines
